"""Tests for the GenomeScope reliability gate.

The fixtures here are real output from NOVA_260909_LA, because every one of these checks
exists to catch something that actually shipped: OG3043's Inf bound published as a 2.5 Mb
genome, OG2983 graded EXCELLENT on a degenerate fit, and a peak detector that invented a
mode on a monotone histogram and then flagged the disagreement.
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import genomescope_reliability as gate


def summary(haploid="358,925,911", repeat="90,607,254", unique="268,318,657",
            fit=("80.9337", "98.6498"), het=("0.375176", "0.421044")):
    return (
        "GenomeScope version 2.0\n"
        "k = 21\n\n"
        "property                      min               max\n"
        "Homozygous (aa)               99.579%           99.6248%\n"
        f"Heterozygous (ab)             {het[0]}%         {het[1]}%\n"
        f"Genome Haploid Length         {haploid.split('|')[0]} bp    {haploid.split('|')[-1]} bp\n"
        f"Genome Repeat Length          {repeat} bp     {repeat} bp\n"
        f"Genome Unique Length          {unique} bp    {unique} bp\n"
        f"Model Fit                     {fit[0]}%          {fit[1]}%\n"
        "Read Error Rate               0.220766%         0.220766%\n"
    )


def write(text, suffix=".txt"):
    handle = tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False)
    handle.write(text)
    handle.close()
    return handle.name


def peaked_histogram(kmercov=30, span=400):
    """An error spike falling into a trough, then a genomic peak at 2x kmercov."""
    hist = {}
    for cov in range(1, span + 1):
        error = int(8e8 / (cov ** 3))
        genomic = int(3e6 * 2.718 ** (-((cov - 2 * kmercov) ** 2) / (2 * (kmercov / 2) ** 2)))
        hist[cov] = error + genomic + 10
    return hist


def monotone_histogram(span=400):
    """What a coverage-starved library looks like: no mode anywhere."""
    return {cov: int(8e8 / (cov ** 2)) + 10 for cov in range(1, span + 1)}


class ParsingTests(unittest.TestCase):
    def test_reads_both_bounds(self):
        parsed = gate.parse_summary(write(summary()))
        self.assertEqual(parsed["haploid_max"], 358925911)
        self.assertAlmostEqual(parsed["model_fit_full"], 98.6498)
        self.assertAlmostEqual(parsed["model_fit_allkmers"], 80.9337)

    def test_inf_upper_bound_is_none_not_the_lower_bound(self):
        # OG3043. The shipped regex could not match "Inf bp" and silently returned the
        # lower bound, publishing 2.5 Mb as the genome size for a 322 Mb assembly.
        parsed = gate.parse_summary(write(summary(haploid="2,499,788|Inf")))
        self.assertEqual(parsed["haploid_min"], 2499788)
        self.assertIsNone(parsed["haploid_max"])

    def test_inf_bound_flags_rather_than_substitutes(self):
        parsed = gate.parse_summary(write(summary(haploid="2,499,788|Inf")))
        flags, details = gate.evaluate(parsed, {})
        self.assertIn("non_finite_genome_size", flags)
        self.assertIsNone(details["estimated_genome_size"])


class ModelFitTests(unittest.TestCase):
    def test_the_two_fit_values_are_not_bounds(self):
        # GenomeScope prints allscore then fullscore under a "min  max" header. They are
        # unrelated statistics, so the first can exceed the second -- which is what the
        # inverted-bounds flag detects, and why neither is named _min or _max here.
        parsed = gate.parse_summary(write(summary(fit=("81.7476", "25.2118"))))
        self.assertGreater(parsed["model_fit_allkmers"], parsed["model_fit_full"])

    def test_gate_uses_the_full_model_fit_not_the_all_kmers_one(self):
        # model_fit_allkmers falls as --max_kmercov admits more unmodelled tail. OG2906
        # read 69.9% at -m 10000 while model_fit_full stayed 98.6498% at every cutoff.
        # Gating on the all-k-mers figure flagged it anyway.
        parsed = gate.parse_summary(write(summary(fit=("69.8861", "98.6498"))))
        flags, _ = gate.evaluate(parsed, {"kmercov": 29.7}, peaked_histogram())
        self.assertNotIn("model_fit_below_50pc", flags)

    def test_genuinely_bad_fit_is_flagged(self):
        parsed = gate.parse_summary(write(summary(fit=("49.1132", "33.7"))))
        flags, _ = gate.evaluate(parsed, {"kmercov": 29.7}, peaked_histogram())
        self.assertIn("model_fit_below_50pc", flags)

    def test_inverted_bounds_flagged(self):
        parsed = gate.parse_summary(write(summary(fit=("81.7476", "25.2118"))))
        flags, _ = gate.evaluate(parsed, {})
        self.assertIn("inverted_model_fit_bounds", flags)

    def test_a_good_fit_to_few_kmers_is_still_reported(self):
        # OG2644 is stored as an 84% fit while the model accounts for only 17% of its
        # k-mers. The gate passes it on fit, but the number has to reach the caller.
        parsed = gate.parse_summary(write(summary(fit=("17.0", "84.0"))))
        _flags, details = gate.evaluate(parsed, {"kmercov": 30.0}, peaked_histogram())
        self.assertAlmostEqual(details["model_fit_allkmers"], 17.0)
        self.assertAlmostEqual(details["model_fit_full"], 84.0)

    def test_large_bias_flagged(self):
        # OG2983: bias 37 at kmercov 569 for a 10.9 Mb "genome", graded EXCELLENT.
        parsed = gate.parse_summary(write(summary()))
        flags, _ = gate.evaluate(parsed, {"kmercov": 569.1, "bias": 37.0}, peaked_histogram())
        self.assertIn("implausible_model_bias", flags)

    def test_zero_lower_heterozygosity_is_degenerate(self):
        parsed = gate.parse_summary(write(summary(het=("0", "14.6879"))))
        flags, _ = gate.evaluate(parsed, {})
        self.assertIn("degenerate_model_fit", flags)


class PeakDetectionTests(unittest.TestCase):
    def test_finds_the_genomic_mode(self):
        peak, reason, _trough = gate.find_kmer_peak(peaked_histogram(kmercov=30))
        self.assertEqual(reason, "ok")
        self.assertTrue(50 <= peak <= 70, peak)

    def test_monotone_histogram_has_no_peak(self):
        # The old detector invented one out in the sparse tail and then flagged the
        # disagreement it had just created, on 22 of 23 tier-4 samples.
        peak, reason, _trough = gate.find_kmer_peak(monotone_histogram())
        self.assertIsNone(peak)
        self.assertEqual(reason, "no_trough_histogram_monotone")

    def test_no_peak_is_its_own_flag(self):
        parsed = gate.parse_summary(write(summary()))
        flags, _ = gate.evaluate(parsed, {"kmercov": 30.0}, monotone_histogram())
        self.assertIn("no_kmer_peak", flags)
        self.assertNotIn("kmer_peak_disagrees_with_fitted_coverage", flags)

    def test_peak_at_twice_kmercov_is_accepted(self):
        # The homozygous mode sits at 2x the fitted haploid coverage; that is normal.
        parsed = gate.parse_summary(write(summary()))
        flags, _ = gate.evaluate(parsed, {"kmercov": 30.0}, peaked_histogram(kmercov=30))
        self.assertNotIn("kmer_peak_disagrees_with_fitted_coverage", flags)

    def test_fitted_coverage_far_from_the_mode_is_flagged(self):
        parsed = gate.parse_summary(write(summary()))
        flags, _ = gate.evaluate(parsed, {"kmercov": 2.0}, peaked_histogram(kmercov=30))
        self.assertIn("kmer_peak_disagrees_with_fitted_coverage", flags)

    def test_peak_search_is_independent_of_max_kmercov(self):
        # Reading the histogram at two different ceilings must not move the answer; tying
        # the search range to genomescope2_m meant changing the cutoff changed the flag.
        path = write("\n".join(f"{c}\t{n}" for c, n in sorted(peaked_histogram().items())), ".hist")
        narrow = gate.find_kmer_peak(gate.read_histogram(path, 1000))
        wide = gate.find_kmer_peak(gate.read_histogram(path, 10000))
        self.assertEqual(narrow[0], wide[0])


class AssemblyCrossCheckTests(unittest.TestCase):
    def test_estimate_far_below_assembly_is_flagged(self):
        # OG3037: 209 Mb estimated at an 86% fit, against a 1.26 Gb assembly.
        parsed = gate.parse_summary(write(summary(haploid="208,969,898|208,969,898")))
        flags, _ = gate.evaluate(parsed, {"kmercov": 39.1}, peaked_histogram(kmercov=39),
                                 assembly_size=1_258_000_000)
        self.assertIn("genome_size_disagrees_with_assembly", flags)

    def test_estimate_near_assembly_passes(self):
        parsed = gate.parse_summary(write(summary()))
        flags, _ = gate.evaluate(parsed, {"kmercov": 30.0}, peaked_histogram(kmercov=30),
                                 assembly_size=375_000_000)
        self.assertEqual(flags, [])

    def test_recheck_merges_into_an_existing_summary(self):
        coverage = {"estimated_genome_size": 208_969_898, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT"}
        stats = write("file\tformat\ttype\tnum_seqs\tsum_len\n"
                      "a.fa\tFASTA\tDNA\t1000\t1258000000\n", ".tsv")
        merged = gate.recheck_against_assembly(
            coverage, gate.assembly_size_from_seqkit(stats))
        self.assertEqual(merged["assembly_size"], 1258000000)
        self.assertIn("genome_size_disagrees_with_assembly", merged["genome_size_flags"])
        self.assertFalse(merged["genome_size_reliable"])
        self.assertEqual(merged["coverage_status"], "UNRELIABLE_GENOME_SIZE_ESTIMATE")

    def test_recheck_keeps_existing_flags(self):
        coverage = {"estimated_genome_size": 400_000_000,
                    "genome_size_flags": "no_kmer_peak", "genome_size_reliable": False,
                    "coverage_status": "UNRELIABLE_GENOME_SIZE_ESTIMATE"}
        merged = gate.recheck_against_assembly(coverage, 400_000_000)
        self.assertEqual(merged["genome_size_flags"], "no_kmer_peak")

    def test_recheck_without_an_assembly_size_records_that(self):
        coverage = {"estimated_genome_size": 400_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT"}
        merged = gate.recheck_against_assembly(coverage, None)
        self.assertEqual(merged["assembly_cross_check"], "no_assembly_size")
        self.assertTrue(merged["genome_size_reliable"])


class CliTests(unittest.TestCase):
    def test_emits_json(self):
        import io
        import contextlib
        path = write(summary())
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            gate.main(["--summary", path])
        parsed = json.loads(buffer.getvalue())
        self.assertEqual(parsed["estimated_genome_size"], 358925911)


if __name__ == "__main__":
    unittest.main()


class ReseedDecisionTests(unittest.TestCase):
    """Both halves of the rule, each fixed before the sweep was run."""

    def fit(self, model_fit_full, kmercov):
        return {"summary": {"model_fit_full": model_fit_full}, "model": {"kmercov": kmercov}}

    def test_keeps_a_reseed_that_improves_fit_and_moves_toward_the_measurement(self):
        # OG2617: 83.5 -> 96.7 with kmercov 17.9 -> 8.3 against a measured 9.2.
        winner, reason = gate.choose_reseeded_fit(
            self.fit(83.5, 17.9), self.fit(96.7, 8.3), 9.2)
        self.assertEqual(winner, "reseeded")
        self.assertIn("96.7", reason)

    def test_rejects_a_reseed_that_hits_the_measurement_but_degrades_the_fit(self):
        # OG2672: kmercov 16.0 -> 7.7 against a measured 6.3, but fit 80.9 -> 62.8. That
        # is a different optimum, not a better one.
        winner, _reason = gate.choose_reseeded_fit(
            self.fit(80.9, 16.0), self.fit(62.8, 7.7), 6.3)
        self.assertEqual(winner, "original")

    def test_rejects_a_reseed_that_improves_fit_while_moving_away(self):
        winner, reason = gate.choose_reseeded_fit(
            self.fit(80.0, 10.0), self.fit(95.0, 40.0), 9.2)
        self.assertEqual(winner, "original")
        self.assertIn("no closer", reason)

    def test_rejects_a_marginal_gain(self):
        winner, reason = gate.choose_reseeded_fit(
            self.fit(90.0, 20.0), self.fit(90.5, 10.0), 9.2)
        self.assertEqual(winner, "original")
        self.assertIn("below the 2", reason)

    def test_without_a_measured_lambda_there_is_nothing_to_judge_against(self):
        winner, _reason = gate.choose_reseeded_fit(
            self.fit(80.0, 17.9), self.fit(96.0, 8.3), None)
        self.assertEqual(winner, "original")


class DepthCrossCheckTests(unittest.TestCase):
    def test_fitted_coverage_far_from_measured_is_flagged(self):
        coverage = {"estimated_genome_size": 123_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT"}
        merged = gate.recheck_against_assembly(
            coverage, 400_000_000, None,
            {"lambda_depth": 2.5, "genomescope_kmercov": 39.1})
        self.assertIn("fitted_coverage_disagrees_with_read_depth", merged["genome_size_flags"])

    def test_host_minority_assembly_is_flagged(self):
        # OG3037: host single-copy genes at 5.8x against an assembly-wide 31.6x. Not a fit
        # problem -- the assembly is not one organism.
        coverage = {"estimated_genome_size": 123_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT"}
        merged = gate.recheck_against_assembly(
            coverage, 400_000_000, None,
            {"lambda_depth": 30.0, "genomescope_kmercov": 30.0,
             "busco_depth_over_assembly_depth": 0.18})
        self.assertIn("assembly_not_dominated_by_host", merged["genome_size_flags"])

    def test_low_host_coverage_says_so_directly(self):
        coverage = {"estimated_genome_size": 123_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT"}
        merged = gate.recheck_against_assembly(
            coverage, 400_000_000, None,
            {"lambda_depth": 2.5, "genomescope_kmercov": 2.5})
        self.assertIn("host_coverage_too_low", merged["genome_size_flags"])

    def test_rechecking_a_recheck_changes_nothing(self):
        # RECHECK_GENOME_SIZE republishes the summary it was given, so on a resumed run its
        # provisional input can be its own earlier output. Every flag NOVA_260909_LA's
        # tier-4 samples carry was stored twice because the appends were unguarded.
        coverage = {"estimated_genome_size": 114_000_000, "genome_size_flags": "no_kmer_peak",
                    "genome_size_reliable": False,
                    "coverage_status": "UNRELIABLE_GENOME_SIZE_ESTIMATE"}
        depth = {"lambda_depth": 3.8, "genomescope_kmercov": 28.9,
                 "busco_depth_over_assembly_depth": 0.15}
        once = gate.recheck_against_assembly(coverage, 327_693_880, None, depth)
        twice = gate.recheck_against_assembly(once, 327_693_880, None, depth)
        self.assertEqual(once["genome_size_flags"], twice["genome_size_flags"])
        self.assertEqual(
            once["genome_size_flags"].split("; "),
            ["no_kmer_peak", "fitted_coverage_disagrees_with_read_depth",
             "assembly_not_dominated_by_host", "host_coverage_too_low"])

    def test_agreement_adds_no_flags(self):
        coverage = {"estimated_genome_size": 359_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT"}
        merged = gate.recheck_against_assembly(
            coverage, 375_000_000, None,
            {"lambda_depth": 30.2, "genomescope_kmercov": 29.7,
             "busco_depth_over_assembly_depth": 0.75})
        self.assertEqual(merged["genome_size_flags"], "")
        self.assertTrue(merged["genome_size_reliable"])


class HostSizeCrossCheckTests(unittest.TestCase):
    """The ratio's denominator. GenomeScope estimates the HOST genome, so measuring it
    against a total assembly that can be half symbiont puts the wrong number underneath
    the fraction in exactly the samples the check exists for."""

    # Measured with samtools coverage, not the idxstats approximation the design note
    # used: that credits a contig the full length of every read recorded against it, and
    # only 55% of a mapped read's bases align on this sample.
    OG2634 = {"lambda_depth": 4.5, "genomescope_kmercov": 4.5,
              "busco_depth_over_assembly_depth": 0.6,
              "host_assembly_size": 207_000_000, "host_assembly_fraction": 0.632,
              "symbiont_assembly_size": 44_700_000, "partition_status": "ok"}

    def base(self, size):
        return {"estimated_genome_size": size, "genome_size_flags": "",
                "genome_size_reliable": True, "coverage_status": "EXCELLENT",
                "kmercov": 4.5}

    def test_host_mass_is_the_denominator_when_the_partition_worked(self):
        # 103.5 Mb against a 328 Mb total assembly is 0.32 and passes on a technicality;
        # against the 207 Mb that is actually the animal it is 0.50, and the flag that
        # matters is neither -- what changes is that the ratio now means something.
        merged = gate.recheck_against_assembly(
            self.base(103_500_000), 328_000_000, None, self.OG2634)
        self.assertEqual(merged["assembly_ratio_basis"], "host")
        self.assertEqual(merged["assembly_ratio_denominator"], 207_000_000)
        self.assertEqual(merged["host_assembly_size"], 207_000_000)
        self.assertAlmostEqual(merged["host_assembly_fraction"], 0.632)

    def test_total_assembly_is_the_denominator_without_a_partition(self):
        merged = gate.recheck_against_assembly(
            self.base(103_500_000), 328_000_000, None,
            {"lambda_depth": 4.5, "partition_status": "no_host_depth",
             "host_assembly_size": None})
        self.assertEqual(merged["assembly_ratio_basis"], "total")
        self.assertEqual(merged["assembly_ratio_denominator"], 328_000_000)

    def test_a_host_minority_assembly_no_longer_hides_a_disagreement(self):
        # 40 Mb of host mass against a 400 Mb assembly: 0.10 of the host, 0.01 of the
        # total. Both flag, but only the host ratio is the quantity GenomeScope estimated.
        merged = gate.recheck_against_assembly(
            self.base(4_000_000), 400_000_000, None,
            dict(self.OG2634, host_assembly_size=40_000_000))
        self.assertIn("genome_size_disagrees_with_assembly", merged["genome_size_flags"])
        self.assertEqual(merged["assembly_ratio_basis"], "host")


class KmercovMergeTests(unittest.TestCase):
    """kmercov was null for 16 NOVA_260909_LA rows, all from --skip_genome_assembly tiers
    where CALCULATE_SEQUENCING_COVERAGE never runs, while the measurement parsed the same
    model file successfully every time."""

    def test_measured_kmercov_fills_an_absent_summary_value(self):
        coverage = {"estimated_genome_size": 359_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT"}
        merged = gate.recheck_against_assembly(
            coverage, 375_000_000, None,
            {"lambda_depth": 30.2, "genomescope_kmercov": 36.0,
             "busco_depth_over_assembly_depth": 0.75})
        self.assertAlmostEqual(merged["kmercov"], 36.0)

    def test_the_summarys_own_value_wins(self):
        coverage = {"estimated_genome_size": 359_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT",
                    "kmercov": 29.7}
        merged = gate.recheck_against_assembly(
            coverage, 375_000_000, None,
            {"lambda_depth": 30.2, "genomescope_kmercov": 36.0,
             "busco_depth_over_assembly_depth": 0.75})
        self.assertAlmostEqual(merged["kmercov"], 29.7)

    def test_the_legacy_spelling_is_still_read(self):
        coverage = {"estimated_genome_size": 359_000_000, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT",
                    "fitted_kmer_coverage": 29.7}
        merged = gate.recheck_against_assembly(
            coverage, 375_000_000, None,
            {"lambda_depth": 30.2, "busco_depth_over_assembly_depth": 0.75})
        self.assertAlmostEqual(merged["kmercov"], 29.7)


class ReseedPlanTests(unittest.TestCase):
    """Whether a reseed is attempted at all.

    The gate the plan reads is the PROVISIONAL verdict, computed before the assembly
    exists, so it can see neither the assembly cross-check nor the lambda ratio. OG2617
    was recorded as "no reseed attempted" on a provisional pass while its fitted coverage
    was 1.95x the measurement.
    """

    def plan(self, coverage, depth, **kwargs):
        import contextlib
        import io
        import genomescope_reseed_decide as decide
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "cov.json").write_text(json.dumps(coverage))
            (root / "depth.json").write_text(json.dumps(depth))
            argv = ["plan", "--coverage-json", str(root / "cov.json"),
                    "--kmer-depth-json", str(root / "depth.json")]
            for key, value in kwargs.items():
                argv += [f"--{key.replace('_', '-')}", str(value)]
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                decide.main(argv)
        return buf.getvalue().strip()

    def test_a_provisional_pass_with_a_disagreeing_ratio_is_still_reseeded(self):
        # OG2617: provisional verdict reliable, fitted/measured 1.954.
        self.assertEqual(
            self.plan({"genome_size_reliable": True},
                      {"lambda_depth": 9.2, "kmercov_over_lambda_depth": 1.954}),
            "9")

    def test_a_provisional_pass_that_agrees_is_left_alone(self):
        self.assertEqual(
            self.plan({"genome_size_reliable": True},
                      {"lambda_depth": 30.2, "kmercov_over_lambda_depth": 0.983}),
            "skip")

    def test_the_ratio_is_symmetric(self):
        # A fitted coverage at half the measurement is the same factor-of-two error as one
        # at twice it -- GenomeScope's characteristic failure, in either direction.
        self.assertEqual(
            self.plan({"genome_size_reliable": True},
                      {"lambda_depth": 20.0, "kmercov_over_lambda_depth": 0.5}),
            "20")

    def test_a_flagged_fit_below_the_seedable_floor_is_still_skipped(self):
        self.assertEqual(
            self.plan({"genome_size_reliable": False},
                      {"lambda_depth": 3.8, "kmercov_over_lambda_depth": 7.6}),
            "skip")
