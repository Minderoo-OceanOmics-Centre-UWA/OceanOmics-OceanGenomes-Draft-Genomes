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


class BuscoFloorTests(unittest.TestCase):
    """The gene-space floor. A size estimate that agrees with an assembly says nothing
    about the host when the assembly is not the host's."""

    def base(self):
        return {"estimated_genome_size": 15_546_537, "genome_size_flags": "",
                "genome_size_reliable": True, "coverage_status": "EXCELLENT"}

    def test_a_failed_assembly_that_matches_its_estimate_is_flagged(self):
        # OG3080: 15.5 Mb estimate, 8.9 Mb assembly (ratio 1.74, inside the band), C:0.0%,
        # too few Complete BUSCOs for any depth measurement. Graded EXCELLENT before.
        merged = gate.recheck_against_assembly(
            self.base(), 8_914_437, None,
            {"status": "too_few_busco_regions", "partition_status": "no_contig_depth"},
            busco_complete=0.0)
        self.assertEqual(merged["genome_size_flags"], "assembly_lacks_host_gene_space")
        self.assertEqual(merged["coverage_status"], "UNRELIABLE_GENOME_SIZE_ESTIMATE")
        self.assertEqual(merged["busco_complete"], 0.0)

    def test_the_floor_is_inclusive(self):
        at = gate.recheck_against_assembly(self.base(), 15_000_000, None, None, 10.0)
        below = gate.recheck_against_assembly(self.base(), 15_000_000, None, None, 9.9)
        self.assertTrue(at["genome_size_reliable"])
        self.assertIn("assembly_lacks_host_gene_space", below["genome_size_flags"])

    def test_no_busco_result_is_not_a_failure(self):
        merged = gate.recheck_against_assembly(self.base(), 15_000_000, None, None, None)
        self.assertTrue(merged["genome_size_reliable"])
        self.assertIsNone(merged["busco_complete"])

    def test_completeness_is_read_from_the_short_summary_json(self):
        path = write(json.dumps({"results": {
            "one_line_summary": "C:0.9%[S:0.7%,D:0.1%],F:4.2%,M:94.9%,n:672",
            "Complete percentage": 0.9}}), suffix=".json")
        self.assertEqual(gate.busco_complete_from_summary(path), 0.9)

    def test_falls_back_to_the_one_line_summary(self):
        path = write(json.dumps({"results": {"one_line_summary": "C:73.8%[S:72.0%,D:1.8%]"}}),
                     suffix=".json")
        self.assertEqual(gate.busco_complete_from_summary(path), 73.8)

    def test_an_unreadable_summary_is_none(self):
        self.assertIsNone(gate.busco_complete_from_summary("/nonexistent/summary.json"))

    def test_cli_flags_from_the_summary_file(self):
        import contextlib
        import io
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "cov.json").write_text(json.dumps(self.base()))
            (root / "seqkit.tsv").write_text(
                "file\tformat\ttype\tnum_seqs\tsum_len\nx.fa\tFASTA\tDNA\t13692\t8914437\n")
            (root / "busco.json").write_text(json.dumps({"results": {"Complete percentage": 0.0}}))
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                gate.main(["--recheck-assembly", "--coverage-json", str(root / "cov.json"),
                           "--seqkit-stats", str(root / "seqkit.tsv"),
                           "--busco-summary", str(root / "busco.json")])
        merged = json.loads(buffer.getvalue())
        self.assertIn("assembly_lacks_host_gene_space", merged["genome_size_flags"])


class PartitionConfidenceGateTests(unittest.TestCase):
    """Which host size the ratio is taken against, now that a partition can say it is
    not to be trusted."""

    def base(self, size):
        return {"estimated_genome_size": size, "genome_size_flags": "",
                "genome_size_reliable": True, "coverage_status": "EXCELLENT"}

    def test_a_low_confidence_partition_falls_back_to_the_total_assembly(self):
        # OG2942: 60 BUSCO regions, anchor 4.4x the assembly's depth mode, "host" 190 of
        # 1,458 Mb.
        merged = gate.recheck_against_assembly(
            self.base(114_000_000), 1_457_700_000, None,
            {"partition_status": "ok", "host_assembly_size": 190_223_452,
             "host_haploid_size": 147_000_000, "partition_confidence": "low",
             "partition_note": "anchor_off_mode",
             "host_assembly_size_mode_anchored": 1_165_000_000})
        self.assertEqual(merged["assembly_ratio_basis"], "total")
        self.assertEqual(merged["assembly_ratio_denominator"], 1_457_700_000)
        self.assertEqual(merged["partition_note"], "anchor_off_mode")
        self.assertEqual(merged["host_assembly_size_mode_anchored"], 1_165_000_000)

    def test_the_haploid_host_size_is_the_denominator_when_present(self):
        # OG3074: 494 Mb counting both haplotypes of split loci, 428 Mb counting them once.
        merged = gate.recheck_against_assembly(
            self.base(338_000_000), 567_000_000, None,
            {"partition_status": "ok", "host_assembly_size": 493_804_448,
             "host_haploid_size": 428_000_000, "partition_confidence": "high"})
        self.assertEqual(merged["assembly_ratio_basis"], "host")
        self.assertEqual(merged["assembly_ratio_denominator"], 428_000_000)
        self.assertTrue(merged["genome_size_reliable"])

    def test_a_partition_from_before_confidence_existed_is_trusted_as_before(self):
        merged = gate.recheck_against_assembly(
            self.base(103_500_000), 328_000_000, None,
            {"partition_status": "ok", "host_assembly_size": 207_000_000})
        self.assertEqual(merged["assembly_ratio_basis"], "host")
        self.assertEqual(merged["assembly_ratio_denominator"], 207_000_000)


class StaleSizeTests(unittest.TestCase):
    """A provisional summary from an older parse can reach the recheck as a precomputed
    input. GenomeScope's failure value is "-1 bp", which that parser read as 1 bp."""

    def test_a_one_bp_genome_is_no_estimate(self):
        # OG2653 as published: 1 bp, 18,007,370,181x, flagged only for disagreeing with
        # the assembly.
        coverage = {"estimated_genome_size": 1, "genome_size_flags": "",
                    "genome_size_reliable": True, "coverage_status": "EXCELLENT",
                    "theoretical_coverage": 18007370181.0, "coverage_before_filtering": 2e10}
        merged = gate.recheck_against_assembly(coverage, 185_000_000)
        self.assertEqual(merged["genome_size_flags"], "no_genome_size_estimate")
        self.assertEqual(merged["estimated_genome_size"], 0)
        self.assertEqual(merged["theoretical_coverage"], 0)
        self.assertEqual(merged["coverage_before_filtering"], 0)
        self.assertEqual(merged["coverage_status"], "UNRELIABLE_GENOME_SIZE_ESTIMATE")

    def test_a_non_finite_size_keeps_its_own_flag(self):
        coverage = {"estimated_genome_size": 0, "genome_size_flags": "non_finite_genome_size",
                    "genome_size_reliable": False}
        merged = gate.recheck_against_assembly(coverage, 185_000_000)
        self.assertEqual(merged["genome_size_flags"], "non_finite_genome_size")

    def test_the_current_parser_reads_minus_one_as_a_failure(self):
        parsed = gate.parse_summary(write(summary(haploid="-1|-1")))
        flags, details = gate.evaluate(parsed, {})
        self.assertIn("no_genome_size_estimate", flags)


class TextReportTests(unittest.TestCase):
    """The text report is rendered from the summary, so it cannot disagree with it."""

    def provisional(self):
        return {"sample_id": "OG2617", "estimated_genome_size": 363_338_117,
                "genome_size_min": 351_373_409, "genome_size_max": 363_338_117,
                "model_fit_allkmers": 48.4, "model_fit_full": 83.5,
                "heterozygosity_min": 3.9, "heterozygosity_max": 12.4,
                "total_reads_after_filtering": 163_999_458,
                "total_bases_after_filtering": 24_362_957_504,
                "theoretical_coverage": 67.05, "coverage_before_filtering": 77.8,
                "data_retention_percent": 86.2, "kmer_peak_coverage": 13, "kmercov": 17.92,
                "genome_size_reliable": True, "genome_size_flags": "",
                "coverage_status": "EXCELLENT",
                "coverage_recommendation": "Suitable for high-quality genome assembly.",
                "assembly_cross_check": "pending"}

    def test_provisional_is_labelled(self):
        text = gate.render_text_report(self.provisional())
        self.assertIn("PROVISIONAL", text)
        self.assertIn("Status: EXCELLENT coverage (67.0x)", text)
        self.assertNotIn("ASSEMBLY CROSS-CHECK", text)

    def test_rechecked_report_follows_the_flags(self):
        # OG2617: k-mer checks passed, then the read-depth check flagged it. The old text
        # file still said "reliable: yes / EXCELLENT".
        merged = gate.recheck_against_assembly(
            self.provisional(), 829_293_888, None,
            {"lambda_depth": 9.17, "partition_status": "ok",
             "host_assembly_size": 709_148_311, "host_assembly_fraction": 0.855})
        text = gate.render_text_report(merged)
        self.assertIn("final (assembly cross-check done)", text)
        self.assertIn("Genome size estimate reliable: no", text)
        self.assertIn("fitted_coverage_disagrees_with_read_depth", text)
        self.assertIn("UNRELIABLE_GENOME_SIZE_ESTIMATE (theoretical coverage 67.0x, not trustworthy)", text)
        self.assertNotIn("EXCELLENT", text)
        self.assertIn("Host fraction of assembly (lower bound): 85.5%", text)

    def test_renders_a_summary_with_fields_missing(self):
        text = gate.render_text_report({"sample_id": "X", "assembly_cross_check": "done"})
        self.assertIn("Estimated genome size: none", text)

    def test_cli_renders_from_a_file(self):
        import contextlib
        import io
        path = write(json.dumps(self.provisional()), suffix=".json")
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            gate.main(["--render-report", path])
        self.assertEqual(buffer.getvalue(), gate.render_text_report(self.provisional()))


class ReseededFitTests(unittest.TestCase):
    """Where the refit won, the verdict has to describe the refit, not the fit it replaced."""

    def provisional(self):
        # OG3057's original fit: 1,858 Mb at kmercov 10.2, flagged at the k-mer stage.
        return {"sample_id": "OG3057", "estimated_genome_size": 1_858_105_218,
                "genome_size_min": 1_858_105_218, "genome_size_max": 1_858_105_218,
                "kmercov": 10.2, "genome_size_flags": "degenerate_model_fit; implausible_model_bias",
                "genome_size_reliable": False, "coverage_status": "UNRELIABLE_GENOME_SIZE_ESTIMATE",
                "total_bases_after_filtering": 30_000_000_000, "data_retention_percent": 80.0,
                "theoretical_coverage": 16.1, "assembly_cross_check": "pending"}

    def refit(self):
        # The refit that won: 992 Mb at kmercov 16.7, against a measured 16.6.
        return (gate.parse_summary(write(summary(haploid="991,669,221|991,669,221",
                                                 fit=("60.0", "97.6")))),
                {"kmercov": 16.7, "bias": 1.2}, peaked_histogram(kmercov=17))

    def test_refit_replaces_the_fit_fields(self):
        merged = gate.apply_reseeded_fit(self.provisional(), *self.refit())
        self.assertEqual(merged["genomescope_fit"], "reseeded")
        self.assertEqual(merged["estimated_genome_size"], 991_669_221)
        self.assertEqual(merged["kmercov"], 16.7)
        self.assertEqual(merged["fitted_kmer_coverage"], 16.7)
        self.assertEqual(merged["model_fit_full"], 97.6)
        self.assertEqual(merged["genome_size_flags"], "")
        self.assertAlmostEqual(merged["theoretical_coverage"], 30_000_000_000 / 991_669_221)
        self.assertAlmostEqual(merged["coverage_before_filtering"], 37_500_000_000 / 991_669_221)
        self.assertEqual(merged["coverage_status"], "GOOD")

    def test_read_depth_check_is_against_the_refit(self):
        kmer_depth = {"lambda_depth": 16.6}
        original = gate.recheck_against_assembly(self.provisional(), 1_361_000_000, None, kmer_depth)
        self.assertIn("fitted_coverage_disagrees_with_read_depth", original["genome_size_flags"])
        refit = gate.recheck_against_assembly(
            gate.apply_reseeded_fit(self.provisional(), *self.refit()), 1_361_000_000, None, kmer_depth)
        self.assertNotIn("fitted_coverage_disagrees_with_read_depth", refit["genome_size_flags"])
        self.assertTrue(refit["genome_size_reliable"])
        self.assertEqual(refit["genomescope_fit"], "reseeded")

    def test_winning_the_reseed_does_not_bypass_the_gate(self):
        parsed = gate.parse_summary(write(summary(haploid="991,669,221|991,669,221",
                                                  het=("0", "4.1"))))
        merged = gate.apply_reseeded_fit(self.provisional(), parsed, {"kmercov": 16.7, "bias": 9.0},
                                         peaked_histogram(kmercov=17))
        self.assertFalse(merged["genome_size_reliable"])
        self.assertIn("degenerate_model_fit", merged["genome_size_flags"])
        self.assertIn("implausible_model_bias", merged["genome_size_flags"])

    def run_cli(self, winner):
        import contextlib
        import io
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "cov.json").write_text(json.dumps(self.provisional()))
            (root / "seqkit.tsv").write_text(
                "file\tformat\ttype\tnum_seqs\tsum_len\nx.fa\tFASTA\tDNA\t838185\t1361000000\n")
            (root / "decision.json").write_text(json.dumps({"winner": winner}))
            (root / "chosen_summary.txt").write_text(summary(haploid="991,669,221|991,669,221",
                                                             fit=("60.0", "97.6")))
            (root / "chosen_model.txt").write_text("kmercov 1.670e+01 1.0 1.0 0.1\nbias 1.2 0.1 1 0.1\n")
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                gate.main(["--recheck-assembly", "--coverage-json", str(root / "cov.json"),
                           "--seqkit-stats", str(root / "seqkit.tsv"),
                           "--reseed-decision", str(root / "decision.json"),
                           "--chosen-summary", str(root / "chosen_summary.txt"),
                           "--chosen-model", str(root / "chosen_model.txt")])
        return json.loads(buffer.getvalue())

    def test_cli_applies_a_winning_refit(self):
        merged = self.run_cli("reseeded")
        self.assertEqual(merged["genomescope_fit"], "reseeded")
        self.assertEqual(merged["estimated_genome_size"], 991_669_221)
        self.assertIn("Fit: reseeded", gate.render_text_report(merged))

    def test_cli_leaves_the_original_alone(self):
        merged = self.run_cli("original")
        self.assertEqual(merged["genomescope_fit"], "original")
        self.assertEqual(merged["estimated_genome_size"], 1_858_105_218)
        self.assertEqual(merged["kmercov"], 10.2)
