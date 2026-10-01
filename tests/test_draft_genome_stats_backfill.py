import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import backfill_draft_genome_stats as backfill
import draft_genome_stats as stats


def write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def make_archive(root, og_id="OG1", date="250101"):
    sample = root / og_id
    fastp = {
        "filtering_result": {"passed_filter_reads": 90, "low_quality_reads": 10},
        "summary": {
            "before_filtering": {"total_reads": 100, "total_bases": 1000, "q20_bases": 900,
                                 "q30_bases": 800, "q20_rate": 0.9, "q30_rate": 0.8,
                                 "read1_mean_length": 100, "read2_mean_length": 100, "gc_content": 0.4},
            "after_filtering": {"total_reads": 90, "total_bases": 900, "q20_bases": 850,
                                "q30_bases": 750, "q20_rate": 0.94, "q30_rate": 0.83,
                                "read1_mean_length": 100, "read2_mean_length": 100, "gc_content": 0.41},
        },
    }
    write(sample / "fastp" / f"{og_id}.ilmn.NOVA_{date}_AD.fastp.json", json.dumps(fastp))
    write(sample / "kmers" / f"{og_id}.ilmn.{date}_genomescope" / f"{og_id}.ilmn.{date}_summary.txt", f"""
name prefix = {og_id}.ilmn.{date}
Homozygous (aa)  90%  91%
Heterozygous (ab)  9%  10%
Genome Haploid Length  999 bp  1000 bp
Genome Repeat Length  99 bp  100 bp
Genome Unique Length  899 bp  900 bp
Model Fit  98%  99%
Read Error Rate  1%  2%
""")
    write(sample / "assemblies/genome/tiara" / f"{og_id}.ilmn.{date}.tiara_filter_summary.txt",
          "Category\tnum_contigs\tbp\nMitochondrion\t1\t10\nPlastid\t2\t20\nProkarya\t3\t30\n")
    write(sample / "assemblies/genome/NCBI" / f"{og_id}.ilmn.{date}.filter_report.txt",
          "EXCLUDE 4 40\nTRIM 5 50\nREVIEW 6 60\n")
    write(sample / "assemblies/genome/NCBI" / f"{og_id}.ilmn.{date}.contig_count_500bp.txt",
          "Number of contigs less than 500bp: 7\n")
    busco = {"parameters": {"in": f"{og_id}.ilmn.{date}.assembly.fna"}, "results": {
        "Complete percentage": "95.0", "Single copy percentage": "90.0",
        "Multi copy percentage": "5.0", "Fragmented percentage": "2.0",
        "Missing percentage": "3.0", "n_markers": 100, "domain": "eukaryota",
        "number_of_scaffolds": 12, "number_of_contigs": 13, "total_length": 1000,
        "percent_gaps": 0.1, "scaffold_n50": 500, "contigs_n50": 400,
        "internal_stop_codon_count": 0, "internal_stop_codon_percent": 0.0,
    }}
    write(sample / "assemblies/genome/busco" / f"{og_id}.ilmn.{date}.busco.short_summary.json", json.dumps(busco))
    write(sample / "kmers" / f"{og_id}.ilmn.{date}.merqury.completeness.stats",
          f"{og_id}.ilmn.{date}\tread_set\t80\t100\t80.0\n")
    write(sample / "kmers" / f"{og_id}.ilmn.{date}.v129mh.merqury.qv",
          f"{og_id}.ilmn.{date}\t70\t100\t30.0\t0.001\n")
    write(sample / "assemblies/genome/gfastats" / f"{og_id}.ilmn.{date}.assembly_summary.txt", """
# scaffolds: 12
Total scaffold length: 1000
Scaffold N50: 500
Largest scaffold: 600
# contigs: 13
Contig N50: 400
GC content %: 41.5
""")
    return sample


class ParserAndDiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        make_archive(self.root)

    def tearDown(self):
        self.temp.cleanup()

    def test_discovers_and_parses_all_six_families(self):
        sample = backfill.discover_sample(self.root, "OG1", "250101")
        self.assertTrue(sample.valid, sample.errors)
        self.assertEqual(set(stats.FAMILY_COLUMNS), set(sample.records))
        self.assertEqual(sample.records["fastp"]["total_reads"], 90)
        self.assertEqual(sample.records["assembly"]["genomesize"], 1000)
        self.assertEqual(sample.records["decontamination"]["bp_review"], 60)
        self.assertEqual(sample.records["busco"]["complete"], 95.0)
        self.assertEqual(sample.records["merqury"]["qv"], 30.0)
        self.assertEqual(sample.records["gfastats"]["gfa_num_contigs"], 13)

    def test_missing_component_fails_validation(self):
        next((self.root / "OG1/kmers").glob("*.merqury.qv")).unlink()
        sample = backfill.discover_sample(self.root, "OG1", "250101")
        self.assertFalse(sample.valid)
        self.assertTrue(any("merqury_qv" in error for error in sample.errors))

    def test_separate_flat_fastp_root_is_supported(self):
        fastp_root = self.root / "s3"
        source = next((self.root / "OG1/fastp").glob("*.json"))
        target = fastp_root / "OG1" / source.name
        write(target, source.read_text(encoding="utf-8"))
        source.unlink()
        sample = backfill.discover_sample(self.root, "OG1", "250101", fastp_root)
        self.assertTrue(sample.valid, sample.errors)
        self.assertEqual(sample.paths["fastp"], target)

    def test_validation_cli_writes_all_reports_without_database_connection(self):
        manifest = write(self.root / "manifest.tsv", "og_id\tseq_date\nOG1\t250101\n")
        config = write(self.root / "database.cfg",
                       "dbname=test\nuser=test\npassword=test\nhost=localhost\nport=5432\n")
        reports = self.root / "reports"
        result = backfill.main(["--archive-root", str(self.root), "--manifest", str(manifest),
                                "--db-config", str(config), "--report-dir", str(reports)])
        self.assertEqual(result, 0)
        for name in ("inventory.tsv", "expected.jsonl", "upload.tsv", "verification.tsv", "summary.txt"):
            self.assertTrue((reports / name).is_file(), name)

    def test_duplicate_component_is_rejected(self):
        source = next((self.root / "OG1/fastp").glob("*.json"))
        duplicate = source.with_name("OG1.other.NOVA_250101_AD.fastp.json")
        duplicate.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
        sample = backfill.discover_sample(self.root, "OG1", "250101")
        self.assertFalse(sample.valid)
        self.assertTrue(any("fastp: ambiguous" in error for error in sample.errors))

    def test_malformed_and_mixed_date_inputs_are_rejected(self):
        tiara = next((self.root / "OG1/assemblies/genome/tiara").glob("*.txt"))
        tiara.write_text("wrong\theader\n", encoding="utf-8")
        sample = backfill.discover_sample(self.root, "OG1", "250101")
        self.assertFalse(sample.valid)
        self.assertTrue(any("tiara: malformed" in error for error in sample.errors))
        comp = next((self.root / "OG1/kmers").glob("*completeness.stats"))
        comp.write_text("OG1.ilmn.250102\tset\t1\t2\t50\n", encoding="utf-8")
        sample = backfill.discover_sample(self.root, "OG1", "250101")
        self.assertFalse(sample.valid)


class FakeCursor:
    def __init__(self):
        self.result = None

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, sql, _params=None):
        if sql == "SELECT 1":
            self.result = (1,)

    def fetchone(self):
        return self.result


class FakeConnection:
    def __init__(self):
        self.commits = 0
        self.rollbacks = 0
        self.closed = False

    def cursor(self):
        return FakeCursor()

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1

    def close(self):
        self.closed = True


class TransactionAndVerificationTests(unittest.TestCase):
    def sample(self):
        records = {family: {"og_id": "OG1", "seq_date": "250101",
                            **{column: 1.0 if column in stats.FLOAT_COLUMNS else 1
                               for column in columns}}
                   for family, columns in stats.FAMILY_COLUMNS.items()}
        records["fastp"]["mach"] = "NOVA"
        records["fastp"]["initial"] = "AD"
        records["busco"]["domain"] = "eukaryota"
        records["merqury"]["k_mer_set"] = "reads"
        return backfill.Sample("OG1", "250101", records=records)

    def test_transaction_commits_after_source_comparison(self):
        sample, connection = self.sample(), FakeConnection()
        actual = {column: sample.records[family][column]
                  for family, columns in stats.FAMILY_COLUMNS.items() for column in columns}
        with mock.patch.object(backfill, "runtime_smoke_check"), \
             mock.patch.object(backfill, "connect", return_value=connection), \
             mock.patch.object(backfill, "upsert_record"), \
             mock.patch.object(backfill, "fetch_database_record", return_value=actual):
            uploads, verification = backfill.apply_samples([sample], Path("config"))
        self.assertEqual(uploads[0]["status"], "COMMITTED")
        self.assertEqual(connection.commits, 1)
        self.assertTrue(all(row["status"] == "PASS" for row in verification))

    def test_failure_rolls_back_entire_sample(self):
        sample, connection = self.sample(), FakeConnection()
        def fail_on_busco(_cursor, family, _record):
            if family == "busco":
                raise RuntimeError("forced failure")
        with mock.patch.object(backfill, "runtime_smoke_check"), \
             mock.patch.object(backfill, "connect", return_value=connection), \
             mock.patch.object(backfill, "upsert_record", side_effect=fail_on_busco):
            uploads, _verification = backfill.apply_samples([sample], Path("config"))
        self.assertEqual(uploads[0]["status"], "ROLLED_BACK")
        self.assertEqual(connection.commits, 0)
        self.assertGreaterEqual(connection.rollbacks, 2)

    def test_idempotent_upserts_use_conflict_update(self):
        cursor = mock.Mock()
        record = self.sample().records["assembly"]
        stats.upsert_record(cursor, "assembly", record)
        stats.upsert_record(cursor, "assembly", record)
        self.assertEqual(cursor.execute.call_count, 2)
        self.assertIn("ON CONFLICT (og_id, seq_date) DO UPDATE", cursor.execute.call_args.args[0])

    def test_verification_detects_integer_float_text_and_null_changes(self):
        sample = self.sample()
        expected = {"fastp": sample.records["fastp"], "busco": sample.records["busco"]}
        actual = {**sample.records["fastp"], **sample.records["busco"]}
        actual.update({"total_reads": 2, "q20_rate": 2.0, "domain": "bacteria", "missing": None})
        result = stats.verify_families(expected, actual)
        self.assertTrue(all(row["status"] == "FAIL" for row in result))

    def test_verification_accepts_database_two_decimal_rounding(self):
        expected = {"assembly": {
            "og_id": "OG1", "seq_date": "250101",
            "homozygosity": 98.187, "heterozygosity": 1.98011,
            "genomesize": 1000, "repeatsize": 100, "uniquesize": 900,
            "modelfit": 98.546, "errorrate": 0.523909,
        }}
        actual = {
            "homozygosity": 98.19, "heterozygosity": 1.98,
            "genomesize": 1000, "repeatsize": 100, "uniquesize": 900,
            "modelfit": 98.55, "errorrate": 0.52,
        }
        result = stats.verify_families(expected, actual)
        self.assertEqual(result[0]["status"], "PASS")


if __name__ == "__main__":
    unittest.main()


def genomescope_summary(og_id="OG1", date="250101", haploid="999 bp  1000 bp",
                        model_fit="98%  99%", het="9%  10%"):
    return f"""
name prefix = {og_id}.ilmn.{date}
Homozygous (aa)  90%  91%
Heterozygous (ab)  {het}
Genome Haploid Length  {haploid}
Genome Repeat Length  99 bp  100 bp
Genome Unique Length  899 bp  900 bp
Model Fit  {model_fit}
Read Error Rate  1%  2%
"""


class GenomescopeBothColumnsTests(unittest.TestCase):
    """Until now only the right-hand column of each GenomeScope row reached the database."""

    def parse(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "OG1.ilmn.250101_summary.txt"
            write(path, text)
            return stats.parse_genomescope(path)

    def test_keeps_both_ends_of_every_bound(self):
        record = self.parse(genomescope_summary())
        self.assertEqual(record["genomesize"], 1000)
        self.assertEqual(record["genomesize_min"], 999)
        self.assertEqual(record["repeatsize_min"], 99)
        self.assertEqual(record["uniquesize_min"], 899)
        self.assertAlmostEqual(record["heterozygosity_min"], 9.0)
        self.assertAlmostEqual(record["homozygosity_min"], 90.0)

    def test_model_fit_columns_are_named_for_what_they_are(self):
        # Not min/max: the first is the share of ALL k-mers modelled, the second the fit
        # within the model's own region. OG2644 is stored as an 84% fit while the model
        # accounts for 17% of its k-mers, and that second number has to survive.
        record = self.parse(genomescope_summary(model_fit="17%  84%"))
        self.assertAlmostEqual(record["modelfit"], 84.0)
        self.assertAlmostEqual(record["modelfit_allkmers"], 17.0)

    def test_inverted_model_fit_is_stored_as_written(self):
        # OG2647 reports 82 then 25. They are unrelated statistics, so the parser must not
        # "correct" the order into a min and a max.
        record = self.parse(genomescope_summary(model_fit="81.7476%  25.2118%"))
        self.assertAlmostEqual(record["modelfit_allkmers"], 81.7476)
        self.assertAlmostEqual(record["modelfit"], 25.2118)

    def test_diverged_fit_keeps_the_finite_bound(self):
        # OG3043 reports "2,499,788 bp  Inf bp". This used to raise, so the sample got no
        # GenomeScope row at all; before that it silently stored 2.5 Mb as the estimate.
        record = self.parse(genomescope_summary(haploid="2,499,788 bp  Inf bp"))
        self.assertEqual(record["genomesize_min"], 2499788)
        self.assertIsNone(record["genomesize"])

    def test_no_finite_size_in_either_column_still_raises(self):
        with self.assertRaises(ValueError):
            self.parse(genomescope_summary(haploid="Inf bp  Inf bp"))


class CoverageSummaryTests(unittest.TestCase):
    def test_reads_the_verdict_and_both_coverages(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "OG1.ilmn.250101_coverage_summary.json"
            write(path, json.dumps({
                "genome_size_reliable": False,
                "genome_size_flags": "no_kmer_peak; host_coverage_too_low",
                "fitted_kmer_coverage": 39.1,
                "lambda_depth": 2.5,
            }))
            record = stats.parse_coverage_summary(path)
        self.assertEqual(record["og_id"], "OG1")
        self.assertFalse(record["genome_size_reliable"])
        self.assertIn("no_kmer_peak", record["genome_size_flags"])
        self.assertAlmostEqual(record["kmercov"], 39.1)
        self.assertAlmostEqual(record["lambda_depth"], 2.5)
