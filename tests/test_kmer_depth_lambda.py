"""Tests for the read-depth lambda measurement.

The calibration that makes this worth trusting is on real data: OG2906's fit is 98.6% with
its k-mer peak at 2x kmercov, so the measurement has to reproduce its fitted 29.7. It
returns 30.2. These tests cover the parts that can fail silently around that.
"""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import kmer_depth_lambda as kdl

AWK_BED = ROOT / "bin" / "busco_complete_bed.awk"


def run_awk_bed(contigs, full_table, count_out=None):
    """Invoke the BED builder the pipeline actually uses."""
    cmd = ["awk"]
    if count_out:
        cmd += ["-v", f"count_out={count_out}"]
    cmd += ["-f", str(AWK_BED), str(contigs), str(full_table)]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in out.stdout.splitlines():
        if not line.strip():
            continue
        contig, start, end, name = line.split("\t")
        rows.append((contig, int(start), int(end), name))
    return rows, out.stderr


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


FULL_TABLE = (
    "# BUSCO version is: 6.0.0\n"
    "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\n"
    "b1\tComplete\tctg1\t100\t200\t+\n"
    "b2\tComplete\tctg1\t800\t500\t-\n"        # minus strand: start > end
    "b3\tComplete\tctg_missing\t10\t20\t+\n"   # contig the BAM does not have
    "b4\tFragmented\tctg1\t10\t20\t+\n"        # not Complete
    "b5\tComplete\tctg2\t900\t1200\t+\n"       # runs past the contig end
)


class RegionTests(unittest.TestCase):
    """The BED builder is bin/busco_complete_bed.awk, driven here as the pipeline drives it.

    It is awk because it runs in the samtools container, which has no python. These are the
    same assertions the python implementation carried before it was removed -- one
    implementation of the coordinate handling, not two that can drift.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.table = write(root / "full_table.tsv", FULL_TABLE)
        self.contigs = write(root / "contigs.tsv", "ctg1\t5000\nctg2\t1000\n")
        self.rows, self.stderr = run_awk_bed(self.contigs, self.table)

    def tearDown(self):
        self.tmp.cleanup()

    def test_only_complete_buscos_are_used(self):
        self.assertNotIn("b4", [r[3] for r in self.rows])

    def test_minus_strand_coordinates_are_normalised(self):
        b2 = next(r for r in self.rows if r[3] == "b2")
        self.assertEqual((b2[1], b2[2]), (499, 800))

    def test_one_based_inclusive_becomes_zero_based_half_open(self):
        b1 = next(r for r in self.rows if r[3] == "b1")
        self.assertEqual((b1[1], b1[2]), (99, 200))

    def test_regions_on_unknown_contigs_are_dropped(self):
        # A stale BUSCO table names contigs from an assembly that no longer exists.
        # samtools bedcov fails the whole sample on the first such line.
        self.assertNotIn("b3", [r[3] for r in self.rows])
        self.assertIn("1 unusable", self.stderr)

    def test_regions_are_clipped_to_the_contig_end(self):
        b5 = next(r for r in self.rows if r[3] == "b5")
        self.assertEqual(b5[2], 1000)

    def test_writes_the_usable_region_count(self):
        count = Path(self.tmp.name) / "n.txt"
        run_awk_bed(self.contigs, self.table, count_out=count)
        self.assertEqual(count.read_text().strip(), "3")


class LambdaTests(unittest.TestCase):
    def test_reproduces_the_calibration_sample(self):
        # OG2906: 70.0x mean BUSCO depth at 146 bp reads and k=21 gives 30.2, against a
        # fitted kmercov of 29.7. If this drifts, nothing downstream is trustworthy.
        self.assertAlmostEqual(kdl.lambda_from_depth(70.0, 146), 30.2, places=1)

    def test_halved_for_haplotype_collapse(self):
        # Both homologs pile onto one collapsed locus, so observed depth is the homozygous
        # peak at 2 x lambda.
        full = 100.0 * (146 - 21 + 1) / 146
        self.assertAlmostEqual(kdl.lambda_from_depth(100.0, 146) * 2, full)

    def test_read_length_shorter_than_k_yields_nothing(self):
        self.assertIsNone(kdl.lambda_from_depth(70.0, 15))

    def test_missing_depth_yields_nothing(self):
        self.assertIsNone(kdl.lambda_from_depth(None, 146))


class MeanDepthTests(unittest.TestCase):
    def test_weights_by_region_length(self):
        # A 10 bp region at depth 100 and a 90 bp region at depth 10 average to 19, not 55.
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "bedcov.tsv",
                         "ctg1\t0\t10\tb1\t1000\nctg1\t100\t190\tb2\t900\n")
            depth, regions = kdl.mean_depth(str(path))
        self.assertEqual(regions, 2)
        self.assertAlmostEqual(depth, 1900 / 100)


if __name__ == "__main__":
    unittest.main()


class UnmeasuredSampleTests(unittest.TestCase):
    """A sample with too few Complete BUSCOs must still produce a row.

    12 of the 80 samples on NOVA_260909_LA have fewer than 20, because their assemblies
    are too fragmented to recover them. Failing the task for those would drop them out of
    the RECHECK join, and so out of the upload join, and so out of the database entirely --
    strictly worse than the unmeasured row they get instead.
    """

    def test_awk_reports_a_low_count_without_failing(self):
        # Too few regions is a finding, not an error: 12 of 80 samples are in that
        # position and failing them would remove them from the database.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = write(root / "ft.tsv", "# h\nb1\tComplete\tctg1\t100\t200\t+\n")
            contigs = write(root / "c.tsv", "ctg1\t5000\n")
            count = root / "n.txt"
            rows, _stderr = run_awk_bed(contigs, table, count_out=count)
            self.assertEqual(len(rows), 1)
            self.assertEqual(count.read_text().strip(), "1")

    def test_summarise_without_bedcov_emits_a_row_with_a_reason(self):
        import contextlib
        import io
        import json
        with tempfile.TemporaryDirectory() as tmp:
            fastp = write(Path(tmp) / "f.json", json.dumps(
                {"summary": {"after_filtering": {"read1_mean_length": 146,
                                                 "read2_mean_length": 146,
                                                 "total_bases": 1000}}}))
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                kdl.main(["--fastp-json", str(fastp),
                          "--status", "too_few_busco_regions", "--sample", "OG2639"])
        result = json.loads(buf.getvalue())
        self.assertEqual(result["sample_id"], "OG2639")
        self.assertEqual(result["status"], "too_few_busco_regions")
        self.assertIsNone(result["lambda_depth"])


def coverage_tsv(rows):
    """A `samtools coverage` table from (name, length, meandepth) triples.

    Columns are #rname startpos endpos numreads covbases coverage meandepth meanbaseq
    meanmapq; only the first, third and seventh are read, but the shape has to be right
    or a format change goes unnoticed.
    """
    out = ["#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"]
    for name, length, depth in rows:
        out.append(f"{name}\t1\t{length}\t0\t{length}\t100.0\t{depth}\t36\t60")
    return "\n".join(out) + "\n"


def run_main(tmp, contig_rows=None, bedcov=None, extra=None):
    """kmer_depth_lambda.main() over a minimal set of inputs, as the module drives it."""
    import contextlib
    import io
    import json as _json
    root = Path(tmp)
    fastp = write(root / "f.json", _json.dumps(
        {"summary": {"after_filtering": {"read1_mean_length": 146,
                                         "read2_mean_length": 146,
                                         "total_bases": 1_000_000}}}))
    argv = ["--fastp-json", str(fastp), "--sample", "OGTEST"]
    if bedcov is not None:
        argv += ["--bedcov", str(write(root / "bedcov.tsv", bedcov))]
    if contig_rows is not None:
        argv += ["--contig-depth", str(write(root / "cd.tsv", coverage_tsv(contig_rows)))]
    argv += extra or []
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        kdl.main(argv)
    return _json.loads(buf.getvalue())


class PartitionTests(unittest.TestCase):
    """The depth partition, which is the only route to a host genome size below about 10x.

    The k-mer route collapses there: of the 32 NOVA_260909_LA samples with measured host
    lambda under 10, not one produced an estimate within 2x of its own assembly.
    """

    def test_a_clean_assembly_is_almost_all_host(self):
        # OG2906's shape: everything at the host's own depth.
        rows = [(f"c{i}", 1_000_000, 94.0) for i in range(10)]
        totals, counts, band = kdl.partition_by_depth(rows, 92.0)
        self.assertEqual(totals["host"], 10_000_000)
        self.assertEqual((totals["above"], totals["below"]), (0, 0))
        self.assertEqual(counts["host"], 10)
        self.assertAlmostEqual(band[0], 46.0)

    def test_symbiont_mass_is_separated_from_host_mass(self):
        # OG2634's shape: BUSCO depth 8.9x, with about half the assembly's length sitting
        # at roughly 6x that -- 46% of it by depth is not the animal.
        rows = ([(f"h{i}", 1_000_000, 9.0) for i in range(17)]
                + [(f"s{i}", 1_000_000, 55.0) for i in range(15)])
        totals, _counts, _band = kdl.partition_by_depth(rows, 8.9)
        self.assertEqual(totals["host"], 17_000_000)
        self.assertEqual(totals["above"], 15_000_000)
        self.assertAlmostEqual(totals["host"] / (totals["host"] + totals["above"]),
                               0.53, places=2)

    def test_short_contigs_are_excluded_whatever_their_depth(self):
        # At a few hundred bases a mean depth is noise. Without the length floor this
        # 3000x fragment would be counted as symbiont mass.
        rows = [("long", 1_000_000, 50.0), ("short", 200, 3000.0)]
        totals, counts, _band = kdl.partition_by_depth(rows, 50.0)
        self.assertEqual(totals["above"], 0)
        self.assertEqual(counts["excluded"], 1)
        self.assertEqual(totals["excluded"], 200)

    def test_json_carries_the_host_size_and_fraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            # 146 bp reads at k=21 halve into lambda, so a BUSCO depth of 100 gives a host
            # depth of 100 -- the partition keys on the depth, not on lambda.
            result = run_main(
                tmp,
                bedcov="ctg1\t0\t1000\tb1\t100000\n",
                contig_rows=[("h1", 2_000_000, 100.0), ("h2", 1_000_000, 90.0),
                             ("s1", 1_000_000, 600.0)])
        self.assertEqual(result["partition_status"], "ok")
        self.assertEqual(result["host_assembly_size"], 3_000_000)
        self.assertEqual(result["symbiont_assembly_size"], 1_000_000)
        self.assertEqual(result["host_contig_count"], 2)
        self.assertAlmostEqual(result["host_assembly_fraction"], 0.75)

    def test_no_host_depth_reports_itself_rather_than_guessing(self):
        # No bedcov means no BUSCO depth, so there is nothing to key the partition on.
        # The row still has to be emitted: a sample that drops out of the join drops out
        # of the database.
        with tempfile.TemporaryDirectory() as tmp:
            result = run_main(tmp, contig_rows=[("c1", 1_000_000, 50.0)])
        self.assertEqual(result["partition_status"], "no_host_depth")
        self.assertIsNone(result["host_assembly_size"])
        self.assertIsNone(result["host_assembly_fraction"])

    def test_no_contig_depth_file_is_not_a_partition(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_main(tmp, bedcov="ctg1\t0\t1000\tb1\t100000\n")
        self.assertEqual(result["partition_status"], "no_contig_depth")
        self.assertIsNone(result["host_assembly_size"])
