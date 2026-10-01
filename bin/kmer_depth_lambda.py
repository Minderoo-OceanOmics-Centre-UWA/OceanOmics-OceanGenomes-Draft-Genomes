#!/usr/bin/env python3
"""Haploid k-mer coverage (GenomeScope's lambda) measured from mapped read depth.

GenomeScope reports a fitted lambda, and every genome size the pipeline publishes is
that lambda divided into the k-mer mass. Nothing checked it. Read depth over single-copy
BUSCO genes estimates the same quantity from the alignments instead, so the two can be
compared -- and where they disagree by a factor of two, so does the genome size.

The check that made this worth building: on NOVA_260909_LA it agreed with GenomeScope
within 8% on all eight samples whose fits were independently trustworthy (ratio
0.92-1.00), then showed that samples with no k-mer peak have their host single-copy genes
at a fifth to an eighth of the assembly's own mean depth. Those fits were not modelling
the host at all.

This is the second half of the job. The first half -- BAM header, BED, index, bedcov --
belongs to samtools, and runs in the samtools container, which has awk and no python at
all. The BED is built there by bin/busco_complete_bed.awk; this reads what comes back.
"""
import argparse
import glob
import json
import os
import sys

# Reads of length L yield (L - k + 1) k-mers, so that fraction of per-base read depth
# survives into k-mer space.
DEFAULT_K = 21


def mean_depth(bedcov_path):
    """Length-weighted mean per-base depth over every region.

    bedcov columns are chrom, start, end, name, summed depth. Summing depth and length
    separately rather than averaging per-region means keeps short genes from carrying the
    same weight as long ones.
    """
    total_depth = total_bases = regions = 0
    with open(bedcov_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                start, end, summed = int(fields[1]), int(fields[2]), float(fields[-1])
            except ValueError:
                continue
            if end <= start:
                continue
            total_depth += summed
            total_bases += end - start
            regions += 1
    if not total_bases:
        return None, 0
    return total_depth / total_bases, regions


def seqkit_total_length(path):
    """Assembly size from a seqkit stats TSV (the sum_len column)."""
    with open(path) as handle:
        rows = [line.split("\t") for line in handle.read().splitlines() if line.strip()]
    if len(rows) < 2:
        return None
    try:
        index = rows[0].index("sum_len")
    except ValueError:
        index = 4  # seqkit's fixed order: file format type num_seqs sum_len ...
    try:
        return int(rows[1][index].replace(",", ""))
    except (IndexError, ValueError):
        return None


def read_contig_depth(path):
    """Per-contig length and mean depth from `samtools coverage`.

    Columns are #rname startpos endpos numreads covbases coverage meandepth meanbaseq
    meanmapq, and meandepth is exact: it is computed from the pileup rather than from
    reads x read length, which over-counts soft-clipped and partially-mapped reads. The
    hand validation that motivated this module used the idxstats approximation; this does
    not.

    Returns a list of (name, length, meandepth).
    """
    rows = []
    try:
        with open(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7:
                    continue
                try:
                    start, end, depth = int(fields[1]), int(fields[2]), float(fields[6])
                except ValueError:
                    continue
                length = end - start + 1
                if length <= 0:
                    continue
                rows.append((fields[0], length, depth))
    except OSError:
        return []
    return rows


def partition_by_depth(rows, host_depth, window=2.0, min_len=500):
    """Split an assembly into host, symbiont and low-depth mass by per-contig depth.

    The host's single-copy genes define one depth; contigs sitting at it are the host's,
    contigs far above it belong to something present in more copies (a symbiont, or a
    high-copy repeat) and contigs far below it are partial or foreign. That is enough to
    say how big the host genome is without a k-mer peak, which is what fails below about
    10x host coverage.

    The window is multiplicative and symmetric: host_depth / window <= meandepth <=
    host_depth * window. Contigs shorter than min_len are excluded entirely -- at a few
    hundred bases a mean depth is noise, and they would otherwise be assigned by it.

    This is a LOWER bound on host genome size. A symbiont that happens to sit at the host's
    depth counts as host, and a high-copy repeat that is genuinely the host's counts as
    symbiont. Describe it as a lower bound everywhere it surfaces.
    """
    totals = {"host": 0, "below": 0, "above": 0, "excluded": 0}
    counts = {"host": 0, "below": 0, "above": 0, "excluded": 0}
    if not rows:
        return totals, counts, None
    if not host_depth or host_depth <= 0 or not window or window <= 0:
        return totals, counts, None

    low, high = host_depth / window, host_depth * window
    for _name, length, depth in rows:
        if length < min_len:
            totals["excluded"] += length
            counts["excluded"] += 1
            continue
        if depth < low:
            bucket = "below"
        elif depth > high:
            bucket = "above"
        else:
            bucket = "host"
        totals[bucket] += length
        counts[bucket] += 1
    return totals, counts, (low, high)


def fitted_kmercov(path):
    """GenomeScope's fitted kmercov, from its model file."""
    try:
        with open(path) as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "kmercov":
                    return float(parts[1])
    except (OSError, ValueError):
        pass
    return None


def total_bases(fastp_json):
    """Bases after filtering. Missing is survivable: it only costs the assembly-wide depth
    comparison, and this module must always emit a row."""
    try:
        with open(fastp_json) as handle:
            return int(json.load(handle)["summary"]["after_filtering"]["total_bases"])
    except (OSError, ValueError, KeyError, TypeError):
        return None


def read_length(fastp_json):
    """Mean read length after filtering. Read from the report, never assumed: the
    k-mer conversion is linear in it."""
    with open(fastp_json) as handle:
        after = json.load(handle)["summary"]["after_filtering"]
    lengths = [after.get("read1_mean_length"), after.get("read2_mean_length")]
    lengths = [float(x) for x in lengths if x]
    return sum(lengths) / len(lengths) if lengths else None


def lambda_from_depth(depth, length, kmer=DEFAULT_K):
    """Read depth -> haploid k-mer coverage.

    Halved because MEGAHIT collapses haplotypes: both homologs of a diploid pile onto one
    locus, so the observed depth corresponds to GenomeScope's homozygous peak at
    2 x lambda rather than to lambda itself.
    """
    if not depth or not length or length <= kmer:
        return None
    return depth * (length - kmer + 1) / length / 2


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bedcov", help="samtools bedcov output; absent means no measurement")
    parser.add_argument("--status", default="ok",
                        help="why there is no measurement, when there is none")
    parser.add_argument("--fastp-json", required=True)
    parser.add_argument("--seqkit-stats", help="stats for the decontaminated assembly")
    parser.add_argument("--model", help="GenomeScope model file, for the fitted kmercov")
    parser.add_argument("--contig-depth",
                        help="samtools coverage output; absent means no depth partition")
    parser.add_argument("--host-depth-window", type=float, default=2.0,
                        help="multiplicative half-width of the host depth band")
    parser.add_argument("--host-depth-min-contig", type=int, default=500,
                        help="contigs shorter than this are too short for a stable depth")
    parser.add_argument("--kmer", type=int, default=DEFAULT_K)
    parser.add_argument("--sample", default="")
    args = parser.parse_args(argv)

    # Every file is parsed here rather than in the shell, because the shell this used to
    # run in had no python and the awk substitutes were where the bugs were.
    assembly_size = seqkit_total_length(args.seqkit_stats) if args.seqkit_stats else None
    kmercov = fitted_kmercov(args.model) if args.model else None
    bases = total_bases(args.fastp_json)

    # No bedcov means no measurement -- too few usable BUSCO regions, or samtools failed.
    # The row is still emitted, with nulls and a status saying why, so the sample stays in
    # every downstream join.
    depth, regions = mean_depth(args.bedcov) if args.bedcov and os.path.exists(args.bedcov) else (None, 0)
    length = read_length(args.fastp_json)
    lam = lambda_from_depth(depth, length, args.kmer)

    # Mean depth across the whole assembly, for comparison. Where the host's conserved
    # single-copy genes sit far below it, the assembly is not one organism and the fit
    # GenomeScope found is not the host's -- that distinction is the point of this module,
    # and it is invisible from lambda alone.
    assembly_depth = (bases / assembly_size if bases and assembly_size else None)

    # Depth partitioning: how much of this assembly sits at the host's own depth.
    # Independent of the k-mer histogram, so it still answers below 10x host coverage
    # where the k-mer route collapses -- on NOVA_260909_LA not one of the 32 samples with
    # measured host lambda under 10 produced an estimate within 2x of its own assembly.
    host_size = host_fraction = above_size = below_size = None
    host_contigs = None
    partition_status = "no_contig_depth"
    if args.contig_depth and os.path.exists(args.contig_depth):
        contig_rows = read_contig_depth(args.contig_depth)
        if not depth:
            # No BUSCO depth means no host depth to key the partition on. The row still
            # carries the sizes as nulls rather than a guess.
            partition_status = "no_host_depth"
        elif not contig_rows:
            partition_status = "no_contig_depth"
        else:
            totals, counts, _band = partition_by_depth(
                contig_rows, depth, args.host_depth_window, args.host_depth_min_contig)
            host_size = totals["host"]
            above_size = totals["above"]
            below_size = totals["below"]
            host_contigs = counts["host"]
            partitioned = totals["host"] + totals["above"] + totals["below"]
            host_fraction = host_size / partitioned if partitioned else None
            partition_status = "ok"

    result = {
        "sample_id": args.sample,
        "status": args.status if lam is None and args.status == "ok" else args.status,
        "busco_regions": regions,
        "busco_mean_depth": depth,
        "read_length": length,
        "lambda_depth": lam,
        "assembly_mean_depth": assembly_depth,
        "busco_depth_over_assembly_depth": (depth / assembly_depth
                                            if depth and assembly_depth else None),
        "assembly_size": assembly_size,
        "genomescope_kmercov": kmercov,
        "kmercov_over_lambda_depth": (kmercov / lam if kmercov and lam else None),
        # A lower bound on the host genome size, and on how much of the assembly is the
        # animal. OG2906 (clean control) gives 344.6 Mb at 91.8% against GenomeScope's
        # 358.9 Mb; OG2634 gives 207.0 Mb at 63.2% with 13.6% above host depth.
        "host_assembly_size": host_size,
        "host_assembly_fraction": host_fraction,
        "symbiont_assembly_size": above_size,
        "low_depth_assembly_size": below_size,
        "host_contig_count": host_contigs,
        "partition_status": partition_status,
    }
    if lam is None and result["status"] == "ok":
        result["status"] = "no_measurement"
    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
