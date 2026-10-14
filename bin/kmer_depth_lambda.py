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
import math
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


def partition_by_depth(rows, host_depth, window=2.0, min_len=500, collapsed_floor=0.75):
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

    Within the host band, contigs below collapsed_floor x host_depth are reported apart as
    host_haplotig. Where heterozygosity is high MEGAHIT keeps both haplotypes of a locus as
    separate contigs, each at about half the collapsed depth, so the plain host total
    counts that locus twice. Across the NOVA_260909_LA fishes the half-depth mass tracks
    heterozygosity (OG3078 0.70% het, 29 Mb; OG3074 2.57% het, 169 Mb), and it is what
    put OG3074's host at 494 Mb against GenomeScope's 338.
    """
    totals = {"host": 0, "host_collapsed": 0, "host_haplotig": 0,
              "below": 0, "above": 0, "excluded": 0}
    counts = dict.fromkeys(totals, 0)
    if not rows:
        return totals, counts, None
    if not host_depth or host_depth <= 0 or not window or window <= 0:
        return totals, counts, None

    low, high = host_depth / window, host_depth * window
    split = max(low, host_depth * collapsed_floor)
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
            sub = "host_haplotig" if depth < split else "host_collapsed"
            totals[sub] += length
            counts[sub] += 1
        totals[bucket] += length
        counts[bucket] += 1
    return totals, counts, (low, high)


def find_depth_mode(rows, min_len=500, bin_width=0.25):
    """The length-weighted modal per-contig depth, or None when there are no usable rows.

    Binned on log2 depth, because the masses being told apart -- host, haplotig, symbiont
    -- differ by factors, not by fixed amounts. Returns the geometric centre of the
    heaviest bin.
    """
    bins = {}
    for _name, length, depth in rows:
        if length < min_len or depth <= 0:
            continue
        key = math.floor(math.log2(depth) / bin_width)
        bins[key] = bins.get(key, 0) + length
    if not bins:
        return None
    best = max(bins, key=lambda k: bins[k])
    return 2 ** ((best + 0.5) * bin_width)


def partition_confidence(regions, host_depth, mode, mode_fraction, window=2.0,
                         min_regions=50, min_mode_fraction=0.4):
    """Whether the BUSCO-anchored partition can be trusted, and if not, why.

    The anchor is depth over the Complete BUSCOs that were found. With few of them, they
    are found where the assembly is best covered, so the anchor sits above the host's
    typical depth: OG2942's 60 regions put it at 27x while 80% of the assembly formed
    one mode near 6x. Returns (confidence, note).

    The anchor-to-mode test does the work. On NOVA_260909_LA it caught every partition
    known to be wrong (OG2942 4.4x, OG2629 8.0x, OG2645 16.8x off the mode), and on
    OG2948 the mode-anchored size (414 Mb) matched GenomeScope's 410 where the BUSCO
    anchor gave 341. The region floor is only a noise guard: at 200 it marked sound
    partitions low (OG2966, OG3024, OG2989, anchors within 20% of the mode) without
    catching anything the mode test missed.
    """
    if not regions or regions < min_regions:
        return "low", "few_busco_regions"
    if mode is None or mode_fraction is None or mode_fraction < min_mode_fraction:
        return "low", "no_dominant_mode"
    if host_depth / mode > window or mode / host_depth > window:
        return "low", "anchor_off_mode"
    return "high", None


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
    parser.add_argument("--host-depth-min-regions", type=int, default=50,
                        help="fewer Complete BUSCO regions than this and the anchor is low confidence")
    parser.add_argument("--host-mode-min-fraction", type=float, default=0.4,
                        help="share of partitioned mass the depth mode's band must hold to count as dominant")
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
    collapsed_size = haplotig_size = haploid_size = None
    mode = anchor_over_mode = mode_anchored_size = None
    confidence = note = None
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
            collapsed_size = totals["host_collapsed"]
            haplotig_size = totals["host_haplotig"]
            haploid_size = collapsed_size + haplotig_size // 2
            partition_status = "ok"

            # Where does the assembly's own mass sit, independently of the BUSCO anchor?
            mode = find_depth_mode(contig_rows, args.host_depth_min_contig)
            mode_fraction = None
            if mode:
                anchor_over_mode = depth / mode
                mode_totals, _c, _b = partition_by_depth(
                    contig_rows, mode, args.host_depth_window, args.host_depth_min_contig)
                mode_fraction = mode_totals["host"] / partitioned if partitioned else None
            confidence, note = partition_confidence(
                regions, depth, mode, mode_fraction, args.host_depth_window,
                args.host_depth_min_regions, args.host_mode_min_fraction)
            # Reported beside the BUSCO-anchored size, never in place of it: a dominant
            # mode is usually the host, but in a symbiont-dominated sponge it need not be.
            if (mode_fraction is not None and mode_fraction >= args.host_mode_min_fraction
                    and (anchor_over_mode > args.host_depth_window
                         or anchor_over_mode < 1 / args.host_depth_window)):
                mode_anchored_size = mode_totals["host"]

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
        # Haplotype-split pairs counted once. host_assembly_size counts both copies.
        "host_collapsed_size": collapsed_size,
        "host_haplotig_size": haplotig_size,
        "host_haploid_size": haploid_size,
        "depth_mode": mode,
        "anchor_over_mode": anchor_over_mode,
        "partition_confidence": confidence,
        "partition_note": note,
        "host_assembly_size_mode_anchored": mode_anchored_size,
    }
    if lam is None and result["status"] == "ok":
        result["status"] = "no_measurement"
    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
