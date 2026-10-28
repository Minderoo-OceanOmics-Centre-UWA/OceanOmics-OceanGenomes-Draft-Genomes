#!/usr/bin/env python3
"""Decide whether a GenomeScope2 fit can be trusted as a genome size.

Everything downstream of the size estimate divides by it -- theoretical coverage, the
coverage grade, gfastats NG50 -- so a bad fit silently becomes a confident-looking
verdict. This module turns "did GenomeScope actually fit this sample" into an explicit
answer with a named reason.

Design notes, each of which is a bug this replaces:

* **GenomeScope's two "Model Fit" values are not min and max.** They are printed under a
  "min   max" header but its own R source emits two unrelated statistics:

      sprintf(format_column_2, percentage_format(model_fit_allscore[1])),
      sprintf(format_column_3, percentage_format(model_fit_fullscore[1]))

  The first is "Percent Kmers Modeled (All Kmers)" -- the share of every k-mer left in the
  histogram that the model accounts for, so it falls mechanically as --max_kmercov admits
  more unmodelled tail. The second is "Percent Kmers Modeled (Full Model)", the fit within
  the model's own region, which does not move with -m. Here they are model_fit_allkmers
  and model_fit_full.

  That is why they can invert (OG2647 82/25, OG2624 70/15) -- unrelated statistics, not
  bounds -- and why gating on the first meant the gate measured its own -m setting: 57 of
  64 re-run samples were flagged, including OG2906 whose model_fit_full was 98.6498% at
  every cutoff tested. Gate on model_fit_full; report model_fit_allkmers, because a sample
  can be a 98% fit to 17% of its k-mers and that is worth knowing.

* **Non-finite bounds are a failure, not a number.** When the model diverges GenomeScope
  writes "Inf bp" (OG3043 does). A regex of `([0-9,]+) bp` cannot match it, so taking the
  last match silently returned the LOWER bound and published 2.5 Mb as that sample's
  genome size.

* **A k-mer peak has to actually exist.** The previous detector walked the histogram
  looking for a trough, re-seating it every time the count fell; on a monotonically
  decreasing histogram -- which is what a coverage-starved library looks like -- it never
  terminated properly and invented a "peak" out in the sparse tail, then flagged the
  disagreement it had just manufactured. That fired on 22 of 23 tier-4 samples. The
  honest answer for those libraries is no_kmer_peak, which means resequence.

* **Degeneracy is a property of the model parameters.** The old check compared the two
  Model Fit values for equality, which stopped firing for OG2983 the moment -m changed,
  letting a 10.9 Mb "genome" at kmercov 569 and bias 37 be graded EXCELLENT.
"""
import argparse
import json
import re
import sys

# A GenomeScope bound is an integer with thousands separators, or a non-finite marker.
NON_FINITE = re.compile(r"-?(inf|na|nan)$", re.I)

DEFAULTS = {
    "min_model_fit": 50.0,          # against model_fit_full; recalibrate from the sweep
    "max_bias": 5.0,                # a diploid fit does not need a 37x bias term
    "min_peak_kmercov_ratio": 0.5,  # peak below half the fitted coverage: wrong mode
    "max_peak_kmercov_ratio": 3.0,  # peak at 1x or 2x kmercov is expected; beyond is not
    "peak_search_max": 1000,        # deliberately NOT params.genomescope2_m
    "min_assembly_ratio": 0.3,      # size estimate vs decontaminated assembly size
    "max_assembly_ratio": 2.0,
    "max_lambda_ratio": 1.5,        # fitted kmercov vs lambda measured from read depth
    "min_host_depth_ratio": 0.4,    # BUSCO-gene depth over assembly-wide depth
    "min_seedable_lambda": 5.0,     # below this there is no host signal to fit
    "min_reseed_fit_gain": 2.0,     # points of model_fit_full a reseed must win to be kept
    "min_busco_complete": 10.0,     # BUSCO C% below which the assembly is not the host's
}


# --------------------------------------------------------------------- parsing
def _bound(token):
    """One bound from a GenomeScope summary line. None when non-finite."""
    token = token.strip()
    if NON_FINITE.match(token):
        return None
    try:
        return int(token.replace(",", ""))
    except ValueError:
        return None


def parse_summary(path):
    """Both reported values for each property. A non-finite one comes back as None and the
    caller must treat it as a failed fit, never fall back to the other.

    For the length and heterozygosity rows the pair really is a lower and upper bound. For
    Model Fit it is not -- see the module docstring -- so those two come back named
    model_fit_allkmers and model_fit_full rather than _min and _max."""
    with open(path) as handle:
        text = handle.read()

    out = {}
    for key, label in (("haploid", "Genome Haploid Length"),
                       ("repeat", "Genome Repeat Length"),
                       ("unique", "Genome Unique Length")):
        match = re.search(rf"^{label}\s+(\S+) bp\s+(\S+) bp", text, flags=re.M)
        lo, hi = (_bound(match.group(1)), _bound(match.group(2))) if match else (None, None)
        out[f"{key}_min"], out[f"{key}_max"] = lo, hi

    for label, first, second in (
            ("Model Fit", "model_fit_allkmers", "model_fit_full"),
            (r"Heterozyg(?:osity|ous \(ab\))", "heterozygosity_min", "heterozygosity_max"),
            ("Homozygous \\(a+\\)", "homozygosity_min", "homozygosity_max"),
            ("Read Error Rate", "error_rate_min", "error_rate_max")):
        match = re.search(rf"^{label}\s+(\S+)%\s+(\S+)%", text, flags=re.M)
        if match:
            lo, hi = match.group(1), match.group(2)
            out[first] = None if NON_FINITE.match(lo) else float(lo)
            out[second] = None if NON_FINITE.match(hi) else float(hi)
        else:
            out[first] = out[second] = None
    return out


def parse_model(path):
    """kmercov, bias and the other nls coefficients. Missing file is not fatal: the
    summary alone still supports most of the gate."""
    out = {"kmercov": None, "bias": None, "d": None, "r1": None}
    try:
        with open(path) as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 2 and parts[0] in out:
                    try:
                        out[parts[0]] = float(parts[1])
                    except ValueError:
                        pass
    except OSError:
        pass
    return out


def read_histogram(path, max_cov):
    hist = {}
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                cov, count = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            if 1 <= cov <= max_cov:
                hist[cov] = count
    return hist


# ----------------------------------------------------------------- peak finding
def find_kmer_peak(hist, lookahead=5, turn=1.05, prominence=1.3, tail_fraction=1e-3):
    """The dominant genomic mode in a k-mer histogram, or None when there isn't one.

    An Illumina histogram with a real genome in it falls steeply from cov=1 (errors),
    troughs, then rises into the genomic peak. Missing that rise is diagnostic, not an
    edge case: it is what a library sequenced below the level where the genomic peak
    separates from the error curve looks like, and the only useful response is more data.

    Returns (peak_coverage, reason, trough_coverage). peak_coverage is None unless
    reason == "ok".
    """
    covs = sorted(hist)
    if len(covs) < 20:
        return None, "histogram_too_short", None

    # Confine the search to where real k-mer mass still is. Past that the counts are a
    # sparse decaying tail, and its noise will happily present as a trough followed by a
    # "peak" -- which is exactly how the old detector ended up reporting a peak at 1958x
    # for a sample whose fitted coverage was 33.7.
    floor = hist[covs[0]] * tail_fraction
    covs = [cov for cov in covs if hist[cov] >= floor]
    if len(covs) < 20:
        return None, "histogram_too_short", None

    trough = None
    for i in range(len(covs) - 1):
        window = range(i + 1, min(len(covs), i + 1 + lookahead))
        if any(hist[covs[j]] > hist[covs[i]] * turn for j in window):
            trough = i
            break
    if trough is None:
        return None, "no_trough_histogram_monotone", None

    best = max(range(trough, len(covs)), key=lambda j: hist[covs[j]])
    trough_count, peak_count = hist[covs[trough]], hist[covs[best]]
    if trough_count <= 0 or peak_count < trough_count * prominence:
        return None, "no_mode_above_trough", covs[trough]
    return covs[best], "ok", covs[trough]


# ------------------------------------------------------------------- the gate
def evaluate(summary, model, hist=None, assembly_size=None, thresholds=None):
    """Flags explaining why this genome size estimate cannot be trusted. Empty means it
    can. Returns (flags, details)."""
    cfg = dict(DEFAULTS, **(thresholds or {}))
    flags = []

    size = summary.get("haploid_max")
    fit_all = summary.get("model_fit_allkmers")
    fit_full = summary.get("model_fit_full")
    het_min, het_max = summary.get("heterozygosity_min"), summary.get("heterozygosity_max")

    # A diverged model reports an infinite upper bound. Never silently substitute the
    # lower one -- for OG3043 that meant publishing 2.5 Mb for a 322 Mb assembly.
    if size is None:
        flags.append("non_finite_genome_size" if summary.get("haploid_min") is not None
                     else "no_genome_size_estimate")
    elif size <= 0:
        flags.append("no_genome_size_estimate")

    if fit_all is not None and fit_full is not None and fit_all > fit_full:
        # The model accounts for a larger share of ALL k-mers than of its own modelled
        # region, which a coherent fit cannot do (OG2647 82/25, OG2624 70/15).
        flags.append("inverted_model_fit_bounds")

    # Gate on the fit within the model's own region. The all-k-mers figure moves with
    # --max_kmercov, so gating on it would gate on the parameter.
    if fit_full is not None and fit_full < cfg["min_model_fit"]:
        flags.append(f"model_fit_below_{cfg['min_model_fit']:g}pc")

    if het_max is None or het_max <= 0 or het_max > 100:
        flags.append("implausible_heterozygosity")
    elif het_min == 0 and het_max > 0:
        # The lower solution collapsed onto zero heterozygosity while the upper did not:
        # the two bounds are not bracketing the same model.
        flags.append("degenerate_model_fit")

    bias, kmercov = model.get("bias"), model.get("kmercov")
    if bias is not None and bias > cfg["max_bias"]:
        # A large bias term is the model buying agreement with a shape that is not a
        # diploid k-mer spectrum (OG2983: bias 37 at kmercov 569 for a 10.9 Mb "genome").
        flags.append("implausible_model_bias")

    peak = None
    if hist is not None:
        peak, reason, _trough = find_kmer_peak(hist)
        if peak is None:
            flags.append("no_kmer_peak")
        elif kmercov and kmercov > 0:
            ratio = peak / kmercov
            # The dominant mode sits at 1x kmercov (heterozygous) or 2x (homozygous), so
            # a ratio anywhere in between is expected. Outside that band GenomeScope
            # fitted something that is not the mode the histogram actually has.
            if ratio > cfg["max_peak_kmercov_ratio"] or ratio < cfg["min_peak_kmercov_ratio"]:
                flags.append("kmer_peak_disagrees_with_fitted_coverage")

    # The check no fit statistic can make. A model can converge beautifully onto the
    # wrong thing: OG3037 fitted at 86% and estimated 209 Mb against a 1.26 Gb assembly.
    if size and assembly_size and assembly_size > 0:
        ratio = size / assembly_size
        if ratio < cfg["min_assembly_ratio"] or ratio > cfg["max_assembly_ratio"]:
            flags.append("genome_size_disagrees_with_assembly")

    details = {
        "estimated_genome_size": size,
        "genome_size_min": summary.get("haploid_min"),
        "genome_size_max": summary.get("haploid_max"),
        "model_fit_allkmers": fit_all,
        "model_fit_full": fit_full,
        # Deprecated aliases, kept for one release so existing readers of the published
        # coverage JSON keep working. model_fit_min was never a minimum.
        "model_fit_min": fit_all,
        "model_fit_max": fit_full,
        "heterozygosity_min": het_min,
        "heterozygosity_max": het_max,
        "kmercov": kmercov,
        "bias": bias,
        "kmer_peak_coverage": peak,
        "assembly_size": assembly_size,
        "genome_size_flags": flags,
        "genome_size_reliable": not flags,
    }
    return flags, details


def assembly_size_from_seqkit(path):
    """Total bases from a seqkit stats TSV (the `sum_len` column)."""
    with open(path) as handle:
        rows = [line.split("\t") for line in handle.read().splitlines() if line.strip()]
    if len(rows) < 2:
        return None
    header, values = rows[0], rows[1]
    try:
        index = header.index("sum_len")
    except ValueError:
        index = 4  # seqkit's fixed column order: file format type num_seqs sum_len ...
    try:
        return int(values[index].replace(",", ""))
    except (IndexError, ValueError):
        return None


def busco_complete_from_summary(path):
    """BUSCO completeness (C%) from a short_summary JSON, or None when it cannot be read.

    Unreadable is not a failure here: a missing completeness costs one check, and the
    rest of the recheck still stands."""
    try:
        with open(path) as handle:
            results = json.load(handle).get("results", {})
    except (OSError, ValueError, AttributeError):
        return None
    value = results.get("Complete percentage")
    if value is None:
        match = re.search(r"C:([0-9.]+)%", results.get("one_line_summary") or "")
        value = match.group(1) if match else None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# Same grades as CALCULATE_SEQUENCING_COVERAGE, for a coverage recomputed from a refit.
COVERAGE_GRADES = (
    (50, "EXCELLENT", "Suitable for high-quality genome assembly and variant calling."),
    (30, "GOOD", "Suitable for genome assembly and most downstream analyses."),
    (20, "ADEQUATE", "Suitable for basic genome assembly."),
    (10, "LOW", "Draft assembly may be possible, but additional sequencing may help."),
    (0, "INSUFFICIENT", "Additional sequencing is strongly recommended."),
)


def apply_reseeded_fit(coverage, summary, model, hist=None, thresholds=None):
    """Swap the provisional summary over to the GenomeScope refit that won.

    GENOMESCOPE_RESEED keeps a refit only when it beats the original on its own terms, but
    the provisional summary was written from the original before the refit existed. Left
    alone, the uploaded verdict described the rejected fit: OG3057 stayed flagged
    fitted_coverage_disagrees_with_read_depth for a kmercov of 10.2 that the refit had
    replaced with 16.7, against a measured 16.6.

    The refit goes through the same k-mer gate as any fit, so winning the reseed does not
    by itself make a sample reliable. Its flags replace the provisional ones, which are all
    k-mer-stage flags; recheck_against_assembly() then adds the assembly checks as usual."""
    flags, details = evaluate(summary, model, hist, None, thresholds)
    merged = dict(coverage)
    size = details["estimated_genome_size"] or 0
    size = size if size > 0 else 0
    for key in ("genome_size_min", "genome_size_max", "model_fit_allkmers", "model_fit_full",
                "heterozygosity_min", "heterozygosity_max", "kmercov", "kmer_peak_coverage"):
        merged[key] = details[key]
    merged["estimated_genome_size"] = size
    # Deprecated aliases, kept in step with the fields they alias.
    merged["model_fit_min"], merged["model_fit_max"] = details["model_fit_allkmers"], details["model_fit_full"]
    merged["fitted_kmer_coverage"] = details["kmercov"]

    bases_after = coverage.get("total_bases_after_filtering") or 0
    retention = coverage.get("data_retention_percent") or 0
    merged["theoretical_coverage"] = bases_after / size if size else 0
    merged["coverage_before_filtering"] = (
        bases_after / (retention / 100) / size if size and retention else 0)
    for floor, status, recommendation in COVERAGE_GRADES:
        if merged["theoretical_coverage"] >= floor:
            merged["coverage_status"], merged["coverage_recommendation"] = status, recommendation
            break

    merged["genome_size_flags"] = "; ".join(flags)
    merged["genome_size_reliable"] = not flags
    merged["genomescope_fit"] = "reseeded"
    return merged


def recheck_against_assembly(coverage, assembly_size, thresholds=None, kmer_depth=None,
                             busco_complete=None):
    """Add the checks that need the finished assembly to an existing coverage summary.

    Kept separate from evaluate() because it runs at a different point in the pipeline:
    the k-mer checks happen before MEGAHIT, these need the assembly and the alignments
    against it. Returns an updated copy of the coverage summary.

    `kmer_depth` is the MEASURE_KMER_COVERAGE result. It carries the one measurement that
    can contradict GenomeScope on its own terms: lambda estimated from read depth over
    single-copy BUSCO genes. A factor of two between that and the fitted kmercov is a
    factor of two in the genome size, in the opposite direction.

    `busco_complete` is the assembly's BUSCO C%. Below the floor the assembly does not
    contain the host's gene space, so a size estimate that agrees with it is agreeing
    with something other than the host: OG3080 was graded EXCELLENT on an 8.9 Mb
    assembly at C:0.0% because its 15.5 Mb estimate matched that assembly, and with too
    few Complete BUSCOs to measure read depth none of the depth checks could fire.
    """
    cfg = dict(DEFAULTS, **(thresholds or {}))
    merged = dict(coverage)

    flags = [f for f in (coverage.get("genome_size_flags") or "").split("; ") if f]
    size = coverage.get("estimated_genome_size") or 0

    # The summary passed in may already carry these flags: on a resumed run the provisional
    # input can itself be a previous recheck's output, so every flag this function raises
    # has to be raised at most once. Appending blind put each k-mer-depth flag in the
    # database twice for all 23 tier-4 samples of NOVA_260909_LA.
    def add(flag):
        if flag not in flags:
            flags.append(flag)

    # GenomeScope writes "-1 bp" when it fails outright, and the parser before this one
    # matched `([0-9,]+) bp` and dropped the sign: OG2653 reached this point as a 1 bp
    # genome at 18,007,370,181x. A summary from a run that old can still arrive here as a
    # precomputed input, so no size this small is taken as an estimate, and the coverage
    # figures derived from it go back to the 0 the coverage step writes for no estimate.
    if size <= 1:
        if "non_finite_genome_size" not in flags:
            add("no_genome_size_estimate")
        size = 0
        merged["estimated_genome_size"] = 0
        merged["theoretical_coverage"] = 0
        merged["coverage_before_filtering"] = 0

    # Which assembly size the ratio is taken against. GenomeScope estimates the HOST
    # genome, so comparing it to a total assembly that can be half symbiont puts the wrong
    # number in the denominator -- in exactly the samples this check exists for. Where the
    # depth partition succeeded, use the host mass instead, and record which was used so
    # the ratio can always be reproduced.
    #
    # Only a confident partition earns that. The partition is keyed on BUSCO-gene depth,
    # and with few Complete BUSCOs those genes sit in the best-covered corner of the
    # assembly rather than at the host's typical depth: OG2942's 60 regions put the anchor
    # at 27x while 80% of its assembly formed one mode near 6x, and the "host" came out
    # as 190 of 1,458 Mb. A partition written before confidence existed is trusted as
    # before, so republishing old rows changes nothing.
    #
    # The haploid host size counts each haplotype-split pair once. The plain host size
    # counts both copies, which is what put OG3074 at 494 Mb against GenomeScope's 338.
    kd = kmer_depth or {}
    host_size = kd.get("host_haploid_size") or kd.get("host_assembly_size")
    confident = kd.get("partition_confidence", "high") == "high"
    if kd.get("partition_status") == "ok" and host_size and confident:
        ratio_basis, ratio_denominator = "host", host_size
    else:
        ratio_basis, ratio_denominator = "total", assembly_size

    if size and ratio_denominator:
        ratio = size / ratio_denominator
        if ratio < cfg["min_assembly_ratio"] or ratio > cfg["max_assembly_ratio"]:
            add("genome_size_disagrees_with_assembly")
    merged["assembly_ratio_basis"] = ratio_basis
    merged["assembly_ratio_denominator"] = ratio_denominator

    if kmer_depth:
        lam = kmer_depth.get("lambda_depth")
        # The provisional summary can carry no kmercov at all: on --skip_genome_assembly
        # tiers CALCULATE_SEQUENCING_COVERAGE never runs, so the summary comes from a
        # published file that predates the field, and 16 rows of NOVA_260909_LA reached the
        # database with kmercov null beside a perfectly good genomescope_kmercov. The
        # measurement parses the same model file every time, so take it when the summary
        # has nothing.
        kmercov = (coverage.get("kmercov") or coverage.get("fitted_kmer_coverage")
                   or kmer_depth.get("genomescope_kmercov"))
        host_ratio = kmer_depth.get("busco_depth_over_assembly_depth")

        merged["lambda_depth"] = lam
        merged["kmercov"] = kmercov
        merged["busco_mean_depth"] = kmer_depth.get("busco_mean_depth")
        merged["busco_regions"] = kmer_depth.get("busco_regions")
        merged["assembly_mean_depth"] = kmer_depth.get("assembly_mean_depth")
        merged["busco_depth_over_assembly_depth"] = host_ratio

        # A lower bound on the host genome size, and on how much of this assembly is the
        # animal. It needs no k-mer peak, so it still answers where GenomeScope cannot.
        for key in ("host_assembly_size", "host_assembly_fraction",
                    "symbiont_assembly_size", "low_depth_assembly_size",
                    "host_contig_count", "partition_status",
                    "host_collapsed_size", "host_haplotig_size", "host_haploid_size",
                    "depth_mode", "anchor_over_mode", "partition_confidence",
                    "partition_note", "host_assembly_size_mode_anchored"):
            merged[key] = kmer_depth.get(key)

        if lam and kmercov:
            ratio = max(kmercov / lam, lam / kmercov)
            merged["kmercov_over_lambda_depth"] = kmercov / lam
            if ratio > cfg["max_lambda_ratio"]:
                add("fitted_coverage_disagrees_with_read_depth")

        # Not a fit problem at all. The host's conserved single-copy genes sitting far
        # below the assembly's own mean depth means the assembly is not one organism, so
        # whatever GenomeScope converged on is not the host's k-mer spectrum.
        if host_ratio is not None and host_ratio < cfg["min_host_depth_ratio"]:
            add("assembly_not_dominated_by_host")

        if lam is not None and lam < cfg["min_seedable_lambda"]:
            # Reported rather than inferred from a bad fit: at this coverage there is no
            # host peak to model and no starting value can produce one.
            add("host_coverage_too_low")

    merged["busco_complete"] = busco_complete
    if busco_complete is not None and busco_complete < cfg["min_busco_complete"]:
        add("assembly_lacks_host_gene_space")

    merged["assembly_size"] = assembly_size
    merged["assembly_cross_check"] = "done" if assembly_size else "no_assembly_size"
    merged["genome_size_flags"] = "; ".join(flags)
    merged["genome_size_reliable"] = not flags

    if flags:
        merged["coverage_status"] = "UNRELIABLE_GENOME_SIZE_ESTIMATE"
        merged["coverage_recommendation"] = (
            "GenomeScope could not fit this sample (" + ", ".join(flags) + "), so the "
            "genome size estimate and every coverage figure derived from it are unreliable. "
            "Assess coverage against the decontaminated assembly size instead."
        )
    return merged


def render_text_report(summary):
    """The human-readable *_sequencing_coverage.txt for a coverage summary.

    One renderer for both the provisional summary and the rechecked one, because the text
    file used to be written only by the provisional step: RECHECK_GENOME_SIZE replaced the
    JSON beside it and left the text saying "reliable: yes / EXCELLENT" for OG2617 after
    the read-depth check had flagged it. Reads nothing but the summary, so the two files
    cannot disagree."""
    s = summary

    def num(value, fmt="{:,}"):
        return "NA" if value is None else fmt.format(value)

    def pct(value):
        return "NA" if value is None else f"{value:.4f}%"

    final = s.get("assembly_cross_check") not in (None, "pending")
    flags = [f for f in (s.get("genome_size_flags") or "").split("; ") if f]
    size = s.get("estimated_genome_size") or 0
    coverage = s.get("theoretical_coverage") or 0

    lines = ["SEQUENCING COVERAGE ANALYSIS", "=" * 50, ""]
    lines.append(f"Sample ID: {s.get('sample_id', 'NA')}")
    lines.append("Verdict: " + ("final (assembly cross-check done)" if final else
                                "PROVISIONAL (assembly cross-check not run yet)"))
    lines.append("")

    lines.append("GENOME SIZE ESTIMATION (GenomeScope):")
    if s.get("genomescope_fit") == "reseeded":
        lines.append("  Fit: reseeded (refit seeded with the measured k-mer coverage, kept over "
                     "the original)")
    elif s.get("genomescope_fit") == "original":
        lines.append("  Fit: original")
    lines.append(f"  Estimated genome size: {num(size) + ' bp' if size else 'none'}")
    if s.get("genome_size_min") and s.get("genome_size_max"):
        lines.append(f"  Genome size range: {s['genome_size_min']:,} - {s['genome_size_max']:,} bp")
    lines.append(f"  Estimated heterozygosity: {pct(s.get('heterozygosity_min'))} - "
                 f"{pct(s.get('heterozygosity_max'))}")
    lines.append(f"  Model fit (full model): {pct(s.get('model_fit_full', s.get('model_fit_max')))}")
    lines.append(f"  K-mers modelled (all k-mers): "
                 f"{pct(s.get('model_fit_allkmers', s.get('model_fit_min')))}")
    lines.append("")

    lines.append("SEQUENCING STATISTICS (FastP):")
    lines.append(f"  Reads after filtering: {num(s.get('total_reads_after_filtering'))}")
    lines.append(f"  Bases after filtering: {num(s.get('total_bases_after_filtering'))} bp")
    lines.append(f"  Data retention: {num(s.get('data_retention_percent'), '{:.1f}')}%")
    lines.append(f"  Theoretical coverage before filtering: "
                 f"{num(s.get('coverage_before_filtering'), '{:.1f}')}x")
    lines.append(f"  Theoretical coverage after filtering: {coverage:.1f}x")
    lines.append("")

    peak, kmercov, lam = s.get("kmer_peak_coverage"), s.get("kmercov"), s.get("lambda_depth")
    lines.append("K-MER COVERAGE:")
    lines.append(f"  K-mer peak coverage: {f'{peak}x' if peak else 'none detected'}")
    lines.append(f"  Fitted k-mer coverage (GenomeScope kmercov): "
                 f"{f'{kmercov:.1f}x' if kmercov else 'NA'}")
    if final:
        lines.append(f"  Measured k-mer coverage (BUSCO read depth): "
                     f"{f'{lam:.1f}x' if lam else 'NA'}")
    lines.append("")

    if final:
        lines.append("ASSEMBLY CROSS-CHECK:")
        lines.append(f"  Decontaminated assembly size: {num(s.get('assembly_size'))} bp")
        lines.append(f"  Compared against: {s.get('assembly_ratio_basis') or 'NA'} assembly "
                     f"({num(s.get('assembly_ratio_denominator'))} bp)")
        fraction = s.get("host_assembly_fraction")
        lines.append(f"  Host fraction of assembly (lower bound): "
                     f"{f'{fraction * 100:.1f}%' if fraction is not None else 'NA'}")
        lines.append(f"  BUSCO complete: {num(s.get('busco_complete'), '{:.1f}')}%")
        lines.append("")

    lines.append("COVERAGE ASSESSMENT:")
    lines.append(f"  Genome size estimate reliable: {'yes' if s.get('genome_size_reliable') else 'no'}")
    if flags:
        lines.append(f"  Genome size warnings: {', '.join(flags)}")
    status = s.get("coverage_status") or "NA"
    if status == "UNRELIABLE_GENOME_SIZE_ESTIMATE":
        # The number is still reported, but never next to a grade it does not support.
        lines.append(f"  Status: {status} (theoretical coverage {coverage:.1f}x, not trustworthy)")
    else:
        lines.append(f"  Status: {status} coverage ({coverage:.1f}x)")
    lines.append(f"  Recommendation: {s.get('coverage_recommendation') or 'NA'}")
    return "\n".join(lines) + "\n"


def choose_reseeded_fit(original, reseeded, lambda_depth, thresholds=None):
    """Which of two GenomeScope fits to keep after reseeding -l with the measured lambda.

    The rule was fixed before the sweep was run, and it is deliberately strict in both
    directions:

      keep the reseed only if the fit within the model's own region improves AND the
      fitted kmercov moves closer to the independently measured lambda.

    Requiring both matters. A reseed that moves kmercov onto the measured value while the
    fit degrades has found a different optimum, not a better one (OG2672 did exactly that:
    kmercov 16.0 -> 7.7 against a measured 6.3, but fit 80.9 -> 62.8). And a size estimate
    that merely lands closer to the assembly proves nothing, because the assembly is what
    the gate's own cross-check measures against -- agreement bought by tuning a starting
    value would be circular. So assembly size is not part of this decision.

    Returns (winner, reason) where winner is "reseeded" or "original".
    """
    cfg = dict(DEFAULTS, **(thresholds or {}))

    base_fit = original.get("summary", {}).get("model_fit_full")
    new_fit = reseeded.get("summary", {}).get("model_fit_full")
    base_kcov = original.get("model", {}).get("kmercov")
    new_kcov = reseeded.get("model", {}).get("kmercov")

    if new_fit is None or base_fit is None:
        return "original", "reseeded fit did not produce a usable model"
    if not lambda_depth:
        return "original", "no measured lambda to judge against"
    if new_kcov is None or base_kcov is None:
        return "original", "missing fitted kmercov"

    fit_gain = new_fit - base_fit
    closer = abs(new_kcov - lambda_depth) < abs(base_kcov - lambda_depth)

    if fit_gain >= cfg["min_reseed_fit_gain"] and closer:
        return "reseeded", (f"model_fit_full {base_fit:.1f} -> {new_fit:.1f} "
                            f"(+{fit_gain:.1f}) and kmercov {base_kcov:.1f} -> {new_kcov:.1f} "
                            f"moved toward the measured {lambda_depth:.1f}")
    if not closer:
        return "original", (f"reseeded kmercov {new_kcov:.1f} is no closer to the measured "
                            f"{lambda_depth:.1f} than {base_kcov:.1f}")
    return "original", (f"model_fit_full gained only {fit_gain:+.1f} points, "
                        f"below the {cfg['min_reseed_fit_gain']:g} required")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary")
    parser.add_argument("--model")
    parser.add_argument("--histogram")
    parser.add_argument("--assembly-size", type=int)
    parser.add_argument("--recheck-assembly", action="store_true",
                        help="add the assembly cross-check to an existing coverage summary")
    parser.add_argument("--coverage-json", help="with --recheck-assembly: the summary to update")
    parser.add_argument("--seqkit-stats", help="with --recheck-assembly: seqkit stats for the assembly")
    parser.add_argument("--kmer-depth-json", help="with --recheck-assembly: MEASURE_KMER_COVERAGE output")
    parser.add_argument("--busco-summary", help="with --recheck-assembly: BUSCO short_summary JSON")
    parser.add_argument("--reseed-decision",
                        help="with --recheck-assembly: GENOMESCOPE_RESEED's decision JSON")
    parser.add_argument("--chosen-summary",
                        help="with --recheck-assembly: the GenomeScope summary the reseed chose")
    parser.add_argument("--chosen-model",
                        help="with --recheck-assembly: the GenomeScope model the reseed chose")
    parser.add_argument("--render-report", metavar="COVERAGE_JSON",
                        help="print the *_sequencing_coverage.txt for a coverage summary")
    for name, value in DEFAULTS.items():
        parser.add_argument(f"--{name.replace('_', '-')}", type=type(value), default=value)
    args = parser.parse_args(argv)
    thresholds = {name: getattr(args, name) for name in DEFAULTS}

    if args.render_report:
        with open(args.render_report) as handle:
            sys.stdout.write(render_text_report(json.load(handle)))
        return 0

    if args.recheck_assembly:
        if not args.coverage_json or not args.seqkit_stats:
            parser.error("--recheck-assembly needs --coverage-json and --seqkit-stats")
        with open(args.coverage_json) as handle:
            coverage = json.load(handle)
        winner = None
        if args.reseed_decision:
            with open(args.reseed_decision) as handle:
                winner = json.load(handle).get("winner")
        if winner == "reseeded":
            if not args.chosen_summary:
                parser.error("the reseed decision chose the refit but no --chosen-summary was given")
            coverage = apply_reseeded_fit(
                coverage, parse_summary(args.chosen_summary),
                parse_model(args.chosen_model) if args.chosen_model else {},
                read_histogram(args.histogram, args.peak_search_max) if args.histogram else None,
                thresholds)
        elif winner:
            coverage = dict(coverage, genomescope_fit="original")
        kmer_depth = None
        if args.kmer_depth_json:
            with open(args.kmer_depth_json) as handle:
                kmer_depth = json.load(handle)
        busco_complete = (busco_complete_from_summary(args.busco_summary)
                          if args.busco_summary else None)
        merged = recheck_against_assembly(
            coverage, assembly_size_from_seqkit(args.seqkit_stats), thresholds, kmer_depth,
            busco_complete)
        json.dump(merged, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    if not args.summary:
        parser.error("--summary is required")
    summary = parse_summary(args.summary)
    model = parse_model(args.model) if args.model else {}
    hist = read_histogram(args.histogram, args.peak_search_max) if args.histogram else None
    _flags, details = evaluate(summary, model, hist, args.assembly_size, thresholds)
    json.dump(details, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
