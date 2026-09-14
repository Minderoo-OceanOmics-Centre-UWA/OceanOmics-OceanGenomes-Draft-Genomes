#!/usr/bin/env bash
# Snapshot the current NOVA_260909_LA QC metrics before anything is re-run.
#
# The tiers publish into the same outdir as the original run, so the numbers this review
# is based on are overwritten the moment tier 3 starts. This copies only the small
# summary files (a few MB in total, no FASTA, no FASTQ), which is enough for
# 99_compare_to_baseline.sh to produce a before/after table.
#
# Safe to re-run: it refuses to overwrite an existing snapshot.
set -euo pipefail

RUN=NOVA_260909_LA
OUT="/scratch/pawsey1348/$USER/${RUN}"
DEST="${OUT}/rerun_baseline"

if [ -d "$DEST" ]; then
    echo "Baseline already exists at $DEST -- refusing to overwrite it."
    echo "That snapshot is the only record of the pre-fix numbers. Move it aside if you"
    echo "really want a fresh one."
    exit 1
fi
if [ ! -d "$OUT/draftgenomes" ]; then
    echo "ERROR: $OUT/draftgenomes not found." >&2
    exit 1
fi

mkdir -p "$DEST"
n=0
for d in "$OUT"/draftgenomes/OG*; do
    [ -d "$d" ] || continue
    s=$(basename "$d")
    mkdir -p "$DEST/$s"
    # Small summary files only.
    for pat in \
        "kmers/*genomescope/*_summary.txt" \
        "kmers/*genomescope/*_model.txt" \
        "kmers/*merqury.qv" \
        "kmers/*merqury.completeness.stats" \
        "coverage/*_coverage_summary.json" \
        "assemblies/genome/seqkit/*seqkit_stats.tsv" \
        "assemblies/genome/gfastats/*assembly_summary" \
        "assemblies/genome/busco/*short_summary.txt" \
        "assemblies/genome/tiara/*tiara_filter_summary.txt" \
        "assemblies/genome/NCBI/*filter_report.txt" \
        "assemblies/genome/NCBI/*contig_count_500bp.txt" \
        "assemblies/genome/NCBI/*summary.txt"
    do
        for f in $d/$pat; do
            [ -f "$f" ] && cp -p "$f" "$DEST/$s/" 2>/dev/null || true
        done
    done
    n=$((n+1))
done

# The FCS-GX action reports are the evidence for what REVIEW held. They are larger, so
# keep just the action column summary rather than the whole report.
mkdir -p "$DEST/_fcsgx_action_summary"
for d in "$OUT"/draftgenomes/OG*; do
    s=$(basename "$d")
    for r in $d/assemblies/genome/NCBI/*fcs_gx_report.txt; do
        [ -f "$r" ] || continue
        awk -F'\t' '!/^#/{n[$5]++; b[$5]+=$4} END{for(k in n) printf "%s\t%d\t%d\n", k, n[k], b[k]}' \
            "$r" > "$DEST/_fcsgx_action_summary/${s}.tsv"
    done
done

echo "Snapshotted $n samples to $DEST"
du -sh "$DEST"
