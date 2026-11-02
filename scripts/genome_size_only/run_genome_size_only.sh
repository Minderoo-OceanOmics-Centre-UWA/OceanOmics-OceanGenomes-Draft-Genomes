#!/usr/bin/env bash
# Re-check the GenomeScope genome size of drafts the depth-only run measured, and upload it
# (-entry GENOME_SIZE_ONLY).
#
# Usage, from the repository root on Setonix:
#   scripts/genome_size_only/run_genome_size_only.sh <list.tsv> [--skip_upload_results]
#
# <list.tsv> has og_id and seq_date per line, tab-separated, no header. Every draft must
# already have a depth-only result (scripts/depth_only/run_depth_only.sh) with status ok;
# one without is written to skipped.tsv.
#
# For each draft this stages, into the depth-only stage directory
# $DEPTH/stage/<og>.<seq_date>/, the small files the re-check needs. Files the depth-only run
# already staged are skipped by rclone --checksum. No BAM or assembly is fetched.
#   Acacia  pawsey0964:oceanomics-draftgenomes/genomes.v2/<og>/
#           meryl histogram, GenomeScope summary and model, BUSCO full table and short
#           summary, seqkit stats
#   S3      s3:ocom.oceangenomes/analysed-data/draft-genomes/<og>/   the fastp report
# and takes the *_kmer_depth.json, and seqkit stats where Acacia has none, from
# $DEPTH/results/draftgenomes/<og>/. A draft whose depth was measured by a full run (260909
# on) has its *_kmer_depth.json on Acacia instead, and that copy is used when the depth-only
# results have none.
#
# Results publish to $BASE/results/draftgenomes/<og>/{coverage,kmers}/. Set GSIZE_BASE to
# use another directory, e.g. to keep one batch's results apart while another is archived.

set -euo pipefail

LIST="${1:?usage: $0 <list.tsv> [extra nextflow args]}"
shift
module load nextflow/25.04.6
module load singularity/4.1.0-nompi

BASE="${GSIZE_BASE:-/scratch/pawsey1348/$USER/genome_size_only}"
DEPTH="/scratch/pawsey1348/$USER/depth_only"
ACACIA="pawsey0964:oceanomics-draftgenomes/genomes.v2"
S3="s3:ocom.oceangenomes/analysed-data/draft-genomes"
RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
BATCH="$(basename "$LIST" .tsv)"
SHEET="$BASE/samplesheets/${BATCH}.csv"
SKIPPED="$BASE/samplesheets/${BATCH}.skipped.tsv"
mkdir -p "$BASE/samplesheets" "$BASE/results"

# The one file matching a pattern under a directory, or nothing when there are none or several.
one() {
    local hits
    hits=$(find "$1" -path "$2" -type f 2>/dev/null)
    if [ -n "$hits" ] && [ "$(printf '%s\n' "$hits" | wc -l)" -eq 1 ]; then printf '%s' "$hits"; fi
}

# The first match in sorted order: the same choice run_depth_only.sh made of the BUSCO table.
first() {
    find "$1" -path "$2" -type f 2>/dev/null | sort | head -n 1
}

echo "og_id,seq_date,fastp_json,genomescope_summary,genomescope_model,meryl_hist,kmer_depth_json,seqkit_stats,busco_short_summary" > "$SHEET"
: > "$SKIPPED"

while IFS=$'\t' read -r og sd; do
    [ -n "$og" ] || continue
    key="${og}.ilmn.${sd}"
    d="$DEPTH/stage/${og}.${sd}"
    r="$DEPTH/results/draftgenomes/${og}"
    mkdir -p "$d"
    rclone copy "$ACACIA/$og/" "$d/" \
        --include "kmers/${key}.meryl.hist" \
        --include "kmers/${key}.genomescope/${key}_genomescope2_summary.txt" \
        --include "kmers/${key}.genomescope/${key}_genomescope2_model.txt" \
        --include "assemblies/genome/busco/${key}.*full_table.tsv" \
        --include "assemblies/genome/busco/${key}.*short_summary.json" \
        --include "assemblies/genome/seqkit/*${key}*" \
        --include "coverage/${key}_kmer_depth.json" \
        --checksum
    rclone copy "$S3/$og/" "$d/fastp/" --include "*${sd}*fastp.json" --checksum

    hist=$(one "$d" "*/kmers/${key}.meryl.hist")
    summary=$(one "$d" "*/${key}_genomescope2_summary.txt")
    model=$(one "$d" "*/${key}_genomescope2_model.txt")
    fastp=$(one "$d" "*/fastp/*${sd}*fastp.json")
    depth=$(one "$r" "*/coverage/${key}_kmer_depth.json")
    [ -n "$depth" ] || depth=$(one "$d" "*/coverage/${key}_kmer_depth.json")
    stats=$(one "$d" "*/seqkit/*${key}*seqkit_stats.tsv")
    [ -n "$stats" ] || stats=$(one "$r" "*${key}*seqkit_stats.tsv")
    # The depth-only run published its stats as null.seqkit_stats.tsv until 2026-11-02. Take
    # that file only if it is about this draft's assembly.
    if [ -z "$stats" ]; then
        stats=$(one "$r" "*/seqkit/null.seqkit_stats.tsv")
        if [ -n "$stats" ] && ! grep -q "${key}" "$stats"; then stats=""; fi
    fi
    # The short summary of the lineage whose full table anchored the depth.
    table=$(first "$d" "*/busco/${key}.*full_table.tsv")
    busco=""
    if [ -n "$table" ] && [ -f "${table%full_table.tsv}short_summary.json" ]; then
        busco="${table%full_table.tsv}short_summary.json"
    fi

    missing=""
    for pair in "meryl_hist:$hist" "genomescope_summary:$summary" "genomescope_model:$model" \
                "fastp_json:$fastp" "kmer_depth_json:$depth" "seqkit_stats:$stats" "busco_short_summary:$busco"; do
        [ -n "${pair#*:}" ] || missing="${missing}${pair%%:*} "
    done
    if [ -n "$missing" ]; then
        printf '%s\t%s\tmissing or ambiguous: %s\n' "$og" "$sd" "$missing" >> "$SKIPPED"
        continue
    fi
    if ! grep -q '"status": "ok"' "$depth"; then
        printf '%s\t%s\tdepth status is not ok\n' "$og" "$sd" >> "$SKIPPED"
        continue
    fi
    echo "$og,$sd,$fastp,$summary,$model,$hist,$depth,$stats,$busco" >> "$SHEET"
done < "$LIST"

echo "samplesheet: $SHEET ($(($(wc -l < "$SHEET") - 1)) drafts)"
echo "skipped:     $SKIPPED ($(wc -l < "$SKIPPED") drafts)"

cd "$BASE"
nextflow -log ".nextflow_genome_size_${BATCH}.log" \
    run "$RUN_DIR/main.nf" \
    -entry GENOME_SIZE_ONLY \
    -resume \
    -name "gsize_${BATCH}_$(date +%Y%m%d_%H%M%S)" \
    -work-dir "./work/$BATCH" \
    -c "$RUN_DIR/pawsey_profile.config" \
    -profile singularity \
    --genome_size_samplesheet "$SHEET" \
    --outdir "$BASE/results" \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    "$@"
