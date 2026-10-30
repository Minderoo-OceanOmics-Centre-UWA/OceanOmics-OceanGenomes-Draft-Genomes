#!/usr/bin/env bash
# Measure read depth for archived drafts that have none, and upload it (-entry DEPTH_ONLY).
#
# Usage, from the repository root on Setonix:
#   scripts/depth_only/run_depth_only.sh <list.tsv> [--skip_upload_results]
#
# <list.tsv> has og_id and seq_date per line, tab-separated, no header: one batch, e.g. the
# 251210 lines of the list in genomes_db_review docs/r6_gate_a0_record.md. Try two drafts first
# and check their coverage/*_kmer_depth.json status before running a whole batch.
#
# For each draft this stages, into $BASE/depth_only/stage/<og>.<seq_date>/, the files the
# measurement needs and nothing else:
#   Acacia  pawsey0964:oceanomics-draftgenomes/genomes.v2/<og>/
#           assembly .fna, bwa/*.sorted.bam, BUSCO full_table, GenomeScope model, seqkit stats
#   S3      s3:ocom.oceangenomes/analysed-data/draft-genomes/<og>/   the fastp report
# A draft missing any required file, or with more than one candidate for it, is written to
# skipped.tsv with the reason and left out of the samplesheet; nothing is guessed.
#
# Results publish to $BASE/depth_only/results/draftgenomes/<og>/coverage/. Archive the
# *_kmer_depth.json files back to Acacia beside the draft when the batch is done.

set -euo pipefail

LIST="${1:?usage: $0 <list.tsv> [extra nextflow args]}"
shift
module load nextflow/25.04.6
module load singularity/4.1.0-nompi

BASE="/scratch/pawsey1348/$USER/depth_only"
ACACIA="pawsey0964:oceanomics-draftgenomes/genomes.v2"
S3="s3:ocom.oceangenomes/analysed-data/draft-genomes"
RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
BATCH="$(basename "$LIST" .tsv)"
SHEET="$BASE/samplesheets/${BATCH}.csv"
SKIPPED="$BASE/samplesheets/${BATCH}.skipped.tsv"
mkdir -p "$BASE/stage" "$BASE/samplesheets" "$BASE/results"

# The one file matching a pattern under a directory, or nothing when there are none or several.
one() {
    local hits
    hits=$(find "$1" -path "$2" -type f 2>/dev/null)
    if [ -n "$hits" ] && [ "$(printf '%s\n' "$hits" | wc -l)" -eq 1 ]; then printf '%s' "$hits"; fi
}

# The first match in sorted order. Older drafts can carry BUSCO tables for two lineages (e.g.
# actinopterygii and vertebrata); either one's Complete genes anchor the depth.
first() {
    find "$1" -path "$2" -type f 2>/dev/null | sort | head -n 1
}

echo "og_id,seq_date,assembly,bam,busco_full_table,genomescope_model,fastp_json,seqkit_stats" > "$SHEET"
: > "$SKIPPED"

while IFS=$'\t' read -r og sd; do
    [ -n "$og" ] || continue
    key="${og}.ilmn.${sd}"
    d="$BASE/stage/${og}.${sd}"
    mkdir -p "$d"
    rclone copy "$ACACIA/$og/" "$d/" \
        --include "assemblies/genome/${key}.*.fna" \
        --include "assemblies/genome/bwa/*.sorted.bam" \
        --include "assemblies/genome/busco/${key}.*full_table.tsv" \
        --include "kmers/${key}.genomescope/${key}_genomescope2_model.txt" \
        --include "assemblies/genome/seqkit/*${key}*" \
        --checksum
    rclone copy "$S3/$og/" "$d/fastp/" --include "*${sd}*fastp.json" --checksum

    fna=$(one "$d" "*/assemblies/genome/${key}.*.fna")
    bam=$(one "$d" "*/bwa/*.sorted.bam")
    table=$(first "$d" "*/busco/${key}.*full_table.tsv")
    model=$(one "$d" "*/${key}_genomescope2_model.txt")
    fastp=$(one "$d" "*/fastp/*${sd}*fastp.json")
    stats=$(one "$d" "*/seqkit/*${key}*seqkit_stats.tsv")

    missing=""
    for pair in "assembly:$fna" "bam:$bam" "busco_full_table:$table" "genomescope_model:$model" "fastp_json:$fastp"; do
        [ -n "${pair#*:}" ] || missing="${missing}${pair%%:*} "
    done
    if [ -n "$missing" ]; then
        printf '%s\t%s\tmissing or ambiguous: %s\n' "$og" "$sd" "$missing" >> "$SKIPPED"
        continue
    fi
    echo "$og,$sd,$fna,$bam,$table,$model,$fastp,$stats" >> "$SHEET"
done < "$LIST"

echo "samplesheet: $SHEET ($(($(wc -l < "$SHEET") - 1)) drafts)"
echo "skipped:     $SKIPPED ($(wc -l < "$SKIPPED") drafts)"

cd "$BASE"
nextflow -log ".nextflow_depth_${BATCH}.log" \
    run "$RUN_DIR/main.nf" \
    -entry DEPTH_ONLY \
    -resume \
    -name "depth_${BATCH}_$(date +%Y%m%d_%H%M%S)" \
    -work-dir "./work/$BATCH" \
    -c "$RUN_DIR/pawsey_profile.config" \
    -profile singularity \
    --depth_samplesheet "$SHEET" \
    --outdir "$BASE/results" \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    "$@"
