#!/bin/bash

# Load in the configfile
. ../configfile.txt

## Use these variables to override configfile settings if needed, otherwise comment them out
# RUN=NOVA_260108_AD
# rundir=/scratch/pawsey0964/tpeirce/NOVA_260108_AD/draftgenomes

S3=s3:oceanomics/OceanGenomes/analysed-data/draft-genomes

size_row() {
  local ogid=$1
  local label=$2
  shift 2

  local output
  if ! output=$(rclone size "$@" 2>&1); then
    echo -e "${ogid}\t${label}_ERROR\tNA\tNA"
    echo "ERROR running rclone size for ${label} ${ogid}: ${output}" >&2
    return 0
  fi

  local objects human bytes
  objects=$(echo "$output" | awk -F: '/Total objects/ {gsub(/^[ \t]+/, "", $2); print $2}')
  human=$(echo "$output" | awk -F: '/Total size/ {sub(/^[ \t]+/, "", $2); sub(/ *\(.*/, "", $2); print $2}')
  bytes=$(echo "$output" | awk -F'[()]' '/Total size/ {print $2}' | awk '{print $1}')

  echo -e "${ogid}\t${objects:-0}\t${human:-0 B}\t${bytes:-0}"
}

mkdir -p "$results"

OGLIST_FILE="$results/OGLIST.$RUN.txt"
find "$rundir" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort > "$OGLIST_FILE"

TSV_ACACIA_LOCAL="$results/draftcheck-workflow-local-acacia-subset.$RUN.tsv"
TSV_S3_LOCAL="$results/draftcheck-workflow-local-s3-subset.$RUN.tsv"

echo -e "OGID\tLocalAcaciaNum\tLocalAcaciaSize\tLocalAcaciaBytes" > "$TSV_ACACIA_LOCAL"
echo -e "OGID\tLocalS3Num\tLocalS3Size\tLocalS3Bytes" > "$TSV_S3_LOCAL"

while IFS= read -r OGID; do
  [ -n "$OGID" ] || continue

  # Acacia receives the final draftgenomes tree, excluding fastp outputs.
  size_row "$OGID" "LocalAcacia" "$rundir/$OGID" \
    --exclude "fastp/**" \
    >> "$TSV_ACACIA_LOCAL"

  # S3 receives only the selected final deliverables copied in 01_final_genome_backup.sh.
  size_row "$OGID" "LocalS3" "$rundir/$OGID" \
    --filter "+ fastp/*.fastq.gz" \
    --filter "+ fastp/*fastp.json" \
    --filter "+ assemblies/genome/*.fna" \
    --filter "+ multiqc/*multiqc_report.html" \
    --filter "- **" \
    >> "$TSV_S3_LOCAL"
done < "$OGLIST_FILE"

echo "Wrote local Acacia subset audit: $TSV_ACACIA_LOCAL"
echo "Wrote local S3 subset audit: $TSV_S3_LOCAL"
