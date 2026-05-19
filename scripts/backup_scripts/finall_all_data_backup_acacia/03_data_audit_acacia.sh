#!/bin/bash

# Load in the configfile
. ../configfile.txt

## Use these variables to override configfile settings if needed, otherwise comment them out
# RUN=NOVA_260108_AD
# rundir=/scratch/pawsey0964/tpeirce/NOVA_260108_AD/draftgenomes

ACACIA=pawsey0964:oceanomics-draftgenomes/genomes.v2
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
TSV_ACACIA_REMOTE="$results/draftcheck-workflow-acacia.$RUN.tsv"
TSV_S3_REMOTE="$results/draftcheck-workflow-s3.$RUN.tsv"

echo -e "OGID\tAcaciaNum\tAcaciaSize\tAcaciaBytes" > "$TSV_ACACIA_REMOTE"
echo -e "OGID\tS3Num\tS3Size\tS3Bytes" > "$TSV_S3_REMOTE"

while IFS= read -r OGID; do
  [ -n "$OGID" ] || continue

  size_row "$OGID" "Acacia" "$ACACIA/$OGID" >> "$TSV_ACACIA_REMOTE"

  size_row "$OGID" "S3" "$S3/$OGID" \
    --filter "+ *.fastq.gz" \
    --filter "+ *fastp.json" \
    --filter "+ *.fna" \
    --filter "+ *multiqc_report.html" \
    --filter "- **" \
    >> "$TSV_S3_REMOTE"
done < "$OGLIST_FILE"

join -t $'\t' \
  "$results/draftcheck-workflow-local-acacia-subset.$RUN.tsv" \
  "$TSV_ACACIA_REMOTE" \
  | awk -F'\t' 'BEGIN {OFS=FS}
      NR==1 {print $0, "ObjectCountMatch", "ByteCountMatch"; next}
      {print $0, ($2==$5 ? "YES" : "NO"), ($4==$7 ? "YES" : "NO")}' \
  > "$results/_draftcheck-workflow-acacia-join.$RUN.tsv"

join -t $'\t' \
  "$results/draftcheck-workflow-local-s3-subset.$RUN.tsv" \
  "$TSV_S3_REMOTE" \
  | awk -F'\t' 'BEGIN {OFS=FS}
      NR==1 {print $0, "ObjectCountMatch", "ByteCountMatch"; next}
      {print $0, ($2==$5 ? "YES" : "NO"), ($4==$7 ? "YES" : "NO")}' \
  > "$results/_draftcheck-workflow-s3-join.$RUN.tsv"

# Backwards-compatible name for the main Acacia audit.
cp "$results/_draftcheck-workflow-acacia-join.$RUN.tsv" "$results/_draftcheck-workflow-join.$RUN.tsv"

echo "Wrote Acacia audit: $results/_draftcheck-workflow-acacia-join.$RUN.tsv"
echo "Wrote S3 audit: $results/_draftcheck-workflow-s3-join.$RUN.tsv"
