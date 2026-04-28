#!/bin/bash --login

#SBATCH --account=pawsey0812
#SBATCH --job-name=final-genome-upload
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=60:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=

. ../configfile.txt

S3=s3://oceanomics/OceanGenomes/analysed-data/draft-genomes
ARCHIVE_FAILED=0

for i in "$rundir"/*; do
    [ -d "$i" ] || continue

    OG=$(basename "$i")
    echo "Processing $OG"

    rclone copy "$i/fastp/" "$S3/$OG/" \
        --include "*.fastq.gz" \
        --include "*fastp.json" \
        --checksum -P

    rclone copy "$i/assemblies/genome/" "$S3/$OG/" \
        --include "*.fna" \
        --checksum -P

    if [ ! -d "$i/kmers" ]; then
        echo "No kmers directory found for $OG"
        continue
    fi

    MERYL_DIR=$(find "$i/kmers" -maxdepth 1 -type d -name "*.meryl" -print -quit)

    if [ -z "$MERYL_DIR" ]; then
        EXISTING_TAR=$(find "$i/kmers" -maxdepth 1 -type f -name "*.meryl.tar.gz" -print -quit)

        if [ -n "$EXISTING_TAR" ] && tar -tzf "$EXISTING_TAR" > /dev/null 2>&1; then
            echo "Valid existing meryl archive found for $OG"
        else
            echo "No .meryl directory found for $OG"
        fi

        continue
    fi

    MERYL_BASENAME=$(basename "$MERYL_DIR")
    TAR_FILE="${MERYL_DIR}.tar.gz"
    TMP_TAR_FILE="${TAR_FILE}.tmp"
    ARCHIVE_NEEDED=true

    if [ -f "$TAR_FILE" ]; then
        echo "Checking existing archive: $TAR_FILE"

        if tar -tzf "$TAR_FILE" > /dev/null 2>&1; then
            echo "Existing archive is valid. Removing original .meryl directory."
            rm -rf "$MERYL_DIR"
            ARCHIVE_NEEDED=false
        else
            echo "Existing archive is invalid. Recreating."
            rm -f "$TAR_FILE"
        fi
    fi

    if [ "$ARCHIVE_NEEDED" = true ]; then
        rm -f "$TMP_TAR_FILE"

        echo "Creating archive: $TAR_FILE"

        if tar -czf "$TMP_TAR_FILE" -C "$(dirname "$MERYL_DIR")" "$MERYL_BASENAME" &&
           tar -tzf "$TMP_TAR_FILE" > /dev/null 2>&1 &&
           mv -f "$TMP_TAR_FILE" "$TAR_FILE" &&
           tar -tzf "$TAR_FILE" > /dev/null 2>&1; then
            rm -rf "$MERYL_DIR"
            echo "Archive complete and .meryl directory removed."
        else
            echo "ERROR: archive failed for $OG. Keeping .meryl directory."
            rm -f "$TMP_TAR_FILE"
            ARCHIVE_FAILED=1
        fi
    fi
done

if [ "$ARCHIVE_FAILED" -ne 0 ]; then
    echo "ERROR: one or more meryl archives failed. Final rclone move not run."
    exit 1
fi

bash 02_data_audit_local.sh

rclone move "$rundir" pawsey0964:oceanomics-draftgenomes/genomes.v2/ \
    --exclude "**/fastp/**" \
    --checksum -P

bash 03_data_audit_acacia.sh