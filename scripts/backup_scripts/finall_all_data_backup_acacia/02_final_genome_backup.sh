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
#-----------------
#Loading the required modules

. ../configfile.txt

# Zip the meryl database before backing up.
for i in  "$rundir"/*; do
    if [ -d "$i" ]; then  # Check if the item in $rundir is a directory
        OG=$(basename "$i")

        # Copy the fastq files and final assembly to the s3 bucket
        S3=s3://oceanomics/OceanGenomes/analysed-data/draft-genomes
        rclone copy "$i/fastp/" "$S3/$OG"/ --include "*.fastq.gz" --checksum -P
        rclone copy "$i/assemblies/genome/" "$S3/$OG"/ --include "*.fna" --checksum -P

        # Find the .meryl directory using 'find'
        MERYL_DIR=$(find "$i/kmers" -maxdepth 1 -type d -name "*.meryl" | head -n 1)
        echo "meryl $MERYL_DIR"
        
        if [ -d "$MERYL_DIR" ]; then  # Now checking if it's a directory
            MERYL_BASENAME=$(basename "$MERYL_DIR")  # Just the name without path
            
            echo "Archiving $MERYL_DIR into ${MERYL_DIR}.tar.gz"
            
            # tar the directory
            tar -czvf "${MERYL_DIR}.tar.gz" -C "$(dirname "$MERYL_DIR")" "$MERYL_BASENAME"
        else
            echo "No .meryl directory found in $i/kmers/"
        fi

        wait
        
        rm -rf $MERYL_DIR

        wait

        bash 01_data_audit_local.sh
    fi
done

rclone move $rundir pawsey0964:oceanomics-draftgenomes/genomes.v2/  --checksum -P