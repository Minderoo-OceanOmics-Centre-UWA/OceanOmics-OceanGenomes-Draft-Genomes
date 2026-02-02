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

## Use these variables to override configfile settings if needed, otherwise comment them out
## If overriding ensure to update the data audit scripts also
# RUN=NOVA_260108_AD
# rundir=/scratch/pawsey0964/tpeirce/NOVA_260108_AD/draftgenomes

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

    fi
done

bash 02_data_audit_local.sh

wait 

rclone move $rundir pawsey0964:oceanomics-draftgenomes/genomes.v2/  --checksum -P

wait

bash 03_data_audit_acacia.sh