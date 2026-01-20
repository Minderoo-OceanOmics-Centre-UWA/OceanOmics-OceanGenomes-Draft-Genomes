#!/bin/bash --login
#SBATCH --account=pawsey0964
#SBATCH --job-name=aws-raw-backup
#SBATCH --partition=work
#SBATCH --mem=15GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --export=NONE

#-----------------
# Update the configfile.txt with your AWS and rclone details
# . ../configfile.txt

# Or overwrite variables here
RUN=NOVA_260108_AD
download=/scratch/pawsey0964/tpeirce/NOVA_260108_AD/basespace/NOVA_260108_AD

#rclone copy $download s3:oceanomics/OceanGenomes/illumina-raw/$RUN --checksum

rclone copy $download s3:oceanomics/OceanGenomes/illumina-raw/$RUN  --checksum --progress