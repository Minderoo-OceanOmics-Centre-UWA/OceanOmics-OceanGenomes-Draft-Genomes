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
. ../configfile.txt

# Or overwrite variables here
# RUN=NOVA_260108_AD2
# download=/scratch/pawsey0964/tpeirce/NOVA_260108_AD2/basespace/NOVA_260108_AD2

rclone copy $download s3://ocom-oceangenomes/illumina-raw/$RUN  --checksum --progress

wait

bash 02_data_audit.sh