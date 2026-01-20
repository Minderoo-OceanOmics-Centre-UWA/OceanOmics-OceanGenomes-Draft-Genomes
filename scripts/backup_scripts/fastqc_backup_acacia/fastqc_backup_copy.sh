#!/bin/bash --login
#SBATCH --account=pawsey0964
#SBATCH --job-name=fastqc-backup
#SBATCH --partition=work
#SBATCH --mem=15GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --export=NONE

#Loading the required modules
#module load rclone/1.63.1

. ../configfile.txt
rundir=/scratch/pawsey0964/tpeirce/NOVA_251215_AD/draftgenomes
SOURCE="$rundir"
DEST="/scratch/pawsey0964/tpeirce/BIGRUN/mitogenomes/NOVA_251215_AD/fastp"

for d in "$SOURCE"/OG*/fastp/; do
    [ -d "$d" ] || continue  # skip if no match
    relpath="${d#$SOURCE/}"  # e.g. OG001/fastp
    echo "Copying $d to $DEST/$relpath"
    rclone copy "$d" "$DEST/$relpath" --progress --checksum
done