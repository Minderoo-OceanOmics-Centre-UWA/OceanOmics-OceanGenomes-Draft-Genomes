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

## Use this variable to override configfile settings if needed, otherwise comment them out
# rundir=/scratch/pawsey0964/tpeirce/NOVA_260108_AD2/draftgenomes


SOURCE="$rundir"
DEST="pawsey0964:oceanomics-filtered-reads"

for d in "$SOURCE"/OG*/fastp/; do
    [ -d "$d" ] || continue  # skip if no match
    relpath="${d#$SOURCE/}"  # e.g. OG001/fastp
    echo "Copying $d to $DEST/$relpath"
    rclone copy "$d" "$DEST/$relpath" --progress --checksum
done

bash 02_fastq_backup_audit.sh