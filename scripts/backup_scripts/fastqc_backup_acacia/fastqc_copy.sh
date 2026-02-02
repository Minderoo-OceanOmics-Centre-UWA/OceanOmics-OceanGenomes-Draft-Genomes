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

## This code is for copying fastq files to a location to run the mitogenome nf-core pipeline on them.

. ../configfile.txt
rundir=/scratch/pawsey0964/tpeirce/__mount
SOURCE="$rundir"
DEST=/scratch/pawsey0964/tpeirce/BIGRUN/mitogenomes/rerun_fastp

## Use this chunk for when you are copying files from a specific directory
# for d in "$SOURCE"/OG*/fastp/*.fastq.gz; do
#     [ -d "$d" ] || continue  # skip if no match
#     relpath="${d#$SOURCE/}"  # e.g. OG001/fastp
#     echo "Copying $d to $DEST"
#     rclone copy "$d" "$DEST" --progress --checksum
# done


## Use this chunk when you want to use a specific list of directories to copy
OG_LIST=og_list.txt

while read -r og; do
    # skip empty lines or comments
    [[ -z "$og" || "$og" =~ ^# ]] && continue

    for f in "$SOURCE/$og"/fastp/*.fastq.gz; do
        [ -e "$f" ] || continue  # skip if no fastqs
        relpath="${f#$SOURCE/}"
        echo "Copying $f to $DEST"
        rclone copy "$f" "$DEST" --progress --checksum
    done
done < "$OG_LIST"