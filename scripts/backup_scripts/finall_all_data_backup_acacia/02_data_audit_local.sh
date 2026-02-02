#!/bin/bash

# Load in the configfile
. ../configfile.txt

## Use these variables to override configfile settings if needed, otherwise comment them out
# RUN=NOVA_260108_AD
# rundir=/scratch/pawsey0964/tpeirce/NOVA_260108_AD/draftgenomes

#1. Generate txt file of OG numbers in directory called OGLIST.txt
# Set the CSV file name

TSV_WF=draftcheck-workflow-local.$RUN.tsv
# Print the TSV header
echo -e "OGID\tLocNum\tLocSize\tLocBytes" | tee -a $TSV_WF

# Create the OGLIST from the data audit
OGLIST_FILE=$results/OGLIST.$RUN.txt
ls $rundir > $OGLIST_FILE

# Loop through each OGID in the OGLIST file
while IFS= read -r OGID; do
  # Get the sizes from rclone for each location
  SIZES=$(echo $(rclone size $rundir/$OGID | sed 's/Total/|-- Local/g'))
  
  # Format and append the results to the TSV
  echo $OGID $SIZES | sed -E 's/(\(|\))//g' | awk '{print $1,$6,$10 $11,$12;}' | sed 's/ /\t/g' | tee -a $TSV_WF
done < "$OGLIST_FILE"