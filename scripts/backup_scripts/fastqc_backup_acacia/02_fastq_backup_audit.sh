#!/bin/bash

# Load in the configfile
. ../configfile.txt


## Use these variables to override configfile settings if needed, otherwise comment them out
# RUN=NOVA_260108_AD2
# rundir=/scratch/pawsey0964/tpeirce/NOVA_260108_AD2/draftgenomes

# Generate txt file with list of OGID into a txt file called OGLIST_$RUN.txt
mkdir -p $results
OGLIST_FILE=OGLIST_$RUN.txt
echo $rundir
ls $rundir > $OGLIST_FILE

# Set the TSV file name  
TSV_L=$results/draftfastpcheck-local.$DATE.TSV

# Print the TSV header
echo -e "OGID\tLocNum\tLocSize\tLocBytes" > $TSV_L

# Loop through each OGID in the OGLIST file
while IFS= read -r OGID; do
  # Get the sizes from rclone for each location
  SIZES_LOCAL=$(echo $(rclone size $rundir/$OGID/fastp | sed 's/Total/|-- Local/g'))
  
  # Format and append the results to the TSV
  echo $OGID $SIZES_LOCAL | sed -E 's/(\(|\))//g' | awk '{OFS="\t"; print $1,$6,$10 $11,$12}' | sed 's/ /\t/g' >> $TSV_L
done < "$OGLIST_FILE"


## 2. Perform check on acacia

# Set the TSV file name  
TSV_A=$results/draftfastpcheck-acacia.$DATE.TSV

# Print the TSV header
echo -e "OGID\tAwsNum\tacaciaSize\tacaciaBytes" > $TSV_A

# Loop through each OGID in the OGLIST file
while IFS= read -r OGID; do
  # Get the sizes from rclone from AWS
  SIZES_ACACIA=$(echo $(rclone size pawsey0964:oceanomics-filtered-reads/$OGID/fastp | sed 's/Total/|-- acacia/g'))
  
  # Format and append the results to the TSV
  echo $OGID $SIZES_ACACIA | sed -E 's/(\(|\))//g' | awk '{print $1,$6,$10 $11,$12;}' | sed 's/ /\t/g' >> $TSV_A
done < "$OGLIST_FILE"

# Join the two for comparison
paste -d '\t' $TSV_L $TSV_A > $results/_draftfastpcheck-join.$DATE.TSV