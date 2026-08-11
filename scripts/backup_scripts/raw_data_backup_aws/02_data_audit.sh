#!/bin/bash

# Load in the configfile
. ../configfile.txt

## Use these variables to override configfile settings if needed, otherwise comment them out
# RUN=NOVA_260108_AD
# download=/scratch/pawsey0964/tpeirce/NOVA_260108_AD/basespace/NOVA_260108_AD

S3=s3://ocom-oceangenomes/illumina-raw
## 1. Perform check on local machine 

# Generate txt file with list of OGID into a txt file called OGLIST_RAW.txt
mkdir -p $results
OGLIST_FILE=$results/OGLIST_$RUN.txt
ls $download > $OGLIST_FILE

# Set the TSV file name  
TSV_L=$results/draftRAWcheck-local.$RUN.TSV

# Print the TSV header
echo -e "OGID\tLocNum\tLocSize\tLocBytes" | tee -a $TSV_L

# Loop through each OGID in the OGLIST file
while IFS= read -r OGID; do
  # Get the sizes from rclone for each location
  SIZES_LOCAL=$(echo $(rclone size $download/$OGID | sed 's/Total/|-- Local/g'))
  
  # Format and append the results to the TSV
  echo $OGID $SIZES_LOCAL | sed -E 's/(\(|\))//g' | awk '{OFS="\t"; print $1,$6,$10 $11,$12}' | sed 's/ /\t/g' | tee -a $TSV_L
done < "$OGLIST_FILE"


## 2. Perform check on aws

# Set the TSV file name  
TSV_A=$results/draftRAWcheck-AWS.$RUN.TSV

# Print the TSV header
echo -e "OGID\tAwsNum\tAwsSize\tAwsBytes" | tee -a $TSV_A

# Loop through each OGID in the OGLIST file
while IFS= read -r OGID; do
  # Get the sizes from rclone from AWS
  SIZES_AWS=$(echo $(rclone size $S3/$RUN/$OGID | sed 's/Total/|-- AWS/g'))
  
  # Format and append the results to the TSV
  echo $OGID $SIZES_AWS | sed -E 's/(\(|\))//g' | awk '{OFS="\t"; print $1,$6,$10 $11,$12}' | sed 's/ /\t/g' | tee -a $TSV_A
done < "$OGLIST_FILE"

# join the two for comparison
paste -d '\t' $TSV_L $TSV_A > $results/_draftRAWcheck-join.$RUN.TSV