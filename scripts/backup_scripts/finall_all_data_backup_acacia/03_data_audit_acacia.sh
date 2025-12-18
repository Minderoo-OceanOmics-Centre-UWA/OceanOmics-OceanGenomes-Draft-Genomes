#!/bin/bash

# Load in the configfile
. ../configfile.txt

###########
# Perform check on acacia to compare with the local workflow check
##########
# Set the TSV file name  
TSV=$results/draftcheck-workflow-acacia.$DATE.tsv

#Add headings
echo -e "OGID\tLocNum\tLocSize\tLocBytes" | tee -a $TSV

# Read the OGLIST from the data audit
OGLIST_FILE=$results/OGLIST.txt


# Loop through each OGID in the OGLIST file
while IFS= read -r OGID; do
  # Get the sizes from rclone for each location
  SIZES=$(echo $(rclone size pawsey0964:oceanomics-draftgenomes/genomes.v2/$OGID | sed 's/Total/|-- Acacia/g'))
  
  # Format and append the results to the TSV
  echo $OGID $SIZES | sed -E 's/(\(|\))//g' | awk '{print $1,$6,$10 $11,$12;}' | sed 's/ /\t/g' | tee -a $TSV
done < "$OGLIST_FILE"

###########
# Join the two tables together to make comparison easier
###########
join -t $'\t' draftcheck-workflow-local.$DATE.tsv $TSV > $results/draftcheck-workflow-join.$DATE.tsv