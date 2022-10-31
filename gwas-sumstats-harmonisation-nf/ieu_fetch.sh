#!/usr/bin/env bash

################################################################################
# Script to collect GWAS summary statistics from IEU OpenGWAS
#
# Usage: bash ieu_fetch.sh <study_id>
#   <study_id> : a string with the study ID to download
#
# Author: David Piñeyro
# Date: 2022-02-23
################################################################################

STUDY="$1"

if [ -z "$STUDY" ]; then
    echo "[ERROR] No study ID detected."
    exit 1
fi

# Collecting all studies in JSON format using the API
curl -X GET "https://gwas-api.mrcieu.ac.uk/gwasinfo" -H  "accept: application/json" -H  "X-Api-Token: null" > all_meta.json

# convert JSON to TSV
# cat all_meta.json | jq --arg gwas_id "$gwas_id" -r '([.[] | keys] | flatten | unique) as $k | ($k | @tsv), (.[] as $obj | $k | map($obj[.]) | @tsv)' > all_meta.tsv

# Query all_studies.json
# Check for the study Id to directly collect it
if [[ $STUDY != "NA" ]] || [[ $STUDY != "" ]]; then
    echo $STUDY > selected_studies.txt
else
    echo "[ERROR] Please, include the study Id in the input file"
    exit 1
fi

# Use the obtained ids to download their GWAS summary statistics VCF files.
while IFS="" read -r gwas_id || [ -n "$gwas_id" ]; do
    echo "Downloading study $gwas_id..."
    wget https://gwas.mrcieu.ac.uk/files/"$gwas_id"/"$gwas_id".vcf.gz
    echo "Study $gwas_id downloaded"
done < selected_studies.txt

# Collecting metadata for the downloaded studies
while IFS="" read -r gwas_id || [ -n "$gwas_id" ]; do
    cat all_meta.json | jq --arg gwas_id "$gwas_id" -r '([.[] | keys] | flatten | unique) as $k | .[] | select(.id == $gwas_id) as $obj | $k | (@tsv), (map($obj[.]) | @tsv)' > "$STUDY"_metadata_ieu.tsv
done < selected_studies.txt

# Print reported trait as a separate file to send to omop trait mapping separately
awk -F "\t" 'NR==1{for(i=1;i<=NF;i++){f[$i]=i};next}{print $f["trait"]}' "$STUDY"_metadata_ieu.tsv > reported_traits.tsv
