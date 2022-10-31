#!/usr/bin/env bash

################################################################################
# Script to collect GWAS summary statistics from EBI GWAS catalog
#
# Usage: bash ebi_fetch.sh <study_id>
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

# Collecting all studies in TSV format
curl -k "https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative" | tr -d '\r' > all_studies.tsv

# Query all_studies.json
# Check for the study Id to directly collect it
if [[ $STUDY != "NA" ]] || [[ $STUDY != "" ]]; then
    echo $STUDY > selected_studies.txt
else
    echo "[ERROR] Please, include the study Id in the input file"
    exit 1
fi

# Use the obtained ids to download their GWAS summary statistics VCF files.
touch downloaded_studies.txt
while IFS="" read -r gwas_id || [ -n "$gwas_id" ]; do
    echo "
    Downloading study $gwas_id...
    "
    # Collecting numeric part of the id
    NUMBER=$(echo $gwas_id | sed 's/[^0-9]*//g')
    DOWN_ID="GCST"${NUMBER::-3}"001"
    THOUSANDS=$(echo $gwas_id | rev | cut -c4)
    if [ $THOUSANDS -eq 9 ]; then
        UP_ID="GCST"${NUMBER::-5}"10000"
    else
        THOUSANDS=$((THOUSANDS+1))
        UP_ID="GCST"${NUMBER::-4}"$THOUSANDS""000"
    fi
    wget -r -l1 -np -nH --cut-dirs 7 "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/$DOWN_ID-$UP_ID/$gwas_id/harmonised/" -A ".h.tsv.gz"
    # Some studies does not have harmonised data yet.
    N_ROWS=$(zcat *$gwas_id*.h.tsv.gz | wc -l)
    if [ $N_ROWS -le 1 ]; then
        echo "
        [WARNING] Study $gwas_id does not have harmonised data. Deleting downloaded file...
        "
        rm *$gwas_id*.h.tsv.gz
    else
        echo "
        [MSG] Study $gwas_id harmonised data downloaded correctly.
        "
        echo $gwas_id >> downloaded_studies.txt
    fi
done < selected_studies.txt

# Transform the INITIAL SAMPLE SIZE column to contain a single number
# The following line performs:
#   1. Take only INITIAL SAMPLE SIZE column (#9)
#   2. Performs the following substitutions:
#       a. Any non-numeric character or non-space char -> "" (removed)
#       b. Any consecutive number of space chars -> " " (to single space char)
#       c. Strip the trailing space char
#   3. Use AWK to print the column name and to sum all the fields per row
cut -f 9 all_studies.tsv | sed -e 's/[^0-9| ]*//g' -e 's/[ ]\+/ /g' -e 's/ $//g' | awk '{if (NR==1) print "sample_size"; else {for (i=1; i<=NF; i++) {a+=$i} print a; a=0}}' > sample_size_col.txt
paste all_studies.tsv sample_size_col.txt > all_studies_ss.tsv

# Collecting metadata for the downloaded studies
while IFS="" read -r gwas_id || [ -n "$gwas_id" ]; do
    awk -v gwas_id="$gwas_id" 'BEGIN{OFS="\t";FS="\t"}NR==1{print}$15==gwas_id{print;exit}' all_studies_ss.tsv > "$STUDY"_metadata_ebi.tsv
done < downloaded_studies.txt

# Print reported trait as a separate file to send to omop trait mapping separately
awk -F "\t" 'NR==1{for(i=1;i<=NF;i++){f[$i]=i};next}{print $f["MAPPED_TRAIT"]}' "$STUDY"_metadata_ebi.tsv > reported_traits.tsv
