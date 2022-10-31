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

# Collecting all studies in JSON format using the API
curl -X GET 'https://www.ebi.ac.uk/gwas/api/search/summaryStatistics?q=fullPvalueSet%3Atrue&max=50000&fl=accessionId%2Cauthor_s%2CauthorAscii_s%2CpubmedId%2Ctitle%2Cpublication%2CpublicationDate%2CmappedLabel%2CmappedUri%2CtraitName_s%2CassociationCount%2CagreedToCc0' -H 'Accept: application/json' > all_studies.json

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

# Collecting metadata for the downloaded studies
printf "pubmedId\ttitle\tauthor_s\tauthorAscii_s\tpublication\tpublicationDate\taccessionId\tassociationCount\ttraitName_s\tmappedLabel\tmappedUri\n" > selected_studies_metadata.tsv
while IFS="" read -r gwas_id || [ -n "$gwas_id" ]; do
    cat all_studies.json | jq -c ".response | .docs | .[] | select(.accessionId == \"$gwas_id\") | \"\(.pubmedId);;\(.title);;\(.author_s);;\(.authorAscii_s);;\(.publication);;\(.publicationDate);;\(.accessionId);;\(.associationCount);;\(.traitName_s);;\(.mappedLabel);;\(.mappedUri)\"" >> selected_studies_metadata.tsv
done < downloaded_studies.txt
# Remove quotes, brackets, backslash and change ";;" to \t
sed -e 's/"//g' -e 's/\[//g' -e 's/\]//g' -e 's/\\//g' -e 's/;;/\t/g' selected_studies_metadata.tsv > "$STUDY"_metadata_ebi.tsv
