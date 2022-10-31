#!/bin/bash -eu

display_usage() {
    echo "<gwas_table> file needs to be provided as argument."
    echo -e "\nUsage: $0 [arguments] \n"
}

if [[ ( $@ == "--help" ) || $@ == "-h" || $# -ne 1 ]]; then
		display_usage
		exit 1
fi

# Take gwas_table file as input
gwas_table="$1"

method=$(grep -h "#Method" $gwas_table | awk -F "\t" '{print $2}')

# Each table should be transformed to match these conditions:
# GWAS_SS_COLUMNS <- c(
#     "SNPID",
#     "SNP",
#     "SI",
#     "INFO",
#     "SE",
#     "POS",
#     "P",
#     "N",
#     "NC",
#     "N_Cases",
#     "ID",
#     "GENPOS",
#     "CHROM",
#     "CHR",
#     "BP",
#     "BETA",
#     "Allele2",
#     "Allele1",
#     "A2FREQ",
#     "A1FREQ",
#     "AF_Allele2",
#     "AF1",
#     "AC_Allele2",
#     "AC"
# )

if [[ ( $method == "regenie" ) ]]; then
    # Rename GENPOS -> POS, remove last column, add col "P" calculated from col "LOG10P"
    sed '/^##/d' $gwas_table | sed -e '1 s/GENPOS/POS/' \
    -e '1 s/ALLELE0/Allele1/' \
    -e '1 s/ALLELE1/Allele2/' |
    awk '{NF-=1} NR==1{print $0, "P"; next} {print $0, 10^(-$13)}'

elif [[ ( $method == "SAIGE" ) ]]; then
    # Rename columns & take just the standard columns
    sed '/^##/d' $gwas_table | sed -e '1 s/AF_Allele2/A2FREQ/' \
    -e '1 s/N.Cases/NC/' \
    -e '1 s/N.Controls/N_Controls/' \
    -e '1 s/imputationInfo/INFO/' \
    -e '1 s/A0/Allele1/' \
    -e '1 s/A1/Allele2/' \
    -e '1 s/AC_Allele2/AC/' \
    -e '1 s/p.value/P/' -e 's/ /\t/g'

elif [[ ( $method == "PLINK2-GLM" ) ]]; then
    # Rename columns, only take rows where TEST == "ADD", replace ":" -> "_" in snp IDs if present
    # Cut out confusing columns
    sed '/^##/d' $gwas_table |
    sed -e '1 s/#//' \
    -e '1 s/A1_FREQ/A1FREQ/' \
    -e '1 s/ID/SNP/' \
    -e '1 s/REF/Allele1/' \
    -e '1 s/ALT/Allele2/' \
    -e '1 s/OBS_CT/N/' |
    awk 'BEGIN{FS=OFS="\t"} \
        NR==1{for(i=1;i<=NF;i++){c[$i]=i};print $0;next} \
        $c["TEST"]=="ADD"{gsub(":","_",$3); print $0}' |
    cut --complement -d$'\t' -f 6,8,9

elif [[ ( $method == "BOLT-LMM" ) ]]; then
    # Rename P and chisq columns, duplicate final column to P
    sed '/^##/d' $gwas_table |
    sed -r -e '1 s/P_([^\t]+)/\1_P/g' \
    -e '1 s/CHISQ_([^\t]+)/\1_CHISQ/g' \
    -e '1 s/ALLELE0/Allele1/' \
    -e '1 s/ALLELE1/Allele2/' \
    -e '1 s/BP/POS/' |
    awk 'BEGIN{FS=OFS="\t"} \
        NR==1{print $0, "P"; next}{print $0, $NF}'

elif [[ ( $method == "fastGWA-GLMM" ) ]]; then
    # Rename columns
    sed '/^##/d' $gwas_table |
    sed -e '1 s/A2/Allele1/' -e '1 s/A1/Allele2/' -e '1 s/AF1/A2FREQ/'

elif [[ ( $method == "METAL" ) ]]; then
    # Rename columns
    sed '/^##/d' $gwas_table |
    sed -e '1 s/MarkerName/SNP/' \
    -e '1 s/Freq1/A1FREQ/' \
    -e '1 s/P.value/P/' \
    -e '1 s/StdErr/SE/' \
    -e '1 s/BP/POS/' \
    -e '1 s/HetISq/IS/' \
    -e '1 s/Effect/BETA/' |
    cut --complement -d$'\t' -f 5

elif [[ ( $method == "Hail" ) ]]; then
    # Split and rename columns
    echo -e 'CHR\tPOS\tAllele1\tAllele2\tN\tsum_x\ttranspose_x\tBETA\tSE\tt_stat\tP'
    sed '/^##/d' $gwas_table |
    sed '/^locus/d' |
    sed -e 's/:/\t/g' |
    sed -E 's/\[\"([A-Z])\",\"([A-Z])\"\]/\1\t\2/g'

else
    echo "$gwas_table format for $method is unsupported at the moment"
fi
