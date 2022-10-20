#!/usr/bin/env Rscript

################################################################################
# Script for splitting input GWAS summary statistics by chromosome
#
# Date: 2022-03-16
################################################################################


############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) ) {
  cat("
     R Script for splitting GWAS summary statistics by chromosome: split_sumstats_sbayesr.R
      Mandatory arguments:
        --input_sumstats         - Input GWAS summary statistics.
        --p_max            - P-value threshold for filtering.
        --chromosome            - Chromosome.
         --help            - helpful documentation.
     Usage:
          The typical command for running the script is as follows:
          Rscript split_sumstats_sbayesr.R --input_sumstats=public-data3-sumstats.txt --p_max=0.5  --chromosome 22 
     Output:
      Returns a space-delimited file of COJO-format containing summary statistics after filtering by user-specifed p-value and chromosome.
      See ./split_sumstats_sbayesr.R --help for more details.
      \n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

sumstats <- as.character(args$input_sumstats)
chrom <- as.character(args$chromosome)
p_max <- as.numeric(args$p_max)


library(data.table)


GWAS<-fread(sumstats)
GWAS<-GWAS[complete.cases(GWAS),]

    
GWAS_tmp<-GWAS[GWAS$CHR == chrom,]
GWAS_tmp<-GWAS_tmp[order(GWAS_tmp$BP),]

###
# Change to COJO format
###

# If OR present, calculate BETA
if( 'OR' %in% names(GWAS)){
  GWAS$BETA<-log(GWAS$OR)
}

# Rename allele frequency column
if( "FREQ" %in% names(GWAS)){
  GWAS$MAF<-GWAS$FREQ
} else {
  GWAS$MAF<-GWAS$REF.FREQ
}

GWAS<-GWAS[,c('SNP','A1','A2','MAF','BETA','SE','P','N'),with=F]
names(GWAS)<-c('SNP','A1','A2','freq','b','se','p','N')

# Check whether per variant sample size is available
if(length(unique(GWAS$N)) == 1){
  per_var_N<-F
}
  

# Set maximum p-value threshold
if(!is.na(p_max) == T){
  GWAS<-GWAS[GWAS$p <= p_max,]
}

# Remove variants with SE of 0
if(sum(GWAS$se == 0) > 0){
  GWAS<-GWAS[GWAS$se!= 0,]
}

# Write out cojo format sumstats
fwrite(GWAS, paste0('GWAS_sumstats_COJO.chr',chrom,'.txt'), sep=' ', na = "NA", quote=F)
