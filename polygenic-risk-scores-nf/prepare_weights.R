#/usr/bin/env Rscript

################################################################################
# Script for reformatting MEGAPRS output variant weights into a format compliant
# plink --score function
#
#
# Date: 2022-03-14
################################################################################


############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) ) {
  cat("
     R Script for reformatting MEGAPRS output: 
      Mandatory arguments:
        --variant_weights         - Variant weights (raw output of MEGAPRS)
        --output                  - File name for output weights.
        --bim_file                - Path to .bim file of target data.
        --help                    - helpful documentation.
     Usage:
          The typical command for running the script is as follows:
          Rscript prepare_weights.R --variant_weights='megaprs_effects.tsv' --bim_file='target_dataset.bim' --output='clean_megaprs_weights' 
     Output:
      Returns a single gzipped file containing variant weights - ready to be used as input for plink --score.
      See ./prepare_weights --help for more details.
      \n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

variant_weights <- as.character(args$variant_weights)
bim_file <- as.character(args$bim_file)
output <- as.character(args$output)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

# load in variant weights
weights <- fread(variant_weights) 
bim <- fread(bim_file)
bim$Predictor<-paste0(bim$V1,':',bim$V4)
score<-merge(weights, bim[,c('Predictor','V2'),with=F], by='Predictor')
score<-score[,c('V2','A1',names(score)[grepl('Model', names(score))]), with=F]
names(score)[1]<-'SNP'
names(score)[grepl('Model', names(score))]<-paste0('SCORE_ldak_',names(score)[grepl('Model', names(score))])

fwrite(score, paste0(output,'.weights'), col.names=T, sep=' ', quote=F)

system(paste0('gzip ',output,'.weights'))