#!/usr/bin/env Rscript

# Import libraries
# ----------------

suppressPackageStartupMessages({
library(optparse)
library(readr)
})

options(warn=-1)

# Parse arguments                                        
# ---------------

option_list = list(
  make_option(c("--pheno"), action="store", type='character',
              help="Input phenofile to be transformed."),
  make_option(c("--transform"), action="store", type='character',
              help="Input transformation table specifying transformations and association tests."),
  make_option(c("--out_prefix"), action="store", type='character',
              help="Prefix for output files."),
  make_option(c("--remove_pheno_na"), action="store_true", default=FALSE,
              help="If specified, do not include samples with NA in the phenotype column in the output.")
)

args = parse_args(OptionParser(option_list=option_list))

pheno_file = args$pheno
transform_file = args$transform
out_prefix = args$out_prefix

# Custom Transform Functions
# --------------------------

ranknorm <- function(u, k = 0.375) {
  # based on RNOmni package https://rdrr.io/cran/RNOmni/src/R/RN.R
  
  # Observations.
  n <- sum(!is.na(u))
  # Ranks.
  r <- rank(u, na.last="keep")
  # Apply transformation.
  out <- qnorm((r - k) / (n - 2 * k + 1))
  return(out)
}

binarise <- function(u, cases) {
  # if the value is in cases: 2, else: 1
  binned <- (u %in% cases)*1 + 1
  return(binned)
}


# Input Files
# -----------

pheno_df <- readr::read_tsv(pheno_file, show_col_types = F)
names(pheno_df)[names(pheno_df) == 'FID'] <- '#FID'

transform_df <- readr::read_tsv(transform_file, show_col_types = F)


# Do Transformations and Output
# -----------------------------

for (i in 1:nrow(transform_df)) {
  r <- transform_df[i,]
  
  out_df <- pheno_df[,c("#FID","IID")]
  ignore <- c()
  
  # apply transforms
  for (col in colnames(r)[-(1:2)]){
    tr <- r[[col]]
    
    if (is.na(tr)) {
      # no transform
      out_df[[col]] <- pheno_df[[col]]
    
    } else if (tr == "IGNORE") {
      # don't include in output df
      ignore <- c(ignore, col)
      
    } else if (tr == "INT") {
      # apply inverse normal transform (ranknorm)
      out_df[[col]] <- ranknorm(pheno_df[[col]])
  
    } else if (tr == "LOG2") {
      # apply log2 transform
      out_df[[col]] <- log2(pheno_df[[col]])
    
    } else {
      # binarise categorical variable
      cases <- strsplit(as.character(tr),',')[[1]]
      out_df[[col]] <- binarise(pheno_df[[col]], cases)
      
    }
  }
  
  # reorder cols (3rd column is always the phenotype to test)
  col_order <- colnames(r)[-(1:2)]
  col_order <- col_order[!(col_order %in% ignore)]
  col_order <- col_order[col_order != r$test]
  col_order <- c("#FID", "IID", r$test, col_order)
  
  out_df <- out_df[, col_order]

  # remove samples with NA in the phenotype column
  if (args$remove_pheno_na){
    out_df <- out_df[!is.na(out_df[[3]]), col_order]
  }

  # recode phenotype column to 1/2 if currently 0/1
  phenovals <- sort(unique(out_df[[3]]))
  if (identical(phenovals,c(0,1))) {
    out_df[[3]] <- binarise(out_df[[3]], 1)
  }
  
  # write output
  write.table(out_df, file = paste0(out_prefix, r$run_id, ".phe" ), 
              sep = '\t', quote = F, row.names = F) 
  
}
