#!/usr/bin/env Rscript

################################################################################
# Script for running LDpred2
#
# Based on tutorial by the authors:
# https://privefl.github.io/bigsnpr/articles/LDpred2.html
#
# Date: 2022-03-04
################################################################################


############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) ) {
  cat("
     R Script for running LDpred2: ldpred2.R
      Mandatory arguments:
        --plink_prefix         - Plink prefix.
        --gt_format        - Format of genotypic data. Options: hard_called or dosages.
        --validation_test_split     - Split ratio of validation/test data.
        --summary_statistics    - Path to input summary statistics file.
        --output_file           - Path to resultant file with polygenic risk scores.
        --random_seed_split     - Random seed set for test/validation data split.
        --ld_reference          - .rds object containing LD reference data. 
         --help                 - helpful documentation.
     Usage:
          The typical command for running the script is as follows:
          Rscript ldpred2.R --summary_statistics=public-data3-sumstats.txt --plink_prefix=target_dataset --gt_format='hard_called' --validation_test_split='50:50' --output_file='PRS_scores.csv' --random_seed_split=123 --ld_reference='map_hm3_ldpred2.rds' 
     Output:
      Returns a single .csv file with individual polygenic risk scores corresponding to testing dataset.
      See ./ldpred2.R --help for more details.
      \n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

plink_prefix <- as.character(args$plink_prefix)
gt_format <- as.character(args$gt_format)
bed_file <- paste0(plink_prefix, ".bed")
input_sumstats <- as.character(args$summary_statistics)
random_seed_split <- as.numeric(args$random_seed_split)
validation_test_split <- as.character(args$validation_test_split)
output_file <- as.character(args$output_file)


library(runonce)
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(stringr)

# Read from hard genotype calls from bed/bim/fam or dosages from bgen, it generates .bk and .rds files.

if (gt_format == "hard_called") { 
    snp_readBed(bed_file)
    } else if (gt_format == "dosages") {
    snp_readBGEN(bed_file)
    } else {
    stop("LDPred2: input format of individual-level genotype data is not recognised. Please select from hard_called/dosages.")
}

# Attach the "bigSNP" object in R session
rds_file <- paste0(plink_prefix, ".rds")
obj.bigSNP <- snp_attach(rds_file)

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
phenotype   <- obj.bigSNP$fam$affection
NCORES <- nb_cores()

# Read external summary statistics

sumstats <- bigreadr::fread2(input_sumstats)

# Split the genotype data
set.seed(random_seed_split)

validation_test = str_split(validation_test_split, ":")
validation_samples = as.integer(length(phenotype)/100 * as.integer(unlist(validation_test[1])))

ind.val <- sample(nrow(G), validation_samples)
ind.test <- setdiff(rows_along(G), ind.val)

# Matching variants between genotype data and summary statistics
sumstats$n_eff <- sumstats$N
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos

# Computing LDpred2 scores genome-wide

# Use genetic maps available at https://github.com/joepickrell/1000-genomes-genetic-maps/ to interpolate physical positions (in bp) to genetic positions (in cM).
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".", ncores = NCORES) 

# Create the on-disk sparse genome-wide correlation matrix on-the-fly
tmp <- tempfile(tmpdir = ".")

for (chr in 1:22) {
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

# LDpred2-inf: infinitesimal model
(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
cor(pred_inf, phenotype[ind.test])

# LDpred2(-grid): grid of models
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(lm(phenotype[ind.val] ~ x))$coef["x", 3]
})


params %>%
  mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
  arrange(desc(score)) %>%
  mutate_at(c("score", "sparsity"), round, digits = 3) %>%
  slice(1:10)

# Choose the best model regardless sparsity
best_beta_grid <- params %>%
  mutate(id = row_number()) %>%
  # filter(sparse) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(id) %>%
  beta_grid[, .]

pred <- big_prodVec(G, best_beta_grid, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
cor(pred, phenotype[ind.test])

#  LDpred2-auto: automatic model
# Recommended to run many chains in parallel with different initial values for p.
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
                               ncores = NCORES)

auto <- multi_auto[[1]]

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])
final_pred_auto <- rowMeans(pred_auto[, keep])
cor(final_pred_auto, phenotype[ind.test])

# compute predictions for all samples and create result dataframe
betas <- cbind(beta_inf, best_beta_grid,
             beta_auto = final_beta_auto)

pred_all <- big_prodMat(G, betas, ind.col = df_beta[["_NUM_ID_"]])

pred_df <- cbind(obj.bigSNP$fam[,c("family.ID","sample.ID")],
                 setNames(as.data.frame(pred_all), colnames(betas)))
pred_df$test_set <- T
pred_df[ind.val,"test_set"] <- F
pred_df <- pred_df %>%
  rename(FID = family.ID, IID = sample.ID) %>%
  select(FID, IID, test_set, everything())

# save results
res <- list(pred = pred_df, params = params, auto = multi_auto[keep])
saveRDS(res, "res_file.rds")
write.csv(pred_df, output_file, row.names=F, quote=F)

