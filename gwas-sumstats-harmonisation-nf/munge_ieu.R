#!/usr/bin/env Rscript

################################################################################
# IEU OpenGWAS harmonisation using MungeSumstats
#
# Usage: Rscript munge_ieu.R <input.vcf> <output.vcf> <N>
#   <input_tsv.> An IEU OpenGWAS summary statistics file (harmonised
#       using their protocol).
#   <output.vcf> The output file. Will be a harmonised GWAS summary statistics
#       VCF with a fix set of columns.
#   <N> Sample size of the study.
################################################################################

suppressPackageStartupMessages({
library(optparse)
library(data.table)
library(MungeSumstats)
library(VariantAnnotation)
library(dplyr)
})

source('harmonise_metadata_funcs.R')
source('munge_funcs.R')

option_list = list(
      make_option(c("-s", "--gwas_vcf"), type="character", default=NULL,
            help="", metavar="character"),
      make_option(c("-m", "--meta_table"), type="character", default=NULL,
            help="", metavar="character"),
      make_option(c("-I", "--study_id"), type="character", default=NULL,
            help="", metavar="character"),
      make_option(c("-o", "--output"), type="character", default=NULL,
            help="", metavar="character"),
      make_option(c("-M", "--col_miss"), type="numeric", default="0.1",
            help="", metavar="numeric"),
      make_option(c("-t", "--method"), type="character", default="MSS",
            help="Sumstats munging method. Valid options: MSS, simple.", metavar="MSS | simple"),
      make_option(c("--dbsnp"), type="numeric", default="155",
            help="", metavar="numeric"),
      make_option(c("-d", "--field_descriptions"), type="character", default="field_descriptions.tsv",
            help="", metavar="character"),
      make_option(c("-N", "--ncpus"), type="numeric", default="1",
            help="", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (!(opt$method %in% c('MSS', 'simple'))) stop(paste("Invalid value for --method :", opt$method, ". Value must be 'MSS' or 'simple'."))

miss_percent_allow <- opt$col_miss
study_id <- opt$study_id

# read IEU study metadata
study_metadata <- data.table::fread(opt$meta_table)

# get vcf header
orig_vcf_header <- VariantAnnotation::scanVcfHeader(opt$gwas_vcf)
orig_metadata <- meta(orig_vcf_header)$SAMPLE

# get genome before harmonising metadata
IEU_genome <- strsplit(study_metadata$build, "/")[[1]][2]

# harmonise IEU study metadata
study_metadata <- harmonise_ieu_metadata(study_metadata)
row.names(study_metadata) <- c(study_id)

# Add studytype to final metadata
study_metadata$StudyType <- orig_metadata[1,"StudyType"]

# read and prepare sumstats table
Rsamtools::indexTabix(opt$gwas_vcf, format="vcf")
read_vcf_threads <- min(opt$ncpus, length(Rsamtools::headerTabix(opt$gwas_vcf)$seqnames))
vcf_dt <- MungeSumstats::read_vcf(opt$gwas_vcf, nThread=read_vcf_threads)
vcf_dt$end <- NULL

# convert to data.frame for compatibility
vcf_dt <- as.data.frame(vcf_dt)

# Get number of rows in input
nrows_gwas <- dim(vcf_dt)[1]
study_metadata$TotalVariants <- nrows_gwas

if ('n_total' %in% colnames(vcf_dt)) {
    vcf_dt$N <- vcf_dt$n_total
    vcf_dt$n_total <- NULL
} else {
    vcf_dt$N <- study_metadata$TotalSamples
}

# Remove columns with more than threshold fraction NA
allowed_miss <- miss_percent_allow * nrows_gwas
columns_to_remove_na <- sapply(vcf_dt, function(x) sum(is.na(x)) > allowed_miss)
vcf_dt <- vcf_dt[, !columns_to_remove_na]

if (opt$method == "MSS") {
    # Trim alleles
    trimmed <- trim_variants(vcf_dt, a1_col = "REF", a2_col = "ALT", bp_col = "BP")
    vcf_dt$REF <- trimmed$a1
    vcf_dt$ALT <- trimmed$a2
    vcf_dt$BP <- trimmed$bp

    # impute missing P vals from beta and SE
    vcf_dt$P <- impute_p(vcf_dt, beta_col='ES', se_col='SE', p_col='P')

    #SI -> INFO
    if("SI" %in% colnames(vcf_dt)){
        vcf_dt$INFO <- vcf_dt$SI
        vcf_dt$SI <- NULL
    }

    # Run sumstats harmonisation
    MSS_results <- MungeSumstats::format_sumstats(path=vcf_dt,
                                                ref_genome = IEU_genome,
                                                INFO_filter = 0,
                                                pos_se = FALSE,
                                                N_dropNA = FALSE,
                                                impute_beta = TRUE,
                                                impute_se = TRUE,
                                                bi_allelic_filter = FALSE, # this unhelpful filter removes anything non-biallelic in dbSNP so keep FALSE
                                                allele_flip_frq = FALSE, # this needs to be FALSE when above filter is FALSE
                                                rmv_chr = NULL,
                                                rmv_chrPrefix = TRUE,
                                                dbSNP = opt$dbsnp,
                                                return_data = TRUE,
                                                return_format = "data.table",
                                                imputation_ind = TRUE,
                                                nThread = opt$ncpus,
                                                log_folder='mss_log',
                                                log_mungesumstats_msgs=TRUE,
                                                log_folder_ind=TRUE)
    rm(vcf_dt)
    sumstats_dt <- MSS_results$sumstats

    # workaround MT chromosome being replaced by NA (fixed in MungeSumstas >= 5.12)
    sumstats_dt <- sumstats_dt %>%
        mutate(CHR = replace(CHR, is.na(CHR) & BP < 16570, "MT"))

    # find multiallelic sites
    sumstats_dt <- sumstats_dt %>%
        group_by(CHR, BP) %>%
        mutate(multiallelic = n() > 1) %>%
        ungroup()

    # Get inverse frequency when flipping allele except for multiallelic
    if ('flipped' %in% colnames(sumstats_dt)) {
    sumstats_dt <- sumstats_dt %>%
        mutate(flipped = replace(flipped, is.na(flipped), FALSE)) %>%
        mutate(FRQ = replace(FRQ, flipped == TRUE, 1-FRQ[flipped == TRUE])) %>%
        mutate(FRQ = replace(FRQ, multiallelic == TRUE && flipped == TRUE, NA))
    }

    # remove intermediate cols
    sumstats_dt <- sumstats_dt %>%
      select(-any_of(c('IMPUTATION_SNP', 'IMPUTATION_BETA', 'IMPUTATION_SE', 'flipped', 'convert_n_integer', 'multiallelic')))

} else if (opt$method == "simple") {
    #rename columns
    cols_renamed <- c(SNP = 'SNP', chr = 'CHR', BP = 'BP',
        REF = 'A1', ALT = 'A2',
        AF = 'FRQ', ES = 'BETA', OR = 'OR', SE = 'SE',
        P = 'P', SS = 'N', N_CAS = 'N_CAS',
        SI = 'INFO', EZ = 'Z')
    vcf_dt <- vcf_dt %>%
        rename_with(~recode(., !!!cols_renamed))

    sumstats_dt <- simple_sumstats_munge(vcf_dt)
    rm(vcf_dt)

    # remove duplicates
    sumstats_dt <- sumstats_dt %>%
        distinct(CHR, BP, A1, A2, .keep_all = TRUE)
}

# remove FILTER col
sumstats_dt <- sumstats_dt %>%
    select(-any_of(c('FILTER')))

sumstats_dt <- as(sumstats_dt, 'data.table')

# rename chromosomes: 1-22, X->23, Y->24, M/MT->26
sumstats_dt$CHR <- recode(sumstats_dt$CHR, X = '23', Y = '24', MT = '26', M = '26')

nrows_harmonised <- dim(sumstats_dt)[1]
study_metadata$HarmonisedVariants <- nrows_harmonised

# Write sumstats as GWAS VCF with metadata
sumstats_dt$study_id <- study_id
write_gwas_vcf(sumstats_dt, study_metadata, IEU_genome,
                output_filename=opt$output, field_descriptions=opt$field_descriptions)
