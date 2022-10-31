#!/usr/bin/env Rscript

################################################################################
# GWAS Sumstats table harmonisation using MungeSumstats
#
# Usage: Rscript munge_gwastable.R \
#            --gwas_table regenie_gwas_table.txt \
#            --meta_table metadata.tsv \
#            --outfile harmonised_gwas_sumstats.vcf \
#            --study_id allancs_regenie
#            --genome_build=GRCh38 \
#            --total_samples=1234 \
#            --ncpus=2
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

option_list=list(
    make_option(c("-s", "--gwas_table"), type="character", default=NULL,
        help="GWAS table input file", metavar="character"),
    make_option(c("-I", "--study_id"), type="character", default=NULL,
        help="Study ID of the GWAS table data", metavar="character"),
    make_option(c("-m", "--meta_table"), type="character", default=NULL,
        help="Input TSV metadata table", metavar="character"),
    make_option(c("-b", "--genome_build"), type="character", default=NULL,
        help="Build of Human genome in study e.g GRCh37|GRCh38", metavar="character"),
    make_option(c("-C", "--total_samples"), type="numeric", default=NULL,
        help="Sample size of the study", metavar="numeric"),
    make_option(c("-o", "--outfile"), type="character", default=NULL,
        help="Name of output file.", metavar="character"),
    make_option(c("-M", "--col_miss"), type="numeric", default="0.1",
        help="Missingness percentage allowed for a column.", metavar="numeric"),
    make_option(c("-t", "--method"), type="character", default="MSS",
        help="Sumstats munging method. Valid options: MSS, simple.", metavar="MSS | simple"),
    make_option(c("--dbsnp"), type="numeric", default="155",
        help="Version of dbSNP to use for munging. Valid options: 144, 155", metavar="144 | 155"),
    make_option(c("-d", "--field_descriptions"), type="character", default="field_descriptions.tsv",
        help="TSV file of descriptions for VCF header fields", metavar="character"),
    make_option(c("-N", "--ncpus"), type="numeric",
        default="1", help="Number of cpus to use in any multithreaded operations within the script", metavar="numeric"),
    make_option(c("--metagwas"), default=FALSE,
        help="Whether the data run has come from a metagwas study")
)

opt_parser=OptionParser(option_list=option_list)
opt=parse_args(opt_parser)

if (!(opt$method %in% c('MSS', 'simple'))) stop(paste("Invalid value for --method :", opt$method, ". Value must be 'MSS' or 'simple'."))

GWAS_SS_COLUMNS <- c(
    "SNPID",
    "SNP",
    "SI",
    "SE",
    "POS",
    "P",
    "N",
    "NC",
    "N_Counts",
    "IS",
    "INFO",
    "ID",
    "GENPOS",
    "CHROM",
    "CHR",
    "BP",
    "BETA",
    "Allele2",
    "Allele1",
    "A2FREQ",
    "A1FREQ",
    "AF_Allele2",
    "AF1",
    "AC",
    "AC_Allele2"
)

miss_percent_allow <- opt$col_miss

# Create meta_table
study_metadata <- as.data.frame(data.table::fread(file=opt$meta_table, sep="\t", sep2=",", blank.lines.skip=TRUE, header=TRUE))
TotalSamples <- c(opt$total_samples)
study_metadata$TotalSamples <- TotalSamples
study_metadata <- harmonise_gwas_table_metadata(study_metadata)
row.names(study_metadata) <- c(opt$study_id)

# Read and prepare sumstats table
gwas_data <- fread(opt$gwas_table, data.table=FALSE)
gwas_data <- gwas_data %>% select(any_of(GWAS_SS_COLUMNS))

# Get number of rows in input
nrows_gwas <- dim(gwas_data)[1]
study_metadata$TotalVariants <- nrows_gwas

# Remove columns with more than threshold fraction NA
allowed_miss <- miss_percent_allow * nrows_gwas
columns_to_remove_na <- sapply(gwas_data, function (x) sum(is.na(x)) > allowed_miss)
gwas_data <- gwas_data[, !columns_to_remove_na]

if (!('N' %in% colnames(gwas_data))) {
    gwas_data$N <- opt$total_samples
}

genome_build=opt$genome_build
# Detect genome build if genome_build is not supplied
if (is.null(genome_build) || genome_build == '') {
    EBI_genome <- MungeSumstats::get_genome_builds(sumstats_list=list(ebi1=gwas_data),
                                                   sampled_snps = 10000,
                                                   dbSNP=opt$dbsnp,
                                                   nThread = opt$ncpus)[[1]]

    genome_build <- EBI_genome
    if (EBI_genome == "GRCH38") genome_build <- "GRCh38"
    if (EBI_genome == "GRCH37") genome_build <- "GRCh37"
}

if (opt$method == "MSS") {
    # Trim alleles
    trimmed <- trim_variants(gwas_data, a1_col="Allele1", a2_col="Allele2", bp_col="POS")
    gwas_data$Allele1 <- trimmed$a1
    gwas_data$Allele2 <- trimmed$a2
    gwas_data$POS <- trimmed$bp

    # impute missing P vals from beta and SE
    gwas_data$P <- impute_p(gwas_data, beta_col='BETA', se_col='SE', p_col='P')

    # Run sumstats harmonisation
    MSS_results <- MungeSumstats::format_sumstats(path=gwas_data,
        ref_genome=genome_build,
        INFO_filter=0,
        pos_se=FALSE,
        N_dropNA=FALSE,
        impute_beta=TRUE,
        impute_se=TRUE,
        bi_allelic_filter=FALSE, # This unhelpful filter removes anything non - biallelic in dbSNP so keep FALSE
        allele_flip_frq=FALSE, # This needs to be FALSE when above filter is FALSE
        rmv_chr=NULL,
        rmv_chrPrefix=TRUE,
        dbSNP=opt$dbsnp,
        return_data=TRUE,
        return_format="data.table",
        imputation_ind=TRUE,
        nThread=opt$ncpus,
        log_folder='mss_log',
        log_mungesumstats_msgs=FALSE,
        log_folder_ind=TRUE)
    rm(gwas_data)
    sumstats_dt <- MSS_results$sumstats

    # Workaround MT chromosome being replaced by NA(fixed in MungeSumstas >= 5.12)
    sumstats_dt <- sumstats_dt %>%
        mutate(CHR=replace(CHR, is.na(CHR) & BP < 16570, "MT"))

    # Fill FRQ with NA if not present
    if (is.null(sumstats_dt$FRQ)) {
        sumstats_dt$FRQ <- NA
    }

    # Find multiallelic sites
    sumstats_dt <- sumstats_dt %>%
        group_by(CHR, BP) %>%
        mutate(multiallelic=n() > 1) %>%
        ungroup()

    # Get inverse frequency when flipping allele except for multiallelic
    if ('flipped' %in% colnames(sumstats_dt)) {
        sumstats_dt <- sumstats_dt %>%
            mutate(flipped=replace(flipped, is.na(flipped), FALSE)) %>%
            mutate(FRQ=replace(FRQ, flipped == TRUE, 1 - FRQ[flipped == TRUE])) %>%
            mutate(FRQ=replace(FRQ, multiallelic == TRUE && flipped == TRUE, NA))
    }

    # Remove intermediate columns
    sumstats_dt <- sumstats_dt %>%
        select(-any_of(c('IMPUTATION_SNP', 'IMPUTATION_BETA', 'IMPUTATION_SE', 'flipped', 'convert_n_integer', 'multiallelic')))

} else if (opt$method == "simple") {
    # Rename columns
    cols_renamed <- c(SNPID='SNP', ID='SNP', CHROM='CHR', POS='BP',
        Allele1='A1', Allele2='A2', INFO='INFO',
        A2FREQ='FRQ', A1='FRQ', beta='BETA', odds_ratio='OR', standard_error='SE',
        p_value='P', n_total='N')
    gwas_data <- gwas_data %>%
        rename_with(~recode(., !!!cols_renamed))

    gwas_data$CHR <- as.character(gwas_data$CHR)

    if (is.null(gwas_data$BETA)) {
        gwas_data$BETA <- log(gwas_data$OR)
    }

    if (is.null(gwas_data$FRQ)) {
        gwas_data$FRQ <- NA
    }

    sumstats_dt <- simple_sumstats_munge(gwas_data)
    rm(gwas_data)

    # remove duplicates
    sumstats_dt <- sumstats_dt %>%
        distinct(CHR, BP, A1, A2, .keep_all=TRUE)
}

sumstats_dt <- as(sumstats_dt, 'data.table')

# rename chromosomes: 1 - 22, X -> 23, Y -> 24, M / MT -> 26
sumstats_dt$CHR <- recode(sumstats_dt$CHR, X='23', Y='24', MT='26', M='26')

nrows_harmonised <- dim(sumstats_dt)[1]
study_metadata$HarmonisedVariants <- nrows_harmonised

# Write sumstats as GWAS VCF with metadata
sumstats_dt$study_id <- opt$study_id
write_gwas_vcf(sumstats_dt, study_metadata, genome_build,
    output_filename=opt$outfile, field_descriptions=opt$field_descriptions,
    metagwas=opt$metagwas)
