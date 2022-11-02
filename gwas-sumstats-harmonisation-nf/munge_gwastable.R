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

harmonise_gwas_table_metadata <- function(input_df) {
    if (!('GeneticModel' %in% colnames(input_df))) {
        input_df$GeneticModel <- NA
    }
    if (!('CovarNames' %in% colnames(input_df))) {
        input_df$CovarNames <- NA
    }
    if (!('PhenoName' %in% colnames(input_df))) {
        input_df$PhenoName <- NA
    }
    out_df <- data.frame(
        TotalCases=as.numeric(input_df[["TotalCases"]]),
        TotalControls=as.numeric(input_df[["TotalControls"]]),
        Population=input_df[["Population"]],
        StudyType=input_df[["StudyType"]],
        Method=input_df[["Method"]],
        TotalSamples=input_df[["TotalSamples"]],
        GeneticModel=input_df[["GeneticModel"]],
        Covariates=input_df[["CovarNames"]],
        TraitReported=input_df[["PhenoName"]])

    return(out_df)
}

# inputs:
# sumstats_dt: dataframe of summary stats with MungeSumStats standard col headers.
#              Required cols - SNP, P, FRQ, BETA, SE, study_id. Optional cols - N, N_CAS, INFO, Z.
#              NB: study_id values should match with rownames in study_metadata.
# study_metadata: dataframe of study-level metadata with 1 row per study and rowname 
#                 corresponding to study ID. Cols will each become a field in '##SAMPLE=<...>' lines.
# genome_build: genome build to add as header line '##genome_build=...'
write_gwas_vcf <- function(sumstats_dt,
                           study_metadata,
                           genome_build,
                           output_filename=NULL,
                           field_descriptions_tsv="field_descriptions.tsv",
                           metagwas=FALSE){


    # set rownames to SNP ID
    if( !"study_id" %in% colnames(sumstats_dt)){
        stop("study_id col is required.")
    }

    #### Convert susmtats to VRanges
    gr <- MungeSumstats:::to_granges(sumstats_dt)
    vr <- VariantAnnotation::makeVRangesFromGRanges(gr,
                                                    ref.field = "A1", 
                                                    alt.field = "A2",
                                                    keep.extra.columns = TRUE,
                                                    sampleNames.field = "study_id")


    # set rownames to SNP ID
    if("SNP" %in% names(mcols(vr))){
        names(vr) <- vr$SNP
        vr$SNP <- NULL
    }else{
        stop("SNP col is required.")
    }

    # transform P to -log10(p)
    if("P" %in% names(mcols(vr))){
        #P -> LP
        vr$LP <- -log(vr$P,base=10)
        vr$P <- NULL
    }else if("LP" %in% names(mcols(vr))) {
        if("P" %in% names(mcols(vr))) vr$P <- NULL
    } else {
        stop("P col or LP col is required.")
    }
    #FRQ -> AF
    if("FRQ" %in% names(mcols(vr))){
        vr$AF <- vr$FRQ
        vr$FRQ <- NULL
    }else{
        stop("FRQ col is required.")
    }
    #BETA -> ES
    if("BETA" %in% names(mcols(vr))){
        vr$ES <- vr$BETA
        vr$BETA <- NULL
    }else{
        stop("BETA col is required.")
    }
    #SE -> SE
    if(!"SE" %in% names(mcols(vr))){
        stop("Standard Error (of effect size) is required.")
    }
    #N -> NS
    if("N" %in% names(mcols(vr))){
        vr$SS <- vr$N
        vr$N <- NULL
    }
    #N_CAS -> NC
    if("N_CAS" %in% names(mcols(vr))){
        vr$NC <- vr$N_CAS
        vr$N_CAS <- NULL
    }
    #INFO -> SI
    if("INFO" %in% names(mcols(vr))){
        vr$SI <- vr$INFO
        vr$INFO <- NULL
    }
    #Z -> EZ
    if("Z" %in% names(mcols(vr))){
        vr$EZ <- vr$Z
        vr$Z <- NULL
    }


    # get field descriptions
    field_descriptions <- as.data.frame(data.table::fread(field_descriptions_tsv))
    row.names(field_descriptions) <- field_descriptions$ID

    outvcf <- VariantAnnotation::asVCF(vr)

    ## change header stuff
    outvcf_header <- header(outvcf)

    # add the FORMAT field header section
    field_header_df <- data.frame(ID=names(mcols(vr)), Number="A", Type="String", Description=".", row.names=names(mcols(vr)))
    specific_field_header_df <- field_descriptions %>% 
        filter(key == "FORMAT" & ID %in% names(mcols(vr))) %>% 
        select(ID, Number, Type, Description)
    
    field_header_df[row.names(specific_field_header_df), ] <- specific_field_header_df

    geno(outvcf_header) <- field_header_df %>% 
        select(Number, Type, Description) %>%
        as("DataFrame")

    vcf_mdata <- meta(outvcf_header)

    vcf_mdata$source <- DataFrame(Value=c('gwas-sumstats-harmonisation-nf'), row.names=c("source"))

    vcf_mdata$GenomeBuild <- DataFrame(Value=c(genome_build), row.names=c("GenomeBuild"))

    if (metagwas){
        vcf_mdata$MetaAnalysis <- DataFrame(Value=c("True"), row.names=c("MetaAnalysis"))
    }

    # add lines describing study metadata fields
    vcf_mdata$META <- field_descriptions %>% 
        filter(key == "META" & ID %in% names(study_metadata) ) %>% 
        select(Number, Type, Description) %>%
        as("DataFrame")

    # add SAMPLE lines containing study metadata for each study
    vcf_mdata$SAMPLE <- study_metadata %>%
        mutate(across(where(is.character), function(x){paste0('"',x,'"')}))

    if (is.null(output_filename)) output_filename <- paste0(row.names(study_metadata)[1], ".gwas.vcf")
    # set the header info and write VCF
    meta(outvcf_header) <- vcf_mdata
    header(outvcf) <- outvcf_header
    suppressWarnings(
        VariantAnnotation::writeVcf(obj = outvcf, filename = output_filename, index = TRUE)
    )
    message(paste("Harmonised GWAS VCF saved to", output_filename, ".bgz"))
}


simple_sumstats_munge <- function(sumstats_df,
                                    a1_col = "A1",
                                    a2_col = "A2",
                                    chr_col = "CHR",
                                    bp_col = "BP",
                                    snp_col = "SNP",
                                    beta_col = "BETA",
                                    n_col = "N",
                                    n_cas_col = "N_CAS") {

    # remove if missing alleles
    filter <- !is.na(sumstats_df[[a1_col]]) & !is.na(sumstats_df[[a1_col]])
    sumstats_df <- sumstats_df[filter, ]

    # remove if missing location
    filter <- !is.na(sumstats_df[[chr_col]]) & !is.na(sumstats_df[[bp_col]])
    sumstats_df <- sumstats_df[filter, ]

    # replace any non-ATCG code with N
    sumstats_df[[a1_col]] <- gsub('[^ATCG]', 'N', sumstats_df[[a1_col]])
    sumstats_df[[a2_col]] <- gsub('[^ATCG]', 'N', sumstats_df[[a2_col]])

    # trim variants
    trimmed <- trim_variants(sumstats_df, a1_col = a1_col, a2_col = a2_col, bp_col = bp_col)
    sumstats_df[[a1_col]] <- trimmed$a1
    sumstats_df[[a2_col]] <- trimmed$a2
    sumstats_df[[bp_col]] <- trimmed$bp

    # fill ID with chr:pos:ref:alt if rsID missing
    fallback_ids <- paste(sumstats_df[[chr_col]], sumstats_df[[bp_col]], sumstats_df[[a1_col]], sumstats_df[[a2_col]], sep=":")
    sumstats_df[is.na(sumstats_df[[snp_col]]), snp_col] <- fallback_ids[is.na(sumstats_df[[snp_col]])]

    # remove if missing beta
    filter <- !is.na(sumstats_df[[beta_col]])
    sumstats_df <- sumstats_df[filter, ]

    # round sample sizes
    if (n_col %in% colnames(sumstats_df)) sumstats_df[[n_col]] <- round(sumstats_df[[n_col]])
    if (n_cas_col %in% colnames(sumstats_df)) sumstats_df[[n_cas_col]] <- round(sumstats_df[[n_cas_col]])

    return(sumstats_df)
}


trim_variants <- function(sumstats_df,
                          a1_col = "A1",
                          a2_col = "A2",
                          bp_col = "BP") {

    df <- data.frame(a1 = sumstats_df[[a1_col]],
                     a2 = sumstats_df[[a2_col]],
                     nca1 = nchar(sumstats_df[[a1_col]]),
                     nca2 = nchar(sumstats_df[[a2_col]]),
                     bp = sumstats_df[[bp_col]])
    
    # right trim
    trimmable <- df$nca1 > 1 &
                 df$nca2 > 1 & 
                 substr(df$a1, df$nca1, df$nca1) == substr(df$a2, df$nca2, df$nca2)
    trimmable[is.na(trimmable)] <- FALSE

    while (any(trimmable, na.rm=TRUE)){
        df[trimmable, 'a1'] <- substr(df[trimmable, 'a1'], 1, df[trimmable, 'nca1'] - 1)
        df[trimmable, 'a2'] <- substr(df[trimmable, 'a2'], 1, df[trimmable, 'nca2'] - 1)
        df[trimmable, 'nca1'] <- nchar(df[trimmable, 'a1'])
        df[trimmable, 'nca2'] <- nchar(df[trimmable, 'a2'])

        trimmable <- df$nca1 > 1 &
                     df$nca2 > 1 & 
                     substr(df$a1, df$nca1, df$nca1) == substr(df$a2, df$nca2, df$nca2)
        trimmable[is.na(trimmable)] <- FALSE
    }

    # left trim
    trimmable <- df$nca1 > 1 &
                 df$nca2 > 1 & 
                 substr(df$a1, 1, 1) == substr(df$a2, 1, 1)
    trimmable[is.na(trimmable)] <- FALSE

    while (any(trimmable, na.rm=TRUE)){
        df[trimmable, 'a1'] <- substr(df[trimmable, 'a1'], 2, df[trimmable, 'nca1'])
        df[trimmable, 'a2'] <- substr(df[trimmable, 'a2'], 2, df[trimmable, 'nca2'])
        df[trimmable, 'nca1'] <- nchar(df[trimmable, 'a1'])
        df[trimmable, 'nca2'] <- nchar(df[trimmable, 'a2'])
        df[trimmable, 'bp'] <- df[trimmable, 'bp'] + 1

        trimmable <- df$nca1 > 1 &
                     df$nca2 > 1 & 
                     substr(df$a1, 1, 1) == substr(df$a2, 1, 1) 
        trimmable[is.na(trimmable)] <- FALSE
    }

    return(df)
}

# impute missing P vals from beta and SE
impute_p <- function(df, beta_col='BETA', se_col='SE', p_col='P') {
  BETA = df[[beta_col]]
  SE = df[[se_col]]
  P = df[[p_col]]
  
  new_P <- if_else(!is.na(BETA) & !is.na(SE) & is.na(P), 
                   pchisq((BETA/SE)^2,  df = 1, lower.tail = F),
                   P)
  return(new_P)
}


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
