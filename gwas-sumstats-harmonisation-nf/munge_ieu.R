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

harmonise_ieu_metadata <- function(input_df) {
    out_df <- data.frame(
        TotalSamples=as.numeric(input_df[["sample_size"]]),
        TotalCases=as.numeric(input_df[["ncase"]]),
        TotalControls=as.numeric(input_df[["ncontrol"]]),
        Population=input_df[["population"]],
        StudySample=input_df[["sex"]],
        TraitReported=input_df[["trait"]],
        TraitEfos=input_df[["ontology"]],
        TraitTags=input_df[["subcategory"]],
        PubAuthor=paste0(input_df[["author"]], "; ", input_df[["consortium"]]),
        PubDate=input_df[["year"]],
        Comments=input_df[["note"]])
        
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
