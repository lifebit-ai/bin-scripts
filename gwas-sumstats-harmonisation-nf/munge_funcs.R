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
