#!/usr/bin/env Rscript

# This script adds trait metadata from OMOP trait mapping to the study metadata
# fields in the input GWAS VCF header. 'trait_table' should have a column for
# 'ReportedTrait' which will be used as the JOIN field for the study metadata
# field 'ReportedTrait' in the GWAS VCF header.

suppressPackageStartupMessages({
library(optparse)
library(VariantAnnotation)
library(dplyr)
})

option_list = list(
    make_option(c("-g", "--gwas_vcf_head"), type="character", default=NULL,
        help="", metavar="character"),
    make_option(c("-m", "--trait_table"), type="character", default=NULL,
        help="", metavar="character"),
    make_option(c("-f", "--field_descriptions"), type="character", default=NULL,
        help="", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

traits <- read.csv(opt$trait_table, colClasses = "character") %>% distinct(TraitReported, .keep_all = TRUE)
field_descriptions <- read.csv(opt$field_descriptions, sep='\t')
row.names(field_descriptions) <- field_descriptions$ID

vcf <- readVcf('header.txt')
vcf_header <- header(vcf)
vcf_mdata <- meta(vcf_header)

orig_study_metadata <- vcf_mdata$SAMPLE
orig_study_metadata[] <- lapply(orig_study_metadata, function(x) gsub('^"|"$', '', x)) # remove surrounding quotes

study_metadata <- merge(orig_study_metadata, traits, by='TraitReported', all.x=TRUE)
rownames(study_metadata) <- rownames(orig_study_metadata)

# Add background traits if present
if (sum(is.na(study_metadata$TraitBackgroundReported)) > 0) {
    bg_traits <- traits
    colnames(bg_traits) <- sub('Trait', 'TraitBackground', colnames(bg_traits))

    study_metadata <- merge(study_metadata, traits, by='TraitBackgroundReported', all.x=TRUE)
    rownames(study_metadata) <- rownames(orig_study_metadata)
}


study_metadata <- study_metadata %>%
    as('data.frame') %>%
    mutate(across(any_of(field_descriptions[field_descriptions$Type=="Integer", "ID"]), as.integer)) %>%
    mutate(across(any_of(field_descriptions[field_descriptions$Type=="Float", "ID"]), as.double))

# add lines describing study metadata fields
vcf_mdata$META <- field_descriptions %>% 
    filter(key == "META" & ID %in% names(study_metadata) ) %>% 
    select(Number, Type, Description) %>%
    as("DataFrame")

# add SAMPLE lines containing study metadata for each study
vcf_mdata$SAMPLE <- study_metadata %>%  
    mutate(across(where(is.character), function(x){paste0('"',x,'"')}))

# set the header info and write VCF
meta(vcf_header) <- vcf_mdata
header(vcf) <- vcf_header
writeVcf(vcf, filename = 'updated_header.txt')
