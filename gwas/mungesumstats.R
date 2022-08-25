#!/usr/bin/env Rscript

##############################################################################
# Wrapper to Script to run MungeSumstats and fix quirks
#
# Usage: Rscript mungesumstats.R --gwas_table=<input> \
#          --outfile=<output_harmonised.vcf> \
#          --map_pos=<map_pos> \
#          --build=<genome_build> \
#          --drop=<drop_cols> \
#          --ncpus=<cpu_threads>
#   <input> [Required] Path to the GWAS summary statistics table.
#   <output_harmonised.vcf> [Optional] Path to the output harmonised GWAS
#       VCF. (Default: <input>.vcf )
#   <map_pos> [Optional] Logical value whether to perform SNP mapping to rsIDs using genomic
#        coordinates. (Default: FALSE)
#   <build> [Optional] Name of genome build (e.g. "GRCh38"). If not set, genome build
#        will be inferred from the data.
#   <drop_cols> [Optional] comma seperated list of columns to drop from sumstat file
#   <cpu_threads> [Optional] Number of CPU threads to use.
##############################################################################

suppressPackageStartupMessages({
  library(MungeSumstats)
  library(data.table)
})

## Parse arguments (we expect the form --arg=value)
args <- commandArgs(TRUE)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

if(is.null(args$outfile)) {args$outfile <- gsub('\\.[^\\.]+$', "", args$gwas_table)} # default is input file basename
if(is.null(args$ncpu)) {args$ncpu <- 1} else {args$ncpu <- as.numeric(args$ncpu)}
if(is.null(args$map_pos)) {do_map_pos <- FALSE} else {do_map_pos <- as.logical(args$map_pos)}

if (do_map_pos && is.null(args$build)) stop("A genome build must be provided with --build when using --map_pos")

# read in the sumstats table and drop cols if requested
drop <- if(is.null(args$drop)) NULL else strsplit(args$drop,",")[[1]]
gwas_df <- data.table::fread(args$gwas_table, drop=drop)

# do the munge
if (is.null(args$build)){
  gbuild <- MungeSumstats::get_genome_builds(sumstats_list=list(gwas_table=gwas_df))[[1]]
} else {
  gbuild <- args$build
}

if (do_map_pos){
  use_rsids <- FALSE
} else {
  use_rsids <- TRUE
}

# VCF conversion does not handle commas and colons in snp IDs well, these are instead substituted to underscores

if ("ID" %in% colnames(gwas_df)) {
    gwas_df$ID <- gsub(':', '_', gwas_df$ID)
    gwas_df$ID <- gsub(',', '_', gwas_df$ID)
} else if ("SNP" %in% colnames(gwas_df)) {
    gwas_df$SNP <- gsub(':', '_', gwas_df$SNP)
    gwas_df$SNP <- gsub(',', '_', gwas_df$SNP)
} else if ("SNPID" %in% colnames(gwas_df)) {
    gwas_df$SNPID <- gsub(':', '_', gwas_df$SNPID)
    gwas_df$SNPID <- gsub(',', '_', gwas_df$SNPID)
}

outfile_gz <- MungeSumstats::format_sumstats(path=gwas_df,
                                ref_genome=gbuild,
                                snp_ids_are_rs_ids=use_rsids,
                                save_path=args$outfile,
                                write_vcf=TRUE,
                                nThread=args$ncpu)

# fix mungesumstats filenaming
outfile <- gsub("\\.gz$", "", outfile_gz)
system(paste("mv", outfile_gz, outfile))

# add genome build to the vcf header
system(paste0("sed -i '4i ##genome_build=", gsub("GRCH", "GRCh", gbuild), "' ", outfile))
