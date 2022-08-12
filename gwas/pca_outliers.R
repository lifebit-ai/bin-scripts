#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)


# required
input_file <- as.character(args$input)
output_file <- as.character(args$out)
# optional
sigma <- ifelse(is.null(args$sigma), 6, as.numeric(args$sigma))
##############################################################################

## Run the outlier removal process
df <- read.table(input_file, header = T, sep="\t", comment.char = "")
all_outlier_idx = rep(FALSE,nrow(df))

for (n in 3:ncol(df)) {
  p <- df[,n]
  
  outlier_idx <- abs(p - mean(p)) > sd(p) * sigma
  
  print(paste(colnames(df)[n], ":", sum(outlier_idx), "outliers"))
  
  all_outlier_idx <- all_outlier_idx | outlier_idx
}

print(paste("Total :", sum(all_outlier_idx ), "outliers"))

plink_remove <- df[all_outlier_idx,1:2]
write.table(plink_remove, output_file, sep="\t", quote=F, row.names=F, col.names=F)


