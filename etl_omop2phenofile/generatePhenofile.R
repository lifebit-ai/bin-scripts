#!/usr/bin/env Rscript

####################################################################################################################################
# This script generates a phenofile for all cohorts in the cohort_table specified using the covariate_spec provided. 
####################################################################################################################################

######################
# Importing packages #
######################

suppressPackageStartupMessages({
  library(tidyverse)
  library(CohortGenerator)
  library(ROhdsiWebApi)
  library(FeatureExtraction)
  library(jsonlite)
  library(RJSONIO)
  library(rlang)
})

options(scipen = 99999)

# Collecting arguments
args <- commandArgs(TRUE)
args

# Default setting when not all arguments passed or help is needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) ) {
  
  cat("
	This script generates a phenofile for all cohorts in the cohort_table specified using the covariate_spec provided. 
    Mandatory arguments:
    --connection_details     	A json containing database connection details
    --cohort_counts           A csv containing individual cohort counts
    --cohort_table						A text file containing the database location of the cohorts
    --covariate_spec          A json containing covariate specifications
    --phenofile_name          The name for the output phenofile
	  Optional arguments:
    --help          This help message.
    Usage:
    ./generatePhenofile.R --connection_details=connectionDetails.json --cohort_counts=cohort_counts.csv --cohort_table=cohort_table_name.txt --covariate_spec=covariate_spec.json
    \n")
  
  q(save="no")
}

#' Title Create a phenofile from cohorts based on covariate_spec
#'
#' @param connectionDetails Connection details to connect to the database via DatabaseConnector
#' @param cohort_counts A table of cohort counts
#' @param cohort_table THe name of the cohort table
#' @param covariate_spec A covariate specification file
#' @param cdmDatabaseSchema The name of the schema with clinical data
#' @param cohortDatabaseSchema The name of the schema to insert the cohort table
#'
#' @return A phenofile for all cohorts in cohort_table
#' @export
#'
get_covariate_data <- function(connectionDetails, cohort_counts, cohort_table, covariate_spec, cdmDatabaseSchema, cohortDatabaseSchema){
  
  covariate_settings <- FeatureExtraction::createCovariateSettings(
    useDemographicsAge = T,
    useDemographicsGender = T,
    useDemographicsRace = T,
    useConditionOccurrenceLongTerm = F,
    useMeasurementValueLongTerm = F,
    longTermStartDays = -99999,
    endDays = 99999
  )
  
  covariate_data_split <- map(cohort_counts$cohortId, ~FeatureExtraction::getDbCovariateData(connectionDetails, cdmDatabaseSchema = cdmDatabaseSchema, cohortDatabaseSchema = cohortDatabaseSchema, cohortTable = cohort_table, cohortId = .x, covariateSettings = covariate_settings))
  
  names(covariate_data_split) <- cohort_counts$cohortId
  
  covariate_data <- list()
  covariate_data$covariateRef <- map_df(covariate_data_split, ~collect(.x$covariateRef)) %>% distinct
  covariate_data$covariates <- bind_rows(map(covariate_data_split, ~collect(.x$covariates)), .id = "cohortDefinitionId")
  
  connection <- DatabaseConnector::connect(connectionDetails)
  
  full_cohort <- DatabaseConnector::querySql(connection, glue::glue("SELECT CAST(COHORT_DEFINITION_ID as character) COHORT_DEFINITION_ID, SUBJECT_ID, PERSON_SOURCE_VALUE FROM {cohort_table} c LEFT JOIN {cdmDatabaseSchema}.person p on c.subject_id = p.person_id", cohort_table = cohort_table), snakeCaseToCamelCase = T)
  
  all_covs <- map(covariate_spec, function(cov){
    
    covariates_to_report <- filter(covariate_data$covariateRef, str_detect(covariateName, cov$covariate))
    
    if(!is.null(cov$concept_ids)) covariates_to_report <- filter(covariates_to_report, conceptId %in% cov$concept_ids)
    
    covariates <- filter(covariate_data$covariates , covariateId %in% !! covariates_to_report$covariateId) %>% collect

    if(!is.null(cov$valid_range)){
      if(!identical(cov$transformation, "value")) stop("valid_range only supported when covariate transformation is value")
      covariates <- filter(covariates, covariateValue >= cov$valid_range[1], covariateValue <= cov$valid_range[2])
    }
  
    if(is.null(cov$transformation)){
      tmp <- left_join(full_cohort, rename(covariates, subjectId = rowId)) %>%
        left_join(select(covariates_to_report, covariateId, covariateName)) %>%
        group_by(cohortDefinitionId, subjectId, personSourceValue) %>%
        summarise(!! cov$covariate_name := str_c(sort(unique(str_remove(covariateName[covariateValue==1], "^[^=]*= ")))))
    }
    
    if(!is.null(cov$transformation)){
      if(cov$transformation == "value"){
        tmp <- left_join(full_cohort, rename(covariates, subjectId = rowId)) %>%
          group_by(cohortDefinitionId, subjectId, personSourceValue) %>%
          summarise(!! cov$covariate_name := covariateValue) %>%
          return()
      }
      
      if(cov$transformation == "binary"){
        tmp <- left_join(full_cohort, rename(covariates, subjectId = rowId)) %>%
          group_by(cohortDefinitionId, subjectId, personSourceValue) %>%
          summarise(!! cov$covariate_name := coalesce(max(covariateValue), 0)) %>%
          return()
      }
      
      if(cov$transformation == "encoded"){
        tmp <- left_join(full_cohort, rename(covariates, subjectId = rowId)) %>%
          left_join(select(covariates_to_report, covariateId, covariateName)) %>%
          select(-covariateId) %>%
          mutate(covariateName = str_c(cov$covariate_name, ": ", coalesce(covariateName, "Unknown"))) %>%
          pivot_wider(names_from = covariateName, values_from = covariateValue, values_fill = 0)
      }
      
      if(cov$transformation == "categorical"){
        tmp <- left_join(full_cohort, rename(covariates, subjectId = rowId)) %>%
          left_join(select(covariates_to_report, covariateId, covariateName)) %>%
          group_by(cohortDefinitionId, subjectId, personSourceValue) %>%
          summarise(!! cov$covariate_name := str_c(sort(unique(str_remove(covariateName[covariateValue==1], "^[^=]*= ")))))
      }
      
      if(str_detect(cov$transformation, "^function")){
        
        transformation_fn <- eval(parse_expr(cov$transformation))
        
        tmp <- left_join(full_cohort, rename(covariates, subjectId = rowId)) %>%
          group_by(cohortDefinitionId, subjectId, personSourceValue) %>%
          summarise(!! cov$covariate_name := transformation_fn(covariateValue)) %>%
          return()
      }
    }

    if(!is.null(cov$imputation)){   
      if(!identical(cov$transformation, "value")) stop("imputation supported when covariate transformation is value")
      if(cov$imputation == "median") tmp <- group_by(tmp, cohortDefinitionId) %>% mutate(!! cov$covariate_name := coalesce(!! sym(cov$covariate_name), median(!! sym(cov$covariate_name), na.rm = T)))
      if(cov$imputation == "mean") tmp <- group_by(tmp, cohortDefinitionId) %>% mutate(!! cov$covariate_name := coalesce(!! sym(cov$covariate_name), mean(!! sym(cov$covariate_name), na.rm = T)))
    }
    
    return(ungroup(tmp))
    
  })
  
  reduce(all_covs, left_join)
  
}

# Parsing arguments
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)


connectionDetailsFull <- jsonlite::read_json(args$connection_details)
cohortCounts <- read_csv(args$cohort_counts)
cohortTable <- readLines(args$cohort_table)
covariateSpec <- jsonlite::fromJSON(args$covariate_spec, simplifyVector = F)
pheno_label <- args$pheno_label
convert_plink <- as.logical(args$convert_plink)
phenofile_name <- args$phenofile_name

if(nrow(cohortCounts) == 0) writeLines("Cohorts contain no patients", "empty_phenofile.tsv")

if(nrow(cohortCounts) > 0){
## Database Connection Jars
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER = getwd())
unzip(args$db_jars)

## Database Connection
connectionDetails <- discard(connectionDetailsFull, names(connectionDetailsFull) %in% c('cohortDatabaseSchema', 'cdmDatabaseSchema'))
connectionDetails <- exec(DatabaseConnector::createConnectionDetails, !!! connectionDetails)

phenofile <- get_covariate_data(connectionDetails, cohortCounts, cohortTable, covariateSpec, cdmDatabaseSchema = connectionDetailsFull$cdmDatabaseSchema, cohortDatabaseSchema = connectionDetailsFull$cohortDatabaseSchema)

## Tidy up tables
connection <- DatabaseConnector::connect(connectionDetails)
walk(str_c("DROP TABLE IF EXISTS ", str_c(cohortTable, c("","_inclusion_result", "_inclusion_stats", "_summary_stats", "_censor_stats"))),
	~DatabaseConnector::dbExecute(connection, .x))
DatabaseConnector::disconnect(connection)

if(convert_plink){  
  phenofile <- select(phenofile, `#FID` = personSourceValue, IID = personSourceValue, everything(), !! pheno_label := cohortDefinitionId, -subjectId) %>%
    mutate(across(matches("SEX|GENDER"), ~case_when(.x == "MALE" ~ "1", .x == "FEMALE" ~ "2", TRUE ~ NA_character_)))

  mutate(phenofile, group = if_else(!! sym(pheno_label) == 2, "case", "control")) %>%
    mutate(index=0) %>% 
    split(.$group) %>%
    list2env(envir = .GlobalEnv)

  set.seed(12345)

  for(i in seq_len(nrow(case))){
    x <- which(between(control$AGE, case$AGE[i] -2, case$AGE[i] +2) & 
               between(control$SEX, as.integer(case$SEX[i]) -0, as.integer(case$SEX[i]) + 0) & 
               control$index==0)
    control$index[sample(x, min(4, length(x)))] <- i
    case$index[i] <- i 
  }

phenofile <- rbind(case, control) %>% filter(index >0) %>%  arrange(index) %>% select(-any_of(c("group", "index")))

}

write_tsv(phenofile, str_c(phenofile_name, ".phe"))

}