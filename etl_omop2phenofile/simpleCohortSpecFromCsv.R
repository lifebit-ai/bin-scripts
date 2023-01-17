#!/usr/bin/env Rscript

################################################################################################################
# This script generates simple cohort_specs from a codelist.
################################################################################################################

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
	library(Capr)
})

options(scipen = 99999)

# Collecting arguments
args <- commandArgs(TRUE)
args

# Default setting when not all arguments passed or help is needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) ) {

	cat("
	This script generates OHDSI JSON cohort definitions based on the simple cohort_specs provided.
    
	Mandatory arguments:
    --codelist
    --connection_details
	
	Optional arguments:
    --help          This help message.

    Usage:
    simpleCohortSpecFromCsv.R --codelist=P000001.csv
    \n")

	q(save="no")
}

# Parsing arguments
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

codelist <- args$codelist
connectionDetailsFull <- jsonlite::read_json(args$connection_details)
vocabularyDatabaseSchema <- connectionDetailsFull$cdmDatabaseSchema
domain <- args$domain
source_or_standard <- str_c(snakecase::to_snake_case(args$concept_type),"s")

## Database Connection Jars
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER = getwd())
unzip(args$db_jars)

## Database Connection
connectionDetails <- discard(connectionDetailsFull, names(connectionDetailsFull) %in% c('cohortDatabaseSchema', 'cdmDatabaseSchema'))
connectionDetails <- exec(DatabaseConnector::createConnectionDetails, !!! connectionDetails)
connection <- connect(connectionDetails)

## Read input file and find codes in OMOP
input_file <- read_csv(codelist, col_types = cols(.default = "c")) %>%
  select(Phenotype_Short, Criteria_Ontology, Criteria, Is_Inclusion_Criteria, Is_Exclusion_Criteria) %>%
  mutate(across(c(Is_Inclusion_Criteria, Is_Exclusion_Criteria), as.logical))

input_file_split <- input_file %>%
  split(.$Criteria_Ontology)

input_file_split$ICD10CM <- input_file_split$ICD10

input_file_full <- map2_df(input_file_split, names(input_file_split), ~ mutate(getConceptCodeDetails(conceptCode = .x$Criteria, vocabulary = .y, connection = connection, vocabularyDatabaseSchema = vocabularyDatabaseSchema, mapToStandard = FALSE), Criteria_Ontology = unique(.x$Criteria_Ontology))) %>%
  select(Criteria_Ontology, Criteria = conceptCode, conceptId) %>%
  mutate(across(everything(), as.character)) %>%
  left_join(input_file) %>%
  split({.}$Phenotype_Short) 

if(length(input_file_full) > 1) stop("Only one phenotype per run support in this version")

## Build simple cohort definition based on codeset
cohort_definition_codes <- input_file_full[[1]]

inclusions <- as.integer(cohort_definition_codes$conceptId[cohort_definition_codes$Is_Inclusion_Criteria])
exclusions <- as.integer(cohort_definition_codes$conceptId[cohort_definition_codes$Is_Exclusion_Criteria])

cohort_definition <- list(
    list(
        cohort_name = names(input_file_full), 
        cohort_identifier = 2, 
        index_event = list(domain = domain, occurrence = "first"),
        control_group = list(cohort_identifier = 1, domain = domain, occurrence = args$control_group_occurrence)
)) 
  
cohort_definition[[1]]$index_event[[source_or_standard]] <- inclusions

if(length(exclusions) > 0){
  cohort_definition[[1]]$qualifying_events <-  list(list(domain = domain, time_period = "all", count = list(logic = "exactly", amount = 0)))
  cohort_definition[[1]]$qualifying_events[[1]][[source_or_standard]] <- exclusions
}

jsonlite::write_json(cohort_definition, str_c("user_",names(input_file_full),".json"), auto_unbox = T, pretty = T)
