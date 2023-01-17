#!/usr/bin/env Rscript

######################
# Importing packages #
######################

##############################################################################################################
# This script generates cohorts in the database using all json cohort definitions in the working directory. 
##############################################################################################################

suppressPackageStartupMessages({
	library(tidyverse)
	library(CohortGenerator)
	library(ROhdsiWebApi)
	library(FeatureExtraction)
	library(jsonlite)
	library(RJSONIO)
})

options(scipen = 99999)

# Collecting arguments
args <- commandArgs(TRUE)
args

# Default setting when not all arguments passed or help is needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) ) {

	cat("
	This script generates cohorts in the database using all json cohort definitions in the working directory. 
    
	Mandatory arguments:
    --connection_details     	A json containing database connection details
	
	Optional arguments:
    --help          This help message.
    
	Usage:
    generateCohorts.R --connection_details=connectionDetails.json
    \n")

	q(save="no")
}

#' Title Insert cohorts into the database
#'
#' @param connectionDetails Connection details to connect to the database via DatabaseConnector
#' @param cohortTableName The table name to insert cohorts into
#' @param cohortJsonFiles The OHDSI cohort definition json files to populate
#' @param cdmDatabaseSchema The name of the schema with clinical data
#' @param cohortDatabaseSchema The name of the schema to insert the cohort table
#'
#' @return
#' @export
#'
generateCohortsFromJson <- function(connectionDetails, cohortTableName = stringr::str_c("cohort_", stringi::stri_rand_strings(n=1, length = 10)), cohortJsonFiles, cdmDatabaseSchema, cohortDatabaseSchema = cdmDatabaseSchema){
	
	cohortsToCreate <- CohortGenerator::createEmptyCohortDefinitionSet()
	
	for (i in 1:length(cohortJsonFiles)) {
		cohortJsonFileName <- cohortJsonFiles[i]
		cohortName <- tools::file_path_sans_ext(basename(cohortJsonFileName))
		cohort_id <- keep(cohort_definitions_spec_input, ~.x$cohort_name == cohortName)[[1]]$cohort_identifier
		cohortJson <- readChar(cohortJsonFileName, file.info(cohortJsonFileName)$size)
		cohortExpression <- CirceR::cohortExpressionFromJson(cohortJson)
		cohortSql <- CirceR::buildCohortQuery(cohortExpression, options = CirceR::createGenerateOptions(generateStats = FALSE))
		cohortsToCreate <- rbind(cohortsToCreate, data.frame(cohortId = cohort_id, cohortName = cohortName,  sql = cohortSql, stringsAsFactors = FALSE))
	}
	
	cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable = cohortTableName)
	
	tables <- CohortGenerator::createCohortTables(connectionDetails = connectionDetails,
																								cohortDatabaseSchema = cohortDatabaseSchema,
																								cohortTableNames = cohortTableNames)
	
	cohortsGenerated <- CohortGenerator::generateCohortSet(
		connectionDetails = connectionDetails,
		cdmDatabaseSchema = cdmDatabaseSchema,
		cohortDatabaseSchema = cohortDatabaseSchema,
		cohortTableNames = cohortTableNames,
		cohortDefinitionSet = cohortsToCreate)
	
	cohortCounts <- CohortGenerator::getCohortCounts(
		connectionDetails = connectionDetails,
		cohortDatabaseSchema = cohortDatabaseSchema,
		cohortTable = cohortTableNames$cohortTable
		) %>% left_join(select(cohortsToCreate, cohortId, cohortName))
	
	return(list(cohort_table = cohortTableName, cohort_counts = cohortCounts))
	
}

# Parsing arguments
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

jsonFiles <- list.files(pattern = "\\.json")
cohortJsonFiles <- jsonFiles[jsonFiles != args$connection_details & jsonFiles != args$cohort_specs]
connectionDetailsFull <- jsonlite::read_json(args$connection_details)

## Database Connection Jars
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER = getwd())
unzip(args$db_jars)

## Database Connection Jars
connectionDetails <- discard(connectionDetailsFull, names(connectionDetailsFull) %in% c('cohortDatabaseSchema', 'cdmDatabaseSchema'))
connectionDetails <- exec(DatabaseConnector::createConnectionDetails, !!! connectionDetails)

cohort_definitions_spec_input <- jsonlite::read_json(args$cohort_specs, simplifyVector = F)

## Write cases
cohorts_generated <- generateCohortsFromJson(connectionDetails, cohortJsonFiles = cohortJsonFiles, cdmDatabaseSchema = connectionDetailsFull$cdmDatabaseSchema, cohortDatabaseSchema = connectionDetailsFull$cohortDatabaseSchema)

## Write controls
connection <- DatabaseConnector::connect(connectionDetails)

for(i in length(cohort_definitions_spec_input)){

	control_cohort_definition <- cohort_definitions_spec_input[[i]]$control_group

	if(identical(NULL, control_cohort_definition)) next

	domain <- snakecase::to_snake_case(control_cohort_definition$domain)
	case_cohort_id <- cohort_definitions_spec_input[[i]]$cohort_identifier
	control_cohort_id <- control_cohort_definition$cohort_identifier
	domain_date <- case_when(domain %in% c("condition_occurrence", "visit_occurrence") ~ str_c(str_remove(domain, "_occurrence"), "_start_date"), domain %in% c("observation", "procedure", "measurement") ~ str_c(domain, "_date"), TRUE ~ str_c(domain, "_start_date"))
	order <- case_when(control_cohort_definition$occurrence == "first" ~ domain_date, control_cohort_definition$occurrence == "last" ~ str_c(domain_date, " DESC"), TRUE ~ "RANDOM()")

	sql <- glue::glue("
		{set_seed}
		INSERT INTO {cohort_table}
		WITH CTE AS (
  			SELECT {control_cohort_id} as cohort_definition_id, d.person_id as subject_id, d.{domain_date} as cohort_start_date, d.{domain_date} as cohort_end_date, ROW_NUMBER() OVER (PARTITION BY d.person_id ORDER BY {order}) AS r
  			FROM {cdmDatabaseSchema}.{domain} d 
  			LEFT JOIN {cohortDatabaseSchema}.{cohort_table} c on d.person_id = c.subject_id and c.cohort_definition_id = {case_cohort_id}
  			WHERE c.subject_id is NULL
			ORDER BY {domain}_id
		)
		SELECT cohort_definition_id, subject_id, cohort_start_date, cohort_end_date FROM CTE where r = 1
		", cohort_table = cohorts_generated$cohort_table, control_cohort_id = control_cohort_id, domain = domain, case_cohort_id = case_cohort_id, cdmDatabaseSchema = connectionDetailsFull$cdmDatabaseSchema, cohortDatabaseSchema = connectionDetailsFull$cohortDatabaseSchema, order = order, domain_date = domain_date, set_seed = ifelse(connectionDetailsFull$dbms == "postgresql", "SELECT setseed(0.5);", ""))

	writeLines(sql, file("query.sql"))

	DatabaseConnector::dbExecute(connection, sql)

	## Count controls
	sql <- glue::glue("
	SELECT cohort_definition_id AS cohort_id, count(*) as cohort_entries, count(distinct subject_id) as cohort_subjects
	FROM {cohortDatabaseSchema}.{cohort_table} 
	WHERE cohort_definition_id = {control_cohort_id}
	GROUP BY cohort_definition_id
	", cohort_table = cohorts_generated$cohort_table, control_cohort_id = control_cohort_id, domain = domain, case_cohort_id = case_cohort_id, cohortDatabaseSchema = connectionDetailsFull$cohortDatabaseSchema)

	counts <- DatabaseConnector::querySql(connection, sql, snakeCaseToCamelCase = T) %>%
		mutate(cohortName = str_c(cohort_definitions_spec_input[[i]]$cohort_name, ": Control"))

	if(nrow(counts) > 0) cohorts_generated$cohort_counts <- bind_rows(cohorts_generated$cohort_counts, counts)
}

DatabaseConnector::disconnect(connection)

writeLines(cohorts_generated$cohort_table, con = "cohort_table_name.txt")
write_csv(cohorts_generated$cohort_counts, "cohort_counts.csv")
