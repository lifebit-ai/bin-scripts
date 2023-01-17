#!/usr/bin/env Rscript

################################################################################################################
# This script generates OHDSI JSON cohort definitions based on the simple cohort_specs provided.
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
    --cohort_specs     	  A json containing one or more cohort specifications
	  --connection_details  A json containing connection details
	  --db_jars             A zipped file containg database connection jars

	Optional arguments:
    --help          This help message.

    Usage:
    createCohortJsonFromSpec.R --cohort_specs=cohorts_spec.json --connection_details=connection_details.json --db_jars=pgjars.zip
    \n")

	q(save="no")
}

#' Creates an OHDSI compliant cohort JSON from a user specification
#'
#' @param cohort_definition The user cohort definition, in list format read from JSON
#' @param connection Database Connection
#' @param vocabularyDatabaseSchema The schema where the concept table can be found in the database
#'
#' @return Writes a cohort JSON file
#' @export
#'
writeCohortJson <- function(cohort_definition, connection, vocabularyDatabaseSchema){
  
  if(cohort_definition$index_event$occurrence == "first") cohort_definition$index_event$occurrence <- "fist" # bug in Capr package

  index_concept_ids <- if(is.null(cohort_definition$index_event$concept_ids)) 0 else cohort_definition$index_event$concept_ids
  index_conceptset_name <-  if(is.null(cohort_definition$index_event$concept_ids)) "Null Codeset" else "Index Event"
    
  index_components <- list()
  attribute_list <- list()
  
  if(!identical(NULL, cohort_definition$index_event$source_concept_ids)){
    source_index_codeset <- getConceptIdDetails(conceptIds = cohort_definition$index_event$source_concept_ids, connection = connection, vocabularyDatabaseSchema = vocabularyDatabaseSchema, mapToStandard = F)  %>%
      createConceptSetExpression(Name = "Source Index Event", includeDescendants = F)
    # Bug in package. Need to add source concept attribute as an index component otherwise the concept set will not register
    attribute_list[[1]] <- createSourceConceptAttribute(cohort_definition$index_event$domain, source_index_codeset)
    index_components[[1]] <-  Capr:::createQuery(Component = source_index_codeset, Domain = cohort_definition$index_event$domain)
  }  
  
  ### Compile indexing event
  index_codeset <- getConceptIdDetails(conceptIds = index_concept_ids, connection = connection, vocabularyDatabaseSchema = vocabularyDatabaseSchema, mapToStandard = F)  %>%
    createConceptSetExpression(Name = index_conceptset_name, includeDescendants = F)

  index_components <-  append(index_components, Capr:::createQuery(Component = index_codeset, Domain = cohort_definition$index_event$domain, attributeList = attribute_list))

  primary <- createPrimaryCriteria(Name = "Primary Criteria", ComponentList = index_components, ObservationWindow = createObservationWindow(0L,0L), Limit =  cohort_definition$index_event$occurrence)
    
  ### Add all inclusion rules. Map because there may be multiple
  inclusion_rules <- map2(cohort_definition$qualifying_events, c(1:length(cohort_definition$qualifying_events)), function(qualifying_event, event_number){
    
    inclusion_concept_ids <- if(is.null(qualifying_event$concept_ids)) 0 else qualifying_event$concept_ids
    inclusion_conceptset_name <-  if(is.null(qualifying_event$concept_ids)) "Null Codeset"  else str_c("Inclusion", event_number)
        
    inclusion_attribute_list <- list()
    inclusion_components <- list()
    
    if(!identical(NULL, qualifying_event$source_concept_ids)){
      source_inclusion_codeset <- getConceptIdDetails(conceptIds = qualifying_event$source_concept_ids, connection = connection, vocabularyDatabaseSchema = vocabularyDatabaseSchema, mapToStandard = F)  %>%
        createConceptSetExpression(Name = str_c("Source Inclusion Event", event_number), includeDescendants = F)
      inclusion_attribute_list[[1]] <- createSourceConceptAttribute(qualifying_event$domain, source_inclusion_codeset)
    }	  

    inclusion_codeset <- getConceptIdDetails(conceptIds = inclusion_concept_ids, connection = connection, vocabularyDatabaseSchema = vocabularyDatabaseSchema, mapToStandard = F)  %>%
      createConceptSetExpression(Name =  inclusion_conceptset_name, includeDescendants = F)
    
    inclusion_components <-  Capr:::createQuery(Component = inclusion_codeset, Domain = qualifying_event$domain,  attributeList = inclusion_attribute_list)

    timeline <- createTimeline(StartWindow = createWindow(StartDays = "All", StartCoeff = "Before", EndDays = "All", EndCoeff = "After"))
    InclusionRule1 <- createCount(Query = inclusion_components, Logic = qualifying_event$count$logic, Count = qualifying_event$count$amount, Timeline = timeline)
    createGroup(Name = str_c("Inclusion", event_number), type = "ALL", criteriaList = list(InclusionRule1), demographicCriteriaList = NULL, Groups = NULL)
  })
  
  Inculsion_rules <- createInclusionRules(Name = "Inclusion Rules", Contents = inclusion_rules, Limit = "First")
  
  ## Compile full cohort definition
  cd <- createCohortDefinition(Name = "My cohort", Description = "My cohoft", PrimaryCriteria = primary, InclusionRules = Inculsion_rules)
  cohort_json <- jsonlite::fromJSON(Capr::compileCohortDefinition(cd), simplifyVector = F)
  
  ## Fix workaround needed for bugs in package
  null_cs_ids <- unique(map_int(cohort_json$ConceptSets, function(cs) if(cs$name == "Null Codeset"){ return(as.integer(cs$id)) }else{ return(as.integer(-1)) }))
  removeCodeSetId <- function(cohortDefinition, null_cs_ids){
    if(class(cohortDefinition) != "list") return(cohortDefinition)
    if("ConditionOccurrenceSourceConcept" %in% names(cohortDefinition)) names(cohortDefinition)[names(cohortDefinition) == "ConditionOccurrenceSourceConcept"] <- "ConditionSourceConcept"
    if("CodesetId" %in% names(cohortDefinition))
      if(cohortDefinition[["CodesetId"]] %in% null_cs_ids) cohortDefinition[["CodesetId"]] <- NULL
    map(cohortDefinition, removeCodeSetId, null_cs_ids)
  }
    
  cohort_json <- removeCodeSetId(cohort_json, null_cs_ids)
  
  if(!identical(NULL, cohort_definition$index_event$source_concept_ids)) cohort_json$PrimaryCriteria$CriteriaList[[1]] <- NULL
  codesets_to_remove <- which(map_lgl(cohort_json$ConceptSets, ~.x$name == "Null Codeset"))
  
  if(length(codesets_to_remove) > 0) cohort_json$ConceptSets[codesets_to_remove] <- NULL
  
  ## Write cohort json
  jsonlite::write_json(cohort_json, stringr::str_c(cohort_definition$cohort_name, ".json"), auto_unbox = T, pretty = T)
  
}

## Parsing arguments
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

connectionDetailsFull <- jsonlite::read_json(args$connection_details)
cdmDatabaseSchema <- connectionDetailsFull$cdmDatabaseSchema

## Database Connection Jars
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER = getwd())
unzip(args$db_jars)

## Database Connection
connectionDetails <- discard(connectionDetailsFull, names(connectionDetailsFull) %in% c('cohortDatabaseSchema', 'cdmDatabaseSchema'))
connectionDetails <- exec(DatabaseConnector::createConnectionDetails, !!! connectionDetails)
connection <- connect(connectionDetails)

## Read cohort specs and generate OHDSI JSON
cohort_definitions_spec_input <- jsonlite::read_json(args$cohort_specs, simplifyVector = F)
walk(cohort_definitions_spec_input, ~writeCohortJson(.x, connection, vocabularyDatabaseSchema = cdmDatabaseSchema))
