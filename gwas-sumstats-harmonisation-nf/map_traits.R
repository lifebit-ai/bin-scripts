#!/usr/bin/env Rscript

# This script takes a list of trait names or codes in csv/tsv format and for each item
# searches for a match in the OMOP concept table (field to search is defined by
# --search_field). If a match is found, the term is mapped to the SNOMED vocabulary
# and a disease based grouping term is provided based on the desired grouping level.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(dtplyr)
})

snake_to_camel <- function(x){
  gsub("(?:\\b|_)([[:alpha:]])", "\\U\\1", x, perl=TRUE)
}

### Remove scientific notation
options(scipen = 99999)

# Collecting arguments
args <- commandArgs(TRUE)
args

# Default setting when not all arguments passed or help is needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) ) {
  
  cat("

    Mandatory arguments:
	  --traits_csv: the location of the traits csv list
    --vocabulary_folder: the path to OMOP vocabularies files.

    Optional arguments:
    --search_field: the field in the OMOP concept table to look for the supplied trait. Defaults to concept_name but could also look for codes in concept_code
    --vocabulary_id: a vocabulary (e.g. ICD10CM) to limit the search.
    --search_priority: location to a csv with vocabularies ordered by priorty in cases where the trait is matched in multiple vocabularies. Defaults to: SNOMED, ICD10CM, UK Biobank, ICD9CM, Read, OXMIS, Mesh
    --snomed_grouping_level: The level at which to group output SNOMED concepts. Defaults to 4.

	Optional arguments:
    --help          This help message.

    Usage:
    ./mapTraits.R --traits_csv=testdata/traits.csv --vocabulary_folder=omop-vocab-files

    \n")
  
  q(save="no")
}



# Parsing arguments
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")
argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

##################################
# Variables reassignment section #
##################################
traits_csv <- args$traits_csv
vocabulary_folder <- args$vocabulary_folder
search_field <- args$search_field
vocabulary_id <- args$vocabulary_id
desired_snomed_group <- as.integer(args$snomed_grouping_level)

message("Vocab folder")
message(vocabulary_folder)

search_terms <- read_csv(traits_csv)[[1]]
traits_input <- data.table(trait_reported = search_terms, vocabulary_id = vocabulary_id)

if(identical(desired_snomed_group, as.integer(NULL))) desired_snomed_group <- 4
if(identical(search_field, NULL)) search_field <- "concept_name"
if(identical(args$search_priority, NULL)) search_priority <- c("SNOMED", "ICD10CM", "UK Biobank", "ICD9CM", "Read", "OXMIS", "Mesh")
if(!identical(args$search_priority, NULL)) search_priority <- read_csv(search_priority)[[1]]


concept <- fread(file.path(vocabulary_folder, "CONCEPT.csv"),na.strings = c("NA", ""))
cr <- fread(file.path(vocabulary_folder, "CONCEPT_RELATIONSHIP.csv"),na.strings = c("NA", ""))
ca <- fread(file.path(vocabulary_folder, "CONCEPT_ANCESTOR.csv"),na.strings = c("NA", ""))
categories <- filter(ca, ancestor_concept_id %in% c(441840), min_levels_of_separation >= desired_snomed_group) %>%
  select(ancestor_concept_id = descendant_concept_id, distance_from_top = min_levels_of_separation)

omop_vocabulary_fields <- c("concept_id", "concept_name", "concept_code", "vocabulary_id", "domain_id", "concept_class_id", "standard_concept")
required_output_cols <- c("trait_reported", "trait_source_origin",  str_c("trait_source_", omop_vocabulary_fields), str_c("trait_", omop_vocabulary_fields), "trait_group_concept_name", "trait_group_concept_id")

maps_from <- mutate(traits_input, join_col = tolower(trait_reported)) %>%
  left_join(mutate(concept, join_col = tolower(!!sym(search_field)))) %>%
  mutate(priority = coalesce(match(vocabulary_id, search_priority), as.integer(999999999))) %>%
  group_by(trait_reported) %>%
  arrange(priority) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.table()

maps_to <- select(maps_from, trait_reported, concept_id_1 = concept_id, standard_concept) %>%
  left_join(filter(cr, relationship_id == "Maps to")) %>%
  select(trait_reported, concept_id = concept_id_2) %>%
  inner_join(filter(concept, standard_concept == "S")) %>%
  as.data.table()

grouping_snomed <- select(maps_to, trait_reported, descendant_concept_id = concept_id) %>%
  filter(!is.na(descendant_concept_id)) %>%
  inner_join(ca) %>%
  inner_join(categories) %>%
  group_by(trait_reported, descendant_concept_id) %>%
  filter(distance_from_top == min(distance_from_top)) %>%
  ungroup %>%
  select(trait_reported, descendant_concept_id, concept_id = ancestor_concept_id) %>%
  left_join(select(concept, concept_name, concept_id)) %>%
  as.data.table()

grouping_other <- filter(maps_to, !concept_id %in% grouping_snomed$descendant_concept_id) %>%
  select(trait_reported, concept_name = domain_id) %>%
  mutate(concept_id = NA_integer_) %>%
  as.data.table()

grouping <- bind_rows(grouping_snomed, grouping_other) %>% 
  select(-descendant_concept_id)

traits_reported <- mutate(traits_input, trait_source_origin =  case_when(is.null(vocabulary_id) ~ "Inferred", TRUE ~ "Provided")) %>%
  as_tibble()

traits_output <- map(list(trait_source = maps_from, trait = maps_to, trait_group = grouping),
    function(x){
      x %>%
        select(trait_reported, any_of(omop_vocabulary_fields)) %>%
        mutate(across(everything(), as.character)) %>%
        as_tibble()
    }) %>%
  map2(., names(.), 
    function(x,y){
      rename_with(x, function(c){paste0(y, "_", c)}, -trait_reported)
    })

plyr::join_all(c(list(traits_reported), traits_output), "trait_reported") %>%
  .[required_output_cols] %>%
  rename_with(snake_to_camel) %>%
  write_csv("mapped.csv")
