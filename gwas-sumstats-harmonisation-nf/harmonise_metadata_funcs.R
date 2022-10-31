harmonise_ebi_metadata <- function(input_df) {
    out_df <- data.frame(
        TotalSamples=as.numeric(input_df[["sample_size"]]),
        TraitReported=input_df[["MAPPED_TRAIT"]],
        TraitEfos=input_df[["MAPPED_TRAIT_URI"]],
        PubTitle=input_df[["STUDY"]],
        PubAuthor=input_df[["FIRST AUTHOR"]],
        PubDate=input_df[["DATE"]],
        PubJournal=input_df[["JOURNAL"]])
    
    return(out_df)
}

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

harmonise_gwas_table_metadata <- function(input_df) {
    if (!('GeneticModel' %in% colnames(input_df))) {
        input_df$GeneticModel <- NA
    }
    if (!('CovarNames' %in% colnames(input_df))) {
        input_df$CovarNames <- NA
    }
    if (!('PhenoName' %in% colnames(input_df))) {
        input_df$PhenoName <- NA
    }
    out_df <- data.frame(
        TotalCases=as.numeric(input_df[["TotalCases"]]),
        TotalControls=as.numeric(input_df[["TotalControls"]]),
        Population=input_df[["Population"]],
        StudyType=input_df[["StudyType"]],
        Method=input_df[["Method"]],
        TotalSamples=input_df[["TotalSamples"]],
        GeneticModel=input_df[["GeneticModel"]],
        Covariates=input_df[["CovarNames"]],
        TraitReported=input_df[["PhenoName"]])

    return(out_df)
}
