#' Title
#' Each recount3 project has a different data source, hence the metadata is structured differently
#' This function extracts the metadata (clusters, RIN, etc) from each recount3 project depending on its data source
#' @param project.metadata raw metadata object as provided by recount3 (https://rdrr.io/bioc/recount3/man/read_metadata.html)
#' @param data.source source of the project in recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' an standardised metadata object 
#' @export
#'
#' @examples
separate_clusters <- function(project.metadata,
                              data.source) {
  
  
  if ( data.source == "data_sources/gtex" ) {
    
    project_metadata_tidy <- project.metadata %>%
      as_tibble()  %>%
      mutate(gtex.smrin = gtex.smrin %>% as.double()) %>%
      filter(gtex.smafrze != "EXCLUDE",
             #!(gtex.smtsd %in% c("Brain - Cortex", "Brain - Cerebellum")),
             gtex.smrin >= 6.0) ## Only fresh-frozen preserved tissues
    
    
    stopifnot(
      "Still there are samples with RIN < 6 and labelled as EXCLUDE" =
        !any(project_metadata_tidy$gtex.smrin < 6)
    )
    stopifnot(
      "Still there are samples labelled as EXCLUDE" =
        !any(project_metadata_tidy$gtex.smafrze == "EXCLUDE")
    )
    
    
    ## TODO same standard column names as with data_sources/sra projects
    project_metadata_tidy <- data.frame(external_id =  project_metadata_tidy$external_id %>% as.character(),
                                        sample_id =  project_metadata_tidy$gtex.sampid %>% as.character(),
                                        age = project_metadata_tidy$gtex.age %>% as.character(),
                                        rin = project_metadata_tidy$gtex.smrin %>% as.double(),
                                        gender = project_metadata_tidy$gtex.sex %>% as.character(),
                                        cluster = project_metadata_tidy$gtex.smtsd,
                                        smnabtcht = project_metadata_tidy$gtex.smnabtcht,
                                        smafrze = project_metadata_tidy$gtex.smafrze,
                                        avg_read_length = project_metadata_tidy$recount_seq_qc.avg_len,
                                        all_mapped_reads = project_metadata_tidy$recount_qc.star.all_mapped_reads,
                                        SRA_project = project_metadata_tidy$recount_project.project) 
    
    
    
    
  } else if (data.source == "data_sources/sra") {
    
    
    
    
    project_metadata_tidy <- project.metadata %>%
      as_tibble() %>%
      dplyr::select(external_id, 
                    sra.experiment_title, 
                    sra.sample_attributes, 
                    all_mapped_reads = recount_qc.star.all_mapped_reads) %>%
      mutate(rn = row_number()) %>%
      separate_rows(sra.sample_attributes, sep = "\\|\\s*") %>%
      separate(sra.sample_attributes, into = c('col1', 'col2'), sep = ";;") %>% 
      pivot_wider(names_from = col1, values_from = col2) %>% 
      dplyr::select(-rn) %>%
      as.data.frame() %>%
      mutate(cluster = ifelse(test = str_detect(sra.experiment_title, pattern="AD"),
                              yes = "AD",
                              no = ifelse( test = str_detect(sra.experiment_title, pattern="P"),
                                           yes = "PD",
                                           no = "control")),
             cluster = cluster %>%as.factor(),
             #rin_score = rin_score %>% as.double(),
             avg_read_length = project.metadata$recount_seq_qc.avg_len,
             SRA_project = project.metadata$recount_project.project) 
    
    project_metadata_tidy <- project_metadata_tidy %>%
      mutate(rin_score = if (exists('rin_score', where = project_metadata_tidy)) rin_score else NA) %>%
      mutate_at(vars(one_of('rin_score')), as.double) %>%
      dplyr::rename(age = "age at death",
                    rin = "rin_score",
                    gender = if (exists('Sex', where = project_metadata_tidy)) "Sex" else "gender") %>%
      mutate(age = age %>% as.integer()) %>%
      as_tibble()
    
    
    if ( !all(is.na(project_metadata_tidy$rin)) ) {
      project_metadata_tidy <- project_metadata_tidy %>%
        filter(rin >= 6.0) ## Only fresh-frozen preserved tissues
    }
  
    project_metadata_tidy
  
    } else {
    ## TODO "data_sources/tcga"
  }
  
  
  return(project_metadata_tidy)
  
}
