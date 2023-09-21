#' Title
#' Performs the split read pairing using the samples of each age group
#' @param age.groups 
#' @param project.id 
#' @param project.name 
#' @param gtf.version 
#' @param age.samples.clusters 
#' @param results.folder 
#' @param replace 
#'
#' @return
#' @export
#'
#' @examples
age_stratification_junction_pairing <- function (age.groups,
                                                 project.id,
                                                 project.name,
                                                 gtf.version,
                                                 age.samples.clusters,
                                                 results.folder,
                                                 replace) {
  
  
  
  ## Loop per age cluster
  doParallel::registerDoParallel(3)
  
  foreach(i = seq(length(age.groups)), .combine = "rbind") %dopar%{
    
    age_cluster <- age.groups[i]
    
    # age_cluster <- age.groups[1]
    # age_cluster <- age.groups[2]
    
    print(paste0(Sys.time(), " - Starting analysis for '", project.id, "' and '", age_cluster, "' samples..."))
    
    #######################################################################
    ## GENERATE THE IDB
    #######################################################################
    
    
    ## Load data ----------------------------------------------------------------
    
    samples_age <- readRDS(file = paste0(results.folder, "/base_data/", 
                                         project.id, "_", age_cluster, "_samples_used.rds"))
    
    split_read_counts_age <- readRDS(file = paste0(results.folder, "/base_data/", 
                                                   project.id, "_", age_cluster, "_split_read_counts.rds")) 
    
    if ( is.null(names(split_read_counts_age)) ) {
      split_read_counts_age <- split_read_counts_age %>%
        as_tibble("junID")
    }
    
    all_split_reads_age <- readRDS(file = paste0(results.folder, "/base_data/",
                                                 project.id, "_", age_cluster, "_all_split_reads.rds")) %>% 
      as_tibble()
    
    if (!identical(split_read_counts_age$junID,
                   all_split_reads_age$junID)) {
      paste0(project.id, " - some jxn have not been considered!") %>% print()
    }
    
    ## Call functions ----------------------------------------------------------------
    
    folder_results <- paste0(results.folder, "/junction_pairing/", age_cluster, "/")
    dir.create(file.path(folder_results), showWarnings = T, recursive = T)
    
    get_distances(cluster = age_cluster,
                  samples = samples_age,
                  split_read_counts = split_read_counts_age,
                  all_split_reads_details = all_split_reads_age,
                  folder_name = folder_results,
                  replace = replace)
    
    
    extract_distances(cluster = age_cluster,
                      samples = samples_age,
                      folder_name = folder_results,
                      replace = replace)
    
    
    get_never_misspliced(cluster = age_cluster,
                         samples = samples_age,
                         split_read_counts = split_read_counts_age,
                         all_split_reads_details = all_split_reads_age,
                         folder_name = folder_results,
                         replace = replace)
    
  }
  
  rm(split_read_counts_age)
  rm(all_split_reads_age)
  rm(samples_age)
  gc()
}