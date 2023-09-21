#' Title
#' Per age group, it creates the split reads annotation files
#' @param age.groups 
#' @param project.name 
#' @param gtf.version 
#' @param project.id 
#' @param age.samples.clusters 
#' @param results.folder 
#' @param supportive.reads 
#'
#' @return
#' @export
#'
#' @examples
age_stratification_annotate <- function (age.groups,
                                         project.name,
                                         gtf.version,
                                         project.id,
                                         age.samples.clusters,
                                         results.folder,
                                         supportive.reads) {
  
  
  all_samples_used <- NULL
  local_results_folder <- file.path(results.folder, "/base_data/")
  
  ## Loop per age cluster
  doParallel::registerDoParallel(3)
  all_samples_age_used <- foreach(i = seq(length(age.groups))) %dopar%{
    
    age_cluster <- age.groups[i]
    message(age_cluster)
    
    # age_cluster <- age.groups[1]
    # age_cluster <- age.groups[2]
    
    ## THE AGE STRATIFICATION SHOULD BE DONE BY BODY SITE, AND STORED IN THE PATH OF THE BODY SITE.
    
    split_read_counts_age <- NULL
    all_split_reads_age <- NULL
    samples_age <- NULL
    j <- 1
    
    ## GENERATE THE COUNTS AND ANNOTATED SPLIT READ COUNTS BY AGE CLUSTER ---------------------------------------------------
    
    for ( cluster_id in age.samples.clusters$cluster %>% unique() ) {
      
      # cluster_id <- (age.samples.clusters$cluster %>% unique())[1]
      
      print(paste0(Sys.time(), " - ", cluster_id, ", ", age_cluster))
      
      ## Get the samples from the current tissue according to the current age cluster
      samples <- age.samples.clusters %>% 
        filter(age_group == age_cluster, cluster == cluster_id) %>% 
        pull(individual)
      
      
      if ( samples %>% length() > 0 ) {
        
        samples_age <- c(samples_age, samples)
        print(paste0(samples %>% length(), " samples... "))
        
        ## Load the annotated split reads (the splicing project contains all samples QC'ed and checked with recount3 origin)
        all_split_reads_local <- readRDS(file = paste0(local_results_folder, "/", project.id,  
                                                       "_", cluster_id, "_all_split_reads.rds"))
        
        # Load split read counts
        split_read_counts_local <- readRDS(file = paste0(local_results_folder, "/", project.id,  
                                                         "_", cluster_id, "_split_read_counts.rds"))
        
        
        if ( is.null(names(split_read_counts_local)) ) {
          split_read_counts_local <- split_read_counts_local %>%
            as_tibble(rownames = "junID")
        }
        
        if ( split_read_counts_local %>%
             mutate(total=rowSums(select_if(., is.numeric))) %>% 
             filter(total < supportive.reads) %>% 
             nrow() != 0 ) {
          paste0(project.id, " - There are split reads with less than the number of supportive reads indicated by parameter!") %>% print
        }
        
        split_read_counts_local <- split_read_counts_local %>%
          dplyr::select(c("junID", any_of(samples)))
        
        ## Only split reads with at least 2 supportive reads after filtering using the local cluster of samples
        split_read_counts_local <- split_read_counts_local %>%
          mutate(total = rowSums(select_if(., is.numeric))) %>% 
          filter(total >= supportive.reads) %>%
          dplyr::select(-total)
        
        
        all_split_reads_local <- all_split_reads_local %>%
          filter(junID %in% split_read_counts_local$junID)
        
        
        if ( j == 1 ) {
          split_read_counts_age <- split_read_counts_local
          all_split_reads_age <- all_split_reads_local
          j <- 2
          
        } else {
          
          ## The age stratification analysis is done per body site category,
          ## thus, only annotated introns overlapping the clusters from each
          ## body site are considered
          
          all_split_reads_age <- data.table::rbindlist(list(all_split_reads_local,
                                                            all_split_reads_age)) %>%
            distinct(junID, .keep_all = T)
          
          split_read_counts_age <- split_read_counts_local %>%
            inner_join(y = split_read_counts_age,
                       by = "junID")
        }
        
        
        #print(paste0(split_read_counts_age %>% names() %>% length(), " samples processed in total."))
        
        rm(all_split_reads_local)
        rm(split_read_counts_local)
        rm(samples)
        gc()
        
      }
      
    }
    
    
    
    print(paste0(Sys.time(), " - SAVING RESULTS ... "))
    
    split_read_counts_age <- split_read_counts_age[rowSums(is.na(split_read_counts_age[,2:ncol(split_read_counts_age)])) !=
                                                     (ncol(split_read_counts_age)-1),]
    
    all_split_reads_age <- all_split_reads_age %>%
      distinct(junID, .keep_all = T) %>%
      filter( junID %in% (split_read_counts_age$junID %>% unique()) )
    
    
    
    ############################################
    ## SAVE RESULTS FOR THE CURRENT AGE CLUSTER
    ############################################
    
    saveRDS(object = samples_age,
            file = paste0(local_results_folder, "/", project.id, "_", age_cluster, "_samples_used.rds"))
    
    ## Split read counts for age
    saveRDS(object = split_read_counts_age %>% data.table::as.data.table(),
            file = paste0(local_results_folder, "/", project.id, "_", age_cluster, "_split_read_counts.rds"))
    
    ## All split reads for age
    saveRDS(object = all_split_reads_age %>% data.table::as.data.table(),
            file = paste0(local_results_folder, "/", project.id, "_", age_cluster, "_all_split_reads.rds"))
    
    print(paste0(Sys.time(), " - ", age_cluster, " - ", project.id, " - RESULTS SAVED!"))
    
    
    return(samples_age)
    
  }
  
  all_samples_age_used %>% head()
  saveRDS(object = all_samples_age_used %>% unlist,
          file = paste0(local_results_folder, "/", project.id, "_age_all_samples_used.rds"))
  
  
  metadata_project <- readRDS(file = paste0(local_results_folder, "/", project.id, "_samples_metadata.rds") )
  
  ## Save age-related metadata
  saveRDS(object = metadata_project %>% 
            filter(external_id %in% (all_samples_age_used %>% unlist)) %>%
            inner_join(y = age.samples.clusters %>% dplyr::select(individual, age_group),
                       by = c("external_id" = "individual")) %>%
            dplyr::select(-cluster) %>%
            dplyr::rename(cluster = age_group),
          file = paste0(local_results_folder, "/", project.id, "_age_samples_metadata.rds"))
  
  
  ## Save age-related clusters
  saveRDS(object = age.groups %>% unique(),
          file = paste0(local_results_folder, "/", project.id, "_age_clusters_used.rds"))
  
  
  gc()
  
}