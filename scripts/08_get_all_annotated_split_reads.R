#' Title
#' Loops through the projects from the current recount3 project to obtain all original unique split reads found across their samples
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param all.clusters Clusters of samples. In GTEx projects, samples were clustered by tissue (eg. 'Puituitary', 'Thyroid', etc)
#' @param database.folder Local path to the folder that will contain the database and results files related with the database
#' @param results.folder Local path to the folder that contains the results of the analyses performed
#'
#' @return
#' @export
#'
#' @examples
get_all_annotated_split_reads <- function(recount3.project.IDs,
                                          all.clusters = NULL,
                                          database.folder,
                                          results.folder) {
  
  
  
  #############################################
  ## These are all split reads from all tissues
  ## obtained from 'USE ME' SAMPLES and samples with more than 6 RIN
  #############################################
  
  ## The sample filter in IntroVerse corresponded to 
  
  doParallel::registerDoParallel(10)
  
  all_split_reads_details_all_sample_clusters <- foreach(i = seq(length(recount3.project.IDs)), .combine = "rbind") %dopar%{
    
    project_id <- recount3.project.IDs[i]
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[2]
    
    folder_results_root <- paste0(results.folder, "/", project_id, "/")
    
      
    if ( is.null(all.clusters) && 
         file.exists(paste0(folder_results_root, "/base_data/", 
                            project_id, "_clusters_used.rds")) ) {
      
      all.clusters <-  readRDS(file = paste0(folder_results_root, "/base_data/", 
                                             project_id, "_clusters_used.rds"))
    } 
    
    if ( !is.null(all.clusters) ) {
      
      all_jxn_qc <- map_df(all.clusters, function(cluster) {
        
        # cluster <- all.clusters[1]
        message(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ...")
        
        if ( file.exists(paste0(folder_results_root, "/base_data/", 
                                project_id, "_", cluster, "_all_split_reads.rds")) ) {
          
          all_split_reads_details_105 <- readRDS(file = paste0(folder_results_root, "/base_data/", project_id, "_",
                                                               cluster, "_all_split_reads.rds"))
          
          all_split_reads_details_tidy <- all_split_reads_details_105 %>% 
            distinct(junID, .keep_all = T) %>% 
            as_tibble() %>%
            return()
          
        } else {
          return(NULL)
        }
        
      })
      
      if (all_jxn_qc %>% nrow() > 0 ) {
        
        all_jxn_qc %>%
          distinct(junID, .keep_all = T) %>% 
          return()
        
      } else {
        return(NULL)
      }
      
    }
    
  }
  
  
  message(Sys.time(), " - saving 'all_annotated_split_reads' for the database!")
  
  all_split_reads_details_all_sample_clusters <- all_split_reads_details_all_sample_clusters %>%
    distinct(junID, .keep_all = T)
  
  dir.create(path = database.folder, recursive = T, showWarnings = T)
  
  
  ## Save data
  saveRDS(object = all_split_reads_details_all_sample_clusters,
          file = paste0(database.folder, "/all_split_reads_qc_level2.rds") )
  
}