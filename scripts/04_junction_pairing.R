
#' Title
#' Function to pair novel junctions with annotated introns across the samples of each sample cluster 
#' (i.e. case/control, tissue, etc)
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615" 
#' @param project.name name given locally to the recount3 project.
#' A given recount3 project can be separated in multiple independent projects in recount3.
#' For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), 
#' but all of them belong to GTEx. Hence, it is useulf to have a folder named "project.name" that will contain
#' multiple subfolders named as the elements contained within the object 'recount3.project.IDs'
#' @param gtf.version Ensembl version (tested using Ensembl v105)
#' e.g. "105"
#'
#' @return
#' @export
#'
#' @examples
junction_pairing <- function(recount3.project.IDs, 
                             project.name,
                             gtf.version,
                             database.folder,
                             results.folder,
                             supporting.reads,
                             replace) {
  
  
  doParallel::registerDoParallel(10)
  
  foreach(i = seq(length(recount3.project.IDs)), .combine = "rbind") %dopar%{
    
    project_id <- recount3.project.IDs[i]
    
    # project_id <- recount3.project.IDs[1]
    
    folder_root <- paste0(results.folder, "/", project_id)
    folder_base_data <- paste0(folder_root, "/base_data/")
    
    print(paste0(Sys.time(), " - getting data from '", project_id, "' project..."))
    
    ## Load clusters
    
    if ( file.exists(paste0(folder_base_data, "/", project_id, "_samples_metadata.rds")) && 
         file.exists(paste0(folder_base_data, "/", project_id, "_clusters_used.rds")) ) {
      
      metadata.info <- readRDS(file = paste0(folder_base_data, "/", project_id, "_samples_metadata.rds"))
      clusters_ID <- readRDS(file = paste0(folder_base_data, "/", project_id, "_clusters_used.rds"))
      
      for (cluster_id in clusters_ID) {
        
        # cluster_id <- clusters_ID[1]
        # cluster_id <- clusters_ID[2]
        
        print(paste0(Sys.time(), " - loading '", cluster_id, "' source data ..."))
        
        ############################################
        ## LOAD DATA FOR THE CURRENT PROJECT
        ############################################
        
        ## Load samples
        samples_used <- readRDS(file = paste0(folder_base_data, "/", 
                                              project_id, "_", cluster_id, "_samples_used.rds"))
        
        ## Load split read data
        all_split_reads_details <- readRDS(file = paste0(folder_base_data, "/", project_id, "_", cluster_id, 
                                                         "_all_split_reads.rds")) %>% as_tibble()
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(folder_base_data, "/", project_id, "_", cluster_id, "_",
                                                   "split_read_counts.rds")) %>% as_tibble()
        
        stopifnot(
            split_read_counts %>% 
            mutate(sumCounts = rowSums(select(., !contains("junID")))) %>%
            filter(sumCounts < supporting.reads) %>% 
            nrow() == 0
        )
        
        
        if ( !identical(all_split_reads_details$junID, split_read_counts$junID) ) {
          print("ERROR! The number of junctions considered is not correct.")
          break;
        }
        
        ############################################
        ## DISTANCES SUITE OF FUNCTIONS
        ############################################
        
        folder_pairing_results <- paste0(folder_root, "/junction_pairing/", cluster_id, "/")
        dir.create(file.path(folder_pairing_results), recursive = TRUE, showWarnings = T)
        
        get_distances(cluster = cluster_id,
                      samples = samples_used,
                      split_read_counts = split_read_counts,
                      all_split_reads_details = all_split_reads_details,
                      folder_name = folder_pairing_results,
                      replace = replace)
        gc()
        
        
        extract_distances(cluster = cluster_id,
                          samples = samples_used,
                          folder_name = folder_pairing_results,
                          replace = replace)
        gc()
        
        
        get_never_misspliced(cluster = cluster_id,
                             samples = samples_used,
                             split_read_counts = split_read_counts,
                             all_split_reads_details = all_split_reads_details,
                             folder_name = folder_pairing_results,
                             replace = replace)
        
        rm(all_split_reads_details)
        rm(split_read_counts)
        gc()
        
        #}
      }
    }
    
  }
}


