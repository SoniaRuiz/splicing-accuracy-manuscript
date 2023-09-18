library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(tidyverse)
library(protr)

#' Title
#' Obtains the raw counts data from recount3 for a given projectID. 
#' It transforms the raw counts and calculates the TPM of the genes across all samples from the projectID.
#' As the raw counts represent the total coverage for a given gene within a given sample, raw counts need to be scaled by library size 
#' as longer genes may obtain higher number of counts, which does not mean higher expression.
#'
#' Then, it separates the TPM values across the samples of each sample cluster.
#' @param recount3.project.IDs Vector with all the recount3 projectIDs to process
#' @param data.source source of the data within recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#' @param results.folder Path to the local folder storing the results produced
#'
#' @return
#' @export
#'
#' @examples
generate_recount3_tpm <- function(recount3.project.IDs,
                                  data.source,
                                  results.folder) {
  
  
  ## Loop through the recount3 projects received by parameter
  doParallel::registerDoParallel(5)
  foreach(i = seq(length(recount3.project.IDs)) ) %dopar%{
    
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[18]
    project_id <- recount3.project.IDs[i]
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    results_folder_local_tpm <- paste0(results_folder_local, "/tpm/")
    dir.create(file.path(results_folder_local_tpm), recursive = TRUE, showWarnings = T)
    
    if ( !file.exists(paste0(here::here(), "/results/tpm/", project_id, "_tpm.rds")) ) {
    
      message("Downloading raw counts from ", project_id, "...")
      
      ## 1. Get expression data from recount3, transform raw counts and calculate TPM
      rse <- recount3::create_rse_manual(
        project = project_id,
        project_home = data.source,
        organism = "human",
        annotation = "gencode_v29",
        type = "gene")
      
      SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
      SummarizedExperiment::assayNames(rse) %>% print()
      
      message("Computing TPM for genes found in ", project_id)
      recount_tpm <- recount::getTPM(rse)
      
      #recount_tpm <- generate_tpm(rse = rse, ref_tidy)
      #recount_tpm %>% head()
      
      ## Save tpm values for all genes across all samples
      
      saveRDS(object = recount_tpm,
              file = paste0(here::here(), "/results/", project_id, "_tpm.rds"))
      
      ## Release some memory
      rm(rse)
      gc()
      
    } else {
      
      message("Loading TPM data for ", project_id, "...")
      recount_tpm <- readRDS(file = paste0(here::here(), "/results/tpm/", project_id, "_tpm.rds"))
      
    }
    
    
    ## 2. For each tissue within the current project, filter the RSE by its samples

    if ( file.exists(paste0(results_folder_local, "/base_data/", project_id, "_samples_metadata.rds")) &&
         file.exists(paste0(results_folder_local, "/base_data/", project_id, "_clusters_used.rds")) ) {

      metadata.info <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_samples_metadata.rds"))
      clusters_ID <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_clusters_used.rds")) %>% unique()

      for ( cluster_id in clusters_ID ) {

        # cluster_id <- clusters_ID[1]
        message("Getting results for ", cluster_id)

        cluster_samples <- readRDS(file = paste0(results_folder_local, "/base_data/",
                                                 project_id, "_", cluster_id, "_samples_used.rds"))

        ## Filter the object for the set of samples corresponding to the current cluster
        recount_tpm_local <- recount_tpm %>% #head %>%
          as_tibble(rownames = "gene_id") %>%
          #rename_with(~stringr::str_replace_all(., '\\.1', '')) %>%
          mutate(gene_id = stringr::str_replace_all(gene_id, '\\.*', '')) %>% 
          dplyr::select(c("gene_id", all_of(cluster_samples)))

        ## Save results
        saveRDS(object = recount_tpm_local,
                file = paste0(results_folder_local_tpm, project_id, "_", cluster_id, "_tpm.rds"))

        rm(recount_tpm_local)
        rm(cluster_samples)
        gc()

      }
    }


    ## 3. We do the same for the ageing analysis

    if ( file.exists(paste0(results_folder_local, "/base_data/", project_id, "_age_samples_metadata.rds")) &&
         file.exists(paste0(results_folder_local, "/base_data/", project_id, "_age_clusters_used.rds")) ) {

      metadata.info <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_age_samples_metadata.rds"))
      clusters_ID <- readRDS(file = paste0(results_folder_local, "/base_data/", project_id, "_age_clusters_used.rds")) %>% unique()

      for ( cluster_id in clusters_ID ) {

        # cluster_id <- clusters_ID[1]
        message("Getting results for ", cluster_id)

        cluster_samples <- readRDS(file = paste0(results_folder_local, "/base_data/",
                                                 project_id, "_", cluster_id, "_samples_used.rds"))

        ## Filter the object for the set of samples corresponding to the current cluster
        recount_tpm_local <- recount_tpm %>% #head %>% 
          as_tibble(rownames = "gene_id") %>%
          #rename_with(~stringr::str_replace_all(., '\\.1', '')) %>%
          mutate(gene_id = stringr::str_replace_all(gene_id, '\\.*', '')) %>%
          dplyr::select(c("gene_id", all_of(cluster_samples)))

        recount_tpm_local %>% head %>% print()

        ## Save results
        saveRDS(object = recount_tpm_local,
                file = paste0(results_folder_local_tpm, project_id, "_", cluster_id, "_tpm.rds"))

        rm(recount_tpm_local)
        rm(cluster_samples)
        gc()

      }
    }

    message(Sys.time(), " - ", project_id, " finished!")
    
    rm(folder_root)
    rm(clusters_ID)
    rm(recount_tpm)
    gc()
    
  }
}


# --------------------------------------------------------------------------------------------------------------
# source("/mnt/PROJECTS/splicing-accuracy-manuscript/scripts/10_generate_recount3_tpm.R")



# recount3_project_IDs <- c( "ADIPOSE_TISSUE", "ADRENAL_GLAND", "BLOOD", "BLOOD_VESSEL", "BONE_MARROW",
#                            "BRAIN",           "COLON",           "ESOPHAGUS",       "HEART",           "KIDNEY",
#                            "LIVER",           "LUNG",            "MUSCLE",          "NERVE",           "PANCREAS",
#                            "PITUITARY",       "SALIVARY_GLAND",  "SKIN",            "SMALL_INTESTINE", "SPLEEN",
#                            "STOMACH",         "THYROID" )
# 
# ## This is the name of the project producing the database
# dependencies_folder <- paste0("/mnt/PROJECTS/splicing-accuracy-manuscript/dependencies/")
# project_name <- "splicing_1read"
# data_source <- "data_sources/gtex"
# gtf_version <- 105
# results_folder <- file.path("/mnt/PROJECTS/splicing-accuracy-manuscript/results", paste0(project_name, "/", gtf_version, "/"))
# #file.path(here::here("results"), paste0(project_name, "/", gtf_version, "/"))
# 
# generate_recount3_tpm(recount3.project.IDs = recount3_project_IDs,
#                       data.source = data_source,
#                       results.folder = results_folder)

