library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(Biostrings)
library(tidyverse)
library(protr)

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/init.R")

base_folder <- here::here()
dependencies_folder <- paste0(base_folder, "/dependencies/")

#####################################
## LOAD SOURCE SCRIPTS
#####################################

setwd(file.path(base_folder,"scripts"))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(base_folder))

get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}

#####################################
## CALLS - PREPARE RECOUNT3 DATA
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_versions <- c(105)

## This is the name of the project producing the database
supportive_reads <- 1
project_name <- paste0("splicing_",supportive_reads,"read")
data_source <- "data_sources/gtex"

## Can be checked here: https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
all_projects <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLADDER",         "BLOOD",           "BLOOD_VESSEL",    "BONE_MARROW",
                   "BRAIN",           "BREAST",          "CERVIX_UTERI",    "COLON",           "ESOPHAGUS",       "FALLOPIAN_TUBE",
                   "HEART",           "KIDNEY",          "LIVER",           "LUNG",            "MUSCLE",          "NERVE",
                   "OVARY",           "PANCREAS",        "PITUITARY",       "PROSTATE",        "SALIVARY_GLAND",  "SKIN",
                   "SMALL_INTESTINE", "SPLEEN",          "STOMACH",         "TESTIS",          "THYROID",         "UTERUS",
                   "VAGINA"  )


recount3_project_IDs <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLOOD",           "BLOOD_VESSEL",    "BONE_MARROW",    
                           "BRAIN",           "COLON",           "ESOPHAGUS",       "HEART",           "KIDNEY",          
                           "LIVER",           "LUNG",            "MUSCLE",          "NERVE",           "PANCREAS",        
                           "PITUITARY",       "SALIVARY_GLAND",  "SKIN",            "SMALL_INTESTINE", "SPLEEN",          
                           "STOMACH",         "THYROID" )





#####################################
## INIT - DATABASE RECOUNT3 PROJECT
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_versions <- c(105)


for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  database_folder <- paste0(base_folder, "/database/", project_name, "/", gtf_version, "/")
  results_folder <- file.path(here::here("results"), paste0(project_name, "/", gtf_version, "/"))
  levelqc1_folder <- paste0(base_folder, "/database/")

  
  # download_recount3_data(recount3.project.IDs = all_projects,
  #                        project.name = project_name,
  #                        gtf.version = gtf_version,
  #                        data.source = data_source,
  #                        database.folder = database_folder,
  #                       levelqc1.folder = levelqc1_folder,
  #                        results.folder = results_folder)
   
  
  # prepare_recount3_data(recount3.project.IDs = recount3_project_IDs,
  #                       project.name = project_name,
  #                       data.source = data_source,
  #                       gtf.version = gtf_version,
  #                       database.folder = database_folder,
  #                       results.folder = results_folder,
  #                       subsampling = F,
  #                       levelqc1.folder = levelqc1_folder,
  #                       supporting.reads = supportive_reads,
  #                       num.cores = 5)

  junction_pairing(recount3.project.IDs = recount3_project_IDs,
                   project.name = project_name,
                   gtf.version = gtf_version,
                   database.folder = database_folder,
                   results.folder = results_folder,
                   supporting.reads = supportive_reads,
                   replace = T)


  get_all_annotated_split_reads(recount3.project.IDs = recount3_project_IDs,
                                database.folder = database_folder,
                                results.folder = results_folder)


  get_all_raw_jxn_pairings(recount3.project.IDs = recount3_project_IDs,
                           database.folder = database_folder,
                           results.folder = results_folder)

  
  ## Use the projects passing the filtering criteria established across functions above
  # recount3_project_IDs <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))
  
   



  # generate_transcript_biotype_percentage(recount3.project.IDs = recount3_project_IDs,
  #                                        project.name = project_name,
  #                                        gtf.version = gtf_version,
  #                                        database.folder = database_folder,
  #                                        results.folder = results_folder)


  # generate_recount3_tpm(recount3.project.IDs = recount3_project_IDs,
  #                       data.source = data_source,
  #                       results.folder = results_folder)
  
  
  # tidy_data_pior_sql(recount3.project.IDs = recount3_project_IDs,
  #                    project.name = project_name,
  #                    database.folder = database_folder,
  #                    results.folder = results_folder)
  # 
  
  # database_path <- paste0(database_folder,  "/", project_name, ".sqlite")
  # 
  # sql_database_generation(database.path = database_path,
  #                         recount3.project.IDs = recount3_project_IDs,
  #                         project.name = project_name,
  #                         gtf.version = gtf_version,
  #                         remove.all = T,
  #                         database.folder = database_folder,
  #                         results.folder = results_folder,
  #                         supportive.reads = supportive_reads)
  
  gc()
}
