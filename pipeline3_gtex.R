
#####################################
## 0 - Load libraries and source data
#####################################

library(stringr)
library(data.table)
library(RSQLite)
library(rtracklayer)
library(tidyverse)
library(ggforce)
library(GenomicRanges)
library(diffloop)
library(ggpubr)
library(matrixStats)
library(xlsx)

# source("/home/sruiz/PROJECTS/splicing-project/pipeline3_gtex.R")

source("/home/sruiz/PROJECTS/splicing-project/pipeline3_methods.R")
# source("/home/sruiz/PROJECTS/splicing-project/pipeline3_analyse_IDB.R")

folder_path <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/"


gtex_samples <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_used.rda")
gtex_tissues <-  readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")

get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}


# tissue <- "Adipose-Subcutaneous"
# tissue <- "Muscle-Skeletal" 
# tissue <- "WholeBlood"
# tissue <- "Liver"  
# tissue <- "Brain-CerebellarHemisphere"
# tissue <- "Cells-Transformedfibroblasts"
# tissue <- "Brain-Hippocampus"
# tissue <- "Brain-Amygdala"


#########################
## 2 - CALLS
#########################


## GENERAL GTEx CALL ------------------------------------------------------------------------

# if (!exists("CNC_CDTS_CONS_gr")) {
#   print(paste0(Sys.time(), " - loading conservation scores ... "))
#   load(file = "/home/egust/Projects/Alu_exonisation/results/CNC_CDTS_CONS_gr.rda")
#   
#   CNC_CDTS_CONS_gr <- CNC_CDTS_CONS_gr %>%
#     diffloop::rmchr()
#   
#   print(paste0(Sys.time(), " - conservation & CDTS scores loaded!"))
#   
# } else {
#   print(paste0(Sys.time(), " - 'CNC_CDTS_CONS_gr.rda' file already loaded!"))
# }

## Connect to the split-read counts DB
con <- dbConnect(RSQLite::SQLite(), "~/PROJECTS/splicing-project/data/splicing_tolerance.sqlite")
dbListTables(con)

versions <- c("76", "81", "90", "97", "104")

for (tissue in gtex_tissues) {
  
  for (version in versions[5]) {
    
    # tissue <- gtex_tissues[11]
    # tissue <- gtex_tissues[40]
    # tissue <- (gtex_tissues[7:17])[2]
    # version <- versions[5]
    
    print(tissue)
    print(version)
    
    #################################
    ## LOAD TISSUE DATA
    #################################
    
    
    samples <- gtex_samples[[tissue]]
    
    
    ## SET RESULTS FOLDER
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/distances/", tissue, "/v", version)
    dir.create(file.path(folder_name), showWarnings = F, recursive = T)
    
    
    # query <- paste0("SELECT ", paste(dbQuoteIdentifier(con, c("juncID", samples %>% as.character())), collapse = ","), " FROM '", tissue, "'")
    # split_read_counts <- dbGetQuery(con, query) %>%
    #   dplyr::rename(junID = juncID) %>%
    #   dplyr::as_tibble()
    
    # ## Load the counts file
    split_read_counts <- fread(paste0("/data/recount/GTEx_SRP012682/gtex_full_split_read_count_table/", tissue, "/", tissue, ".csv"))
    if (any(split_read_counts %>% names() == "juncID")) {
      split_read_counts <- split_read_counts %>% dplyr::rename(junID = juncID) %>% data.table::as.data.table()
    }
    
    # ## LOAD THE ANNOTATED SPLIT READS FILE
    all_split_reads_details_104 <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/base_data/", tissue, "/",
                                                         tissue, "_annotated_SR_details_length_v", version, ".rds")) %>% data.table::as.data.table()
    
    # all_split_reads_details_104[1,]
    # all_split_reads_details_104$seqnames %>% unique()
    
    
    #################################
    ## GET THE DISTANCES
    #################################
    
    # get_distances(cluster = tissue,
    #               samples = samples,
    #               split_read_counts = split_read_counts,
    #               all_split_reads_details_104 = all_split_reads_details_104,
    #               folder_name = folder_name)
    
    extract_distances(cluster = tissue,
                      samples = samples,
                      split_read_counts = split_read_counts,
                      folder_name = folder_name)
    
    # get_never_misspliced(cluster = tissue,
    #                      samples = samples,
    #                      split_read_counts = split_read_counts,
    #                      all_split_reads_details_104 = all_split_reads_details_104,
    #                      folder_name = folder_name)
    
    extract_never_misspliced(cluster = tissue,
                             samples = samples,
                             split_read_counts = split_read_counts,
                             folder_name = folder_name)
    
    add_never_misspliced_to_df(cluster = tissue,
                               samples = samples,
                               all_split_reads_details_104 = all_split_reads_details_104,
                               folder_name = folder_name)
    
    
    #################################
    ## GET MIS-SPLICING RATIO
    #################################
    
    folder_idb_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue, "/v", version)
    dir.create(file.path(folder_idb_name), showWarnings = F, recursive = T)
    
    get_missplicing_ratio(cluster = tissue,
                          samples = samples,
                          split_read_counts = split_read_counts,
                          folder_name = folder_name,
                          folder_save_name = folder_idb_name)
    
    
    
    add_missplicing_class_to_df(cluster = tissue,
                                folder_name = folder_idb_name)
    
    
    
    get_missplicing_QC(cluster = tissue,
                       all_split_reads_details_104 = all_split_reads_details_104,
                       samples = samples,
                       split_read_counts = split_read_counts,
                       folder_name = folder_idb_name)
    
    
    
    #################################
    ## ADD FEATURES TO THE IDB
    #################################
    
    add_intron_type(cluster = tissue,
                    folder_name = folder_idb_name)
    
    remove_MT_genes(cluster = tissue,
                    folder_name = folder_idb_name)
    
    add_cdts_cons_scores(cluster = tissue,
                         CNC_CDTS_CONS_gr = CNC_CDTS_CONS_gr,
                         folder_name = folder_idb_name)
    
    
    clinvar_analysis(cluster = tissue,
                     folder_name = folder_idb_name)
    
    
    # add_gene_tpm_length(cluster = tissue,
    #                     folder_name = folder_idb_name,
    #                     GTEx = T)
    
    # add_MANE_info(cluster = tissue,
    #               folder_name = folder_idb_name)
    # 
    # 
    # #################################
    # ## ADD PROTEIN PERCENTAGE
    # #################################
    # 
    # 
    # protein_folder <- "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_annotated_SR_details_length_104_biotype.rds"
    # folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/"
    # 
    # 
    # add_protein_percentage_to_database(cluster = tissue,
    #                                    protein_folder = protein_folder,
    #                                    folder_root = folder_root)
    
    
  }
  
  rm(all_split_reads_details_104)
  rm(split_read_counts)
  gc()
  
}


