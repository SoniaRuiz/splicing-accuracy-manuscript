library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/init.R")

#biomaRt::biomartCacheClear()

.libPaths( c( "/home/sruiz/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()) )
Sys.setenv("AWS_ACCESS_KEY_ID"="AKIAUKEFAIIIETDHJZOI",
           "AWS_SECRET_ACCESS_KEY"="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",
           "AWS_DEFAULT_REGION"="eu-west-2")


projects_used <- readRDS(file = "~/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
clusters_used <- readRDS(file = "~/PROJECTS/splicing-project/splicing-recount3-projects/all_clusters_used.rds")


source("/home/sruiz/PROJECTS/splicing-project-recount3/sql_helper.R")
source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline2_annotate_from_recount.R")
#source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_sql_generation.R")
source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_SQL_generation_extra.R")

#####################################
## FUNCTIONS - PREPARE RECOUNT3 DATA
#####################################

init_recount3_gtex_data <- function (gtf_version) {
  
  # project_id <- projects_used[6]
  all_projects <- readRDS(file = "~/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  if (file.exists(paste0("~/PROJECTS/splicing-project-recount3/database/all_split_reads_raw.rds"))) {
    
    print(paste0(Sys.time(), " - loading 'all_split_reads_raw.rds' file..."))
    all_split_reads_raw_tidy <- readRDS(file = paste0("~/PROJECTS/splicing-project-recount3/database/all_split_reads_raw.rds"))
    all_split_reads_raw_tidy %>% distinct(junID) %>% nrow() %>% print()
    
  } else {
    
    all_split_reads_raw <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
      
      folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                            project_id, "/")
      dir.create(file.path(folder_root), recursive = TRUE, showWarnings = T)
      
      rse <- recount3::create_rse_manual(
        project = project_id,
        project_home = "data_sources/gtex",
        organism = "human",
        annotation = "gencode_v29",
        type = "jxn"
      )
      
      all_split_reads <- data.frame(chr = rse %>% SummarizedExperiment::seqnames(),
                                    start = rse %>% SummarizedExperiment::start(),
                                    end = rse %>% SummarizedExperiment::end(),
                                    strand = rse %>% SummarizedExperiment::strand(),
                                    width = rse %>% SummarizedExperiment::width(),
                                    junID = rse %>% rownames())
      all_split_reads %>% nrow()
      saveRDS(object = all_split_reads %>% data.table::as.data.table(),
              file = paste0(folder_root, "/all_split_reads_raw.rds"))
      
      all_split_reads %>% data.table::as.data.table() %>%
        return()
    })
    
    all_split_reads_raw_tidy <- all_split_reads_raw %>%
      mutate(junID_tidy = paste0(chr, ":", start, "-", end, ":", strand)) %>%
      distinct(junID_tidy, .keep_all = T)
    
    all_split_reads_raw_tidy %>% nrow()
    
    folder_path <- paste0("~/PROJECTS/splicing-project-recount3/database/")
    dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
    saveRDS(object = all_split_reads_raw_tidy %>% data.table::as.data.table(),
            file = paste0(folder_path, "/all_split_reads_raw.rds"))
  }
  
  #######################################################################
  ## Remove split reads located in unplaced sequences in the chromosomes
  #######################################################################
  
  print(paste0(Sys.time(), " - removing unplaced sequences genome..."))
  
  all_split_reads_raw_tidy_gr <- all_split_reads_raw_tidy %>%
    GenomicRanges::GRanges() %>%
    diffloop::rmchr()
  
  
  all_split_reads_raw_tidy_gr %>% length()
  all_split_reads_raw_tidy_gr <- GenomeInfoDb::keepSeqlevels(x = all_split_reads_raw_tidy_gr,
                                                             value = intersect(all_split_reads_raw_tidy_gr %>% 
                                                                                 GenomeInfoDb::seqnames() %>% 
                                                                                 levels(), 
                                                                               c( "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                                                  "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X", "Y")), 
                                                             pruning.mode = "tidy")
  all_split_reads_raw_tidy_gr %>% head()
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T) %>% 
    nrow() %>%
    print()
  
  
  #######################################################
  ## Remove split reads overlapping the ENCODE backlist
  #######################################################
  
  print(paste0(Sys.time(), " - removing blacklist sequences..."))
  
  
  blacklist_path <- "/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed"
  all_split_reads_raw_tidy_gr <- remove_encode_blacklist_regions(GRdata = all_split_reads_raw_tidy_gr,
                                                                 blacklist_path = blacklist_path)
  
  
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T) %>% 
    nrow() %>% 
    print()
  
  
  #######################################################
  ## Anotate using 'dasper'
  #######################################################
  
  print(paste0(Sys.time(), " - annotating dasper..."))
  
  gtf_path <- paste0("/data/references/ensembl/gtf/v", gtf_version, "/Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")
  edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  
  all_split_reads_details_w_symbol <- dasper::junction_annot(junctions = all_split_reads_raw_tidy_gr %>% GRanges(), 
                                                                 ref = edb,
                                                                 ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"), 
                                                                 ref_cols_to_merge = c("gene_id", "gene_name", "tx_id"))
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol %>% 
    as_tibble()
  
  all_split_reads_details_w_symbol_reduced <- all_split_reads_details_w_symbol %>% 
    dplyr::select(seqnames, start, end, strand, junID, 
                  gene_id = gene_id_junction, in_ref, type,
                  tx_id_junction )
  
  folder_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/")
  dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
  saveRDS(all_split_reads_details_w_symbol_reduced,
          file = paste0(folder_path, "/all_split_reads_details_w_symbol.rds"))
  
  
  
  
  ################################################################################
  ## Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  print(paste0(Sys.time(), " - removing split reads not classified as 'annotated', 'novel_donor' or 'novel_acceptor'..."))
  
  all_split_reads_details_w_symbol_reduced <- readRDS(file = paste0(folder_path, "/all_split_reads_details_w_symbol.rds"))
  
  all_split_reads_details_w_symbol_reduced_discard <- all_split_reads_details_w_symbol_reduced %>%
    dplyr::filter(!(type %in% c("annotated", "novel_donor", "novel_acceptor"))) %>%
    data.table::as.data.table()
  all_split_reads_details_w_symbol_reduced_discard %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type) %>%
    print()
  all_split_reads_details_w_symbol_reduced_discard %>%
    distinct(junID) %>%
    nrow() %>%
    print()
  
  
  all_split_reads_details_w_symbol_reduced_keep <- all_split_reads_details_w_symbol_reduced %>%
    data.table::as.data.table() %>%
    dplyr::filter(!(junID %in% all_split_reads_details_w_symbol_reduced_discard$junID))
  
  all_split_reads_details_w_symbol_reduced_keep %>% nrow() %>% print()
  
  ############################################
  ## Discard all junctions shorter than 25bp
  ############################################
  
  print(paste0(Sys.time(), " - removing split reads shorter than 25bp..."))
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- all_split_reads_details_w_symbol_reduced_keep %>%
    GRanges() %>%
    as_tibble()
  
  all_split_reads_details_w_symbol_reduced_keep_gr  %>% 
    dplyr::filter(width < 25) %>%
    nrow() %>% 
    print()
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- all_split_reads_details_w_symbol_reduced_keep_gr  %>% 
    dplyr::filter(width >=25)
  
  all_split_reads_details_w_symbol_reduced_keep_gr %>% nrow() %>% print()
  
  ############################################################################
  ## Discard all introns assigned to multiple genes (i.e. ambiguous introns)
  ############################################################################
  
  print(paste0(Sys.time(), " - discarding ambiguous introns..."))
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- all_split_reads_details_w_symbol_reduced_keep_gr %>%
    distinct(junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  ambiguous_introns <- all_split_reads_details_w_symbol_reduced_keep_gr %>%
    dplyr::filter(ambiguous == T)
  
  ambiguous_introns %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type) %>% 
    print()
  
  ambiguous_introns %>% 
    distinct(junID, .keep_all = T) %>%
    print()
  
  saveRDS(object = ambiguous_introns,
          file = paste0(folder_path, "/all_ambiguous_introns.rds"))
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- all_split_reads_details_w_symbol_reduced_keep_gr %>%
    dplyr::filter(ambiguous == F) %>%
    dplyr::select(-ambiguous)
  
  all_split_reads_details_w_symbol_reduced_keep_gr
  
  saveRDS(object = all_split_reads_details_w_symbol_reduced_keep_gr,
          file = paste0(folder_path, "/all_split_reads_details_w_symbol_reduced_keep.rds"))
  
  
  
}


tidy_recount3_data_per_tissue <- function(projects_used, 
                                          main_project,
                                          gtf_version) {
  
  folder_main <- paste0("~/PROJECTS/splicing-project-recount3/database/v", 
                        gtf_version, "/")
  all_split_reads_details_w_symbol_reduced_keep_gr <- readRDS(file = paste0(folder_main, "/all_split_reads_details_",
                                                                            gtf_version, "_w_symbol_reduced_keep.rds"))
  
  # ## Generate the raw file
  for (project_id in projects_used) {
    
    # project_id <- projects_used[10]
    
    folder_root <- paste0("~/PROJECTS/splicing-project/splicing-recount3-projects/",
                          project_id, "/v", gtf_version, "/", main_project, "_project/")
    folder_path <- paste0(folder_root, "raw_data/")
    dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
    
    print(paste0(Sys.time(), " - getting junction data from recount3 - '", project_id, "' tissue..."))
    
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      type = "jxn"
    )
    
    metadata.info <- rse %>% 
      SummarizedExperiment::colData() %>%
      as_tibble() %>%
      filter(gtex.smrin >= 6.0,
             gtex.smafrze != "EXCLUDE")
    
    saveRDS(object = metadata.info, file = paste0(folder_path, "/samples_metadata.rds"))
    
    clusters_ID <- rse$gtex.smtsd %>% unique()
    
    for (cluster_id in clusters_ID) {
      
      # cluster_id <- clusters_ID[2]
      print(paste0(Sys.time(), " - filtering junction data by cluster - '", cluster_id, "' tissue..."))
      
      cluster_samples <- metadata.info %>% 
        as_tibble() %>%
        filter(gtex.smtsd == cluster_id,
               gtex.smrin >= 6.0,
               gtex.smafrze != "EXCLUDE") %>%
        distinct(external_id) %>% 
        pull()
      
      if (cluster_samples %>% length() >= 70) {
        
        saveRDS(object = cluster_samples, 
                file = paste0(folder_path, "/", project_id, "_", cluster_id, "_samples_used.rds"))
        
        local_rse <- rse[, rse$external_id %in% cluster_samples]
        
        print(paste0(Sys.time(), " - getting split read counts..."))
        counts <- (local_rse %>% SummarizedExperiment::assays())[[1]]
        counts <- counts[(rownames(counts) %in% all_split_reads_details_w_symbol_reduced_keep_gr$junID),]
        counts <- counts[rowSums(counts) > 0, ]
        counts <- counts %>% as.matrix()
        counts %>% rownames()
        counts %>% nrow()
        saveRDS(object = counts,
                file = paste0(folder_path, "/", project_id, "_", cluster_id, "_split_read_counts_sample_tidy.rds"))
        
        
        print(paste0(Sys.time(), " - getting split read IDs..."))
        all_split_reads <- data.frame(junID = local_rse %>% rownames()) %>%
          filter(junID %in% (counts %>% rownames()))
        all_split_reads %>% nrow()
        
        all_split_reads_tidy <- all_split_reads %>%
          data.table::as.data.table() %>%
          inner_join(y = all_split_reads_details_w_symbol_reduced_keep_gr,
                     by = "junID")
        
        all_split_reads_tidy %>% nrow()
        
        if (any(all_split_reads_tidy$width < 25) | 
            any(str_detect(string = all_split_reads_tidy$chr, pattern = "random")) |
            any(str_detect(string = str_to_lower(all_split_reads_tidy$chr), pattern = "u")) |
            any(!(all_split_reads_tidy$type %in% c("annotated",
                                                   "novel_acceptor",
                                                   "novel_donor")))) {
          print("ERROR! The split reads do not meet the minimum filtering criteria.")
          break;
        }
        
        if (setdiff(rownames(counts), all_split_reads_tidy$junID) %>% length() > 0) {
          print("ERROR! The split reads are not identical.")
          break;
        }
        
        saveRDS(object = all_split_reads_tidy %>% data.table::as.data.table(),
                file = paste0(folder_path, "/", project_id, "_", cluster_id, "_all_split_reads_sample_tidy.rds"))
        
      }
    } 
    
    rm(rse)
  }
  
  ## Save clusters_used, as it will be part of the database
  clusters_used <- data.frame(project_id = as.character(),
                              cluster_id = as.character(),
                              sample_id = as.character())
  
  for (project_id in projects_used) {
    
    # project_id <- projects_used[10]
    
    folder_root <- paste0("~/PROJECTS/splicing-project/splicing-recount3-projects/",
                          project_id, "/v", gtf_version, "/", main_project, "_project/")
    folder_path <- paste0(folder_root, "/raw_data/")
    
    print(paste0(Sys.time(), " - getting junction data from recount3 - '", project_id, "' tissue..."))
    
    metadata.info <- readRDS(file = paste0(folder_path, "/samples_metadata.rds"))
    
    for (cluster_id in (metadata.info$gtex.smtsd %>% unique())) {
      
      # cluster_id <- clusters_ID[1]
      print(paste0(Sys.time(), " - filtering junction data by cluster - '", cluster_id, "' tissue..."))
      
      cluster_samples <- metadata.info %>% 
        as_tibble() %>%
        filter(gtex.smtsd == cluster_id,
               gtex.smrin >= 6.0,
               gtex.smafrze != "EXCLUDE") %>%
        distinct(external_id) %>% 
        pull()
      
      if (cluster_samples %>% length() >= 70) {
        
        clusters_used <- rbind(clusters_used,
                               data.frame(project = project_id,
                                          cluster = cluster_id,
                                          sample_id = cluster_samples))
        
        
      }
    }  
  }
  
  dir.create(file.path(paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", 
                              gtf_version, "/", main_project)), recursive = TRUE, showWarnings = T)
  saveRDS(clusters_used,
          file =  paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", 
                         gtf_version, "/", main_project, "/metadata.rds"))
}


junction_pairing <- function(projects_used, 
                             main_project) {
  source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_junction_pairing.R")
  
  for (project_id in projects_used) {
    
    # project_id <- projects_used[9]
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          project_id, "/", main_project, "_project/")
    folder_path <- paste0(folder_root, "/raw_data/")
    
    print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
    
    ## Load clusters
    metadata.info <- readRDS(file = paste0(folder_path, "/samples_metadata.rds"))
    clusters_ID <- metadata.info$gtex.smtsd %>% unique()
    
    
    for (cluster_id in clusters_ID) {
    
      # cluster_id <- clusters_ID[1]
      # cluster <- clusters_ID[2]
      
      print(paste0(Sys.time(), " - loading '", cluster_id, "' source data ..."))
      
      ############################################
      ## LOAD DATA FOR THE CURRENT PROJECT
      ############################################
      
      ## Load samples
      samples <- metadata.info %>% 
        as_tibble() %>%
        filter(gtex.smtsd == cluster_id,
               gtex.smrin >= 6.0,
               gtex.smafrze != "EXCLUDE") %>%
        distinct(external_id) %>% 
        pull()
      
      if (samples %>% length() >= 70) {
        
        samples_used <- readRDS(file = paste0(folder_path, "/", project_id, "_", cluster_id, "_samples_used.rds"))
        
        if (!identical(samples %>% sort(), samples_used %>% sort())) {
          print("ERROR - samples are not identical.")
          break;
        }
      
        folder_name <- paste0(folder_root, "/results/", cluster_id, "/distances/v", gtf_version, "/")
        dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)
        
        
        ## Load split read data
        all_split_reads_details <- readRDS(file = paste0(folder_path, "/", project_id, "_", cluster_id, 
                                                         "_all_split_reads_sample_tidy.rds"))
        
        
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(folder_path, "/", project_id, "_", cluster_id, "_",
                                                   "split_read_counts_sample_tidy.rds")) %>% 
          as_tibble(rownames = "junID")
       
        
        ############################################
        ## DISTANCES SUITE OF FUNCTIONS
        ############################################
        
        get_distances(cluster = cluster_id,
                      samples = samples,
                      split_read_counts = split_read_counts,
                      all_split_reads_details = all_split_reads_details,
                      folder_name)
        
        
        extract_distances(cluster = cluster_id,
                          samples,
                          split_read_counts = split_read_counts,
                          folder_name = folder_name)
        
        
        get_never_misspliced(cluster = cluster_id,
                             samples = samples,
                             split_read_counts = split_read_counts,
                             all_split_reads_details = all_split_reads_details,
                             folder_name = folder_name,
                             save_results = T)
        
        rm(all_split_reads_details)
        rm(split_read_counts)
        gc()
      }
    }
  }
}


tidy_data_pior_sql <- function (projects_used, 
                                gtf_version,
                                all_clusters,
                                main_project) {
  
  ## Load main object
  all_split_reads_details_w_symbol_reduced_keep_gr <- readRDS(file = paste0("~/PROJECTS/splicing-project-recount3/database/v",
                                                                            gtf_version, "/all_split_reads_details_",
                                                                            gtf_version, "_w_symbol_reduced_keep.rds"))
  
  ############################################
  ## Discard all junctions from EXCLUDE ME samples
  ## Sex-related tissues
  ## Samples with RIN < 6
  ## Tissues with less than 70 samples
  ############################################
  
  all_split_reads_details <- readRDS(file = paste0("~/PROJECTS/splicing-project-recount3/database/v",
                                                   gtf_version, "/", main_project, 
                                                   "/all_split_reads_105_length_all_tissues.rds"))
  
  
  
  ## This should be zero
  setdiff(all_split_reads_details$junID, all_split_reads_details_w_symbol_reduced_keep_gr$junID)
  
  ## These are the junctions from EXCLUDE ME samples, sex-related tissues, and tissues with less than 70 samples
  setdiff(all_split_reads_details_w_symbol_reduced_keep_gr$junID, all_split_reads_details$junID) %>% unique %>% length()
  
  
  all_split_reads_details %>%
    as.data.table() %>%
    dplyr::count(type)
  all_split_reads_details %>% distinct(junID) %>% nrow()
  
  ############################################
  ## QC
  ############################################ 
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_w_symbol_reduced_keep_gr$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_w_symbol_reduced_keep_gr[ind, "junID"] <- str_replace(string = all_split_reads_details_w_symbol_reduced_keep_gr[ind, "junID"]$junID, 
                                                                                  pattern = "\\*", 
                                                                                  replacement = all_split_reads_details_w_symbol_reduced_keep_gr[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_w_symbol_reduced_keep_gr$junID, pattern = "\\*"))
  }
  
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details[ind, "junID"] <- str_replace(string = all_split_reads_details[ind, "junID"]$junID, 
                                                         pattern = "\\*", 
                                                         replacement = all_split_reads_details[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details$junID, pattern = "\\*"))
  }
  
  
  ##########################################
  ## Load all distances pairings
  ##########################################
  
  df_all_distances_pairings_raw <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v",
                                                         gtf_version, "/",
                                                         main_project, "/df_all_distances_pairings_raw.rds"))
  
  ## QC
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = df_all_distances_pairings_raw$ref_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_distances_pairings_raw[ind, "ref_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "ref_junID"]$ref_junID, 
                                                                   pattern = "\\*", 
                                                                   replacement = df_all_distances_pairings_raw[ind, "ref_strand"]$ref_strand %>% as.character())
  }
  ## Remove potential * in the junID of the novel junctions
  ind <- which(str_detect(string = df_all_distances_pairings_raw$novel_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_distances_pairings_raw[ind, "novel_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "novel_junID"]$novel_junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = df_all_distances_pairings_raw[ind, "novel_strand"]$novel_strand  %>% as.character())
  }
  any(str_detect(df_all_distances_pairings_raw$ref_junID, pattern = "\\*"))
  any(str_detect(df_all_distances_pairings_raw$novel_junID, pattern = "\\*"))
  
  
  ##########################################
  ## Get never mis-spliced
  #########################################
  
  df_never_misspliced <- get_intron_never_misspliced(all_projects = projects_used,
                                                     all_clusters = all_clusters,
                                                     main_project)
  
  ## Remove the introns paired with novel junctions
  df_never_misspliced_tidy <- df_never_misspliced %>%
    dplyr::filter(!(ref_junID %in% df_all_distances_pairings_raw$ref_junID)) %>%
    as_tibble()
  
  df_never_misspliced_tidy %>% distinct(ref_junID) %>% as_tibble()
  
  
  
  ############################################
  ## GET all not paired
  ############################################ 
  
  # df_all_novel_raw_tidy
  df_not_paired <- all_split_reads_details %>%
    as.data.table() %>%
    dplyr::filter(!(junID %in% c(df_all_distances_pairings_raw$ref_junID,
                                 df_all_distances_pairings_raw$novel_junID)))
  
  ## These are all the non-paired, including the never mis-spliced. Thus, [768,646 - 38,521 = 730125]
  df_not_paired %>% 
    distinct(junID)
  df_not_paired %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  
  if (any(str_detect(df_never_misspliced_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(df_not_paired$junID, pattern = "\\*"))) {
    print("ERROR!")
  }
  
  ## All never mis-spliced should be categorised as not paired.
  ## Thus, this should be zero
  setdiff(df_never_misspliced_tidy$ref_junID, df_not_paired %>% dplyr::filter(type == "annotated") %>% pull(junID))
  
  df_not_paired_tidy <- df_not_paired %>%
    dplyr::filter(!(junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    distinct(junID, .keep_all = T) %>%
    as_tibble()
  
  df_not_paired_tidy %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  df_not_paired_tidy %>% distinct(junID)
  
  ## This should be zero
  intersect(df_not_paired_tidy$junID, df_all_distances_pairings_raw$novel_junID) %>% length()
  
  
  df_all_distances_pairings_raw %>% distinct(novel_junID) %>% nrow() +
    df_all_distances_pairings_raw %>% distinct(ref_junID) %>% nrow() +
    df_never_misspliced_tidy %>% distinct(ref_junID) %>% nrow()
  
  ##########################################
  ## Remove ambiguous junctions
  ##########################################
  
  ## All these should be zero
  if( intersect(df_not_paired_tidy$junID, df_all_distances_pairings_raw$novel_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_all_distances_pairings_raw$ref_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_never_misspliced_tidy$ref_junID) %>% length() > 0) {
    print("ERROR!")
  }
  
  
  df_ambiguous_novel <- df_all_distances_pairings_raw %>%
    dplyr::filter(!(novel_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    group_by(novel_junID) %>%
    mutate(distances_sd = distance %>% sd()) %>%
    dplyr::filter(distances_sd > 0)
  
  
  df_ambiguous_novel %>% 
    ungroup() %>%
    distinct(novel_junID)
  
  df_ambiguous_novel %>% 
    ungroup() %>%
    distinct(ref_junID)
  
  
  saveRDS(df_ambiguous_novel,
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", 
                        gtf_version, "/", main_project, "/df_all_tissues_raw_distances_ambiguous.rds"))
  
  
  ## Remove ambiguous junctions
  df_all_distances_pairings_raw_tidy <- df_all_distances_pairings_raw %>%
    dplyr::filter(!(novel_junID %in% df_ambiguous_novel$novel_junID)) %>%
    as.data.table() %>%
    distinct(novel_junID, ref_junID, .keep_all = T) %>%
    mutate(ref_strand = ref_strand %>% as.character(),
           novel_strand = novel_strand %>% as.character()) 
  
  
  
  df_all_distances_pairings_raw_tidy %>%
    dplyr::distinct(novel_junID)
  df_all_distances_pairings_raw_tidy %>%
    dplyr::distinct(ref_junID)
  
  
  
  
  (df_ambiguous_novel %>% 
      ungroup() %>%
      distinct(ref_junID) %>% nrow()) - (intersect(c(df_all_distances_pairings_raw_tidy$novel_junID,
                                                     df_all_distances_pairings_raw_tidy$ref_junID),
                                                   df_ambiguous_novel %>% 
                                                     ungroup() %>%
                                                     distinct(ref_junID) %>% pull) %>% length())
  
  
  df_all_distances_pairings_raw_tidy %>%
    dplyr::distinct(novel_junID) %>% nrow() +
    df_all_distances_pairings_raw_tidy %>%
    dplyr::distinct(ref_junID) %>% nrow() +
    df_never_misspliced_tidy %>% distinct(ref_junID) %>% nrow()
  
  
  ##############################################################################
  ## SAVE FINAL OBJECT
  ##############################################################################
  
  ## 1. DISTANCES PAIRINGS
  
  if (any(str_detect(string = df_all_distances_pairings_raw_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_all_distances_pairings_raw_tidy$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  
  df_all_distances_pairings_raw_tidy <- df_all_distances_pairings_raw_tidy %>%
    inner_join(y = all_split_reads_details %>% dplyr::select(junID, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  
  saveRDS(object = df_all_distances_pairings_raw_tidy,
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", 
                        gtf_version, "/",
                        main_project, "/all_paired_intron_novel_tidy.rds"))
  
  
  ## 2. NEVER MIS-SPLICED
  
  if (any(str_detect(string = df_never_misspliced_tidy$ref_junID, pattern = "\\*")) ) {
    print("ERROR! Some NEVER MIS-SPLICED junctions still have a * in their IDs!")
  }
  df_never_misspliced_tidy <- df_never_misspliced_tidy %>%
    inner_join(y = all_split_reads_details %>% 
                 dplyr::select(junID, 
                               seqnames, start, end, width, strand,
                               gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  
  saveRDS(object = df_never_misspliced_tidy,
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", 
                        gtf_version, "/",
                        main_project, "/df_all_nevermisspliced_introns.rds"))
  
  
  ## AMBIGUOUS JUNCTIONS
  if (any(str_detect(string = df_ambiguous_novel$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_ambiguous_novel$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  saveRDS(df_ambiguous_novel,
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", 
                        gtf_version, "/",
                        main_project, "/df_all_tissues_raw_distances_ambiguous.rds"))
}


sql_database_generation <- function(database_path,
                                    projects_used, 
                                    main_project,
                                    gtf_version = gtf_version,
                                    remove_all) {
  

  remove_tables(database_path, remove_all)
  
  if (remove_all) {
    create_metadata_table(database_path,
                          main_project = main_project,
                          gtf_version = gtf_version,
                          SRA_projects = projects_used)
    
    create_mane_table(database_path)
  
    create_master_tables(database_path,
                         main_project = main_project,
                         gtf_version = gtf_version)
  } 
  
  create_cluster_tables(database_path,
                        gtf_version = gtf_version,
                        main_project = main_project)

}

#####################################
## CALLS - PREPARE RECOUNT3 DATA
#####################################

gtf_version <- 105
main_project <- "splicing"
database_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v",
                        gtf_version, "/", main_project, "/", main_project, ".sqlite")

sql_database_generation(database_path,
                        projects_used,
                        main_project,
                        gtf_version,
                        remove_all = T)

# tidy_recount3_data_per_tissue(projects_used,
#                               gtf_version = 97)
