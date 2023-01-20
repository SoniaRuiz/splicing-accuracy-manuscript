library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/init.R")

#biomaRt::biomartCacheClear()

.libPaths( c( "/home/sruiz/R/x86_64-pc-linux-gnu-library/4.0/", .libPaths()) )

setwd("~/PROJECTS/splicing-project-recount3/")

dependencies_folder <- paste0(getwd(), "/../introverse-app/database_generation/dependencies/")

source(paste0(getwd(), "/pipeline3_junction_pairing.R"))
source(paste0(getwd(), "/pipeline3_database_SQL_helper.R"))
source(paste0(getwd(), "/pipeline3_database_SQL_generation.R"))



#####################################
## FUNCTIONS - PREPARE RECOUNT3 DATA
#####################################

init_recount3_gtex_data <- function (projects_used,
                                     gtf_version) {
  
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  if ( file.exists(paste0(getwd(), "/database/all_split_reads_raw.rds")) ) {
    
    print(paste0(Sys.time(), " - loading 'all_split_reads_raw.rds' file..."))
    all_split_reads_raw_tidy <- readRDS(file = paste0(getwd(), "/database/all_split_reads_raw.rds"))
    
  } else {
    
    all_split_reads_raw <- map_df(projects_used, function(project_id) {
      
      # project_id <- projects_used[31]
      print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
      
      folder_root <- paste0(getwd(), "/results/", project_id, "/")
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
      rm(rse)
      gc()
      
      all_split_reads %>% nrow()
      saveRDS(object = all_split_reads %>% data.table::as.data.table(),
              file = paste0(folder_root, "/all_split_reads_raw.rds"))
      
      all_split_reads %>% data.table::as.data.table() %>%
        return()
    })
    
    all_split_reads_raw_tidy <- all_split_reads_raw %>%
      distinct(junID, .keep_all = T)
    
    all_split_reads_raw_tidy %>% nrow()
    
    folder_path <- paste0(getwd(), "/database/")
    dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
    saveRDS(object = all_split_reads_raw_tidy %>% data.table::as.data.table(),
            file = paste0(folder_path, "/all_split_reads_raw.rds"))
  }
  
  
  all_split_reads_raw_tidy %>% 
    distinct(junID) %>%
    nrow()
  
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
  ## Data check
  all_split_reads_raw_tidy_gr %>% head()
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>%
    nrow() %>%
    print()
  all_split_reads_raw_tidy %>% nrow() -
    all_split_reads_raw_tidy_gr %>% length()
  
  #######################################################
  ## Remove split reads overlapping the ENCODE backlist
  #######################################################
  
  print(paste0(Sys.time(), " - removing blacklist sequences..."))
  
  
  blacklist_path <- paste0(getwd(), "/dependencies/hg38-blacklist.v2.bed")
  all_split_reads_raw_tidy_gr <- remove_encode_blacklist_regions(GRdata = all_split_reads_raw_tidy_gr,
                                                                 blacklist_path = blacklist_path)
  
  ## Data check
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>%
    nrow() %>%
    print()
  
  
  #######################################################
  ## Anotate using 'dasper'
  #######################################################
  
  print(paste0(Sys.time(), " - annotating dasper..."))
  
  if ( file.exists(paste0(getwd(), "/dependencies/Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")) ) {
    gtf_path <- paste0(getwd(), "/dependencies/Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")
  } else {
    gtf_path <- paste0(getwd(), "/dependencies/Homo_sapiens.GRCh38.", gtf_version, ".gtf")
  }
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
  
  folder_path <- paste0(getwd(), "/database/v", gtf_version, "/")
  dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
  saveRDS(all_split_reads_details_w_symbol_reduced,
          file = paste0(folder_path, "/all_split_reads_gtex_recount3_", gtf_version, ".rds"))
  
  
  
  
  ################################################################################
  ## Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  print(paste0(Sys.time(), " - removing split reads not classified as 'annotated', 'novel_donor' or 'novel_acceptor'..."))
  
  if ( !exists("all_split_reads_details_w_symbol_reduced") ) {
    all_split_reads_details_w_symbol_reduced <- readRDS(file = paste0(folder_path, "/all_split_reads_gtex_recount3_", gtf_version, ".rds"))
  }
  
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
    dplyr::filter(width >= 25)
  
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
          file = paste0(folder_path, "/all_ambiguous_jxn_gtex_recount3_", gtf_version, ".rds"))
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- all_split_reads_details_w_symbol_reduced_keep_gr %>%
    dplyr::filter(ambiguous == F) %>%
    dplyr::select(-ambiguous)
  
  all_split_reads_details_w_symbol_reduced_keep_gr
  
  saveRDS(object = all_split_reads_details_w_symbol_reduced_keep_gr,
          file = paste0(folder_path, "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))
  
  
  ## FREE UP SOME MEMORY
  rm(edb)
  rm(ambiguous_introns)
  rm(all_split_reads_raw_tidy_gr)
  rm(all_split_reads_details_w_symbol)
  rm(all_split_reads_details_w_symbol_reduced)
  rm(all_split_reads_details_w_symbol_reduced_keep_gr)
  rm(all_split_reads_details_w_symbol_reduced_keep)
  rm(all_split_reads_details_w_symbol_reduced_discard)
  gc()
}


tidy_recount3_data_per_tissue <- function(projects_used, 
                                          main_project,
                                          gtf_version) {
  
  folder_database <- paste0(getwd(), "/database/v", gtf_version, "/")
  
  print("Loading tidy split reads from recount3....")
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- 
    readRDS(file = paste0(folder_database, "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))
  
  all_projects_used <- NULL
  
  # ## Generate the raw file
  for (project_id in projects_used) {
    
    # project_id <- projects_used[1]
    # project_id <- projects_used[6]
    # project_id <- projects_used[31]
    
    folder_results <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", 
                             main_project, "/base_data/")
    dir.create(file.path(folder_results), recursive = TRUE, showWarnings = T)
    
    print(paste0(Sys.time(), " - getting junction data from recount3 - '", project_id, "' tissue..."))
    
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      type = "jxn"
    )
    
    #################################
    ## GET METADATA
    #################################
    
    metadata.info <- rse %>% 
      SummarizedExperiment::colData() %>%
      as_tibble() %>%
      filter(gtex.smafrze != "EXCLUDE") %>%
      filter(gtex.smrin >= 6.0)

    saveRDS(object = metadata.info, 
            file = paste0(folder_results, "/", project_id, "_samples_metadata.rds"))
    
    
    #################################
    ## GET SPLIT READS AND COUNTS  
    ## PER TISSUE
    #################################
    
    clusters_ID <- rse$gtex.smtsd %>% unique()
    clusters_used <- NULL
    gc()
    
    for (cluster_id in clusters_ID) {
      
      # cluster_id <- clusters_ID[1]
      print(paste0(Sys.time(), " - filtering junction data by cluster - '", cluster_id, "' tissue..."))
      
      ################
      ## Get clusters
      ################
      
      minimum_samples <- 70
      cluster_samples <- metadata.info %>% 
        filter(gtex.smtsd == cluster_id) %>%
        distinct(external_id) %>% 
        pull()
        
      
      saveRDS(object = cluster_samples, 
              file = paste0(folder_results, "/", project_id, "_", cluster_id, "_samples_used.rds"))
      
      
      if ( (cluster_samples %>% length() >= minimum_samples ) ) {
        
        clusters_used <- c(clusters_used, cluster_id)
        all_projects_used <- c(all_projects_used, project_id)
        
        ################
        ## Get counts
        ################
        
        print(paste0(Sys.time(), " - getting split read counts..."))
        
        counts <- (rse[, rse$external_id %in% cluster_samples] %>% SummarizedExperiment::assays())[[1]]
        counts <- counts[(rownames(counts) %in%
                            all_split_reads_details_w_symbol_reduced_keep_gr$junID),]
        counts <- counts[rowSums(counts) > 0, ]
        counts <- counts %>% as.matrix()
        counts <- counts %>% as_tibble(rownames = "junID")
        
        # counts %>% head
        print(object.size(counts), units = "GB")
        saveRDS(object = counts,
                file = paste0(folder_results, "/", project_id, "_", cluster_id, "_split_read_counts.rds"))
        gc()
        
        #######################
        ## Get split reads ID
        #######################
        
        print(paste0(Sys.time(), " - getting split read IDs..."))
        all_split_reads <- counts %>%
          dplyr::select(junID) %>%
          data.table::as.data.table() %>%
          inner_join(y = all_split_reads_details_w_symbol_reduced_keep_gr,
                     by = "junID")
        
        #######################
        ## qc and save results
        #######################
        
        if (any(all_split_reads$width < 25) |
            any(str_detect(string = all_split_reads$chr, pattern = "random")) |
            any(str_detect(string = str_to_lower(all_split_reads$chr), pattern = "u")) |
            any(!(all_split_reads$type %in% c("annotated", "novel_acceptor", "novel_donor")))) {
          print("ERROR! The split reads do not meet the minimum filtering criteria.")
          break;
        }
        # if (setdiff(counts$i, all_split_reads$index) %>% length() > 0) {
        #   print("ERROR! The split reads are not identical.")
        #   break;
        # }
        
        ## Check how much memory this object uses
        print(object.size(all_split_reads %>% data.table::as.data.table()), units = "GB")
        
        ## Save the object
        saveRDS(object = all_split_reads %>% data.table::as.data.table(),
                file = paste0(folder_results, "/", project_id, "_", cluster_id, "_all_split_reads.rds"))
        
        
        ## Free up some memory
        rm(all_split_reads)
        rm(counts)
        gc()
      }
      
    } 
    
    if ( !is.null(clusters_used) ) {
      saveRDS(object = clusters_used, 
              file = paste0(folder_results, "/", project_id, "_clusters_used.rds"))
    }
    
    rm(rse)
    gc()
  }
  
  
  saveRDS(object = all_projects_used, 
          file = paste0(folder_results, "/all_projects_used.rds"))
  
}


junction_pairing <- function(projects_used, 
                             gtf_version,
                             main_project) {
  
  
  
  for (project_id in projects_used) {
    
    # project_id <- projects_used[1]
    
    folder_root <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", 
                          main_project, "/")
    folder_path <- paste0(folder_root, "/base_data/")
    
    print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
    
    ## Load clusters
    metadata.info <- readRDS(file = paste0(folder_path, "/", project_id, "_samples_metadata.rds"))
    clusters_ID <- metadata.info$gtex.smtsd %>% unique()
    
    
    for (cluster_id in clusters_ID) {
      
      # cluster_id <- clusters_ID[1]
      # cluster <- clusters_ID[2]
      
      print(paste0(Sys.time(), " - loading '", cluster_id, "' source data ..."))
      
      ############################################
      ## LOAD DATA FOR THE CURRENT PROJECT
      ############################################
      
      ## Load samples
      samples_used <- readRDS(file = paste0(folder_path, "/", 
                                            project_id, "_", cluster_id, "_samples_used.rds"))
      
      ## IntroVerse project considers all samples regardless of their RIN numbers. It also includes all
      ## tissues, regardless of their number of samples. However, this is configurable:
      if ( main_project == "splicing" ) {
        minimum_samples <- 70
      }
      
      
      if ( samples_used %>% length() >= minimum_samples ) {
        
        
        folder_name <- paste0(folder_root, "/results/", cluster_id, "/")
        dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)
        
        ## Load split read data
        all_split_reads_details <- readRDS(file = paste0(folder_path, "/", project_id, "_", cluster_id, 
                                                         "_all_split_reads.rds")) %>% as_tibble()
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(folder_path, "/", project_id, "_", cluster_id, "_",
                                                   "split_read_counts.rds")) %>% as_tibble()
        
        if ( !identical(all_split_reads_details$junID, split_read_counts$junID) ) {
          print("ERROR! The number of junctions considered is not correct.")
          break;
        }
        
        ############################################
        ## DISTANCES SUITE OF FUNCTIONS
        ############################################
        
        get_distances(cluster = cluster_id,
                      samples = samples_used,
                      split_read_counts = split_read_counts,
                      all_split_reads_details = all_split_reads_details,
                      folder_name)
        gc()
        
        
        extract_distances(cluster = cluster_id,
                          samples = samples_used,
                          folder_name = folder_name)
        gc()
        
        
        get_never_misspliced(cluster = cluster_id,
                             samples = samples_used,
                             split_read_counts = split_read_counts,
                             all_split_reads_details = all_split_reads_details,
                             folder_name = folder_name)
        
        rm(all_split_reads_details)
        rm(split_read_counts)
        gc()
        
      }
    }
  }
}


tidy_data_pior_sql <- function (projects_used, 
                                gtf_version,
                                all_clusters = NULL,
                                main_project) {
  
  
  print(paste0(Sys.time(), " - loading base data ..."))
  
  ## Load base recount3 object containing the GTEx v8 split reads passing the QC criteria
  all_split_reads_details_w_symbol_reduced_keep_gr <- readRDS(file = paste0(getwd(), "/database/v", gtf_version, 
                                                                            "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))
  all_split_reads_details_w_symbol_reduced_keep_gr %>% nrow()
  
  ############################################
  ## Discard all junctions from 
  ## EXCLUDE ME samples
  ############################################
  
  all_split_reads_details_all_tissues <- readRDS(file = paste0(getwd(), "/database/v", gtf_version, "/", main_project, 
                                                               "/all_split_reads_details_all_tissues.rds"))
  
  all_split_reads_details_all_tissues %>% nrow()
  all_split_reads_details_all_tissues %>% head()
  
  ## This should be zero
  if ( setdiff(all_split_reads_details_all_tissues$junID, 
               all_split_reads_details_w_symbol_reduced_keep_gr$junID) %>% 
       length() > 0) {
    print("ERROR! Some of the annotated split reads are not found within the base tidied recount3 object.")
    break;
  }
  
  
  # ## These are the junctions from EXCLUDE ME samples
  setdiff(all_split_reads_details_w_symbol_reduced_keep_gr$junID, 
          all_split_reads_details_all_tissues$junID) %>% unique %>% length()
  
  # ## Stats
  # all_split_reads_details %>% as.data.table() %>% dplyr::count(type)
  # all_split_reads_details %>% distinct(junID) %>% nrow()
  
  ############################################
  ## QC
  ############################################ 
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_w_symbol_reduced_keep_gr$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_w_symbol_reduced_keep_gr[ind, "junID"] <- str_replace(string = all_split_reads_details_w_symbol_reduced_keep_gr[ind, "junID"]$junID, 
                                                                                  pattern = "\\*", 
                                                                                  replacement = all_split_reads_details_w_symbol_reduced_keep_gr[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_w_symbol_reduced_keep_gr$junID, pattern = "\\*")) %>% print()
  }
  
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_all_tissues$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_all_tissues[ind, "junID"] <- str_replace(string = all_split_reads_details_all_tissues[ind, "junID"]$junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = all_split_reads_details_all_tissues[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_all_tissues$junID, pattern = "\\*")) %>% print()
  }
  
  
  ##########################################
  ## Load all distances pairings
  ##########################################
  
  df_all_distances_pairings_raw <- readRDS(file = paste0(getwd(), "/database/v", gtf_version, "/",
                                                         main_project, "/all_distances_pairings_all_tissues.rds"))
  
  ## QC
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = df_all_distances_pairings_raw$ref_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_distances_pairings_raw[ind, "ref_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "ref_junID"]$ref_junID, 
                                                                   pattern = "\\*", 
                                                                   replacement = df_all_distances_pairings_raw[ind, "ref_strand"]$ref_strand %>% as.character())
    any(str_detect(df_all_distances_pairings_raw$ref_junID, pattern = "\\*")) %>% print()
  }
  ## Remove potential * in the junID of the novel junctions
  ind <- which(str_detect(string = df_all_distances_pairings_raw$novel_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_distances_pairings_raw[ind, "novel_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "novel_junID"]$novel_junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = df_all_distances_pairings_raw[ind, "novel_strand"]$novel_strand  %>% as.character())
    any(str_detect(df_all_distances_pairings_raw$novel_junID, pattern = "\\*")) %>% print()
  }
  
  
  ##########################################
  ## Get never mis-spliced
  #########################################
  
  df_never_misspliced <- get_intron_never_misspliced(projects_used = projects_used,
                                                     all_clusters = all_clusters,
                                                     main_project = main_project)
  
  if ( any(str_detect(df_never_misspliced$ref_junID, pattern = "\\*")) ) {
    
    print("ERROR! Some junctions contain a '*' in their IDs")
    break;
  }
  
  ## Remove the introns paired with novel junctions (i.e. mis-spliced)
  df_never_misspliced_tidy <- df_never_misspliced %>%
    dplyr::filter(!(ref_junID %in% df_all_distances_pairings_raw$ref_junID)) %>%
    as_tibble()
  
  df_never_misspliced_tidy %>% distinct(ref_junID) %>% as_tibble()
  
  
  ############################################
  ## GET all not paired
  ############################################ 
  
  # df_all_novel_raw_tidy
  df_not_paired <- all_split_reads_details_all_tissues %>%
    as.data.table() %>%
    dplyr::filter(!(junID %in% c(df_all_distances_pairings_raw$ref_junID,
                                 df_all_distances_pairings_raw$novel_junID)))
  
  ## These are all the non-paired, including the never mis-spliced. Thus, GTEx 'splicing' data: [768,646 - 56,090 = 730125]
  df_not_paired %>% 
    distinct(junID)
  
  ## Check them by category
  df_not_paired %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  
  if (any(str_detect(df_never_misspliced_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(df_not_paired$junID, pattern = "\\*"))) {
    print("ERROR!")
  }
  
  ## All never mis-spliced should be categorised as not paired.
  ## Thus, this should be zero
  setdiff(df_never_misspliced_tidy$ref_junID, 
          df_not_paired %>% dplyr::filter(type == "annotated") %>% pull(junID))
  
  df_not_paired_tidy <- df_not_paired %>%
    dplyr::filter(!(junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    distinct(junID, .keep_all = T) %>%
    as_tibble()
  
  df_not_paired_tidy %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  df_not_paired_tidy %>% distinct(junID)
  
  
  ## This should be zero
  if ( intersect(df_not_paired_tidy$junID, 
                 df_all_distances_pairings_raw$novel_junID) %>% length() > 0 ) {
    print("ERROR!")
  }
  
  
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
  
  
  ## 1. Obtain the ambiguous junctions
  
  df_ambiguous_novel <- df_all_distances_pairings_raw %>%
    dplyr::filter(!(novel_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    dplyr::group_by(novel_junID) %>%
    mutate(distances_sd = distance %>% sd()) %>%
    dplyr::filter(distances_sd > 0)
  
  
  df_ambiguous_novel %>% ungroup() %>% distinct(novel_junID)
  df_ambiguous_novel %>% ungroup() %>% distinct(ref_junID)
  
  
  ## 2. Remove ambiguous junctions
  
  df_all_distances_pairings_raw_tidy <- df_all_distances_pairings_raw %>%
    dplyr::filter(!(novel_junID %in% df_ambiguous_novel$novel_junID)) %>%
    as.data.table() %>%
    distinct(novel_junID, ref_junID, .keep_all = T) %>%
    mutate(ref_strand = ref_strand %>% as.character(),
           novel_strand = novel_strand %>% as.character()) 
  
  
  ## Some numbers for the graph
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
  
  (df_all_distances_pairings_raw_tidy %>%
      dplyr::distinct(novel_junID) %>% 
      nrow()) + 
    (df_all_distances_pairings_raw_tidy %>%
       dplyr::distinct(ref_junID) %>% 
       nrow()) + 
    (df_never_misspliced_tidy %>% 
       distinct(ref_junID) %>% 
       nrow())
  
  
  ##############################################################################
  ## SAVE FINAL OBJECT
  ##############################################################################
  
  ## 1. DISTANCES PAIRINGS
  if (any(str_detect(string = df_all_distances_pairings_raw_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_all_distances_pairings_raw_tidy$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  df_all_distances_pairings_raw_tidy <- df_all_distances_pairings_raw_tidy %>%
    inner_join(y = all_split_reads_details_all_tissues %>% dplyr::select(junID, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  saveRDS(object = df_all_distances_pairings_raw_tidy,
          file = paste0(getwd(), "/database/v", gtf_version, "/",
                        main_project, "/all_distances_correct_pairings_all_tissues.rds"))
  
  
  ## 2. NEVER MIS-SPLICED
  if (any(str_detect(string = df_never_misspliced_tidy$ref_junID, pattern = "\\*")) ) {
    print("ERROR! Some NEVER MIS-SPLICED junctions still have a '*' in their IDs!")
  }
  df_never_misspliced_tidy <- df_never_misspliced_tidy %>%
    inner_join(y = all_split_reads_details_all_tissues %>% 
                 dplyr::select(junID, seqnames, start, end, width, strand, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  
  saveRDS(object = df_never_misspliced_tidy,
          file = paste0(getwd(), "/database/v", gtf_version, "/",
                        main_project, "/all_nevermisspliced_introns_all_tissues.rds"))
  
  
  ## 3. AMBIGUOUS JUNCTIONS
  if (any(str_detect(string = df_ambiguous_novel$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_ambiguous_novel$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  saveRDS(df_ambiguous_novel,
          file = paste0(getwd(), "/database/v", gtf_version, "/", main_project, 
                        "/all_distances_ambiguous_pairings_all_tissues.rds"))
}


sql_database_generation <- function(database_path,
                                    projects_used, 
                                    main_project,
                                    gtf_version,
                                    remove_all = NULL) {
  
  
  if ( !is.null(remove_all) ) {
    
    remove_tables(database_path, remove_all)
  } 
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  tables <- DBI::dbListTables(conn = con)
  
  
  if ( !any(tables == 'master') ) {
    create_metadata_table(database_path,
                          main_project = main_project,
                          gtf_version = gtf_version,
                          all_projects = projects_used)
  }
  
  
  if ( ! any(tables %in% c('intron', 'novel', 'gene', 'transcript')) ) {
    create_master_tables(database_path,
                         main_project = main_project,
                         gtf_version = gtf_version)
  }
  
  print(paste0(Sys.time(), " - creating cluster tables ..."))
  
  rm(CNC_CDTS_CONS_gr)
  rm(hg_mane_transcripts)
  rm(hg_MANE_tidy)
  rm(hg_MANE)
  rm(hg38_transcripts)
  rm(hg38)
  gc()
  
  create_cluster_tables(database_path = database_path,
                        gtf_version = gtf_version,
                        all_projects = projects_used,
                        main_project = main_project)
  
}

#####################################
## CALLS - PREPARE RECOUNT3 DATA
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_versions <- c(105)

## This is the name of the project producing the database
main_project <- "splicing"

## Can be checked here: https://jhubiostatistics.shinyapps.io/recount3-study-explorer/
all_projects <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLADDER",         "BLOOD",           "BLOOD_VESSEL",    "BONE_MARROW",    
                   "BRAIN",           "BREAST",          "CERVIX_UTERI",    "COLON",           "ESOPHAGUS",       "FALLOPIAN_TUBE", 
                   "HEART",           "KIDNEY",          "LIVER",           "LUNG",            "MUSCLE",          "NERVE",          
                   "OVARY",           "PANCREAS",        "PITUITARY",       "PROSTATE",        "SALIVARY_GLAND",  "SKIN",           
                   "SMALL_INTESTINE", "SPLEEN",          "STOMACH",         "TESTIS",          "THYROID",         "UTERUS",         
                   "VAGINA"  )


splicing_projects <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLOOD",           "BLOOD_VESSEL",    "BONE_MARROW",    
                        "BRAIN",           "COLON",           "ESOPHAGUS",       "HEART",           "KIDNEY",          
                        "LIVER",           "LUNG",            "MUSCLE",          "NERVE",           "PANCREAS",        
                        "PITUITARY",       "SALIVARY_GLAND",  "SKIN",            "SMALL_INTESTINE", "SPLEEN",          
                        "STOMACH",         "THYROID")

for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  # init_recount3_gtex_data(projects_used = all_projects,
  #                         gtf_version = gtf_version)
  
  
  # tidy_recount3_data_per_tissue(projects_used = splicing_projects,
  #                               main_project,
  #                               gtf_version = gtf_version)
  
  
  
  
  
  # junction_pairing(projects_used = splicing_projects[-c(1:2)],
  #                  main_project,
  #                  gtf_version = gtf_version)
  
  
  # get_all_annotated_split_reads(projects_used = splicing_projects,
  #                               gtf_version = gtf_version,
  #                               main_project = main_project)
  # 
  # 
  # get_all_raw_distances_pairings(projects_used = splicing_projects,
  #                                gtf_version = gtf_version,
  #                                main_project = main_project)
  
  
  # tidy_data_pior_sql(projects_used = splicing_projects,
  #                    gtf_version = gtf_version,
  #                    main_project = main_project)
  
  
  # generate_transcript_biotype_percentage(projects_used = splicing_projects,
  #                                        homo_sapiens_v105_path = paste0(dependencies_folder, 
  #                                                                        "/Homo_sapiens.GRCh38.105.chr.gtf"),
  #                                        main_project,
  #                                        gtf_version = gtf_version)
  
  
  all_final_projects_used <- readRDS(file = paste0(getwd(),"/results/all_final_projects_used.rds"))
  
  # generate_recount3_tpm(projects_used = all_final_projects_used,
  #                       gtf_version = gtf_version,
  #                       main_project = main_project)
  
  database_folder <- paste0(getwd(), "/database/v", gtf_version, "/", main_project)
  dir.create(file.path(database_folder), recursive = TRUE, showWarnings = T)
  database_path <- paste0(database_folder,  "/", main_project, ".sqlite")

  sql_database_generation(database_path = database_path,
                          projects_used = all_final_projects_used,
                          main_project = main_project,
                          gtf_version = gtf_version,
                          remove_all = F)
  
  gc()
}

