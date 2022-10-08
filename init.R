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


get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}


database_path <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/splicing/splicing.sqlite"
projects_used <- readRDS(file = "~/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
clusters_used <- readRDS(file = "~/PROJECTS/splicing-project/splicing-recount3-projects/all_clusters_used.rds")
gtf_version <- 105

source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_database_generation.R")
source("/home/sruiz/PROJECTS/splicing-project-recount3/sql_helper.R")
source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_SQL_generation_extra.R")

#####################################
## FUNCTIONS - PREPARE RECOUNT3 DATA
#####################################

download_recount_split_reads <- function(projects_used, 
                                         main_project = "splicing") {
  
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- readRDS(file = "~/PROJECTS/splicing-project-recount3/database/introverse/all_split_reads_details_105_w_symbol_reduced_keep.rds")
  
  # ## Generate the raw file
  for (project_id in projects_used) {
    
    # project_id <- projects_used[1]
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          project_id, "/", main_project, "_project/")
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
    
    clusters_ID = rse$gtex.smtsd %>% unique()
    
    for (cluster_id in clusters_ID) {
      
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
        
        saveRDS(object = cluster_samples, 
                file = paste0(folder_path, "/", project_id, "_", cluster_id, "_samples_used.rds"))
        
        local_rse <- rse[, rse$external_id %in% cluster_samples]
        
        
        print(paste0(Sys.time(), " - getting split read counts..."))
        counts <- (local_rse %>% SummarizedExperiment::assays())[[1]]
        counts <- counts[(rownames(counts) %in% all_split_reads_details_105_w_symbol_reduced_keep_gr$junID),]
        counts <- counts[rowSums(counts) > 0, ]
        counts <- counts %>% as.matrix()
        counts %>% rownames()
        counts %>% nrow()
        saveRDS(object = counts,
                file = paste0(folder_path, "/", project_id, "_", cluster_id, "_split_read_counts_sample_tidy.rds"))
        
        
        print(paste0(Sys.time(), " - getting split read IDs..."))
        all_split_reads <- data.frame(junID = local_rse %>% rownames())
        all_split_reads %>% nrow()
        
        all_split_reads_tidy <- all_split_reads %>%
          data.table::as.data.table() %>%
          inner_join(y = all_split_reads_details_105_w_symbol_reduced_keep_gr,
                     by = "junID")
        
        
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
  }
  
  ## Save clusters_used, as it will be part of the database
  clusters_used <- data.frame(project_id = as.character(),
                              cluster_id = as.character(),
                              sample_id = as.character())
  
  for (project_id in projects_used) {
    
    # project_id <- projects_used[10]
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          project_id, "/", main_project, "_project/")
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
  
  
  saveRDS(clusters_used,
          file =  "/home/sruiz/PROJECTS/splicing-project-recount3/database/splicing/metadata.rds")
}

prepare_recount_split_reads <- function(projects_used) {
  
  # project_id <- projects_used[6]
  
  ##########################################################
  ## Load original recount3 data pre-filtered for IntroVerse
  ##########################################################
  
  intersect(all_split_reads_splicing_project$junID,
            all_split_reads_details_105_w_symbol_reduced_keep_gr$junID) %>% 
    unique() %>% 
    length()
  
  ##########################################################
  ## Read all the split reads filtered for the splicing project
  ##########################################################
  
  all_split_reads_splicing_project <- map_df(projects_used[1:5], function(project_id) {
    
    # project_id <- projects_used[1]
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
    folder_path <- paste0(folder_root, "/raw_data/")
    
    
    metadata.info <- readRDS(file = paste0(folder_path, "/samples_metadata.rds"))
    clusters <- metadata.info %>% 
      as_tibble() %>%
      filter(gtex.smrin >= 6.0,
             gtex.smafrze != "EXCLUDE") %>%
      distinct(gtex.smtsd) %>% 
      pull()
    
    map_df(clusters, function(cluster_id) {
      
      print(paste0(Sys.time(), " - getting data from '", cluster_id, "' tissue..."))
      
      readRDS(file = paste0(folder_path, "/", project_id, "_", cluster_id, "_all_split_reads_sample_tidy.rds")) %>%
        return()
      
    })
    
  })
  
  
  all_split_reads_splicing_project %>% nrow()
  
  all_split_reads_splicing_project_tidy <- all_split_reads_splicing_project %>%
    distinct(junID, .keep_all = T) %>%
    mutate(junID_tidy = paste0(chr, ":", start, "-", end, ":", strand)) %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T)
  
  all_split_reads_splicing_project_tidy %>% nrow()
  
  
  ############################################
  ## QC
  ############################################ 
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_105_w_symbol_reduced_keep_gr$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_105_w_symbol_reduced_keep_gr[ind, "junID"] <- str_replace(string = all_split_reads_details_105_w_symbol_reduced_keep_gr[ind, "junID"]$junID, 
                                                                                      pattern = "\\*", 
                                                                                      replacement = all_split_reads_details_105_w_symbol_reduced_keep_gr[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_105_w_symbol_reduced_keep_gr$junID, pattern = "\\*"))
  }
  
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_splicing_project_tidy$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_splicing_project_tidy[ind, "junID"] <- str_replace(string = all_split_reads_splicing_project_tidy[ind, "junID"]$junID, 
                                                             pattern = "\\*", 
                                                             replacement = all_split_reads_splicing_project_tidy[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_splicing_project_tidy$junID, pattern = "\\*"))
  }
  
  
  if (!(identical(all_split_reads_splicing_project_tidy$junID, 
                  all_split_reads_splicing_project_tidy$junID_tidy))) {
    print("ERROR! some IDs are not equal.")
    break;
  }
  
  
  ## This should be zero
  intersect(all_split_reads_splicing_project_tidy$junID,
            all_split_reads_details_105_w_symbol_reduced_keep_gr$junID) %>% 
    unique() %>% 
    length()
  
  
  ## This should correspond to the number of junctions discarded as they did not pass the splicing-project criteria
  setdiff(all_split_reads_details_105_w_symbol_reduced_keep_gr$junID, 
          all_split_reads_splicing_project_tidy$junID)
  
  
  ############################################
  ## GET all not paired
  ############################################ 
  ## The object is now ready - distances should start
  
  ############################################
  ## GET all not paired
  ############################################ 
  
  df_all_distances_pairings_raw <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/df_all_distances_pairings_raw.rds")
  
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
  
  ## Get all non-paired
  
  # df_all_novel_raw_tidy
  df_not_paired <- all_split_reads_details_105 %>%
    as.data.table() %>%
    dplyr::filter(!(junID %in% c(df_all_distances_pairings_raw$ref_junID,
                                 df_all_distances_pairings_raw$novel_junID)))
  
  ## These are all the non-paired, including the never mis-spliced. Thus, [768,646 - 38,521 = 730125]
  df_not_paired %>% distinct(junID)
  df_not_paired %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  
  ##########################################
  ## Get never mis-spliced
  ##########################################
  
  df_never_misspliced <- get_intron_never_misspliced(database_path = database_path)
  
  ## Remove the introns paired with novel junctions
  df_never_misspliced_tidy <- df_never_misspliced %>%
    dplyr::filter(!(ref_junID %in% df_all_distances_pairings_raw$ref_junID)) %>%
    as_tibble()
  
  df_never_misspliced_tidy %>% distinct(ref_junID) %>% as_tibble()
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
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/df_all_tissues_raw_distances_ambiguous.rds")
  
  
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
  
  ##############################################################################
  ## SAVE FINAL OBJECT
  ##############################################################################
  
  ## 1. DISTANCES PAIRINGS
  
  if (any(str_detect(string = df_all_distances_pairings_raw_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_all_distances_pairings_raw_tidy$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  
  df_all_distances_pairings_raw_tidy <- df_all_distances_pairings_raw_tidy %>%
    inner_join(y = all_split_reads_details_105 %>% dplyr::select(junID, gene_id = gene_id_junction, tx_id_junction),
               by = c("ref_junID" = "junID"))
  
  saveRDS(object = df_all_distances_pairings_raw_tidy,
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/all_paired_intron_novel_tidy.rds")
  
  
  ## 2. NEVER MIS-SPLICED
  
  if (any(str_detect(string = df_never_misspliced_tidy$ref_junID, pattern = "\\*")) ) {
    print("ERROR! Some NEVER MIS-SPLICED junctions still have a * in their IDs!")
  }
  df_never_misspliced_tidy <- df_never_misspliced_tidy %>%
    inner_join(y = all_split_reads_details_105 %>% 
                 dplyr::select(junID, 
                               seqnames,start,end, width, strand,
                               gene_id = gene_id_junction, tx_id_junction),
               by = c("ref_junID" = "junID"))
  saveRDS(object = df_never_misspliced_tidy,
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/df_all_nevermisspliced_introns.rds")
  
  
  ## AMBIGUOUS JUNCTIONS
  if (any(str_detect(string = df_ambiguous_novel$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_ambiguous_novel$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  saveRDS(df_ambiguous_novel,
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/df_all_tissues_raw_distances_ambiguous.rds")
}

junction_pairing <- function(projects_used, 
                             main_project = "splicing") {
  
  
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
      }
    }
  }
}

sql_database_generation <- function(database_path,
                                    projects_used, 
                                    main_project = "splicing") {
  
  create_metadata_table(database_path,
                        main_project = "splicing",
                        SRA_projects = projects_used)
  
  get_all_raw_distances_pairings(database_path,
                                 main_project = "splicing")
 
  ## Tengo que coger todos los files all_split_reads_annotated_105 de cada tissue
  ## y dejar únicamente los IDs únicos
  ## eso compararlo con el fichero común de 6,105... split reads
  ## la differencia son los split reas que pertenecen a EXCLUDE ME samples,
  ## sex-related tissues, and tissues with less than 70 samples
  ## 2. comparo el fichero 'all_split_reads_annotated_105' con los paired distances,
  ## y la diferencia son los junctions que no se han emparejado.
  ## a estos últimos tengo que descontarle los introns never mis-spliced
  ## 3. después elimino los novel junctions que sean ambiguous
  ## y los introns que únicamente los estén parenting a ellos.
  ## 4. Al final, me tendrá que salir la cuenta de los introns and junctions 
  ## almacenados en spliceverse
  
  ## cuando pueda reproducir los análisis que tengo dise;ados, entonces
  ## hago el poster para el ASHG
  
  
  create_master_tables(database_path, 
                       file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                     main_project, "_project_all_distances_pairings_raw.rds"))
}

#####################################
## CALLS - PREPARE RECOUNT3 DATA
#####################################

#download_recount_split_reads(projects_used[9])
#junction_pairing(projects_used[-c(1:8)])
sql_database_generation(database_path,
                        projects_used, 
                        main_project = "splicing")

