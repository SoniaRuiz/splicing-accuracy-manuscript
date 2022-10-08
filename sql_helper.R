library(tidyverse)
library(data.table)
library(GenomicRanges)
library(DBI)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/helper.R")
# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline1_download_from_recount.R")
# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline2_annotate_from_recount.R")
# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_generation.R")


database_path <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/splicing-intronextra-v2.sqlite"

#################################################################

get_all_raw_distances_pairings <- function(database_path,
                                           main_project = "splicing") {
  
  
  gtf_version <-  105
  
  ## Connect to the DB
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  # GET INFO FROM MASTER
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  

  
  ## LOOP THROUGH PROJECTS
  df_all_distances_pairings_raw <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          db, "/", main_project, "_project/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(clusters, function(cluster) {
      
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      
      ## Load samples
      samples <- readRDS(file = paste0(base_folder, "/raw_data/", db, "_", cluster,  "_samples_used.rds"))
      
      if (samples %>% length() > 0) {
        
        folder_name <- paste0(base_folder, "/results/", cluster, "/distances/v", gtf_version, "/")
        
        ## Obtain the distances across all samples
        df_all <- map_df(samples, function(sample) { 
          
          # sample <- samples[1]
          print(paste0(cluster, " - ", sample))
          file_name <- paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")
          
          
          if (file.exists(file_name)) {
            
            df <- readRDS(file = file_name)
            
            return(df)
          } 
          
        })
        
        saveRDS(object = df_all %>%
                  distinct(novel_junID, ref_junID, .keep_all = T) %>%
                  mutate(tissue = cluster),
                file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
        
        return(df_all %>%
                 distinct(novel_junID, ref_junID, .keep_all = T))
      } else {
        return(NULL)
      }
    })  
  })  
  
  
  
  
  saveRDS(object = df_all_distances_pairings_raw %>%
            distinct(novel_junID, ref_junID, .keep_all = T) %>% 
            as.data.table(),
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                        main_project, "_project_all_distances_pairings_raw.rds"))

}



#################################################################

get_mean_coverage <- function(split_read_counts,
                              samples,
                              junID) {
  
  split_read_counts_intron <- split_read_counts %>%
    dplyr::filter(junID %in% junID) %>%
    dplyr::select(junID, all_of(samples %>% as.character())) 
  
  split_read_counts_intron[,"n_individuals"] <- (matrixStats::rowCounts(split_read_counts_intron[, -c(1)] > 0, na.rm = T)) 
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron[,"sum_counts"] <- Matrix::rowSums(split_read_counts_intron[,-c(split_read_counts_intron %>% ncol(),1)], na.rm = T)
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron <- split_read_counts_intron[, c(1,(split_read_counts_intron %>% ncol() - 1),(split_read_counts_intron %>% ncol()))]
  
  
  if (any(split_read_counts_intron[, "n_individuals"] < 1)) {
    print("Error: some ref junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_intron %>% return()
}

################################################################

get_intron_never_misspliced <- function (database_path) {
  
  gtf_version <- 105
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbListTables(conn = con)
  
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  

  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  df_never <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[6]
    
    # print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(clusters, function(cluster) { 
      
      # cluster <- clusters[1]
      
      # print(paste0(Sys.time(), " --> ", cluster))
      if (file.exists(paste0(base_folder, "results/pipeline3/distances/", 
                             cluster, "/v", gtf_version, "/not-misspliced/", cluster, "_all_notmisspliced.rds"))) {
        df_introns_never <- readRDS(file = paste0(base_folder, "results/pipeline3/distances/", 
                                                  cluster, "/v", gtf_version, "/not-misspliced/", cluster, "_all_notmisspliced.rds")) %>% as_tibble()
        return(data.frame(ref_junID = df_introns_never$value))
      } else {
        return(NULL)
      }
      
    })
  })
  
  
  

  df_never %>%
    distinct(ref_junID) %>%
    return()
  
  
}




#################################################################

get_all_annotated_split_reads_from_use_me_samples <- function() {
  
  
  database_path <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/splicing-intronextra-v2.sqlite"
  
  getwd()
  setwd("~/PROJECTS/splicing-project-recount3/")
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  # all_split_reads_details_105 <- readRDS(file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_105_length_all_tissues.rds")
  
  #############################################
  ## GET ALL SPLIT READS FROM 'USE ME' SAMPLES
  #############################################
  
  all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
  
  all_split_reads_details_105 <- map_df(all_projects, function(project_id) {
      
    # project_id <- all_projects[1]
    # project_id <- "BONE_MARROW"
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
    
    
    jxn_qc <- map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[1]
      print(paste0(Sys.time(), " - loading '", cluster, "'  data ..."))
      
      if (file.exists(paste0(folder_root, "/results/base_data/", cluster, "/", cluster, "_annotated_SR_details_length_105.rds"))) {
        
        all_split_reads_details_105 <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/",
                                                             cluster, "_annotated_SR_details_length_105.rds"))
        
        
        
        ## Remove split reads annotated to multiple genes
        all_split_reads_details_tidy <- all_split_reads_details_105 %>%
          distinct(junID, .keep_all = T) %>% 
          rowwise() %>%
          mutate(ambiguous = ifelse(gene_id_junction %>% unlist() %>% length() > 1, T, F))
        all_split_reads_details_tidy <- all_split_reads_details_tidy %>%
          filter(ambiguous == F)
        
        saveRDS(all_split_reads_details_tidy %>% dplyr::select(-ambiguous),
                file =  paste0(folder_root, "/results/base_data/", cluster, "/",
                               cluster, "_annotated_SR_details_length_105.rds"))
        
        ## Print message
        print(paste0(all_split_reads_details_105 %>% nrow, " - ", all_split_reads_details_tidy %>% nrow()))
        
        all_split_reads_details_tidy %>%
          distinct(junID, .keep_all = T) %>% 
          as_tibble() %>%
          return()
        
      } else {
        return(NULL)
      }
      
    })
    if (jxn_qc %>% nrow() > 0 ) {
      jxn_qc %>%
        distinct(junID, .keep_all = T) %>% return()
    } else {
      return(NULL)
    }
    
  })
  
  all_split_reads_details_105 <- all_split_reads_details_105 %>%
    distinct(junID, .keep_all = T)
  
  saveRDS(object = all_split_reads_details_105,
          file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_105_length_all_tissues.rds")
  
    
  

  
}

# get_all_annotated_split_reads_from_use_me_samples()


##################################################################

init_data <- function(projects_used) {
  
  # project_id <- projects_used[6]
  
  
  for (project_id in projects_used) {
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
            file = paste0(folder_path, "/", "all_split_reads_raw.rds"))
  }
  
  
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  all_split_reads_raw <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
    folder_path <- paste0(folder_root, "/raw_data/")
    
    
    print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
    
    readRDS(file = paste0(folder_path, "/all_split_reads_raw.rds")) %>%
      return()
  })
  
  
  all_split_reads_raw_tidy <- all_split_reads_raw %>%
    mutate(junID_tidy = paste0(chr, ":", start, "-", end, ":", strand)) %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T)
  
  all_split_reads_raw_tidy %>% nrow()
  
  
  
  
  #######################################################################
  ## Remove split reads located in unplaced sequences in the chromosomes
  #######################################################################
  
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
    nrow()
  
  
  #######################################################
  ## Remove split reads overlapping the ENCODE backlist
  #######################################################
  
  source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline2_annotate_from_recount.R")
  blacklist_path <- "/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed"
  all_split_reads_raw_tidy_gr <- remove_encode_blacklist_regions(GRdata = all_split_reads_raw_tidy_gr,
                                                                 blacklist_path = blacklist_path)
  
  
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T) %>% 
    nrow()
  
  
  #######################################################
  ## Anotate using 'dasper'
  #######################################################
  
  gtf_path <- paste0("/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  
  all_split_reads_details_105_w_symbol <- dasper::junction_annot(junctions = all_split_reads_raw_tidy_gr %>% GRanges(), 
                                                                 ref = edb)
  
  all_split_reads_details_105_w_symbol <- all_split_reads_details_105_w_symbol %>% 
    as_tibble()
  
  all_split_reads_details_105_w_symbol_reduced <- all_split_reads_details_105_w_symbol %>% 
    dplyr::select(seqnames, start, end, strand, junID, gene_id = gene_id_junction, in_ref, type )
  
  saveRDS(all_split_reads_details_105_w_symbol_reduced,
          file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_details_105_w_symbol.rds")
  
  
  
  
  ################################################################################
  ## Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  all_split_reads_details_105_w_symbol_reduced <- readRDS(file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_details_105_w_symbol.rds")
  
  all_split_reads_details_105_w_symbol_reduced_discard <- all_split_reads_details_105_w_symbol_reduced %>%
    dplyr::filter(!(type %in% c("annotated", "novel_donor", "novel_acceptor"))) %>%
    as.data.table()
  all_split_reads_details_105_w_symbol_reduced_discard %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  all_split_reads_details_105_w_symbol_reduced_discard %>%
    distinct(junID) %>%
    nrow()
  
  
  all_split_reads_details_105_w_symbol_reduced_keep <- all_split_reads_details_105_w_symbol_reduced %>%
    as.data.table() %>%
    dplyr::filter(!(junID %in% all_split_reads_details_105_w_symbol_reduced_discard$junID))
  
  ############################################
  ## Discard all junctions shorter than 25bp
  ############################################
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep %>%
    GRanges() %>%
    as_tibble()
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr  %>% 
    dplyr::filter(width < 25) %>%
    nrow()
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep_gr  %>% 
    dplyr::filter(width >=25)
  
  
  ############################################################################
  ## Discard all introns assigned to multiple genes (i.e. ambiguous introns)
  ############################################################################
  
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep_gr %>%
    distinct(junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  ambiguous_introns <- all_split_reads_details_105_w_symbol_reduced_keep_gr %>%
    dplyr::filter(ambiguous == T)
  
  ambiguous_introns %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  
  saveRDS(object = ambiguous_introns,
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/all_ambiguous_introns.rds")
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep_gr %>%
    dplyr::filter(ambiguous == F) %>%
    dplyr::select(-ambiguous)
  
  saveRDS(object = all_split_reads_details_105_w_symbol_reduced_keep_gr,
          file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_details_105_w_symbol_reduced_keep.rds")
  
  
  ############################################
  ## Discard all junctions from EXCLUDE ME samples
  ############################################
  
  
  all_split_reads_details_105 <- readRDS(file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_105_length_all_tissues.rds")
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- readRDS(file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_details_105_w_symbol_reduced_keep.rds")
  
  
  ## This should be zero
  setdiff(all_split_reads_details_105$junID, all_split_reads_details_105_w_symbol_reduced_keep_gr$junID)
  
  ## These are the junctions from EXCLUDE ME samples
  setdiff(all_split_reads_details_105_w_symbol_reduced_keep_gr$junID, all_split_reads_details_105$junID) %>% unique %>% length()
  
  
  all_split_reads_details_105 %>%
    as.data.table() %>%
    dplyr::count(type)
  
  
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
  ind <- which(str_detect(string = all_split_reads_details_105$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_105[ind, "junID"] <- str_replace(string = all_split_reads_details_105[ind, "junID"]$junID, 
                                                             pattern = "\\*", 
                                                             replacement = all_split_reads_details_105[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_105$junID, pattern = "\\*"))
  }
  
  
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

