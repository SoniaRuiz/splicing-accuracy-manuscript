library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/init.R")

setwd(normalizePath("."))

dependencies_folder <- paste0(getwd(), "/dependencies/")

#source(paste0(getwd(), "/database_junction_pairing.R"))
source(paste0(getwd(), "/database_SQL_helper.R"))
source(paste0(getwd(), "/database_SQL_generation.R"))



#####################################
## FUNCTIONS - PREPARE RECOUNT3 DATA
#####################################

download_recount3_data <- function (recount3_project_IDs,
                                    gtf_version) {
  
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  if ( file.exists(paste0(getwd(), "/database/all_split_reads_raw.rds")) ) {
    
    print(paste0(Sys.time(), " - loading 'all_split_reads_raw.rds' file..."))
    all_split_reads_raw_tidy <- readRDS(file = paste0(getwd(), "/database/all_split_reads_raw.rds"))
    
  } else {
    
    all_split_reads_raw <- map_df(recount3_project_IDs, function(project_id) {
      
      # project_id <- recount3_project_IDs[31]
      print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
      
      folder_root <- paste0(getwd(), "/results/", project_id, "/")
      dir.create(file.path(folder_root), recursive = TRUE, showWarnings = T)
      
      ## Build the RSE object from recount3
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
  
  
  ############################################
  ## QC STEP 1. 
  ## Discard all junctions shorter than 25bp
  ############################################
  
  
  all_split_reads_raw_tidy %>% 
    distinct(junID) %>%
    nrow()
  
  all_split_reads_raw_tidy <- all_split_reads_raw_tidy %>%
    filter(width >= 25)
  
  all_split_reads_raw_tidy %>% 
    distinct(junID) %>%
    nrow()
  
  
  ############################################
  ## QC STEP 2. 
  ## Remove split reads located in
  ## unplaced sequences in the chromosomes
  ############################################
  
  
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
  
  
  
  ############################################
  ## QC STEP 3. 
  ## Remove split reads overlapping 
  ## the ENCODE backlist
  ############################################
  
  
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
  
  
  
  ############################################
  ## QC STEP 4. 
  ## Anotate split reads using 'dasper'
  ############################################
  
  
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
  
  
  ############################################
  ## QC STEP 5. 
  ## Discard all junctions that are not 
  ## annotated, novel donor or novel acceptor
  ############################################
  
  
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
  ## QC Step 6. 
  ## Discard all introns assigned to
  ## multiple genes (i.e. ambiguous introns)
  ############################################
  
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


tidy_recount3_data <- function(recount3_project_IDs, 
                               main_project,
                               gtf_version) {
  
  folder_database <- paste0(getwd(), "/database/v", gtf_version, "/")
  
  print("Loading tidy split reads from recount3....")
  
  all_split_reads_details_w_symbol_reduced_keep_gr <- 
    readRDS(file = paste0(folder_database, "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))
  
  all_projects_used <- NULL
  
  # ## Generate the raw file
  for (project_id in recount3_project_IDs) {
    
    # project_id <- recount3_project_IDs[1]
    # project_id <- recount3_project_IDs[6]
    # project_id <- recount3_project_IDs[31]
    
    folder_results <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", 
                             main_project, "/base_data/unique/")
    dir.create(file.path(folder_results), recursive = TRUE, showWarnings = T)
    
    print(paste0(Sys.time(), " - getting junction data from recount3 - '", project_id, "' tissue..."))
    
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      type = "jxn"
    )
    gc()
    
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
        
        ## At least two supportive reads per junction
        counts <- counts[rowSums(counts) >= 2, ]
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


junction_pairing <- function(recount3_project_IDs, 
                             gtf_version,
                             main_project) {
  
  
  
  for (project_id in recount3_project_IDs) {
    
    # project_id <- recount3_project_IDs[1]
    
    folder_root <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", 
                          main_project, "/")
    folder_path <- paste0(folder_root, "/base_data/")
    
    print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
    
    ## Load clusters
    
    if ( file.exists(paste0(folder_path, "/", project_id, "_samples_metadata.rds")) ) {
      
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
          
          ## Only split reads with at least 2 supportive reads
          split_read_counts_tidy <- split_read_counts %>%
            mutate(total=rowSums(select_if(., is.numeric))) %>% 
            filter(total >= 2) %>%
            dplyr::select(-total)
          
  
          all_split_reads_details_tidy <- all_split_reads_details %>%
            filter(junID %in% split_read_counts_tidy$junID)
            
          
          if ( !identical(all_split_reads_details_tidy$junID, split_read_counts_tidy$junID) ) {
            print("ERROR! The number of junctions considered is not correct.")
            break;
          }
          
          ############################################
          ## DISTANCES SUITE OF FUNCTIONS
          ############################################
          
          get_distances(cluster = cluster_id,
                        samples = samples_used,
                        split_read_counts = split_read_counts_tidy,
                        all_split_reads_details = all_split_reads_details_tidy,
                        folder_name)
          gc()
          
          
          extract_distances(cluster = cluster_id,
                            samples = samples_used,
                            folder_name = folder_name)
          gc()
          
          
          get_never_misspliced(cluster = cluster_id,
                               samples = samples_used,
                               split_read_counts = split_read_counts_tidy,
                               all_split_reads_details = all_split_reads_details_tidy,
                               folder_name = folder_name)
          
          rm(all_split_reads_details_tidy)
          rm(split_read_counts_tidy)
          rm(all_split_reads_details)
          rm(split_read_counts)
          gc()
          
        }
      }
    } else {
      message(Sys.time(), " there is no local metadata downloaded for '", project_id, "' project from recount3. ")
    }
  }
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


splicing_projects <- c( "ADIPOSE_TISSUE",  "ADRENAL_GLAND",   "BLOOD",           "BLOOD_VESSEL",    
                        "BONE_MARROW",    
                        "BRAIN",           "COLON",           "ESOPHAGUS",       "HEART",           "KIDNEY",          
                        "LIVER",           "LUNG",            "MUSCLE",          "NERVE",           "PANCREAS",        
                        "PITUITARY",       "SALIVARY_GLAND",  "SKIN",            "SMALL_INTESTINE", "SPLEEN",          
                        "STOMACH",         "THYROID")



for (gtf_version in gtf_versions) {
  
  # gtf_version <- gtf_versions[1]
  
  # download_recount3_data(recount3_project_IDs = all_projects,
  #                         gtf_version = gtf_version)
  # 
  # 
  # tidy_recount3_data(recount3_project_IDs = splicing_projects,
  #                               main_project,
  #                               gtf_version = gtf_version)
  
  
  # junction_pairing(recount3_project_IDs = splicing_projects,
  #                  main_project,
  #                  gtf_version = gtf_version)
  
  
  # get_all_annotated_split_reads(recount3_project_IDs = splicing_projects,
  #                               gtf_version = gtf_version,
  #                               main_project = main_project)


  # get_all_raw_distances_pairings(recount3_project_IDs = splicing_projects,
  #                                gtf_version = gtf_version,
  #                                main_project = main_project)
  # 
  # 
  # tidy_data_pior_sql(recount3_project_IDs = splicing_projects,
  #                    gtf_version = gtf_version,
  #                    main_project = main_project)
  #  
  #  
  #  generate_transcript_biotype_percentage(recount3_project_IDs = splicing_projects,
  #                                         homo_sapiens_v105_path = paste0(dependencies_folder,
  #                                                                         "/Homo_sapiens.GRCh38.105.chr.gtf"),
  #                                         main_project,
  #                                         gtf_version = gtf_version)
   
   
   all_final_projects_used <- readRDS(file = paste0(getwd(),"/results/all_final_projects_used.rds"))
   
   # homo_sapiens_v105 <- rtracklayer::import(con = paste0(dependencies_folder,"/Homo_sapiens.GRCh38.105.chr.gtf"))
   # generate_recount3_tpm(recount3_project_IDs = all_final_projects_used,
   #                       gtf_version = gtf_version,
   #                       ref = homo_sapiens_v105,
   #                       main_project = main_project)
  
   database_folder <- paste0(getwd(), "/database/v", gtf_version, "/", main_project)
   dir.create(file.path(database_folder), recursive = TRUE, showWarnings = T)
   database_path <- paste0(database_folder,  "/", main_project, ".sqlite")
   print("Starting SQL generation...")
   sql_database_generation(database_path = database_path,
                           recount3_project_IDs = all_final_projects_used,
                           main_project = main_project,
                           gtf_version = gtf_version,
                           remove_all = F)
  
  
  # rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  # gc()
}

