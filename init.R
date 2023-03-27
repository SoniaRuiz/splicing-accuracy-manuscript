library(tidyverse)
library(GenomicRanges)

setwd(normalizePath("."))
dependencies_folder <- paste0(getwd(), "/dependencies/")

source(paste0(getwd(), "/database_junction_pairing.R"))
source(paste0(getwd(), "/database_SQL_helper.R"))
source(paste0(getwd(), "/database_SQL_generation.R"))



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
    
    ## We access recount3 files directly to reduce memory usage
    for ( project_id in projects_used ) {
      
      # project_id <- projects_used[1]
      # project_id <- "BRAIN"
      print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
      
      folder_root <- paste0(getwd(), "/results/", project_id, "/")
      dir.create(file.path(folder_root), recursive = TRUE, showWarnings = T)
      
      #############################################################################
      
      # metadata <- recount3::read_metadata(recount3::file_retrieve(
      #   url = recount3::locate_url(
      #     project = project_id,
      #     project_home = "data_sources/gtex",
      #     type = "metadata",
      #     organism = "human",
      #     annotation = "gencode_v29",
      #     recount3_url = getOption("recount3_url", "http://duffel.rail.bio/recount3"),
      #   ),
      #   bfc = recount3::recount3_cache(),
      #   verbose = getOption("recount3_verbose", TRUE)
      # ))
      
      jxn_files <- recount3::locate_url(
        project = project_id,
        project_home = "data_sources/gtex",
        type = "jxn",
        organism = "human",
        annotation = "gencode_v29",
        jxn_format = c("UNIQUE"),
        recount3_url = getOption("recount3_url", "http://duffel.rail.bio/recount3")
      )
      
      feature_info <- utils::read.delim(recount3::file_retrieve(
        url = jxn_files[grep("\\.RR\\.gz$", jxn_files)],
        bfc = recount3::recount3_cache(),
        verbose = getOption("recount3_verbose", TRUE)
      ))
      
      feature_info$strand[feature_info$strand == "?"] <- "*"
      #feature_info <- GenomicRanges::GRanges(feature_info)
      
      all_split_reads <- data.frame(chr = feature_info$chromosome,
                                    start = feature_info$start,
                                    end = feature_info$end,
                                    strand = feature_info$strand,
                                    width = feature_info$length,
                                    annotated = feature_info$annotated,
                                    junID = paste0(feature_info$chromosome, ":",
                                                   feature_info$start, "-", feature_info$end, ":",
                                                   feature_info$strand))
      
      
      # counts <- Matrix::readMM(recount3::file_retrieve(
      #   url = jxn_files[grep("\\.MM\\.gz$", jxn_files)],
      #   bfc = recount3::recount3_cache(),
      #   verbose = getOption("recount3_verbose", TRUE)
      # ))
      # gc()
      # 
      # counts %>% head()
      # 
      # if (verbose) {
      #   message(
      #     Sys.time(),
      #     " matching exon-exon junction counts with the metadata."
      #   )
      # }
      # ## The samples in the MM jxn table are not in the same order as the
      # ## metadata!
      # jxn_rail <- read.delim(recount3::file_retrieve(
      #   url = jxn_files[grep("\\.ID\\.gz$", jxn_files)],
      #   bfc = recount3::recount3_cache(),
      #   verbose = getOption("recount3_verbose", TRUE)
      # ))
      # m <- match(metadata$rail_id, jxn_rail$rail_id)
      # stopifnot(
      #   "Metadata rail_id and exon-exon junctions rail_id are not matching." =
      #     !all(is.na(m))
      # )
      # counts <- counts[, m, drop = FALSE]
      # colnames(counts) <- metadata$external_id
      
      ######################################################
      
      # rse <- recount3::create_rse_manual (
      #   project = project_id,
      #   project_home = "data_sources/gtex",
      #   organism = "human",
      #   annotation = "gencode_v29",
      #   jxn_format = c("UNIQUE"),
      #   type = "jxn"
      # )
      # 
      # all_split_reads <- data.frame(chr = rse %>% SummarizedExperiment::seqnames(),
      #                               start = rse %>% SummarizedExperiment::start(),
      #                               end = rse %>% SummarizedExperiment::end(),
      #                               strand = rse %>% SummarizedExperiment::strand(),
      #                               width = rse %>% SummarizedExperiment::width(),
      #                               junID = rse %>% rownames())
      # rm(rse)
      # gc()
      
      ######################################################
      
      all_split_reads %>% nrow()
      saveRDS(object = all_split_reads %>% data.table::as.data.table(),
              file = paste0(folder_root, "/all_split_reads_raw.rds"))
      
      rm(all_split_reads)
      gc()
    }
    
    gc()
    
    all_split_reads_raw <- map_df(projects_used, function(project_id) {
      # project_id <- projects_used[31]
      print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
      
      folder_root <- paste0(getwd(), "/results/", project_id, "/")
      all_split_reads <- readRDS(file = paste0(folder_root, "/all_split_reads_raw.rds"))
      
      return(all_split_reads %>%
               distinct(junID, .keep_all = T))
    })
    
    all_split_reads_raw_tidy <- all_split_reads_raw %>%
      distinct(junID, .keep_all = T)
    
    all_split_reads_raw_tidy %>% nrow() %>% print()
    gc()
    
    ## Save data
    folder_path <- paste0(getwd(), "/database/")
    dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
    saveRDS(object = all_split_reads_raw_tidy %>% data.table::as.data.table(),
            file = paste0(folder_path, "/all_split_reads_raw.rds"))
    
  }
  
  # all_split_reads_raw_tidy %>% 
  #   distinct(junID) %>%
  #   nrow()
  
  #######################################################################
  ## Remove split reads located in unplaced sequences in the chromosomes
  #######################################################################
  
  print(paste0(Sys.time(), " - removing unplaced sequences genome..."))
  
  all_split_reads_raw_tidy_gr <- all_split_reads_raw_tidy %>%
    GenomicRanges::GRanges() %>%
    diffloop::rmchr()
  rm(all_split_reads_raw_tidy)
  gc()
  
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
  #all_split_reads_raw_tidy %>% nrow() - all_split_reads_raw_tidy_gr %>% length()
  
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
  ## Remove split reads shorter than 25 bp
  #######################################################
  
  print(paste0(Sys.time(), " - removing split reads shorter than 25bp..."))
  
  indx <- which(all_split_reads_raw_tidy_gr %>% width() < 25)
  all_split_reads_raw_tidy_gr <- all_split_reads_raw_tidy_gr[-indx,]
  
  any(all_split_reads_raw_tidy_gr %>% width() < 25)
  
  rm(indx)
  gc()
  
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
  
  rm(all_split_reads_raw_tidy_gr)
  rm(edb)
  gc()
  
  
  ################################################################################
  ## Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  print(paste0(Sys.time(), " - removing split reads not classified as 'annotated', 'novel_donor' or 'novel_acceptor'..."))
  
  ## Subset columns
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[,c("junID", "gene_id_junction", "in_ref", "type", "tx_id_junction")]
  
  ## Only use annotated introns, novel donor and novel acceptor junctions
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[(elementMetadata(all_split_reads_details_w_symbol)[,"type"] %in% c("annotated", "novel_donor", "novel_acceptor"))]
  
  
  all_split_reads_details_w_symbol %>% length()
  
  
  ############################################
  ## Discard all junctions shorter than 25bp
  ############################################
  
  print(paste0(Sys.time(), " - removing split reads shorter than 25bp..."))
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[all_split_reads_details_w_symbol %>% width() >= 25,]
  
  
  
  folder_path <- paste0(getwd(), "/database/v", gtf_version, "/")
  dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
  saveRDS(all_split_reads_details_w_symbol,
          file = paste0(folder_path, "/all_split_reads_gtex_recount3_", gtf_version, ".rds"))
  
  
  
  
  ############################################################################
  ## Discard all introns assigned to multiple genes (i.e. ambiguous introns)
  ############################################################################
  
  if ( !exists("all_split_reads_details_w_symbol") ) {
    folder_path <- paste0(getwd(), "/database/v", gtf_version, "/")
    all_split_reads_details_w_symbol <- readRDS(file = paste0(folder_path, "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))
  }
  
  print(paste0(Sys.time(), " - discarding ambiguous introns..."))
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id_junction %>% unlist() %>% length() > 1, T, F))
  
  ambiguous_introns <- all_split_reads_details_w_symbol %>%
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
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol %>%
    dplyr::filter(ambiguous == F) %>%
    dplyr::select(-ambiguous)
  
  all_split_reads_details_w_symbol %>% nrow()
  
  saveRDS(object = all_split_reads_details_w_symbol %>% dplyr::rename(gene_id = gene_id_junction),
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
  all_projects_used <- NULL
  
  message(Sys.time()," loading tidy split reads ID from recount3.")
  all_split_reads_details_w_symbol_reduced_keep_gr <- 
    readRDS(file = paste0(folder_database, "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))
  gc()
  
  # ## Generate the raw file
  for (project_id in projects_used) {
    
    # project_id <- projects_used[1]
    # project_id <- "BRAIN"
    
    folder_results <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", 
                             main_project, "/base_data/")
    dir.create(file.path(folder_results), recursive = TRUE, showWarnings = T)
    print(paste0(Sys.time(), " - getting junction data from recount3 - '", project_id, "' tissue ..."))
    
    ###########################################################################
    
    bfc <- recount3::recount3_cache()
    recount3_url <- getOption("recount3_url", "http://duffel.rail.bio/recount3")
    verbose <- getOption("recount3_verbose", TRUE)
    
    message(Sys.time()," loading metadata info.")
    metadata.info <- recount3::read_metadata(recount3::file_retrieve(
      url = recount3::locate_url(
        project = project_id,
        project_home = "data_sources/gtex",
        organism = "human",
        annotation = "gencode_v29",
        type = "metadata",
        recount3_url = recount3_url
      ),
      bfc = bfc,
      verbose = verbose
    ))
    
    jxn_files <- recount3::locate_url(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      jxn_format = c("UNIQUE"),
      type = "jxn",
      recount3_url = recount3_url
    )
    
    # feature_info <- utils::read.delim(recount3::file_retrieve(
    #   url = jxn_files[grep("\\.RR\\.gz$", jxn_files)],
    #   bfc = bfc,
    #   verbose = verbose
    # ))
    # feature_info$strand[feature_info$strand == "?"] <- "*"
    # feature_info <- as.character(GenomicRanges::GRanges(feature_info))
    
    feature_info <- 
      readRDS(file = paste0(getwd(), "/results/", project_id, "/all_split_reads_raw.rds")) %>%
      pull(junID)
    gc()
    
    
    
    
    if (verbose) {
      message(
        Sys.time(),
        " matching exon-exon junction counts with the metadata."
      )
    }
    ## The samples in the MM jxn table are not in the same order as the
    ## metadata!
    jxn_rail <- read.delim(recount3::file_retrieve(
      url = jxn_files[grep("\\.ID\\.gz$", jxn_files)],
      bfc = bfc,
      verbose = verbose
    ))
    m <- match(metadata.info$rail_id, jxn_rail$rail_id)
    stopifnot(
      "Metadata rail_id and exon-exon junctions rail_id are not matching." =
        !all(is.na(m))
    )
    
    
    
    
    ###########################################################################
    
    
    # rse <- recount3::create_rse_manual(
    #   project = project_id,
    #   project_home = "data_sources/gtex",
    #   organism = "human",
    #   annotation = "gencode_v29",
    #   jxn_format = c("UNIQUE"),
    #   type = "jxn"
    # )
    
    
    #################################
    ## GET METADATA
    #################################
    
    # Samples passing the GTEx QC filter and with RNA quality greater or equal than 6
    metadata.info <- metadata.info %>%
      as_tibble() %>%
      filter(gtex.smafrze != "EXCLUDE",
             gtex.smrin >= 6.0)

    
    saveRDS(object = metadata.info, 
            file = paste0(folder_results, "/", project_id, "_samples_metadata.rds"))
    # metadata.info <- readRDS(file = paste0(folder_results, "/", project_id, "_samples_metadata.rds")) 
    
    
    
    #################################
    ## GET SPLIT READS AND COUNTS  
    ## PER TISSUE
    #################################
    
    clusters_ID <- metadata.info$gtex.smtsd %>% unique()
    clusters_used <- NULL
    gc()
    
    for (cluster_id in clusters_ID) {
      
      # cluster_id <- clusters_ID[1]
      print(paste0(Sys.time(), " - filtering junction data by cluster - '", cluster_id, "' tissue"))
      
      message(Sys.time()," loading count matrix.")
      counts <- Matrix::readMM(recount3::file_retrieve(
        url = jxn_files[grep("\\.MM\\.gz$", jxn_files)],
        bfc = bfc,
        verbose = verbose
      ))
      #counts <- counts[, m, drop = FALSE]
      message(Sys.time()," ordering count matrix.")
      colnames(counts) <- metadata.info$external_id[m]
      rownames(counts) <- feature_info
      gc()
      
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
      
      ################
      ## Get counts
      ################
      
      ## Tissues at least with 70 samples
      if ( (cluster_samples %>% length() >= minimum_samples ) ) {
        
        clusters_used <- c(clusters_used, cluster_id)
        all_projects_used <- c(all_projects_used, project_id)
        print(paste0(Sys.time(), " - filtering split read counts matrix for the current matrix."))
        
        local_counts <- counts[,(colnames(counts) %in% cluster_samples)]
        gc()
        local_counts <- local_counts[(rownames(local_counts) %in% 
                                        all_split_reads_details_w_symbol_reduced_keep_gr$junID),]
        gc()
        
        print(paste0(Sys.time(), " - converting the split read counts matrix into tibble."))
        local_counts <- local_counts %>% as.matrix()
        local_counts <- local_counts[rowSums(local_counts) > 0, ]
        gc()
        local_counts <- local_counts %>% as_tibble(rownames = "junID")
        local_counts %>% nrow()
        gc()
        
        print(paste0(Sys.time(), " - saving data."))
        print(object.size(local_counts), units = "GB")
        saveRDS(object = local_counts,
                file = paste0(folder_results, "/", project_id, "_", cluster_id, "_split_read_counts.rds"))
        gc()
        
        #######################
        ## Get split reads ID
        #######################
        
        print(paste0(Sys.time(), " - getting split read IDs"))
        all_split_reads <- local_counts %>%
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
        
        print(paste0(Sys.time(), " - saving data."))
        ## Save the object
        saveRDS(object = all_split_reads %>% data.table::as.data.table(),
                file = paste0(folder_results, "/", project_id, "_", cluster_id, "_all_split_reads.rds"))
        
        
        ## Free up some memory
        rm(all_split_reads)
        rm(local_counts)
        rm(cluster_samples)
        rm(counts)
        gc()
        
      } 
    }
    
    if ( !is.null(clusters_used) ) {
      saveRDS(object = clusters_used, 
              file = paste0(folder_results, "/", project_id, "_clusters_used.rds"))
    }
    
    #rm(rse)
    rm(metadata.info)
    rm(clusters_used)
    rm(clusters_ID)
    
    rm(feature_info)
    rm(jxn_files)
    rm(jxn_rail)
    rm(m)
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
  
  # gtf_version <- 105
  
  # init_recount3_gtex_data(projects_used = splicing_projects,
  #                         gtf_version = gtf_version)
  
  
  # tidy_recount3_data_per_tissue(projects_used = splicing_projects[c(6,17:22)],
  #                               main_project,
  #                               gtf_version = gtf_version)
  # 
  # 
  # junction_pairing(projects_used = splicing_projects[-c(17:22)],
  #                  main_project,
  #                  gtf_version = gtf_version)
  
  
  # get_all_annotated_split_reads(projects_used = splicing_projects[-c(6,17:22)],
  #                               gtf_version = gtf_version,
  #                               main_project = main_project)


  # get_all_raw_distances_pairings(projects_used = splicing_projects[-c(6,17:22)],
  #                                gtf_version = gtf_version,
  #                                main_project = main_project)
  # 
  # 
  # tidy_data_pior_sql(projects_used = splicing_projects[-c(6,17:22)],
  #                    gtf_version = gtf_version,
  #                    main_project = main_project)
  # 
  # 
  # generate_transcript_biotype_percentage(projects_used = splicing_projects[-c(6,17:22)],
  #                                        homo_sapiens_v105_path = paste0(dependencies_folder,
  #                                                                        "/Homo_sapiens.GRCh38.105.chr.gtf"),
  #                                        main_project,
  #                                        gtf_version = gtf_version)


  all_final_projects_used <- readRDS(file = paste0(getwd(),"/results/", main_project, "_final_projects_used.rds"))

  # generate_recount3_tpm(projects_used = all_final_projects_used,
  #                       gtf_version = gtf_version,
  #                       main_project = main_project)

  database_folder <- paste0(getwd(), "/database/v", gtf_version, "/", main_project)
  dir.create(file.path(database_folder), recursive = TRUE, showWarnings = T)
  database_path <- paste0(database_folder,  "/", main_project, ".sqlite")

  sql_database_generation(database_path = database_path,
                          projects_used = all_final_projects_used[9:13],
                          main_project = main_project,
                          gtf_version = gtf_version,
                          remove_all = F)

  
  # rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  # gc()
}

