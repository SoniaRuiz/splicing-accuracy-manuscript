#' Title
#' downloads, quality-control and annotates the split read data from a given recount3 project
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615"
#' @param project.name name given locally to the recount3 project.
#' A given recount3 project can be separated in multiple independent projects in recount3.
#' For instance, GTEx is stored in recount3 in multiple independent projects ID (e.g. BRAIN, SKIN, BLOOD, etc), 
#' but all of them belong to GTEx. Hence, it is useulf to have a folder named "project.name" that will contain
#' multiple subfolders named as the elements contained within the object 'recount3.project.IDs'
#' @param gtf.version Ensembl version (tested using Ensembl v105)
#' e.g. "105"
#' @param data.source source of the data within recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' @export
#'
#' @examples
download_recount3_data <- function (recount3.project.IDs,
                                    project.name,
                                    gtf.version,
                                    data.source,
                                    database.folder,
                                    results.folder) {
  
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  if ( file.exists(paste0(here::here(), "/database/all_split_reads_raw.rds")) ) {
    
    message("Loading 'all_split_reads_raw.rds' file...")
    all_split_reads_raw <- readRDS(file = paste0(database.folder, "/all_split_reads_raw.rds"))
    
  } else {
    
    doParallel::registerDoParallel(10)
    all_split_reads_raw <- foreach(j = seq(length(recount3.project.IDs)), .combine = "rbind") %do% {
      
      project_id <- recount3.project.IDs[j]
      
      # project_id <- recount3.project.IDs[1]
      # project_id <- "KIDNEY"
      message(Sys.time(), " - getting data from '", project_id, "' recount3 project...")
      
      folder_root <- paste0(results.folder, "/", project_id, "/base_data/")
      dir.create(file.path(folder_root), recursive = TRUE, showWarnings = T)
      
      #############################################################################
      
      ## THERE ARE OTHER RECOMMENDED WAYS OF DOWNLOADING DATA FROM RECOUNT3
      ## HOWEVER, THIS ALTERNATIVE METHOD IS FOLLOWED IN ORDER TO REDUCE MEMORY USAGE
      ## PARTICULARLY USEFUL WITH LARGE DATASETS SUCH AS GTEX
      ## FOR OTHER RECOMMENDED METHODS, SEE (https://bioconductor.org/packages/release/bioc/manuals/recount3/man/recount3.pdf) - ACCESSED 08/07/2023
      
      #############################################################################
      
      
      jxn_files <- recount3::locate_url(
        project = project_id,
        project_home = data.source,
        type = "jxn",
        organism = "human",
        annotation = "gencode_v29",
        #jxn_format = c("UNIQUE"),
        recount3_url = getOption("recount3_url", "http://duffel.rail.bio/recount3")
      )
      
      feature_info <- utils::read.delim(recount3::file_retrieve(
        url = jxn_files[grep("\\.RR\\.gz$", jxn_files)],
        bfc = recount3::recount3_cache(),
        verbose = getOption("recount3_verbose", TRUE)
      ))
      
      ## Convert unstranded junctions from '*' to '?' to facilitate later conversion to genome:ranges
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
      
      saveRDS(object = all_split_reads %>% data.table::as.data.table(),
              file = paste0(folder_root, "/all_split_reads_raw.rds"))
      
      return(all_split_reads %>%
               distinct(junID, .keep_all = T))
    }
    
    # ## We access recount3 files directly to save memory
    # all_split_reads_raw <- map_df(recount3.project.IDs, function(project_id) {
    #   
    #   
    # })
    
    
    all_split_reads_raw <- all_split_reads_raw %>%
      distinct(junID, .keep_all = T)
    
    print(paste0(Sys.time(), " - ", all_split_reads_raw %>% nrow(), " initial number of split reads"))
    gc()
    
    ## Save data
    dir.create(file.path(database.folder), recursive = TRUE, showWarnings = T)
    saveRDS(object = all_split_reads_raw %>% data.table::as.data.table(),
            file = paste0(database.folder, "/all_split_reads_raw.rds"))
    
  }
  
  all_split_reads_raw %>% 
    distinct(junID) %>%
    nrow()
  
  #######################################################################
  ## 1. Discard all split reads shorter than 25 bp
  #######################################################################
  
  all_split_reads_raw %>% head()
  all_split_reads_raw %>% nrow()
  
  all_split_reads_raw_tidy <- all_split_reads_raw %>%
    filter(width >= 25)
  
  all_split_reads_raw_tidy %>% head()
  all_split_reads_raw_tidy %>% nrow()
  
  rm(all_split_reads_raw)
  gc()
  
  #######################################################################
  ## 2. Remove split reads located in unplaced sequences in the chromosomes
  #######################################################################
  
  print(paste0(Sys.time(), " - removing unplaced sequences genome..."))
  
  all_split_reads_raw_tidy_gr <- all_split_reads_raw_tidy %>%
    filter(!(strand=="?")) %>%
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
  
  ##########################################################
  ## 3. Remove split reads overlapping the ENCODE backlist
  ##########################################################
  
  print(paste0(Sys.time(), " - removing blacklist sequences..."))
  
  blacklist_path <- paste0(dependencies_folder, "/hg38-blacklist.v2.bed")
  all_split_reads_raw_tidy_gr <- remove_encode_blacklist_regions(GRdata = all_split_reads_raw_tidy_gr,
                                                                 blacklist_path = blacklist_path)
  
  ## Data check
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID, .keep_all = T) %>%
    nrow() %>%
    print()
  
  
  
  #######################################################
  ## 4. Anotate using 'dasper'
  #######################################################
  
  message("Annotating dasper...")
  
  if ( file.exists(paste0(dependencies_folder, "/Homo_sapiens.GRCh38.", gtf.version, ".chr.gtf")) ) {
    gtf_path <- paste0(dependencies_folder, "/Homo_sapiens.GRCh38.", gtf.version, ".chr.gtf")
  } else {
    gtf_path <- paste0(dependencies_folder, "/Homo_sapiens.GRCh38.", gtf.version, ".gtf")
  }
  edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  
  all_split_reads_details_w_symbol <- dasper::junction_annot(junctions = all_split_reads_raw_tidy_gr %>% GRanges(), 
                                                             ref = edb,
                                                             ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"), 
                                                             ref_cols_to_merge = c("gene_id", "gene_name", "tx_id"))
  all_split_reads_details_w_symbol %>% head
  rm(all_split_reads_raw_tidy_gr)
  rm(edb)
  gc()
  
  
  ################################################################################
  ## 5. Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  print(paste0(Sys.time(), " - removing split reads not classified as 'annotated', 'novel_donor' or 'novel_acceptor'..."))
  
  ## Subset columns
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[,c("junID", "gene_id_junction", "in_ref", "type", "tx_id_junction")]
  
  ## Only use annotated introns, novel donor and novel acceptor junctions
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[(elementMetadata(all_split_reads_details_w_symbol)[,"type"] %in% 
                                                                          c("annotated", "novel_donor", "novel_acceptor"))]
  
  
  all_split_reads_details_w_symbol %>% length()
  
  
  # ############################################
  # ## 6. Make sure none of the split reads are shorter than 25bp
  # ############################################
  # 
  # print(paste0(Sys.time(), " - removing split reads shorter than 25bp..."))
  # 
  # all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol[all_split_reads_details_w_symbol %>% width() >= 25,]
  
  
  ############################################################################
  ## 6. Discard all introns assigned to multiple genes (i.e. ambiguous introns)
  ############################################################################
  
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
          file = paste0(database.folder, "/all_ambiguous_jxn.rds"))
  
  all_split_reads_details_w_symbol <- all_split_reads_details_w_symbol %>%
    dplyr::filter(ambiguous == F) %>%
    dplyr::select(-ambiguous)
  
  all_split_reads_details_w_symbol %>% nrow()
  
  saveRDS(object = all_split_reads_details_w_symbol %>% dplyr::rename(gene_id = gene_id_junction),
          file = paste0(database.folder, "/all_split_reads_qc_level1.rds"))
  
  
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