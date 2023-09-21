
#' Title
#' Separates the quality-controlled split-read data by sample cluster for a given recount3 project
#' (e.g. by case/control cluster, by tissue, etc. )
#' @param recount3.project.IDs array with the recount3 projects to download
#' e.g. "SRP009615"
#' @param data.source source of the project in recount3. It can be:
#' - "data_sources/sra"  
#' - "data_sources/gtex" 
#' - "data_sources/tcga"
#'
#' @return
#' @export
#'
#' @examples
prepare_recount3_data <- function(recount3.project.IDs, 
                                  data.source,
                                  results.folder,
                                  subsampling = F,
                                  levelqc1.folder,
                                  supporting.reads,
                                  num.cores) {
  
  
  message(Sys.time()," loading all QC split reads ...")
  
  if ( !file.exists(paste0(levelqc1.folder, "/all_split_reads_qc_level1.rds")) ) {
    message("Error! The file with the split reads passing LEVEL 1 of filtering criteria does not exist!")
    break;
  } else {
    all_split_reads_qc_level1 <- 
      readRDS(file = paste0(levelqc1.folder, "/all_split_reads_qc_level1.rds"))  
  }
  
  message(all_split_reads_qc_level1 %>% nrow(), " initial number of split reads...")
  
  doParallel::registerDoParallel(num.cores)
  foreach(i = seq(length(recount3.project.IDs))) %dopar%{
    
    project_id <- recount3.project.IDs[i]
    
    
    ## If this is GTEx, the recount3.project.IDs object can contain up to 54 different values
    ## e.g. KIDNEY, BRAIN, BLOOD, etc
    # project_id <- recount3.project.IDs[12]
    # project_id <- "LUNG"
    
    
    local_folder_results <- paste0(results.folder, "/", project_id, "/base_data/")
    dir.create(file.path(local_folder_results), recursive = TRUE, showWarnings = T)
    message("Getting junction data from recount3 - '", project_id, "' ...")
    
    
    rse_jxn <- recount3::create_rse_manual(
      project = project_id,
      project_home = data.source,
      organism = "human",
      annotation = "gencode_v29",
      type = "jxn"
    )
    
    #############################################################################
    
    ## THERE ARE OTHER RECOMMENDED WAYS OF DOWNLOADING DATA FROM RECOUNT3
    ## HOWEVER, THIS ALTERNATIVE METHOD IS FOLLOWED IN ORDER TO REDUCE MEMORY USAGE
    ## PARTICULARLY USEFUL WITH LARGE DATASETS SUCH AS GTEX (https://github.com/LieberInstitute/recount3/blob/devel/R/create_rse_manual.R)
    ## FOR OTHER RECOMMENDED METHODS, SEE (https://bioconductor.org/packages/release/bioc/manuals/recount3/man/recount3.pdf) - ACCESSED 08/07/2023
    
    #############################################################################
    
    # bfc <- recount3::recount3_cache()
    # recount3_url <- getOption("recount3_url", "http://duffel.rail.bio/recount3")
    # verbose <- getOption("recount3_verbose", TRUE)
    # 
    # message(Sys.time()," loading metadata info.")
    # metadata.info <- recount3::read_metadata(recount3::file_retrieve(
    #   url = recount3::locate_url(
    #     project = project_id,
    #     project_home = data.source,
    #     organism = "human",
    #     annotation = "gencode_v29",
    #     type = "metadata",
    #     recount3_url = recount3_url
    #   ),
    #   bfc = bfc,
    #   verbose = verbose
    # ))
    # 
    # jxn_files <- recount3::locate_url(
    #   project = project_id,
    #   project_home = data.source,
    #   organism = "human",
    #   annotation = "gencode_v29",
    #   #jxn_format = c("UNIQUE"),
    #   type = "jxn",
    #   recount3_url = recount3_url
    # )
    # 
    # feature_info <- utils::read.delim(recount3::file_retrieve(
    #   url = jxn_files[grep("\\.RR\\.gz$", jxn_files)],
    #   bfc = bfc,
    #   verbose = verbose
    # ))
    # feature_info$strand[feature_info$strand == "?"] <- "*"
    # feature_info <- GenomicRanges::GRanges(feature_info)
    # feature_info %>% length()
    # 
    # # feature_info_prev <- readRDS(file = paste0(local_folder_results, "/all_split_reads_raw.rds")) %>% pull(junID) 
    # # gc()
    # 
    # if (verbose) {
    #   message(
    #     Sys.time(),
    #     " matching exon-exon junction counts with the metadata."
    #   )
    # }
    # 
    # ## The samples in the MM jxn table are not in the same order as the metadata!
    # jxn_rail <- read.delim(recount3::file_retrieve(
    #   url = jxn_files[grep("\\.ID\\.gz$", jxn_files)],
    #   bfc = bfc,
    #   verbose = verbose
    # ))
    # m <- match(metadata.info$rail_id, jxn_rail$rail_id)
    # stopifnot(
    #   "Metadata rail_id and exon-exon junctions rail_id are not matching." =
    #     !all(is.na(m))
    # )
    # 
    # 
    # 
    # #################################
    # ## GET SPLIT READS AND COUNTS  
    # #################################
    # 
    # message(Sys.time()," loading count matrix.")
    # counts <- Matrix::readMM(recount3::file_retrieve(
    #   url = jxn_files[grep("\\.MM\\.gz$", jxn_files)],
    #   bfc = bfc,
    #   verbose = verbose
    # ))
    # counts %>% nrow()
    # 
    # message(Sys.time()," ordering count matrix.")
    # counts <- counts[, m, drop = FALSE]
    # colnames(counts) <- metadata.info$external_id[m]
    # rownames(counts) <- feature_info %>% as.character()
    # colnames(counts) %>% length()
    # 
    # ## Make names consistent
    # names(feature_info) <- rownames(counts)
    # rownames(metadata.info) <- colnames(counts)
    # 
    
    
    metadata.info <- colData(rse_jxn)
    saveRDS(object = metadata.info, 
            file = paste0(local_folder_results, "/", project_id, "_samples_raw_metadata.rds"))
    
    
    #################################
    ## GET SAMPLE CLUSTERS
    #################################
    
    # metadata.info <- readRDS(file = paste0(folder_results, "/", project_id, "_samples_raw_metadata.rds"))
    metadata_tidy <- separate_clusters(project.metadata = metadata.info, data.source)

    if ( metadata_tidy %>% nrow() > 0) {
      
      if (subsampling) {
        set.seed(12)
        ## From the samples obtained, only keep the ones matching similar RIN numbers
        m.out <- MatchIt::matchit(cluster ~ rin_score,
                                  data = metadata_tidy,
                                  distance = metadata_tidy$rin_score,
                                  method = "nearest",
                                  caliper = c(rin_score = 0),
                                  std.caliper = F)
        metadata_tidy_filter <- MatchIt::match.data(m.out) %>%
          dplyr::select(-c("distance", "weights", "subclass"))
        metadata_tidy_filter %>% 
          distinct(external_id, .keep_all = T) %>%
          count(cluster)
        metadata_tidy_filter %>% 
          distinct(external_id, .keep_all = T) %>%
          as_tibble()
        
      } else {
        
        metadata_tidy_filter <- metadata_tidy
      }
      
      
      stopifnot(
        "There are samples with RIN lower than 6!" =
          all(metadata_tidy_filter$rin >= 6)
      )
      
      
      saveRDS(object = metadata_tidy_filter, 
              file = paste0(local_folder_results, "/", project_id, "_samples_metadata.rds"))
      
      
      # metadata_tidy_filter$all_mapped_reads %>% min
      # metadata_tidy_filter %>% as_tibble()
      # metadata_tidy_filter$external_id %>% sort()
      
      ###################################
      ## GET SPLIT READ DATA PER CLUSTER
      ###################################
      
      ## Separate data per cluster
      clusters_ID <- metadata_tidy$cluster %>% unique()
      clusters_used <- NULL
      gc()
      
      for (cluster_id in clusters_ID) {
        
        # cluster_id <- clusters_ID[1]
        print(paste0(Sys.time(), " - '", cluster_id, "' tissue"))
        
        ################
        ## Get clusters
        ################
        
        ## Only consider samples that passed the filtering process performed above
        ## (e.g. RIN score and other filters, in case of the 'gtex' project) 
        cluster_samples <- metadata_tidy_filter %>% 
          filter(cluster == cluster_id) %>%
          distinct(external_id) %>% 
          pull()
        
        ## When working with 'gtex' data, only tissues with at least 70 samples were considered
        if ( (data.source == "data_sources/gtex" && cluster_samples %>% length() >= 70) ||
             (data.source != "data_sources/gtex" && cluster_samples %>% length() >= 1) ) {
          
          saveRDS(object = cluster_samples, 
                  file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_samples_used.rds"))
          
          clusters_used <- c(clusters_used, cluster_id)
          
          ################
          ## Get counts
          ################
          
          message("Getting split read counts matrix for the current cluster...")
          
          ## Only samples passing the filtering criteria
          local_counts <- assay(rse_jxn, "counts")[,(colnames(assay(rse_jxn, "counts")) %in% cluster_samples)]
          local_counts %>% nrow()
          
          ## Only split reads passing the LEVEL1 filtering criteria
          local_counts <- local_counts[(rownames(local_counts) %in% 
                                          all_split_reads_qc_level1$junID),]
          local_counts %>% nrow()
          
          message(cluster_id, " - converting the split read counts matrix into tibble ...")
          
          local_counts %>% nrow()
          local_counts %>% head()
          local_counts <- local_counts %>% as.matrix()
        
          ## At least N number of supporting reads
          local_counts <- local_counts[rowSums(local_counts) >= supporting.reads, ]
          local_counts %>% nrow()

          #local_counts["chr1:43422633-43422820:+", ]
          
          local_counts <- local_counts %>% as_tibble(rownames = "junID")

          print(paste0(Sys.time(), " - saving data."))
          print(object.size(local_counts), units = "GB")
          
          stopifnot(
            "Still there are split reads with a lower number of supporting reads indicated by parameter." =
              local_counts %>% 
              mutate(sumCounts = rowSums(select(., !contains("junID")))) %>%
              filter(sumCounts < supporting.reads) %>% 
              nrow() == 0
          )
          
          saveRDS(object = local_counts,
                  file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_split_read_counts.rds"))
          gc()
          
          #######################
          ## Get split reads ID
          #######################
          
          message(cluster_id, " - getting split read IDs")
          
          all_split_reads <- local_counts %>%
            dplyr::select(junID) %>%
            data.table::as.data.table() %>%
            inner_join(y = all_split_reads_qc_level1,
                       by = "junID")
          
          #######################
          ## qc and save results
          #######################
          
          if (any(all_split_reads$width < 25) |
              any(str_detect(string = all_split_reads$chr, pattern = "random")) |
              any(str_detect(string = str_to_lower(all_split_reads$chr), pattern = "u")) |
              any(!(all_split_reads$type %in% c("annotated", "novel_acceptor", "novel_donor")))) {
            print("ERROR! The split reads do not meet the minimum LEVEL1 filtering criteria.")
            break;
          }
          
          ## Check how much memory this object uses
          #print(object.size(all_split_reads %>% data.table::as.data.table()), units = "GB")
          
          message(cluster_id, " - saving data.")
          
          ## Save the object
          saveRDS(object = all_split_reads %>% data.table::as.data.table(),
                  file = paste0(local_folder_results, "/", project_id, "_", cluster_id, "_all_split_reads.rds"))
          
          
          ## Free up some memory
          rm(all_split_reads)
          rm(local_counts)
          
        }
        
        rm(cluster_samples)
        gc()
      }
      
      
      if ( !is.null(clusters_used) ) {
        saveRDS(object = clusters_used, 
                file = paste0(local_folder_results, "/", project_id, "_clusters_used.rds"))
      }
      
    } else {
      print(paste0(project_id, " does not have any sample that meet the minimum criteria for inclusion."))
    }
    
    
    rm(metadata.info)
    rm(metadata_tidy)
    rm(rse_jxn)
    gc()
  }
  
}