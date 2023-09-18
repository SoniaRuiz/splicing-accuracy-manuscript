source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/conservation/conservation_general/conservation_general_functions_from_bw.R")
source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/constraint/constraint_general/constraint_general_functions_bw.R")

#' Title
#' Calculates the CDTS and PhastCons20 mean scores of the sequences overlapping the -5/+35bp sequence at the donor site of the intron (i.e. the exon-intron junction), and the -5/+35bp sequence at the acceptor splice 
#' site of the intron (i.e. the intron-exon junction)
#' @param cluster 
#' @param db_introns Dataframe of introns to calculate the scores
#' @param folder_name 
#'
#' @return
#' @export
#'
#' @examples
generate_cdts_phastcons_scores <- function(db_introns = NULL,
                                           cluster = NULL,
                                           folder_name = NULL,
                                           exon_size = 5,
                                           intron_size = 35) {
  
  
  if (is.null(db_introns) && !is.null(cluster)) {
    ## Load the IDB 
    db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
  }
  
  db_introns <- db_introns %>%
    GRanges() %>% 
    diffloop::addchr()
  
      
  
  ###########################################
  ## PHASTCONS SCORES
  ###########################################
  
  for (intron_size in c(35, 70, 100)) {
    
    for (phastcons_type in c(4, 7, 20)) {
      
      # phastcons_type <-20
      message("PhastCons", phastcons_type, "...")
      
      bw_path = paste0("/data/conservation/phastCons/hg38.phastCons", phastcons_type, "way.bw")
      
      if (!file.exists(bw_path)) {
        paste0("PhastCons file '", bw_path, "' does not exist!") %>% print()
        break;
      }
        
      
      # intron_size <- 35
      
      ## Calculate donor scores
      message(intron_size, "bp - Calculating PhastCons" , phastcons_type, " scores overlapping donor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                   ranges = IRanges(start = db_introns %>% start() - exon_size, 
                                                    end = db_introns %>% start() + intron_size))
      values(gr) <- DataFrame(ref_junID = (db_introns)$ref_junID )
      phastCons_5ss <- get_conservation_score_for_regions_bw(bw_path = bw_path,
                                                             gr = gr, 
                                                             summaryFun = "mean") %>% 
        as_tibble() %>%
        dplyr::rename_with(.fn = ~paste0(., "5ss_", intron_size), .cols = paste0("mean_phastCons", phastcons_type, "way"))
      
      
      
      
      ## Calculate acceptor scores
      message(intron_size, "bp - Calculating PhastCons" , phastcons_type, " scores overlapping acceptor sequences...")
      gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                   ranges = IRanges(start = db_introns %>% end() - intron_size, 
                                                    end = db_introns %>% end() + exon_size))
      values(gr) <- DataFrame(ref_junID = (db_introns)$ref_junID )
      phastCons_3ss <- get_conservation_score_for_regions_bw(bw_path = bw_path,
                                                             gr = gr, 
                                                             summaryFun = "mean") %>% 
        as_tibble()  %>%
        dplyr::rename_with(.fn = ~paste0(., "3ss_", intron_size), .cols = paste0("mean_phastCons", phastcons_type, "way"))
      
      
      
      ## Add columns to master data
      db_introns <- db_introns %>%
        as_tibble() %>%
        inner_join(phastCons_5ss %>% dplyr::select(ref_junID, paste0("mean_phastCons", phastcons_type, "way5ss_", intron_size)),
                   by = "ref_junID") %>%
        inner_join(phastCons_3ss %>% dplyr::select(ref_junID, paste0("mean_phastCons", phastcons_type, "way3ss_", intron_size)),
                   by = "ref_junID") %>%
        GRanges()
      
      
    }
    
    ###########################################
    ## CDTS SCORES
    ###########################################
    
    bw_path = paste0(here::here(), "/natcomms_review/constraint/CDTS_percentile_N7794_unrelated_all_chrs.bw")
    
    # intron_size <- 35
    
    ## Calculate donor scores
    message(intron_size, "bp - CDTS calculating donor sequences")
    gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                 ranges = IRanges(start = db_introns %>% start() - exon_size, 
                                                  end = db_introns %>% start() + intron_size))
    values(gr) <- DataFrame(ref_junID = (db_introns)$ref_junID )
    CDTS_5ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
                                                    gr = gr,
                                                    summaryFun = "mean") %>% 
      as_tibble() %>%
      dplyr::rename_with(.fn = ~paste0(., "5ss_", intron_size), .cols = mean_CDTS)
    
    
    
    
    ## Calculate acceptor scores
    message(intron_size, "bp - CDTS calculating acceptor sequences")
    gr <- GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                 ranges = IRanges(start = db_introns %>% end() - intron_size, 
                                                  end = db_introns %>% end() + exon_size))
    values(gr) <- DataFrame(ref_junID = (db_introns)$ref_junID )
    CDTS_3ss <- get_constraint_score_for_regions_bw(bw_path = bw_path,
                                                    gr = gr, 
                                                    summaryFun = "mean") %>% 
      as_tibble()  %>%
      dplyr::rename_with(.fn = ~paste0(., "3ss_", intron_size), .cols = paste0("mean_CDTS"))
    
    
    
    ## Add columns to master data
    db_introns <- db_introns %>%
      as_tibble() %>%
      inner_join(CDTS_5ss %>% dplyr::select(ref_junID, paste0("mean_CDTS5ss_", intron_size)),
                 by = "ref_junID") %>%
      inner_join(CDTS_3ss %>% dplyr::select(ref_junID, paste0("mean_CDTS3ss_", intron_size)),
                 by = "ref_junID") %>%
      GRanges()
    
  }
  
  
  ####################
  ## Save results
  ####################
  
  if (!is.null(cluster)) {

    saveRDS(object = db_introns,
            file = paste0(folder_name, "/", cluster, "_db_introns.rds"))

    print(paste0(Sys.time(), " - CDTS and PhastCons scores added! Database updated!"))


    rm(gr)
    rm(db_introns)
    rm(phastCons_5ss)
    rm(phastCons_3ss)

  } else {
    return(db_introns)
  }
      
  
}




  
  
