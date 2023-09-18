#' Title
#' Removes the ambiguous junctions and prepares the data prior generation of the SQL database
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param project.name Name given to the project 
#' @param all.clusters Clusters of samples. In GTEx projects, samples were clustered by tissue (eg. 'Puituitary', 'Thyroid', etc)
#'
#' @return
#' @export
#'
#' @examples
tidy_data_pior_sql <- function (recount3.project.IDs,
                                project.name,
                                all.clusters = NULL,
                                database.folder,
                                results.folder) {
  
  
  message("Loading split reads QC level 1 ...")
  
  ## Load base recount3 object containing the split reads passing the QC criteria
  if (str_detect(database.folder, pattern = "age")) {
    all_split_reads_details_qc_level1 <- readRDS(file = paste0(here::here(), "/database/",project.name,"/105/all_split_reads_qc_level1.rds")) %>%
      as_tibble()
  } else {
    all_split_reads_details_qc_level1 <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level1.rds")) %>%
      as_tibble()
  }
  all_split_reads_details_qc_level1 %>% nrow()
  
  ############################################
  ## Discard all junctions from 
  ## 'EXCLUDE ME' samples, tissues with < 70 samples
  ## samples RIN < 6, Brain Cortex and Cerebellum
  ############################################
  
  message("Loading split reads QC level 2 ...")
  
  all_split_reads_details_qc_level2 <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level2.rds"))  %>%
    as_tibble()
  
  all_split_reads_details_qc_level2 %>% nrow()
  all_split_reads_details_qc_level2 %>% head()
  
  ## This should be zero
  if ( setdiff(all_split_reads_details_qc_level2$junID, 
               all_split_reads_details_qc_level1$junID) %>% 
       length() > 0) {
    print("ERROR! Some of the annotated split reads that passed the 2nd QC level are not found within the split reads from the 1st QC level.")
    break;
  }
  
  # ## These are the number of split reads from the samples excluded
  setdiff(all_split_reads_details_qc_level1$junID, 
          all_split_reads_details_qc_level2$junID) %>% unique %>% length()
  
  
  
  ############################################
  ## QC
  ############################################ 
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_qc_level1$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_qc_level1[ind, "junID"] <- 
      str_replace(string = all_split_reads_details_qc_level1[ind, "junID"]$junID, 
                  pattern = "\\*", 
                  replacement = all_split_reads_details_qc_level1[ind, "strand"]$strand %>% as.character() )
    if (any(str_detect(all_split_reads_details_qc_level1$junID, pattern = "\\*"))) {
      print("ERROR!")
      break;
    }
  }
  
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_qc_level2$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_qc_level2[ind, "junID"] <- str_replace(string = all_split_reads_details_qc_level2[ind, "junID"]$junID, 
                                                                   pattern = "\\*", 
                                                                   replacement = all_split_reads_details_qc_level2[ind, "strand"]$strand %>% as.character() )
    if (any(str_detect(all_split_reads_details_qc_level2$junID, pattern = "\\*")) ) {
      print("ERROR!")
      break;
    }
  }
  
  
  ##########################################
  ## Load all distances pairings
  ##########################################
  
  df_all_jxn_pairings <- readRDS(file = paste0(database.folder, "/all_raw_jxn_pairings.rds"))
  
  ## QC
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = df_all_jxn_pairings$ref_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_jxn_pairings[ind, "ref_junID"] <- str_replace(string = df_all_jxn_pairings[ind, "ref_junID"]$ref_junID, 
                                                         pattern = "\\*", 
                                                         replacement = df_all_jxn_pairings[ind, "ref_strand"]$ref_strand %>% as.character())
    if( any(str_detect(df_all_jxn_pairings$ref_junID, pattern = "\\*")) ) {
      print("ERROR!")
      break;
    }
  }
  ## Remove potential * in the junID of the novel junctions
  ind <- which(str_detect(string = df_all_jxn_pairings$novel_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_jxn_pairings[ind, "novel_junID"] <- str_replace(string = df_all_jxn_pairings[ind, "novel_junID"]$novel_junID, 
                                                           pattern = "\\*", 
                                                           replacement = df_all_jxn_pairings[ind, "novel_strand"]$novel_strand  %>% as.character())
    if (any(str_detect(df_all_jxn_pairings$novel_junID, pattern = "\\*")) ) {
      print("ERROR!")
      break;
    }
  }
  
  
  ##########################################
  ## LOAD NEVER MIS-SPLICED JUNCTIONS
  #########################################
  
  message("Obtaining never mis-spliced junctions ...")
  
  df_never_misspliced <- get_all_intron_never_misspliced(recount3.project.IDs = recount3.project.IDs,
                                                         all.clusters = all.clusters,
                                                         project.name = project.name,
                                                         database.folder,
                                                         results.folder)
  
  if ( any(str_detect(df_never_misspliced$ref_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still contain a '*' in their IDs")
    break;
  }
  
  ## Remove the introns paired with novel junctions (i.e. mis-spliced)
  df_never_misspliced_tidy <- df_never_misspliced %>%
    dplyr::filter(!(ref_junID %in% df_all_jxn_pairings$ref_junID)) %>%
    as_tibble()
  
  df_never_misspliced_tidy %>% distinct(ref_junID) %>% as_tibble()
  
  
  ############################################
  ## GET ALL JUNCTIONS THAT HAVE NOT BEEN PAIRED
  ############################################ 
  
  message("Obtained all split reads that have not been paired ...")
  
  ## These are all the non-paired junctions, including the never mis-spliced 
  df_not_paired <- all_split_reads_details_qc_level2 %>%
    dplyr::filter(!(junID %in% c(df_all_jxn_pairings$ref_junID,
                                 df_all_jxn_pairings$novel_junID)))
  
  ## Hence, this should equal to the number of never-misspliced
  if ( ! identical(intersect(df_not_paired$junID, 
                             df_never_misspliced$ref_junID) %>% sort(), 
                   df_never_misspliced_tidy %>% distinct(ref_junID) %>% pull(ref_junID) %>% sort()) ) {
    print("ERROR! Some of the never mis-spliced junctions have not been found as not paired")
    break;
  }
  
  if (any(str_detect(df_never_misspliced_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(df_not_paired$junID, pattern = "\\*"))) {
    print("ERROR!")
  }
  
  ## All never mis-spliced should be categorised as not paired.
  ## Thus, this should be zero
  stopifnot(
    "Some never mis-spliced junctions have been paired!" =
      all(setdiff(df_never_misspliced_tidy$ref_junID, 
                  df_not_paired %>% dplyr::filter(type == "annotated") %>% pull(junID)) %>% length() == 0)
  )
  
  ## Separate the never mis-spliced from the non-paired junctions
  df_not_paired_tidy <- df_not_paired %>%
    dplyr::filter(!(junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    distinct(junID, .keep_all = T) %>%
    as_tibble()
  
  df_not_paired_tidy %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  ## This is the final number of non-paired junctions
  df_not_paired_tidy %>% distinct(junID)
  
  
  ## Hence, this should be zero
  if ( intersect(df_not_paired_tidy$junID, 
                 df_all_jxn_pairings$novel_junID) %>% length() > 0 ) {
    print("ERROR!")
  }
  
  
  df_all_jxn_pairings %>% distinct(novel_junID) %>% nrow() +
    df_all_jxn_pairings %>% distinct(ref_junID) %>% nrow() +
    df_never_misspliced_tidy %>% distinct(ref_junID) %>% nrow()
  
  ##########################################
  ## Remove ambiguous junctions
  ##########################################
  
  message("Removing ambiguous junctions ...")
  
  ## All these should be zero
  if( intersect(df_not_paired_tidy$junID, df_all_jxn_pairings$novel_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_all_jxn_pairings$ref_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_never_misspliced_tidy$ref_junID) %>% length() > 0) {
    print("ERROR!")
  }
  
  
  
  ## 1. Obtain the ambiguous junctions
  df_ambiguous_novel <- df_all_jxn_pairings %>%
    dplyr::filter(!(novel_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    dplyr::group_by(novel_junID) %>%
    mutate(distances_sd = distance %>% sd()) %>%
    dplyr::filter(distances_sd > 0)
  
  ## This is the number of ambiguous novel junctions to remove
  df_ambiguous_novel %>% ungroup() %>% distinct(novel_junID)

  
  
  ## 2. Remove ambiguous junctions
  df_all_jxn_pairings_tidy <- df_all_jxn_pairings %>%
    dplyr::filter(!(novel_junID %in% df_ambiguous_novel$novel_junID)) %>%
    distinct(novel_junID, ref_junID, .keep_all = T) %>%
    mutate(ref_strand = ref_strand %>% as.character(),
           novel_strand = novel_strand %>% as.character()) 
  
  
  
  ## Introns may parent multiple novel junctions. Hence, the introns thar are left orphaned after
  ## removing the ambiguous novel junctions are:
  (df_ambiguous_novel %>% 
      ungroup() %>%
      distinct(ref_junID) %>% nrow()) - (intersect(c(df_all_jxn_pairings_tidy$novel_junID,
                                                     df_all_jxn_pairings_tidy$ref_junID),
                                                   df_ambiguous_novel %>% 
                                                     ungroup() %>%
                                                     distinct(ref_junID) %>% pull) %>% length())
  
  ## 3. Get ambiguous figures and stats
  
  
  ## This is the number of unique novel junctions to be stored in the DB
  df_all_jxn_pairings_tidy %>%
    dplyr::distinct(novel_junID)
  
  # This is the number of unique introns to be stored inth DB
  df_all_jxn_pairings_tidy %>%
    dplyr::distinct(ref_junID)
  
  ## If we include the number of never mis-spliced junctions, the final number
  ## of junctions to be stored within the DB is:
  (df_all_jxn_pairings_tidy %>%
      dplyr::distinct(novel_junID) %>% 
      nrow()) + 
    (df_all_jxn_pairings_tidy %>%
       dplyr::distinct(ref_junID) %>% 
       nrow()) + 
    (df_never_misspliced_tidy %>% 
       distinct(ref_junID) %>% 
       nrow())
  
  
  
  
  
  
  
  
  
  ##############################################################################
  ## SAVE FINAL OBJECT
  ##############################################################################
  
  message("Saving results ...")
  
  ## 1. DISTANCES PAIRINGS
  
  if (any(str_detect(string = df_all_jxn_pairings_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_all_jxn_pairings_tidy$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a '*' within their IDs!")
  }
  df_all_jxn_pairings_tidy <- df_all_jxn_pairings_tidy %>%
    inner_join(y = all_split_reads_details_qc_level2 %>% 
                 dplyr::select(junID, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  saveRDS(object = df_all_jxn_pairings_tidy,
          file = paste0(database.folder, "/all_jxn_correct_pairings.rds"))
  
  
  ## 2. NEVER MIS-SPLICED
  
  if (any(str_detect(string = df_never_misspliced_tidy$ref_junID, pattern = "\\*")) ) {
    print("ERROR! Some NEVER MIS-SPLICED junctions still have a '*' in their IDs!")
  }
  df_never_misspliced_tidy <- df_never_misspliced_tidy %>%
    inner_join(y = all_split_reads_details_qc_level2 %>% 
                 dplyr::select(junID, seqnames, start, end, width, strand, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  
  saveRDS(object = df_never_misspliced_tidy,
          file = paste0(database.folder, "/all_jxn_never_misspliced.rds"))
  
  
  ## 3. AMBIGUOUS JUNCTIONS
  
  if (any(str_detect(string = df_ambiguous_novel$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_ambiguous_novel$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  saveRDS(df_ambiguous_novel,
          file = paste0(database.folder, "/all_jxn_ambiguous_pairings.rds"))
}
