#' Title
#' Per junction, this function calculates the percentage of transcripts in which this junction may appear, that are protein-coding
#' @param recount3.project.IDs List of recount3 projects to analyse
#' @param project.name Name given to the project 
#' @param gtf.version Version of the reference transcriptome to use. In this case it has been used '105' corresponding
#' to Ensembl v105
#' @param database.folder Path to the local folder that stores the database to be produced and the files needed to produce it 
#' @param results.folder Local path to the folder that contains the results of the analyses performed
#'
#' @return
#' @export
#'
#' @examples
generate_transcript_biotype_percentage <- function(recount3.project.IDs,
                                                   project.name,
                                                   gtf.version,
                                                   database.folder,
                                                   results.folder) {
  
  
  #######################################
  ## GET THE TRANSCRIPT BIOTYPE
  #######################################
  
  print(paste0(Sys.time(), " - loading the human reference transcriptome ... "))
  
  ## Import HUMAN REFERENCE transcriptome
  homo_sapiens_v105 <- rtracklayer::import(con = paste0(dependencies_folder, 
                                                        "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf")) %>% 
    as_tibble()
  
  ## Get v105 transcripts
  transcripts_v105 <- homo_sapiens_v105 %>%
    filter(type == "transcript") %>% 
    dplyr::select(transcript_id, transcript_biotype, gene_id)
  
  transcripts_v105 %>% head()
  
  
  #######################################
  ## LOAD ALL SPLIT READS 
  #######################################
  
  print(paste0(Sys.time(), " - loading the recount3 split reads ... "))
  
  ## LOAD the all split reads from all recount3 GTEx projects
  all_split_reads_details_all_tissues <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level2.rds") )
  
  all_split_reads_details_all_tissues %>% head()
  all_split_reads_details_all_tissues %>% nrow()
  
  
  
  
  #######################################
  ## EXPLORE AND TIDY THE RESULT
  #######################################
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_all_tissues$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_all_tissues[ind, "junID"] <- 
      str_replace(string = all_split_reads_details_all_tissues[ind, "junID"]$junID, 
                  pattern = "\\*", 
                  replacement = all_split_reads_details_all_tissues[ind, "strand"]$strand %>% as.character() )
  }
  
  print(object.size(all_split_reads_details_all_tissues), units = "Gb")
  
  ## Merge datasets to add transcript biotype
  print(paste0(Sys.time(), " - adding reference transcript biotype to the junctions paired..."))
  
  ## Merging using data table structures saves time
  df_all_junctions <- all_split_reads_details_all_tissues %>% unnest(tx_id_junction) %>% data.table::as.data.table()
  transcripts_v105 <- transcripts_v105 %>% data.table::as.data.table()
  
  df_all_junctions %>% head()
  transcripts_v105 %>% head()
  
  df_all_junctions_tx <- df_all_junctions %>% 
    left_join(y = transcripts_v105,
              by = c("tx_id_junction" = "transcript_id"))
  
  # print(object.size(df_all_junctions_tx), units = "Gb")
  
  
  #######################################
  ## CALCULATE THE TRANSCRIPT PERCENTAGE
  #######################################
  
  
  print(paste0(Sys.time(), " - starting protein-coding percentage calculation ..."))
  print(paste0(Sys.time(), " - ", df_all_junctions$junID %>% unique() %>% length(), " total number of junctions."))
  
  
  ## Only keep non ambiguous jxn (i.e. junctions belonging to multiple genes)
  junID_OK <- df_all_junctions %>% 
    dplyr::group_by(junID) %>% 
    distinct(gene_id, .keep_all = T) %>% 
    dplyr::count() %>% 
    filter(n == 1) %>%
    pull(junID)
  
  # junID_OK %>% unique() %>% length()
  
  ## Calculate the biotype percentage
  df_all_percentage <- df_all_junctions_tx %>% 
    filter(junID %in% junID_OK) %>%
    dplyr::group_by(junID, transcript_biotype) %>%
    distinct(tx_id_junction, .keep_all = T) %>% 
    summarise(n = n()) %>% 
    mutate(percent = (n / sum(n)) * 100) %>%
    ungroup()
  
  df_all_percentage %>% head()
  saveRDS(object = df_all_percentage,
          file = paste0(results.folder, "/all_split_reads_qc_level2_biotype.rds"))
  
  ## Only filter by the protein-coding biotype
  df_all_percentage_tidy_PC <- df_all_percentage %>% 
    dplyr::group_by(junID) %>%
    rowwise() %>%
    mutate(percent = ifelse (transcript_biotype == "protein_coding", percent, 0)) %>%
    ungroup() %>%
    dplyr::group_by(junID) %>%
    filter(percent == max(percent)) %>%
    dplyr::select(-transcript_biotype, -n) %>%
    distinct(junID, .keep_all = T) %>%
    ungroup()
  
  ## Only filter by the lncRNA biotype
  df_all_percentage_tidy_lncRNA <- df_all_percentage %>% 
    dplyr::group_by(junID) %>%
    rowwise() %>%
    mutate(percent = ifelse (transcript_biotype == "lncRNA", percent, 0)) %>%
    ungroup() %>%
    dplyr::group_by(junID) %>%
    filter(percent == max(percent)) %>%
    dplyr::select(-transcript_biotype, -n) %>%
    distinct(junID, .keep_all = T) %>%
    ungroup()
  
  
  df_all_percentage_tidy_merged <- merge(x = df_all_percentage_tidy_PC %>% dplyr::rename(percent_PC = percent),
                                         y = df_all_percentage_tidy_lncRNA %>% dplyr::rename(percent_lncRNA = percent),
                                         by = "junID")
  df_all_percentage_tidy_merged %>% nrow()
  
  if (df_all_percentage_tidy_merged %>% filter(percent_PC == 100) %>% distinct(percent_lncRNA) %>% pull() != 0) {
    print("ERROR! some only protein-coding introns have been also classified as lncRNAs!")
  }
  if (df_all_percentage_tidy_merged %>% filter(percent_lncRNA == 100) %>% distinct(percent_PC) %>% pull() != 0) {
    print("ERROR! some only lncRNA introns have been also classified as protein-coding!")
  }
  
  print(object.size(df_all_percentage_tidy_merged), units = "Gb")
  
  saveRDS(object = df_all_percentage_tidy_merged,
          file = paste0(results.folder, "/all_split_reads_qc_level2_PC_biotype.rds"))
  print(paste0(Sys.time(), " - results saved!"))
  
  ##########################################
  ## FREE UP SOME MEMORY 
  ##########################################
  
  rm(all_split_reads_details_all_tissues)
  rm(df_all_junctions)
  rm(df_all_junctions_tx)
  rm(homo_sapiens_v105)
  rm(transcripts_v105)
  gc()
  
}