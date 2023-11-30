
#' Title
#' Creates the 'intron' and 'novel' master tables
#' @param database.path 
#' @param gtf.version 
#' @param database.folder 
#' @param results.folder 
#'
#' @return
#' @export
#'
#' @examples
sql_create_master_tables <- function(database.path,
                                     gtf.version,
                                     database.folder,
                                     results.folder) {
  
  
  message(Sys.time(), " - loading GRCh38 reference...")
  
  hg38 <- rtracklayer::import(con = paste0(dependencies_folder,
                                           "/Homo_sapiens.GRCh38.",gtf.version,".chr.gtf"))
  
  ##########################################
  ## LOAD AND TIDY THE PAIR-WISE DISTANCES
  ##########################################
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  
  if ( file.exists(paste0(database.folder, "/all_jxn_correct_pairings.rds")) ) {
    
    message(Sys.time(), " - loading the pre-generated pair-wise distances data...")
    
    df_all_distances_pairings <- readRDS(file = paste0(database.folder, "/all_jxn_correct_pairings.rds"))
    
    df_ambiguous_novel <- readRDS(file = paste0(database.folder, "/all_jxn_ambiguous_pairings.rds"))
    
    df_introns_never <- readRDS(file = paste0(database.folder, "/all_jxn_never_misspliced.rds"))
    
  } else {
    
    print("ERROR loading file dependencies!")
    break;
    
  }
  
  ####################################
  ## A) GET ALL INTRONS
  ####################################
  
  
  ## Get all annotated introns
  print(paste0(Sys.time(), " - getting the mis-spliced introns..."))
  
  df_all_introns <- df_all_distances_pairings %>%
    distinct(ref_junID, .keep_all = T) %>%
    as_tibble() %>%
    mutate(misspliced = T) %>%
    dplyr::select(ref_junID,seqnames = ref_seq,
                  start = ref_start,
                  strand = ref_strand,
                  end = ref_end,
                  gene_id, 
                  tx_id_junction,
                  misspliced) %>%
    GRanges() %>% ## This is to get the width() as calculated by GRanges
    as_tibble()
  
  
  print(paste0(Sys.time(), " - getting the never mis-spliced introns..."))
  
  ## Remove potential * in the junID of the reference introns
  if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
    print("ERROR! some never mis-spliced junctions still have an *!")
    break;
  }
  
  df_introns_never_tidy <- df_introns_never %>% 
    as_tibble() %>%
    mutate(misspliced = F) %>%
    dplyr::filter(!(ref_junID %in% df_all_introns$ref_junID)) %>%
    dplyr::filter(!(ref_junID %in% df_ambiguous_novel$ref_junID)) 
  
  
  if (any(df_introns_never_tidy$width %>% abs() < 25)) {
    print("ERROR! some mis-spliced introns are shorter than 25bp!")
  }
  if (any(df_introns_never_tidy$width %>% abs() < 25)) {
    print("ERROR! some never mis-spliced introns are shorter than 25bp!")
  }
  if (intersect(df_introns_never_tidy$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    print("ERROR! some never mis-spliced introns are mis-spliced!")
  }
  
  
  ## Merge mis-spliced introns and never mis-spliced introns
  df_introns_introverse <- rbind(df_introns_never_tidy, df_all_introns)
  if (any(str_detect(string = df_introns_introverse$ref_junID, pattern = "\\*" ))) {
    print("ERROR! Some introns still have an ambiguous '*' strand!")
  }
  if (any(df_introns_introverse$ref_junID %>% duplicated())) {
    print("ERROR! Duplicated IDs")
  }
  
  #saveRDS(object = df_introns_introverse %>%
  #          distinct(ref_junID, .keep_all = T),
  #        file = paste0(folder_database, "/df_all_introns_database.rds"))
  
  
  
  ######################################
  ## INTRONS - REMOVE AMBIGUOUS INTRONS
  ######################################
  
  print(paste0(Sys.time(), " - getting the ambiguous introns..."))
  
  # df_introns_introverse <- df_introns_introverse %>%
  #   dplyr::filter(!(ref_junID %in% df_ambiguous_intron$junID)) %>%
  #   as_tibble()
  
  df_introns_introverse_tidy <- df_introns_introverse %>%
    distinct(ref_junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  ## There should not be any ambiguous intron at this point
  if ( any(df_introns_introverse_tidy %>% dplyr::filter(ambiguous == T) %>% nrow() > 0) ) {
    print("ERROR! Still there are some ambiguous introns")
    break;
  } else {
    df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
      dplyr::select(-ambiguous)
  }
  
  message(df_introns_introverse_tidy %>%
            distinct(ref_junID) %>%
            nrow(), " annotated introns")
  
  #saveRDS(object = df_introns_introverse_tidy,
  #        file = paste0(folder_database, "/df_all_introns_database_tidy.rds"))
  
  ######################################
  ## GENES - CREATE GENE TABLE
  ######################################
  
  
  # 
  # df_protein <- readRDS(file = paste0(results.folder, "/all_genes_PC_biotype.rds")) %>%
  #   as_tibble()
  
  
  sql_create_master_table_gene(database.path = database.path,
                               hg38 = hg38,
                               # protein_biotype = df_protein,
                               gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id) )
  
  
  ######################################
  ## TX_JUNCTION - CREATE TX TABLE
  ######################################
  
  
  sql_create_master_table_transcript(database.path = database.path,
                                     gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id),
                                     hg38 = hg38,
                                     tx_ids = df_introns_introverse_tidy %>% unnest(tx_id_junction) %>% distinct(tx_id_junction))
  
  
  ############################################
  ## INTRONS - ADD THE TRANSCRIPT FOREING KEY
  ############################################
  
  print(paste0(Sys.time(), " --> adding TRANSCRIPT foreing key to the introns..."))
  
  query <- paste0("SELECT id, transcript_id FROM 'transcript'")
  df_transcripts <- dbGetQuery(con, query) %>% as_tibble()
  
  
  ## Add the GENE ID for the foreign key
  df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
    unnest(tx_id_junction) %>%
    inner_join(df_transcripts,
               by = c("tx_id_junction" = "transcript_id")) %>%
    dplyr::select(-gene_id, -tx_id_junction) %>%
    dplyr::rename(transcript_id = id) 
  
  
  message(df_introns_introverse_tidy %>%
            distinct(ref_junID) %>% nrow(), " introns to store")
  
  
  ######################################
  ## INTRONS - ADD MAXENTSCAN INFO 
  ######################################
  
  print(paste0(Sys.time(), " --> adding the MaxEntScan info ..."))
  
  wd <- getwd()
  if ( !file.exists(paste0(dependencies_folder, 
                           "/Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
    print(paste0("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
                 does not exist within the specified dependencies folder."))
    break;
  }
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_introns_introverse_tidy %>% dplyr::rename(junID = ref_junID) %>% distinct(junID, .keep_all = T),
                                                 max_ent_tool_path = paste0(dependencies_folder, "/fordownload/"),
                                                 homo_sapiens_fasta_path = paste0(dependencies_folder, 
                                                                                  "/Homo_sapiens.GRCh38.dna.primary_assembly.fa") )
  rm(df_introns_introverse_tidy)
  gc()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-donorSeqStart,
                  -donorSeqStop,
                  -AcceptorSeqStart,
                  -AcceptorSeqStop) %>%
    dplyr::rename(ref_donor_sequence = donor_sequence,
                  ref_acceptor_sequence = acceptor_sequence)
  
  setwd(wd)
  
  df_all_introns <- all_split_reads_tidy %>%
    mutate(ref_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if ((setdiff(df_all_introns$junID, df_all_introns$ref_junID) %>% length()) > 0) {
    print("ERROR! ")
    break;
  }
  
  df_all_introns %>% as_tibble()
  
  df_all_introns <- df_all_introns %>%
    dplyr::select(-one_of("junID", "ref_ss5score", "ref_ss3score")) %>% 
    dplyr::rename(ref_mes5ss = ss5score, 
                  ref_mes3ss = ss3score) %>%
    dplyr::relocate(ref_junID) %>%
    dplyr::relocate(c(ref_mes5ss, ref_mes3ss), .before = transcript_id ) 
  
  
  df_all_introns %>% as_tibble()
  
  
  ################################################
  ## INTRONS - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  df_all_introns %>% as_tibble()
  
  print(paste0(Sys.time(), " - adding CDTS and Conservation scores..."))
  
  df_all_introns <- generate_cdts_phastcons_scores(db_introns = df_all_introns %>% 
                                                     distinct(ref_junID, .keep_all = T) %>%
                                                     as_tibble(),
                                                   intron_size = c(35,50,100),
                                                   phastcons_type = 17) %>%
    as_tibble()
  
  df_all_introns <- df_all_introns %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
  # df_all_introns <- df_all_introns %>% 
  #   dplyr::select(-c("CDTS_5ss_mean", "CDTS_3ss_mean", "phastCons20way_5ss_mean", "phastCons20way_3ss_mean",
  #                    "phastCons7way_5ss_mean", "phastCons7way_3ss_mean", "phastCons4way_5ss_mean", "phastCons4way_3ss_mean"))
  # 
  # df_all_introns <- df_all_introns %>%
  #   dplyr::rename(ref_CDTS_5ss = CDTS_5ss_mean,
  #                 ref_CDTS_3ss = CDTS_3ss_mean,
  #                 ref_phastCons20_5ss = phastCons20way_5ss_mean,
  #                 ref_phastCons20_3ss = phastCons20way_3ss_mean,
  #                 ref_phastCons7_5ss = phastCons7way_5ss_mean,
  #                 ref_phastCons7_3ss = phastCons7way_3ss_mean,
  #                 ref_phastCons4_5ss = phastCons4way_5ss_mean,
  #                 ref_phastCons4_3ss = phastCons4way_3ss_mean)
  
  df_all_introns %>% as_tibble()
  
  
  ######################################
  ## INTRONS - ADD THE CLINVAR DATA
  ######################################
  
  print(paste0(Sys.time(), " - adding the ClinVar data...")) 
  
  clinvar_tidy <- readRDS(file = paste0(dependencies_folder, "/clinvar_intronic_tidy.rds")) %>% 
    GRanges() %>% 
    diffloop::addchr()
  clinvar_tidy %>% head()
  
  df_all_introns_tidy <- df_all_introns %>%
    mutate(clinvar = F) %>%
    GRanges() 
  
  ## Find overlaps between clinvar mutations and mis-splicing ratios
  overlaps <- GenomicRanges::findOverlaps(query = clinvar_tidy,
                                          subject = df_all_introns_tidy,
                                          type = "any")
  
  
  elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar"] <- T
    # ifelse(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"] == "none",
    #        "-", paste0(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"], "-",
    #                    elementMetadata(clinvar_tidy)[queryHits(overlaps), "distance"], "bp"))
  #elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_start"] <- clinvar_tidy[queryHits(overlaps), ] %>% start()
  #elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_end"] <- clinvar_tidy[queryHits(overlaps), ] %>% end()
  
  message(df_all_introns_tidy %>%
            as_tibble() %>%
            filter(clinvar == T) %>%
            nrow(), " introns containing ClinVar variants")
  df_all_introns_tidy %>% head()
  
  
  ######################################
  ## INTRONS - ADD THE INTRON TYPE
  ## This intron type corresponds to whether the intron is 
  ## spliced out by the minor or the major spliceosome
  ######################################
  
  print(paste0(Sys.time(), " - adding the IAOD intron data...")) 
  
  ## Load intron type files
  u12_introns <- readRDS(file = paste0(dependencies_folder,
                                       "/minor_introns_tidy.rds")) %>%
    GRanges() %>%
    diffloop::addchr()
  
  ## Add two new columns to incorporate info about intron type
  df_all_introns_tidy <- df_all_introns_tidy %>%
    as_tibble() %>%
    mutate(u2_intron = T)
  
  
  ## MINOR INTRON
  print(paste0(Sys.time(), " - Getting junctions spliced out by the minor spliceosome."))
  overlaps <- GenomicRanges::findOverlaps(query = u12_introns,
                                          subject = df_all_introns_tidy %>% GRanges(),
                                          ignore.strand = FALSE,
                                          type = "equal")
  
  message(queryHits(overlaps) %>% length(), " introns spliced by the minor spliceosome!")
  df_all_introns_tidy[subjectHits(overlaps),]$u2_intron <- F
  
  
  df_all_introns_tidy <- df_all_introns_tidy %>%
    filter(u2_intron == T) %>%
    dplyr::select(-u2_intron)
  
  
  message(df_all_introns_tidy$ref_junID %>% unique %>% length(), " introns to be stored!")
  df_all_introns_tidy %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(misspliced)
  
  
  ######################################
  ## INTRONS - ADD THE TRANSCRIPT BIOTYPE
  ######################################
  
  
  df_protein_junID <- readRDS(file = paste0(results.folder, "/all_junID_PC_biotype.rds")) %>%
    as_tibble()
  
  df_all_introns_tidy <- df_all_introns_tidy %>%
    inner_join(y = df_protein_junID %>% dplyr::select(junID, protein_coding),
               by = c("ref_junID" = "junID"))
  
  ######################################
  ## INTRONS - POPULATE THE TABLE
  ######################################
  
  # mean_phastCons7way5ss_35 DOUBLE NOT NULL, 
  # mean_phastCons7way3ss_35 DOUBLE NOT NULL, 
  # mean_phastCons20way5ss_35 DOUBLE NOT NULL, 
  # mean_phastCons20way3ss_35 DOUBLE NOT NULL, 
  # mean_phastCons4way5ss_70 DOUBLE NOT NULL, 
  # mean_phastCons4way3ss_70 DOUBLE NOT NULL, 
  # mean_phastCons7way5ss_70 DOUBLE NOT NULL, 
  # mean_phastCons7way3ss_70 DOUBLE NOT NULL, 
  # mean_phastCons20way5ss_70 DOUBLE NOT NULL, 
  # mean_phastCons20way3ss_70 DOUBLE NOT NULL, 
  # mean_CDTS5ss_70 DOUBLE NOT NULL, 
  # mean_CDTS3ss_70 DOUBLE NOT NULL, 
  # mean_phastCons4way5ss_100 DOUBLE NOT NULL, 
  # mean_phastCons4way3ss_100 DOUBLE NOT NULL, 
  # mean_phastCons7way5ss_100 DOUBLE NOT NULL, 
  # mean_phastCons7way3ss_100 DOUBLE NOT NULL, 
  # mean_phastCons20way5ss_100 DOUBLE NOT NULL, 
  # mean_phastCons20way3ss_100 DOUBLE NOT NULL, 
  # mean_CDTS5ss_100 DOUBLE NOT NULL, 
  # mean_CDTS3ss_100 DOUBLE NOT NULL, 
  
  df_all_introns_tidy %>% names()
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'intron'",
                  "(ref_junID NUMERIC PRIMARY KEY NOT NULL,
                  ref_coordinates TEXT NOT NULL, 

                  seqnames NUMERIC NOT NULL,
                  start NUMERIC NOT NULL,
                  end NUMERIC NOT NULL,
                  strand TEXT NOT NULL, 

                  ref_length INTEGER NOT NULL, 
                  ref_mes5ss DOUBLE NOT NULL, 
                  ref_mes3ss DOUBLE NOT NULL, 

                  mean_phastCons17way5ss_100 DOUBLE NOT NULL, 
                  mean_phastCons17way3ss_100 DOUBLE NOT NULL, 
                  mean_phastCons17way5ss_50 DOUBLE NOT NULL, 
                  mean_phastCons17way3ss_50 DOUBLE NOT NULL, 
                  mean_phastCons17way5ss_35 DOUBLE NOT NULL, 
                  mean_phastCons17way3ss_35 DOUBLE NOT NULL, 
                  
                  mean_CDTS5ss_100 DOUBLE NOT NULL, 
                  mean_CDTS3ss_100 DOUBLE NOT NULL, 
                  mean_CDTS5ss_50 DOUBLE NOT NULL, 
                  mean_CDTS3ss_50 DOUBLE NOT NULL, 
                  mean_CDTS5ss_35 DOUBLE NOT NULL, 
                  mean_CDTS3ss_35 DOUBLE NOT NULL, 
                  
                  ref_donor_sequence TEXT NOT NULL,
                  ref_acceptor_sequence TEXT NOT NULL,

                  u2_intron BOOL,
                  clinvar BOOL NOT NULL, 
                  protein_coding DOUBLE NOT NULL, 
                  
                  misspliced BOOL NOT NULL,
                  transcript_id INTEGER NOT NULL,
                  FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database.path)
  dbListTables(con)
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  print("'Intron' table created!")
  
  
  ## POPULATE INTRON TABLE ----------------------------------------------------
  df_all_introns_tidy_final <- df_all_introns_tidy %>% 
    dplyr::rename(ref_length = width,
                  ref_coordinates = ref_junID) %>%
    distinct(ref_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("ref_junID")
  
  if (all(df_all_introns_tidy_final$misspliced == T)) {
    print("ERROR! all introns classified as misspliced!")
  }
  # if (any(is.na(df_all_introns_tidy_final$lncRNA))) {
  #   print("ERROR! some introns do not have a biotype assigned")
  # }
  if (any(df_all_introns_tidy_final$transcript_id %>% is.na())) {
    print("ERROR! some introns do not have a gene assigned")
  }
  
  DBI::dbAppendTable(conn = con,
                     name = "intron", 
                     value = df_all_introns_tidy_final)
  
  message("'Intron' master table populated! ", 
               df_all_introns_tidy_final %>% distinct(ref_junID) %>% nrow(), " annotated introns stored!" )
  
  ## CREATE INDEX TO SPEED UP QUERIES ------------------------------------------
  query <- paste0("CREATE UNIQUE INDEX 'index_intron' ON 'intron'(ref_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  ########################################
  ## B) NOVEL JUNCTION - ADD MAXENTSCAN INFO 
  ########################################
  
  print(paste0(Sys.time(), " - adding MaxEntScan scores to the NOVEL JUNCTONS..."))
  
  df_all_novel_raw_tidy <- df_all_distances_pairings %>%
    mutate( start = novel_start %>% as.integer(),
            end = novel_end %>% as.integer()) %>%
    dplyr::select(seqnames = novel_seq,
                  start, end, 
                  strand = novel_strand,
                  novel_junID, ref_junID, 
                  novel_type = type, distance) %>%
    distinct(novel_junID, .keep_all = T) %>%
    filter(ref_junID %in% df_all_introns_tidy_final$ref_coordinates)
  
  wd <- getwd()
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_all_novel_raw_tidy %>% dplyr::rename(junID = novel_junID) %>% distinct(junID, .keep_all = T),
                                                 max_ent_tool_path = paste0(dependencies_folder,"/fordownload/"),
                                                 homo_sapiens_fasta_path = paste0(dependencies_folder,
                                                                                  "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
  
  all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-donorSeqStart,
                  -donorSeqStop,
                  -AcceptorSeqStart,
                  -AcceptorSeqStop) %>%
    dplyr::rename(novel_donor_sequence = donor_sequence,
                  novel_acceptor_sequence = acceptor_sequence)
  
  setwd(wd)
  
  
  df_all_novels_tidy <- all_split_reads_tidy %>%
    mutate(novel_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if ((setdiff(df_all_novels_tidy$junID, df_all_novels_tidy$novel_junID) %>% length()) > 0) {
    print("ERROR! Novel junctions have been analysed under different sorting.")
    break;
  }
  
  rm(all_split_reads_tidy)
  gc()
  
  df_all_novels_tidy <- df_all_novels_tidy %>%
    dplyr::select(-one_of("junID","novel_ss5score","novel_ss3score")) %>% 
    dplyr::rename(novel_mes5ss = ss5score, 
                  novel_mes3ss = ss3score) %>%
    dplyr::relocate(ref_junID, novel_junID) %>%
    dplyr::relocate(c(novel_mes5ss, novel_mes3ss), .before = novel_type ) %>% 
    as_tibble()
  
  ################################################
  ## NOVEL - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  # df_all_novels_tidy %>% as_tibble()
  # 
  # print(paste0(Sys.time(), " - adding CDTS and Conservation scores..."))
  # 
  # df_all_novels_tidy <- generate_cdts_phastcons_scores(db_introns = df_all_novels_tidy %>% 
  #                                                        distinct(novel_junID, .keep_all = T) %>%
  #                                                        as_tibble()) %>%
  #   as_tibble()
  # 
  # # df_all_novels_tidy <- df_all_novels_tidy %>%
  # #   dplyr::rename(novel_CDTS5score = CDTS_5ss_mean,
  # #                 novel_CDTS3score = CDTS_3ss_mean,
  # #                 novel_cons5score = phastCons20way_5ss_mean,
  # #                 novel_cons3score = phastCons20way_3ss_mean)
  # 
  # df_all_novels_tidy %>% as_tibble()
  
  ####################################
  ## NOVEL - ADD INTRON FOREIGN KEY REFERENCE
  ####################################
  
  message(Sys.time(), " - adding the INTRON foreign key to the NOVEL JUNCTONS...")
  
  df_all_novels_tidy_final <- df_all_novels_tidy %>% 
    inner_join(y = df_all_introns_tidy_final %>%
                 dplyr::select(ref_junID, ref_coordinates),
               by = c("ref_junID" = "ref_coordinates" )) %>%
    dplyr::select(-ref_junID) %>%
    dplyr::rename(ref_junID = ref_junID.y) %>% 
    as_tibble() %>%
    dplyr::relocate(ref_junID)
  
  rm(df_all_novels_tidy)
  
  df_all_novels_tidy_final
  ####################################
  ## CREATE NOVEL JUNCTION TABLE
  ####################################
  
  # dbRemoveTable(conn = con, "novel")
  query <- paste0("CREATE TABLE IF NOT EXISTS 'novel'",
                  "(novel_junID NUMERIC NOT NULL,
                  ref_junID NUMERIC NOT NULL,

                  seqnames NUMERIC NOT NULL,
                  start NUMERIC NOT NULL,
                  end NUMERIC NOT NULL,
                  strand TEXT NOT NULL, 

                  novel_coordinates TEXT NOT NULL, 
                  novel_mes5ss DOUBLE NOT NULL, 
                  novel_mes3ss DOUBLE NOT NULL,

                  novel_donor_sequence TEXT NOT NULL,
                  novel_acceptor_sequence TEXT NOT NULL,
              
                  novel_length INTEGER NOT NULL, 
                  novel_type TEXT NOT NULL, 
                  distance INTEGER NOT NULL,
                  PRIMARY KEY (ref_junID, novel_junID),
                  FOREIGN KEY (ref_junID) REFERENCES 'intron'(ref_junID))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  message("'Novel' master table created!")
  
  df_all_novels_tidy_final <- df_all_novels_tidy_final %>% 
    as_tibble() %>%
    dplyr::rename(novel_coordinates = novel_junID ) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    GRanges() %>%
    as_tibble() %>%
    tibble::rowid_to_column("novel_junID") %>%
    dplyr::rename(novel_length = width)
  
  df_all_novels_tidy_final
  
  if (any(duplicated(df_all_novels_tidy_final$novel_coordinates))) {
    print("ERROR! some novel junctions are duplicated")
    break;
  }
  
  DBI::dbAppendTable(conn = con,
                     name = "novel", 
                     value = df_all_novels_tidy_final)
  
  message("'Novel' master table populated! ", 
          df_all_novels_tidy_final %>% distinct(novel_junID) %>% nrow(), 
          " novel junctions stored!" )
  
  DBI::dbDisconnect(conn = con)
  
}




