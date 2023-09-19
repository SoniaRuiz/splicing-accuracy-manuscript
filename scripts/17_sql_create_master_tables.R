

## LOAD DATA DEPENDENCIES ------------------------------------------------------------------------

## MANE transcripts
if ( !exists("hg_mane_transcripts") ) {
  
  print("Loading the 'hg_MANE' file...")
  hg_MANE <- rtracklayer::import(con = paste0(dependencies_folder,
                                              "/MANE.GRCh38.v1.0.ensembl_genomic.gtf"))
  hg_MANE_tidy <- hg_MANE %>%
    as_tibble() %>%
    dplyr::select(-source, -score, -phase, -gene_id, -gene_type, -tag, -protein_id,
                  -db_xref,-transcript_type,-exon_id,-exon_number, -width ) %>%
    mutate(transcript_id = transcript_id %>% str_sub(start = 1, end = 15)) %>%
    drop_na()
  
  hg_mane_transcripts <- hg_MANE_tidy %>%
    dplyr::filter(type == "transcript") %>%
    distinct(transcript_id) %>%
    mutate(MANE = T)
  
} else {
  print("'hg_mane_transcripts' file already loaded!")
}

# GET INFO FROM THE HG38
if ( !exists("hg38_transcripts") ) {
  
  print("Loading the 'hg38' file...")
  
  if ( !exists("hg38") ) {
    hg38 <- rtracklayer::import(con = paste0(dependencies_folder,
                                             "/Homo_sapiens.GRCh38.105.chr.gtf"))
  }
  
  hg38_transcripts <- hg38 %>%
    as_tibble() %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(transcript_id, transcript_support_level) %>%
    mutate(transcript_support_level = str_sub(string = transcript_support_level,
                                              start = 1,
                                              end = 1))
} else {
  print("'hg38_transcripts' file already loaded!")
}



## Create 'intron' and 'novel' master tables ----------------------------------------------------

sql_create_master_tables <- function(database.path,
                                     project.name,
                                     gtf.version,
                                     database.folder,
                                     results.folder) {
  
  
  ##########################################
  ## LOAD AND TIDY THE PAIR-WISE DISTANCES
  ##########################################
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  
  if ( file.exists(paste0(database.folder, "/all_jxn_correct_pairings.rds")) ) {
    
    print(paste0(Sys.time(), " - loading the pre-generated pair-wise distances data..."))
    
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
  
  #saveRDS(object = df_introns_introverse_tidy,
  #        file = paste0(folder_database, "/df_all_introns_database_tidy.rds"))
  
  ######################################
  ## GENES - CREATE GENE TABLE
  ######################################
  
  sql_create_master_table_gene(database.path = database.path,
                               hg38,
                               gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id) )
  
  
  ######################################
  ## TX_JUNCTION - CREATE TX TABLE
  ######################################
  
  
  sql_create_master_table_transcript(database.path = database.path,
                                     gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id),
                                     hg38,
                                     tx_ids = df_introns_introverse_tidy %>% unnest(tx_id_junction) %>% distinct(tx_id_junction))
  
  
  ######################################
  ## INTRONS - ADD THE TRANSCRIPT FOREING KEY
  ######################################
  
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
  
  
  df_introns_introverse_tidy %>%
    distinct(ref_junID)
  
  
  ######################################
  ## INTRONS - ADD MAXENTSCAN INFO 
  ######################################
  
  print(paste0(Sys.time(), " --> adding the MaxEntScan info ..."))
  
  wd <- getwd()
  if ( !file.exists(paste0(dependencies_folder, 
                           "/Homo_sapiens.GRCh38.dna.primary_assembly.fa")) ) {
    print(paste0("ERROR! File dependency 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
                 does not exist within the specified dependencies folder."))
  }
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_introns_introverse_tidy %>% dplyr::rename(junID = ref_junID),
                                                 max_ent_tool_path = paste0(dependencies_folder, "/fordownload/"),
                                                 homo_sapiens_fasta_path = paste0(dependencies_folder, 
                                                                                  "/Homo_sapiens.GRCh38.dna.primary_assembly.fa") )
  
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
    dplyr::select(-one_of("junID","ref_ss5score","ref_ss3score")) %>% 
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
                                                     as_tibble()) %>%
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
  ## INTRONS - ADD THE BIOTYPE INFO
  ######################################
  
  print(paste0(Sys.time(), " - adding the transcript biotype ... "))
  
  df_protein <- readRDS(file = paste0(results.folder, "/all_split_reads_qc_level2_PC_biotype.rds")) %>%
    as_tibble()
  
  print(paste0(Sys.time(), " - ", intersect(df_all_introns$ref_junID, df_protein$junID) %>% length(), " total number of introns found!."))
  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns."))
  
  df_all_introns <- df_all_introns %>%
    left_join(y = df_protein, by = c("ref_junID" = "junID")) %>%
    dplyr::rename(protein_coding = percent_PC,
                  lncRNA = percent_lncRNA)
  
  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns after adding transcript biotype."))
  
  ######################################
  ## INTRONS - ADD THE CLINVAR DATA
  ######################################
  
  print(paste0(Sys.time(), " - adding the ClinVar data...")) 
  
  clinvar_tidy <- readRDS(file = paste0(dependencies_folder, "/clinvar_intronic_tidy.rds")) %>%
    GRanges() 
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
  
  df_all_introns_tidy %>% head()
  
  
  ######################################
  ## INTRONS - ADD THE INTRON TYPE
  ## This intron type corresponds to whether the intron is spliced out by the minor or the major spliceosome
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
  overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u12_introns %>% seqnames(),
                                                                         ranges = IRanges(start = u12_introns %>% start(), 
                                                                                          end = u12_introns %>% end()),
                                                                         strand = u12_introns %>% strand()),
                                          subject = GenomicRanges::GRanges(seqnames = df_all_introns_tidy$seqnames,
                                                                           ranges = IRanges(start = df_all_introns_tidy$start, 
                                                                                            end = df_all_introns_tidy$end),
                                                                           strand = df_all_introns_tidy$strand),
                                          ignore.strand = FALSE,
                                          type = "equal")
  
  queryHits(overlaps) %>% length()
  df_all_introns_tidy[subjectHits(overlaps),]$u2_intron <- F
  
  
  
  
  ## MAJOR INTRON
  #print(paste0(Sys.time(), " - Getting junctions spliced out by the major spliceosome."))
  
  #overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u2_introns$seqnames %>% as.character(),
  #                                                                       ranges = IRanges(start = u2_introns$start,
  #                                                                                        end = u2_introns$end),
  #                                                                       strand = u2_introns$strand),
  #                                        subject = GenomicRanges::GRanges(seqnames = df_all_introns_tidy$seqnames,
  #                                                                         ranges = IRanges(start = df_all_introns_tidy$start,
  #                                                                                          end = df_all_introns_tidy$end),
  #                                                                         strand = df_all_introns_tidy$strand),
  #                                        ignore.strand = FALSE,
  #                                        type = "equal")
  
  # print(overlaps)
  #subjectHits(overlaps) %>% length()
  #df_all_introns_tidy[subjectHits(overlaps),]$u2_intron <- T
  
  
  df_all_introns_tidy$ref_junID %>% unique %>% length()
  df_all_introns_tidy %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(misspliced)
  
  
  ######################################
  ## INTRONS - POPULATE THE TABLE
  ######################################
  
  df_all_introns_tidy %>% names()
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'intron'",
                  "(ref_junID NUMERIC PRIMARY KEY NOT NULL,
                  ref_coordinates TEXT NOT NULL, 
                  ref_length INTEGER NOT NULL, 
                  ref_mes5ss DOUBLE NOT NULL, 
                  ref_mes3ss DOUBLE NOT NULL, 

                  mean_phastCons4way5ss_35 DOUBLE NOT NULL, 
                  mean_phastCons4way3ss_35 DOUBLE NOT NULL, 
                  mean_phastCons7way5ss_35 DOUBLE NOT NULL, 
                  mean_phastCons7way3ss_35 DOUBLE NOT NULL, 
                  mean_phastCons20way5ss_35 DOUBLE NOT NULL, 
                  mean_phastCons20way3ss_35 DOUBLE NOT NULL, 
                  mean_CDTS5ss_35 DOUBLE NOT NULL, 
                  mean_CDTS3ss_35 DOUBLE NOT NULL, 
                  mean_phastCons4way5ss_70 DOUBLE NOT NULL, 
                  mean_phastCons4way3ss_70 DOUBLE NOT NULL, 
                  mean_phastCons7way5ss_70 DOUBLE NOT NULL, 
                  mean_phastCons7way3ss_70 DOUBLE NOT NULL, 
                  mean_phastCons20way5ss_70 DOUBLE NOT NULL, 
                  mean_phastCons20way3ss_70 DOUBLE NOT NULL, 
                  mean_CDTS5ss_70 DOUBLE NOT NULL, 
                  mean_CDTS3ss_70 DOUBLE NOT NULL, 
                  mean_phastCons4way5ss_100 DOUBLE NOT NULL, 
                  mean_phastCons4way3ss_100 DOUBLE NOT NULL, 
                  mean_phastCons7way5ss_100 DOUBLE NOT NULL, 
                  mean_phastCons7way3ss_100 DOUBLE NOT NULL, 
                  mean_phastCons20way5ss_100 DOUBLE NOT NULL, 
                  mean_phastCons20way3ss_100 DOUBLE NOT NULL, 
                  mean_CDTS5ss_100 DOUBLE NOT NULL, 
                  mean_CDTS3ss_100 DOUBLE NOT NULL, 
                  
                  ref_donor_sequence TEXT NOT NULL,
                  ref_acceptor_sequence TEXT NOT NULL,

                  u2_intron BOOL,
                  clinvar BOOL NOT NULL, 
                  lncRNA INTEGER NOT NULL,
                  protein_coding INTEGER NOT NULL,
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
    #as_tibble() %>%
    #unnest(gene_id) %>%
    dplyr::rename(ref_length = width,
                  ref_coordinates = ref_junID) %>%
                  #ref_mes5ss = ref_ss5score, 
                  #ref_mes3ss = ref_ss3score,
                  #clinvar = clinvar_type) %>%
    distinct(ref_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("ref_junID") %>%
    dplyr::select(-c(seqnames, start, end, strand) ) 
                     #clinvar_start, clinvar_end) )
  
  if (all(df_all_introns_tidy_final$misspliced == T)) {
    print("ERROR! all introns classified as misspliced!")
  }
  if (any(is.na(df_all_introns_tidy_final$lncRNA))) {
    print("ERROR! some introns do not have a biotype assigned")
  }
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
    distinct(novel_junID, .keep_all = T)
  
  wd <- getwd()
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_all_novel_raw_tidy %>% dplyr::rename(junID = novel_junID),
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
  
  df_all_novels_tidy <- df_all_novels_tidy %>%
    dplyr::select(-one_of("junID","novel_ss5score","novel_ss3score")) %>% 
    dplyr::rename(novel_mes5ss = ss5score, 
                  novel_mes3ss = ss3score) %>%
    dplyr::relocate(ref_junID, novel_junID) %>%
    dplyr::relocate(c(novel_mes5ss, novel_mes3ss), .before = novel_type ) 
  
  
  df_all_novels_tidy <- df_all_novels_tidy %>% as_tibble()
  df_all_novels_tidy
  
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
  
  
  ####################################
  ## CREATE NOVEL JUNCTION TABLE
  ####################################
  
  # dbRemoveTable(conn = con, "novel")
  query <- paste0("CREATE TABLE IF NOT EXISTS 'novel'",
                  "(novel_junID NUMERIC NOT NULL,
                  ref_junID NUMERIC NOT NULL,

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
    dplyr::select(-seqnames, -start, -end, -strand) %>%
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




