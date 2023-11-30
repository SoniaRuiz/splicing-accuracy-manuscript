#' Title
#' Create master table 'transcript'
#' @param database.path Local path to the .sqlite database
#' @param gene_ids List of gene Ensembl IDs
#' @param hg38 Human genome. Needed to extract info about the transcripts to be stored.
#' @param tx_ids List of transcript Ensembl IDs
#'
#' @return
#' @export
#'
#' @examples
sql_create_master_table_transcript  <- function(database.path,
                                                gene_ids,
                                                hg38,
                                                tx_ids) {
  
  
  ## Load MANE transcripts
  hg_mane_transcripts <- rtracklayer::import(con = paste0(dependencies_folder,
                                              "/MANE.GRCh38.v1.0.ensembl_genomic.gtf")) %>%
    as_tibble() %>%
    dplyr::select(-source, -score, -phase, -gene_id, -gene_type, -tag, -protein_id,
                  -db_xref,-transcript_type,-exon_id,-exon_number, -width ) %>%
    mutate(transcript_id = transcript_id %>% str_sub(start = 1, end = 15)) %>%
    drop_na() %>%
    dplyr::filter(type == "transcript") %>%
    distinct(transcript_id) %>%
    mutate(MANE = T)

  
  ## Create 'Transcript' MANE data
  
  message(Sys.time(), " --> creating 'transcript' master table...")
  con <- dbConnect(RSQLite::SQLite(), database.path)
  
  ## GET TRANSCRIPT ID
  hg38_transcripts_gene <- hg38 %>%
    as_tibble() %>%
    dplyr::filter(transcript_id %in% tx_ids$tx_id_junction,
                  type == "transcript") %>%
    mutate(gene_id = str_sub(gene_id, start = 1, end = 15)) %>%
    #distinct(gene_id, .keep_all = T) %>%
    dplyr::filter(gene_id %in% gene_ids$gene_id) %>%
    dplyr::select(transcript_id, TSL = transcript_support_level, gene_id) %>%
    mutate(TSL = str_sub(TSL, start = 1, end = 2) %>% as.integer())
  
  
  hg38_transcripts_gene[is.na(hg38_transcripts_gene[,"TSL"]),"TSL"] <- 10
  
  ## ADD MANE INFO
  hg38_transcripts_gene_mane <- hg38_transcripts_gene %>%
    left_join(y = hg_mane_transcripts,
              by = "transcript_id" )
  hg38_transcripts_gene_mane[is.na(hg38_transcripts_gene_mane[,"MANE"]),"MANE"] <- F
  
  
  
  ## ADD THE GENE FOREING KEY
  print(paste0(Sys.time(), " --> adding GENE foreing key to the transcripts..."))
  
  query <- paste0("SELECT id, gene_id FROM 'gene'")
  df_genes <- dbGetQuery(con, query) %>% as_tibble()
  
  ## Add the GENE ID for the foreign key
  hg38_transcripts_gene_mane <- hg38_transcripts_gene_mane %>%
    inner_join(df_genes %>% dplyr::select(id, gene_id),
               by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    dplyr::rename(gene_id = id) 
  
  ## Add the Id to the table
  hg38_transcripts_final <- hg38_transcripts_gene_mane %>% 
    as_tibble() %>%
    tibble::rowid_to_column("id") 
  
  
  ####################################
  ## CREATE TRANSCRIPT TABLE AND POPULATE
  ####################################
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'transcript'",
                  "(id NUMERIC NOT NULL,
                  transcript_id TEXT NOT NULL,
                  TSL NUMERIC NOT NULL, 
                  MANE BOOL NOT NULL, 
                  gene_id INTEGER NOT NULL,
                  PRIMARY KEY (id),
                  FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  print(paste0(Sys.time(), " --> 'Transcript' table created!"))
  
  
  if ( any(duplicated(hg38_transcripts_final$transcript_id)) ) {
    print("ERROR! some novel junctions are duplicated")
  }
  DBI::dbAppendTable(conn = con,
                     name = "transcript", 
                     value = hg38_transcripts_final)
  
  
  print( paste0(Sys.time(), " --> 'Transcript' table populated!") )
  print(paste0(Sys.time(), " --> ", hg38_transcripts_final %>% nrow(), " transcripts stored!"))
  
  DBI::dbDisconnect(con)
}