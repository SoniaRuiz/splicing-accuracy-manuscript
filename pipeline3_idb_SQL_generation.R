###################################
## CONNECTION TO THE DB  
###################################

library(DBI)
library(tidyverse)
library(data.table)

# Create an ephemeral in-memory RSQLite database
# setwd("/home/sruiz/PROJECTS/splicing-project-app/introverse/")

gtf_version <- 105
###################################
## REMOVE ALL TABLES
###################################

remove_tables <- function(database_path,
                          all,
                          age = F) {
  
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  tables <- dbListTables(con)
  tables
  for (table in tables) {
    #if (str_detect(table %>% tolower(), pattern = "kidney")) {
    if (!all) {
      if (!(table %in% c("gene","intron", "mane", "master", "novel"))) {
        dbRemoveTable(conn = con, table)
        print(paste0("Table: '", table, "' has been removed!"))
      }
    } else {
      dbRemoveTable(conn = con, table)
      print(paste0("Table: '", table, "' has been removed!"))
    }
  }
  dbListTables(con)
}


###################################
## CREATE MASTER TABLE
###################################

create_master_table <- function(database_path,
                                age = F, QC = F)  {
  
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  if (age) {
    
    df_metadata <- map_df(SRA_projects, function(project) { 
      
      print(paste0(Sys.time(), " - getting info from ", project))
      
      df_project_metadata <- data.frame(cluster = c("60-79", "40-59", "20-39"),
                                        SRA_project = project) %>%
        return()
      
    })
    
  } else {
    
    df_metadata <- map_df(SRA_projects, function(project) { 
      
      print(paste0(Sys.time(), " - getting info from ", project))
      # project <- SRA_projects[1]
      # project <- SRA_projects[2]
      # project <- SRA_projects[6]
      metadata <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                        project, "/raw_data/samples_metadata.rds"))
      if (metadata %>% nrow() >= 1) {
        
        print(metadata$gtex.smtsd %>% unique())
        print(metadata %>% distinct(gtex.sampid) %>% nrow())
        
        df_project_metadata <- data.frame(age = metadata$gtex.age %>% as.character(),
                                          rin = metadata$gtex.smrin %>% as.character(),
                                          gender = metadata$gtex.sex %>% as.character(),
                                          tissue = metadata$gtex.smtsd,
                                          cluster = metadata$gtex.smtsd, #str_remove_all(metadata$gtex.smtsd, pattern = " ") %>% tolower(),
                                          cluster_tidy = metadata$gtex.smtsd,
                                          smnabtcht = metadata$gtex.smnabtcht, 
                                          sample_id = metadata %>% distinct(gtex.sampid) %>% nrow(),
                                          smafrze = metadata$gtex.smafrze,
                                          avg_read_length = metadata$recount_seq_qc.avg_len,
                                          mapped_read_count = metadata$recount_qc.star.all_mapped_reads,
                                          SRA_project_tidy = metadata$recount_project.project,
                                          SRA_project = metadata$recount_project.project) %>% 
          arrange(SRA_project, cluster) %>% 
          return()
        
      } else {
        
        return(NULL)
      }
      
    })
  }
  
  if (QC) {
    sex_tissues <- c("BREAST", "CERVIX_UTERI", "FALLOPIAN_TUBE", "OVARY", "PROSTATE", "TESTIS", "UTERUS", "VAGINA")
    
    df_metadata <- df_metadata %>% 
      dplyr::filter(!(SRA_project %in% sex_tissues),
             !(cluster %in% c("Brain - Cortex", "Brain - Cerebellum"))) %>%
      dplyr::filter(cluster %in% (df_metadata %>% dplyr::count(cluster) %>% dplyr::filter(n >= 70) %>% pull(cluster)))
  }
  
  
  # saveRDS(object = df_all_projects_metadata,
  #         file = "./dependencies/df_all_projects_metadata.rds")
  
  
  DBI::dbWriteTable(conn = con,
                    name = "master",
                    value = df_metadata,
                    overwrite = T)
  
  DBI::dbDisconnect(conn = con)
  
  print(paste0("Table: 'master' created!"))
}
# create_master_table()

###################################
## CREATE MANE TABLE
###################################

create_mane_table <- function(database_path,
                              age = F) {
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  ## LOAD THIS FILE
  #hg_MANE <- rtracklayer::import(con = "/data/references/MANE/MANE.GRCh38.v1.0.ensembl_genomic.gtf")
  
  hg_MANE_tidy <- hg_MANE %>%
    as_tibble() %>%
    dplyr::select(-source, -score, -phase, -gene_id, -gene_type, -tag, -protein_id, 
                  -db_xref,-transcript_type,-exon_id,-exon_number, -width ) %>%
    mutate(transcript_id = transcript_id %>% str_sub(start = 1, end = 15)) %>%
    drop_na()
  
  
  DBI::dbWriteTable(conn = con,
                    name = "mane",
                    value = hg_MANE_tidy,
                    overwrite = T)
  
  # dbRemoveTable(conn = con, "mane")
  DBI::dbDisconnect(conn = con)
  
  print(paste0("Table: 'mane' created!"))
}
# hg_MANE <- rtracklayer::import(con = "/data/references/MANE/MANE.GRCh38.v1.0.ensembl_genomic.gtf")
# create_mane_table()

###################################
## CREATE GENE TABLE
###################################

create_gene_table <- function(database_path,
                              age = F) {
  
  # DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  # dbRemoveTable(conn = con, "gene")
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  ## GET ALL DATA FROM MASTER
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  ## Loop through the cluster to obtain all the genes
  gene_ids <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(clusters, function(cluster) { 
      
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      
      if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                             cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))) {
        
        
        df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                            cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds")) %>% 
          as_tibble() %>%
          dplyr::select(gene_id)
        
        df_novel <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                          cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds")) %>% 
          as_tibble() %>%
          dplyr::select(gene_id)
        
        
        
        return( rbind(df_introns %>% dplyr::select(gene_id),
                      df_novel %>% dplyr::select(gene_id)) %>%
                  unnest(gene_id) %>%
                  distinct(gene_id))
        
      } else {
        
        return(NULL)
        
      }
    })
    
  })
  

  gene_ids <- gene_ids %>% distinct(gene_id)
  gene_ids %>% nrow()
  
  # hg38 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  # hg38 %>% as_tibble %>% distinct(gene_id, .keep_all = T) %>% nrow()
  
  hg38_transcripts <- hg38 %>%
    as_tibble() %>%
    #dplyr::select(gene_id, gene_name) %>%
    mutate(gene_id = str_sub(gene_id, start = 1, end = 15)) %>%
    #distinct(gene_id, .keep_all = T) %>%
    dplyr::filter(gene_id %in% gene_ids$gene_id) %>%
    dplyr::count(gene_id, type) %>%
    dplyr::filter(type == "transcript") %>%
    unnest(gene_id) %>%
    dplyr::select(-type) %>%
    dplyr::rename(n_transcripts = n)
  
  hg38_genes <- hg38 %>%
    as_tibble() %>%
    mutate(gene_id = str_sub(gene_id, start = 1,end = 15)) %>%
    dplyr::filter(gene_id %in% gene_ids$gene_id) %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name, gene_width = width)
  
  
  
  hg38_tidy <- hg38_transcripts %>%
    inner_join(hg38_genes, by = "gene_id" ) %>% 
    as_tibble()
  
  
  ## CREATE GENE_NAME TABLE ---------------------------------------------------
  sql_statement <- paste0("CREATE TABLE IF NOT EXISTS 'gene'", 
                                 "(id INTEGER PRIMARY KEY NOT NULL,
                                 gene_id TEXT NOT NULL,
                                 gene_name TEXT,
                                 n_transcripts INTEGER NOT NULL,
                          gene_width INTEGER NOT NULL)")
  res <- DBI::dbSendQuery(conn = con, statement = sql_statement)
  DBI::dbClearResult(res)
  
  #sql_statement <- paste0("CREATE UNIQUE INDEX gene_index ON gene(id)");
  #res <- DBI::dbSendQuery(conn = con, statement = sql_statement)
  #DBI::dbClearResult(res)
  
  
  ## POPULATE GENE_NAME TABLE -------------------------------------------------
  
  
  hg38_tidy <- hg38_tidy %>% 
    as_tibble() %>%
    tibble::rowid_to_column("id")
  
  print(Sys.time(), " - the gene set is identical to the hg: ",
        identical(hg38_tidy$gene_id %>% sort(),
                  gene_ids$gene_id %>% sort()))
  
  # hg38_genes %>%
  #   dplyr::filter(str_detect(gene_name, pattern = "c\\("))
  # any(str_detect(hg38_genes$gene_name, pattern = "c\\("))
  
  DBI::dbAppendTable(conn = con,
                     name = "gene", 
                     value = hg38_tidy)
  
  
  
  DBI::dbDisconnect(conn = con)
  print(paste0("Table: 'gene' created!"))
}
# hg38 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
# create_gene_table()

###################################
## CREATE INTRON TABLE
##################################



create_intron_table <- function(database_path,
                                age = F) {
  
  ## Connect to the DB
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  # dbRemoveTable(conn = con, "intron")
  
  ## GET ALL GENES
  query = paste0("SELECT * FROM 'gene'")
  all_genes <- dbGetQuery(con, query)
  
  # GET INFO FROM MASTER
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  # GET INFO FROM MANE
  query = paste0("SELECT * FROM 'mane'")
  mane_transcripts <- dbGetQuery(con, query) %>% 
    dplyr::filter(type == "transcript") %>%
    distinct(transcript_id) %>% pull()
  
  hg38_transcripts <- hg38 %>%
    as_tibble() %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(transcript_id, transcript_support_level) %>%
    mutate(transcript_support_level = str_sub(string = transcript_support_level,
                                              start = 1,
                                              end = 1))
  
  
  if (!exists("CNC_CDTS_CONS_gr")) {
    print("Loading the 'CNC_CDTS_CONS_gr' file...")
    aws.s3::s3load(object = "CNC_CDTS_CONS_gr.rda", bucket = "data-references", region = "eu-west-2")
    print("'CNC_CDTS_CONS_gr' file loaded!")
  } else {
    print("Loading the 'CNC_CDTS_CONS_gr' file...")
    aws.s3::s3load(object = "CNC_CDTS_CONS_gr.rda", bucket = "data-references", region = "eu-west-2")
    print("'CNC_CDTS_CONS_gr' file already loaded!")
  }
  
  ###################################
  ## CREATE INTRONS TABLE
  ###################################
  
  # Get all GTEx projects
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  ## LOOP THROUGH PROJECTS
  df_all_introns <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[6]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()

    map_df(clusters, function(cluster) {
      
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      
      if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                             cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))) {
        
        df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                            cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))
        

        
        df_introns_tidy <- df_introns %>% 
          as_tibble() %>%
          #dplyr::mutate(ref_mean_counts = round((ref_sum_counts / ref_n_individuals), digits = 2)) %>% 
          #dplyr::filter(u2_intron == T) %>%
          distinct(ref_junID, .keep_all = T) %>% 
          as_tibble() %>%
          dplyr::select(ref_junID,
                        seqnames,
                        start,
                        end,
                        strand,
                        ref_length = width,
                        #ref_ss5score, 
                        #ref_ss3score,
                        #ref_cons5score = phastCons20way_5ss_mean,
                        #ref_cons3score = phastCons20way_3ss_mean,
                        #ref_CDTS5score = CDTS_5ss_mean,
                        #ref_CDTS3score = CDTS_3ss_mean,
                        clinvar_type, 
                        u2_intron,
                        u12_intron,
                        #protein_coding,
                        tx_id_junction,
                        gene_id)
        
        return(df_introns_tidy)
        
      } else {
        return(NULL)
      }
    })  
  })  
 
  
  ######################################
  ## TIDY THE JUNCTION ID
  ######################################
  
  ## QC in case some coordinates were stored with a * instead of the strand
  df_all_introns <- df_all_introns %>%
    distinct(ref_junID, .keep_all = T) %>%
    mutate(strand = strand %>% as.character()) %>%
    rowwise() %>%
    mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                              str_replace(string = ref_junID, pattern = "\\*", strand ),
                              ref_junID))
  
  ######################################
  ## ADD MANE INFO - ROWISE
  ######################################
  
  print(paste0(Sys.time(), " --> adding MANE info ..."))
  
  df_all_introns <- df_all_introns %>%
    distinct(ref_junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(MANE = ifelse(any((tx_id_junction %>% unlist() %>% unname()) %in% mane_transcripts), T, F)) 
  
  
  ######################################
  ## ADD TSL INFO 
  ######################################
  
  print(paste0(Sys.time(), " --> adding TSL info ..."))
 
  ## Store a 10 value into the TSL field so it can be an integer field and reduce disk space
  df_all_introns <- df_all_introns %>%
    unnest(tx_id_junction) %>%
    dplyr::left_join(y = hg38_transcripts,
                     by = c("tx_id_junction" = "transcript_id")) %>%
    mutate(TSL = ifelse(is.na(transcript_support_level) | 
                          transcript_support_level == "N", 
                        10, transcript_support_level),
           TSL = TSL %>% as.integer())

  ## Remove the columns we do not need
  df_all_introns <- df_all_introns %>% 
    group_by(ref_junID) %>%
    mutate(TSL = TSL %>% min) %>%
    dplyr::select(-tx_id_junction,
                  -transcript_support_level) %>%
    distinct(ref_junID, .keep_all = T)
  
  
  ######################################
  ## ADD MAXENTSCAN INFO 
  ######################################
  
  print(paste0(Sys.time(), " --> adding the MaxEntScan info ..."))
  
  wd <- getwd()
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_all_introns %>% dplyr::rename(junID = ref_junID),
                                                 max_ent_tool_path = "/home/sruiz/fordownload/",
                                                 homo_sapiens_fasta_path = "/data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
  
  all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-donorSeqStart,
                  -donorSeqStop,
                  -AcceptorSeqStart,
                  -AcceptorSeqStop,
                  -donor_sequence,
                  -acceptor_sequence)
  
  setwd(wd)
  
  df_all_introns <- all_split_reads_tidy %>%
    mutate(ref_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if ((setdiff(df_all_introns$junID, df_all_introns$ref_junID) %>% length()) > 0) {
    print("ERROR!")
    break;
  }
  
  df_all_introns <- df_all_introns %>%
    dplyr::select(-one_of("junID","ref_ss5score","ref_ss3score")) %>% 
    dplyr::rename(ref_ss5score = ss5score, 
                  ref_ss3score = ss3score) %>%
    dplyr::relocate(ref_junID) %>%
    dplyr::relocate(c(ref_ss5score, ref_ss3score), .before = ref_length ) 
    
  
  df_all_introns %>% as_tibble()
  
  
  ######################################
  ## ADD THE CONSERVATION AND CDTS INFO
  ######################################
  
  df_all_introns %>% as_tibble()
  
  print(paste0(Sys.time(), " - adding CDTS and Conservation scores..."))
  
  df_all_introns <- add_cdts_cons_scores(db_introns = df_all_introns %>% as_tibble()) %>%
    as_tibble()

  df_all_introns <- df_all_introns %>%
    dplyr::rename(ref_CDTS5score = CDTS_5ss_mean,
                  ref_CDTS3score = CDTS_3ss_mean,
                  ref_cons5score = phastCons20way_5ss_mean,
                  ref_cons3score = phastCons20way_3ss_mean)
  
  ######################################
  ## ADD THE BIOTYPE INFO
  ######################################
  
  print(paste0(Sys.time(), " - adding the transcript biotype ... "))
  
  df_protein <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_annotated_SR_details_length_105_biotype.rds") %>%
    as_tibble()
  
  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns."))
  
  df_all_introns <- df_all_introns %>%
    left_join(y = df_protein,
              by = c("ref_junID" = "junID")) %>%
    dplyr::rename(protein_coding = percent_PC,
                  lncRNA = percent_lncRNA)

  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns after adding transcript biotype."))
  
  #######################################
  ## QC
  ######################################
  
  print(paste0(Sys.time(), " --> starting QC ..."))
  
  df_all_introns %>% nrow()
  df_all_introns %>% head()
  
  if ((df_all_introns %>% 
       unnest(gene_id) %>%
       group_by(ref_junID) %>% 
       distinct(gene_id, .keep_all = T) %>%
       dplyr::count(gene_id) %>%
       dplyr::filter(n > 1) %>% 
       nrow()) != 0) {
    print(paste0("Error! There are introns that have been assigned to multiple genes across tissues."))
  } 
  if ((df_all_introns %>% 
       group_by(ref_junID) %>% 
       distinct(MANE, .keep_all = T) %>%
       dplyr::count(MANE) %>%
       dplyr::filter(n > 1) %>% nrow()) != 0) {
    print(paste0("Error! There are introns that have been assigned multiple MANE values across tissues."))
  } 
  if ((df_all_introns %>%
    rowwise() %>%
    dplyr::filter(ref_junID %>% str_detect( pattern = "\\*")) %>%
    dplyr::select(ref_junID)) %>% nrow() > 1) {
    print(paste0("Error! There are introns with '*' strand."))
  }


  #######################################
  ## TIDY THE DATAFRAME
  ######################################
  
  print(paste0(Sys.time(), " --> tidying the dataframe ..."))
  

  ## Flatten each gene list element internally
  ## Only keep introns belonging to one single gene
  df_all_introns_tidy <- df_all_introns %>%
    distinct(ref_junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(gene = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  
  df_ambiguous <- df_all_introns_tidy %>%
    dplyr::filter(gene == T)

  
  ## If there are ambigous junctions, we discard them
  if (df_ambiguous %>% nrow() > 1) {
    
    print(paste0(Sys.time(), " --> there are ambiguous introns: ", df_ambiguous %>% nrow()))
    
    df_all_introns_tidy <- df_all_introns_tidy %>% 
      dplyr::filter(!(ref_junID %in% df_ambiguous$ref_junID))
    
    if (age) {
      file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/ambiguous_introns_age.rds"
    } else {
      file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/ambiguous_introns.rds"
    }
    saveRDS(object = df_ambiguous, file = file_name)
    
  }
  
  df_all_introns_tidy <- df_all_introns_tidy %>%
    dplyr::select(-gene) %>%
    mutate_if(is.list, simplify_all) %>%    
    unnest(gene_id)
  
  ## Add the GENE ID for the foreign key
  df_all_introns_tidy <- df_all_introns_tidy %>%
    inner_join(all_genes %>% data.table::as.data.table(),
               by = "gene_id") %>%
    dplyr::select(-gene_name, -gene_id, -n_transcripts, -gene_width) %>%
    dplyr::rename(gene_id = id)
  
  
 
  
  ## Extra QC
  ## Same introns should have been assigned same scores across tissues
  # df_all_introns %>%
  #   dplyr::group_by(ref_junID, gene_id) %>%
  #   distinct(ref_ss5score, .keep_all = T) %>% nrow()
  # df_all_introns %>%
  #   dplyr::group_by(ref_junID, gene_id) %>%
  #   distinct(ref_ss3score, .keep_all = T) %>% nrow()
  # df_all_introns %>%
  #   dplyr::group_by(ref_junID) %>%
  #   distinct(ref_ss3score, .keep_all = T) %>% nrow()
  # ## Each intron should only belong to the same gene
  # df_all_introns %>%
  #   dplyr::count(ref_junID, gene_id) 
  # 
  # which(duplicated(df_all_introns_tidy$ref_junID) == T)
  # df_all_introns_tidy %>%
  #   dplyr::filter(ref_junID == (df_all_introns_tidy[651,1]$ref_junID))
  
  #######################################
  ## CREATE INTRON TABLE
  ######################################
  
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'intron'",
  "(ref_junID NUMERIC PRIMARY KEY NOT NULL,
  ref_coordinates TEXT NOT NULL, 
  ref_length INTEGER NOT NULL, 
  ref_ss5score DOUBLE NOT NULL, 
  ref_ss3score DOUBLE NOT NULL, 
  ref_cons5score DOUBLE NOT NULL,
  ref_cons3score DOUBLE NOT NULL,
  ref_CDTS5score DOUBLE NOT NULL,
  ref_CDTS3score DOUBLE NOT NULL,
  u2_intron BOOL,
  u12_intron BOOL,
  clinvar TEXT NOT NULL, 
  MANE BOOL NOT NULL,
  TSL NUMERIC NOT NULL,
  lncRNA INTEGER NOT NULL,
  protein_coding INTEGER NOT NULL,
  gene_id INTEGER NOT NULL,
  FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  print("'Intron' table created!")
  
  
  ## POPULATE INTRON TABLE ----------------------------------------------------
  df_all_introns_tidy <- df_all_introns_tidy %>% 
    as_tibble() %>%
    mutate(ref_cons5score = ifelse(ref_cons5score %>% is.na(), 0, ref_cons5score),
           ref_cons3score = ifelse(ref_cons3score %>% is.na(), 0, ref_cons3score)) %>%
    dplyr::rename(ref_coordinates = ref_junID,
                  clinvar = clinvar_type) %>%
    distinct(ref_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("ref_junID")
  
  df_all_introns_tidy <- df_all_introns_tidy %>%
    dplyr::select(-c(seqnames, start, end, strand,  width) )
  

  if (any(df_all_introns_tidy$gene_id %>% is.na())) {
    print("ERROR! some introns do not have a gene assigned")
  }
  
  DBI::dbAppendTable(conn = con,
                     name = "intron", 
                     value = df_all_introns_tidy)
  
  print("'Intron' table populated!")
  ## CREATE INDEX TO SPEED UP QUERIES ------------------------------------------
  query <- paste0("CREATE UNIQUE INDEX 'index_intron' ON 'intron'(ref_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  
  DBI::dbDisconnect(conn = con)
}
# create_intron_table()

###################################
## CREATE NOVEL JUNCTION TABLE
##################################




create_novel_table <- function(database_path,
                               age = F) {
  
  ## Establish connection to the DB
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  # DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  # dbRemoveTable(conn = con, "novel")
  
  ## GET GENES
  query = paste0("SELECT id, gene_id FROM 'gene'")
  all_genes <- dbGetQuery(con, query)
  all_genes %>% nrow()
  
  ## GET REFERENCE INTRONS
  query = paste0("SELECT ref_junID, ref_coordinates FROM 'intron'")
  df_all_introns <- dbGetQuery(con, query) 
  df_all_introns %>% nrow()
  
  
  # GET MASTER INFO
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  
  ###################################
  ## GET NOVEL JUNCTION DATA 
  ###################################
  
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  gtf_version <- "105"
  
  
  df_all_novels <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[1]
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    df_all_novels <- map_df(clusters, function(cluster) {
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      
      if (file.exists( paste0(base_folder, "results/pipeline3/missplicing-ratio/",
                              cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))){
        
        df_novel <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/",
                                          cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))
        
        df_novel_tidy <- df_novel %>% 
          as_tibble() %>%
          dplyr::mutate(novel_type = novel_type %>% as.character()) %>%
          dplyr::select(ref_junID,
                        seqnames, start, end, strand,
                        novel_junID,
                        novel_ss5score, 
                        novel_ss3score, 
                        novel_type,
                        distance, 
                        gene_id)
        
        return(df_novel_tidy)
        
      } else {
        return(NULL)
      }
    })  
  })
  
  
  #############################################
  ## PREPARE DATA BEFORE DATABASE 
  #############################################
  
  
  
  df_all_novels_tidy <- df_all_novels %>%
    distinct(ref_junID, novel_junID, .keep_all = T) %>%
    rowwise() %>%
    mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                              str_replace(string = ref_junID, pattern = "\\*", strand ),
                              ref_junID))
  
  any(str_detect(string = df_all_novels_tidy$ref_junID, pattern = "\\*"))
  
  
  #############################################
  ## ADD INTRON ID info
  #############################################
  
  ## There are some introns that are ambiguous,thus novel junctions depending on them should be discarded.
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/ambiguous_introns.rds"
  df_ambiguous_introns <- readRDS(file = file_name)

  df_all_novels_tidy <- df_all_novels_tidy %>% 
    distinct(ref_junID, novel_junID, .keep_all = T) %>%
    data.table::as.data.table() %>%
    dplyr::filter(!(ref_junID %in% df_ambiguous_introns$ref_junID))
  
 
  
  ## This should be zero
  if (setdiff(df_all_novels_tidy$ref_junID, df_all_introns$ref_coordinates) %>% length() > 0) {
    print("ERROR!")
  }
  
  
  
  df_all_novels_tidy <- merge(x = df_all_novels_tidy,
                              y = df_all_introns,
                              by.x = "ref_junID",
                              by.y = "ref_coordinates",
                              all.x = T) %>%
    dplyr::select(-ref_junID) %>%
    dplyr::rename(ref_junID = ref_junID.y) %>% 
    as_tibble()
  
  
  
  # df_all_novels_tidy <- df_all_novels_tidy %>%
  #   inner_join(df_all_introns,
  #              by = c("ref_junID" = "ref_coordinates")) %>%
  #   dplyr::select(-ref_junID) %>%
  #   dplyr::rename(ref_junID = ref_junID.y) %>% 
  #   as_tibble()
  
  df_all_novels_tidy <- df_all_novels_tidy %>% 
    unnest(gene_id)
  

  df_all_novels_tidy %>% nrow()
  
  
  ######################################
  ## MERGE WITH GENE TABLE - FOREIGN KEY
  ######################################
  
  print(paste0(Sys.time(), " - merging novel junctions with other lists..."))
  
  ## Add the GENE ID for the foreign key
  # df_all_novels_tidy <- df_all_novels_tidy %>%
  #   left_join(all_genes %>% as.data.table(),
  #              by = "gene_id") %>%
  #   dplyr::select(-gene_id) %>%
  #   dplyr::rename(gene_id = id)
  
  df_all_novels_tidy <- merge(x = df_all_novels_tidy %>% as.data.table(),
                              y = all_genes %>% as.data.table(),
                              by.x = "gene_id",
                              by.y = "gene_id",
                              all.x = T) %>%
    dplyr::select(-gene_id) %>%
    dplyr::rename(gene_id = id)
  
  
  ####################################
  ## DISCARD AMBIGUOUS JUNCTIONS
  ####################################
  
  print("Removing ambiguous novel junctions across tissues ...")
  
  # Same novels should have assigned the same intron across tissues
  
  df_ambiguous <- df_all_novels_tidy %>%
    dplyr::group_by(novel_junID) %>%
    distinct(ref_junID) %>%
    dplyr::count() %>% 
    dplyr::filter(n > 1)
  
  ## If there are ambigous junctions, we remove them
  if (df_ambiguous %>% nrow() > 1) {
    
    if (age) {
      file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/ambiguous_novel_junctions_age.rds"
    } else {
      file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/ambiguous_novel_junctions.rds"
    }
    
    
    saveRDS(object = df_all_novels_tidy %>% 
              dplyr::filter(novel_junID %in% df_ambiguous$novel_junID), 
            file = file_name)
    
    df_all_novels_tidy <- df_all_novels_tidy %>% 
      dplyr::filter(!(novel_junID %in% df_ambiguous$novel_junID))
    
    # any((df_all_novels_tidy %>% 
    #        dplyr::filter(novel_junID %in% df_ambiguous$novel_junID))$ref_junID %in% db_missing_intron$ref_junID)
  }
  

  ######################################
  ## ADD MAXENTSCAN INFO 
  ######################################
  
  print(paste0(Sys.time(), " - adding MaxEntScan scores ..."))
  
  df_all_novels_tidy <- df_all_novels_tidy %>%
    distinct(novel_junID, .keep_all = T)
  
  wd <- getwd()
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_all_novels_tidy %>% dplyr::rename(junID = novel_junID),
                                                 max_ent_tool_path = "/home/sruiz/fordownload/",
                                                 homo_sapiens_fasta_path = "/data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
  
  all_split_reads_tidy %>% as_tibble()
  
  all_split_reads_tidy <- all_split_reads_tidy %>% 
    dplyr::select(-donorSeqStart,
                  -donorSeqStop,
                  -AcceptorSeqStart,
                  -AcceptorSeqStop,
                  -donor_sequence,
                  -acceptor_sequence)
  
  setwd(wd)
  
  
  df_all_novels_tidy <- all_split_reads_tidy %>%
    mutate(novel_junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) 
  
  if (setdiff(all_split_reads_tidy$junID, df_all_novels_tidy$novel_junID) %>% length() > 0) {
    print("ERROR!")
  }
  
  df_all_novels_tidy <- df_all_novels_tidy %>%
    dplyr::select(-one_of("junID","novel_ss5score","novel_ss3score")) %>% 
    dplyr::rename(novel_ss5score = ss5score, 
                  novel_ss3score = ss3score) %>%
    dplyr::relocate(ref_junID, novel_junID) %>%
    dplyr::relocate(c(novel_ss5score, novel_ss3score), .before = novel_type ) 
  
  
  df_all_novels_tidy %>% as_tibble()
  
  ####################################
  ## CREATE NOVEL JUNCTION TABLE
  ####################################

  query <- paste0("CREATE TABLE IF NOT EXISTS 'novel'",
                  "(novel_junID NUMERIC NOT NULL,
                  ref_junID NUMERIC NOT NULL,
                  novel_coordinates TEXT NOT NULL, 
                  novel_ss5score DOUBLE NOT NULL, 
                  novel_ss3score DOUBLE NOT NULL,
                  novel_type TEXT NOT NULL, 
                  distance INTEGER NOT NULL,
                  PRIMARY KEY (ref_junID, novel_junID),
                  FOREIGN KEY (ref_junID) REFERENCES 'intron'(ref_junID))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  print("'Novel' table created!")
  
  
  ## POPULATE NOVEL JUNCTION TABLE  -----------------------------------------
  
  df_all_novels_tidy <- df_all_novels_tidy %>% 
    as_tibble() %>%
    dplyr::rename(novel_coordinates = novel_junID ) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("novel_junID") %>%
    dplyr::select(-gene_id, -seqnames, 
                  -start,
                  -end, -strand)
  
  DBI::dbAppendTable(conn = con,
                     name = "novel", 
                     value = df_all_novels_tidy)
  
  
  print("'Novel' table populated!")
  
  ## ADD INDEX TO THE TABLE TO SPEED UP QUERIES -------------------------------
  query <- paste0("CREATE UNIQUE INDEX 'index_novel' ON 'novel'(novel_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  
  ## DISCONNECT ---------------------------------------------------------------
  DBI::dbDisconnect(conn = con)
  
}
# create_novel_table()

###################################
## CREATE AND POPULATE TABLES 
###################################

create_cluster_tables <- function(database_path,
                                  age = F) {
    
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbListTables(conn = con)
  
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## GET FROM INTRON TABLE
  query = paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query) 
  df_intron %>% nrow()
  
  ## GET FROM NOVEL JUNCTION TABLE
  query = paste0("SELECT * FROM 'novel'")
  df_novel <- dbGetQuery(con, query) 
  df_novel %>% nrow()
  
  ## GET FROM GENE TABLE
  query = paste0("SELECT * FROM 'gene'")
  df_gene <- dbGetQuery(con, query)
  
  if (age) {
    file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/ambiguous_novel_junctions_age.rds"
  } else {
    file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/ambiguous_novel_junctions.rds"
  }
  ambiguous_novel_junc <- readRDS(file = file_name)
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  for (db in SRA_projects) {
    
    # db <- SRA_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
  
    for (cluster in clusters) { 
      
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
     
      
      ###############################
      ## CREATE INTRON TABLE
      ###############################
      
      # dbRemoveTable(conn = con, paste0(cluster, "_", db))
      query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_misspliced"), "'", 
                                       "(ref_junID INTEGER NOT NULL,
                                       novel_junID INTEGER NOT NULL,
                                       ref_n_individuals INTEGER NOT NULL,
                                       ref_sum_counts INTEGER NOT NULL,
                                       ref_type TEXT NOT NULL, 
                                       novel_n_individuals INTEGER NOT NULL, 
                                       novel_sum_counts INTEGER NOT NULL, 
                                        
                                       MSR_D DOUBLE NOT NULL, 
                                       MSR_A DOUBLE NOT NULL, 
                                       
                                       gene_tpm DOUBLE,
                                       gene_id INTEGER NOT NULL,
                                       FOREIGN KEY (ref_junID, novel_junID) REFERENCES novel (ref_junID, novel_junID),
                      FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
      
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
      query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_nevermisspliced"), "'", 
                      "(ref_junID INTEGER NOT NULL,
                                       ref_n_individuals INTEGER NOT NULL, 
                                       ref_sum_counts INTEGER NOT NULL,
                                       MSR_D DOUBLE NOT NULL, 
                                       MSR_A DOUBLE NOT NULL, 
                                       ref_type TEXT NOT NULL, 
                                       gene_tpm DOUBLE,
                                       gene_id INTEGER NOT NULL,
                                       FOREIGN KEY (ref_junID) REFERENCES intron (ref_junID),
                      FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
      
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
      
      ###############################
      ## INSERT DATA
      ###############################
      
      
      
      if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                        cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds")) && 
          file.exists( paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                              cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))) {
        
        
        ## LOAD INTRONS AND NOVEL JUNCTIONS ------------------------------------
        
        df_introns_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                            cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds")) %>% as_tibble()
        
        df_novel_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                             cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds")) %>% as_tibble()
        
        
     
        # df_introns_gr %>%
        #   dplyr::filter(ref_junID == "chrX:100629987-100630758:-")
        # df_novel_gr %>%
        #   dplyr::filter(ref_junID == "chrX:100629987-100630758:-")
        
        ## TIDY DATA ----------------------------------------------------------- 
        
        df_introns_tidy <- df_introns_gr %>%
          distinct(ref_junID, .keep_all = T) %>%
          mutate(strand = strand %>% as.character()) %>%
          rowwise() %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID))%>%
          dplyr::select(ref_junID, 
                        ref_n_individuals, 
                        ref_sum_counts,
                        MSR_D = ref_missplicing_ratio_tissue_ND,
                        MSR_A = ref_missplicing_ratio_tissue_NA,
                        #gene_tpm = tpm_mean_rct,
                        ref_type)
        
       
        
        df_novel_tidy <- df_novel_gr %>%
          distinct(novel_junID, .keep_all = T) %>%
          mutate(strand = strand %>% as.character()) %>%
          rowwise() %>%
          mutate(novel_junID = ifelse(str_detect(string = novel_junID, pattern = "\\*"), 
                                      str_replace(string = novel_junID, pattern = "\\*", strand ),
                                      novel_junID)) %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID)) %>% 
          dplyr::select(ref_junID,
                        novel_junID,
                        novel_n_individuals, 
                        novel_sum_counts)
        

        
        ####################################
        ## QC
        ####################################
        
        if (any(str_detect(df_introns_tidy$ref_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(str_detect(df_novel_tidy$ref_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(str_detect(df_novel_tidy$novel_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(setdiff(df_novel_tidy$ref_junID, df_introns_tidy$ref_junID) > 0)) {
          print("Error! Some novel junctions have different reference introns!")
        }
        if (any(df_introns_tidy %>%
                dplyr::filter(ref_type == "never") %>%
                pull(ref_junID) %in% df_novel_tidy$ref_junID)) {
          print("Error! Some 'never misspliced' juctions are linked to novel junctions")
        }
        if ((df_introns_tidy %>% dplyr::filter(ref_type != "never") %>% distinct(ref_junID) %>% nrow() == 
             df_novel_tidy$ref_junID %>% unique %>% length()) && 
            (df_novel_tidy$ref_junID %>% unique %>% length() ==
             intersect(df_introns_tidy$ref_junID, df_novel_tidy$ref_junID) %>% unique %>% length())) {
          print("OK!")
        }
        
        if (age) {
          ambiguous_introns <- readRDS(file = "database/ambiguous_introns_age.rds")
          df_introns_tidy <- df_introns_tidy %>%
            dplyr::filter(!(ref_junID %in% ambiguous_introns))
          df_novel_tidy <- df_novel_tidy %>%
            dplyr::filter(!(ref_junID %in% ambiguous_introns))
        }
        
        
        
        ## JOIN LOCAL INTRON AND NOVELS
        df_all_misspliced <- df_introns_tidy %>% 
          inner_join(df_novel_tidy %>% as.data.table(),
                     by = "ref_junID") 
        
        
        if (any(df_all_misspliced$ref_junID %>% is.na())) {
          print("ERROR!!")
        }
  
        
        ## JOIN data with MASTER NOVEL table
        df_all_misspliced <- df_all_misspliced %>% 
          inner_join(y = df_novel %>% 
                       dplyr::select(novel_junID, novel_coordinates, novel_type) %>% 
                       as.data.table(),
                     by = c("novel_junID" = "novel_coordinates")) %>%
          dplyr::rename(novel_coordinates = novel_junID) %>%
          dplyr::rename(novel_junID = novel_junID.y)

        
        ## JOIN data with MASTER INTRON table
        df_all_misspliced <- df_all_misspliced %>% 
          inner_join(y = df_intron %>% 
                       dplyr::select(ref_junID, ref_coordinates, gene_id) %>% 
                       as.data.table(),
                     by = c("ref_junID" = "ref_coordinates")) %>%
          dplyr::select(-ref_junID) %>%
          dplyr::rename(ref_junID = ref_junID.y)
        
        
        ## PREPARE DATA PRIOR POPULATING THE TABLE
        df_all_misspliced <- df_all_misspliced %>%
          relocate(ref_junID, novel_junID)
        
       
        
        #######################################
        ## CHECK INTEGRITY WITH PARENT TABLE
        #######################################
        
        master_novel <- df_novel %>%
          dplyr::filter(novel_junID %in% 
                   (df_all_misspliced %>%
                      pull(novel_junID))) %>% 
          dplyr::select(novel_coordinates) %>% 
          as.data.table()
        
        if (!(identical(df_all_misspliced$novel_coordinates, 
                        master_novel$novel_coordinates))) {
          
          if (all(intersect(setdiff(df_all_misspliced$novel_coordinates, 
                                    master_novel$novel_coordinates),
                            ambiguous_novel_junc$novel_junID) == setdiff(df_all_misspliced$novel_coordinates, 
                                                             master_novel$novel_coordinates)) == T) {
            df_all_misspliced <- df_all_misspliced %>%
              as.data.table() %>%
              dplyr::filter(!(novel_coordinates %in% ambiguous_novel_junc$novel_junID))
          }
        }
        
        if (identical(df_all_misspliced$novel_coordinates %>% sort(), 
                      master_novel$novel_coordinates %>% sort())) {
          
          df_all_misspliced <- df_all_misspliced %>%
            dplyr::select(-novel_coordinates)
          
          ###################################################################
          
          df <- df_novel %>%
            dplyr::select(novel_junID, ref_junID) %>%
            arrange(novel_junID) %>% 
            as.data.table() %>%
            inner_join(df_all_misspliced %>%
                         dplyr::select(novel_junID, ref_junID) %>%
                         arrange(novel_junID) %>% 
                         as.data.table(),
                       by = "novel_junID")
          
          
          
          diff <- df %>% 
            dplyr::filter(ref_junID.x != ref_junID.y)
          
          
          if (diff %>% nrow() > 0) {
            df_all_misspliced <- df_all_misspliced %>%
              dplyr::filter(!(ref_junID %in% diff$ref_junID.y)) 
            print(diff)
          } else {
            print("good")
          }
          
          
          #####################################
          ## CALCULATE THE NEW MSR
          #####################################
          
          # df_all_misspliced %>% dplyr::filter(ref_junID == "17") %>% as.data.frame()
          
          db_introns <- df_all_misspliced %>%
            group_by(ref_junID, novel_type) %>%
            mutate(MSR = sum(novel_sum_counts)/(sum(novel_sum_counts) + ref_sum_counts))
          
          # db_introns %>% dplyr::filter(ref_junID == "17") %>% as.data.frame()
          
          db_introns <- db_introns %>% 
            spread(key = novel_type, value = MSR, fill = 0) %>%
            group_by(ref_junID) %>% 
            dplyr::mutate(MSR_Donor = max(novel_donor, na.rm = T)) %>%
            dplyr::mutate(MSR_Acceptor = max(novel_acceptor, na.rm = T))  %>%
            #distinct(ref_junID, .keep_all = T)  %>%
            ungroup()
          
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
          tpm <- readRDS(file = paste0(base_folder, "results/base_data/", cluster, 
                                       "/", db, "_", cluster, "_tpm.rds")) %>% 
            dplyr::select(gene_id = gene, 
                          tpm_median = TPM_median,
                          tpm_mean = TPM_mean) 
          
          ## In case there are any duplicates, take the genes with the maximum tpms
          tpm <- tpm %>% 
            group_by(gene_id) %>% 
            summarize_all(max)
            
          tpm <- tpm %>%
            inner_join(y = df_gene %>% 
                        as_tibble(),
                      by = c("gene_id" = "gene_id"))
          
         
          
          
          db_introns <- db_introns %>%
            left_join(y = tpm,
                       by = c("gene_id" = "id")) %>% 
            dplyr::rename(gene_tpm = tpm_median) %>%
            dplyr::select(-tpm_mean, -gene_id.y, 
                          -gene_name, -n_transcripts, 
                          -gene_width)
          
          
          
          # db_introns %>% dplyr::filter(ref_junID == "43423") %>% as.data.frame()
          #####################################
          ## POPULATE THE TABLE
          #####################################
          
          ## To check if there are any duplicated
          db_introns %>% 
            group_by(ref_junID, novel_junID) %>% 
            mutate(dupe = n()>1) %>%
            dplyr::filter(dupe == T)
          
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster, "_", db, "_misspliced"), 
                             value = db_introns %>%
                               dplyr::select(-novel_acceptor,-novel_donor,
                                             -MSR_D,-MSR_A)%>%
                               dplyr::rename(MSR_D = MSR_Donor, MSR_A = MSR_Acceptor))
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", paste0(cluster, "_", db, "_misspliced"), "' ON '",
                          paste0(cluster, "_", db, "_misspliced"), "'(ref_junID,novel_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          
          print(paste0(Sys.time(), ". '", paste0(cluster, "_", db, "_misspliced"), "' table populated!"))
          
          ####################################
          ## NEVER MISSPLICED
          ####################################
          
          df_never <- df_introns_tidy %>% 
            left_join(df_novel_tidy %>% as.data.table(),
                      by = "ref_junID") %>%
            dplyr::filter(is.na(novel_junID)) %>%
            #relocate(ref_mean_counts, .after = ref_n_individuals) %>%
            dplyr::select(-novel_junID,
                          -novel_n_individuals,
                          -novel_sum_counts)
          
           
          df_never <- df_never %>%
            inner_join(df_intron %>% 
                         dplyr::select(ref_junID, ref_coordinates, gene_id) %>% 
                         as.data.table(),
                       by = c("ref_junID" = "ref_coordinates")) %>%
            dplyr::filter(!is.na(ref_junID)) %>%
            dplyr::select(-ref_junID) %>% 
            dplyr::rename(ref_junID = ref_junID.y) %>%
            relocate(ref_junID)
          
          df_never %>% dplyr::filter(!is.na(ref_junID))
          
          if (any(df_never$ref_type != "never")) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
          }
          if ((intersect(df_never$ref_junID, db_introns$ref_junID) %>%
              length()) > 0) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
          }
          if (any(duplicated(df_never$ref_junID))) {
            print("Error! Some never mis-spliced junctions are duplicated.")
          }
          
          ## Remove introns and novel junctions with not TPM data
          # df_never <- df_never %>%
          #   dplyr::filter(!(gene_tpm %>% is.na()))
          
          #####################################
          ## GET THE GENE TPM
          #####################################
        
          
          
          df_never <- df_never %>%
            left_join(y = tpm,
                      by = c("gene_id" = "id")) %>% 
            dplyr::rename(gene_tpm = tpm_median) %>%
            dplyr::select(-tpm_mean, -gene_id.y, 
                          -gene_name, -n_transcripts, 
                          -gene_width)
          
          
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster, "_", db, "_nevermisspliced"), 
                             value = df_never)
          
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", paste0(cluster, "_", db, "_nevermisspliced"), "' ON '",
                          paste0(cluster, "_", db, "_nevermisspliced"),"'(ref_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          print(paste0(Sys.time(), ". '", paste0(cluster, "_", db, "_nevermisspliced"), "' table populated!"))
          
        } else {
          print("Error: novel junctions are distinct!")
          break;
        }
        
      }
      # ## QC
      # df_novel_tidy_test <- df_novel_tidy[1,]
      # df_novel_tidy_test[, "ref_junID"] <- NULL
      # df_novel_tidy_test[, "novel_junID"] <- NA
      # DBI::dbAppendTable(conn = con,
      #                    name = paste0(cluster, "_db_novel"), 
      #                    value = df_novel_tidy_test)
      # 
      # DBI::dbReadTable(conn = con, name = paste0(cluster, "_", db, "_db_novel")) %>%
      #   dplyr::filter(is.na(ref_junID))
      
    }
  }
  
  #query <- paste0("SELECT * FROM 'Brain - Putamen (basal ganglia)_BRAIN' LIMIT 10")
  #dbGetQuery(con, query)

  
  
             
  
  dbListTables(con)
  DBI::dbDisconnect(conn = con)

}
# create_cluster_tables()

###################################
## QC
###################################



# con <- dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")
# query <- paste0("SELECT * FROM 'Brain - Hippocampus_BRAIN_db_intron'")
# all_introns <- dbGetQuery(con, query)
# 
# query <- paste0("SELECT * FROM 'Brain - Hippocampus_BRAIN_db_novel'")
# all_novel <- dbGetQuery(con, query)
# 
# all_novel %>% nrow()
# all_novel %>%
#   dplyr::filter(ref_junID %in% all_introns$ref_junID) %>% nrow()
# all_introns
