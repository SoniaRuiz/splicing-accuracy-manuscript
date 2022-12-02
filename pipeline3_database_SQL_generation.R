###################################
## CONNECTION TO THE DB  
###################################
# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_SQL_generation_extra.R")
library(DBI)
library(tidyverse)
library(data.table)


# database_path <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/splicing-intronextra-v2.sqlite"
# SRA_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects.rds")
# 
# 
# 
# SRA_projects %>% length() %>% print()
# 
# con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database_path)
# dbListTables(con)

source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_junction_pairing.R")
##########################################
## LOAD REFERNCES
##########################################

if (!exists("CNC_CDTS_CONS_gr")) {

  Sys.setenv("AWS_ACCESS_KEY_ID"="AKIAUKEFAIIIAL7JUG3W",
             "AWS_SECRET_ACCESS_KEY"="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",
             "AWS_DEFAULT_REGION"="eu-west-2")
  print("Loading the 'CNC_CDTS_CONS_gr' file...")
  aws.s3::s3load(object = "CNC_CDTS_CONS_gr.rda", bucket = "data-references", region = "eu-west-2")
  print("'CNC_CDTS_CONS_gr' file loaded!")
} else {
  print("'CNC_CDTS_CONS_gr' file already loaded!")
}

## GET INFO FROM MANE
if ( !exists("hg_mane_transcripts") ) {
  
  print("Loading the 'hg_MANE' file...")
  hg_MANE <- rtracklayer::import(con = "/data/references/MANE/MANE.GRCh38.v1.0.ensembl_genomic.gtf")
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
  hg38 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
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


## Remove the tables -----------------------------------------------------------

remove_tables <- function(database_path,
                          all) {
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database_path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  tables <- dbListTables(con)
  tables
  for (table in tables) {
    #if (str_detect(table %>% tolower(), pattern = "kidney")) {
    if (!all) {
      if (!(table %in% c("gene", "transcript", "intron", 
                         "master", "novel"))) {
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



## METADATA table --------------------------------------------------------------

create_metadata_table <- function(database_path,
                                  SRA_projects,
                                  gtf_version,
                                  main_project,
                                  QC = F)  {
  
  print(paste0(Sys.time(), " - creating metadata table ... "))
  
  
  if ( str_detect(string = main_project,
                  pattern = "age") ) {
    
    df_metadata <- map_df(SRA_projects, function(project) { 
      
      # project <- SRA_projects[1]
      
      print(paste0(Sys.time(), " - getting info from ", project))
      
      project_init <- age_stratification_init_data(projects_id = project,
                                                   gtf_version = gtf_version,
                                                   main_project = main_project)
      
      project_init %>%
        dplyr::select(-subclass,-weights) %>%
        dplyr::rename(cluster = age_group,
                      SRA_project = project) %>%
        return()
    })
    
  } else {
    
    df_metadata <- map_df(SRA_projects, function(project) { 
      
      print(paste0(Sys.time(), " - getting info from ", project))
      # project <- SRA_projects[1]
      # project <- SRA_projects[2]
      # project <- SRA_projects[6]
      
      if (file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                             project, "/v", gtf_version, "/",
                             main_project, "_project/raw_data/samples_metadata.rds"))) {
        
      metadata <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                                        project, "/v", gtf_version, "/",
                                        main_project, "_project/raw_data/samples_metadata.rds"))
      
      metadata <- metadata %>% 
        as_tibble() %>%
        filter(gtex.smafrze != "EXCLUDE") %>%
        distinct(external_id, .keep_all = T) 
      
      if (main_project == "introverse") {
        
        min_samples <- 1
        
      } else {
        
        min_samples <- 70
        metadata <- metadata %>% 
          filter(gtex.smrin >= 6.0) 
        
      }
      
      if (metadata %>% nrow() >= min_samples) {
        
        
        print(metadata$gtex.smtsd %>% unique())
        print(metadata %>% distinct(gtex.sampid) %>% nrow())
        
        df_project_metadata <- data.frame(age = metadata$gtex.age %>% as.character(),
                                          rin = metadata$gtex.smrin %>% as.character(),
                                          gender = metadata$gtex.sex %>% as.character(),
                                          cluster = metadata$gtex.smtsd, 
                                          smnabtcht = metadata$gtex.smnabtcht, 
                                          sample_id = metadata %>% distinct(gtex.sampid) %>% nrow(),
                                          smafrze = metadata$gtex.smafrze,
                                          avg_read_length = metadata$recount_seq_qc.avg_len,
                                          mapped_read_count = metadata$recount_qc.star.all_mapped_reads,
                                          SRA_project = metadata$recount_project.project) %>% 
          arrange(SRA_project, cluster) %>% 
          return()
      }else {
        
        return(NULL)
      }
        
      } else {
        
        return(NULL)
      }
      
    })
  }
  
  df_metadata %>% as_tibble() %>% dplyr::count(cluster)
  
  # if (QC) {
  #   sex_tissues <- c("BREAST", "CERVIX_UTERI", "FALLOPIAN_TUBE", "OVARY", "PROSTATE", "TESTIS", "UTERUS", "VAGINA")
  #   
  #   df_metadata <- df_metadata %>% 
  #     dplyr::filter(!(SRA_project %in% sex_tissues),
  #                   !(cluster %in% c("Brain - Cortex", "Brain - Cerebellum"))) %>%
  #     dplyr::filter(cluster %in% (df_metadata %>% dplyr::count(cluster) %>% dplyr::filter(n >= 70) %>% pull(cluster)))
  # }
  
  
  # saveRDS(object = df_all_projects_metadata,
  #         file = "./dependencies/df_all_projects_metadata.rds")
  
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbWriteTable(conn = con,
                    name = "master",
                    value = df_metadata %>% as_tibble(),
                    overwrite = T)
  
  DBI::dbDisconnect(conn = con)
  
  print(paste0("Table: 'master' created!"))
}



# create_mane_table <- function(database_path) {
#   
#   print(paste0(Sys.time(), " - creating master table ... "))
#   
#   con <- dbConnect(RSQLite::SQLite(), database_path)
#   
#   DBI::dbWriteTable(conn = con,
#                     name = "mane",
#                     value = hg_MANE_tidy,
#                     overwrite = T)
#   
#   # dbRemoveTable(conn = con, "mane")
#   DBI::dbDisconnect(conn = con)
#   
#   print(paste0("Table: 'mane' created!"))
# }


## Create the master tables ----------------------------------------------------

create_master_tables <- function(database_path,
                                 main_project,
                                 gtf_version) {

  
  
  
  ##########################################
  ## LOAD AND TIDY THE PAIR-WISE DISTANCES
  ##########################################
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  if ( file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
                          main_project, "/all_paired_intron_novel_tidy.rds")) ) {
    
    print(paste0(Sys.time(), " - loading the pre-generated pair-wise distances data..."))
    
    df_all_distances_pairings_raw <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
                                                           main_project, "/all_paired_intron_novel_tidy.rds")) #"/home/sruiz/PROJECTS/splicing-project-recount3/database/all_paired_intron_novel_tidy.rds")
    df_ambiguous_novel <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
                                                main_project, "/df_all_tissues_raw_distances_ambiguous.rds"))
    df_introns_never <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
                                              main_project, "/df_all_nevermisspliced_introns.rds"))
    
  } else {
    
    print("ERROR loading dependencies!")
    # df_all_novel_raw <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/df_all_distances_pairings_raw.rds")
    # 
    # df_ambiguous_novel <- df_all_novel_raw %>%
    #   group_by(novel_junID) %>%
    #   mutate(distances_sd = distance %>% sd()) %>%
    #   dplyr::filter(distances_sd > 0)
    # 
    # saveRDS(df_ambiguous_novel,
    #           file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/df_all_tissues_raw_distances_ambiguous.rds")
    # df_ambiguous_novel %>% distinct(novel_junID)
    # 
    # print(paste0(Sys.time(), " - tidying the pair-wise distances data..."))
    # 
    # ## Remove ambiguous junctions
    # df_all_distances_pairings_raw <- df_all_novel_raw %>%
    #   distinct(novel_junID, ref_junID, .keep_all = T) %>%
    #   mutate(ref_strand = ref_strand %>% as.character(),
    #          novel_strand = novel_strand %>% as.character()) %>%
    #   dplyr::filter(!(novel_junID %in% df_ambiguous_novel$novel_junID))
    # 
    # 
    # ## Remove potential * in the junID of the reference introns
    # ind <- which(str_detect(string = df_all_distances_pairings_raw$ref_junID,
    #                  pattern = "\\*"))
    # df_all_distances_pairings_raw[ind, "ref_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "ref_junID"]$ref_junID, 
    #                                                        pattern = "\\*", 
    #                                                        replacement = df_all_distances_pairings_raw[ind, "ref_strand"]$ref_strand )
    # 
    # ## Remove potential * in the junID of the novel junctions
    # ind <- which(str_detect(string = df_all_distances_pairings_raw$novel_junID,
    #                         pattern = "\\*"))
    # df_all_distances_pairings_raw[ind, "novel_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "novel_junID"]$novel_junID, 
    #                                                          pattern = "\\*", 
    #                                                          replacement = df_all_distances_pairings_raw[ind, "novel_strand"]$novel_strand )
    # 
    # 
    # if (any(str_detect(string = df_all_distances_pairings_raw$ref_junID, pattern = "\\*")) |
    #     any(str_detect(string = df_all_distances_pairings_raw$novel_junID, pattern = "\\*")) ) {
    #   print("ERROR! Some junctions still have a * in their IDs!")
    # }
    # 
    # saveRDS(object = df_all_distances_pairings_raw,
    #         file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/all_paired_intron_novel_tidy.rds")
    # 
    # 
    # aux_introns_tidy <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/df_all_introns.rds")  %>%
    #   distinct(ref_junID, .keep_all = T)
    # 
    # 
    # ## Get the never-misspliced
    # df_introns_never <- get_intron_never_misspliced(database_path = database_path) 
    # df_introns_never <- df_introns_never %>%
    #   distinct(junID, .keep_all = T) %>%
    #   dplyr::rename(ref_junID = junID)
    # 
    # df_introns_never_tidy <- df_introns_never  %>%
    #   as_tibble() %>%
    #   left_join(y = aux_introns_tidy %>% dplyr::select(-tissue),
    #             by = c("ref_junID" = "ref_junID")) %>%
    #   distinct(ref_junID, .keep_all = T) %>%
    #   mutate(misspliced = F)  %>%
    #   mutate(strand = strand %>% as.character()) %>%
    #   rowwise() %>%
    #   mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
    #                             str_replace(string = ref_junID, pattern = "\\*", strand ),
    #                             ref_junID))
    # 
    # if (any(str_detect(string = df_introns_never_tidy$ref_junID, pattern = "\\*"))) {
    #   print("ERROR! Some never mis-spliced introns still have a * in their IDs!")
    # }
    # saveRDS(object = df_introns_never_tidy,
    #         file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/df_all_nevermisspliced_introns.rds")
    # 
    # 
    # ## Remove potential * in the junID of the reference introns
    # aux_introns_tidy <- aux_introns_tidy %>%
    #   mutate(strand = strand %>% as.character())
    # ind <- which(str_detect(string = aux_introns_tidy$ref_junID,
    #                         pattern = "\\*"))
    # aux_introns_tidy[ind, "ref_junID"] <- str_replace(string = aux_introns_tidy[ind, "ref_junID"]$ref_junID, 
    #                                              pattern = "\\*", 
    #                                              replacement = aux_introns_tidy[ind, "strand"]$strand )
    # if (any(str_detect(string = aux_introns_tidy$ref_junID, pattern = "\\*"))) {
    #   print("ERROR! Some introns still have a * in their IDs!")
    # }
    # saveRDS(object = aux_introns_tidy,
    #         file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/df_all_introns_tidy.rds")
  }
  
  


  ####################################
  ## A) GET ALL INTRONS
  ####################################
  
  
  ## Get all annotated introns
  print(paste0(Sys.time(), " - getting the mis-spliced introns..."))
  
  df_all_introns <- df_all_distances_pairings_raw %>%
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
  }
  
  df_introns_never <- df_introns_never %>% 
    as_tibble() %>%
    mutate(misspliced = F) %>%
    dplyr::filter(!(ref_junID %in% df_all_introns$ref_junID)) %>%
    dplyr::filter(!(ref_junID %in% df_ambiguous_novel$ref_junID)) 

  
  if (any(df_all_introns$width %>% abs() < 25)) {
    print("ERROR! some mis-spliced introns are shorter than 25bp!")
  }
  if (any(df_introns_never$width %>% abs() < 25)) {
    print("ERROR! some never mis-spliced introns are shorter than 25bp!")
  }
  if (intersect(df_introns_never$ref_junID, df_all_introns$ref_junID) %>% length() > 0) {
    print("ERROR! some never mis-spliced introns are mis-spliced!")
  }
  
  
  ## Merge
  df_introns_introverse <- rbind(df_introns_never, df_all_introns)
  if (any(str_detect(string = df_introns_introverse$ref_junID, pattern = "\\*" ))) {
    print("ERROR! IDs with an star!")
  }
  if (any(df_introns_introverse$ref_junID %>% duplicated())) {
    print("ERROR! Duplicated IDs")
  }

  saveRDS(object = df_introns_introverse,
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
                        main_project, "/df_all_introns_introverse.rds"))
  
  
  
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
    saveRDS(object = df_introns_introverse_tidy %>% dplyr::filter(ambiguous == T),
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/all_ambiguous_introns.rds")
  } else {
    df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
      dplyr::select(-ambiguous)
  }
  
  saveRDS(object = df_introns_introverse_tidy,
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
                        main_project, "/df_all_introns_introverse_tidy.rds"))
  
  ######################################
  ## GENES - CREATE GENE TABLE
  ######################################
  
  # df_introns_introverse_tidy <- readpaste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", main_project, "/df_all_introns_introverse_tidy.rds")
  create_gene_table(database_path = database_path,
                    gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id) )
  
  
  ######################################
  ## TX_JUNCTION - CREATE TX TABLE
  ######################################
  
  # df_introns_introverse_tidy <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", main_project, "/df_all_introns_introverse_tidy.rds"))
  create_transcript_table(database_path = database_path,
                          gene_ids = df_introns_introverse_tidy %>% unnest(gene_id) %>% distinct(gene_id),
                          tx_ids = df_introns_introverse_tidy %>% unnest(tx_id_junction) %>% distinct(tx_id_junction))
  
  
  ######################################
  ## INTRONS - ADD THE GENE FOREING KEY
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
  

  
  # ######################################
  # ## INTRONS - ADD MANE INFO
  # ######################################
  # 
  # print(paste0(Sys.time(), " --> adding MANE info ..."))
  # 
  # df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
  #   as.data.table() %>%
  #   distinct(ref_junID, .keep_all = T) %>% 
  #   rowwise() %>%
  #   mutate(MANE = ifelse(any((tx_id_junction %>% unlist() %>% unname()) %in% hg_mane_transcripts), T, F)) 
  # 
  # 
  # ######################################
  # ## INTRONS - ADD TSL INFO 
  # ######################################
  # 
  # print(paste0(Sys.time(), " --> adding TSL info ..."))
  # 
  # ## Store a 10 value into the TSL field so it can be an integer field and reduce disk space
  # df_introns_introverse_tidy <- df_introns_introverse_tidy %>%
  #   unnest(tx_id_junction) %>%
  #   dplyr::left_join(y = hg38_transcripts,
  #                    by = c("tx_id_junction" = "transcript_id")) %>%
  #   mutate(TSL = ifelse(is.na(transcript_support_level) | 
  #                         transcript_support_level == "N", 
  #                       10, transcript_support_level),
  #          TSL = TSL %>% as.integer())
  # 
  # ## Remove the columns we do not need
  # df_introns_introverse_tidy <- df_introns_introverse_tidy %>% 
  #   group_by(ref_junID) %>%
  #   mutate(TSL = TSL %>% min) %>%
  #   dplyr::select(-tx_id_junction,
  #                 -transcript_support_level) %>%
  #   distinct(ref_junID, .keep_all = T)
  # 
  # if (any(str_detect(df_introns_introverse_tidy$ref_junID, pattern = "\\*"))) {
  #   print("ERROR!")
  # }
  
  ######################################
  ## INTRONS - ADD MAXENTSCAN INFO 
  ######################################
  
  print(paste0(Sys.time(), " --> adding the MaxEntScan info ..."))
  
  wd <- getwd()
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_introns_introverse_tidy %>% dplyr::rename(junID = ref_junID),
                                                 max_ent_tool_path = "/home/sruiz/fordownload/",
                                                 homo_sapiens_fasta_path = "/data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
  
  all_split_reads_tidy <- all_split_reads_tidy %>% as_tibble()
  
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
    print("ERROR! ")
    break;
  }
  
  df_all_introns %>% as_tibble()
  
  df_all_introns <- df_all_introns %>%
    dplyr::select(-one_of("junID","ref_ss5score","ref_ss3score")) %>% 
    dplyr::rename(ref_ss5score = ss5score, 
                  ref_ss3score = ss3score) %>%
    dplyr::relocate(ref_junID) %>%
    dplyr::relocate(c(ref_ss5score, ref_ss3score), .before = transcript_id ) 
  
  
  df_all_introns %>% as_tibble()
  
  
  ################################################
  ## INTRONS - ADD THE CONSERVATION AND CDTS INFO
  ################################################
  
  df_all_introns %>% as_tibble()
  
  print(paste0(Sys.time(), " - adding CDTS and Conservation scores..."))
  
  df_all_introns <- add_cdts_cons_scores(db_introns = df_all_introns %>% as_tibble()) %>%
    as_tibble()
  
  df_all_introns <- df_all_introns %>%
    dplyr::rename(ref_CDTS5score = CDTS_5ss_mean,
                  ref_CDTS3score = CDTS_3ss_mean,
                  ref_cons5score = phastCons20way_5ss_mean,
                  ref_cons3score = phastCons20way_3ss_mean)
  
  df_all_introns %>% as_tibble()
  
  ######################################
  ## INTRONS - ADD THE BIOTYPE INFO
  ######################################
  
  print(paste0(Sys.time(), " - adding the transcript biotype ... "))
  
  df_protein <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_annotated_SR_details_length_105_biotype.rds") %>%
    as_tibble()
  
  print(paste0(Sys.time(), " - ", intersect(df_all_introns$ref_junID, df_protein$junID) %>% length(), " total number of introns found!."))
  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns."))
  
  df_all_introns <- df_all_introns %>%
    left_join(y = df_protein,
              by = c("ref_junID" = "junID")) %>%
    dplyr::rename(protein_coding = percent_PC,
                  lncRNA = percent_lncRNA)
  
  print(paste0(Sys.time(), " - ", df_all_introns %>% distinct(ref_junID) %>% nrow(), " total number of introns after adding transcript biotype."))
  
  
  
  ######################################
  ## INTRONS - ADD THE CLINVAR DATA
  ######################################
  
  print(paste0(Sys.time(), " - adding the ClinVar data...")) 
  
  clinvar_tidy <- readRDS(file = "/data/references/clinvar/clinvar_intronic_tidy.rda") %>%
    GRanges() %>%
    diffloop::rmchr()
  clinvar_tidy %>% head()
  
  df_all_introns_tidy <- df_all_introns%>%
    mutate(clinvar_type = "-") %>%
    mutate(clinvar_start = 0) %>%
    mutate(clinvar_end = 0) %>%
    GRanges() 
  
  ## Find overlaps between clinvar mutations and mis-splicing ratios
  overlaps <- GenomicRanges::findOverlaps(query = clinvar_tidy,
                                          subject = df_all_introns_tidy,
                                          ignore.strand = F,
                                          type = "any")
  
  
  elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_type"] <- ifelse(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"] == "none",
                                                                               "-", paste0(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"], "-",
                                                                                           elementMetadata(clinvar_tidy)[queryHits(overlaps), "distance"], "bp"))
  elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_start"] <- clinvar_tidy[queryHits(overlaps), ] %>% start()
  elementMetadata(df_all_introns_tidy)[subjectHits(overlaps), "clinvar_end"] <- clinvar_tidy[queryHits(overlaps), ] %>% end()
  
  
  df_all_introns_tidy
  
  
  ######################################
  ## INTRONS - ADD THE INTRON TYPE
  ######################################
  
  print(paste0(Sys.time(), " - adding the IAOD intron data...")) 
  
  ## Load intron type files
  u12_introns <- readRDS(file = "/data/references/IAOD/minor_introns_tidy.rds")
  u2_introns <- readRDS(file = "/data/references/IAOD/major_introns_tidy.rds")
  
  
  ## Add two new columns to incorporate info about intron type
  df_all_introns_tidy <- df_all_introns_tidy %>%
    as_tibble() %>%
    mutate(u12_intron = F) %>%
    mutate(u2_intron = F)
  

  ## MINOR INTRON
  print(paste0(Sys.time(), " - Getting junctions spliced out by the minor spliceosome."))
  overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u12_introns$seqnames,
                                                                         ranges = IRanges(start = u12_introns$start, 
                                                                                          end = u12_introns$end),
                                                                         strand = u12_introns$strand),
                                          subject = GenomicRanges::GRanges(seqnames = df_all_introns_tidy$seqnames,
                                                                           ranges = IRanges(start = df_all_introns_tidy$start, 
                                                                                            end = df_all_introns_tidy$end),
                                                                           strand = df_all_introns_tidy$strand),
                                          ignore.strand = FALSE,
                                          type = "equal")
  
  queryHits(overlaps) %>% length()
  df_all_introns_tidy[subjectHits(overlaps),]$u12_intron <- T
  
  
  

  ## MAJOR INTRON
  print(paste0(Sys.time(), " - Getting junctions spliced out by the major spliceosome."))
  
  overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u2_introns$seqnames %>% as.character(),
                                                                         ranges = IRanges(start = u2_introns$start,
                                                                                          end = u2_introns$end),
                                                                         strand = u2_introns$strand),
                                          subject = GenomicRanges::GRanges(seqnames = df_all_introns_tidy$seqnames,
                                                                           ranges = IRanges(start = df_all_introns_tidy$start,
                                                                                            end = df_all_introns_tidy$end),
                                                                           strand = df_all_introns_tidy$strand),
                                          ignore.strand = FALSE,
                                          type = "equal")
  
  # print(overlaps)
  subjectHits(overlaps) %>% length()
  df_all_introns_tidy[subjectHits(overlaps),]$u2_intron <- T
  
  ######################################
  ## INTRONS - POPULATE THE TABLE
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
  lncRNA INTEGER NOT NULL,
  protein_coding INTEGER NOT NULL,
  misspliced BOOL NOT NULL,
  transcript_id INTEGER NOT NULL,
  FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database_path)
  dbListTables(con)
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  print("'Intron' table created!")
  
  
  ## POPULATE INTRON TABLE ----------------------------------------------------
  df_all_introns_tidy_final <- df_all_introns_tidy %>% 
    #as_tibble() %>%
    #unnest(gene_id) %>%
    dplyr::rename(ref_length = width) %>%
    mutate(ref_cons5score = ifelse(ref_cons5score %>% is.na(), 0, ref_cons5score),
           ref_cons3score = ifelse(ref_cons3score %>% is.na(), 0, ref_cons3score)) %>%
    dplyr::rename(ref_coordinates = ref_junID,
                  clinvar = clinvar_type) %>%
    distinct(ref_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("ref_junID") %>%
    dplyr::select(-c(seqnames, start, end, strand, 
                     clinvar_start, clinvar_end) )
  
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
  
 
  
  print("'Intron' table populated!")
  
  ## CREATE INDEX TO SPEED UP QUERIES ------------------------------------------
  query <- paste0("CREATE UNIQUE INDEX 'index_intron' ON 'intron'(ref_junID)");
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)

  #setdiff(df_never_misspliced$ref_coordinates, df_all_introns_tidy_final$ref_coordinates) %>% length()
  
  
  ########################################
  ## B) NOVEL JUNCTION - ADD MAXENTSCAN INFO 
  ########################################
  
  print(paste0(Sys.time(), " - adding MaxEntScan scores to the NOVEL JUNCTONS..."))
  
  df_all_novel_raw_tidy <- df_all_distances_pairings_raw %>%
    mutate( start = novel_start %>% as.integer(),
            end = novel_end %>% as.integer()) %>%
    dplyr::select(seqnames = novel_seq,
                  start, end, strand = novel_strand,
                  novel_junID, ref_junID, novel_type = type, distance) %>%
    distinct(novel_junID, .keep_all = T)
  
  wd <- getwd()
  
  ## Add MaxEntScan score to the split reads
  all_split_reads_tidy <- generate_max_ent_score(junc_tidy = df_all_novel_raw_tidy %>% dplyr::rename(junID = novel_junID),
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
  
  
  
  df_all_novels_tidy <- df_all_novels_tidy %>%
    dplyr::select(-one_of("junID","novel_ss5score","novel_ss3score")) %>% 
    dplyr::rename(novel_ss5score = ss5score, 
                  novel_ss3score = ss3score) %>%
    dplyr::relocate(ref_junID, novel_junID) %>%
    dplyr::relocate(c(novel_ss5score, novel_ss3score), .before = novel_type ) 
  
  
  df_all_novels_tidy <- df_all_novels_tidy %>% as_tibble()
  df_all_novels_tidy
  
  ####################################
  ## NOVEL - ADD INTRON REFERENCE
  ####################################

  print(paste0(Sys.time(), " - adding the INTRON foreign key to the NOVEL JUNCTONS..."))
  
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
  
  
  df_all_novels_tidy_final <- df_all_novels_tidy_final %>% 
    as_tibble() %>%
    dplyr::rename(novel_coordinates = novel_junID ) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("novel_junID") %>%
    dplyr::select(-seqnames, 
                  -start,
                  -end, -strand)
  
  df_all_novels_tidy_final
  if (any(duplicated(df_all_novels_tidy_final$novel_coordinates))) {
    print("ERROR! some novel junctions are duplicated")
  }
  DBI::dbAppendTable(conn = con,
                     name = "novel", 
                     value = df_all_novels_tidy_final)
  
  
  print("'Novel' table populated!")
  
  
}


create_gene_table <- function(database_path,
                              gene_ids) {
  
  
  print(paste0(Sys.time(), " --> creating the gene table..."))
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  
  hg38_transcripts_gene <- hg38 %>%
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
  
  hg38_tidy <- hg38_transcripts_gene %>%
    inner_join(hg38_genes, by = "gene_id" ) %>% 
    as_tibble()
  
  
  ## CREATE GENE_NAME TABLE
  sql_statement <- paste0("CREATE TABLE IF NOT EXISTS 'gene'", 
                          "(id INTEGER PRIMARY KEY NOT NULL,
                          gene_id TEXT NOT NULL,
                          gene_name TEXT,
                          n_transcripts INTEGER NOT NULL,
                          gene_width INTEGER NOT NULL)")
  res <- DBI::dbSendQuery(conn = con, statement = sql_statement)
  DBI::dbClearResult(res)
  
  
  ## POPULATE GENE TABLE
  hg38_tidy_final <- hg38_tidy %>% 
    as_tibble() %>%
    tibble::rowid_to_column("id")
  
  DBI::dbAppendTable(conn = con,
                     name = "gene", 
                     value = hg38_tidy_final)
  
  print(paste0(Sys.time(), " --> gene table created!"))
}


create_transcript_table  <- function(database_path,
                                     gene_ids,
                                     tx_ids) {
  
  print(paste0(Sys.time(), " --> creating the 'transcript' table..."))
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
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
}


## Create the cluster tables ----------------------------------------------------

create_cluster_tables <- function(database_path,
                                  main_project,
                                  SRA_projects = NULL,
                                  gtf_version) {
  
  all_split_reads_details_105 <- readRDS(file = paste0("~/PROJECTS/splicing-project-recount3/database/v", 
                                                       gtf_version, "/",
                                                       main_project, "/all_split_reads_105_length_all_tissues.rds"))
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbListTables(conn = con)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  print( paste0(Sys.time(), " --> SQL connection stablished!") )
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## GET FROM INTRON TABLE
  query = paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query) %>% as_tibble()
  df_intron %>% nrow()
  df_intron %>% 
    dplyr::count(misspliced)
  
  ## GET FROM NOVEL JUNCTION TABLE
  query = paste0("SELECT * FROM 'novel'")
  df_novel <- dbGetQuery(con, query) %>% as_tibble()
  df_novel %>% nrow()
  df_novel %>% 
    dplyr::count(novel_type) %>%
    print()
  
  
  ## GET FROM GENE TABLE
  query = paste0("SELECT * FROM 'gene'")
  df_gene <- dbGetQuery(con, query) %>% as_tibble()
  df_gene %>% nrow()
  
  ## GET FROM TRANSCRIPT TABLE
  query = paste0("SELECT * FROM 'transcript'")
  df_transcript <- dbGetQuery(con, query) %>% as_tibble()
  df_transcript %>% nrow()
  
  if ( is.null(SRA_projects) ){
    SRA_projects <- (df_metadata$SRA_project %>% unique())
  }
  
  
  
  for (db in SRA_projects) {
    
    # db <- SRA_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                          db, "/v", gtf_version, "/", main_project, "_project/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    for (cluster in clusters) { 
      
      # cluster <- clusters[1]
      print(paste0(Sys.time(), " --> ", cluster))
      
      ###############################
      ## PREPARE DATA
      ###############################
      
      if ( file.exists(paste0(base_folder, "results/", 
                              cluster, "/", cluster, "_raw_distances_tidy.rds")) ) {
        
        
        ## BASE DATA -----------------------------------------------------------
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(base_folder, "/raw_data/", 
                                                   db, "_", cluster, 
                                                   "_split_read_counts_sample_tidy.rds")) 
        if ( is.null(names(split_read_counts)) ) {
          split_read_counts <- split_read_counts %>%
            as_tibble(rownames = "junID")
        }
        
        ## Load samples
        samples <- readRDS(file = paste0(base_folder, "/raw_data/", 
                                         db, "_", cluster, "_samples_used.rds"))
        
        
        
        ## LOAD INTRONS AND NOVEL JUNCTIONS ------------------------------------
        df_cluster_distances <- readRDS(file = paste0(base_folder, "results/", 
                                                      cluster, "/", cluster, 
                                                      "_raw_distances_tidy.rds")) %>% as_tibble()
        
        
        
        ## INTRONS ---------------------------------------------------
        df_introns_gr <- df_cluster_distances %>%
          distinct(ref_junID, .keep_all = T) %>%
          dplyr::select(ref_junID, seqnames = ref_seq, start = ref_start,
                        end = ref_end, strand = ref_strand)
        
        
        
        ## Add reads detected for the introns in the current tissue
        split_read_counts_intron <- get_mean_coverage(split_read_counts = split_read_counts,
                                                      samples = samples,
                                                      junIDs = df_introns_gr$ref_junID) %>%
          dplyr::rename(ref_n_individuals = n_individuals,
                        ref_sum_counts = sum_counts)
        split_read_counts_intron %>% head()
        
        df_introns_merged <- df_introns_gr %>%
          left_join(y = split_read_counts_intron,
                    by = c("ref_junID" = "junID"))
        
        
        ## NOVEL JUNCTIONS -----------------------------------------------------
        df_novel_gr <- df_cluster_distances %>%
          distinct(novel_junID, .keep_all = T) %>%
          dplyr::select(novel_junID, ref_junID, seqnames = novel_seq, 
                        start = novel_start, end = novel_end, strand = novel_strand)
        
       
        split_read_counts_novel <- get_mean_coverage(split_read_counts = split_read_counts,
                                                     samples = samples,
                                                     junID = df_novel_gr$novel_junID) %>%
          dplyr::rename(novel_n_individuals = n_individuals,
                        novel_sum_counts = sum_counts)
        
        
        df_novel_merged <- df_novel_gr %>%
          left_join(y = split_read_counts_novel,
                    by = c("novel_junID" = "junID"))
        
        df_novel_merged %>% as_tibble()
        
        
        ## QC data -------------------------------------------------------------
        df_novel_merged <- df_novel_merged %>% 
          rowwise() %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID)) 
        
        df_novel_merged <- df_novel_merged %>% 
          rowwise() %>%
          mutate(novel_junID = ifelse(str_detect(string = novel_junID, pattern = "\\*"),
                                      str_replace(string = novel_junID, pattern = "\\*", strand ),
                                      novel_junID)) 
        
        df_introns_merged <- df_introns_merged %>% 
          rowwise() %>%
          mutate(strand = strand %>% as.character()) %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"),
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID))
        
        
        ## Merge tidy data -----------------------------------------------------
        df_intron_novel_merged <- df_novel_merged %>%
          dplyr::select(-c(seqnames, start, end, strand)) %>%
          inner_join(y = df_introns_merged %>% 
                       dplyr::select(-c(seqnames, start, end, strand)),
                     by = "ref_junID")
        
        
        
        
        ####################################
        ## QC
        ####################################
        
        if (any(str_detect(df_intron_novel_merged$ref_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(str_detect(df_intron_novel_merged$novel_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(df_intron_novel_merged$ref_junID %>% is.na())) {
          print("ERROR!!")
        }
        
     
        
        ## JOIN data with MASTER NOVEL table
        df_all_misspliced <- df_intron_novel_merged %>% 
          inner_join(y = df_novel %>% 
                       dplyr::select(novel_junID, novel_coordinates, novel_type) %>% 
                       as.data.table(),
                     by = c("novel_junID" = "novel_coordinates")) %>%
          dplyr::rename(novel_coordinates = novel_junID) %>%
          dplyr::rename(novel_junID = novel_junID.y)
        
       
        if(setdiff(df_all_misspliced$ref_junID, df_intron$ref_coordinates) %>% length() > 0) {
          print("ERROR! introns detected in this tissue are not stored in the master intron table.")
          break;
        }
        
        ## JOIN data with MASTER INTRON table
        df_all_misspliced <- df_all_misspliced %>% 
          left_join(y = df_intron %>%
                      dplyr::filter(misspliced == T) %>%
                      dplyr::select(ref_junID, ref_coordinates, transcript_id) %>% 
                      as.data.table(),
                    by = c("ref_junID" = "ref_coordinates")) %>%
          dplyr::select(-ref_junID) %>%
          dplyr::rename(ref_junID = ref_junID.y)
        
        
        if (which(str_detect(df_all_misspliced$novel_coordinates,pattern = "//*")) %>% length() > 0) {
          print("ERROR: some IDs contain '*'")
        }
        
        
        ## PREPARE DATA PRIOR POPULATING THE TABLE
        df_all_misspliced <- df_all_misspliced %>%
          relocate(ref_junID, novel_junID)
        
        
        
        #######################################
        ## CHECK INTEGRITY WITH PARENT TABLES
        #######################################
        
        master_novel <- df_novel %>%
          dplyr::filter(novel_junID %in% 
                          (df_all_misspliced %>%
                             pull(novel_junID))) %>% 
          dplyr::select(novel_coordinates) %>% 
          as.data.table()
        
        if ( !(identical(df_all_misspliced$novel_coordinates %>% sort(), 
                         master_novel$novel_coordinates %>% sort())) ) {
          
          if ( all(intersect(setdiff(df_all_misspliced$novel_coordinates, 
                                     master_novel$novel_coordinates),
                             ambiguous_novel_junc$novel_junID) == setdiff(df_all_misspliced$novel_coordinates, 
                                                                          master_novel$novel_coordinates)) == T ) {
            df_all_misspliced <- df_all_misspliced %>%
              as.data.table() %>%
              dplyr::filter(!(novel_coordinates %in% ambiguous_novel_junc$novel_junID))
          }
        }
        
        if ( identical(df_all_misspliced$novel_coordinates %>% sort(), 
                       master_novel$novel_coordinates %>% sort()) ) {
          
          df_all_misspliced <- df_all_misspliced %>%
            dplyr::select(-novel_coordinates)
          
          
          ## CHECK PARENT INTEGRITY
          
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
            print(paste0("ERROR!: ", diff, " --> mismatch junctions."))
            break;
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
            ungroup()
          
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
          
          
          if ( str_detect(string = main_project, pattern = "age") ) {
            
            body_sites <- df_metadata %>%
              filter(SRA_project == db) %>%
              distinct(region) %>%
              pull()
            
            tpm <- map_df(body_sites, function(site) {
              # site <- body_sites[1]
              tpm <- readRDS(file = 
                               paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", db, 
                                      "/v", gtf_version, "/splicing_project/results/tpm/",
                                      db, "_", site, "_tpm.rds")) %>% 
                dplyr::select(gene_id = gene, any_of(samples))
              
            })
            
            if ( !identical(names(tpm)[-1] , samples) ){
              print("ERROR! Samples obtained for the TPM analysis are not identical")
              break;
            }
            tpm[is.na(tpm)] <- 0
            
          } else {
            
            tpm <- readRDS(file = paste0(base_folder, "results/tpm/",
                                         "/", db, "_", cluster, "_tpm.rds")) %>% 
              dplyr::select(gene_id = gene, all_of(samples))
          }
          
          tpm <- tpm  %>%
            dplyr::mutate(tpm_median = matrixStats::rowMedians(x = as.matrix(.[2:(ncol(tpm))]))) %>%
            dplyr::select(gene_id, tpm_median) 
          
          
          ## In case there are any duplicates, take the genes with the maximum tpms
          tpm <- tpm %>% 
            group_by(gene_id) %>% 
            summarize_all(max)
          
          tpm_tidy <- tpm %>%
            inner_join(y = df_gene %>% 
                         as_tibble(),
                       by = c("gene_id" = "gene_id")) %>%
            inner_join(y = df_transcript %>% 
                         as_tibble(),
                       by = c("id" = "gene_id")) %>%
            dplyr::select(transcript_id = id.y, 
                          tpm_median)
          
          db_introns <- db_introns %>%
            left_join(y = tpm_tidy,
                      by = c("transcript_id" = "transcript_id")) %>% 
            dplyr::rename(gene_tpm = tpm_median)
          
          
          
          #####################################
          ## ADD THE REFERENCE TYPE
          #####################################
          
          
          db_introns[, "ref_type"] <- ""
          db_introns %>% head()
          
          ## TYPE 'BOTH'
          indx <- which(db_introns$MSR_Donor > 0 & db_introns$MSR_Acceptor > 0)
          db_introns[indx, "ref_type"] <- "both"
          print(paste0("Introns type 'both': ", indx %>% length()))
          
          
          
          ## TYPE 'DONOR'
          indx <- which(db_introns$MSR_Donor > 0 & db_introns$MSR_Acceptor == 0)
          db_introns[indx, "ref_type"] <- "donor"
          print(paste0("Junctions type 'donor': ", indx %>% length()))
          
          
          ## TYPE 'ACCEPTOR'
          indx <- which(db_introns$MSR_Donor == 0 & db_introns$MSR_Acceptor > 0)
          db_introns[indx, "ref_type"] <- "acceptor"
          print(paste0("Junctions type 'acceptor': ", indx %>% length()))
          
         
      

          
          db_introns %>% nrow()
          db_introns %>% distinct(ref_junID) %>% nrow()
          
          
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
                          transcript_id INTEGER NOT NULL, 
                          FOREIGN KEY (ref_junID, novel_junID) REFERENCES novel (ref_junID, novel_junID),
                          FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
          
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          
          #####################################
          ## POPULATE THE TABLE
          #####################################
          
   
          db_introns_final <- db_introns %>%
            dplyr::select(-novel_acceptor,-novel_donor)%>%
            dplyr::rename(MSR_D = MSR_Donor, MSR_A = MSR_Acceptor)
          
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster, "_", db, "_misspliced"), 
                             value = db_introns_final)
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", paste0(cluster, "_", db, "_misspliced"), "' ON '",
                          paste0(cluster, "_", db, "_misspliced"), "'(ref_junID,novel_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          
          print(paste0(Sys.time(), ". '", paste0(cluster, "_", db, "_misspliced"), 
                       "' table populated!"))
          
          ####################################
          ## NEVER MISSPLICED
          ####################################
          
          ## TYPE 'NONE'
          introns_never <- readRDS(file = paste0(base_folder, "results/", 
                                                cluster, "/not-misspliced/", 
                                                cluster, "_all_notmisspliced.rds"))
          ## The introns not misspliced in this tissue, should have not been detected as spliced.
          ## Thus, this should be zero
          if (intersect(introns_never, df_intron %>%
                        dplyr::filter(ref_junID %in% db_introns_final$ref_junID) %>%
                        pull(ref_coordinates) %>% 
                        unique()) %>% length() > 0 ){
            print("ERROR! The introns not misspliced in this tissue, should have not been detected as spliced.")
            break;
          }
          df_introns_never <- data.frame(ref_junID = introns_never)
          
          
        
          if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
            print("ERROR! * in the IDs")
          }
          
          split_read_counts <- split_read_counts %>%
            left_join(y = all_split_reads_details_105 %>% 
                        dplyr::select(junID, seqnames, start, end, strand),
                      by = "junID") %>%
            dplyr::select(-junID) %>%
            mutate(junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>%
            dplyr::relocate(junID)
          
          split_read_counts_intron_never <- get_mean_coverage(split_read_counts = split_read_counts,
                                                              samples = samples,
                                                              junID = df_introns_never$ref_junID) %>%
            dplyr::rename(ref_n_individuals = n_individuals,
                          ref_sum_counts = sum_counts) 
          
          if (any(str_detect(split_read_counts_intron_never$junID,pattern = "\\*"))) {
            print("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          
          
          df_never_merged <- df_introns_never %>%
            inner_join(y = split_read_counts_intron_never,
                      by = c("ref_junID" = "junID")) %>%
            mutate(MSR_D = 0, MSR_A = 0) %>%
            mutate(ref_type = "never")
          df_never_merged <- df_never_merged %>% as_tibble()
          
          
          
          if (any(df_never_merged$ref_n_individuals %>% is.na)) {
            print("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          
          
          ## QC
          if (any(str_detect(df_never_merged$ref_junID, pattern = "\\*"))) {
            print("ERROR! * in the IDs")
            break;
          }
          
          
          ## TYPE 'NONE'
          print(paste0("Junctions type 'never': ", df_never_merged %>% distinct(ref_junID) %>% nrow()))
          
          ##################################################
          ## ADD THE INTRON REFERENCE KEY
          ##################################################
          
          df_never_merged <- df_never_merged %>%
            inner_join(df_intron %>% 
                         dplyr::select(ref_junID, ref_coordinates, transcript_id) %>% 
                         as.data.table(),
                       by = c("ref_junID" = "ref_coordinates")) %>%
            dplyr::filter(!is.na(ref_junID)) %>%
            dplyr::select(-ref_junID) %>% 
            dplyr::rename(ref_junID = ref_junID.y) %>%
            relocate(ref_junID)
          
          
          if (df_never_merged %>% dplyr::filter(is.na(ref_junID)) %>% nrow > 0) {
            print("ERROR! IDs are NA")
          }
          if (any(df_never_merged$ref_type != "never")) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
          }
          if ((intersect(df_never_merged$ref_junID, db_introns$ref_junID) %>% length()) > 0) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
          }
          if (any(duplicated(df_never_merged$ref_junID))) {
            print("Error! Some never mis-spliced junctions are duplicated.")
          }
          
          ## Remove introns and novel junctions with not TPM data
          # df_never <- df_never %>%
          #   dplyr::filter(!(gene_tpm %>% is.na()))    
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
      
          df_never_merged <- df_never_merged %>%
            left_join(y = tpm_tidy,
                      by = "transcript_id") %>% 
            dplyr::rename(gene_tpm = tpm_median) 
          
          
          
          ## QC --------------------------------------------------------
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "donor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "aceptor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "both") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          if (intersect(db_introns %>%
                        dplyr::filter(ref_type == "acceptor") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "donor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
          }
          
          
          ###############################
          ## CREATE NEVER MIS-SPLICED TABLE
          ###############################
          
          query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_nevermisspliced"), "'", 
                          "(ref_junID INTEGER NOT NULL,
                          ref_n_individuals INTEGER NOT NULL, 
                          ref_sum_counts INTEGER NOT NULL,
                          MSR_D DOUBLE NOT NULL, 
                          MSR_A DOUBLE NOT NULL, 
                          ref_type TEXT NOT NULL, 
                          gene_tpm DOUBLE,
                          transcript_id INTEGER NOT NULL,
                          FOREIGN KEY (ref_junID) REFERENCES intron (ref_junID),
                          FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
          
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          
          ###############################
          ## POPULATE TABLE
          ###############################
          
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster, "_", db, "_nevermisspliced"), 
                             value = df_never_merged)
          
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", 
                          paste0(cluster, "_", db, "_nevermisspliced"), "' ON '",
                          paste0(cluster, "_", db, "_nevermisspliced"),"'(ref_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          print(paste0(Sys.time(), ". '", 
                       paste0(cluster, "_", db, "_nevermisspliced"), 
                       "' table populated!"))
          
        } else {
          print("Error: novel junctions are distinct!")
          break;
        }
        
      }
   
    }
  }
  
  
  dbListTables(con)
  DBI::dbDisconnect(conn = con)
  
}



###################################################
## CALLS
###################################################

#remove_tables(database_path, all = T)
#create_metadata_table(database_path)
#create_master_tables(database_path)
#create_cluster_tables(database_path)
