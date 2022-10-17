library(tidyverse)
library(data.table)
library(GenomicRanges)
library(DBI)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/sql_helper.R")
# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline1_download_from_recount.R")
# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline2_annotate_from_recount.R")
# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_generation.R")

################################################################
## UTILS
#################################################################


get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}

get_mean_coverage <- function(split_read_counts,
                              samples,
                              junIDs) {
  
  split_read_counts_intron <- split_read_counts %>%
    dplyr::filter(junID %in% junIDs) %>%
    dplyr::select(junID, all_of(samples %>% as.character())) 
  
  split_read_counts_intron[,"n_individuals"] <- (matrixStats::rowCounts(split_read_counts_intron[, -c(1)] > 0, na.rm = T)) 
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron[,"sum_counts"] <- Matrix::rowSums(split_read_counts_intron[,-c(split_read_counts_intron %>% ncol(),1)], na.rm = T)
  split_read_counts_intron <- split_read_counts_intron %>% as.data.frame()
  
  split_read_counts_intron <- split_read_counts_intron[, c(1,(split_read_counts_intron %>% ncol() - 1),(split_read_counts_intron %>% ncol()))]
  
  
  if (any(split_read_counts_intron[, "n_individuals"] < 1)) {
    print("Error: some ref junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_intron %>% return()
}

# all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
get_all_annotated_split_reads <- function(all_projects,
                                          gtf_version,
                                          all_clusters = NULL,
                                          main_project = "splicing") {
  
  
 
  
  # all_split_reads_details_105 <- readRDS(file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_105_length_all_tissues.rds")
  
  #############################################
  ## GET ALL SPLIT READS FROM 'USE ME' SAMPLES
  #############################################
  
  
  
  all_split_reads_details_105 <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    # project_id <- "BONE_MARROW"
    
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          project_id, "/v", gtf_version, "/", 
                          main_project, "_project/")
    
    if (is.null(all_clusters)) {
      
      metadata.info <- readRDS(file = paste0(folder_root, "/raw_data/samples_metadata.rds"))
      all_clusters <-  metadata.info %>% 
        as_tibble() %>%
        filter(gtex.smrin >= 6.0,
               gtex.smafrze != "EXCLUDE") %>%
        distinct(gtex.smtsd) %>% 
        pull()
    }
    
    
    jxn_qc <- map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[1]
      print(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
      
      if (file.exists(paste0(folder_root, "/raw_data/", project_id, "_",
                             cluster, "_all_split_reads_sample_tidy.rds"))) {
        
        all_split_reads_details_105 <- readRDS(file = paste0(folder_root, "/raw_data/", project_id, "_",
                                                             cluster, "_all_split_reads_sample_tidy.rds"))
        
        ## Remove split reads annotated to multiple genes
        all_split_reads_details_tidy <- all_split_reads_details_105 %>%
          distinct(junID, .keep_all = T) %>% 
          rowwise() %>%
          mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
        
        all_split_reads_details_tidy <- all_split_reads_details_tidy %>%
          filter(ambiguous == F)
        
        saveRDS(all_split_reads_details_tidy %>% dplyr::select(-ambiguous),
                file =  paste0(folder_root, "/raw_data/", project_id, "_",
                               cluster, "_annotated_SR_details_length_105.rds"))
        
        ## Print message
        print(paste0(all_split_reads_details_105 %>% nrow, " - ", all_split_reads_details_tidy %>% nrow()))
        
        all_split_reads_details_tidy %>%
          distinct(junID, .keep_all = T) %>% 
          as_tibble() %>%
          return()
        
      } else {
        return(NULL)
      }
      
    })
    if (jxn_qc %>% nrow() > 0 ) {
      jxn_qc %>%
        distinct(junID, .keep_all = T) %>% return()
    } else {
      return(NULL)
    }
    
  })
  
  
  print(paste0(Sys.time(), " - saving 'all_annotated_split_reads' for the database!"))
  
  all_split_reads_details_105 <- all_split_reads_details_105 %>%
    distinct(junID, .keep_all = T)
  
  database_folder <- paste0("~/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", main_project, "/")
  dir.create(file.path(database_folder), showWarnings = T, recursive = T)
  
  saveRDS(object = all_split_reads_details_105,
          file = paste0(database_folder, "/all_split_reads_105_length_all_tissues.rds"))
  
  
  
  
  
}


get_all_raw_distances_pairings <- function(all_projects,
                                           gtf_version,
                                           all_clusters = NULL,
                                           main_project) {
  
  
  ## LOOP THROUGH PROJECTS
  df_all_distances_pairings_raw <- map_df(all_projects, function(db) {
    
    # db <- all_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          db, "/v", gtf_version, "/", main_project, "_project/")
    
    if (is.null(all_clusters)) {
      
      metadata.info <- readRDS(file = paste0(folder_root, "/raw_data/samples_metadata.rds"))
      all_clusters <-  metadata.info %>% 
        as_tibble() %>%
        filter(gtex.smrin >= 6.0,
               gtex.smafrze != "EXCLUDE") %>%
        distinct(gtex.smtsd) %>% 
        pull()
      
    }
    
    map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[2]
      
      print(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
      
      ## Load samples
      samples <- readRDS(file = paste0(folder_root, "/raw_data/", db, "_", cluster,  "_samples_used.rds"))
      
      if (samples %>% length() > 0) {
        
        folder_name <- paste0(folder_root, "/results/", cluster, "/")
        
        ## Obtain the distances across all samples
        df_all <- map_df(samples, function(sample) { 
          
          # sample <- samples[1]
          print(paste0(cluster, " - ", sample))
          file_name <- paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")
          
          
          if (file.exists(file_name)) {
            
            df <- readRDS(file = file_name)
            
            return(df)
          } 
          
        })
        
        if (!file.exists(paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))) {
          
          saveRDS(object = df_all %>%
                    distinct(novel_junID, ref_junID, .keep_all = T) %>%
                    mutate(tissue = cluster),
                  file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
        }
        
        
        # df_all2 <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
        return(df_all %>%
                 distinct(novel_junID, ref_junID, .keep_all = T))
      } else {
        return(NULL)
      }
    })  
  })  
  
  print(paste0(Sys.time(), " - saving 'df_all_distances_pairings_raw' for the database!"))
  
  
  saveRDS(object = df_all_distances_pairings_raw %>%
            distinct(novel_junID, ref_junID, .keep_all = T) %>% 
            as.data.table(),
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
                        main_project, "/df_all_distances_pairings_raw.rds"))
  
}

get_intron_never_misspliced <- function (projects_used,
                                         all_clusters = NULL,
                                         main_project) {
  
 
  df_never <- map_df(projects_used, function(db) {
    
    # db <- projects_used[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          db, "/v", gtf_version, "/", main_project, "_project/")
    
    if (is.null(all_clusters)) {
      
      if (file.exists(paste0(base_folder, "/raw_data/samples_metadata.rds"))) {
        metadata.info <- readRDS(file = paste0(base_folder, "/raw_data/samples_metadata.rds"))
        all_clusters <-  metadata.info %>% 
          as_tibble() %>%
          filter(gtex.smrin >= 6.0,
                 gtex.smafrze != "EXCLUDE") %>%
          distinct(gtex.smtsd) %>% 
          pull()  
      }
      
    }
    
    map_df(all_clusters, function(cluster) { 
      
      # cluster <- all_clusters[1]
      
      # print(paste0(Sys.time(), " --> ", cluster))
      if ( file.exists(paste0(base_folder, "results/", 
                              cluster, "/not-misspliced/", cluster, "_all_notmisspliced.rds")) ) {
        df_introns_never <- readRDS(file = paste0(base_folder, "results/", 
                                                  cluster, "/not-misspliced/", cluster, "_all_notmisspliced.rds")) %>% as_tibble()
        return(data.frame(ref_junID = df_introns_never$value))
      } else {
        return(NULL)
      }
      
    })
  })

  df_never %>%
    distinct(ref_junID) %>%
    return()
}



generate_max_ent_score <- function(junc_tidy,
                                   # bedtools_path = NULL,
                                   max_ent_tool_path,
                                   homo_sapiens_fasta_path){
  
  library(Biostrings)
  library(tidyverse)
  library(protr)
  
  # ## Load the annotated split reads data (both for PD and control samples)
  # junc_tidy <-readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/data/",
  #                                   project_id, "/results/base_data/", project_id, "_split_reads.rda"))
  
  junc_tidy <- junc_tidy %>% dplyr::as_tibble()
  
  junc_tidy$seqnames <- junc_tidy$seqnames %>% as.character()
  
  if (any(junc_tidy$seqnames == "M")) {
    print("Error! There's data for chr-MT!")
  } 
  
  ## 0. Prepare the object ---------------------------------------------------
  
  ## get the ranges for the donor and acceptor sequences needed for the MaxEntScan
  junc_tidy <- junc_tidy %>%  mutate(donorSeqStart = 
                                       ifelse(strand == "-",
                                              end - 6, start - 4),
                                     donorSeqStop =
                                       ifelse(strand == "-",
                                              end + 3, start + 5),
                                     AcceptorSeqStart =
                                       ifelse(strand == "-",
                                              start - 4, end - 20),
                                     AcceptorSeqStop =
                                       ifelse(strand == "-",
                                              start + 19, end + 3)) 
  
  junc_tidy[1,]
  
  to.BED <- data.frame(seqnames = junc_tidy$seqnames,
                       starts = as.integer(junc_tidy$donorSeqStart),
                       ends = as.integer(junc_tidy$donorSeqStop),
                       names = as.character(junc_tidy$junID),
                       scores = c(rep(".", nrow(junc_tidy))),
                       strands = junc_tidy$strand)
  to.BED[1,]
  
  ## 1. Obtain the genomic sequence for splice sites ------------------------------------------------
  
  ## Get the donor genomic sequence
  
  tmp.file <- tempfile()
  
  ## get the maxentscan for the 5' splice site
  write.table(to.BED, file = tmp.file, quote = F, sep = "\t", row.names = F, col.names = F)
  tmp.file_seq <- tempfile()
  system(paste0("bedtools getfasta -name -s -fi ", homo_sapiens_fasta_path, " -bed ",
                tmp.file, " -tab -fo ", tmp.file_seq))
  donor_sequences_input <- read.delim(tmp.file_seq, header = F)
  head(donor_sequences_input)
  head(junc_tidy)
  
  
  stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(donor_sequences_input$V1)))),
                      junc_tidy$junID %>% as.character()))
  junc_tidy <- cbind(junc_tidy, 
                     donor_sequence = as.character(donor_sequences_input$V2))
  
  junc_tidy %>% head()
  
  
  ## Get the acceptor genomic sequence
  
  to.BED <- data.frame(seqnames = junc_tidy$seqnames,
                       starts = as.integer(junc_tidy$AcceptorSeqStart),
                       ends = as.integer(junc_tidy$AcceptorSeqStop),
                       names = as.character(junc_tidy$junID),
                       scores = c(rep(".", nrow(junc_tidy))),
                       strands = junc_tidy$strand)
  
  
  tmp.file <- tempfile()
  
  write.table(to.BED, file = tmp.file, quote = F, sep = "\t", row.names = F, col.names = F)
  tmp.file_seq <- tempfile()
  system(paste0("bedtools getfasta -name -s -fi ", homo_sapiens_fasta_path, " -bed ",
                tmp.file, " -tab -fo ", tmp.file_seq))
  acceptor_sequences_input <- read.delim(tmp.file_seq, header = F)
  
  head(acceptor_sequences_input)
  head(donor_sequences_input)
  
  stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(acceptor_sequences_input$V1)))),
                      junc_tidy$junID %>% as.character()))
  junc_tidy <- cbind(junc_tidy,
                     acceptor_sequence = as.character(acceptor_sequences_input$V2))
  
  junc_tidy %>% head()
  
  
  ## Remove temporary files
  rm(to.BED, tmp.file, tmp.file_seq)
  
  
  
  ## 2. Generate the MaxEntScore --------------------------------------------------------------------
  
  
  ## get the sequences
  tmp.file <- tempfile()
  ## get the maxentscan for the 5' splice site
  
  ## check how many sequences contain "N"
  length(grep("N",as.character(junc_tidy$donor_sequence)))
  
  write.table(gsub("N","A",as.character(junc_tidy$donor_sequence)),file=tmp.file,row.names=F,col.names=F,quote=F)
  setwd(max_ent_tool_path)
  ss5score <- read.delim(pipe(paste0("perl ", max_ent_tool_path, "score5.pl ", tmp.file)),header = F)
  identical(as.character(ss5score$V1),gsub("N","A",as.character(junc_tidy$donor_sequence)))
  junc_tidy <- cbind(junc_tidy, ss5score = ss5score$V2)
  
  print("MaxEntScan score generated for the donor sequences!")
  
  
  ## get the maxentscan for the 3' splice site
  length(grep("N",as.character(junc_tidy$acceptor_sequence)))
  
  write.table(gsub("N","A",as.character(acceptor_sequences_input$V2)),file=tmp.file,row.names=F,col.names=F,quote=F)
  ss3score <- read.delim(pipe(paste0("perl ", max_ent_tool_path, "/score3.pl ", tmp.file)),header = F)
  identical(as.character(ss3score$V1),gsub("N","A",as.character(junc_tidy$acceptor_sequence)))
  junc_tidy <- cbind(junc_tidy, ss3score = ss3score$V2)
  
  print("MaxEntScan score generated for the acceptor sequences!")
  
  rm(ss5score, ss3score, tmp.file)
  
  junc_tidy[1,]
  
  
  
  return(junc_tidy)
}

#################################################################

get_and_tidy_recount3_raw_GTEx_split_reads <- function(projects_used) {
  
  # project_id <- projects_used[6]
  
  
  for (project_id in projects_used) {
    
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      type = "jxn"
    )
    
    all_split_reads <- data.frame(chr = rse %>% SummarizedExperiment::seqnames(),
                                  start = rse %>% SummarizedExperiment::start(),
                                  end = rse %>% SummarizedExperiment::end(),
                                  strand = rse %>% SummarizedExperiment::strand(),
                                  width = rse %>% SummarizedExperiment::width(),
                                  junID = rse %>% rownames())
    all_split_reads %>% nrow()
    saveRDS(object = all_split_reads %>% data.table::as.data.table(),
            file = paste0(folder_path, "/", "all_split_reads_raw.rds"))
  }
  
  
  
  ##########################################################
  ## Read all the split reads and return them by tissue
  ##########################################################
  
  all_split_reads_raw <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
    folder_path <- paste0(folder_root, "/raw_data/")
    
    
    print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))
    
    readRDS(file = paste0(folder_path, "/all_split_reads_raw.rds")) %>%
      return()
  })
  
  
  all_split_reads_raw_tidy <- all_split_reads_raw %>%
    mutate(junID_tidy = paste0(chr, ":", start, "-", end, ":", strand)) %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T)
  
  all_split_reads_raw_tidy %>% nrow()
  
  
  
  
  #######################################################################
  ## Remove split reads located in unplaced sequences in the chromosomes
  #######################################################################
  
  all_split_reads_raw_tidy_gr <- all_split_reads_raw_tidy %>%
    GenomicRanges::GRanges() %>%
    diffloop::rmchr()
  
  
  all_split_reads_raw_tidy_gr %>% length()
  all_split_reads_raw_tidy_gr <- GenomeInfoDb::keepSeqlevels(x = all_split_reads_raw_tidy_gr,
                                                             value = intersect(all_split_reads_raw_tidy_gr %>% 
                                                                                 GenomeInfoDb::seqnames() %>% 
                                                                                 levels(), 
                                                                               c( "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                                                  "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X", "Y")), 
                                                             pruning.mode = "tidy")
  all_split_reads_raw_tidy_gr %>% head()
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T) %>% 
    nrow()
  
  
  #######################################################
  ## Remove split reads overlapping the ENCODE backlist
  #######################################################
  
  source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline2_annotate_from_recount.R")
  blacklist_path <- "/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed"
  all_split_reads_raw_tidy_gr <- remove_encode_blacklist_regions(GRdata = all_split_reads_raw_tidy_gr,
                                                                 blacklist_path = blacklist_path)
  
  
  all_split_reads_raw_tidy_gr %>%
    as_tibble() %>%
    distinct(junID_tidy, .keep_all = T) %>% 
    nrow()
  
  
  #######################################################
  ## Anotate using 'dasper'
  #######################################################
  
  gtf_path <- paste0("/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  
  all_split_reads_details_105_w_symbol <- dasper::junction_annot(junctions = all_split_reads_raw_tidy_gr %>% GRanges(), 
                                                                 ref = edb)
  
  all_split_reads_details_105_w_symbol <- all_split_reads_details_105_w_symbol %>% 
    as_tibble()
  
  all_split_reads_details_105_w_symbol_reduced <- all_split_reads_details_105_w_symbol %>% 
    dplyr::select(seqnames, start, end, strand, junID, gene_id = gene_id_junction, in_ref, type )
  
  saveRDS(all_split_reads_details_105_w_symbol_reduced,
          file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_details_105_w_symbol.rds")
  
  
  
  
  ################################################################################
  ## Discard all junctions that are not annotated, novel donor or novel acceptor
  ################################################################################
  
  all_split_reads_details_105_w_symbol_reduced <- readRDS(file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_details_105_w_symbol.rds")
  
  all_split_reads_details_105_w_symbol_reduced_discard <- all_split_reads_details_105_w_symbol_reduced %>%
    dplyr::filter(!(type %in% c("annotated", "novel_donor", "novel_acceptor"))) %>%
    as.data.table()
  all_split_reads_details_105_w_symbol_reduced_discard %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  all_split_reads_details_105_w_symbol_reduced_discard %>%
    distinct(junID) %>%
    nrow()
  
  
  all_split_reads_details_105_w_symbol_reduced_keep <- all_split_reads_details_105_w_symbol_reduced %>%
    as.data.table() %>%
    dplyr::filter(!(junID %in% all_split_reads_details_105_w_symbol_reduced_discard$junID))
  
  ############################################
  ## Discard all junctions shorter than 25bp
  ############################################
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep %>%
    GRanges() %>%
    as_tibble()
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr  %>% 
    dplyr::filter(width < 25) %>%
    nrow()
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep_gr  %>% 
    dplyr::filter(width >=25)
  
  
  ############################################################################
  ## Discard all introns assigned to multiple genes (i.e. ambiguous introns)
  ############################################################################
  
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep_gr %>%
    distinct(junID, .keep_all = T) %>% 
    rowwise() %>%
    mutate(ambiguous = ifelse(gene_id %>% unlist() %>% length() > 1, T, F))
  
  ambiguous_introns <- all_split_reads_details_105_w_symbol_reduced_keep_gr %>%
    dplyr::filter(ambiguous == T)
  
  ambiguous_introns %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  
  saveRDS(object = ambiguous_introns,
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/database/all_ambiguous_introns.rds")
  
  all_split_reads_details_105_w_symbol_reduced_keep_gr <- all_split_reads_details_105_w_symbol_reduced_keep_gr %>%
    dplyr::filter(ambiguous == F) %>%
    dplyr::select(-ambiguous)
  
  saveRDS(object = all_split_reads_details_105_w_symbol_reduced_keep_gr,
          file = "~/PROJECTS/splicing-project-recount3/database/all_split_reads_details_105_w_symbol_reduced_keep.rds")
  
}


##################################################################





