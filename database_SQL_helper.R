################################################################
## UTILS
## functions to help in the QC data, add intronic features
#################################################################


get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}

remove_encode_blacklist_regions <- function(GRdata,
                                            blacklist_path) {
  
  
  if (!exists("encode_blacklist_hg38")) {
    encode_blacklist_hg38 <- rtracklayer::import(con = blacklist_path) %>% diffloop::rmchr()
  } else {
    print("'encode_blacklist_hg38' file already loaded!")
  }
  
  overlaped_junctions <- GenomicRanges::findOverlaps(query = encode_blacklist_hg38, 
                                                     subject = GRdata %>% diffloop::rmchr(),
                                                     #type = "any",
                                                     ignore.strand = F)
  
  ## JuncID indexes to be removed: they overlap with a black region
  indexes <- S4Vectors::subjectHits(overlaped_junctions)
  
  if (length(indexes) > 0) {
    print(paste0(length(unique(indexes)), " junctions overlap with a ENCODE blacklist region"))
    GRdata <- GRdata[-indexes, ]
    print(paste0(length(GRdata), " junctions after removing overlaps with ENCODE BlackList regions!"))
  }else{
    print("No junctions overlapping with an ENCODE blacklist region")
  }
  
  ## Return tidied data
  return(GRdata)
}

generate_max_ent_score <- function(junc_tidy,
                                   max_ent_tool_path,
                                   homo_sapiens_fasta_path){
  
  library(Biostrings)
  library(tidyverse)
  library(protr)
  
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
  
  
  #stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(donor_sequences_input$V1) ))),
  #                    junc_tidy$junID %>% as.character()))
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
  
  #stopifnot(identical(gsub("\\(\\+\\)", "", gsub("\\(\\*\\)", "", gsub("\\(-\\)", "", as.character(acceptor_sequences_input$V1)))),
  #                    junc_tidy$junID %>% as.character()))
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

add_cdts_cons_scores <- function(cluster = NULL,
                                 db_introns = NULL,
                                 folder_name = NULL) {
  
  
  if (is.null(db_introns) && !is.null(cluster)) {
    ## Load the IDB 
    db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
  }
  
  db_introns <- db_introns %>%
    mutate(CDTS_5ss_mean = 0.0,
           CDTS_3ss_mean = 0.0,
           phastCons20way_5ss_mean = 0.0,
           phastCons20way_3ss_mean = 0.0) %>%
    GRanges()
  
  
  print(paste0(Sys.time(), " - getting 5' scores assigned to introns from IntroVerse ..."))
  
  ## https://www.nature.com/articles/nature09000
  ## Scores from the 5'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% start() - 5,
                                                                                            end = db_introns %>% start() + 35),
                                                                           strand = db_introns %>% strand()),
                                          ignore.strand = FALSE,
                                          type = "any")
  
  overlaps_tidy <- overlaps %>%
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_5ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_5ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  
  
  print(paste0(Sys.time(), " - getting 3' scores assigned to introns from IntroVerse ..."))
  
  ## Scores from the 3'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% end() - 35,
                                                                                            end = db_introns %>% end() + 5),
                                                                           strand = db_introns %>% strand()),
                                          
                                          ignore.strand = FALSE,
                                          type = "any")
  
  overlaps_tidy <- overlaps %>% 
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    as_tibble() %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_3ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_3ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  #####################
  ## SAVE RESULTS
  #####################
  
  if (!is.null(cluster)) {
    file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
    saveRDS(object = db_introns %>% data.table::as.data.table(),
            file = file_name)
    
    print(paste0(Sys.time(), " - CDTS and Conservation scores added! IDB updated!"))
    
    
    rm(overlaps_tidy)
    rm(db_introns)
    rm(file_name)
    rm(overlaps)
    
  } else {
    return(db_introns)
  }
  
  
  
  
  #gc()
}

remove_MT_genes <- function(cluster,
                            folder_name) {
  
  
  
  ## Load mis-splicing ratios df
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  df_introns %>% head()
  df_introns %>% nrow()
  
  
  ## Load MT genes
  
  print(paste0(Sys.time(), " - checking the existance of MT genes..."))
  
  MT_genes <- readRDS(file = paste0(dependencies_folder,
                                    "/MT_genes/MT_genes.rds"))
  MT_geneID <- MT_genes %>% distinct(gene_id) %>% pull(gene_id)
  
  if (any(df_introns$gene_id %in% MT_geneID)) {
    
    ## Remove MT genes
    df_introns <- df_introns %>%
      dplyr::filter(!(gene_id %in% MT_geneID)) 
    
    print(paste0(Sys.time(), " - MT genes removed! ", df_introns %>% nrow(), " final number of junctions."))
    
    
    ## SAVE RESULTS
    saveRDS(object = df_introns,
            file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
    
    print(paste0(Sys.time(), " - results saved!"))
    
  } else {
    print(paste0(Sys.time(), " - the dataset doesn't contain any MT gene!"))
  }
  
  
  ## FREE SOME MEMORY
  rm(df_introns)
  rm(MT_geneID)
  rm(MT_genes)
  gc()
  
}

annotate_contamination_rates <- function(gtf_versions,
                                         splicing_projects) {
  
  gtf_versions <- c("76", "81", "90", "104")
  splicing_projects <- "BRAIN"
  # gtf_versions <- c("97")
  
  for (gtf_version in gtf_versions) {
    
    # gtf_version <- gtf_versions[1]
    
    ## re-annotate them
    if ( file.exists(paste0("/data/references/ensembl/gtf_gff3/v", gtf_version, "/Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")) ) {
      gtf_path <- paste0("/data/references/ensembl/gtf_gff3/v", gtf_version, "/Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")
    } else {
      gtf_path <- paste0("/data/references/ensembl/gtf_gff3/v", gtf_version, "/Homo_sapiens.GRCh38.", gtf_version, ".gtf")
    }
    
    edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    
    ## Per each project
    for (project_id in splicing_projects) {
      
      # project_id <- splicing_projects[1]
      
      folder_results <- paste0(getwd(), "/results/", project_id, "/v105/", main_project, "/base_data/")
      
      all_clusters <- readRDS(file = paste0(folder_results, "/", project_id, "_clusters_used.rds"))
      
      
      ## Per each cluster
      for ( cluster_id in all_clusters[6] ) {
        
        # cluster_id <- all_clusters[1]
        
        ## get the QC split reads
        all_split_reads <- readRDS(file = paste0(folder_results, "/", project_id, "_", cluster_id, "_all_split_reads.rds"))
        
        
        
        
        all_split_reads_details_w_symbol <- dasper::junction_annot(junctions = all_split_reads %>%
                                                                     dplyr::select(junID, seqnames, start, end, strand) %>%
                                                                     GRanges(), 
                                                                   ref = edb)
        
        
        folder_path <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/")
        dir.create(file.path(folder_path), recursive = TRUE, showWarnings = T)
        saveRDS(object = all_split_reads_details_w_symbol,
                file = paste0(folder_path, "/", project_id, "_", cluster_id, "_all_split_reads.rds") )
      }
      
    }
    
    rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    gc()
    
  }
  
}


#################################################################
## SQL HELPER
## auxiliary functions to help in the generation of a database
#################################################################



get_all_annotated_split_reads <- function(projects_used,
                                          gtf_version,
                                          all_clusters = NULL,
                                          main_project) {
  
  
  
  #############################################
  ## These are all split reads from all tissues
  ## obtained from 'USE ME' SAMPLES and samples with more than 6 RIN
  #############################################
  
  ## The sample filter in IntroVerse corresponded to 
  
  all_split_reads_details_all_tissues <- map_df(projects_used, function(project_id) {
    
    # project_id <- projects_used[1]
    # project_id <- "BONE_MARROW"
    
    folder_root <- paste0(getwd(), "/results/", project_id, "/v", 
                          gtf_version, "/", main_project, "/")
    
    if ( is.null(all_clusters) ) {
      ## Read all clusters considered from the current tissue
      metadata.info <- readRDS(file = paste0(folder_root, "/base_data/",
                                             project_id, "_samples_metadata.rds"))
      all_clusters <-  metadata.info %>%
        distinct(gtex.smtsd) %>%
        pull()

      if ( any(metadata.info$gtex.smrin < 6) ) {
        print("ERROR! The samples used have not been filtered adequately.")
        break;
      }
    }
   
    
    all_jxn_qc <- map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[1]
      print(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
      
      if ( file.exists(paste0(folder_root, "/base_data/", 
                              project_id, "_", cluster, "_all_split_reads.rds")) ) {
        
        all_split_reads_details_105 <- readRDS(file = paste0(folder_root, "/base_data/", project_id, "_",
                                                             cluster, "_all_split_reads.rds"))
        
        all_split_reads_details_tidy <- all_split_reads_details_105 %>% 
          distinct(junID, .keep_all = T) %>% 
          as_tibble() %>%
          return()
        
      } else {
        return(NULL)
      }
      
    })
    
    if (all_jxn_qc %>% nrow() > 0 ) {
      
      all_jxn_qc %>%
        distinct(junID, .keep_all = T) %>% 
        return()
      
    } else {
      return(NULL)
    }
    
  })
  
  
  print(paste0(Sys.time(), " - saving 'all_annotated_split_reads' for the database!"))
  
  all_split_reads_details_all_tissues <- all_split_reads_details_all_tissues %>%
    distinct(junID, .keep_all = T)
  
  database_folder <- paste0(getwd(), "/database/v", gtf_version, "/", main_project, "/")
  dir.create(file.path(database_folder), showWarnings = T, recursive = T)
  
  ## Save data
  saveRDS(object = all_split_reads_details_all_tissues,
          file = paste0(database_folder, "/all_split_reads_details_all_tissues.rds") )
  
}


get_all_raw_distances_pairings <- function(projects_used,
                                           gtf_version,
                                           all_clusters = NULL,
                                           main_project) {
  
  
  ## LOOP THROUGH PROJECTS
  df_all_distances_pairings_raw <- map_df(projects_used, function(project_id) {
    
    # project_id <- projects_used[1]
    # project_id <- "KIDNEY"
    
    print(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    folder_root <- paste0(getwd(), "/results/", 
                          project_id, "/v", gtf_version, "/", main_project, "/")
    
    if ( is.null(all_clusters) ) {
      ## Read all clusters considered from the current tissue
      metadata.info <- readRDS(file = paste0(folder_root, "/base_data/",
                                             project_id, "_samples_metadata.rds"))
      all_clusters <-  metadata.info %>%
        distinct(gtex.smtsd) %>%
        pull()
  
      if ( any(metadata.info$gtex.smrin < 6) ) {
        print("ERROR! The samples used have not been filtered adequately.")
        break;
      }
    }
    
    map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[1]
      
      print(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
      
      ## Load samples
      if ( file.exists(paste0(folder_root, "/base_data/", project_id, "_", cluster, "_samples_used.rds")) ) {
        
        samples <- readRDS(file = paste0(folder_root, "/base_data/", project_id, "_", cluster,  "_samples_used.rds"))
        
        if ( samples %>% length() > 0 ) {
          
          folder_name <- paste0(folder_root, "/results/", cluster, "/")
          
          if ( !file.exists(paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) ) {
            
            ## Obtain the distances across all samples
            df_all <- map_df(samples, function(sample) { 
              
              # sample <- samples[1]
              
              file_name <- paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")
              
              
              if (file.exists(file_name)) {
                print(paste0(cluster, " - ", sample))
                df <- readRDS(file = file_name)
                
                return(df)
              } else {
                return(NULL)
              }
              
            })
            
            if ( nrow(df_all) > 0 ) {
              saveRDS(object = df_all %>%
                        distinct(novel_junID, ref_junID, .keep_all = T) %>%
                        mutate(tissue = cluster),
                      file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
            }
            
          } else {
            print(paste0("File '", cluster, "_raw_distances_tidy.rds' already exists!"))
            df_all <- readRDS( file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds") )
          }
          
          
          if ( nrow(df_all) > 0 ) {
            
            df_all %>%
              distinct(novel_junID, ref_junID, .keep_all = T) %>%
              mutate(project = project_id) %>%
              return()
            
          } else {
            return(NULL)
          }          
          # df_all2 <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
          
        } else {
          print(paste0("ERROR: no samples available for the tissue: ", project_id))
          return(NULL)
        }
        
      } else {
        return(NULL)
      }
    })  
  })  
  
  
  print(paste0(Sys.time(), " - saving 'df_all_distances_pairings' for the database!"))
  
  
  saveRDS(object = df_all_distances_pairings_raw %>% distinct(project) %>% pull(),
          file = paste0(getwd(),"/results/",main_project,"_final_projects_used.rds"))
  
  saveRDS(object = df_all_distances_pairings_raw %>%
            distinct(novel_junID, ref_junID, .keep_all = T),
          file = paste0(getwd(),"/database/v", gtf_version, "/",
                        main_project, "/all_distances_pairings_all_tissues.rds"))
  
}

get_intron_never_misspliced <- function (projects_used,
                                         all_clusters,
                                         main_project) {
  
  
  df_never <- map_df(projects_used, function(project_id) {
    
    # project_id <- projects_used[1]
    
    print(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    base_folder <- paste0(getwd(), "/results/", 
                          project_id, "/v", gtf_version, "/", main_project, "/")
    
    
    
    if ( file.exists(paste0(base_folder, "/base_data/", 
                            project_id, "_samples_metadata.rds")) ) {
      metadata.info <- readRDS(file = paste0(base_folder, "/base_data/", 
                                             project_id, "_samples_metadata.rds"))
      all_clusters <-  metadata.info %>% 
        distinct(gtex.smtsd) %>% 
        pull()
      
    }
    
    
    map_df(all_clusters, function(cluster) { 
      
      # cluster <- all_clusters[1]
      
      # print(paste0(Sys.time(), " --> ", cluster))
      if ( file.exists(paste0(base_folder, "results/", cluster, 
                              "/not-misspliced/", cluster, "_all_notmisspliced.rds")) ) {
        df_introns_never <- readRDS(file = paste0(base_folder, "results/", cluster, 
                                                  "/not-misspliced/", 
                                                  cluster, "_all_notmisspliced.rds")) %>% as_tibble()
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


generate_transcript_biotype_percentage <- function(projects_used,
                                                   homo_sapiens_v105_path,
                                                   main_project,
                                                   gtf_version) {
  
  
  #######################################
  ## GET THE TRANSCRIPT BIOTYPE
  #######################################
  
  print(paste0(Sys.time(), " - loading the human reference transcriptome ... "))
  
  ## Import HUMAN REFERENCE transcriptome
  homo_sapiens_v105 <- rtracklayer::import(con = homo_sapiens_v105_path) %>% 
    as_tibble()
  
  ## Get v105 transcripts
  transcripts_v105 <- homo_sapiens_v105 %>%
    filter(type == "transcript") %>% 
    dplyr::select(transcript_id, transcript_biotype, gene_id)
  
  transcripts_v105 %>% head()
  
  
  #######################################
  ## LOAD ALL SPLIT READS ALL TISSUES
  #######################################
  
  print(paste0(Sys.time(), " - loading the recount3 GTEx split reads ... "))
  
  ## LOAD the all split reads from all recount3 GTEx projects
  database_folder <- paste0(getwd(), "/database/v", gtf_version, "/", main_project, "/")
  all_split_reads_details_all_tissues <- readRDS(file = paste0(database_folder, "/all_split_reads_details_all_tissues.rds") )
  
  all_split_reads_details_all_tissues %>% head()
  all_split_reads_details_all_tissues %>% nrow()
  
  
  
  
  #######################################
  ## EXPLORE AND TIDY THE RESULT
  #######################################
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_all_tissues$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_all_tissues[ind, "junID"] <- str_replace(string = all_split_reads_details_all_tissues[ind, "junID"]$junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = all_split_reads_details_all_tissues[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_all_tissues$junID, pattern = "\\*")) %>% print()
  }
  
  print(object.size(all_split_reads_details_all_tissues), units = "Gb")
  
  ## Merge datasets to add transcript biotype
  print(paste0(Sys.time(), " --> adding transcript biotype..."))
  
  ## Merging using data table structures saves time
  df_all_junctions <- all_split_reads_details_all_tissues %>% unnest(tx_id_junction) %>% data.table::as.data.table()
  transcripts_v105 <- transcripts_v105 %>% data.table::as.data.table()
  
  df_all_junctions %>% head()
  transcripts_v105 %>% head()
  
  df_all_junctions_tx <- df_all_junctions %>% 
    left_join(y = transcripts_v105,
              by = c("tx_id_junction" = "transcript_id"))
  
  print(object.size(df_all_junctions_tx), units = "Gb")
  
  
  #######################################
  ## CALCULATE THE TRANSCRIPT PERCENTAGE
  #######################################
  
  
  print(paste0(Sys.time(), " --> starting protein-coding percentage calculation ..."))
  print(paste0(Sys.time(), " --> ", df_all_junctions$junID %>% unique() %>% length(), " total number of junctions."))
  
  
  ## Only keep not ambiguous jxn (i.e. junctions belonging to multiple genes)
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
          file = paste0(getwd(), "/results/all_split_reads_all_tissues_all_biotypes.rds"))
  
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
  
  
  df_all_percentage_tidy_merged <- merge(x = df_all_percentage_tidy_PC %>% dplyr::rename(percent_PC = percent) %>% data.table::as.data.table(),
                                         y = df_all_percentage_tidy_lncRNA %>% dplyr::rename(percent_lncRNA = percent) %>% data.table::as.data.table(),
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
          file = paste0(getwd(), "/results/all_split_reads_all_tissues_PC_biotype.rds"))
  print(paste0(Sys.time(), " - files saved!"))
  
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

calculate_tpm <- function(rse, ref_tidy) {
  
  
  # Remove anything after . in ensembl id
  rownames(rse) <- rownames(rse) %>% 
    str_remove("\\..*")
  
  
  # Convert to tpm, which is calculated by:
  # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
  # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
  # 3. Divide the RPK values by the “per million” scaling factor.
  
  
  srp_rpk <- 
    rse %>% 
    SummarizedExperiment::assay() %>%
    as_tibble(rownames = "gene") %>% 
    tidyr::pivot_longer( ## equivalent to gather
      cols = -c("gene"),
      names_to = "recount_id",
      values_to = "counts"
    ) %>% 
    dplyr::inner_join(
      ref_tidy %>% 
        as_tibble() %>% 
        dplyr::select(gene_id, width),
      by = c("gene" = "gene_id")
    ) %>% # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
    dplyr::mutate(
      rpk = counts/width
    ) 
  srp_rpk %>% head()
  
  
  tpm <- 
    srp_rpk %>% 
    # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
    dplyr::inner_join(
      srp_rpk %>% 
        dplyr::group_by(recount_id) %>% 
        dplyr::summarise(
          scaling_factor = sum(rpk)/1e6
        ),
      by = "recount_id"
    ) %>% # 3. Divide the RPK values by the “per million” scaling factor.
    dplyr::mutate(
      tpm = rpk/scaling_factor
    )   
  
  tpm <- tpm %>%
    dplyr::select(gene, recount_id, tpm) %>% 
    spread(key = recount_id, value = tpm)
  
  return(tpm)
}


generate_recount3_tpm <- function(projects_used,
                                  ref = NULL,
                                  gtf_version,
                                  main_project) {
  
  # Reference gtf
  if (is.null(ref)) {
    ref <- rtracklayer::import(con = paste0(dependencies_folder, "/Homo_sapiens.GRCh38.",gtf_version,".chr.gtf"))
  }
  ref <- ref %>% GenomeInfoDb::keepSeqlevels(c(1:22), pruning.mode = "coarse") 
  ref %>% head()
  
  ## Per each gene, get its length
  ref_tidy <- ref %>%
    as_tibble() %>% 
    filter(type == "gene") %>%
    distinct(gene_id, .keep_all = T) %>%
    dplyr::select(gene_id, gene_name, width) 
  
  for (project_id in projects_used) {
    
    # project_id <- projects_used[1]
    
    ## 1. Get expression data from recount3 
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      type = "gene")
    
    ## 2. Transform raw counts
    SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
    
    ## 3. Calculate TPM using the length of the coding sequence
    tpm <- calculate_tpm(rse, ref_tidy)
    tpm %>% head() %>% print()
    
    ## 4. Free up memory and save
    print(paste0(Sys.time(), " - ", project_id, " finished!"))
    dir_path <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/",
                       main_project, "/tpm/")
    dir.create(path = dir_path, showWarnings = T, recursive = T)
    saveRDS(object = tpm,
            file = paste0(dir_path, "/",  project_id, "_tpm.rds"))
    
    
    rm(rse)
    gc()
    
  }
  
  rm(ref_tidy)
  gc()
  
}



filter_recount3_tpm <- function(projects_used,
                                gtf_version,
                                all_clusters = NULL,
                                main_project) {
  
  
  for (project_id in projects_used) {
    
    
    # project_id <- projects_used[1]
    
    print(paste0(Sys.time(), " - calculating TPM values for '", project_id, "'"))
    
    
    folder_root <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/")
    recount_tpm <- readRDS(file = paste0(folder_root, "/tpm_all_samples.rds"))
    
    
    ## 2. For each tissue within the current project, filter the RSE by its samples
    folder_root <- paste0(folder_root, "/", main_project, "/")
    
    
    
    if ( is.null(all_clusters) && file.exists(paste0(folder_root, "/base_data/", project_id, "_clusters_used.rds")) ) {
      all_clusters <- readRDS(paste0(folder_root, "/base_data/", project_id, "_clusters_used.rds"))
    } 
    
    if ( !is.null(all_clusters) ) {
      
      for ( cluster_id in all_clusters ) {

        # cluster_id <- all_clusters[1]
        
        print(paste0(Sys.time(), " - '", cluster_id, "' ..."))

        cluster_samples <- NULL
        
        if ( file.exists(paste0(folder_root, "/base_data/", project_id, "_", cluster_id, "_samples_used.rds"))) {
            
            cluster_samples <- readRDS(file = paste0(folder_root, "/base_data/",
                                                     project_id, "_", cluster_id, "_samples_used.rds"))  
          
        }

        if ( !is.null(cluster_samples) ) {


          ## Filter the object for the set of samples corresponding to the current cluster
          recount_tpm_local <- recount_tpm %>%
            dplyr::select(c("gene", all_of(cluster_samples)))

          ## Save results
          folder_name <- paste0(folder_root, "/tpm/")
          dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)
          saveRDS(object = recount_tpm_local,
                  file = paste0(folder_name, project_id, "_", cluster_id, "_tpm.rds"))


          rm(recount_tpm_local)
          rm(cluster_samples)
          gc()

        }

      }
    }
  }
  
}
  
update_table_tpm <- function(database_path,
                             projects_used,
                             gtf_version,
                             main_project) {
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  
  for (project_id in all_projects) {
    
    # project_id <- all_projects[15]
    
    print(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    base_folder <- paste0(getwd(),"/results/", project_id, "/v", gtf_version, "/", main_project, "/")
    
    clusters <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version,
                                      "/", main_project, "/base_data/", project_id, "_clusters_used.rds"))
    
    for (cluster_id in clusters) { 
      
      # cluster_id <- clusters[1]
      print(paste0(Sys.time(), " --> ", cluster_id))
      
      if ( file.exists( paste0(base_folder, "results/tpm/",
                               "/", project_id, "_", cluster_id, "_tpm.rds")) ) {
        
        samples_used <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version, 
                                              "/", main_project, "/base_data/", project_id, "_", cluster_id, "_samples_used.rds"))
        
        
        query <- paste0("SELECT tissue.*, gene.gene_id FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue
                        INNER JOIN 'transcript' ON transcript.id = tissue.transcript_id 
                        INNER JOIN 'gene' ON gene.id = transcript.gene_id;")
        df_misspliced <- dbGetQuery(con, query) %>% as_tibble()
        
        
        tpm <- readRDS(file = paste0(base_folder, "results/tpm/",
                                     "/", project_id, "_", cluster_id, "_tpm.rds")) %>% 
          dplyr::select(gene_id = gene, all_of(samples_used))
        
        
        tpm <- tpm  %>%
          dplyr::mutate(tpm_median = matrixStats::rowMedians(x = as.matrix(.[2:(ncol(tpm))]))) %>%
          dplyr::select(gene_id, tpm_median) 
        
        
        ## In case there are any duplicates, take the genes with the maximum tpms
        tpm <- tpm %>% 
          group_by(gene_id) %>% 
          summarize_all(max) %>%
          ungroup()
        
        tpm_tidy <- tpm %>%
          inner_join(y = df_misspliced %>% 
                       as_tibble(),
                     by = c("gene_id" = "gene_id")) %>%
          inner_join(y = df_transcript %>% 
                       as_tibble(),
                     by = c("id" = "gene_id")) %>%
          dplyr::select(transcript_id = id.y, 
                        tpm_median)
        
        df_misspliced_tidy <- df_misspliced %>%
          left_join(y = tpm,
                    by = c("gene_id")) %>% 
          dplyr::select(-gene_id) %>%
          dplyr::rename(gene_tpm = tpm_median)
      }
    }
  }
  
}


##################################################################
## DATABASE GENERATION
##################################################################

tidy_data_pior_sql <- function (projects_used, 
                                gtf_version,
                                all_clusters = NULL,
                                main_project) {
  
  
  print(paste0(Sys.time(), " - loading base data ..."))
  
  ## Load base recount3 object containing the GTEx v8 split reads passing the QC criteria
  all_split_reads_details_w_symbol_reduced_keep_gr <- readRDS(file = paste0(getwd(), "/database/v", gtf_version, 
                                                                            "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))
  all_split_reads_details_w_symbol_reduced_keep_gr %>% nrow()
  
  ############################################
  ## Discard all junctions from 
  ## EXCLUDE ME samples
  ############################################
  
  all_split_reads_details_all_tissues <- readRDS(file = paste0(getwd(), "/database/v", gtf_version, "/", main_project, 
                                                               "/all_split_reads_details_all_tissues.rds"))
  
  all_split_reads_details_all_tissues %>% nrow()
  all_split_reads_details_all_tissues %>% head()
  
  ## This should be zero
  if ( setdiff(all_split_reads_details_all_tissues$junID, 
               all_split_reads_details_w_symbol_reduced_keep_gr$junID) %>% 
       length() > 0) {
    print("ERROR! Some of the annotated split reads are not found within the base tidied recount3 object.")
    break;
  }
  
  
  # ## These are the junctions from EXCLUDE ME samples
  setdiff(all_split_reads_details_w_symbol_reduced_keep_gr$junID, 
          all_split_reads_details_all_tissues$junID) %>% unique %>% length()
  
  # ## Stats
  # all_split_reads_details %>% as.data.table() %>% dplyr::count(type)
  # all_split_reads_details %>% distinct(junID) %>% nrow()
  
  ############################################
  ## QC
  ############################################ 
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_w_symbol_reduced_keep_gr$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_w_symbol_reduced_keep_gr[ind, "junID"] <- str_replace(string = all_split_reads_details_w_symbol_reduced_keep_gr[ind, "junID"]$junID, 
                                                                                  pattern = "\\*", 
                                                                                  replacement = all_split_reads_details_w_symbol_reduced_keep_gr[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_w_symbol_reduced_keep_gr$junID, pattern = "\\*")) %>% print()
  }
  
  
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = all_split_reads_details_all_tissues$junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    all_split_reads_details_all_tissues[ind, "junID"] <- str_replace(string = all_split_reads_details_all_tissues[ind, "junID"]$junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = all_split_reads_details_all_tissues[ind, "strand"]$strand %>% as.character() )
    any(str_detect(all_split_reads_details_all_tissues$junID, pattern = "\\*")) %>% print()
  }
  
  
  ##########################################
  ## Load all distances pairings
  ##########################################
  
  df_all_distances_pairings_raw <- readRDS(file = paste0(getwd(), "/database/v", gtf_version, "/",
                                                         main_project, "/all_distances_pairings_all_tissues.rds"))
  
  ## QC
  ## Remove potential * in the junID of the reference introns
  ind <- which(str_detect(string = df_all_distances_pairings_raw$ref_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_distances_pairings_raw[ind, "ref_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "ref_junID"]$ref_junID, 
                                                                   pattern = "\\*", 
                                                                   replacement = df_all_distances_pairings_raw[ind, "ref_strand"]$ref_strand %>% as.character())
    any(str_detect(df_all_distances_pairings_raw$ref_junID, pattern = "\\*")) %>% print()
  }
  ## Remove potential * in the junID of the novel junctions
  ind <- which(str_detect(string = df_all_distances_pairings_raw$novel_junID, pattern = "\\*"))
  if (ind %>% length() > 0) {
    df_all_distances_pairings_raw[ind, "novel_junID"] <- str_replace(string = df_all_distances_pairings_raw[ind, "novel_junID"]$novel_junID, 
                                                                     pattern = "\\*", 
                                                                     replacement = df_all_distances_pairings_raw[ind, "novel_strand"]$novel_strand  %>% as.character())
    any(str_detect(df_all_distances_pairings_raw$novel_junID, pattern = "\\*")) %>% print()
  }
  
  
  ##########################################
  ## Get never mis-spliced
  #########################################
  
  df_never_misspliced <- get_intron_never_misspliced(projects_used = projects_used,
                                                     all_clusters = all_clusters,
                                                     main_project = main_project)
  
  if ( any(str_detect(df_never_misspliced$ref_junID, pattern = "\\*")) ) {
    
    print("ERROR! Some junctions contain a '*' in their IDs")
    break;
  }
  
  ## Remove the introns paired with novel junctions (i.e. mis-spliced)
  df_never_misspliced_tidy <- df_never_misspliced %>%
    dplyr::filter(!(ref_junID %in% df_all_distances_pairings_raw$ref_junID)) %>%
    as_tibble()
  
  df_never_misspliced_tidy %>% distinct(ref_junID) %>% as_tibble()
  
  
  ############################################
  ## GET all not paired
  ############################################ 
  
  # df_all_novel_raw_tidy
  df_not_paired <- all_split_reads_details_all_tissues %>%
    data.table::as.data.table() %>%
    dplyr::filter(!(junID %in% c(df_all_distances_pairings_raw$ref_junID,
                                 df_all_distances_pairings_raw$novel_junID)))
  
  ## These are all the non-paired, including the never mis-spliced. Thus, GTEx 'splicing' data: [768,646 - 56,090 = 730125]
  df_not_paired %>% 
    distinct(junID)
  
  ## Check them by category
  df_not_paired %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  
  
  if (any(str_detect(df_never_misspliced_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(df_not_paired$junID, pattern = "\\*"))) {
    print("ERROR!")
  }
  
  ## All never mis-spliced should be categorised as not paired.
  ## Thus, this should be zero
  setdiff(df_never_misspliced_tidy$ref_junID, 
          df_not_paired %>% dplyr::filter(type == "annotated") %>% pull(junID))
  
  df_not_paired_tidy <- df_not_paired %>%
    dplyr::filter(!(junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    distinct(junID, .keep_all = T) %>%
    as_tibble()
  
  df_not_paired_tidy %>%
    distinct(junID, .keep_all = T) %>%
    dplyr::count(type)
  df_not_paired_tidy %>% distinct(junID)
  
  
  ## This should be zero
  if ( intersect(df_not_paired_tidy$junID, 
                 df_all_distances_pairings_raw$novel_junID) %>% length() > 0 ) {
    print("ERROR!")
  }
  
  
  df_all_distances_pairings_raw %>% distinct(novel_junID) %>% nrow() +
    df_all_distances_pairings_raw %>% distinct(ref_junID) %>% nrow() +
    df_never_misspliced_tidy %>% distinct(ref_junID) %>% nrow()
  
  ##########################################
  ## Remove ambiguous junctions
  ##########################################
  
  
  ## All these should be zero
  if( intersect(df_not_paired_tidy$junID, df_all_distances_pairings_raw$novel_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_all_distances_pairings_raw$ref_junID) %>% length() > 0 |
      intersect(df_not_paired_tidy$junID, df_never_misspliced_tidy$ref_junID) %>% length() > 0) {
    print("ERROR!")
  }
  
  
  ## 1. Obtain the ambiguous junctions
  
  df_ambiguous_novel <- df_all_distances_pairings_raw %>%
    dplyr::filter(!(novel_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_not_paired_tidy$junID),
                  !(ref_junID %in% df_never_misspliced_tidy$ref_junID)) %>%
    dplyr::group_by(novel_junID) %>%
    mutate(distances_sd = distance %>% sd()) %>%
    dplyr::filter(distances_sd > 0)
  
  
  df_ambiguous_novel %>% ungroup() %>% distinct(novel_junID)
  df_ambiguous_novel %>% ungroup() %>% distinct(ref_junID)
  
  
  ## 2. Remove ambiguous junctions
  
  df_all_distances_pairings_raw_tidy <- df_all_distances_pairings_raw %>%
    dplyr::filter(!(novel_junID %in% df_ambiguous_novel$novel_junID)) %>%
    data.table::as.data.table() %>%
    distinct(novel_junID, ref_junID, .keep_all = T) %>%
    mutate(ref_strand = ref_strand %>% as.character(),
           novel_strand = novel_strand %>% as.character()) 
  
  
  ## Some numbers for the graph
  df_all_distances_pairings_raw_tidy %>%
    dplyr::distinct(novel_junID)
  df_all_distances_pairings_raw_tidy %>%
    dplyr::distinct(ref_junID)
  
  
  (df_ambiguous_novel %>% 
      ungroup() %>%
      distinct(ref_junID) %>% nrow()) - (intersect(c(df_all_distances_pairings_raw_tidy$novel_junID,
                                                     df_all_distances_pairings_raw_tidy$ref_junID),
                                                   df_ambiguous_novel %>% 
                                                     ungroup() %>%
                                                     distinct(ref_junID) %>% pull) %>% length())
  
  (df_all_distances_pairings_raw_tidy %>%
      dplyr::distinct(novel_junID) %>% 
      nrow()) + 
    (df_all_distances_pairings_raw_tidy %>%
       dplyr::distinct(ref_junID) %>% 
       nrow()) + 
    (df_never_misspliced_tidy %>% 
       distinct(ref_junID) %>% 
       nrow())
  
  
  ##############################################################################
  ## SAVE FINAL OBJECT
  ##############################################################################
  
  ## 1. DISTANCES PAIRINGS
  if (any(str_detect(string = df_all_distances_pairings_raw_tidy$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_all_distances_pairings_raw_tidy$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  df_all_distances_pairings_raw_tidy <- df_all_distances_pairings_raw_tidy %>%
    inner_join(y = all_split_reads_details_all_tissues %>% dplyr::select(junID, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  saveRDS(object = df_all_distances_pairings_raw_tidy,
          file = paste0(getwd(), "/database/v", gtf_version, "/",
                        main_project, "/all_distances_correct_pairings_all_tissues.rds"))
  
  
  ## 2. NEVER MIS-SPLICED
  if (any(str_detect(string = df_never_misspliced_tidy$ref_junID, pattern = "\\*")) ) {
    print("ERROR! Some NEVER MIS-SPLICED junctions still have a '*' in their IDs!")
  }
  df_never_misspliced_tidy <- df_never_misspliced_tidy %>%
    inner_join(y = all_split_reads_details_all_tissues %>% 
                 dplyr::select(junID, seqnames, start, end, width, strand, gene_id, tx_id_junction),
               by = c("ref_junID" = "junID"))
  
  saveRDS(object = df_never_misspliced_tidy,
          file = paste0(getwd(), "/database/v", gtf_version, "/",
                        main_project, "/all_nevermisspliced_introns_all_tissues.rds"))
  
  
  ## 3. AMBIGUOUS JUNCTIONS
  if (any(str_detect(string = df_ambiguous_novel$ref_junID, pattern = "\\*")) |
      any(str_detect(string = df_ambiguous_novel$novel_junID, pattern = "\\*")) ) {
    print("ERROR! Some junctions still have a * in their IDs!")
  }
  saveRDS(df_ambiguous_novel,
          file = paste0(getwd(), "/database/v", gtf_version, "/", main_project, 
                        "/all_distances_ambiguous_pairings_all_tissues.rds"))
}


sql_database_generation <- function(database_path,
                                    projects_used, 
                                    main_project,
                                    gtf_version,
                                    remove_all = NULL) {
  
  
  if ( !is.null(remove_all) ) {
    
    remove_tables(database_path, remove_all)
  } 
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  tables <- DBI::dbListTables(conn = con)
  
  
  if ( !any(tables == 'master') ) {
    create_metadata_table(database_path,
                          main_project = main_project,
                          gtf_version = gtf_version,
                          all_projects = projects_used)
  }
  
  
  if ( ! any(tables %in% c('intron', 'novel', 'gene', 'transcript')) ) {
    create_master_tables(database_path,
                         main_project = main_project,
                         gtf_version = gtf_version)
  }
  
  print(paste0(Sys.time(), " - creating cluster tables ..."))
  
  rm(CNC_CDTS_CONS_gr)
  rm(hg_mane_transcripts)
  rm(hg_MANE_tidy)
  rm(hg_MANE)
  rm(hg38_transcripts)
  rm(hg38)
  gc()
  
  create_cluster_tables(database_path = database_path,
                        gtf_version = gtf_version,
                        all_projects = projects_used,
                        main_project = main_project)
  
}
