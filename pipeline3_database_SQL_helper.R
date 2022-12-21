library(tidyverse)
library(data.table)
library(GenomicRanges)
library(DBI)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/sql_helper.R")


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
                                   # bedtools_path = NULL,
                                   max_ent_tool_path,
                                   homo_sapiens_fasta_path){
  
  library(Biostrings)
  library(tidyverse)
  library(protr)
  
  # ## Load the annotated split reads data (both for PD and control samples)
  # junc_tidy <-readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount2-projects/data/",
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
  
  MT_genes <- readRDS(file = "/data/references/MT_genes/MT_genes.rds")
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

generate_recount3_median_tpm <- function(all_projects,
                                         gtf_version,
                                         main_project) {
  
  
  
  
  for (project_id in all_projects) {
    
    # project_id <- all_projects[1]
    # project_id <- "KIDNEY"
    
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      type = "gene")
    
    SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
    
    SummarizedExperiment::assays(rse)$TPM <- recount::getTPM(rse)
    
    recount_tpm <- SummarizedExperiment::assays(rse)$TPM %>%
      as_tibble(rownames = "gene") %>% 
      mutate(gene = gene %>% str_remove("\\..*"))
    
    rm(rse)
    gc()
    
    ## 1. For each tissue within the current project, filter the RSE by its samples
    ## 2. Calculate median TPM value of each gene across samples of the current tissue
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          project_id, "/v", gtf_version, "/", main_project, "_project/")
    
    if (file.exists(paste0(folder_root, "/raw_data/samples_metadata.rds"))) {
      
      metadata.info <- readRDS(file = paste0(folder_root, "/raw_data/samples_metadata.rds"))
      
      clusters_ID <- metadata.info$gtex.smtsd %>% unique()
      
      for (cluster in clusters_ID) {
        
        # cluster <- clusters_ID[1]
        
        
        
        if ( main_project == "introverse" ) {
          samples <- metadata.info %>% 
            as_tibble() %>%
            filter(gtex.smtsd == cluster,
                   #gtex.smrin >= 6.0,
                   gtex.smafrze != "EXCLUDE") %>%
            distinct(external_id) %>% 
            pull()
        } else {
          samples <- metadata.info %>% 
            as_tibble() %>%
            filter(gtex.smtsd == cluster,
                   gtex.smrin >= 6.0,
                   gtex.smafrze != "EXCLUDE") %>%
            distinct(external_id) %>% 
            pull()
        }
        
        
        
        
        if (((main_project == "splicing" ||
              str_detect(main_project, pattern = "age")) && length(samples) >= 70) || 
            (main_project == "introverse" && length(samples) >= 1)) {
          
          cluster_samples <- readRDS(file = paste0(folder_root, "/raw_data/", 
                                                   project_id, "_", cluster, "_samples_used.rds"))
          
          if (!identical(samples, cluster_samples)) {
            print("ERROR - samples are not identical")  
          }
          
          n_cols <- length(cluster_samples) # length(names(recount_tpm))-1
          
          recount_tpm_local <- recount_tpm %>%
            dplyr::select(c("gene", all_of(cluster_samples))) 
          
          # recount_tpm_local$TPM_median %>% unique()
          
          folder_name <- paste0(folder_root, "/results/tpm/")
          dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)
          saveRDS(object = recount_tpm_local,
                  file = paste0(folder_name, 
                                project_id, "_", cluster, "_tpm.rds"))
          
          
          rm(recount_tpm_local)
          rm(cluster_samples)
          gc()
          
        }
        
      }
    }
    
    print(paste0(Sys.time(), " - ", project_id, " finished!"))
    
    rm(folder_root)
    rm(clusters_ID)
    rm(recount_tpm)
    gc()
    
  }
  
}

tidy_gtex_tpm <- function() {
  
  ## LOAD THE TPM DATA
  tpm <- aws.s3::s3read_using(FUN = CePa::read.gct,
                              object = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                              bucket = "data-references") %>%
    as.data.frame()
  
  colnames(tpm) <- all_clusters_used[-24]
  tpm <- tpm %>%
    tibble::rownames_to_column("gene_id") %>%
    mutate(gene_id = gene_id %>% str_remove("\\..*"))
  
  tpm %>% head()
  
  ## Save
  saveRDS(object = tpm,
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/tpm_tidy.rds")
  
  
}

generate_biotype_percentage <- function() {
  
  
  #######################################
  ## GET THE TRANSCRIPT BIOTYPE
  #######################################
  
  
  print(paste0(Sys.time(), " - loading the human reference transcriptome ... "))
  
  ## Import HUMAN REFERENCE transcriptome
  homo_sapiens_v105 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf") %>% 
    as_tibble()
  
  ## Get v105 transcripts
  transcripts_v105 <- homo_sapiens_v105 %>%
    filter(type == "transcript") %>% 
    dplyr::select(transcript_id, transcript_biotype, gene_id)
  
  transcripts_v105 %>% head()
  
  
  ## LOAD ALL ANNOTATED RECOUNT3 JUNCTIONS
  ## We load all projects, one by one to load the SR data 
  
  print(paste0(Sys.time(), " - loading the recount3 split reads ... "))
  
  all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
  
  #######################################
  ## LOOP THROUGH THE TISSUES
  #######################################
  
  df_all <- map_df(all_projects, function(project) {
    
    # project <- all_projects[7]
    print(project)
    
    all_clusters <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                          project, "/raw_data/all_clusters.rds"))
    
    df_all <- map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[1]
      print(cluster)
      df_all <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                      project, "/results/base_data/", cluster, "/", 
                                      cluster, "_annotated_SR_details_length_105.rds")) %>%
        dplyr::select(junID, seqnames, start, end, strand, tx_id_junction)
      
      return(df_all)
      
    })
    
    return(df_all)
    
  })
  
  ## Check if the junction "chr1:100007157-100011364:+" belongs to multiple genes across the tissues
  df_test_1 <- merge(x = df_all %>% filter(junID == "chr1:100007157-100011364:+") %>% unnest(tx_id_junction),
                     y = transcripts_v105,
                     by.x = "tx_id_junction",
                     by.y = "transcript_id",
                     all.x = T)
  df_test_1 %>%
    group_by(junID) %>% 
    distinct(gene_id, .keep_all = T) %>% 
    dplyr::count() %>% 
    distinct(n, .keep_all = T)
  
  ## Check if the junction "chr7:117604990-117616228:-", which ultimately has
  ## been classified as lncRNA = 100, has only lncRNA transcripts
  df_biotype %>%
    filter(junID == "chr7:117604990-117616228:-")
  
  saveRDS(object = df_all %>%
            distinct(junID, .keep_all = T) %>% 
            mutate(ref_coord = paste0("chr", seqnames, ":", 
                                      start, "-", end, ":", strand)) %>% 
            dplyr::select(-junID) %>%
            dplyr::rename(junID = ref_coord) %>% 
            data.table::as.data.table(),
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_introns_raw.rds")
  
  #######################################
  ## EXPLORE AND TIDY THE RESULT
  #######################################
  
  df_all_introns <- df_all %>%
    mutate(ref_coord = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>% 
    dplyr::select(-junID) %>%
    dplyr::rename(junID = ref_coord)
  
  # indx <- which(str_detect(df_all_introns %>% distinct(junID) %>% pull(), pattern = "\\*"))
  # if (indx %>% length > 0) {
  #   print("ERROR!")
  # }
  
  #df_all_introns %>% head()
  #df_all_introns %>% nrow()
  #print(object.size(df_all_introns), units = "Gb")
  
  
  ## Merge datasets to add transcript biotype
  print(paste0(Sys.time(), " --> adding transcript biotype..."))
  
  df_all_introns <- df_all_introns %>% data.table::as.data.table()
  transcripts_v105 <- transcripts_v105 %>% data.table::as.data.table()
  
  df_all_introns <- merge(x = df_all_introns,
                          y = transcripts_v105,
                          by.x = "tx_id_junction",
                          by.y = "transcript_id",
                          all.x = T)
  
  # df_all_introns %>% head()
  # df_all_introns %>% nrow()
  
  #print(object.size(df_all_introns), units = "Gb")
  
  
  
  #######################################
  ## CALCULATE THE TRANSCRIPT PERCENTAGE
  #######################################
  
  
  print(paste0(Sys.time(), " --> starting protein-coding percentage calculation ..."))
  #print(paste0(Sys.time(), " --> ", df_all_introns$junID %>% unique() %>% length(), " total number of junctions."))
  
  
  
  ## Remove ambiguous jxn (belonging to multiple genes)
  junID_OK <- df_all_introns %>% 
    group_by(junID) %>% 
    distinct(gene_id, .keep_all = T) %>% 
    dplyr::count() %>% 
    filter(n == 1) %>%
    pull(junID)
  
  # junID %>% unique() %>% length()
  
  ## Calculate the biotype percentage
  df_all_percentage <- df_all_introns %>% 
    filter(junID %in% junID_OK) %>%
    group_by(junID, transcript_biotype) %>%
    distinct(tx_id_junction, .keep_all = T) %>% 
    summarise(n = n()) %>% 
    mutate(percent = (n / sum(n)) * 100) %>%
    ungroup()
  
  saveRDS(object = df_all_percentage,
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_annotated_SR_details_length_105_raw_biotype.rds")
  
  ## Only filter by the protein-coding biotype
  df_all_percentage_tidy_PC <- df_all_percentage %>% 
    group_by(junID) %>%
    rowwise() %>%
    mutate(percent = ifelse (transcript_biotype == "protein_coding", percent, 0)) %>%
    ungroup() %>%
    group_by(junID) %>%
    filter(percent == max(percent)) %>%
    dplyr::select(-transcript_biotype, -n) %>%
    distinct(junID, .keep_all = T) %>%
    ungroup()
  
  ## Only filter by the lncRNA biotype
  df_all_percentage_tidy_lncRNA <- df_all_percentage %>% 
    group_by(junID) %>%
    rowwise() %>%
    mutate(percent = ifelse (transcript_biotype == "lncRNA", percent, 0)) %>%
    ungroup() %>%
    group_by(junID) %>%
    filter(percent == max(percent)) %>%
    dplyr::select(-transcript_biotype, -n) %>%
    distinct(junID, .keep_all = T) %>%
    ungroup()
  
  
  df_all_percentage_tidy_merged <- merge(x = df_all_percentage_tidy_PC %>% dplyr::rename(percent_PC = percent) %>% as.data.table(),
                                         y = df_all_percentage_tidy_lncRNA %>% dplyr::rename(percent_lncRNA = percent) %>% as.data.table(),
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
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_annotated_SR_details_length_105_biotype.rds")
  
  
  ##########################################
  ## FREE UP SOME MEMORY 
  ##########################################
  
  rm(df_all)
  rm(df_all_introns)
  rm(homo_sapiens_v105)
  rm(transcripts_v105)
  print(paste0(Sys.time(), " - file saved!"))
}


#################################################################
## SQL HELPER
## auxiliary functions to help in the generation of a database
#################################################################

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

# all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects.rds")
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
    
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                          project_id, "/v", gtf_version, "/", 
                          main_project, "_project/")
    
    if (is.null(all_clusters)) {
      
      metadata.info <- readRDS(file = paste0(folder_root, "/raw_data/samples_metadata.rds"))
      
      if ( main_project == "introverse" ) {
        all_clusters <-  metadata.info %>% 
          as_tibble() %>%
          filter(#gtex.smrin >= 6.0,
            gtex.smafrze != "EXCLUDE") %>%
          distinct(gtex.smtsd) %>% 
          pull()
      } else {
        all_clusters <-  metadata.info %>% 
          as_tibble() %>%
          filter(gtex.smrin >= 6.0,
                 gtex.smafrze != "EXCLUDE") %>%
          distinct(gtex.smtsd) %>% 
          pull()
      }
      
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
  df_all_distances_pairings_raw <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    # project_id <- "KIDNEY"
    
    print(paste0(Sys.time(), " --> Working with '", project_id, "' DataBase..."))
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                          project_id, "/v", gtf_version, "/", main_project, "_project/")
    
    if (is.null(all_clusters)) {
      
      metadata.info <- readRDS(file = paste0(folder_root, "/raw_data/samples_metadata.rds"))
      
      
      if ( main_project == "introverse" ) {
        all_clusters <-  metadata.info %>% 
          as_tibble() %>%
          filter(#gtex.smrin >= 6.0,
            gtex.smafrze != "EXCLUDE") %>%
          distinct(gtex.smtsd) %>% 
          pull()
      } else {
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
      
      print(paste0(Sys.time(), " - ", project_id, " loading '", cluster, "'  data ..."))
      
      ## Load samples
      if ( file.exists( paste0(folder_root, "/raw_data/", project_id, "_", cluster,  "_samples_used.rds") ) ) {
        
        samples <- readRDS(file = paste0(folder_root, "/raw_data/", project_id, "_", cluster,  "_samples_used.rds"))
        
        if (samples %>% length() > 0) {
          
          folder_name <- paste0(folder_root, "/results/", cluster, "/")
          
          if ( !file.exists(paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) ) {
            
            ## Obtain the distances across all samples
            df_all <- map_df(samples, function(sample) { 
              
              # sample <- samples[1]
              print(paste0(cluster, " - ", sample))
              file_name <- paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")
              
              
              if (file.exists(file_name)) {
                
                df <- readRDS(file = file_name)
                
                return(df)
              } else {
                return(NULL)
              }
              
            })
            
            if ( !is.null(df_all) ) {
              saveRDS(object = df_all %>%
                        distinct(novel_junID, ref_junID, .keep_all = T) %>%
                        mutate(tissue = cluster),
                      file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
            } 
          } else {
            df_all <- readRDS( file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds") )
          }
          
          
          # df_all2 <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
          return(df_all %>%
                   distinct(novel_junID, ref_junID, .keep_all = T))
        } else {
          return(NULL)
        }
        
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
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                          db, "/v", gtf_version, "/", main_project, "_project/")
    
    if ( is.null(all_clusters) ) {
      
      if (file.exists(paste0(base_folder, "/raw_data/samples_metadata.rds"))) {
        metadata.info <- readRDS(file = paste0(base_folder, "/raw_data/samples_metadata.rds"))
        
        if ( main_project == "introverse" ) {
          all_clusters <-  metadata.info %>% 
            as_tibble() %>%
            filter(#gtex.smrin >= 6.0,
              gtex.smafrze != "EXCLUDE") %>%
            distinct(gtex.smtsd) %>% 
            pull()
        } else {
          all_clusters <-  metadata.info %>% 
            as_tibble() %>%
            filter(gtex.smrin >= 6.0,
                   gtex.smafrze != "EXCLUDE") %>%
            distinct(gtex.smtsd) %>% 
            pull()
        }
        
      }
      
    }
    
    map_df(all_clusters, function(cluster) { 
      
      # cluster <- all_clusters[1]
      
      # print(paste0(Sys.time(), " --> ", cluster))
      if ( file.exists(paste0(base_folder, "results/", cluster, 
                              "/not-misspliced/", cluster, "_all_notmisspliced.rds")) ) {
        df_introns_never <- readRDS(file = paste0(base_folder, "results/", cluster, 
                                                  "/not-misspliced/", cluster, "_all_notmisspliced.rds")) %>% as_tibble()
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
    
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", project_id, "/")
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





