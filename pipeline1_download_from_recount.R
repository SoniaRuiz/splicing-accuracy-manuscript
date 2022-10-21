################################
## FUNCTIONS
################################

library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline1_download_from_recount.R")

prepare_data_from_rse <- function(rse,
                                  folder_path) {

  
  ## Get samples and store
  samples <- rse %>% 
    SummarizedExperiment::colData() %>%
    as_tibble() %>%
    filter(gtex.smafrze != "EXCLUDE") %>%
    pull(external_id)
  saveRDS(object = samples, file = paste0(folder_path, "/samples.rds"))
  
  ## Save samples metadata info
  metadata.info <- rse %>% 
    SummarizedExperiment::colData() %>%
    as_tibble() %>%
    filter(gtex.smrin >= 6.0,
           gtex.smafrze != "EXCLUDE")
  saveRDS(object = metadata.info, file = paste0(folder_path, "/samples_metadata.rds"))
  
  ## Save junctions metadata info
  junctions.info <- rse %>% 
    SummarizedExperiment::elementMetadata()
  saveRDS(object = junctions.info, file = paste0(folder_path, "/junctions_metadata.rds"))
  
  ## Save all_split_reads data
  all_split_reads <- data.frame(chr = rse %>% SummarizedExperiment::seqnames(),
                                start = rse %>% SummarizedExperiment::start(),
                                end = rse %>% SummarizedExperiment::end(),
                                strand = rse %>% SummarizedExperiment::strand(),
                                width = rse %>% SummarizedExperiment::width(),
                                junID = rse %>% rownames())
  saveRDS(object = all_split_reads, file = paste0(folder_path, "/all_split_reads.rds"))
  # all_split_reads <- readRDS(file = paste0(folder_path, "/all_split_reads.rds")) 

  
  
  
  rm(samples)
  rm(metadata.info)
  rm(junctions.info)
  rm(all_split_reads)
  #rm(counts)
  gc()
}

generate_all_split_reads <- function(folder_path,
                                     folder_destiny) {
  
  ## Junction ID BED
  all_split_reads <- readRDS(file = paste0(folder_path, "/all_split_reads.rds"))
  all_split_reads %>% head()
  
  
  all_split_reads_tidy <- all_split_reads %>%
    GenomicRanges::GRanges() %>%
    diffloop::rmchr()
  
  
  all_split_reads_tidy %>% head()
  all_split_reads_tidy <- GenomeInfoDb::keepSeqlevels(x = all_split_reads_tidy,
                                                      value = intersect(all_split_reads_tidy %>% 
                                                                          seqnames %>% 
                                                                          levels(), 
                                                                        c( "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                                           "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X", "Y")), 
                                                      pruning.mode = "tidy")
  all_split_reads_tidy %>% head()
  
  ## QC
  if (any((GenomeInfoDb::seqlevels(all_split_reads_tidy) %>% unique()) == "MT")) {
    print("Error: some junctions from the MT chromosome have not been removed!")
  }
  
  # test <- dasper::junction_annot(all_split_reads_tidy %>% GRanges(), 
  #                                ref = edb,
  #                                ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"), 
  #                                ref_cols_to_merge = c("gene_id", "gene_name", "tx_id"))
  # test %>%
  #   as_tibble() %>%
  #   dplyr::count(type)
  
  # ranges(all_split_reads_tidy) <- IRanges::IRanges(start = start(all_split_reads_tidy) + 1, 
  #                                                  end = end(all_split_reads_tidy) + 1)
  
  
    
  all_split_reads_tidy %>% head()
  
  if (any(all_split_reads_tidy$junID %>% duplicated())) {
    print(paste0("Split reads duplicated: ", any(all_split_reads_tidy$junID %>% duplicated())))
  }
  
  
  
  ##########################
  ## SAVE RESULTS
  ##########################
  
  
  saveRDS(object = all_split_reads_tidy, file = paste0(folder_destiny, "/all_split_reads.rds"))
  
  
  rm(all_split_reads)
  rm(all_split_reads_tidy)
  gc()
}

separate_samples <- function(clusters_ID,
                             project_id,
                             folder_origin,
                             folder_destiny,
                             GTEx = T) {
  
  ## Load recount2 project metadata
  metadata.info <- readRDS(file = paste0(folder_origin, "/samples_metadata.rds"))
  metadata.info %>% head()
  
  for (cluster in clusters_ID) { 
    # cluster <- clusters_ID[1]
    
    print(paste0(Sys.time(), " - getting samples from '", cluster, "' cluster..."))
    
    if (GTEx) {
      
      cluster_samples <- metadata.info %>% 
        as_tibble() %>%
        filter(gtex.smtsd == cluster,
               gtex.smrin >= 6.0,
               gtex.smafrze != "EXCLUDE") %>%
        distinct(external_id) %>% 
        pull()
      
      print(paste0(Sys.time(), " - ", cluster, ": ", cluster_samples %>% length(), " samples!"))
      
    } else {
      ## Get the control samples and save them
      indx <- which(str_detect(metadata.info$sra.sample_title, cluster) == T) 
      
      if (indx %>% length() == 0) {
        indx <- which(str_detect(metadata.info$sample_attributes, cluster) == T) 
      }
      
      cluster_samples <- metadata.info[indx, ] %>% rownames() %>% as.character()
      
      if (any(cluster_samples == "SRR2015746")) {
        ind <- which(cluster_samples == "SRR2015746")
        cluster_samples <- cluster_samples[-ind]
      }
      rm(indx)
      
    }
    
    folder_save <- paste0(folder_destiny, "/", cluster)
    dir.create(file.path(folder_save), recursive = T, showWarnings = T)
    
    saveRDS(object = cluster_samples,
            file = paste0(folder_save, "/", project_id, "_", cluster, "_samples.rds"))
    
    
    
    rm(cluster_samples)
    rm(folder_save)
    gc()
  }
  
  rm(metadata.info)
  gc()
  
}

separate_split_read_counts_cluster <- function(clusters_ID,
                                               project_id,
                                               folder_origin,
                                               folder_destiny,
                                               rse = NULL,
                                               GTEx = T) {
  
  all_split_reads <- readRDS(file = paste0(folder_destiny, "/all_split_reads.rds"))$junID
  gc()
  
  for (cluster in clusters_ID) { 
    # cluster <- clusters_ID[1]
    cluster_samples <- readRDS(file = paste0(folder_destiny, "/", cluster, "/", project_id, "_", cluster, "_samples.rds"))
    
    if (cluster_samples %>% length() > 0) {
      
      print(paste0(Sys.time(), " - getting split_read_counts from ", cluster_samples %>% length(), 
                   " '", cluster, "' samples..."))
      if ("raw_counts" %in% assayNames(rse)) {
        SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
      }
      counts <- (rse %>% SummarizedExperiment::assays())[[1]]
      counts <- counts[, cluster_samples %>% as.character(), drop = FALSE]
      counts <- counts[rownames(counts) %in% all_split_reads, , drop = FALSE]
     
      counts <- counts[rowSums(counts) > 0, ]
     
      counts <- counts %>% as.matrix()
      junID <- counts %>% rownames()
      
      counts <- counts %>%
        as_tibble() %>%
        dplyr::mutate("junID" = junID) %>%
        dplyr::relocate(junID) %>%
        filter(junID %in% all_split_reads)
      
      
      print(paste0(Sys.time(), " - ", counts %>% nrow(), " split reads."))
      # counts[rowSums(counts %>% select(-junID))>0,]
      # 
      # if (GTEx) {
      # 
      #   split_read_counts_cluster <- split_read_counts[,cluster_samples %>% as.character(), drop = FALSE]
      #   gc()
      #   
      #   split_read_counts_cluster <- as.matrix(split_read_counts_cluster)
      #   gc()
      #   
      #   split_read_counts_cluster <- split_read_counts_cluster %>%
      #     as_tibble() %>% 
      #     mutate("junID" = rse %>% rownames()) %>%
      #     relocate("junID")
      #   gc()
      #   
      #  
      #   #split_read_counts_cluster[split_read_counts_cluster == 0] <- NA
      #   
      # } else {
      #   #split_read_counts %>% head()
      #   split_read_counts_cluster <- split_read_counts %>%
      #     select(junID, cluster_samples %>% as.character())
      #   
      # }
      # 
      # split_read_counts_cluster <- split_read_counts_cluster[rowSums(split_read_counts_cluster[, c(2:ncol(split_read_counts_cluster))])>0,]
      gc()
      
      folder_save <- paste0(folder_destiny, "/", cluster)
      dir.create(file.path(folder_save), recursive = T, showWarnings = T)
      
      print(paste0(Sys.time(), " - ", cluster, " saving split reads..."))
      
      saveRDS(object = counts,
              file = paste0(folder_save, "/", project_id, "_", cluster, "_split_read_counts.rds"))
      
      
      ## FREE UP SOME MEMORY
      rm(cluster_samples)
      rm(folder_save)
      rm(counts)
      
    }
    
    gc()
    
  }
  
  rm(all_split_reads)
  rm(counts)
  gc()
  
}





## Requires sam tools

### Requires bed tools
# wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
# tar -zxvf bedtools-2.30.0.tar.gz
# cd bedtools2
# make

### Requires perl

### Requires MaxEntScan score tool
# wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
# gunzip fordownload.tar.gz
# tar -xvf fordownload.tar

## FASTA file
# sudo wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --no-check-certificate
# sudo gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

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





############################################
## EXTRA FUNCTIONS NEEDED TO GENERATE 
## DEPENDENCIES
############################################

generate_recount2_median_tpm <- function() {
  
  source("/home/sruiz/secondary_projects/asap-gba-gbap1/R/file_paths.R")
  
  if (project_id == "SRP058181") {
    
    # SRP058181 data
    dds <- readRDS(
      file =
        file.path(
          path_to_lbd_seq,
          "results/SRP058181/gene_level_quant/SRP058181_DESeqDataSet.Rds"
        )
    )
    dds %>% colData %>% as.data.frame() %>% filter( Disease_Group == "Control") %>% select(Age_at_death) %>% pull(Age_at_death) %>% mean
    dds %>% colData %>% as.data.frame() %>% filter( Disease_Group == "PD") %>% select(Age_at_death) %>% pull(Age_at_death) %>% mean
    # Remove SRR2015746 (C0061), due to uncertainty re. sex
    dds <- dds[, dds$recount_id != "SRR2015746"]
    
  } else if (project_id == "SRP049203") {
    
    load(paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/",
                project_id, "/original_files/rse_fc_", project_id, ".Rdata"))
    
    rse_fc %>% colData %>% names
    rse_fc %>% colData %>% as.data.frame()
    
    library("DESeq2")
    ddsSE <- DESeqDataSet(rse_fc, design = ~ characteristics)
    ddsSE
    
    
  }
  
  # AIBS data
  # load(file="/home/rreynolds/data/scRNAseq_AIBS/MTG/AIBS2018_DataForEWCE.Rda")
  
  # Reference gtf
  ref <- rtracklayer::import(path_to_ref_gtf)
  ref <- ref %>% GenomeInfoDb::keepSeqlevels(c(1:22,"X","Y"), pruning.mode = "coarse")
  ref <- ref[ref$type == "gene"]
  
  
  
  
  # Remove anything after . in ensembl id
  rownames(dds) <-
    rownames(dds) %>%
    str_remove("\\..*")
  
  # Convert to tpm, which is calculated by:
  # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
  # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
  # 3. Divide the RPK values by the “per million” scaling factor.
  srp_rpk <-
    dds %>%
    assay() %>%
    as_tibble(rownames = "gene") %>%
    tidyr::pivot_longer(
      cols = -c("gene"),
      names_to = "recount_id",
      values_to = "counts"
    ) %>%
    dplyr::inner_join(
      ref %>%
        as_tibble() %>%
        dplyr::select(gene_id, gene_name, width),
      by = c("gene" = "gene_id")
    ) %>%
    dplyr::mutate(
      rpk = counts/width
    )
  srp_rpk[1,]
  srp_tpm <-
    srp_rpk %>%
    #dplyr::filter(gene %in% genes$gene_id) %>%
    dplyr::inner_join(
      srp_rpk %>%
        dplyr::group_by(recount_id) %>%
        dplyr::summarise(
          scaling_factor = sum(rpk)/1e6
        ),
      by = "recount_id"
    ) %>%
    dplyr::mutate(
      tpm = rpk/scaling_factor
    ) %>%
    dplyr::inner_join(
      colData(dds) %>%
        as_tibble(),
      by = "recount_id"
    )
  
  srp_tpm[1,]
  
  tpm <- srp_tpm %>%
    dplyr::group_by(gene) %>%
    mutate(tpm_median = tpm %>% median()) %>%
    dplyr::ungroup() %>%
    distinct(gene, .keep_all = T) %>%
    select(gene_id = gene, tpm, tpm_median)
  
  saveRDS(object = tpm,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/", project_id, "/tpm.rds"))
  
  
  
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


############################################
## CALLS
############################################

gtf_version <- 105
main_project <- "introverse"

if ( main_project == "introverse" ) {
  all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
} else {
  all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  #all_projects <- "NERVE"
}

generate_recount3_median_tpm(all_projects = all_projects[7],
                             gtf_version = gtf_version,
                             main_project = main_project)