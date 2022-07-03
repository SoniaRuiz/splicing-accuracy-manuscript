library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/init.R")

.libPaths( c( "/home/sruiz/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()) )
Sys.setenv("AWS_ACCESS_KEY_ID"="AKIAUKEFAIIIH7PS3Z7V",
           "AWS_SECRET_ACCESS_KEY"="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",
           "AWS_DEFAULT_REGION"="eu-west-2")


get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}



##################################
## CALLS - PREPARE RECOUNT3 DATA
##################################

source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline1_download_from_recount.R")

download_from_recount <- function(project_id,
                                  recount_version = 3,
                                  clusters_ID,
                                  rse = NULL,
                                  folder_root = NULL) {
  
  print(paste0(Sys.time(), " - ", project_id, " --> preparing data from RSE object."))
  
  if (is.null(folder_root)) {
    folder_root <- getwd()
    folder_root <- paste0(folder_root, "/", project_id, "/")
    # folder_root <- "PROJECTS/splicing-project/splicing-recount3-projects/"
  }
  
  dir.create(file.path(folder_root), recursive = T, showWarnings = T)
  
  folder_origin <- paste0(folder_root, "/raw_data/")
  dir.create(file.path(folder_origin), recursive = T, showWarnings = T)
  
  folder_destiny <- paste0(folder_root, "/results/base_data/")
  dir.create(file.path(folder_destiny), recursive = T, showWarnings = T)
  
  ## prepare_data_from_rse -----------------------

  if (is.null(rse)) {
    if (recount_version == 2) {
      download.file(url = URLencode(URL = paste0("http://duffel.rail.bio/recount/", project_id, "/rse_jx.Rdata")),
                    destfile = paste0(folder_origin, "/rse_jx.Rdata"), quiet = F)
  
      load(paste0(folder_origin, "/rse_jx.Rdata"))
    } else {
      rse <- recount3::create_rse_manual(
        project = project_id,
        project_home = "data_sources/gtex",
        organism = "human",
        annotation = "gencode_v29",
        type = "jxn"
      )
    }
  }

  prepare_data_from_rse(rse = rse,
                        folder_path = folder_origin)

  ## generate_all_split_reads ---------------------

  print(paste0(Sys.time(), " - ", project_id, " --> generating split reads file."))
  
  generate_all_split_reads(folder_path = folder_origin,
                           folder_destiny = folder_destiny)


  ## separate_samples -------------------------------

  # if (project_id == "SRP051844") {
  #   clusters_ID <- c("C_","H_")#(rse_jx %>% SummarizedExperiment::colData())$sra.sample_title
  # } else if (project_id =="SRP058181") {
  #   clusters_ID <- c("C_","P_")#(rse_jx %>% SummarizedExperiment::colData())$sra.sample_title
  # } else if (project_id == "ADRENAL_GLAND") {
  #   clusters_ID <- readRDS(file = "PROJECTS/splicing-project/splicing-recount3-projects/ADRENAL_GLAND/raw_data/metadata.rds")$gtex.smtsd %>%
  #     unique()
  # }
  
  print(paste0(Sys.time(), " - ", project_id, " --> separating samples into clusters."))
  
  saveRDS(object = clusters_ID, file = paste0(folder_origin, "/all_clusters_used.rds"))
  
  separate_samples(clusters_ID = clusters_ID,
                   folder_origin = folder_origin,
                   project_id = project_id,
                   folder_destiny = folder_destiny,
                   GTEx = T)

  
  ## separate_split_read_counts_cluster --------------------
  
  print(paste0(Sys.time(), " - ", project_id, " --> generating split reads counts file across clusters."))
  
  # clusters_ID <- readRDS(file = paste0(folder_origin, "/all_clusters_used.rds"))
  separate_split_read_counts_cluster(clusters_ID = clusters_ID,
                                     folder_origin = folder_origin,
                                     project_id = project_id,
                                     folder_destiny = folder_destiny,
                                     rse = rse,
                                     GTEx = T)
  
}


############################
## CALLS - ANNOTATION
############################

source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline2_annotate_from_recount.R")


annotate_from_recount <- function(project_id,
                                  gtf_version = NULL,
                                  clusters_ID = NULL,
                                  folder_root = NULL) {
  
  if (is.null(folder_root)) {
    folder_root <- getwd()
    folder_root <- paste0(folder_root, "/", project_id, "/")
    dir.create(file.path(folder_root), recursive = T, showWarnings = T)
    # folder_root <- paste0("PROJECTS/splicing-project/splicing-recount2-projects/", project_id, "/")
  }
  
  if (is.null(clusters_ID)) {
    clusters_ID <- readRDS(file = paste0(folder_root, "/raw_data/all_clusters_used.rds"))
  }
  
  if (is.null(gtf_version)) {
    gtf_version <- 105
  }
    
  if (file.exists(paste0("/data/references/ensembl/gtf/v", gtf_version, "/Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf"))) {
    gtf_path <- paste0("/data/references/ensembl/gtf/v", gtf_version, "/Homo_sapiens.GRCh38.", gtf_version, ".chr.gtf")
  } else {
    gtf_path <- paste0("/data/references/ensembl/gtf/v", gtf_version, "/Homo_sapiens.GRCh38.", gtf_version, ".gtf")
  }
    
  edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  
  ## Load all split reads
  folder_path <- paste0(folder_root, "results/base_data/")
  all_split_reads <- readRDS(file = paste0(folder_path, "/all_split_reads.rds"))
  
  
  ## Annotate split reads per cluster
  for (cluster in clusters_ID) {
    
    print(paste0(Sys.time(), " - Annotating '", cluster, "' samples..."))
    
    # cluster <- clusters_ID[2]
    folder_source <- paste0(folder_path, "/", cluster)
    
    ## Load split read counts
    split_read_counts <- readRDS(file = paste0(folder_source, "/", project_id, "_", cluster, "_split_read_counts.rds"))
    
    ## Load samples
    cluster_samples <- readRDS(file = paste0(folder_source, "/", project_id, "_", cluster, "_samples.rds"))
    
    
    if (cluster_samples %>% length() > 0) {
      
      get_base_data_annotated(cluster = cluster,
                              samples = cluster_samples,
                              all_split_reads = all_split_reads,
                              split_read_counts = split_read_counts,
                              blacklist_path = "/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed",
                              edb = edb,
                              gtf_version = gtf_version,
                              folder_save_path = folder_source)
      
      
      all_split_reads_details_105 <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/",
                                                           cluster, "_annotated_SR_details_length_", gtf_version, ".rds"))
  
      split_read_counts <- split_read_counts %>%
        filter(junID %in% all_split_reads_details_105$junID)
      
      split_read_counts <- split_read_counts[rowSums(split_read_counts[, c(2:ncol(split_read_counts))])>0,]
      
      saveRDS(object = split_read_counts, file = paste0(folder_source, "/", 
                                                        project_id, "_", cluster, "_split_read_counts_", gtf_version, ".rds"))
      
      rm(all_split_reads_details_105)
      rm(split_read_counts)
      rm(cluster_samples)
      rm(folder_source)
      
    }
    
    gc()
  }

  
  rm(all_split_reads)
  rm(edb)
  gc()
}

############################
## IDB GENERATION
############################

source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_generation.R")

idb_generation <- function(project_id,
                           gtf_version,
                           clusters_ID = NULL,
                           folder_root = NULL) {
  
  if (is.null(folder_root)) {
    folder_root <- getwd()
    folder_root <- paste0(folder_root, "/", project_id, "/")
    dir.create(file.path(folder_path), recursive = T, showWarnings = T)
  }
  
  if (is.null(clusters_ID)) {
    clusters_ID <- readRDS(file = paste0(folder_root, "/raw_data/all_clusters_used.rds"))
  }
  
  for (cluster in clusters_ID) {
    
    # cluster <- clusters_ID[1]
    # cluster <- clusters_ID[2]
    
    print(paste0(Sys.time(), " - loading '", cluster, "' source data ..."))
    
    ############################################
    ## LOAD DATA FOR THE CURRENT PROJECT
    ############################################
    
    ## Load samples
    samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster,  "_samples.rds"))
    
    if (samples %>% length() > 0) {

      folder_name <- paste0(folder_root, "/results/pipeline3/distances/", cluster, "/v", gtf_version, "/")
      dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)


      ## Load split read data
      all_split_reads_details <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/",
                                                           cluster, "_annotated_SR_details_length_", gtf_version, ".rds"))
      ## Load split read counts
      split_read_counts <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/",
                                                 project_id, "_", cluster, "_split_read_counts_", gtf_version, ".rds"))


      ############################################
      ## DISTANCES SUITE OF FUNCTIONS
      ############################################

      get_distances(cluster = cluster,
                    samples = samples,
                    split_read_counts = split_read_counts,
                    all_split_reads_details = all_split_reads_details,
                    folder_name)


      extract_distances(cluster = cluster,
                        samples,
                        split_read_counts = split_read_counts,
                        folder_name = folder_name)


      get_never_misspliced(cluster = cluster,
                           samples = samples,
                           split_read_counts = split_read_counts,
                           all_split_reads_details = all_split_reads_details,
                           folder_name = folder_name,
                           save_results = T)


      extract_never_misspliced(cluster = cluster,
                               samples = samples,
                               split_read_counts = split_read_counts,
                               folder_name = folder_name)

      add_never_misspliced_to_df(cluster = cluster,
                                 samples = samples,
                                 all_split_reads_details = all_split_reads_details,
                                 folder_name = folder_name)
  
      ############################################
      ## MIS-SPLICING RATIO SUITE OF FUNCTIONS
      ############################################
  
      folder_idb_name <-paste0(folder_root, "/results/pipeline3/missplicing-ratio/", cluster, "/v", gtf_version, "/")
      dir.create(file.path(folder_idb_name), recursive = T, showWarnings = F)
  
      get_missplicing_ratio(cluster = cluster,
                            samples = samples,
                            split_read_counts = split_read_counts,
                            folder_name = folder_name,
                            folder_save_name = folder_idb_name)

      add_missplicing_class_to_df(cluster = cluster,
                                  folder_name = folder_idb_name)

      get_missplicing_QC(cluster = cluster,
                         samples = samples,
                         split_read_counts = split_read_counts,
                         all_split_reads_details = all_split_reads_details,
                         folder_name = folder_idb_name)
  
  
      ############################################
      ## ADD FEATURES TO THE IDB
      ############################################
  
      remove_MT_genes(cluster = cluster,
                      folder_name = folder_idb_name)

      add_intron_type(cluster = cluster,
                      folder_name = folder_idb_name)

      clinvar_analysis(cluster = cluster,
                       folder_name = folder_idb_name)

      add_MANE_info(cluster = cluster,
                    folder_name = folder_idb_name)
      
      
      
      #################################
      ## CONFIGURE TO ADD THESE PROPERTIES TO THE IDB
      ## Note: is the Shiny app the portal in which I
      ## restrict the info that I show there.
      #################################
      
      add_gene_tpm_length(cluster = cluster,
                          tpm_file = paste0(folder_root, "/results/base_data/", cluster, 
                                            "/", project_id, "_", cluster, "_tpm.rds"),
                          folder_name = folder_idb_name)
      
      add_cdts_cons_scores(cluster = cluster,
                           folder_name = folder_idb_name)

      add_protein_percentage_to_database(cluster = cluster,
                                         folder_root = folder_idb_name)
    
    }
    
  }
}



############################
## CALLS 
############################

## Load the dependency needed for the conservation/constraint scores        
# print("Loading the 'CNC_CDTS_CONS_gr' file...")
# aws.s3::s3load(object = "CNC_CDTS_CONS_gr.rda", bucket = "data-references", region = "eu-west-2")
# print("'CNC_CDTS_CONS_gr' file loaded!")

## Loop through each project
all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
gtf_version <- 105

for (project_id in all_projects) {

  # project_id <- all_projects[31]
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
  
  print(paste0(Sys.time(), " - getting data from '", project_id, "' tissue..."))

  rse <- recount3::create_rse_manual(
    project = project_id,
    project_home = "data_sources/gtex",
    organism = "human",
    annotation = "gencode_v29",
    type = "jxn"
  )

  download_from_recount(project_id = project_id,
                        recount_version = 3,
                        clusters_ID = rse$gtex.smtsd %>% unique(),
                        rse = rse,
                        folder_root = folder_root)

  rm(rse)
  gc()
   
  clusters_ID <- readRDS(file = paste0(folder_root, "/raw_data/all_clusters_used.rds"))
   
  annotate_from_recount(project_id = project_id,
                        gtf_version = gtf_version,
                        clusters_ID = clusters_ID,
                        folder_root = folder_root)
  gc()

  # idb_generation(project_id = project_id,
  #                gtf_version = gtf_version,
  #                clusters_ID = clusters_ID,
  #                folder_root = folder_root)
  # gc()
  

  
}






#############################################
## ADDING THE TPM VALUE GENERATED BY RECOUNT3
#############################################

# for (project_id in all_projects) {
#   
#   folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
#   clusters_ID <- readRDS(file = paste0(folder_root, "/raw_data/all_clusters_used.rds"))
#   
#   for (cluster in clusters_ID) {
#   
#     # cluster <- clusters_ID[1]
#   
#     tpm_file <- paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster, "_tpm.rds")
#     folder_name <- paste0(folder_root, "/results/pipeline3/missplicing-ratio/", cluster, "/v", gtf_version, "/")
#   
#     print(cluster)
#   
#     ## LOAD SOURCE DATA
#     file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
#   
#     if (file.exists(file_name) && file.exists(tpm_file)) {
#   
#       df_introns <- readRDS(file = file_name)
#   
#       ## GENE TPM
#       print(paste0(Sys.time(), " - Adding gene TPM..."))
#   
#       tpm <- readRDS(file = tpm_file)
#   
#       tpm %>% head()
#       tpm %>% nrow()
#   
#       tpm <- tpm  %>%
#         dplyr::select(gene_id = gene,
#                       tpm_median = TPM_median,
#                       tpm_mean = TPM_mean)
#   
#       tpm %>% head()
#       tpm %>% nrow()
#   
#   
#       df_introns <- df_introns %>%
#         unnest(gene_id)
#   
#   
#       df_merged <- merge(x = df_introns %>% data.table::as.data.table(),
#                          y = tpm %>% data.table::as.data.table(),
#                          by = "gene_id",
#                          all.x = T) %>%
#         dplyr::rename(tpm_mean_rct3 = tpm_mean)%>%
#         dplyr::rename(tpm_median_rct3 = tpm_median)
#   
#   
#       df_merged %>% head()
#       df_merged %>% nrow()
#   
#       saveRDS(object = df_merged, file = file_name)
#     }
#   
#   }
# }



#############################################
## ADDING ONLY THE CONSERVATION SCORES
#############################################

# for (project_id in all_projects[-c(1:15)]) {
# 
#   
#   folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
#   clusters_ID <- readRDS(file = paste0(folder_root, "/raw_data/all_clusters_used.rds"))
# 
#   
#   for (cluster in clusters_ID) {
# 
#     # cluster <- clusters_ID[1]
#     print(paste0(Sys.time(), " - ", cluster))
# 
#     folder_name <- paste0(folder_root, "/results/pipeline3/missplicing-ratio/", cluster, "/v", gtf_version, "/")
# 
#     add_cdts_cons_scores(cluster = cluster,
#                          folder_name = folder_name)
# 
# 
#   }
#   
# }
