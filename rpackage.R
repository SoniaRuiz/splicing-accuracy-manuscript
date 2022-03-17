library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)

# source("/home/sruiz/PROJECTS/splicing-project/rpackage.R")

.libPaths( c( "/home/sruiz/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()) )


# project_id <- "SRP058181"
# clusters_ID <- c("C_", "P_")


# project_id <- "SRP051844"
# clusters_ID <- c("C_", "H_")
# 
# 
# project_id <- "SRP035988"
# clusters_ID <- c("normal", "lesional")

# project_id <- "BRAIN"
# project_id <- "ADRENAL_GLAND"
# gtf_version <- 105
# folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")

get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}

##################################
## CALLS - PREPARE RECOUNT2 DATA
##################################

source("/home/sruiz/PROJECTS/splicing-project/pipeline0_prepare_data.R")

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

source("/home/sruiz/PROJECTS/splicing-project/pipeline1_QC_split_reads.R")


annotate_from_recount <- function(project_id = "BRAIN",
                                  gtf_path = NULL,
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
  
  
  if (is.null(gtf_path)) {
    gtf_path <- "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf"
    gtf_version <- 105
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
    
    saveRDS(object = split_read_counts, file = paste0(folder_source, "/", project_id, "_", cluster, "_split_read_counts.rds"))
    
    rm(all_split_reads_details_105)
    rm(split_read_counts)
    rm(cluster_samples)
    rm(folder_source)
    
    gc()
  }
  
  
  rm(all_split_reads)
  rm(edb)
  gc()
}

############################
## IDB GENERATION
############################

source("/home/sruiz/PROJECTS/splicing-project/pipeline3_methods.R")

idb_generation <- function(project_id,
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
    
    folder_name <- paste0(folder_root, "/results/pipeline3/distances/", cluster, "/v", gtf_version, "/")
    dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)
    
    
    ## Load split read data
    all_split_reads_details_105 <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/",
                                                         cluster, "_annotated_SR_details_length_", gtf_version, ".rds"))
    ## Load split read counts
    split_read_counts <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/",
                                               project_id, "_", cluster, "_split_read_counts.rds"))

    
    ############################################
    ## DISTANCES SUITE OF FUNCTIONS
    ############################################

    get_distances(cluster = cluster,
                  samples = samples[119:length(samples)],
                  split_read_counts = split_read_counts,
                  all_split_reads_details = all_split_reads_details_105,
                  folder_name)


    extract_distances(cluster = cluster,
                      samples,
                      split_read_counts = split_read_counts,
                      folder_name = folder_name)


    get_never_misspliced(cluster = cluster,
                         samples = samples,
                         split_read_counts = split_read_counts,
                         all_split_reads_details = all_split_reads_details_105,
                         folder_name = folder_name,
                         save_results = T)

     
    extract_never_misspliced(cluster = cluster,
                             samples = samples,
                             split_read_counts = split_read_counts,
                             folder_name = folder_name)

    add_never_misspliced_to_df(cluster = cluster,
                               samples = samples,
                               all_split_reads_details = all_split_reads_details_105,
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

    # get_missplicing_QC(cluster = cluster,
    #                    samples = samples,
    #                    split_read_counts = split_read_counts,
    #                    all_split_reads_details = all_split_reads_details_105,
    #                    folder_name = folder_idb_name)


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
    
    
    # add_cdts_cons_scores(cluster = cluster,
    #                      CNC_CDTS_CONS_gr = CNC_CDTS_CONS_gr,
    #                      folder_name = folder_idb_name)
    
    # add_gene_tpm_length(cluster = cluster,
    #                     folder_name = folder_idb_name,
    #                     GTEx = F)
    
    
    #################################
    ## ADD PROTEIN PERCENTAGE
    #################################
    
    
    # protein_folder <- "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_annotated_SR_details_length_104_biotype.rds"
    # folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/"
    # 
    # 
    # add_protein_percentage_to_database(cluster = cluster,
    #                                    protein_folder = protein_folder,
    #                                    folder_root = folder_root)
    
  }
}




############################
## CALLS 
############################

## DONE:
# project_id <- "BRAIN"
# project_id <- "ADRENAL_GLAND"
# project_id <- "KIDNEY"
# project_id <- "SMALL_INTESTINE"
# project_id <- "SALIVARY_GLAND"


## TODO:
# project_id <- "BLOOD"
projects <- c("SKIN",
              "ESOPHAGUS",
              "BLOOD_VESSEL",
              "ADIPOSE_TISSUE",
              "HEART",
              "MUSCLE",
              "COLON",
              "THYROID",
              "NERVE",
              "LUNG",
              "BREAST",
              "TESTIS",
              "STOMACH",
              "PANCREAS",
              "PITUITARY",
              "PROSTATE",
              "SPLEEN",
              "LIVER",
              "BONE_MARROW",
              "OVARY",
              "VAGINA",
              "UTERUS",
              "BLADDER",
              "CERVIX_UTERI",
              "FALLOPIAN_TUBE")


for (project_id in projects[1]) {
  
  project_id <- "BLOOD"
  # project_id <- projects[25]
  
  gtf_version <- 105
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, "/")
  
  # rse <- recount3::create_rse_manual(
  #   project = project_id,
  #   project_home = "data_sources/gtex",
  #   organism = "human",
  #   annotation = "gencode_v29",
  #   type = "jxn"
  # )
  # clusters_ID <- rse$gtex.smtsd %>% unique()
  # 
  # download_from_recount(project_id = project_id,
  #                       recount_version = 3,
  #                       clusters_ID = clusters_ID,
  #                       rse = rse,
  #                       folder_root = folder_root)
  # 
  # rm(rse)
  # gc()
  
  annotate_from_recount(project_id = project_id,
                        folder_root = folder_root)

  gc()
  
  idb_generation(project_id = project_id,
                 clusters_ID = readRDS(file = paste0(folder_root, "/raw_data/all_clusters_used.rds")),
                 folder_root = folder_root)
  
  
  gc()
  
}



