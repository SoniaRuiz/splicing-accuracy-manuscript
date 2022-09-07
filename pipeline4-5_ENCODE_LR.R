library(tidyverse)

###########################
## LOAD DATA
## Junctions generated from source data: 
## https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=Transcription&assay_title=long+read+RNA-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=GRCh38&files.platform.term_name=Pacific+Biosciences+Sequel+II&biosample_ontology.classification=tissue
###########################

source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline1_download_from_recount.R")
source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline2_annotate_from_recount.R")


folder_root <- "/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/junctions/"


generate_metadata_encode_LR <- function() {
  
  encode_tissue <- c("aorta", "aorta", 
                     "astrocyte",
                     "cardiac_septum",
                     "cortical_neuron",
                     "kidney",
                     "left_cardiac_atrium",
                     "left_ventricle", "left_ventricle", "left_ventricle", "left_ventricle",
                     "left_ventricle_myocard_sup",
                     "left_ventricle_myocardium_inferior",
                     "neural_crest_cell",
                     "posterior_vena_cava", "posterior_vena_cava",
                     "right_cardiac_atrium", "right_cardiac_atrium", "right_cardiac_atrium", "right_cardiac_atrium",
                     "right_lobe_liver", "right_lobe_liver",
                     "right_ventricle", "right_ventricle", "right_ventricle",
                     "right_ventricle_myocardium_inferior", "right_ventricle_myocardium_inferior")
  
  gtexv8_tissue <- c("Artery - Aorta", "Artery - Aorta",
                     "Brain - Frontal Cortex (BA9)",
                     "Artery - Coronary",
                     "Brain - Frontal Cortex (BA9)",
                     "Kidney - Cortex", 
                     "Artery - Coronary",
                     "Heart - Left Ventricle","Heart - Left Ventricle","Heart - Left Ventricle","Heart - Left Ventricle",
                     "Heart - Left Ventricle",
                     "Heart - Left Ventricle",
                     "Brain - Frontal Cortex (BA9)",
                     "Artery - Coronary","Artery - Coronary",
                     "Heart - Atrial Appendage","Heart - Atrial Appendage","Heart - Atrial Appendage","Heart - Atrial Appendage",
                     "Liver","Liver",
                     "Heart - Atrial Appendage","Heart - Atrial Appendage","Heart - Atrial Appendage",
                     "Heart - Atrial Appendage","Heart - Atrial Appendage")
  
  encode_sample <- c("ENCSR425HFS", "ENCSR700EBI", 
                     "ENCSR071IHY",
                     "ENCSR549ELD",
                     "ENCSR257JBF",
                     "ENCSR432YKA",
                     "ENCSR424QFN",
                     "ENCSR194YUY", "ENCSR575LWI", "ENCSR700XDQ", "ENCSR994YZY",
                     "ENCSR777CCI",
                     "ENCSR786FLO",
                     "ENCSR056MYH",
                     "ENCSR138TAS", "ENCSR853YZN",
                     "ENCSR435UUS", "ENCSR514YQN", "ENCSR553SVP", "ENCSR728TXV",
                     "ENCSR293MOX", "ENCSR657HJW",
                     "ENCSR329ZQG", "ENCSR782LGT", "ENCSR984OAE",
                     "ENCSR591OZR", "ENCSR899GAP")
  
  gtexv8_project <- c("BLOOD_VESSEL","BLOOD_VESSEL",
                      "BRAIN",
                      "BLOOD_VESSEL",
                      "BRAIN",
                      "KIDNEY",
                      "BLOOD_VESSEL",
                      "HEART","HEART","HEART","HEART",
                      "HEART",
                      "HEART",
                      "BRAIN",
                      "BLOOD_VESSEL","BLOOD_VESSEL",
                      "HEART","HEART","HEART","HEART",
                      "LIVER", "LIVER",
                      "HEART","HEART","HEART",
                      "HEART","HEART")
  
  
  
  saveRDS(object = data.frame(encode_tissue,
                              encode_sample,
                              gtexv8_project,
                              gtexv8_tissue),
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/metadata.rds")## TODO: build the df & save it
  
  
}

annotate_split_reads_encode_LR <- function() {
  
  bucket_name <- "s3://encode-lr"
  gtf_version <- 105
  edb <-ensembldb::ensDbFromGtf(gtf = URLencode(paste0("http://ftp.ensembl.org/pub/release-", gtf_version, 
                                                       "/gtf/homo_sapiens/Homo_sapiens.GRCh38.", 
                                                       gtf_version, ".gtf.gz")), 
                                outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
  
  if (!aws.s3::bucket_exists(bucket = bucket_name)) {
    print("ERROR! bucket doesn't exist.")
  }
  
  Sys.setenv("AWS_ACCESS_KEY_ID"=Sys.getenv("AWS_ACCESS_KEY_ID"),
             "AWS_SECRET_ACCESS_KEY"=Sys.getenv("AWS_SECRET_ACCESS_KEY"),
             "AWS_DEFAULT_REGION"=Sys.getenv("AWS_DEFAULT_REGION"))
  
  
  for (file in (aws.s3::get_bucket(bucket = bucket_name))) {
    
    # file <- (aws.s3::get_bucket(bucket = bucket_name))[1]
    
    folder_path <- paste0(folder_root, "/", sub("\\_ENC.*", "", file$Key))
    dir.create(path = folder_path, showWarnings = T)
    
    if (!file.exists(paste0(folder_path, "/all_split_reads.rds"))) {
      
      ##############################
      ## EXTRACT THE SPLIT READS
      ##############################
      
      raw_junc <- aws.s3::s3read_using(FUN = read.delim,
                                       object = file$Key,
                                       bucket = bucket_name)
      
      raw_junc <- raw_junc %>% 
        dplyr::mutate(junID = paste0(chrom, ":", 
                                     genomic_start_coord, "-", 
                                     genomic_end_coord, ":", strand))
      
      junc_tidy <- raw_junc %>% 
        select(seqnames = chrom,
               start = genomic_start_coord,
               end = genomic_end_coord,
               strand,
               junction_category,
               junID) %>%
        as_tibble()
    
      saveRDS(object = junc_tidy,
              file = paste0(folder_path, "/all_split_reads.rds"))
      
      generate_all_split_reads(folder_path = folder_path, folder_destiny = folder_path)
      
    } else {
      junc_tidy <- readRDS(file = paste0(folder_path, "/all_split_reads.rds"))
      print(paste0(Sys.time(), " - '", paste0(folder_path, "/all_split_reads.rds"), "' exists!"))
    }
    
    
    ##############################
    ## ANNOTATE THE SPLIT READS
    ##############################
    
    get_base_data_annotated(cluster = sub("\\_ENC.*", "", file$Key),
                            all_split_reads = junc_tidy,
                            blacklist_path = "/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed",
                            edb = edb,
                            gtf_version = gtf_version,
                            folder_save_path = folder_path)
  }
  
}  
  
crossing_SR_LR <- function() {
  
  ##################################
  ## CONNECT TO THE DATABASE
  ##################################
  
  database_path <- "/home/sruiz/PROJECTS/splicing-project-recount3/database/introverse_qc.sqlite"
  con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
  DBI::dbListTables(con)
  
  
  ##################################
  ## CROSSING SR AND LR DATA
  ##################################
  
  df_crossing_metadata <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/metadata.rds")
  
  df_crossing_data <- map_df(df_crossing_metadata$encode_sample, function(sample) {
  
    # sample <- "ENCSR425HFS"
    
    #########################################
    ## LOAD THE ENCODE LR DATA
    #########################################
    
    encode_sample <- df_crossing_metadata %>%
      filter(encode_sample == sample)
    
    
    print(paste0(Sys.time(), " - ", sample))
    all_split_reads_LR_novel <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/junctions/",
                                                      encode_sample$encode_tissue, "/", encode_sample$encode_tissue, 
                                                      "_annotated_SR_details_length_105.rds"))
    
    all_split_reads_LR_novel <- all_split_reads_LR_novel %>%
      GRanges() %>%
      diffloop::addchr() %>%
      as_tibble() %>%
      filter(type != "annotated")
    
    all_split_reads_LR_novel
    
    if (any(all_split_reads_LR_novel$junction_category != "novel")) {
      print("ERROR! Some of the LR novel junctions have been originally classified as not novel.")
    }
    #########################################
    ## LOAD THE GTEx SR DATA
    #########################################
    
    query = paste0("SELECT tissue.*, gene.gene_name, novel.novel_type, novel.novel_coordinates, intron.ref_coordinates
                 FROM '", encode_sample$gtexv8_tissue, "_", encode_sample$gtexv8_project, "_misspliced' as tissue
                 INNER JOIN 'gene' ON tissue.gene_id=gene.id
                 INNER JOIN 'intron' ON tissue.ref_junID=intron.ref_junID
                 INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID")
    db_tissue <- DBI::dbGetQuery(con, query) %>% 
      dplyr::select(-gene_id,-novel_junID, -ref_junID) %>% as_tibble() %>%
      dplyr::relocate(ref_coordinates,novel_coordinates) %>%
      dplyr::relocate(novel_type,.after = "novel_sum_counts")
    db_tissue %>% as_tibble() 
    
    
    #########################################
    ## CROSS SR AND LR DATA
    #########################################
    
    SR_hits <- db_tissue %>%
      filter(novel_coordinates %in% all_split_reads_LR_novel$junID)
    
    LR_hits <- all_split_reads_LR_novel %>%
      filter(junID %in% db_tissue$novel_coordinates) 

    
    
    df_SR_hits <- SR_hits %>%
      mutate(GTEx_tissue = encode_sample$gtexv8_tissue,
             encode_sample_id = sample,
             encode_tissue = encode_sample$encode_tissue) %>%
      as_tibble() 
    
    print(paste0(sample, '_GTEx_', str_remove_all(string = encode_sample$gtexv8_tissue,pattern = " "), 
                 '-ENCODE_', encode_sample$encode_tissue))
    
    if (!file.exists("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/GTExv8_SR-ENCODE_LR.xlsx")) {
      print("it doesn't exist")
      
      xlsx::write.xlsx(x = df_SR_hits %>% as.data.frame(),
                       sheetName = paste0(sample, '_GTEx_', str_remove_all(string = encode_sample$gtexv8_tissue,pattern = " "), 
                                          '-ENCODE_', encode_sample$encode_tissue),
                       row.names = F,
                       file = "/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/GTExv8_SR-ENCODE_LR.xlsx")
      
    } else {
      xlsx::write.xlsx(x = df_SR_hits %>% as.data.frame(),
                       sheetName = paste0(sample, '_GTEx_', str_remove_all(string = encode_sample$gtexv8_tissue,pattern = " "),
                                          '-ENCODE_', encode_sample$encode_tissue),
                       col.names = T,
                       row.names = F,
                       append = T,
                       file = "/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/GTExv8_SR-ENCODE_LR.xlsx")
    }
    
    df_SR_hits %>%
      return()

  })
  
  xlsx::write.xlsx(x = df_crossing_data %>% as.data.frame(),
                   sheetName = "all",
                   row.names = F,
                   append = T,
                   file = "/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/long-read/GTExv8_SR-ENCODE_LR.xlsx")
  
 
}

################################

con <- DBI::dbConnect(RSQLite::SQLite(), "/home/sruiz/PROJECTS/splicing-project-app/intron_db/dependencies/splicing.sqlite")
query = paste0("SELECT * FROM 'Brain - Hippocampus_BRAIN_db_novel'")
db_novel <- DBI::dbGetQuery(con, query) %>% as_tibble() %>%
  dplyr::rename(gene_id = gene_name)
query = paste0("SELECT * FROM 'gene_name'")
db_genes <- DBI::dbGetQuery(con, query) %>% as_tibble()
DBI::dbDisconnect(conn = con)

db_novel <- merge(x = db_novel,
                  y = db_genes,
                  by = "gene_id",
                  all.x = T) %>%
  select(-gene_id)
  
db_novel %>%
  distinct(novel_junID) %>%
  count()

db_novel %>%
  distinct(gene_name)

#####################################
## SAVE RESULTS
#####################################

SR_hits <- db_novel %>%
  filter(novel_junID %in% all_split_reads_LR_novel$junID)

LR_hits <- all_split_reads_LR_novel %>%
  filter(junID %in% db_novel$novel_junID) %>%
  select(
         -gene_id_start,
         -gene_name_start,
         -symbol_start,
         -tx_id_start,
         -exon_id_start,
         -strand_start,
         -exon_width_start,
         -gene_id_end,
         -gene_name_end,
         -symbol_end,
         -tx_id_end,
         -exon_id_end,
         -strand_end,
         -exon_width_end,
         -tx_id_junction) %>% 
  mutate( gene_id_junction =  gene_id_junction %>% as.character())

LR_hits <- merge(x = LR_hits,
                 y = raw_junc %>% 
                   select(junID, isoform),
                 by = "junID", 
                 all.x = T) %>%
  distinct(junID, .keep_all = T)
 

intersect(SR_hits$novel_junID, LR_hits$junID)
xlsx::write.xlsx(x = SR_hits,
                 sheetName = 'GTExv8_SR_Hippocampus',
                 row.names = F,
                 file = "/home/sruiz/PROJECTS/encode_lr/junctions/left_ventricle/GTExv8-Hippocampus_ENCODELR-ArteryAorta.xlsx")
xlsx::write.xlsx(x = LR_hits ,
                 sheetName = 'ENCODE_LR_Artery_Aorta',
                 append = T,
                 row.names = F,
                 showNA = F,
                 file = "/home/sruiz/PROJECTS/encode_lr/junctions/left_ventricle/GTExv8-Hippocampus_ENCODELR-ArteryAorta.xlsx")

xlsx::write.xlsx(x = merge(x = SR_hits,
                           y = LR_hits,
                           by.x = "novel_junID",
                           by.y = "junID"), 
                 sheetName = 'All',
                 append = T,
                 row.names = F,
                 showNA = F,
                 file = "/home/sruiz/PROJECTS/encode_lr/junctions/left_ventricle/GTExv8-Hippocampus_ENCODELR-ArteryAorta.xlsx")

          