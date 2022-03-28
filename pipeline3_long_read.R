library(tidyverse)

###########################
## LOAD DATA
###########################

source("/home/sruiz/PROJECTS/splicing-project/pipeline0_prepare_data.R")
source("/home/sruiz/PROJECTS/splicing-project/pipeline1_QC_split_reads.R")

folder_path <- "/home/sruiz/PROJECTS/encode_lr/junctions/left_ventricle/"


add_encode_lr_data_to_idb <- function() {
  
  
  Sys.setenv("AWS_ACCESS_KEY_ID"=Sys.getenv("AWS_ACCESS_KEY_ID"),
             "AWS_SECRET_ACCESS_KEY"=Sys.getenv("AWS_SECRET_ACCESS_KEY"),
             "AWS_DEFAULT_REGION"=Sys.getenv("AWS_DEFAULT_REGION"))
  
  aws.s3::bucket_exists(bucket = "s3://encode-lr")
  raw_junc <- aws.s3::s3read_using(FUN = read.delim,
                                   object = "left_ventricle_ENCSR575LWI_junctions.txt",
                                   bucket = "s3://encode-lr")
  raw_junc <- raw_junc %>% 
    dplyr::mutate(junID = paste0(chrom, ":", genomic_start_coord, "-", genomic_end_coord, ":", strand))
  
  junc_tidy <- raw_junc %>% 
    select(seqnames = chrom,
           start = genomic_start_coord,
           end = genomic_end_coord,
           strand,
           junction_category) %>%
    as_tibble()
  saveRDS(object = junc_tidy,
          file = paste0(folder_path, "/all_split_reads.rds"))
  
  
  generate_all_split_reads(folder_path = folder_path, folder_destiny = folder_path)
}  

generate_metadata_encode_lr <- function() {
  
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
  
  gtexv8_project <- c("BLOOD_VESSEL",
                      "",
                      "",
                      "",
                      "KIDNEY",
                      "",
                      "HEART",
                      "HEART",
                      "HEART",
                      "",
                      "",
                      "HEART",
                      "LIVER", 
                      "HEART",
                      "HEART")
  
  gtexv8_tissue <- c("Artery - Aorta",
                     "",
                     "",
                     "",
                     "Kidney - Cortex","Kidney - Medulla",
                     "",
                     "Heart - Left Ventricle",
                     "Heart - Left Ventricle",
                     "Heart - Left Ventricle",
                     "",
                     "",
                     "Heart - Atrial Appendage",
                     "Liver",
                     "",
                     "")
  
  ## TODO: build the df & save it
  
  
}
raw_junc <- read.delim(file = paste0(folder_path, "/left_ventricle_ENCSR575LWI_junctions.txt"), header = T,sep = "\t")
raw_junc <- raw_junc %>% 
  dplyr::mutate(junID = paste0(chrom, ":", genomic_start_coord, "-", genomic_end_coord, ":", strand))

junc_tidy <- raw_junc %>% 
  select(seqnames = chrom,
         start = genomic_start_coord,
         end = genomic_end_coord,
         strand,
         junction_category) %>%
  as_tibble()
saveRDS(object = junc_tidy,
        file = paste0(folder_path, "/all_split_reads.rds"))


generate_all_split_reads(folder_path = folder_path, folder_destiny = folder_path)




## Load reference transcriptome Ensembl v105
blacklist_path = "/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed"
gtf_path <- "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf"
edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))

all_split_reads_tidy <- readRDS(file = paste0(folder_path, "/all_split_reads.rds"))


all_split_reads_LR <- annotate_dasper(all_split_reads = all_split_reads_tidy %>% 
                                        GenomicRanges::GRanges() %>% 
                                        diffloop::rmchr(),
                                      edb = edb,
                                      blacklist_path = blacklist_path)

all_split_reads_LR %>%
  dplyr::count(type) 

all_split_reads_LR_novel <- all_split_reads_LR %>%
  filter(type %in% c("novel_donor", "novel_acceptor")) %>%
  dplyr::mutate(junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand))
  


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

          