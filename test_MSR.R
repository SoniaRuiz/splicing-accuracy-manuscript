library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)
library(biomaRt)
library(DBI)
library(doParallel)
library(ggridges)

####################################################
## CONNECT TO THE SPLICING DATABASE ################
####################################################

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/splicing_accuracy_manuscript_figures.R")

project_id <- "BRAIN"
cluster_id <- "Brain - Frontal Cortex (BA9)"
print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))


## CONNECT TO THE OLD DATABASE ------------------------------

gtf_version <- 105

main_project <- paste0("splicing_1read")

database_path <- paste0(getwd(), "/database/", main_project, "/", gtf_version, "/", 
                        "/", main_project, "_backup.sqlite")


con_old <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con_old)

query = paste0("SELECT * FROM 'metadata'")
df_metadata_old <- dbGetQuery(con_old, query) %>% as_tibble()

all_projects_old <- df_metadata_old$SRA_project %>% unique



## Load data from FCTX


query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                FROM '", cluster_id, "_", project_id, "_nevermisspliced' ")
fctx_introns_old <- dbGetQuery(con_old, query) %>% as_tibble()
query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                FROM '", cluster_id, "_", project_id, "_misspliced' ")
fctx_introns_old <- rbind(fctx_introns_old, dbGetQuery(con_old, query) %>% as_tibble())

query <- paste0("SELECT DISTINCT intron.ref_junID, intron.protein_coding
                FROM 'intron' 
                WHERE intron.ref_junID IN (",
                paste(fctx_introns_old$ref_junID, collapse = ","),")")
fctx_introns_old <- fctx_introns_old %>%
  left_join(y = dbGetQuery(con_old, query) %>% as_tibble(),
            by = "ref_junID") %>% 
  as_tibble() 

fctx_introns_old

print(paste0(Sys.time(), " - ", fctx_introns_old %>% nrow(), " - introns!"))



## CONNECT TO THE NEW DATABASE ------------------------------

supporting_reads <- 1
gtf_version <- 105

main_project <- paste0("splicing_",supporting_reads,"read")

database_path <- paste0(getwd(), "/database/", main_project, "/", gtf_version, "/", 
                        "/", main_project, "_backup2.sqlite")


con_new <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con_new)

query = paste0("SELECT * FROM 'metadata'")
df_metadata_new <- dbGetQuery(con_new, query) %>% as_tibble()

all_projects_new <- df_metadata_new$SRA_project %>% unique


## Load data from FCTX

query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                FROM '", cluster_id, "_", project_id, "_nevermisspliced' ")
fctx_introns_new <- dbGetQuery(con_new, query) %>% as_tibble()
query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                FROM '", cluster_id, "_", project_id, "_misspliced' ")
fctx_introns_new <- rbind(fctx_introns_new, dbGetQuery(con_new, query) %>% as_tibble())

query <- paste0("SELECT DISTINCT intron.ref_junID, intron.protein_coding
                FROM 'intron' 
                WHERE intron.ref_junID IN (",
                paste(fctx_introns_new$ref_junID, collapse = ","),")")
fctx_introns_new <- fctx_introns_new %>%
  left_join(y = dbGetQuery(con_new, query) %>% as_tibble(),
            by = "ref_junID") %>% 
  as_tibble() 

fctx_introns_new






## SUBSAMPLE OLD INTRONS ------------------------------------------------------


fctx_introns_old_tidy <- fctx_introns_old %>%
  filter(protein_coding %in% c(0,100)) %>%
  mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC"))  %>%
  dplyr::rename(biotype = type_PC) %>%
  distinct(ref_junID, .keep_all = T) %>%
  group_by(ref_junID) %>%
  mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
  mutate(mean_coverage = mean_coverage %>% log10()) %>%
  ungroup()



  ## Subsampling introns to control by similarity in mean read coverage
  m.out <- MatchIt::matchit(biotype ~ mean_coverage, 
                            data = fctx_introns_old_tidy, 
                            distance = fctx_introns_old_tidy$mean_coverage,
                            method = "nearest", 
                            caliper = c(mean_coverage = 0.005), 
                            std.caliper = FALSE)
  subsample_old <- MatchIt::match.data(m.out)
  subsample_old %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)
  

## SUBSAMPLE NEW INTRONS ------------------------------------------------------
  
  fctx_introns_new_tidy <- fctx_introns_new %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC"))  %>%
    dplyr::rename(biotype = type_PC) %>%
    distinct(ref_junID, .keep_all = T) %>%
    group_by(ref_junID) %>%
    mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
    mutate(mean_coverage = mean_coverage %>% log10()) %>%
    ungroup()
  
  
  
  ## Subsampling introns to control by similarity in mean read coverage
  m.out.new <- MatchIt::matchit(biotype ~ mean_coverage, 
                                data = fctx_introns_new_tidy, 
                                distance = fctx_introns_new_tidy$mean_coverage,
                                method = "nearest", 
                                caliper = c(mean_coverage = 0.005), 
                                std.caliper = FALSE)
  subsample_new <- MatchIt::match.data(m.out.new)
  subsample_new %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)

  
  
  
## COMPARE FCTX NEW AND OLD DATABASE ------------------------------------------
  
  fctx_introns_new %>%
    filter(protein_coding %in% c(0,100)) %>%
    dplyr::count(protein_coding)
  
  fctx_introns_new %>%
    filter(protein_coding == 0)
  
  print(paste0(Sys.time(), " - ", fctx_introns_new %>% nrow(), " - NEW introns before subsampling!"))
  print(paste0(Sys.time(), " - ", subsample_new %>% nrow(), " - NEW introns after subsampling!"))
  
  wilcox.test(x = subsample_new %>% filter(protein_coding == 100) %>% pull(MSR_D),
              y = subsample_new %>% filter(protein_coding == 0) %>% pull(MSR_D),
              alternative = "less",
              paired = T,
              correct = T)
  
  ##
  
  
  fctx_introns_old %>%
    filter(protein_coding %in% c(0,100)) %>%
    dplyr::count(protein_coding)
  
  fctx_introns_old %>%
    filter(ref_junID == 216327)
  
  
  print(paste0(Sys.time(), " - ", fctx_introns_old %>% nrow(), " - OLD introns before subsampling!"))
  print(paste0(Sys.time(), " - ", subsample_old %>% nrow(), " - OLD introns after subsampling!"))
  
  
wilcox.test(x = subsample_old %>% filter(protein_coding == 100) %>% pull(MSR_D),
            y = subsample_old %>% filter(protein_coding == 0) %>% pull(MSR_D),
            alternative = "less",
            paired = T,
            correct = T)







## Study ref_junID == 216327



query <- paste0("SELECT DISTINCT ref_junID, ref_coordinates
                  FROM 'intron' ")
introns <- dbGetQuery(con, query) %>% as_tibble()

fctx_introns_old %>%
  filter(ref_junID == 132991) %>%
  inner_join(y = introns,
             by ="ref_junID")

fctx_introns_new %>%
  filter(ref_junID == 216327)








######################################################################################################
######################################################################################################

################################################################
## DONOR ADAR
################################################################


ENCODE_ENCORI_ADAR_MSR_D <- readRDS("~/PROJECTS/splicing-accuracy-manuscript/natcomms_review/results/ENCODE_ENCORI_ADAR_MSR_D.rds") %>%
  distinct(ref_coordinates, .keep_all =T)

table(ENCODE_ENCORI_ADAR_MSR_D$MSR_direction, ENCODE_ENCORI_ADAR_MSR_D$RBP_motif)

ENCODE_ENCORI_ADAR_MSR_D %>%
  dplyr::relocate(MSR_direction, .after = MSR_type) %>%
  dplyr::count(MSR_direction)




################################################################
## ACCEPTOR ADAR
################################################################

ENCODE_ENCORI_ADAR_MSR_A <- readRDS("~/PROJECTS/splicing-accuracy-manuscript/natcomms_review/results/ENCODE_ENCORI_ADAR_MSR_A.rds") %>%
  distinct(ref_coordinates, .keep_all =T)

table(ENCODE_ENCORI_ADAR_MSR_A$MSR_direction, ENCODE_ENCORI_ADAR_MSR_A$RBP_motif)

ENCODE_ENCORI_ADAR_MSR_A %>%
  dplyr::relocate(MSR_direction, .after = MSR_type) %>%
  dplyr::count(MSR_direction)

ENCODE_ENCORI_ADAR_MSR_A %>%
  dplyr::count(RBP_motif)


ENCORI_iCLIP_ADAR_ChiSquare <- readRDS("~/PROJECTS/splicing-accuracy-manuscript/natcomms_review/results/ENCORI_iCLIP_ADAR_ChiSquare.rds")

ENCORI_iCLIP_ADAR_ChiSquare %>% head()
ENCORI_iCLIP_ADAR_ChiSquare %>% dplyr::count(MSR_direction)

table(ENCORI_iCLIP_ADAR_ChiSquare$MSR_direction, ENCORI_iCLIP_ADAR_ChiSquare$RBP_iCLIP_type)
chisq <- chisq.test(table(ENCORI_iCLIP_ADAR_ChiSquare$MSR_direction, ENCORI_iCLIP_ADAR_ChiSquare$RBP_iCLIP_type))

chisq
