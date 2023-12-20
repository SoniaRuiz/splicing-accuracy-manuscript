library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(Biostrings)
library(tidyverse)
library(protr)

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/init.R")

base_folder <- here::here()
dependencies_folder <- paste0(base_folder, "/dependencies/")

#####################################
## LOAD SOURCE SCRIPTS
#####################################

setwd(file.path(base_folder,"scripts"))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(base_folder))


#####################################
## SET MAIN VARIABLES
#####################################

## This is the Ensembl gtf transcriptome version 
gtf_versions <- c(105)

## This is the name of the project producing the database
supportive_reads <- 1
project_name <- paste0("splicing_", supportive_reads, "read")
data_source <- "data_sources/gtex"

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/gene_expression.R")
source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/scripts/27_age_effect_uncorrected_TPM_lm.R")
source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/scripts/99_utils.R")

#################################
## SET VARIABLES
#################################

# source(file.path(here::here("database_SQL_helper.R")))
# dependencies_folder <- file.path(here::here("dependencies"))
# 
gtf_version <- 105

database_folder <- paste0(base_folder, "/database/", project_name, "/", gtf_version, "/")
results_folder <- paste0(base_folder, "/results/", project_name, "/", gtf_version, "/")
levelqc1_folder <- paste0(base_folder, "/database/")
tpm_folder <- paste0(base_folder, "/results/tpm/")


recount3_project_IDs <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))

# project_name <- "splicing_1read"
# 
# results_folder <- file.path(here::here("results"), paste0(project_name, "/", gtf_version, "/"))

gene_type <- "nmd"

if (gene_type == "nmd") {
  gene_list <- data.frame(id = c("ENSG00000005007", "ENSG00000151461", "ENSG00000157106", "ENSG00000070366", "ENSG00000116698"),
                          name = c("UPF1", "UPF2", "SMG1", "SMG6", "SMG7"))
} else {
  gene_list <- all_RBPs <- xlsx::read.xlsx(file = file.path(here::here(dependencies_folder, '/RBPs_subgroups.xlsx')), 
                                           header = TRUE,
                                           sheetIndex = 1) %>% 
    as_tibble() %>%
    distinct(name, .keep_all = T)
}

# projects_id <- c( "ADIPOSE_TISSUE", "ADRENAL_GLAND", "BLOOD", "BLOOD_VESSEL", "BRAIN", "COLON", 
#                   "ESOPHAGUS", "HEART", "LUNG", "MUSCLE", "NERVE", "PANCREAS", "SALIVARY_GLAND", 
#                   "SKIN", "SMALL_INTESTINE", "SPLEEN", "STOMACH", "THYROID" )

################################
## FUNCTIONS       
################################


## 1. GET RECOUNT 3 EXPRESSION DATA AND METADATA FOR THE INDICATED TISSUE

if ( !file.exists(file.path(results_folder, paste0(gene_type, "_genes_age_lm.rds")))) {
  
  gene_age_lm <- map_df(recount3_project_IDs, function(project_id) {
    
    # project_id <- recount3_project_IDs[5]
    
    message(Sys.time(), " - ", project_id)
    local_results_folder <- file.path(results_folder, project_id)
    
    
    if ( file.exists(paste0(local_results_folder, "/", gene_type, "_expression/", project_id, "_tpm_ensembl",gtf_version,".rds")) &&
         file.exists(paste0(local_results_folder, "/", gene_type, "_expression/", project_id, "_covariates.rds")) ) {
      
      message(Sys.time(), " - loading data for ", project_id, "...")
      
      dds_tpm <- readRDS(file = paste0(local_results_folder, "/", gene_type, "_expression/", project_id, "_tpm_ensembl",gtf_version,".rds"))
      sample_metadata <- readRDS(file = paste0(local_results_folder,  "/", gene_type, "_expression/", project_id, "_covariates.rds"))
      
      
    } else {
      
      dir.create(file.path(local_results_folder, paste0(gene_type, "_expression")), recursive = TRUE, showWarnings = T)
      
      ## Get metadata for local tissue
      metadata <- readRDS(file = paste0(local_results_folder, "/base_data/", project_id, "_samples_raw_metadata.rds"))
      
      metadata$gtex.smafrze %>% unique
      metadata$gtex.smrin %>% unique %>% min
      
      
      ## Get base TPM for local tissue
      recount_tpm <- readRDS(file = file.path(here::here("results/tpm/"), paste0(project_id, "_tpm.rds"))) %>%
        as_tibble(rownames = "gene")
      
      
      recount_tpm %>% head
      recount_tpm %>% nrow()
      recount_tpm %>% ncol()
      
      
      ## Tidy the object
      
      ## Filter by the samples of the current cluster
      dds_tpm <- recount_tpm %>% 
        dplyr::select(gene, all_of(metadata$external_id)) %>%
        mutate(gene = gsub(pattern = "\\..*", replacement = "", x = gene)) %>%
        filter(gene %in% gene_list$id)
      
      
      saveRDS(object = dds_tpm,
              file = paste0(local_results_folder, "/", gene_type, "_expression/", project_id, "_tpm_ensembl",gtf_version,".rds"))
      
      
      # 2. Get the covariates to correct for
      sample_metadata <- tidy_sample_metadata(sample.metadata = metadata, 
                                              samples = metadata$external_id) %>% 
        as_tibble(rownames = "covariates")
      
      saveRDS(object = sample_metadata, 
              file = paste0(local_results_folder,  "/", gene_type, "_expression/", project_id, "_covariates.rds"))
      
      
      rm(recount_tpm)
      #gc()
      
    }
    
    
    ##########################################
    # 2. ANALYSIS
    # Test if TPM values are significantly affected by age
    ##########################################
    
    message(Sys.time(), " - getting linear models for ", project_id, "...")
    
    
    lm_output <- age_effect_uncorrected_TPM_lm(project.id = project_id,
                                               tpm.uncorrected = dds_tpm,
                                               sample.metadata = sample_metadata,
                                               gene.list = gene_list,
                                               results.folder = local_results_folder)
    
    return(lm_output %>%
             mutate(project = project_id))
    
    
  })
  
  
  
  saveRDS(object = gene_age_lm %>% as_tibble(),
          file = file.path(results_folder, paste0(gene_type, "_genes_age_lm.rds")) ) 
} else {
  gene_age_lm <- readRDS(file = file.path(results_folder, paste0(gene_type, "_genes_age_lm.rds")) ) 
}



##################################
## GET STATS FOR TISSUES
##################################

tissues_data_tidy <- gene_age_lm %>%
  filter(covariate == "gtex.age") 

tissues_data_tidy %>%
  filter(q <= 0.05, Estimate < 0) %>%
  group_by(project) %>%
  dplyr::count()

tissues_data_tidy %>%
  filter(q <= 0.05, Estimate < 0) %>%
  pull(q) %>%
  summary()
  
((tissues_data_tidy %>%
    filter(q <= 0.05) %>%
    nrow) * 100) / (tissues_data_tidy %>% nrow)

tissues_data_tidy %>%
  filter(q <= 0.05) %>%
  pull(q) %>%
  summary()

##################################
## GET STATS FOR BRAIN
##################################

brain_data_tidy <- gene_age_lm %>%
  filter(covariate == "gtex.age",
         project == "BRAIN") %>%
  distinct(name, .keep_all =T)

write.csv(x = brain_data_tidy %>%
            dplyr::select(name, Estimate, q, pval) %>%
            as.data.frame,
          file = paste0(results_folder, "/_paper_review/results/brain_age_RPBs_lm.csv"))

((brain_data_tidy %>%
    filter(q <= 0.05) %>%
    nrow) * 100) / (brain_data_tidy %>% nrow)

brain_data_tidy %>%
  filter(q <= 0.05) %>%
  pull(q) %>%
  summary()

##################################
## PLOT
##################################

## Check the NMD factors with TPM values affected by age
gene_age_lm_tidy <- gene_age_lm %>%
  filter(covariate == "gtex.age") %>%
  #filter(Estimate < 0) %>%
  mutate(Estimate = ifelse(Estimate < 0, NA, Estimate)) %>%
  mutate(q = ifelse(q > 0.05, NA, q)) %>%
  mutate(`log10(q)` = q %>% log10()) %>%
  group_by(name) %>%
  mutate(name = factor(name, levels=(gene_age_lm$name)[order(gene_age_lm$name %>% unique %>% dplyr::desc())] %>% unique)) 

#gene_age_lm_tidy$name <- factor(gene_age_lm_tidy$name, levels=(gene_age_lm_tidy$name)[order(gene_age_lm_tidy$Estimate)])

ggplot() + 
  geom_tile(data = gene_age_lm_tidy, 
            mapping = aes(x = name, y = fct_rev(project), fill = `log10(q)`, colour = "q>=0.05")) +
  scale_fill_gradient(low = "red", high = "white", na.value = '#cccccc') +
  scale_colour_manual(values = c( "q>=0.05" = "#cccccc")) +
  xlab("") + 
  ylab("") + 
  theme_light()  +
  custom_ggtheme +
  #theme(legend.box.margin=margin(b = -11, l = -20)) +
  guides(colour = guide_legend(override.aes = list(fill = '#cccccc'),
                               #title = "log10 q",
                               label.position = "bottom",
                               order = 2)) +
  
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        legend.position = "top",
        axis.text.y = element_text(size = 5, colour = "black"))

file_name <- paste0(results_folder, "/_paper_review/figures/",gene_type,"_age.png")
ggplot2::ggsave(filename = file_name, width = 90, height = 90, units = "mm", dpi = 300)
