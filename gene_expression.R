library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(DBI)
library(here)


# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/gene_expression.R")


# source(file.path(here::here("database_SQL_helper.R")))
dependencies_folder <- file.path(here::here("dependencies"))

gtf_version <- 105
project_name <- "splicing_2reads"


################################
## FUNCTIONS       
################################

tidy_sample_metadata <- function(sample.metadata,
                                 samples) {
  
 
  
  age_numeric <- as.numeric(factor(as.matrix(sample.metadata$gtex.age))) 
  sample.metadata$gtex.age <- age_numeric
 
  sample.metadata <- sample.metadata %>%
    tibble::column_to_rownames(var = "external_id")
  
  covariates <- c("gtex.age", "gtex.smcenter", "gtex.smtsd",
                  "gtex.smgebtch", "gtex.smgebtchd",
                  "gtex.smnabtch", "gtex.smnabtchd", "gtex.smnabtcht",
                  "gtex.dthhrdy", "gtex.sex", "gtex.smrin")
  
  
  smcenter <- as.numeric(factor(as.matrix(sample.metadata$gtex.smcenter))) 
  sample.metadata$gtex.smcenter <- smcenter
  
  gtex.smgebtch <- as.numeric(factor(as.matrix(sample.metadata$gtex.smgebtch))) 
  sample.metadata$gtex.smgebtch <- gtex.smgebtch
  
  gtex.smgebtchd <- as.numeric(factor(as.matrix(sample.metadata$gtex.smgebtchd))) 
  sample.metadata$gtex.smgebtchd <- gtex.smgebtchd
  
  gtex.smnabtch <- as.numeric(factor(as.matrix(sample.metadata$gtex.smnabtch))) 
  sample.metadata$gtex.smnabtch <- gtex.smnabtch
  
  gtex.smnabtchd <- as.numeric(factor(as.matrix(sample.metadata$gtex.smnabtchd))) 
  sample.metadata$gtex.smnabtchd <- gtex.smnabtchd
  
  gtex.smnabtcht <- as.numeric(factor(as.matrix(sample.metadata$gtex.smnabtcht))) 
  sample.metadata$gtex.smnabtcht <- gtex.smnabtcht
  
  gtex.smtsd <- as.numeric(factor(as.matrix(sample.metadata$gtex.smtsd))) 
  sample.metadata$gtex.smtsd <- gtex.smtsd
  
  
  # 8. Return covariates ---------------------------
  
  return(t(sample.metadata %>%
           dplyr::select(all_of(covariates))))
}

age_effect_uncorrected_TPM_lm <- function(project.id,
                                          tpm.uncorrected,
                                          sample.metadata,
                                          gene.list) {
  
  ################################
  ## LOAD GENES OF INTEREST
  ################################
  
  ## Only keep the TPM value of the GENES of interest
  tpm.uncorrected <- tpm.uncorrected %>%
    filter(gene %in% gene.list$id)
 
  
  ## Check the samples used are the same
  if ( !( identical(x = colnames(tpm.uncorrected)[-1] %>% unique %>% sort(),
                    y = colnames(sample.metadata)[-1] %>% unique %>% sort()) ) ) {
    print("ERROR!")
    break;
  }
  
  
  genes_of_interest <- (tpm.uncorrected$gene) %>% unique
    

  ## Per each gene, obtain TPM and metadata values and run lm
  df_lm_output <- map_df(genes_of_interest, function(gene_id) {
    
    # gene_id <- genes_of_interest[1]
    # gene_id <- "ENSG00000002079"
    
    print(gene_id)
    
    tpm <- tpm.uncorrected %>%
      filter(gene == gene_id) %>%
      mutate(gene = "tpm") %>%
      as_tibble()
    
    if ( rowSums(tpm[, c(2:ncol(tpm))], na.rm = T) != 0 ) {

      
      ## Transform TPM to log10 so resilduals after modelling in the linear model are normally distributed
      ## https://stats.stackexchange.com/questions/40907/in-regression-analysis-what-does-taking-the-log-of-a-variable-do
      df_gene <- rbind(tpm %>% dplyr::rename(covariates = "gene"),
                  sample.metadata) %>%
        gather(sample, value, -covariates)  %>%
        drop_na() %>%
        spread(key = covariates, value = value) %>%
        drop_na() %>%
        dplyr::mutate(gtex.age = gtex.age %>% as.double(),
                      gtex.dthhrdy = gtex.dthhrdy %>% as.double(),
                      gtex.sex = gtex.sex %>% as.double(),
                      gtex.smcenter = gtex.smcenter %>% as.double(),
                      gtex.smgebtch = gtex.smgebtch %>% as.double(),
                      gtex.smgebtchd = gtex.smgebtchd %>% as.double(),
                      gtex.smnabtch = gtex.smnabtch %>% as.double(),
                      gtex.smnabtchd = gtex.smnabtchd %>% as.double(),
                      gtex.smnabtcht = gtex.smnabtcht %>% as.double(),
                      gtex.smrin = gtex.smrin %>% as.double(),
                      tpm = (tpm %>% as.double()) %>% log10) %>%
        as_tibble() %>%
        filter_all(all_vars(!is.infinite(.)))
      
      # print(df)
      lm_output <- lm(tpm ~ 
                        gtex.smrin + 
                        gtex.smcenter +
                        gtex.smgebtch +
                        gtex.smgebtchd +
                        gtex.smnabtch +
                        gtex.smnabtchd +
                        gtex.smnabtcht +
                        gtex.dthhrdy +
                        gtex.sex +
                        gtex.smtsd +
                        gtex.age,
                      data = df_gene)
      
      lm_output <- lm_output %>% summary()
      df_lm_output <- lm_output$coefficients %>% as_tibble(rownames = "covariate")

    
    df_lm_output <- df_lm_output %>%
      mutate(gene_id = gene_id) %>%
      dplyr::select(covariate, Estimate, pval = `Pr(>|t|)`, gene_id)
    
    
    # plot(density(lm_output$residuals)) 
    
    return(df_lm_output)
    
  } else {
    
    return(NULL)
    
  }
    
  }) 
    
    
  ## Filter by age covariate and order by estimate size
  df_lm_output_tidy <- df_lm_output %>%
    mutate(q = p.adjust(p = pval, method = "fdr"))  %>%
    arrange(Estimate)
  

  ## Add gene SYMBOL info
  df_lm_output_tidy <- left_join(x = df_lm_output_tidy, 
                                 y = gene.list %>% drop_na(), 
                                 by = c("gene_id" = "id")) %>%
    group_by(gene_id) 
  
 
  
  ## Save results
  write_csv(x = df_lm_output_tidy,
            file = file.path(here::here("results/_paper/results/", 
                                        paste0(project.id, "_", main_project,"_tpm_lm_all_ensembl105.csv"))))
  
  
  ###############################################
  ## Get GENE that decrease expression with AGE
  ###############################################
  
  
  ## Genes with expression levels decreasing with age (i.e. with negative effect over their TPM associated to the 'age' covariate)
  genes_age_lm <- df_lm_output_tidy %>%
    filter(str_detect(string = covariate, pattern = "gtex.age")) %>%
    #filter(Estimate < 0) %>%
    #filter(pval <= 0.05) %>%
    arrange(name)
  
  return(genes_age_lm %>%
           mutate(project = project.id))
  
  
}



###########################################
## 1. GET RECOUNT 3 EXPRESSION DATA
##    AND METADATA FOR THE INDICATED TISSUE
###########################################

gene_type <- "NMD"

if (gene_type == "NMD") {
  gene_list <- data.frame(id = c("ENSG00000005007", "ENSG00000151461", "ENSG00000157106", "ENSG00000070366", "ENSG00000116698"),
                          name = c("UPF1", "UPF2", "SMG1", "SMG6", "SMG7"))
} else {
  gene_list <- all_RBPs <- xlsx::read.xlsx(file = file.path(here::here(dependencies_folder, '/RBPs_subgroups.xlsx')), 
                                           header = TRUE,
                                           sheetIndex = 1) %>% 
    as_tibble() %>%
    distinct(name, .keep_all = T)
}

projects_id <- c( "ADIPOSE_TISSUE", "ADRENAL_GLAND", "BLOOD", "BLOOD_VESSEL", "BRAIN", "COLON", 
                  "ESOPHAGUS", "HEART", "LUNG", "MUSCLE", "NERVE", "PANCREAS", "SALIVARY_GLAND", 
                  "SKIN", "SMALL_INTESTINE", "SPLEEN", "STOMACH", "THYROID" )


if ( file.exists(file.path(here::here("results"), paste0("NMD_genes_age_lm.rds"))) ) {
  
  gene_age_lm <- readRDS(file = file.path(here::here("results"), paste0("NMD_genes_age_lm.rds"))) 
  
} else {
  
  gene_age_lm <- map_df(projects_id, function(project_id) {
    
    results_folder <- file.path(here::here("results"), paste0(project_name, "/", gtf_version, "/"))
    
    # project_id <- projects_id[1]
    
    ####################################
    # 1. GENE EXPRESSION PRE-PROCESSING
    ####################################
    
    message(Sys.time(), " - ", project_id)
    folder_name <- file.path(here::here("results/RBP_EXPRESSION/splicing", project_id))
    
    
    if ( file.exists(paste0(folder_name, "/", project_id, "_tpm_ensembl105.csv")) && 
         file.exists(paste0(folder_name, "/", project_id, "_covariates.csv"))) {
      
      dds_tpm <- read.csv(file = paste0(folder_name, "/", project_id, "_tpm_ensembl105.csv"), header = T)
      sample_metadata <- read.csv(file = paste0(folder_name,  "/", project_id, "_covariates.csv"), header = T)
      
      
    } else {
      
      dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)
      
      ## Get metadata for brain tissue
      metadata <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", 
                                        gtf_version, "/",main_project,"/base_data/", 
                                        project_id, "_samples_metadata.rds"))
      
      metadata$gtex.smafrze %>% unique
      metadata$gtex.smrin %>% unique %>% min
      
      
      ## Get TPM for brain tissue
      
      recount_tpm <- generate_recount3_tpm(recount3.project.IDs = project_id,
                                           gtf.version = gtf_version,
                                           data.source = "data_sources/gtex", 
                                           results.folder = results_folder)
      
      recount_tpm %>% head
      recount_tpm %>% nrow()
      recount_tpm %>% ncol()
      
      
      
      ## 1. Get the TPM expression using GTEX data
      
      ## Filter by the samples of the current cluster
      dds_tpm <- recount_tpm %>% 
        dplyr::select(gene, all_of(metadata$external_id))
      
      
      write.csv(x = dds_tpm,
                file = paste0(folder_name, "/", project_id, "_tpm_ensembl105.csv"))
      
      
      # 2. Get the covariates to correct for
      sample_metadata <- tidy_sample_metadata(metadata, 
                                              samples = metadata$external_id) %>% 
        as_tibble(rownames = "covariates")
      
      write_csv(x = sample_metadata %>% as_tibble(rownames = "covariates"), 
                file = paste0(folder_name,  "/", project_id, "_covariates.csv"))
      
      
      rm(recount_tpm)
      gc()
      
    }
    
    
    ##########################################
    # 2. ANALYSIS
    # Test if TPM values are significantly affected by age
    ##########################################
    
    
    lm_output <- age_effect_uncorrected_TPM_lm(project.id = project_id,
                                               tpm.uncorrected = dds_tpm,
                                               sample.metadata = sample_metadata,
                                               gene.list = gene_list)
    
    return(lm_output)
    
  })
  
  
  
  saveRDS(object = gene_age_lm %>% as_tibble(),
          file = file.path(here::here("results"), paste0("NMD_genes_age_lm.rds")) ) 
}


## Check the NMD factors with TPM values affected by age
gene_age_lm %>%
  mutate(q = ifelse(q > 0.05, NA, q)) %>%
  mutate(`log10(q)` = q %>% log10())%>%
  arrange(desc(Estimate)) %>%
  as.data.frame()


ggplot() + 
  geom_tile(data = gene_age_lm %>%
              mutate(q = ifelse(q > 0.05, NA, q)) %>%
              mutate(`log10(q)` = q %>% log10())%>%
              arrange(desc(Estimate)) %>%
              as.data.frame(), 
            mapping = aes(y = project, x = name, fill = `log10(q)`, colour = "q>=0.05")) +
  scale_fill_gradient(low = "red", high = "white", na.value = '#cccccc') +
  scale_colour_manual(values = c( "q>=0.05" = "#cccccc")) +
  xlab("NMD gene") + 
  ylab("GTEx tissue") + 
  theme_light() +
  custom_ggtheme  +
  theme(legend.box = "horizontal",
        legend.box.margin=margin(b = -11,l = -20)) +
  guides(colour = guide_legend(override.aes = list(fill = '#cccccc'),
                               #title = "log10 q",
                               label.position = "bottom",
                               order = 2))

file_name <- paste0(getwd(), "/results/_paper/figures/nmd_age")
ggplot2::ggsave(paste0(file_name, ".png"), width = 110, height = 90, units = "mm", dpi = 300)
