library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(DBI)


setwd(normalizePath("."))

source(paste0(getwd(), "/database_SQL_helper.R"))
dependencies_folder <- paste0(getwd(), "/dependencies/")

gtf_version <- 105
main_project <- "splicing"


# ## Load reference GTF and get only the genes
ensembl105 <- rtracklayer::import(con = paste0(dependencies_folder, "/Homo_sapiens.GRCh38.105.chr.gtf")) 


################################
## FUNCTIONS       
################################

tidy_sample_metadata <- function(sample_metadata, samples) {
  
 
  
  age_numeric <- as.numeric(factor(as.matrix(sample_metadata$gtex.age))) 
  sample_metadata$gtex.age <- age_numeric
 
  sample_metadata <- sample_metadata %>%
    tibble::column_to_rownames(var = "external_id")
  
  covariates <- c("gtex.age", "gtex.smcenter", "gtex.smtsd",
                  "gtex.smgebtch", "gtex.smgebtchd",
                  "gtex.smnabtch", "gtex.smnabtchd", "gtex.smnabtcht",
                  "gtex.dthhrdy", "gtex.sex", "gtex.smrin")
  
  
  smcenter <- as.numeric(factor(as.matrix(sample_metadata$gtex.smcenter))) 
  sample_metadata$gtex.smcenter <- smcenter
  
  gtex.smgebtch <- as.numeric(factor(as.matrix(sample_metadata$gtex.smgebtch))) 
  sample_metadata$gtex.smgebtch <- gtex.smgebtch
  
  gtex.smgebtchd <- as.numeric(factor(as.matrix(sample_metadata$gtex.smgebtchd))) 
  sample_metadata$gtex.smgebtchd <- gtex.smgebtchd
  
  gtex.smnabtch <- as.numeric(factor(as.matrix(sample_metadata$gtex.smnabtch))) 
  sample_metadata$gtex.smnabtch <- gtex.smnabtch
  
  gtex.smnabtchd <- as.numeric(factor(as.matrix(sample_metadata$gtex.smnabtchd))) 
  sample_metadata$gtex.smnabtchd <- gtex.smnabtchd
  
  gtex.smnabtcht <- as.numeric(factor(as.matrix(sample_metadata$gtex.smnabtcht))) 
  sample_metadata$gtex.smnabtcht <- gtex.smnabtcht
  
  gtex.smtsd <- as.numeric(factor(as.matrix(sample_metadata$gtex.smtsd))) 
  sample_metadata$gtex.smtsd <- gtex.smtsd
  
  
  # 8. Return covariates ---------------------------
  
  return(t(sample_metadata %>%
           dplyr::select(all_of(covariates))))
}

RBP_uncorrected_TPM_lm <- function(project_id = "BRAIN",
                                   tpm_uncorrected = dds_tpm,
                                   sample_metadata) {
  
  
  ################################
  ## LOAD RBPs OF INTEREST
  ################################
  
  all_RBPs <- xlsx::read.xlsx(file = paste0(dependencies_folder, '/RBPs_subgroups.xlsx'), 
                              header = TRUE,
                              sheetIndex = 1) %>% 
    as_tibble() %>%
    distinct(name, .keep_all = T) 
  
  
  ## Only keep the TPM value of the RBPs of interest
  tpm_uncorrected <- tpm_uncorrected %>%
    filter(gene %in% all_RBPs$id)
 
  
  ## Check the samples used are the same
  if ( !( identical(x = colnames(tpm_uncorrected)[-1] %>% unique %>% sort(),
                    y = colnames(sample_metadata)[-1] %>% unique %>% sort()) ) ) {
    print("ERROR!")
    break;
  }
  
  
  RBPs <- (tpm_uncorrected$gene) %>% unique
    

  ## Per each RBP, obtain TPM and metadata values and run lm
  df_lm_output <- map_df(RBPs, function(RBP) {
    
    # RBP <- RBPs[1]
    # RBP <- "ENSG00000002079"
    
    print(RBP)
    
    tpm <- tpm_uncorrected %>%
      filter(gene == RBP) %>%
      mutate(gene = "tpm") %>%
      as_tibble()
    
    if ( rowSums(tpm[, c(2:ncol(tpm))], na.rm = T) != 0 ) {

      
      ## Transform TPM to log10 so resilduals after modelling in the linear model are normally distributed
      ## https://stats.stackexchange.com/questions/40907/in-regression-analysis-what-does-taking-the-log-of-a-variable-do
      df_rbp <- rbind(tpm %>%  dplyr::rename(covariates = "gene"),
                  sample_metadata) %>%
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
                        data = df_rbp)
      
      
      
      lm_output <- lm_output %>% summary()
      df_lm_output <- lm_output$coefficients %>% as_tibble(rownames = "covariate")

      
      df_lm_output <- df_lm_output %>%
        mutate(RBP_ID = RBP) %>%
        dplyr::select(covariate, Estimate, pval = `Pr(>|t|)`, RBP_ID)
      
      
      plot(density(lm_output$residuals)) 
      
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
                              y = all_RBPs %>% drop_na(), 
                              by = c("RBP_ID" = "id")) %>%
    group_by(RBP_ID) 
  
 
  
  ## Save results
  write_csv(x = df_lm_output_tidy,
            file = paste0(getwd(), "/results/_paper/results/",project_id, "_",main_project,"_tpm_lm_all_ensembl105.csv"))
  
  
  ###############################################
  ## Get RBPs that decrease expression with age
  ###############################################
  
  ## RBPs with expression levels decreasing with age (i.e. with negative effect over their TPM associated to the 'age' covariate)
  RBP_decrease <- df_lm_output_tidy %>%
    filter(str_detect(string = covariate, pattern = "gtex.age")) %>%
    #filter(Estimate < 0) %>%
    filter(q <= 0.05) %>%
    arrange(name)
  
  RBP_decrease %>% print(n=50)
  RBP_decrease %>% filter(q <= 0.05) %>% print(n=50)
  RBP_decrease %>% filter(q <= 0.05) %>% pull(q) %>% summary()
  
  ## Count them as a percentage
  ((RBP_decrease %>% filter(q <= 0.05) %>%  nrow) * 100) / (all_RBPs %>% nrow)
  
  RBP_decrease$q %>% sort()
  
  ## RBPs with expression levels increasing with age (i.e. with positive effect over their TPM associated to the 'age' covariate)
  RBP_increase <- df_lm_output_tidy %>%
    filter(str_detect(string = covariate, pattern = "gtex.age")) %>%
    filter(Estimate > 0) %>%
    arrange(name)
  
  #RBP_increase %>% print(n=50)
  RBP_increase %>% print(n=50)
  RBP_increase %>% filter(q <= 0.05) %>% print(n=50)
  
  ## Count them as a percentage
  ((RBP_increase %>% filter(q <= 0.05) %>% nrow) * 100) / (all_RBPs %>% nrow)
  
  
  ####
  
  ((df_lm_output_tidy %>% filter(q > 0.05) %>% nrow) * 100) / (all_RBPs %>% nrow)
  
  intersect(RBP_decrease$name,
            RBP_increase$name)

  RBP_decrease$name %>% unique() %>% sort()
  
  
}



###########################################
## 1. GET RECOUNT 3 EXPRESSION DATA
##    AND METADATA FOR BRAIN TISSUE
###########################################


project_id <- "BRAIN"
message(Sys.time(), " - ", project_id)


## Get metadata for brain tissue
metadata <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", 
                                  gtf_version, "/",main_project,"/base_data/", project_id, "_samples_metadata.rds"))

metadata$gtex.smafrze %>% unique
metadata$gtex.smrin %>% unique %>% min


## Get TPM for brain tissue

recount_tpm <- generate_recount3_tpm(recount3_project_IDs = project_id,
                                     ref = ensembl105,
                                     gtf_version,
                                     main_project)

recount_tpm %>% head
recount_tpm %>% nrow()
recount_tpm %>% ncol()

## 1. Get the TPM expression using GTEX data

## Filter by the samples of the current cluster
dds_tpm <- recount_tpm %>% 
  dplyr::select(gene, all_of(metadata$external_id))


dds_tpm %>% head
recount_tpm %>% nrow()
recount_tpm %>% ncol()
  

# folder_name <- paste0(getwd(),"/results/RBP_EXPRESSION/splicing/")
# write.csv(x = dds_tpm_wsymbol,
#           file = paste0(folder_name, "/", project_id, "_tpm_ensembl105.csv"))


# 2. Get the covariates to correct for
sample_metadata <- tidy_sample_metadata(metadata, samples = metadata$external_id) %>% 
  as_tibble(rownames = "covariates")
# write_csv(x = sample_metadata %>% as_tibble(rownames = "covariates"), 
#           file = paste0(folder_name,  "/", project_id, "_covariates.csv"))


##################
# 1. CALLS 
##################

RBP_uncorrected_TPM_lm(project_id = project_id,
                       tpm_uncorrected = dds_tpm,
                       sample_metadata = sample_metadata)