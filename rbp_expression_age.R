library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(DBI)

## source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline4-3_RBP_expression.R")
setwd("~/splicing-accuracy-manuscript/")


gtf_version <- 105
main_project <- "splicing"


dependencies_folder <- paste0(getwd(), "/dependencies/")

## Load reference GTF and get only the genes
ensembl105 <- rtracklayer::import(con = paste0(dependencies_folder, "/Homo_sapiens.GRCh38.105.chr.gtf")) %>% 
  as_tibble() %>% 
  dplyr::select(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = T)


################################
## FUNCTIONS       
################################

get_genes_to_analyse_expression <- function(type, ensembl105) {
  
  if (type == "RBP") {
    
    ## Load the RBPs and add the ensemblID
    genes_annotated <- xlsx::read.xlsx(file = paste0(dependencies_folder, '/RBPs_subgroups.xlsx'), 
                                header = TRUE,
                                sheetIndex = 1) %>%  
      as_tibble() %>%
      dplyr::rename(gene_id = id)
    

  
    
  } else {
    
    ## Load the NMD molecules
    all_NMD <- read.delim(file = paste0(dependencies_folder, '/NMD.txt'),
                          sep = '\t', header = TRUE) %>% as_tibble() %>% 
      dplyr::select(hgnc_symbol = Name) %>%
      distinct(hgnc_symbol)
    
    ## Tidy them
    genes_annotated <- all_NMD %>% 
      left_join( y = ensembl105,
                 by = c("hgnc_symbol" = "gene_name"))

  }
  
  return(genes_annotated)
}

get_GTEx_gene_expression <- function(rse, ensembl, recount3 = T) {
  
  
  ################################################
  ## Compute TPM values using MANUAL approach
  ################################################

  if (!recount3) {
    ## Counts have been transformed, so they are ready to work with them
    ddsSE <- rse #DESeqDataSet(rse, design = ~ gtex.smnabtcht)
    
    # Remove anything after . in ensembl id
    rownames(ddsSE) <- rownames(ddsSE) %>%
      str_remove("\\..*")
    
    # Convert to tpm, which is calculated by:
    # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
    # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
    # 3. Divide the RPK values by the “per million” scaling factor.
    dds_rpk <- ddsSE %>%
      SummarizedExperiment::assay("counts") %>%
      as_tibble(rownames = "gene") %>%
      tidyr::pivot_longer(cols = -c("gene"),
                          names_to = "recount_id",
                          values_to = "counts") %>%
      #filter(recount_id %in% sample_used) %>%
      dplyr::inner_join(ref %>%
                          as_tibble() %>%
                          dplyr::select(gene_id, gene_name, width),
                        by = c("gene" = "gene_id")) %>%
      dplyr::mutate(rpk = counts/width)
    
    dds_tpm <- dds_rpk %>%
      group_by(recount_id) %>%
      mutate(scaling_factor = sum(rpk)/1e6) %>%
      ungroup() %>%
      dplyr::mutate(tpm = rpk/scaling_factor) %>%
      dplyr::mutate(tpm_log10 = log10(tpm)) ## Normalisation log10
    
    # dds_tpm <- dds_tpm %>% dplyr::arrange(gene)
  }
  
  ################################################
  ## Compute TPM values using recount3 approach
  ################################################
  
  SummarizedExperiment::assays(rse)$TPM <- recount::getTPM(rse)
  ## Should all be equal to 1
  indx <- which(colSums(SummarizedExperiment::assay(rse, "TPM")) / 1e6 > 1) 
  
  any(SummarizedExperiment::assays(rse)$TPM[-(indx %>% unlist %>% unname),]/ 1e6 > 1)
  
  ## Tidy the dataframe
  recount_dds_tpm <- SummarizedExperiment::assays(rse)$TPM[-(indx %>% unlist %>% unname),] %>%
    as_tibble(rownames = "gene") %>% 
    mutate(gene = gene %>% str_remove("\\..*")) %>% 
    tidyr::gather(key = sample, value = tpm, -gene)
  
  recount_dds_tpm 
  

  
  ## Filter by 'USE ME' sample
  
  ################################################
  ## Compare both results
  ################################################
  
  ## 1. Filter by only "USE ME" samples
  
  sample_used <- rse %>% 
    SummarizedExperiment::colData() %>%
    as_tibble() %>%
    filter(gtex.smafrze != "EXCLUDE") %>%
    pull(external_id)
  
  if (!recount3) {
    dds_tpm <- dds_tpm %>% 
      filter(recount_id %in% sample_used)
  }

  
  recount_dds_tpm <- recount_dds_tpm %>% 
    filter(sample %in% sample_used) 

  
  
  ## 2. Add gene name
  
  mapping <- biomaRt::getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = recount_dds_tpm$gene %>% unique,
    mart = ensembl
  )
  
  df_return <- recount_dds_tpm %>% 
    mutate(project = project_id) %>%
    dplyr::left_join(mapping,
                     by =  c("gene" = "ensembl_gene_id"))
  
  ## 3. Return recount3 TPM counts
  
  return(df_return)
  
}

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

RBP_uncorrected_TPM_lm <- function(projects_used = c("BRAIN")) {
  
  
  ################################
  ## LOAD RBPs OF INTEREST
  ################################
  
  if ( !exists("all_RBPs") ) {
    all_RBPs <- xlsx::read.xlsx(file = paste0(dependencies_folder, '/RBPs_subgroups.xlsx'), 
                                header = TRUE,
                                sheetIndex = 1) %>%  as_tibble() 
  }
  
  
  all_RBPs_tidy <- all_RBPs %>%
    distinct(name, .keep_all = T) 
  

  
  for (project_id in projects_used) {
    
    ## project_id <- "BRAIN"
    
    ## Load clusters used
    if (str_detect(main_project, pattern = "age")) {
      all_clusters <- c("20-39","40-59","60-79")
    } else {
      all_clusters <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", 
                                            gtf_version, "/", main_project, "/base_data/", project_id, "_clusters_used.rds"))  
    }
    
    
    ## Load and tidy the uncorrected TPMs
    tpm_uncorrected <- map_df(all_clusters, function(cluster_id) {
      
      print(paste0(Sys.time(), " - ", cluster_id, " ... "))
      if ( file.exists(paste0(getwd(), "/results/RBP_EXPRESSION/", main_project, "/RBP/", 
                              project_id, "/", cluster_id, "/tpm_ensembl105.csv")) ) {
        read.csv(file = paste0(getwd(), "/results/RBP_EXPRESSION/", main_project, "/RBP/", 
                               project_id, "/", cluster_id, "/tpm_ensembl105.csv"), 
                 header = T, fileEncoding = "UTF-8") 
      } else {
        return(NULL)
      }   
    })
    
    tpm_uncorrected_tidy <- tpm_uncorrected %>%
      dplyr::rename(gene = "X") %>%
      as_tibble() %>%
      distinct(gene, .keep_all = T)
    
    
    
    ## Load and tidy the sample metadata
    sample_metadata <- map_df(all_clusters, function(cluster_id) {
      
      print(paste0(Sys.time(), " - ", cluster_id, " ... "))
      if (file.exists(paste0(getwd(), "/results/RBP_EXPRESSION/", main_project, "/RBP/", 
                             project_id, "/", cluster_id, "/covariates.csv"))) {
        read.csv(file = paste0(getwd(), "/results/RBP_EXPRESSION/", main_project, "/RBP/", 
                               project_id, "/", cluster_id, "/covariates.csv"), 
                 header = T, fileEncoding = "UTF-8") 
      } else {
        return(NULL)
      }   
    })
    
    
    print(sample_metadata %>% nrow)
    
    ## Check the samples are the same
    if ( !( identical(x = names(tpm_uncorrected_tidy)[-1] %>% sort(),
                      y = names(sample_metadata)[-1] %>% sort()) ) ) {
      print("ERROR!")
      break;
    }
    
    RBPs <- (tpm_uncorrected_tidy$gene)
    
    
    # tpm_uncorrected[rowSums(tpm_uncorrected[, c(2:ncol(tpm_uncorrected))]),]
    
    ## Obtain values and run lm
    df_lm_output <- map_df(RBPs, function(RBP) {
      
      # RBP <- RBPs[1]
      # RBP <- "ENSG00000277564"
      print(RBP)
      
      tpm <- tpm_uncorrected_tidy %>%
        filter(gene == RBP) %>%
        mutate(gene = "tpm") %>%
        as_tibble()
      
      if ( rowSums(tpm[, c(2:ncol(tpm))], na.rm = T) != 0 ) {
        

        
        df <- rbind(tpm %>%  dplyr::rename(covariates = "gene"),
                    sample_metadata) %>%
          gather(sample, value, -covariates)  %>%
          drop_na() %>%
          group_by(sample, covariates) %>%
          mutate(median_tpm = value %>% median) %>%
          distinct(median_tpm, .keep_all = T) %>%
          dplyr::select(-value) %>%
          spread(key = covariates, value = median_tpm) %>%
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
                        tpm = tpm %>% as.double()) %>%
          as_tibble()
        

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
                          data = df)
        
        
        
        lm_output <- lm_output %>% summary()
        df_lm_output <- lm_output$coefficients %>% as_tibble(rownames = "covariate")
        
        # cor(df)
        df_lm_output <- df_lm_output %>%
          mutate(RBP_ID = RBP) %>%
          dplyr::select(covariate, Estimate, pval = `Pr(>|t|)`, RBP_ID)
        
        return(df_lm_output)
        
      } else {
        
        return(NULL)
        
      }
      
    }) 
    
    
    ## Filter by age covariate and order by estimate size
    df_lm_output_age <- df_lm_output %>%
      filter(str_detect(string = covariate, pattern = "gtex.age")) %>%
      arrange(Estimate)
    
    
    # if ( exists("RBPs_annotated") ) {
    #   RBPs_annotated_tidy <- left_join(x = RBPs_annotated %>% as_tibble(), 
    #                                    y = all_RBPs_tidy %>% as_tibble(), 
    #                                    by = c("ensembl_gene_id" = "ensembl_gene_id"))
    # } else {
    #   RBPs_annotated_tidy <- all_RBPs_tidy
    # }
    
    ## Add gene SYMBOL info
    df_lm_age_tidy <- left_join(x = df_lm_output_age, 
                                y = all_RBPs_tidy %>% drop_na(), 
                                by = c("RBP_ID" = "id")) %>%
      group_by(RBP_ID) 
    
    
    write_csv(x = df_lm_age_tidy,
              file = paste0(getwd(), "/results/RBP_EXPRESSION/RBP/",project_id,"/", project_id, "_",main_project,"_tpm_lm_all_ensembl105.csv"))
    
  }
  
  ###############################################
  ## Get RBPs that decrease expression with age
  ###############################################
  
  
  RBP_decrease <- df_lm_age_tidy %>%
    filter(Estimate < 0) %>%
    arrange(desc(abs(Estimate)))
  RBP_decrease %>% print(n=50)
  

  RBP_increase <- df_lm_age_tidy %>%
    filter(Estimate > 0) %>%
    arrange(desc(Estimate))
  RBP_increase$name %>% sort

  intersect(RBP_decrease$name,
            RBP_increase$name)

  RBP_decrease$name %>% unique() %>% sort()
  
  
}







################################################################################
# 0. Prepare the necessary data ------------------------------------------------
################################################################################


## Only for the subsampled samples - age stratification subsampled to correct by RIN

projects_used <- readRDS(file = paste0(getwd(), "/results/all_final_projects_used.rds"))

for (project_id in projects_used) {
  
  # project_id <- projects_used[1]
  # project_id <- "BRAIN"

  
  message(Sys.time(), " - ", project_id)
  
  ## Load clusters used
  if (str_detect(main_project,pattern = "age")) {
    clusterIDs <- c("20-39","40-59","60-79")
  } else {
    clusterIDs <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/", project_id, "_clusters_used.rds"))  
  }
  
  ## Get metadata
  metadata <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", 
                                    gtf_version, "/",main_project,"/base_data/", project_id, "_samples_metadata.rds"))
 
  
  for (cluster_id in clusterIDs) {
    
    # cluster_id <-  clusterIDs[1]

    if (str_detect(main_project,pattern = "age")) {
      ## Get metadata
      samples_used <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/", project_id, "_", cluster_id, "_samples_used.rds"))
      cluster_metadata <- metadata %>%
        filter(external_id %in% samples_used)
    } else {
      cluster_metadata <- metadata %>%
        filter(gtex.smtsd == cluster_id) 
      ## Get metadata
      samples_used <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/", project_id, "_", cluster_id, "_samples_used.rds"))
    }
    
    ################################################################################
    # 1. Get the TPM expression using GTEX data --------------------------------
    ################################################################################
    
    ## 1. Filter by the samples of the current cluster
    tpm <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/tpm_all_samples.rds"))
    
    dds_tpm <- tpm %>% 
      tidyr::gather(key = sample, value = tpm, -gene) %>% 
      filter(tpm > 0, sample %in% samples_used)
    
    ## 2. Add gene name
    dds_tpm_wsymbol <- dds_tpm %>% 
      mutate(project = project_id) %>%
      dplyr::left_join(ensembl105,
                       by =  c("gene" = "gene_id"))
    
    ## 3. Get TPM data per RBP
    for (type in c("RBP", "NMD")) {
      
      # type <- c("RBP", "NMD")[2]
      
      print(paste0(Sys.time(), " - ", cluster_id, "...", type))
      
      genes <- get_genes_to_analyse_expression(type, ensembl105 = ensembl105)
      
      # Tidy the 'dds_tpm' object
      dds_tpm_pv <- dds_tpm %>%
        dplyr::filter(gene %in% genes$gene_id) %>%
        #dplyr::select(gene, sample, tpm) %>%
        filter(sample %in% cluster_metadata$external_id) %>%
        group_by(gene) %>%
        distinct(sample, .keep_all = T) %>%
        pivot_wider(names_from = sample, values_from = tpm, values_fill = 0)
      
      ## Save data
      folder_name <- paste0(getwd(), "/results/RBP_EXPRESSION/", main_project, "/", type, "/", 
                            project_id, "/", cluster_id, "/")
      dir.create(file.path(folder_name), recursive = TRUE, showWarnings = T)
      
      write.csv(x = dds_tpm_pv %>% tibble::column_to_rownames("gene"),
                file = paste0(folder_name, "/tpm_ensembl105.csv"))
      
      ################################################################################
      # 2. Get the covariates to correct for -----------------------------------------
      ################################################################################
      
      sample_metadata <- tidy_sample_metadata(cluster_metadata, samples = cluster_metadata$external_id)
      
      
      write_csv(x = sample_metadata %>% as_tibble(rownames = "covariates"), 
                file = paste0(folder_name, "/covariates.csv"))
      
      if (!identical(names(sample_metadata %>% as_tibble(rownames = "covariates"))[-1] %>% sort(), 
                     names(dds_tpm_pv %>% tibble::column_to_rownames("gene")) %>% sort())) {
        print("ERROR! Covariates and gene expression have different sample IDs associated.")
        break;
      }
      
    }
  
  }
  
  rm(rse)
  rm(dds_tpm)
  rm(dds_tpm_pv)
  rm(sample_metadata)
  gc()
  
}

