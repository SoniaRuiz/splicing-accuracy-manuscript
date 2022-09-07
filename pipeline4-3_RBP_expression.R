library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)

## source("/home/sruiz/PROJECTS/splicing-project/pipeline3_RBP_expression.R")
## source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline4-2_age_stratification.R")

project_id <- "BRAIN"

## Load reference gtf and get only the genes
ensembl105 <- biomaRt::useEnsembl(biomart = 'genes', 
                                  dataset = 'hsapiens_gene_ensembl',
                                  version = 105)


## Load the RBPs and add the ensemblID
all_RBPs <- xlsx::read.xlsx(file = '/home/sruiz/PROJECTS/splicing-project-recount3/markdowns/41586_2020_2077_MOESM3_ESM.xlsx', 
                            sep = '\t', header = TRUE,
                            sheetIndex = 1) %>%  as_tibble() 


################################
## FUNCTIONS       
################################

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

tidy_sample_metadata <- function(rse) {
  
  sample_metadata <- rse %>% 
    SummarizedExperiment::colData() %>%
    as_tibble() %>%
    filter(gtex.smafrze != "EXCLUDE")
  
  age_numeric <- as.numeric(factor(as.matrix(sample_metadata$gtex.age))) 
  sample_metadata$gtex.age <- age_numeric
  
  # sample_metadata <- sample_metadata %>%
  #   mutate(gtex.age = ifelse(gtex.age == 1 | gtex.age == 2, 1, gtex.age)) %>%
  #   mutate(gtex.age = ifelse(gtex.age == 3 | gtex.age == 4, 2, gtex.age)) %>%
  #   mutate(gtex.age = ifelse(gtex.age == 5 | gtex.age == 6, 3, gtex.age))
  # sample_metadata$gtex.age
  
  
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

RBP_corrected_TPM_analysis <- function(project_id) {
  
  # project_id <- "BRAIN"
  # project_id <- "BLOOD"
  # project_id <- "MUSCLE"
  
  ## LOAD TPM CORRECTED VALUES
  tpm_corrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                          project_id, "/results/pipeline3/rbp/tpm_residuals.csv"), header = T)
  
  
  ## INIT AGE SUPERGROUPS
  age_samples_clusters_tidy <- age_stratification_init_data(project_id)
  
  
  
  ## FIlter TPM by age group and merge
  tpm_20_39 <- tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "20-39") %>% 
                     pull(individual))) %>%
    gather(key = "RBP", value = "TPM", -X ) %>%
    mutate(age_text = paste("20-39_", X)) %>%
    mutate(age = "20-39") %>% 
    as_tibble()
  
  
  tpm_40_59 <-tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "40-59") %>% 
                     pull(individual)))  %>%
    gather(key = "RBP", value = "TPM", -X ) %>%
    mutate(age_text = paste("40-59_", X)) %>%
    mutate(age = "40-59") %>% as_tibble()
  
  
  tpm_60_79 <- tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "60-79") %>% 
                     pull(individual))) %>%
    gather(key = "RBP", value = "TPM", -X ) %>%
    mutate(age_text = paste("60-79_", X)) %>%
    mutate(age = "60-79") %>% as_tibble()
  
  
  gattered_matrix <- do.call("rbind", list(tpm_20_39, tpm_40_59, tpm_60_79)) %>%
    dplyr::rename(sample = X) %>% as_tibble() %>%
    group_by(age, RBP) %>%
    mutate(mean_tpm = (TPM %>% mean())) %>%
    distinct(mean_tpm, .keep_all = T)  %>%
    ungroup() %>%
    arrange(age , mean_tpm) %>%
    mutate(RBP = fct_inorder(RBP))
  
  
  ## PLOT HEATMAP
  ggplot(data = gattered_matrix, 
         aes(x = age, y = RBP, fill = mean_tpm)) +
    geom_tile() + 
    scale_fill_gradient(low = "green", high = "red") +
    ggtitle(paste0("Mean RBP level of expression across samples\nfrom each age cluster - ", project_id
                   ,".\nTPM values have been covariate corrected.")) +
    ylab("RBP") +
    xlab("age cluster") +
    theme(axis.text.y = element_blank()) 
  
  ggsave(filename = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                           "/results/pipeline3/rbp/heatmap.png"))
  
  
  
  
  ## GET RBPs decreasing
  spread_matrix <- gattered_matrix %>%
    dplyr::select(RBP, age, mean_tpm) %>%
    spread(key = age, value = mean_tpm)
  spread_matrix
  write_csv(x = spread_matrix,
            file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/results/pipeline3/rbp/tpm_age_spread.csv"))
  
}

RBP_uncorrected_TPM_lm <- function(projects_id = c("BRAIN")) {
  
  
  ################################
  ## LOAD RBPs OF INTEREST
  ################################
  
  if (!exists("all_RBPs")) {
    all_RBPs <- xlsx::read.xlsx(file = '/home/sruiz/PROJECTS/splicing-project-recount3/markdowns/41586_2020_2077_MOESM3_ESM.xlsx', 
                                sep = '\t', header = TRUE,
                                sheetIndex = 1) %>%  as_tibble() 
  }
  
  
  all_RBPs_tidy <- all_RBPs %>%
    mutate(type = ifelse(Splicing.regulation == 1, "splicing regulation", "other")) %>%
    dplyr::select(name, ensembl_gene_id = id, type) %>%
    distinct(name, .keep_all = T) 
  
  
  
  for (project_id in projects_id) {
    ## project_id <- "BRAIN"
    
    ## Load and tidy the uncorrected TPMs
    ## TPMs here should be uncorrected as the covariates are going to be included in the linear models
    tpm_uncorrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                              project_id, "/results/pipeline3/rbp/tpm_ensembl105.csv"))
    
    
    tpm_uncorrected <- tpm_uncorrected %>%
      dplyr::rename(gene = "X") %>%
      as_tibble() %>%
      distinct(gene, .keep_all = T)
    
    
    
    ## Load and tidy the sample metadata
    sample_metadata <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                              project_id, "/results/pipeline3/rbp/covariates.csv"), header = T,
                                fileEncoding = "UTF-8") %>% as_tibble()
    
    print(sample_metadata %>% nrow)
    ## Check the samples are the same
    if ( !( identical(x = names(tpm_uncorrected)[-1],
                      y = names(sample_metadata)[-1]) ) ) {
      print("ERROR!")
      break;
    }
    
    
    
    RBPs <- (tpm_uncorrected$gene)
    
    # tpm_uncorrected[rowSums(tpm_uncorrected[, c(2:ncol(tpm_uncorrected))]),]
    
    ## Obtain values and run lm
    df_lm_output <- map_df(RBPs, function(RBP) {
      
      # RBP <- RBPs[1]
      # RBP <- "ENSG00000277564"
      print(RBP)
      
      tpm <- tpm_uncorrected %>%
        filter(gene == RBP) %>%
        mutate(gene = "tpm") %>%
        as_tibble()
      
      if (rowSums(tpm[, c(2:ncol(tpm))]) != 0) {
        
        df <- rbind(tpm %>% 
                      dplyr::rename(covariates = "gene"),
                    sample_metadata) %>%
          gather(sample, value, -covariates)  %>%
          spread(key = covariates, value) %>%
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
        
        if (project_id == "MUSCLE") {
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
                            #gtex.smtsd +
                            gtex.age,
                          data = df)
        } else {
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
        }
        
        
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
    
    
    # any(df_lm_output$pval > 0.05)
    
    
    ## Filter by age covariate and order
    df_lm_output_age <- df_lm_output %>%
      filter(str_detect(string = covariate, pattern = "gtex.age")) %>%
      arrange(Estimate)
    
    
    if (exists("RBPs_annotated")) {

      RBPs_annotated_tidy <- left_join(x = RBPs_annotated %>% as_tibble(), 
                                       y = all_RBPs_tidy %>% as_tibble(), 
                                       by = c("ensembl_gene_id" = "ensembl_gene_id"))
    } else {
      RBPs_annotated_tidy <- all_RBPs_tidy
    }
    
    ## Add gene SYMBOL info
    df_lm_age_tidy <- left_join(x = df_lm_output_age, 
                                y = RBPs_annotated_tidy %>% drop_na(), 
                                by = c("RBP_ID" = "ensembl_gene_id")) %>%
      group_by(RBP_ID) #%>%
    #distinct(pval, .keep_all = T)
    
    
    write_csv(x = df_lm_age_tidy,
              file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                            "/results/pipeline3/rbp/tpm_lm_all_ensembl105.csv"))
    
  }
  
  ###############################################
  ## Get RBPs that decrease expression with age
  ###############################################
  
  
  # RBP_decrease <- df_lm_age_tidy %>%
  #   filter(Estimate < 0) %>%
  #   arrange(desc(abs(Estimate)))
  # 
  # RBP_increase <- df_lm_age_tidy %>%
  #   filter(Estimate > 0) %>%
  #   arrange(desc(Estimate))
  # 
  # intersect(RBP_decrease$hgnc_symbol,
  #           RBP_increase$hgnc_symbol)
  # 
  # RBP_decrease$hgnc_symbol %>% unique() %>% sort()
  
  
}







################################################################################
# 0. Prepare the necessary data ------------------------------------------------
################################################################################


# RBPs <- read.table(file = '/home/sruiz/PROJECTS/splicing-project-recount3/markdowns/experiment_report_2022_3_28_16h_50m.tsv', sep = '\t', header = TRUE) %>%
#   as_tibble() %>%
#   distinct(Target.gene.symbol, .keep_all = T)

RBPs_annotated <- biomaRt::getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = all_RBPs$name %>% unique,
  mart = ensembl105
)
RBPs_annotated <- RBPs_annotated %>%
  distinct(hgnc_symbol, .keep_all = T)


## Load the recount3 data
rse <- recount3::create_rse_manual(
  project = project_id,
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene")

## Transform counts
SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
## See that now we have two assayNames()
SummarizedExperiment::assayNames(rse)


################################################################################
# 1. Get the RBP TPM expression using GTEX data --------------------------------
################################################################################

dds_tpm <- get_GTEx_gene_expression(rse = rse, ensembl = ensembl105)
dds_tpm %>% head()

dds_tpm <- dds_tpm %>%
  filter(tpm > 0)

# Tidy the 'dds_tpm' object
dds_tpm_pv <- dds_tpm %>%
  dplyr::filter(gene %in% RBPs_annotated$ensembl_gene_id) %>%
  dplyr::select(gene, sample, tpm) %>% 
  group_by(gene) %>%
  distinct(sample, .keep_all = T) %>%
  pivot_wider(names_from = sample, values_from = tpm, values_fill = 0)

write.csv(x = dds_tpm_pv %>%
            tibble::column_to_rownames("gene"),
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                        "/results/pipeline3/rbp/tpm_ensembl105.csv"))

################################################################################
# 2. Get the covariates to correct for -----------------------------------------
################################################################################

sample_metadata <- tidy_sample_metadata(rse)


write_csv(x = sample_metadata %>%
            as_tibble(rownames = "covariates"), 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/covariates.csv"))

write_csv(x = sample_metadata %>%
            as_tibble(rownames = "covariates") %>%
            filter(covariates != "gtex.smtsd"), 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/all_covariates_correct_for.csv"))

write_csv(x = sample_metadata %>%
            as_tibble(rownames = "covariates") %>%
            filter(!(covariates %in% c("gtex.age","gtex.smtsd"))), 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/except_age_covariates_correct_for.csv"))




# ## Save only batch covariates
# 
# sample_metadata_batch <- sample_metadata %>%
#   as_tibble(rownames = "covariates") %>%
#   dplyr::filter(covariates %in% c("gtex.smcenter",
#                                   "gtex.smgebtch", "gtex.smgebtchd",
#                                   "gtex.smnabtch", "gtex.smnabtchd",
#                                   "gtex.smnabtcht", "gtex.dthhrdy",
#                                   "gtex.sex", "gtex.smrin"))
#   
# 
# write_csv(x = sample_metadata_batch, 
#           file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
#                         project_id, "/results/pipeline3/rbp/batch_covariates.csv"))
# 
# 
# 
# 
# 
# 
# ##





















################################################################################
################################################################################
################################################################################

# 2. Get the best covariates to correct for ------------------------------------

# ## 2. input the 'dds_tpm_pv' to the PCA function
# 
# pca <- prcomp(x = dds_tpm_pv %>%
#                 tibble::column_to_rownames(var = "gene"), 
#               center = TRUE, scale. = T)
# 
# 
# # 3. Select PCA1 - ~20 models
# pca_20models <- pca$rotation[,1:20]
# 
# 
# # 4. Tidy the metadata
sample_metadata <- tidy_sample_metadata(rse)


# # 5. Foreach column in 'sample_metadata_tidy' do cor.test with each PCA 1 to 20 (cor function)
# cor_matrix <- cor(x = sample_metadata,
#                   y = pca_20models)
# 
# saveRDS(object = cor_matrix,
#         file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                       "/results/pipeline3/rbp/correlation_matrix.rds"))
# 
# # 6. Point above will return a matrix of correlation metrics
# # (highest correlation values & p-values are the covariates influencing the most)
# 
# gattered_cor_matrix <- reshape2::melt(cor_matrix)
# head(gattered_cor_matrix)
# 
# ggplot(data = gattered_cor_matrix, aes(x=Var1, y=Var2, fill = value)) +
#   geom_tile()+ 
#   theme_minimal() +
#   scale_fill_gradient2(
#     name = "Scaled Value",
#     low = "darkblue",
#     mid = "gray",
#     high = "darkred"
#   ) +
#   theme(
#     axis.text.x = element_text(
#       angle = 90,
#       hjust = 1,
#       vjust = 0.5
#     ))
# 
# # 7. Order covariates by PC 2
# covariates <- gattered_cor_matrix %>% 
#   as_tibble() %>%
#   #filter(Var2 == "PC2") %>%
#   mutate(value = abs(value)) %>%
#   arrange(desc(value)) %>%
#   filter(value > 0.25) %>%
#   #filter(Var1 != "gtex.age") %>%
#   pull(Var1)


# 8. Save the covariates ---------------------------


write_csv(x = sample_metadata %>%
            dplyr::rename(gtex.age = gtex.age.num), 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/covariates.csv"))






########################################################
## 9. Save only batch covariates ---------------------

sample_metadata_batch <- sample_metadata %>%
  dplyr::select(all_of(c("sample", "gtex.smcenter",
                         "gtex.smgebtch", "gtex.smgebtchd",
                         "gtex.smnabtch", "gtex.smnabtchd",
                         "gtex.smnabtcht", "gtex.dthhrdy")))


write_csv(x = t(sample_metadata_batch), 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/batch_covariates.csv"))


###############################
# 9. Call Aine's function
###############################




# # 2. Run the Principal Component Analysis (PCA) --------------------------------
# 
# ## The idea of PCA is simple — reduce the number of variables of a data set, while preserving as much information as possible.
# 
# pc <- run_PCA_analysis(rse)
# 
# ## Visualize eigenvalues (scree plot). 
# ## Show the percentage of variances explained by each principal component.
# factoextra::fviz_eig(pc)
# 
# 
# ## Graph of variables. 
# ## Positive correlated variables point to the same side of the plot. 
# ## Negative correlated variables point to opposite sides of the graph.
# factoextra::fviz_pca_var(pc,
#                          col.var = "contrib", # Color by contributions to the PC
#                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                          repel = TRUE     # Avoid text overlapping
# )
# 
# #compute variance
# pr_var <- (pc$sdev)^2
# #proportion of variance explained
# prop_varex <- pr_var/sum(pr_var)
# prop_varex[1:20] ## This shows that first principal component explains 31.5% variance
# 
# 
# ## This plot shows that ~ 20 components explains around 99% variance in the data set.
# plot(prop_varex, xlab = "Principal Component",
#      ylab = "Proportion of Variance Explained",
#      type = "b")
# 
# 
# #cumulative scree plot
# plot(cumsum(prop_varex), xlab = "Principal Component",
#        ylab = "Cumulative Proportion of Variance Explained",
#        type = "b")
# 
# # Therefore, in this case, we’ll select number of components as 20 [PC1 to PC20] and proceed to the modeling stage. 
# # This completes the steps to implement PCA on train data. 
# # For modeling, we’ll use these 20 components as predictor variables and follow the normal procedures.


# RBP_corrected_TPM_lm <- function(project_id) {
#   
#   
#   ## Load and tidy the uncorrected TPMs
#   ## TPMs here should be uncorrected as the covariates are going to be included in the linear models
#   tpm_uncorrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
#                                             project_id, "/results/pipeline3/rbp/tpm_residuals.csv"))
#   
#   
#   tpm_uncorrected <- tpm_uncorrected %>%
#     as_tibble()
#   
#   tpm_uncorrected_gather <- gather(data = tpm_uncorrected,
#                                    key = gene,
#                                    value = tpm, -X) %>% as_tibble()
#   tpm_uncorrected <- tpm_uncorrected_gather %>%
#     tidyr::spread(key = X, value = tpm) %>% as_tibble()
#   
#   
#   
#   ## Load and tidy the sample metadata
#   sample_metadata <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
#                                             project_id, "/results/pipeline3/rbp/covariates.csv"), header = T, check.names = F,
#                               fileEncoding = "UTF-8") %>% as_tibble()
#   
#   
#   sample_metadata <- sample_metadata %>%
#     dplyr::select("covariates",colnames(tpm_uncorrected)[-1])
#   
#   
#   
#   ## Check the samples are the same
#   identical(x = names(tpm_uncorrected)[-1] %>% sort(),
#             y = names(sample_metadata)[-1] %>% sort())
#   
#   
#   
#   RBPs <- (tpm_uncorrected$gene)
#   
#   # tpm_uncorrected[rowSums(tpm_uncorrected[, c(2:ncol(tpm_uncorrected))]),]
#   
#   ## Obtain values and run lm
#   df_lm_output <- map_df(RBPs, function(RBP) {
#     
#     # RBP <- RBPs[1]
#     
#     print(RBP)
#     
#     tpm <- tpm_uncorrected %>%
#       filter(gene == RBP) %>%
#       mutate(gene = "tpm")
#     
#     if (rowSums(tpm[, c(2:ncol(tpm))]) != 0) {
#       
#       df <- rbind(tpm %>% 
#                     dplyr::rename(covariates = "gene"),
#                   sample_metadata) %>%
#         gather(sample, value, -covariates)  %>%
#         spread(key = covariates, value) %>%
#         dplyr::mutate(gtex.age = gtex.age %>% as.double(),
#                       gtex.dthhrdy = gtex.dthhrdy %>% as.double(),
#                       gtex.sex = gtex.sex %>% as.double(),
#                       gtex.smcenter = gtex.smcenter %>% as.double(),
#                       gtex.smgebtch = gtex.smgebtch %>% as.double(),
#                       gtex.smgebtchd = gtex.smgebtchd %>% as.double(),
#                       gtex.smnabtch = gtex.smnabtch %>% as.double(),
#                       gtex.smnabtchd = gtex.smnabtchd %>% as.double(),
#                       gtex.smnabtcht = gtex.smnabtcht %>% as.double(),
#                       gtex.smrin = gtex.smrin %>% as.double(),
#                       tpm = tpm %>% as.double()) %>%
#         as_tibble()
#       
#       if (project_id != "MUSCLE") {
#         # print(df)
#         lm_output <- lm(tpm ~ 
#                           #gtex.smrin + 
#                           #gtex.smcenter +
#                           #gtex.smgebtch +
#                           #gtex.smgebtchd +
#                           #gtex.smnabtch +
#                           #gtex.smnabtchd +
#                           #gtex.smnabtcht +
#                           #gtex.dthhrdy +
#                           #gtex.sex +
#                           #gtex.smtsd +
#                           gtex.age,
#                         data = df)
#       } else {
#         # print(df)
#         lm_output <- lm(tpm ~ 
#                           #gtex.smrin + 
#                           #gtex.smcenter +
#                           #gtex.smgebtch +
#                           #gtex.smgebtchd +
#                           #gtex.smnabtch +
#                           #gtex.smnabtchd +
#                           #gtex.smnabtcht +
#                           #gtex.dthhrdy +
#                           #gtex.sex +
#                           gtex.smtsd +
#                           gtex.age,
#                         data = df)
#       }
#       
#       
#       lm_output <- lm_output %>% summary()
#       df_lm_output <- lm_output$coefficients %>% as_tibble(rownames = "covariate")
#       
#       # cor(df)
#       df_lm_output <- df_lm_output %>%
#         mutate(RBP_ID = RBP) %>%
#         dplyr::select(covariate, Estimate, pval = `Pr(>|t|)`, RBP_ID)
#       
#       return(df_lm_output)
#     }
#     
#   }) 
#   
#   
#   # any(df_lm_output$pval > 0.05)
#   
#   
#   ## Filter by age covariate and order
#   df_lm_output_age <- df_lm_output %>%
#     filter( pval <= 0.05,
#             str_detect(string = covariate, pattern = "gtex.age")) %>%
#     arrange(Estimate)
#   
#   
#   ## Add gene SYMBOL info
#   df_lm_age_tidy <- left_join(x = df_lm_output_age, 
#                               y = RBPs_tidy, 
#                               by = c("RBP_ID" = "ensembl_gene_id")) %>%
#     group_by(RBP_ID) %>%
#     distinct(pval, .keep_all = T)
#   
#   
#   write_csv(x = df_lm_age_tidy,
#             file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                           "/results/pipeline3/rbp/tpm_lm_corrected.csv"))
#   
#   ###############################################
#   ## Get RBPs that decrease expression with age
#   ###############################################
#   
#   
#   RBP_decrease <- df_lm_age_tidy %>%
#     filter(Estimate < 0) %>%
#     arrange(desc(abs(Estimate)))
#   
#   RBP_increase <- df_lm_age_tidy %>%
#     filter(Estimate > 0) %>%
#     arrange(desc(Estimate))
#   
#   intersect(RBP_decrease$hgnc_symbol,
#             RBP_increase$hgnc_symbol)
#   
#   RBP_decrease$hgnc_symbol %>% unique() %>% sort()
#   RBP_increase$hgnc_symbol %>% unique() %>% sort()
#   
#   
# }
# 
# RBP_analysis_comparison <- function(project_id) {
#   
#   
#   df_analysis_corrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                                                   "/results/pipeline3/rbp/tpm_age_spread.csv"), check.names = F) %>% as_tibble()
#   
#   df_analysis_lm_uncorrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                                                        "/results/pipeline3/rbp/tpm_lm.csv"))
#   
#   df_analysis_lm_corrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                                                      "/results/pipeline3/rbp/tpm_lm_corrected.csv"))
#   
#   
#   
#   RBP_corrected_decreasing <- df_analysis_corrected %>%
#     filter(`20-39` > `40-59`, `20-39` > `60-79`) %>%
#     pull(RBP)
#   
#   RBP_uncorrected_decreasing <- df_analysis_lm_uncorrected %>%
#     filter(Estimate < 0) %>%
#     pull(RBP_ID)
#   
#   
#   RBPs <- biomaRt::getBM(
#     attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = intersect(RBP_corrected_decreasing,
#                        RBP_uncorrected_decreasing),
#     mart = ensembl106
#   )
#   
#   RBPs$hgnc_symbol %>% sort()
#   
#   ##
#   
#   RBP_lm_corrected_decreasing <- df_analysis_lm_corrected %>%
#     filter(Estimate < 0)
#   
#   RBP_lm_corrected_decreasing$hgnc_symbol %>% sort()
#   
#   intersect(RBPs$hgnc_symbol %>% sort(),
#             RBP_lm_corrected_decreasing$hgnc_symbol)
# }
# 
# RBP_lm_analysis_tissue_comparison <- function() {
#   
#   projects_id <- c("BRAIN", "MUSCLE", "BLOOD")
#   
#   df_tissues <- map_df(projects_id, function(project_id) {
#     
#     # project_id <- projects_id[1]
#     
#     df_analysis_lm_uncorrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                                                          "/results/pipeline3/rbp/tpm_lm.csv")) %>%
#       mutate(tissue = project_id)
#     
#     return(df_analysis_lm_uncorrected)
#     
#   })
#   
#   df_tissues %>%
#     dplyr::select(Estimate, hgnc_symbol, tissue) %>%
#     spread(key = tissue, value = Estimate) %>%
#     dplyr::rename(RBP = hgnc_symbol) %>%
#     DT::datatable()
#   
#   
# }

