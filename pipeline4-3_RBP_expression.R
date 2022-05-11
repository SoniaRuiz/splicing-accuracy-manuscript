library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(biomaRt)

## source("/home/sruiz/PROJECTS/splicing-project/pipeline3_RBP_expression.R")

################################
## FUNCTIONS       
################################

get_GTEx_gene_expression <- function(rse, ensembl106, recount3 = T) {
  
  
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
  
  assays(rse)$TPM <- recount::getTPM(rse)
  ## Should all be equal to 1
  colSums(assay(rse, "TPM")) / 1e6
  
  ## Tidy the dataframe
  recount_dds_tpm <- assays(rse)$TPM %>%
    as_tibble(rownames = "gene") %>% 
    mutate(gene = gene %>% str_remove("\\..*")) %>% 
    tidyr::gather(key = sample, value = tpm, -gene)
  
  # recount_dds_tpm <- recount_dds_tpm %>% dplyr::arrange(gene)
  
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
  
  
  ## 2. Compare TPM values
  # dds_tpm %>% filter(gene == "ENSG00000160201")
  
  recount_dds_tpm %>%
    filter(gene == "ENSG00000160201")
  
  # df_tpm_merged <- merge(x = dds_tpm %>% data.table::as.data.table(),
  #                        y = recount_dds_tpm %>% data.table::as.data.table(),
  #                        by = "gene") 
  # 
  # 
  # dds_tpm %>% distinct(gene) %>% nrow()
  # recount_dds_tpm %>% distinct(gene) %>% nrow()
  
  

  
  
  ## Return recount3 TPM counts
  
  
  # listAttributes(ensembl106)
  # searchAttributes(mart = ensembl106, pattern = "hgnc_symbol")
  mapping <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = recount_dds_tpm$gene %>% unique,
    mart = ensembl106
  )
  
  df_return <- recount_dds_tpm %>% 
    mutate(project = project_id) %>%
    dplyr::left_join(mapping,
                     by =  c("gene" = "ensembl_gene_id"))
  

  return(df_return)
  
}

tidy_sample_metadata <- function(rse) {
  
  sample_metadata <- rse %>% 
    SummarizedExperiment::colData() %>%
    as_tibble() %>%
    filter(gtex.smafrze != "EXCLUDE")
  
  age_numeric <- as.numeric(factor(as.matrix(sample_metadata$gtex.age))) 
  sample_metadata$gtex.age <- age_numeric
  
  sample_metadata <- sample_metadata %>%
    mutate(gtex.age = ifelse(gtex.age == 1 | gtex.age == 2, 1, gtex.age)) %>%
    mutate(gtex.age = ifelse(gtex.age == 3 | gtex.age == 4, 2, gtex.age)) %>%
    mutate(gtex.age = ifelse(gtex.age == 5 | gtex.age == 6, 3, gtex.age))
  sample_metadata$gtex.age
  
  
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
  
  
  # 8. Return covariates ---------------------------
  
  return(t(sample_metadata %>%
           dplyr::select(all_of(covariates))))
}

################################
## CALLS       
################################

## RBP EXPRESSION ANALYSIS

################################################################################
# 0. Prepare the necessary data ------------------------------------------------
################################################################################

## Load reference gtf and get only the genes
ensembl106 <- biomaRt::useEnsembl(biomart = 'genes', 
                                  dataset = 'hsapiens_gene_ensembl',
                                  version = 106)


## Load the RBPs and add the ensemblID
RBPs <- read.table(file = 'markdowns/experiment_report_2022_3_28_16h_50m.tsv', sep = '\t', header = TRUE) %>%
  as_tibble() %>%
  distinct(Target.gene.symbol)

RBPs_tidy <- biomaRt::getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = RBPs$Target.gene.symbol %>% unique,
  mart = ensembl106
)

## Load the recount3 data
project_id <- "BRAIN"
rse <- recount3::create_rse_manual(
  project = project_id,
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene")

assays(rse)$counts <- recount3::transform_counts(rse)

## See that now we have two assayNames()
rse
assayNames(rse)





################################################################################
# 1. Get the RBP TPM expression using GTEX data --------------------------------
################################################################################

dds_tpm <- get_GTEx_gene_expression(rse, ensembl106)
dds_tpm %>% head()

# Tidy the 'dds_tpm' object
dds_tpm_pv <- dds_tpm %>%
  dplyr::filter(gene %in% RBPs_tidy$ensembl_gene_id) %>%
  dplyr::select(gene, sample, tpm) %>% 
  group_by(gene) %>%
  distinct(sample, .keep_all = T) %>%
  pivot_wider(names_from = sample, values_from = tpm, values_fill = 0)

is.na(dds_tpm_pv) <- sapply(dds_tpm_pv, is.infinite)
dds_tpm_pv[is.na(dds_tpm_pv)] <- 0

write.csv(x = dds_tpm_pv %>%
            tibble::column_to_rownames("gene"),
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                        "/results/pipeline3/rbp/tpm.csv"))


# 2. Get the best covariates to correct for ------------------------------------

sample_metadata <- tidy_sample_metadata(rse)

write_csv(x = sample_metadata %>%
            as_tibble(rownames = "covariates"), 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/covariates.csv"))






########################################################
## 9. Save only batch covariates ---------------------

sample_metadata_batch <- sample_metadata %>%
  as_tibble(rownames = "covariates") %>%
  dplyr::filter(covariates %in% c("gtex.smcenter",
                                  "gtex.smgebtch", "gtex.smgebtchd",
                                  "gtex.smnabtch", "gtex.smnabtchd",
                                  "gtex.smnabtcht", "gtex.dthhrdy",
                                  "gtex.sex", "gtex.smrin"))
  

write_csv(x = sample_metadata_batch, 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/batch_covariates.csv"))






##





















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

