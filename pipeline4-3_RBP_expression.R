library(tidyverse)
library(GenomicRanges)

## source("/home/sruiz/PROJECTS/splicing-project/pipeline3_RBP_expression.R")

################################
## FUNCTIONS       
################################

get_GTEx_gene_expression <- function(rse, ref) {
  
  library("DESeq2")
  
  
  
  
  
  sample_use <- rse %>% 
    SummarizedExperiment::colData() %>%
    as_tibble() %>%
    filter(gtex.smafrze != "EXCLUDE") %>%
    pull(external_id)
  
  
  #dds_tpm <- map_df(project_ids, function(project_id) {
  
  
  
  ddsSE <- rse#DESeqDataSet(rse, design = ~ 1)
  
  
  # Remove anything after . in ensembl id
  rownames(ddsSE) <- rownames(ddsSE) %>%
    str_remove("\\..*")
  
  
  # Convert to tpm, which is calculated by:
  # 1. Divide the read counts by the length of each gene in kilobases (i.e. RPK)
  # 2. Count up all the RPK values in a sample and divide this number by 1,000,000.
  # 3. Divide the RPK values by the “per million” scaling factor.
  dds_rpk <- ddsSE %>%
    assay() %>%
    as_tibble(rownames = "gene") %>%
    tidyr::pivot_longer(cols = -c("gene"),
                        names_to = "recount_id",
                        values_to = "counts") %>%
    filter(recount_id %in% sample_use) %>%
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
  
  ## Return df
  return(dds_tpm %>% 
           mutate(project = project_id))
  
  # ggplot2::ggplot(data = dds_tpm %>%
  #                   filter(gene == "ENSG00000063244")) +
  #   geom_density(mapping = aes(x = tpm_log10))
  #                 
  # 
  # 
  # ## Add age supercluster info
  # dds_tpm <- merge(x = dds_tpm %>% data.table::as.data.table(),
  #                  y = age_samples_clusters_tidy %>% 
  #                    select(age_group, run) %>% 
  #                    data.table::as.data.table(),
  #                  by.x = "recount_id",
  #                  by.y = "run",
  #                  all.x = T)
  # 
  
  #})
  
  # ## Calculate median values across samples  
  # tpm <- dds_tpm %>%
  #   select(gene, gene_name, age_group, tpm) %>%
  #   #drop_na() %>%
  #   dplyr::group_by(gene, age_group) %>%
  #   mutate(tpm_median = tpm %>% median()) %>%
  #   mutate(tpm_mean = tpm %>% mean()) %>%
  #   distinct(gene, .keep_all = T) %>%
  #   dplyr::ungroup() %>%
  #   select(gene, gene_name, age_group, tpm_median, tpm_mean) 
  # 
  # 
  # # saveRDS(object = dds_tpm,
  # #         file = paste0("/home/sruiz/PROJECTS/splicing-project/markdowns/age_stratification/tpm.rds"))
  # 
  # saveRDS(object = tpm,
  #         file = paste0("/home/sruiz/PROJECTS/splicing-project/markdowns/age_stratification/tpm.rds"))
  # 
  # 
  # ###################################################################
  # ## LOAD THE TMP DATA
  # ###################################################################
  # 
  # tpm <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/markdowns/age_stratification/tpm.rds"))
  # 
  # tpm_tidy <- tpm %>%
  #   filter(gene_name %in% c('SF3B4','U2AF2','SF3B4','SF3A3','SMNDC1','GPKOW','SMNDC1','U2AF2','SF3B1','BUD13','EFTUD2','U2AF2','XRN2','XRN2','LSM11','U2AF1','EFTUD2','CDC40','BUD13','RBM22','DDX42','KHSRP','U2AF1','PRPF8','TDP43','SUPV3L1','TBRG4','RBM22')   )
  # 
  # 
  # wilcox.test(x = tpm_tidy %>% filter(age_group == "20-39") %>% pull(tpm_median),
  #             y = tpm_tidy %>% filter(age_group == "40-59") %>% pull(tpm_median),
  #             alternative = "greater")
  # wilcox.test(x = tpm_tidy %>% filter(age_group == "20-39") %>% pull(tpm_median),
  #             y = tpm_tidy %>% filter(age_group == "60-79") %>% pull(tpm_median),
  #             alternative = "greater")
  # 
  # ggplot(data = tpm_tidy) +
  #   geom_density(mapping = aes(x = tpm_mean, fill = age_group), alpha = 0.6) + 
  #   facet_zoom(xlim = c(0,100))
  # 
  # tpm_tidy %>%
  #   select(-gene, -tpm_mean) %>%
  #   spread(key = age_group, value = tpm_median)
  # 
  # tpm_tidy_hm <- tpm_tidy %>%
  #   group_by(age_group, gene_name) %>%
  #   distinct(tpm_median, .keep_all = T) %>%
  #   select(gene_name, age_group, tpm_median) 
  # 
  # ggplot(data = tpm_tidy_hm, aes(x = age_group, y = gene_name, fill= tpm_median)) + 
  #   geom_tile()
  # 
  # tpm_tidy_hm %>% as.data.frame()
  # 
  # 
  # # filter(gene %in% c("ENSG00000160201", "ENSG00000063244", "ENSG00000136450",
  # #                    "ENSG00000115524", "ENSG00000011304" )) %>%
  # tpm %>%
  #   spread(key = age_group, value = tpm_median)
  
}

tidy_sample_metadata <- function(rse) {
  
  sample_metadata <- rse %>% 
    SummarizedExperiment::colData() %>%
    as_tibble() %>%
    filter(gtex.smafrze != "EXCLUDE")
  
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
  
  
  # 8. Save the covariates ---------------------------
  
  sample_metadata_t <- as.data.frame(t(sample_metadata %>%
                                         dplyr::select(all_of(covariates))))
  
  # ## No N/A
  # sample_metadata_tidy[rowSums(is.na(sample_metadata_tidy)) > 0, ] %>% nrow()
  # 
  # ## No columns with zero or constant values
  # sample_metadata_tidy %>% head %>% as.data.frame()
  
  # pc <- prcomp(x = sample_metadata_tidy, center = TRUE, scale. = T)
  # 
  # attributes(pc)
  # pc$center
  # print(pc)
  # summary(pc)
  # 
  # var_explained <- pc$sdev^2/sum(pc$sdev^2)
  # var_explained[1:5]
  
  return(sample_metadata_t)
}

################################
## CALLS       
################################

## RBP EXPRESSION ANALYSIS


# 0. Prepare the necessary data ------------------------------------------------

## Load reference gtf and get only the genes
ref <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
ref <- ref %>% GenomeInfoDb::keepSeqlevels(c(1:22,"X","Y"), pruning.mode = "coarse")
ref <- ref %>%
  as_tibble() %>%
  filter(type == "gene") %>%
  select(seqnames, start, end, strand, gene_id, gene_name)%>%
  GRanges()


RBPs <- read.table(file = 'markdowns/experiment_report_2022_3_28_16h_50m.tsv', sep = '\t', header = TRUE) %>%
  as_tibble() %>%
  distinct(Target.gene.symbol)
RBPs_tidy <- merge(x = RBPs %>% data.table::as.data.table() %>% dplyr::rename(gene_name = Target.gene.symbol), 
                   y = ref %>% data.table::as.data.table() %>% select(gene_id, gene_name),
                   by = "gene_name")


project_id <- "BLOOD"

rse <- recount3::create_rse_manual(
  project = project_id,
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene")



# 1. Get the RBP TPM expression using GTEX data --------------------------------

dds_tpm <- get_GTEx_gene_expression(rse, ref)

dds_tpm %>% head


# 2. Get the best covariates to correct for ------------------------------------

## 1. pivot longer the 'dds_tpm' object

dds_tpm_pv <- dds_tpm %>%
  dplyr::filter(gene %in% RBPs_tidy$gene_id) %>%
  select(gene, recount_id, tpm_log10) %>% 
  group_by(gene) %>%
  distinct(recount_id, .keep_all = T) %>%
  pivot_wider(names_from = recount_id, values_from = tpm_log10, values_fill = 0)

is.na(dds_tpm_pv) <- sapply(dds_tpm_pv, is.infinite)
dds_tpm_pv[is.na(dds_tpm_pv)]<-0

write.csv(x = dds_tpm_pv %>%
            tibble::column_to_rownames(var = "gene"),
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                        "/results/pipeline3/rbp/tpm.csv"))


## 2. input the 'dds_tpm_pv' to the PCA function

pca <- prcomp(x = dds_tpm_pv %>%
                tibble::column_to_rownames(var = "gene"), 
              center = TRUE, scale. = T)


# 3. Select PCA1 - ~20 models
pca_20models <- pca$rotation[,1:20]


# 4. Tidy the metadata
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
            tibble::rownames_to_column(var = "covariate"), 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                        project_id, "/results/pipeline3/rbp/covariates.csv"))


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
