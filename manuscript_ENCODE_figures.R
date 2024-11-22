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

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/manuscript_ENCODE_figures.R")


## CONNECT TO THE DATABASE ------------------------------
supportive_reads <- 1
gtf_version <- 105
analysis_type = "shRNA"
project_name <- paste0("ENCODE_SR_", supportive_reads, "read_", analysis_type)


database_folder <- here::here(file.path("database/", project_name, gtf_version))


database_path <- paste0(database_folder,  "/", project_name, ".sqlite")
con <- dbConnect(RSQLite::SQLite(), database_path)
tables <- dbListTables(con)

## SET PATHS TO FOLDERS

base_folder <- here::here()

args <-
  list(
    dependencies_folder = file.path(here::here(), "dependencies"),
    results_folder = file.path(here::here(), "results", project_name, gtf_version, "_paper_review", "results"),
    figures_folder = file.path(here::here(), "results", project_name, gtf_version, "_paper_review", "figures"),
    data_folder = file.path(here::here(), "results", project_name, gtf_version, "_paper_review", "data")
  )


dir.create(file.path(args$results_folder), recursive = TRUE, showWarnings = F)
dir.create(file.path(args$figures_folder), recursive = TRUE, showWarnings = F)
dir.create(file.path(args$data_folder), recursive = TRUE, showWarnings = F)


## QUERY MASTER TABLES 

query = paste0("SELECT * FROM 'metadata'")
master_metadata <- dbGetQuery(con, query) %>% as_tibble()
all_projects <- master_metadata$SRA_project %>% unique
all_projects %>% length

query <- paste0("SELECT * FROM 'intron'")
master_introns <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'novel'")
master_novel_junctions <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'transcript'")
master_transcript <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'gene'")
master_gene <- dbGetQuery(con, query) %>% as_tibble()


## UTILS FUNCTION

get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}


custom_ggtheme <-  theme(text = element_text(size = 7, colour = "black"),
                         axis.ticks = element_line(colour = "black", linewidth = 2),
                         axis.text = element_text(size = 7, colour = "black"),
                         axis.line = element_line(colour = "black"),
                         axis.title = element_text(size = 7,  colour = "black"),
                         axis.text.y = element_text(size = 7, colour = "black"),
                         axis.text.x = element_text(size = 7, colour = "black", hjust = 0.5, vjust = 0.5),
                         strip.text = element_text(size = 7, colour = "black"),
                         legend.text = element_text(size = "7", colour = "black"),
                         legend.title = element_blank(),
                         legend.position = "top")


##############################################
## FUNCTIONS - Figures for the Splicing paper
##############################################


## SECTION 1 ---------------------------------------------

get_database_stats <- function() {
  
  tables <- dbListTables(con)
  
  ## Methods: number of samples by RIN number
  master_metadata %>% nrow() 
  
  master_metadata %>%
    filter(rin >= 8) %>% nrow()
  master_metadata %>%
    filter(rin >= 7) %>% nrow()
  master_metadata %>%
    filter(rin >= 6) %>% nrow()
  
  if ( master_metadata %>% filter(rin < 6) %>% nrow() > 1 ) {
    print("ERROR! Some of the samples considered present a RIN number lower than 6!")
    break;
  }
  
  
  master_metadata %>%
    dplyr::count(cluster) %>%
    print(n = 50)
  
  ## We found that 268,988 (82.8%) annotated introns had at least a single associated novel donor or acceptor junction,
  ## with only 55,968 annotated introns appearing to be precisely spliced across all the samples and tissues studied.
  
  master_introns %>% head()
  master_introns %>% distinct(ref_junID) %>% nrow()
  master_introns %>%
    dplyr::count(misspliced)
  
  
  ## Collectively, we detected 3,865,268 unique novel junctions, equating to 14 novel junctions per annotated intron. 
  master_novel_junctions %>% head
  master_novel_junctions %>% nrow() 
  master_novel_junctions %>% 
    dplyr::count(novel_type)
  
  (master_novel_junctions %>% nrow()) / (master_introns %>% filter(misspliced==1) %>% nrow())
  master_novel_junctions$seqnames %>% unique
  
  ## Collectively, we detected 31,811 genes and 199,551 transcripts
  query <- paste0("SELECT * FROM 'transcript'")
  master_transcripts <- dbGetQuery(con, query) %>% as_tibble() 
  master_transcripts %>% nrow()
  
  ## JOIN WITH GENE DATA
  query <- paste0("SELECT * FROM 'gene'")
  master_genes <- dbGetQuery(con, query) %>% as_tibble() 
  master_genes %>% distinct(gene_id) %>% nrow()
  
  ## Novel junctions exceed in X fold to annotated introns
  (master_novel_junctions %>% distinct(novel_junID) %>% nrow()) / (master_introns %>% distinct(ref_junID) %>% nrow())
  
  
  ## Percentage of mis-spliced introns
  (((master_introns %>%
       dplyr::count(misspliced) %>%
       filter(misspliced == 1) %>%
       pull(n)) * 100 )) /
    ((master_introns %>% distinct(ref_junID) %>% nrow()) %>%
       round(digits = 1))
  
  
  ## Collectively, we detected X novel junctions
  master_novel_junctions %>% distinct(novel_junID) %>% nrow()
  
  
  ## equating to 14 novel junctions per an annotated junction.
  ((master_novel_junctions %>% distinct(novel_junID) %>% nrow()) / (master_introns %>%
                                                                      dplyr::count(misspliced) %>%
                                                                      filter(misspliced == "Yes") %>%
                                                                      pull(n))) %>%
    round()
  
  
  ## After accounting for sample number, we found that the highest numbers of unique 
  ## novel junctions were identified in X tissue with the lowest numbers in Y tissue, 
  
  db_metadata_tidy <- master_metadata %>%
    dplyr::group_by(cluster) %>%
    mutate(n=n()) %>%
    dplyr::select(SRA_project,cluster,n) %>%
    distinct( SRA_project, cluster, .keep_all = T)
  
  df_novel_jxn_count <- map_df( (db_metadata_tidy$SRA_project %>% unique()), function (project_id)  {
    
    # project_id <- (db_metadata_tidy$SRA_project %>% unique())[1]
    all_clusters <- db_metadata_tidy %>%
      filter(SRA_project == project_id) %>%
      pull(cluster) %>%
      unique()
    
    map_df(all_clusters, function (cluster_id)  {
      
      # cluster_id <- all_clusters[1]
      print(paste0(cluster_id))
      query <- paste0("SELECT COUNT(DISTINCT novel_junID), AVG(MSR_D), AVG(MSR_A) 
                      FROM '", cluster_id, "_", project_id, "_misspliced'")
      novel_junctions <- dbGetQuery(con, query) %>% as.double()
      
      return(data.frame(cluster = cluster_id,
                        n_novel = novel_junctions[1],
                        mean_msrd = novel_junctions[2],
                        mean_msra = novel_junctions[3]))
    })
    
  })
  
  
  ## The detection of unique novel donor and acceptor junctions was a common finding across all tissues, 
  ## with the highest numbers found in “Cells - EBV-transformed lymphocytes” tissue and the lowest in “Whole Blood”.
  
  df_n_samples_msr <- db_metadata_tidy %>%
    dplyr::rename(n_sample = n) %>%
    inner_join(y = df_novel_jxn_count, by = "cluster") %>%
    mutate(prop_novel = n_novel/n_sample)
  
  df_n_samples_msr %>%
    arrange(desc(prop_novel)) %>%
    as.data.frame()
  
  
  df_n_samples_msr$cluster <- df_n_samples_msr$cluster %>% as.factor() 
  plot_data <- df_n_samples_msr %>%
    mutate(avg_n_novel_sample = n_novel/n_sample) %>%
    ungroup() %>%
    arrange(desc(avg_n_novel_sample))%>%
    mutate(cluster = fct_inorder(cluster))
  
  ggplot(data =  plot_data) +
    geom_bar(mapping = aes(x = cluster, 
                           y = avg_n_novel_sample),
             stat = "identity") +
    ylab("Average number of unique novel\njunctions across samples") +
    xlab("GTEx tissue") +
    theme_light() +
    custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  
  file_name <- paste0(args$figures_folder, "/avg_novel_jnx_per_tissue")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)
  
}


get_common_introns_across_experiments <- function (RBPs = NULL,
                                                   required_clusters = NULL) {
  
  ## Get list of RBPs
  if ( is.null(RBPs) || RBPs == "all" ) {
    
    all_projects <- master_metadata$SRA_project %>% unique()
    RBPs <- "all"
  } else {
    all_projects <- RBPs
  }
  
  if ( is.null(required_clusters) ) {
    required_clusters <- c("case", "control")
  }
  
  
  print(paste0(Sys.time(), " - getting unique and common junctions across ENCODE experiments..."))
  
  ## Getting all introns that are common across all ENCODE experiments -------------------------------------------
  
  all_introns <- list()
  
  for (project_id in all_projects) {
    
    # project_id <- all_projects[1]
    print(paste0(Sys.time(), " - ", project_id))
    
    ## GET THE CLUSTERS
    
    all_clusters <- master_metadata %>%
      filter(SRA_project == project_id,
             cluster %in% required_clusters) %>%
      distinct(cluster) %>%
      pull()
    
    print(all_clusters)
    
    for(cluster_id in all_clusters) { 
      
      query <- paste0("SELECT DISTINCT ref_junID 
                      FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT DISTINCT ref_junID 
                      FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      
      all_introns[[paste(c(project_id,cluster_id), collapse = "_")]] <- introns$ref_junID %>% unique()
      
      print(paste0(Sys.time(), " - ", introns$ref_junID %>% unique() %>% length(), 
                   " unique introns collected from '", project_id, "'"))
      
    }
    
  }
  
  common_introns <- data.frame(ref_junID = Reduce(intersect,  all_introns))
  common_introns %>% head()
  common_introns %>% nrow()
  
  
  ## Getting splicing data from common introns across all ENCODE experiments -------------------------------------------
  
  common_introns_all_experiments <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <- master_metadata %>%
      filter(SRA_project == project_id,
             cluster %in% required_clusters) %>%
      distinct(cluster) %>%
      pull()
    
    
    map_df(all_clusters, function(cluster_id) { 
      
      print(paste0(" --> ", cluster_id))
      
      query <- paste0("SELECT * 
                      FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT * 
                      FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
      
      introns %>%
        filter(ref_junID %in% common_introns$ref_junID) %>%
        mutate(RBP = project_id, 
               cluster = cluster_id) %>% 
        as_tibble() %>%
        return()
      
    })
  })
  
  ## The number of unique introns should be the same across experiments
  common_introns_all_experiments %>%
    dplyr::group_by(RBP, cluster) %>%
    distinct(ref_junID) %>%
    dplyr::count() %>%
    ungroup() %>%
    print
  
  ## The number of novel junctions linked to the unique introns studied can differ across experiments
  common_introns_all_experiments %>%
    dplyr::group_by(RBP, cluster)%>%
    distinct(novel_junID) %>%
    dplyr::count() %>%
    ungroup()%>%
    print
  
  file_name <- paste0(args$results_folder, "/common_introns_details_",
                      paste(RBPs,collapse = "_"), "_experiments_", 
                      paste(required_clusters,collapse = "_"), ".rds")
  saveRDS(object = common_introns_all_experiments, file = file_name)
  
  return(common_introns_all_experiments)
}




##################################
## MAIN FIGURE 5
##################################

get_data_main_figure5_c <- function(replace = T) {
  
  ## LOAD DATABASE METADATA
  
  metadata_RBPs <- if (analysis_type=="shRNA") {
    master_metadata %>% 
      pivot_longer(c("Splicing regulation", 
                     "Spliceosome", 
                     "Exon Junction Complex", 
                     "NMD", 
                     "Novel_RBP", 
                     "RNA modification"), names_to = "Category") %>%
      filter(value == 1) %>%
      dplyr::select(-value) %>%
      distinct(target_gene, sample_id, .keep_all = T)
  } else { master_metadata }
  
  metadata_RBPs$Category %>% unique
  
  
  ## Load common introns ---------------------------------------
  
  target_genes = "all"
  required_clusters <- c("case", "control")
  
  common_introns_path <- paste0(args$results_folder, "/common_introns_details_", 
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  common_introns_all_experiments <- if (!file.exists(common_introns_path)) {
    message("Calculating common introns...")
    get_common_introns_across_experiments(RBPs = target_genes, required_clusters = c("case", "control"))
  } else {
    message("Loading common introns across '", target_genes, "' experiments...")
    readRDS(file = common_introns_path)
  }
  
  
  ## Run the wilcox test on MSR_D --------------------------
  
  if (!file.exists(paste0(args$results_folder, "/ENCODE_effectsize_MSRD.rds")) || replace) {
    
    message("Calculating MSR_D effect sizes...")
    
    # Load the MSR tables
    MSR_D <- common_introns_all_experiments %>%
      group_by(RBP, cluster) %>%
      distinct(ref_junID, .keep_all = T) %>%
      ungroup() %>%
      dplyr::select(ref_junID, MSR_D, RBP, cluster) %>%
      pivot_wider(id_cols = ref_junID, 
                  names_from = c("RBP", "cluster"), 
                  values_from = c("MSR_D") )
    MSR_D %>% head()
    
    MSR_D_tests <- GenerateMSRtests(target_RBPs = all_projects, 
                                    MSR_Table = MSR_D, 
                                    file_output = paste0(args$results_folder, "/ENCODE_effectsize_MSRD.rds"), 
                                    overwrite = T,
                                    num_cores = 8)
    
    ## Add the categories
    if (analysis_type=="shRNA") {
      MSR_D_tests <- MSR_D_tests %>%
        left_join(y = metadata_RBPs %>% distinct(target_gene, Category), by = "target_gene") 
    }
    
    ## Add bonferroni correction
    MSR_D_tests <- MSR_D_tests %>% mutate(FDR = NA, .after = p.value)
    MSR_D_tests$FDR <- p.adjust(MSR_D_tests$p.value, method = "fdr")
    
    ## Save results
    saveRDS(object = MSR_D_tests, file = paste0(args$results_folder, "/ENCODE_effectsize_MSRD.rds"))
    write_csv(x = MSR_D_tests %>%
                distinct(target_gene, .keep_all=T) %>%
                mutate(statistical_test = "Wilcoxon Rank text: rstatix::wilcox_test(data, formula, paired = TRUE, correct = TRUE, alternative = 'greater')",
                       H0 = "The MSR_D observations in case and control samples are symmetric about their median value.",
                       H1 = "The MSR_D observations in case samples are greater at their median value than the MSR_D observations in control samples."),
              file = paste0(args$results_folder, "/ENCODE_effectsize_MSRD.csv"), col_names = T)
    
  } else {
    message("Loading 'ENCODE_effectsize_MSRD.rds' file ...")
    MSR_D_tests <- readRDS(file = paste0(args$results_folder, "/ENCODE_effectsize_MSRD.rds"))
  }
  
  if (!file.exists(paste0(args$results_folder, "/ENCODE_effectsize_MSRA.rds"))) {
    
    message("Calculating MSR_A effect sizes...")
    
    MSR_A <- common_introns_all_experiments %>%
      group_by(RBP, cluster) %>%
      distinct(ref_junID, .keep_all = T) %>%
      ungroup() %>%
      dplyr::select(ref_junID, MSR_A, RBP, cluster) %>%
      pivot_wider(id_cols = ref_junID, 
                  names_from = c("RBP", "cluster"), 
                  values_from = c("MSR_A"))
    MSR_A %>% head()
    
    MSR_A_tests <- GenerateMSRtests(target_RBPs = all_projects, 
                                    MSR_Table = MSR_A, 
                                    file_output = paste0(args$results_folder, "/ENCODE_effectsize_MSRA.rds"),  
                                    overwrite = T, 
                                    num_cores = 8)
    
    ## Add the categories
    if (analysis_type=="shRNA") {
      MSR_A_tests <- MSR_A_tests %>%
        left_join(y = metadata_RBPs %>% distinct(target_gene, Category), by = "target_gene")
    }
    
    ## Add bonferroni correction
    MSR_A_tests <- MSR_A_tests %>% distinct(target_gene, .keep_all=T) %>% mutate(FDR = NA, .after = p.value)
    MSR_A_tests$FDR <- p.adjust(MSR_A_tests$p.value, method = "fdr")
    
    ## Save results
    saveRDS(object = MSR_A_tests, file = paste0(args$results_folder, "/ENCODE_effectsize_MSRA.rds"))
    write_csv(x = MSR_A_tests %>%
                distinct(target_gene, .keep_all=T) %>%
                mutate(statistical_test = "Wilcoxon Rank text: rstatix::wilcox_test(data, formula, paired = TRUE, correct = TRUE, alternative = 'greater')",
                       H0 = "The observations MSR_A in case samples vs control samples are symmetric about their median value.",
                       H1 = "The observations MSR_A in case samples are greater at their median value than the observations MSR_A in control samples."),
              file = paste0(args$results_folder, "/ENCODE_effectsize_MSRA.csv"))
  } else {
    message("Loading 'ENCODE_effectsize_MSRA.rds' file ...")
    MSR_A_tests <- readRDS(file = paste0(args$results_folder, "/ENCODE_effectsize_MSRA.rds"))
  }
  
  
  ##########################
  ## PAPER STATS
  ##########################
  
  ## 1. TEST  ---------------------------------------------------------------------------------------------------------
  
  ## MSR_D
  ## Firstly, it revealed a significant increase in mis-splicing rates in samples with gene knockdowns compared to untreated controls for 90% of the 54 genes considered
  
  (MSR_D_tests %>%
      distinct(target_gene, .keep_all = T) %>%
      filter(FDR <= 0.05) %>%
      distinct(target_gene) %>%
      nrow() * 100) / ( MSR_D_tests %>%
                          distinct(target_gene, .keep_all = T) %>%
                          distinct(target_gene)%>%
                          nrow() )
  
  MSR_D_tests %>% distinct(target_gene, .keep_all = T) %>% filter(FDR <= 0.05) %>% distinct(target_gene, .keep_all=T) 
  MSR_D_tests %>% distinct(target_gene, .keep_all = T) %>% filter(FDR <= 0.05) %>% distinct(target_gene, .keep_all=T) %>% pull(FDR) %>% summary
  
  
  
  ## MSR_A
  ## increase in mis-splicing rates in samples with gene knockdowns compared to untreated controls for 90% of the 54 genes considered
  
  (MSR_A_tests %>%
      filter(FDR <= 0.05) %>%
      distinct(target_gene) %>%
      nrow() * 100) / ( MSR_A_tests %>%
                          distinct(target_gene)%>%
                          nrow() )
  
  MSR_A_tests %>% distinct(target_gene, .keep_all = T) %>% filter(FDR <= 0.05) %>% distinct(target_gene, .keep_all=T)
  MSR_A_tests %>% distinct(target_gene, .keep_all = T) %>% filter(FDR <= 0.05) %>% distinct(target_gene, .keep_all=T) %>% pull(FDR) %>% summary
  
 
  ## 2ND TEST ------------------------------------------------------------------------------------
  
  ## Knockdowns of the splicing machinery components tended to have a greater effect on 3’ss than 5’ss mis-splicing 
  ## MSR_D
  
  MSR_D_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05,
           effect_size==max(effect_size)) %>%
    distinct(target_gene, .keep_all=T) 
  
  MSR_D_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05) %>%
    distinct(target_gene, .keep_all=T) %>%
    pull(effect_size) %>%
    summary
  
  ## MSR_A
  MSR_A_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05,
           Category == "Spliceosome") %>%
    #distinct(target_gene, .keep_all=T) %>%
    pull(effect_size) %>%
    summary
  
  
  
  ## (...) except for 6 genes (including SAFB2 which is not thought to impact on splicing and so was used as a negative control) 
  
  c(MSR_D_tests %>%
      filter(FDR > 0.05)  %>%
      distinct(target_gene, .keep_all=T) %>% pull(target_gene),
    MSR_A_tests %>%
      filter(FDR > 0.05)  %>%
      distinct(target_gene, .keep_all=T) %>% pull(target_gene)) %>% unique
  
  
  ## Notably, AQR, EIF4A3, SF3A3, U2AF1, U2AF2, and MAGOH knockdowns resulted in the highest increases in 5’ss 
  
  MSR_D_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05) %>%
    arrange(desc(effect_size) ) 
  
  
  ## Notably, AQR, EFTUD2, HNRNPC, MAGOH, SF3A3, SF3B4 and U2AF1 knockdowns resulted in the highest increases in 3’ss mis-splicing
  
  MSR_A_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05) %>%
    arrange(desc(effect_size) )
}
  
main_figure5_c <- function() {
  
  ## Run the wilcox test on MSR_A --------------------------
  
  message("Loading 'ENCODE_effectsize_MSRD.rds' file ...")
  MSR_D_tests <- readRDS(file = paste0(args$results_folder, "/ENCODE_effectsize_MSRD.rds")) %>% distinct(target_gene, Category, .keep_all =T)
  
  message("Loading 'ENCODE_effectsize_MSRA.rds' file ...")
  MSR_A_tests <- readRDS(file = paste0(args$results_folder, "/ENCODE_effectsize_MSRA.rds")) %>% distinct(target_gene, Category, .keep_all =T)
  
  RBPs_subgroups <- readxl::read_excel("dependencies/RBPs_subgroups.xlsx")
  
  
  #Tidy
  MSR_D_tests <- MSR_D_tests %>% 
    left_join(y = RBPs_subgroups %>% gather(category, type, -id,-name) %>% filter(type == 1), 
              by = c("target_gene" = "name")) %>%
    mutate(category = ifelse(is.na(category), Category, category)) %>% 
    dplyr::select(-c(id, Category, type)) %>% dplyr::rename(Category=category)
  MSR_A_tests <- MSR_A_tests  %>% 
    left_join(y = RBPs_subgroups %>% gather(category, type, -id,-name) %>% filter(type == 1), 
              by = c("target_gene" = "name")) %>%
    mutate(category = ifelse(is.na(category), Category, category)) %>% 
    dplyr::select(-c(id, Category, type)) %>% dplyr::rename(Category=category) 
  
  
  ##########################
  ## TIDY EFFECT SIZES
  ##########################
  
  # Combine both MSR_A and MSR_D -------------------
  MSR_combined = rbind(MSR_A_tests %>% mutate(MSR_type = "MSR_A"), 
                       MSR_D_tests %>% mutate(MSR_type = "MSR_D")) %>% 
    filter(FDR <= 0.05 | target_gene == "SAFB2")## 'SAFB2' is a novel RBP that will be used as control

  #MSR_combined[is.na(MSR_combined$Category),] %>% print(n=100)
  #MSR_combined$Category[is.na(MSR_combined$Category)] <- "Novel_RBP"
  
  max_genes <- 10
  
  # Filter the Splicing regulation category
  filter_splicing_regulation <- MSR_combined %>%
    distinct(target_gene, .keep_all = T) %>% 
    arrange(-effect_size) %>%
    filter(Category == "Splicing regulation") %>%
    head(max_genes) %>%
    pull(target_gene)
  filter_spliceosome <- MSR_combined %>%
    distinct(target_gene, .keep_all = T) %>% 
    arrange(-effect_size) %>%
    filter(Category == "Spliceosome") %>%
    head(max_genes) %>%
    pull(target_gene)
  
  MSR_graph_data <- MSR_combined %>%
    filter(target_gene %in% c(filter_splicing_regulation, filter_spliceosome, "UPF1", "UPF2", "MAGOH", "SAFB2")) %>% #| Category != c("Splicing regulation") ) %>%
    arrange(-effect_size) %>%
    mutate(target_gene = factor(target_gene, levels = .$target_gene %>% unique)) %>%
    distinct(target_gene, MSR_type, .keep_all = T)
  
  if (analysis_type == "shRNA") {
    MSR_graph_data <- MSR_graph_data %>%
      mutate(# Use factors to sort the graph
        Category = factor(Category, levels = c("Splicing regulation", 
                                                "Spliceosome", 
                                                "NMD", 
                                                "Exon Junction Complex", 
                                                "Novel_RBP")))  # Use factors to sort the graph)
    print(MSR_graph_data$Category %>% unique())
  }
 
  
  ##########################
  ## METADATA KEFF
  ##########################
  
  metadata_kEff_path <- file.path("ENCODE_SR/metadata_WB_kEff.tsv")
  metadata_kEff <- readr::read_delim(metadata_kEff_path, show_col_types = F) %>%
    mutate(kEff_text = ifelse(is.na(kEff_avg), kEff_avg, paste0(round(kEff_avg), "%"))) %>%
    mutate(kEff_text = kEff_text %>% as.factor())

  MSR_graph_data <- MSR_graph_data %>%
    left_join(y = metadata_kEff, by = "target_gene") %>% 
    dplyr::select(-c(statistical_test, H0,H1)) %>% 
    distinct(target_gene, Category, MSR_type, .keep_all = T)
  
  MSR_graph_data %>% head()
  
  ## Save data
  write.csv(x = MSR_graph_data, file = file.path(args$data_folder,"figure5_c.csv"), row.names = T)
  
  ##########################
  ## PLOT
  ##########################
  
  MSR_graph_data <- MSR_graph_data %>%
    mutate(Category = as.character(Category)) %>%
    mutate(Category = ifelse(Category == "Exon Junction Complex", "EJC", Category)) %>%
    mutate(Category = ifelse(Category == "Novel_RBP", "Control", Category)) %>%
    mutate(Category = str_replace(string = Category, pattern = "_", replacement = " ")) %>%
    mutate(Category = ifelse(Category == "Splicing regulation", str_to_title(string = Category), Category)) %>%
    mutate(Category = factor(Category, 
                             levels = c("Splicing Regulation", "Spliceosome",  "NMD", "EJC", "Control" ))) %>%
    arrange(Category, desc(effect_size)) %>%  # First by 'type', then by 'prop' descending
    mutate(target_gene = factor(target_gene, levels = unique(target_gene)))  # Reorder target_gene
  MSR_graph_data$MSR_type = factor(MSR_graph_data$MSR_type, levels = c( "MSR_A","MSR_D"))
  
  # Plot the graph
  plot_effectsize <- ggplot(MSR_graph_data, aes(x = target_gene, y = effect_size)) + 
    geom_bar(aes(fill = MSR_type), 
             stat = "identity", color = "black", 
             linewidth = 0.25, width = 0.80, position = "dodge") + 
    scale_y_continuous(limits = c(0,0.8),
                       expand = expansion(mult = c(0, 0.02)), 
                       breaks = seq(0, 0.8, 0.1),
                       labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "", "")) +
    scale_x_discrete(expand = expansion(add = c(0.7, 0.7))) +
    labs(x = "Target gene shRNA knockdown", 
         y = "Probability of superior MSR in\ngene knockdown vs. untreated samples") + 
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      labels = c("MSR_A" = "MSR Acceptor", "MSR_D" = "MSR Donor"),
                      breaks = c("MSR_D", "MSR_A"))  +
    guides(fill = guide_legend(title = NULL, order = 2, ncol = 2,  nrow = 1 )) +
    theme_light()  +
    ggforce::facet_row(facets = vars(Category), 
                       scales = "free_x", space = "free",
                       #labeller = labeller(Category = category_labels),
                       drop = T,
                       shrink = T) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.box = "horizontal") +
    guides(colour = guide_legend(override.aes = list(fill = '#999999'), title = NULL, label.position = "bottom", order = 1)) +
   ggnewscale::new_scale_fill() +
   geom_tile(stat = "identity",
             aes(y = 0.68, fill = kEff_avg, color = "No data\navailable"),
             linewidth = 0.5, width = 1, height = 0.045) +
   geom_text(aes(y = 0.76, label = kEff_text), color = "black", size = 2.5) +
   viridis::scale_fill_viridis(option = "inferno",
                               # na.value = "#999999",
                               name = "Knockdown\nEfficiency",
                               limits = c(0, 100),
                               breaks = c(NA, seq(0, 100, 25)),
                               labels = c("none", paste0(seq(0, 100, 25), "%")),
                               guide = guide_colourbar(frame.colour = "black",
                                                       frame.linewidth = 0.4,
                                                       order = 1,
                                                       ticks.colour = "black",
                                                       barwidth = 10,
                                                       barheight = 1.5)) +
   scale_colour_manual(values = c( "No data\navailable" = "black")) +
   custom_ggtheme +
     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.box = "horizontal") +
     guides(colour = guide_legend(override.aes = list(fill = '#999999'), title = NULL, label.position = "bottom", order = 1))
  
  
  plot_effectsize
  # Save the graph
  figure_name <- file.path(args$figures_folder, "main_figure5c")
  ggsave(file = paste0(figure_name, ".png"),  plot = plot_effectsize, width = 180, height = 90, units = "mm", dpi = 300)
  ggsave(file = paste0(figure_name, ".svg"),  plot = plot_effectsize, width = 180, height = 90, units = "mm", dpi = 300)
  
  
  
}


##################################
## MAIN FIGURE 6
##################################

main_figure6_a <- function() {
  
  ## LOAD COMMON INTRONS AND NOVEL FOR AQR AND U2AF2
  target_RBPs <- c("AQR", "U2AF2")
  required_clusters <- c("case", "control")
  common_introns_path <- paste0(args$results_folder, "/common_introns_details_", 
                                paste(target_RBPs, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  ## LOAD COMMON INTRONS
  common_introns <- if (file.exists(common_introns_path)) {
    readRDS(file = common_introns_path) %>% as_tibble()
  } else {
    get_common_introns_across_experiments(RBPs = target_RBPs,
                                          required_clusters = required_clusters)
  }
  
  limit_bp = 30
  
  RBP_novel <- common_introns %>%
    filter(RBP %in% target_RBPs) %>%
    inner_join(master_introns %>% dplyr::select(ref_junID), by = "ref_junID") %>%
    left_join(master_novel_junctions %>% dplyr::select(novel_junID, novel_type, distance), by = "novel_junID") %>%
    group_by(RBP, cluster) %>%
    distinct(novel_junID, .keep_all = TRUE) %>%
    ungroup() %>%
    mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " ")) %>%
    filter(abs(distance) < limit_bp) %>%
    mutate(RBP = factor(x = RBP, levels = target_RBPs)) %>%
    dplyr::select(distance, cluster, novel_type, RBP) %>%
    mutate(cluster = factor(ifelse(cluster == "case" & novel_type == "novel donor", "gene knockdown donor",
                                   ifelse(cluster == "case" & novel_type == "novel acceptor", "gene knockdown acceptor", "control")),
                            levels = c("control", "gene knockdown donor", "gene knockdown acceptor"))) %>%
    mutate(novel_type = str_to_title(novel_type)) 
  
  ## Save source data
  
  write_csv(x = RBP_novel, file = file.path(args$data_folder, "figure6_a.csv"), col_names = T)
  
  #################################################
  ## PLOT AQR AND U2AF1 DISTANCES
  #################################################
  
  distances_plot <- ggplot(data = RBP_novel) + 
    geom_histogram(aes(x = distance, fill = fct_rev(cluster)), 
                   bins = 60,  binwidth = 1,  position = "identity", alpha = 1,  color = "black",  linewidth = 0.1) +
    scale_x_continuous(breaks = seq(-limit_bp, limit_bp, length.out = 5)) + 
    scale_fill_manual(values = c("#35B779FF", "#8d03b0", "#666666"),
                      breaks = c("gene knockdown donor", "gene knockdown acceptor", "control"),
                      labels = c("gene knockdown donor", "gene knockdown acceptor", "control")) +
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    facet_grid(fct_rev(novel_type)~RBP) +
    guides(fill = guide_legend(title = "Sample type: ", ncol = 2, nrow = 1 )) +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 60), fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 33),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 30, ymax = 31), fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15, y = 48),  size = 3, label = "intron") +
    theme_void()
  
  distances_plot <- distances_plot / (distance_rectangle +  distance_rectangle) + patchwork::plot_layout(heights = c(8, 1))
  distances_plot
  
  figure_name <- file.path(args$figures_folder, "main_figure6_a")
  ggsave(file = paste0(figure_name, ".png"), width = 180, height = 80, dpi = 300, units = "mm")
  ggsave(file = paste0(figure_name, ".svg"), width = 180, height = 80, dpi = 300, units = "mm")
  
  
}

main_figure6_b <- function() {
  
  
  ## LOAD COMMON INTRONS AND NOVEL FOR AQR AND U2AF2
  
  target_RBPs <- c("AQR", "U2AF2")
  required_clusters <- c("case", "control")
  common_introns_path <- paste0(args$results_folder, "/common_introns_details_", 
                                paste(target_RBPs, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  ## LOAD COMMON INTRONS
  common_introns <- if (file.exists(common_introns_path)) {
    readRDS(file = common_introns_path) %>% as_tibble()
  } else {
    get_common_introns_across_experiments(RBPs = target_RBPs,
                                          required_clusters = required_clusters)
  }
  
  ## Add MES 5' and 3' info to the common introns and join with novel junctions
  RBP_novel <- common_introns %>%
    filter(RBP %in% target_RBPs) %>%
    mutate(cluster = factor(ifelse(cluster == "case", "gene knockdown", "control"), 
                            levels = c("control", "gene knockdown"))) %>%
    inner_join(master_introns %>% dplyr::select(ref_junID, ref_mes5ss, ref_mes3ss), by = "ref_junID") %>%
    left_join(master_novel_junctions %>% 
                dplyr::select(novel_junID, novel_type, distance, novel_mes5ss, novel_mes3ss), 
              by = "novel_junID") %>%
    group_by(RBP, cluster) %>%
    distinct(novel_junID, .keep_all = TRUE) %>%
    ungroup() %>%
    mutate(delta_ss5score = ref_mes5ss - novel_mes5ss, 
           delta_ss3score = ref_mes3ss - novel_mes3ss)
  
  
  ## Statistical Tests
  run_wilcox_test <- function(data, RBP) {
    wilcox.test(x = data %>% filter(RBP == RBP, cluster == "gene knockdown") %>% pull(delta_ss3score),
                y = data %>% filter(RBP == RBP, cluster == "control") %>% pull(delta_ss3score),
                paired = FALSE,
                alternative = "greater")
  }
  RBP_stats <- map(target_RBPs, ~run_wilcox_test(RBP_novel, .))
  
  
  #################################################
  ## PREPARE DATA TO PLOT 
  #################################################

  
  RBP_novel_acceptor_delta <- RBP_novel %>%
    filter(novel_type == "novel_acceptor") %>%
    dplyr::group_by(RBP, cluster) %>%
    mutate(medianMSR = median(delta_ss3score)) %>%
    ungroup() %>%
    mutate(cluster = factor(cluster, levels = c("gene knockdown", "control"))) %>%
    dplyr::group_by(RBP) %>%
    mutate(
      p.value = format(
        wilcox.test(delta_ss3score ~ cluster, 
                    data = cur_data(), 
                    alternative = "greater",
                    subset = cluster %in% c("gene knockdown", "control"),
                    exact = T)$p.value,
        digits = 2, scientific = TRUE
      )
    ) %>%
    ungroup()%>%
    mutate(p.value = ifelse(as.numeric(p.value) == 0, "2.2e-16", p.value)) %>%
    mutate(RBP = factor(RBP, levels = target_RBPs)) %>%
    dplyr::select(delta_ss3score, cluster, medianMSR, RBP, p.value) 
  
  
  # RBP_novel_acceptor_delta <- RBP_novel %>%
  #   filter(novel_type == "novel_acceptor") %>%
  #   
  #   dplyr::group_by(RBP, cluster) %>%
  #   mutate(medianMSR = median(delta_ss3score)) %>%
  #   ungroup() %>%
  #   dplyr::group_by(RBP) %>%
  #   mutate(p.value = format(wilcox.test(delta_ss3score ~ cluster)$p.value, digits = 2, scientific = TRUE)) %>%
  #   ungroup %>%
  #   mutate(p.value = ifelse(as.numeric(p.value) == 0, "2.2e-16", p.value)) %>%
  #   mutate(RBP = factor(RBP, levels = target_RBPs)) %>%
  #   dplyr::select(delta_ss3score, cluster, medianMSR, RBP, p.value) 
    
  
  write.csv(RBP_novel_acceptor_delta, file = file.path(args$data_folder, "figure6_b.csv"), row.names = FALSE)
  
  #################################################
  ## PLOT DELTA MES ONLY ACCEPTOR - AQR and U2AF2
  #################################################
  
  delta_mes_acceptor <- ggplot(data = RBP_novel_acceptor_delta)  +
    geom_density(aes(x = delta_ss3score, fill = cluster), alpha = 0.8, linewidth = 0.3, color = "black") +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = medianMSR, color = cluster), linetype = "dashed", linewidth = 0.9) +
    facet_wrap(vars(RBP)) +
    geom_text(data = distinct(RBP_novel_acceptor_delta, RBP, p.value), 
              aes(label = paste0("P<", p.value)), x = 25, y = 0.08, size = 3, color = "#333333") +
    scale_fill_manual(values = c("#64037d", "#999999"), breaks = c("gene knockdown", "control"), labels = c("shRNA knockdown", "Control")) + 
    scale_colour_manual(values = c("#64037d", "#333333"), breaks = c("gene knockdown", "control"), labels = c("shRNA knockdown", "Control")) + 
    labs(x = "Delta MES Acceptor") +
    theme_light() +
    custom_ggtheme
  
  delta_mes_acceptor
  
  figure_name <- file.path(args$figures_folder, "main_figure6_b")
  ggsave(file = paste0(figure_name, ".png"), width = 180, height = 60, dpi = 300, units = "mm")
  ggsave(file = paste0(figure_name, ".svg"), width = 180, height = 60, dpi = 300, units = "mm")
  
  ## Supplementary Figure for AQR distances
  # supplementary_figure14(RBP_novel)

}

main_figure6_c <- function() {
  

  ################################
  ## LOAD METADATA & RBPs
  ################################
  
  all_projects <- unique(master_metadata$SRA_project)
  target_genes <- "all"
  required_clusters <- c("case", "control")
  
  file_name <- paste0(args$results_folder, "/common_introns_details_", target_genes, "_experiments_", 
                      paste(required_clusters, collapse = "_"), ".rds")
  
  common_introns_all_experiments <- if (file.exists(file_name)) {
    message("Loading common introns across all experiments...")
    readRDS(file = file_name)
  } else {
    get_common_introns_across_experiments(RBPs = target_genes, required_clusters = required_clusters)
  }
  
  # Summary of introns across experiments
  common_introns_summary <- common_introns_all_experiments %>% group_by(RBP, cluster) %>% distinct(ref_junID) %>% dplyr::count() %>% ungroup()
  
  ####################################
  ## LOAD MSR VALUES
  ####################################
  
  MSR_RBPs <- common_introns_all_experiments %>%
    inner_join(master_introns %>% dplyr::select(ref_junID, seqnames, start, end, strand), by = "ref_junID") %>%
    dplyr::select(ref_junID, MSR_D, MSR_A, RBP, cluster, seqnames, start, end, strand) %>%
    group_by(RBP, cluster) %>%
    distinct(ref_junID, .keep_all = TRUE) %>%
    ungroup()
  
  ####################################
  ## MERGE ENCODE and iCLIP DATA
  ####################################
  
  iclip_results_path <- file.path(args$results_folder, "/iCLIP_ENCODE_all_introns_chisq.csv")
  
  if (!file.exists(iclip_results_path)) {
    
    process_iCLIP_data <- function(target_RBP) {
      message("Processing RBP: ", target_RBP)
      
      # Load iCLIP data for HepG2 and K562
      iclip_RBP <- map_df(c("HepG2", "K562"), function(celltype) {
        file_path <- file.path(args$dependencies_folder, "ENCORI", 
                               paste0("ENCORI_hg38_RBPTarget_", target_RBP, "_", celltype, ".txt"))
        
        if (file.exists(file_path)) {
          readr::read_delim(file_path, delim = "\t", col_names = TRUE, skip = 3, show_col_types = FALSE)
        } else {
          NULL
        }
      })
      
      if (nrow(iclip_RBP) > 1) {
        
        process_MSR <- function(MSR_type) {
          MSR_RBP <- MSR_RBPs %>%
            filter(RBP == target_RBP) %>%
            pivot_wider(id_cols = ref_junID, names_from = cluster, values_from = all_of(MSR_type)) %>%
            inner_join(MSR_RBPs %>% filter(RBP == target_RBP) %>% dplyr::select(ref_junID, seqnames, start, end, strand), by = "ref_junID") %>%
            mutate(start = start - 100, end = end + 100)
          
          MSR_RBP_inc <- MSR_RBP %>%
            filter(case > control) %>%
            distinct(ref_junID, .keep_all = TRUE) %>%
            GRanges()
          
          MSR_RBP_notinc <- MSR_RBP %>%
            filter(case <= control) %>%
            distinct(ref_junID, .keep_all = TRUE) %>%
            GRanges()
          
          iclip_RBP_gr <- iclip_RBP %>%
            dplyr::select(seqnames = "chromosome", start = "broadStart", end = "broadEnd", strand) %>%
            GRanges()
          
          # Find overlaps with iCLIP
          MSR_RBP_inc_overlaps <- findOverlaps(query = MSR_RBP_inc, subject = iclip_RBP_gr, type = "any")
          MSR_RBP_notinc_overlaps <- findOverlaps(query = MSR_RBP_notinc, subject = iclip_RBP_gr, type = "any")
          
          list(MSR_RBP_inc_overlaps = MSR_RBP_inc_overlaps, MSR_RBP_notinc_overlaps = MSR_RBP_notinc_overlaps)
        }
        
        # Run the process_MSR function for both MSR_D and MSR_A
        result_MSR_D <- process_MSR("MSR_D")
        result_MSR_A <- process_MSR("MSR_A")
        
        return(list(result_MSR_D = result_MSR_D, result_MSR_A = result_MSR_A))
      } else {
        message("No iCLIP HepG2|K562 data for ", target_RBP)
        return(NULL)
      }
    }
    
    iCLIP_binding_MSR_results <- map_df(all_projects, process_iCLIP_data)
    write.csv(iCLIP_binding_MSR_results, file = iclip_results_path, row.names = FALSE)
    
  } else {
    message("Loading existing iCLIP results...")
    iCLIP_binding_MSR_results <- read.csv(iclip_results_path)
  }
  
  
  ########################################
  ## PAPER STATS & CALCULATE FOLD-CHANGES
  ########################################
  
  
  # Filter out SAFB2 and summarize p-values
  summary_stats <- function(msr_type) {
    iCLIP_binding_MSR_results %>%
      filter(MSR_tested == msr_type, RBP != "SAFB2") %>%
      pull(chisq_pvalue) %>%
      summary()
  }
  
  summary_D <- summary_stats("MSR_D")
  summary_A <- summary_stats("MSR_A")
  
  # Calculate fold-changes
  fold_changes <- iCLIP_binding_MSR_results %>%
    rowwise() %>%
    mutate(MSR_increasing_w_motif = MSR_increasing_w_RBP_motif / (MSR_increasing_w_RBP_motif + MSR_increasing_N_RBP_motif),
           MSR_increasing_n_motif = MSR_increasing_N_RBP_motif / (MSR_increasing_w_RBP_motif + MSR_increasing_N_RBP_motif),
           
           MSR_not_increasing_w_motif = MSR_not_increasing_w_RBP_motif / (MSR_not_increasing_w_RBP_motif + MSR_not_increasing_N_RBP_motif),
           MSR_not_increasing_n_motif = MSR_not_increasing_N_RBP_motif / (MSR_not_increasing_w_RBP_motif + MSR_not_increasing_N_RBP_motif)) %>%
    
    mutate(fold_change_w_motif = log2(MSR_increasing_w_motif) - log2(MSR_not_increasing_w_motif)) %>% 
    mutate(fold_change_not_motif = log2(MSR_increasing_n_motif) - log2(MSR_not_increasing_n_motif)) %>%
    
    dplyr::select(RBP, 
                  MSR_tested, 
                  "Introns containing binding sites for each RBP" = fold_change_w_motif,
                  "Introns not containing binding sites for each RBP" = fold_change_not_motif) %>%
    
    gather(type,  "2fold", -MSR_tested, -RBP) %>% 
    arrange( MSR_tested,  desc("2fold"))  %>%
    ungroup
  
  
  ########################################################
  ## ADD KNOCKDOWN EFFICIENCY
  ########################################################
  
  metadata_kEff_path <- file.path("ENCODE_SR/metadata_WB_kEff.tsv")
  metadata_kEff <- readr::read_delim(metadata_kEff_path, show_col_types = F) %>%
    mutate(kEff_text = ifelse(is.na(kEff_avg), kEff_avg, paste0(round(kEff_avg), "%"))) %>%
    mutate(kEff_text = kEff_text %>% as.factor())
  
  iCLIP_binding_MSR_results_prop <- fold_changes %>%
    left_join(y = metadata_kEff, by = c("RBP" = "target_gene")) %>%
    filter(RBP != "SAFB2")
  
  # Export data
  write.csv(x = iCLIP_binding_MSR_results_prop, file = file.path(args$data_folder, "/figure6_c.csv"), row.names = F)
  
  
  ####################################
  ## PLOT RESULTS
  ####################################
  
  iCLIP_binding_MSR_results_prop <- iCLIP_binding_MSR_results_prop %>%
    mutate(MSR_tested = factor(MSR_tested, levels = c("MSR_D", "MSR_A"), labels = c("MSR Donor", "MSR Acceptor")))
  
  
  
  iCLIP_plot <-  ggplot(data = iCLIP_binding_MSR_results_prop,
                        aes(x = RBP, y = `2fold`, group = factor(MSR_tested)) ) +
    geom_bar(stat = "identity",
             aes(x = tidytext::reorder_within(x = RBP, by = -`2fold`, within = type), y = `2fold`,
                 fill = factor(MSR_tested)), position = position_dodge(width = 0.9)) +
    ylab("Log2 Fold Change of introns with increasing MSR\nvalues in shRNA knockdown vs control sampmles") +
    xlab("") +
    theme_light() +
    ggforce::facet_row(~ type, scales = "free_x") +
    tidytext::scale_x_reordered() +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2,  nrow = 1 ))  +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("MSR Donor","MSR Acceptor")) +
    theme(plot.margin = margin(t = 5, r = 5, b = -5, l = 5), legend.box.margin = margin(l = -10, b = -10, t = -5 )) +
    ggnewscale::new_scale_fill() +
    geom_tile(aes(y = 1.4, x = tidytext::reorder_within(x = RBP, by = -`2fold`, within = type), fill = kEff_avg), 
              linewidth = 0.5, width = 1, height = 0.1, stat = "identity") + 
    viridis::scale_fill_viridis(option = "inferno", 
                                name = "Knockdown\nEfficiency", 
                                limits = c(0, 100), 
                                breaks = c(NA, seq(0, 100, 25)),
                                labels = c("none", paste0(seq(0, 100, 25), "%")),
                                guide = guide_colourbar(frame.colour = "black", 
                                                        frame.linewidth = 0.4,
                                                        order = 1, 
                                                        ticks.colour = "black", 
                                                        barwidth = 10, 
                                                        barheight = 1.5)) +
    scale_colour_manual(values = c( "No data\navailable" = "black")) +
    guides(colour = guide_legend(override.aes = list(fill = '#999999'), title = NULL, label.position = "bottom", order = 1))
  
  iCLIP_plot
  
  
  figure_name <- file.path(args$figures_folder, "/main_figure6c")
  ggsave(file = paste0(figure_name, ".png"), width = 180, height = 75, units = "mm", dpi = 300)
  ggsave(file = paste0(figure_name, ".svg"), width = 180, height = 75, units = "mm", dpi = 300)
  
  
}



##################################
## SUPPLEMENTARY FIGURES
##################################

supplementary_figure7 <- function() {
  
  ## RBP knockdowns were aligned
  ## junctions from shRNA control experiments followed by RNA-sequencing data downloaded from the ENCODE platform
  ## were extracted using a minimum anchor length of 8 bp required to call the presence of a junction 
  ## (regtools junctions extract -a 8, https://regtools.readthedocs.io/en/latest/).
  ## This represented an increase in stringency of 3 bp with respect to the anchor used in the GTEx v8 data (5bp by default).
  
  required_clusters <- "control"
  target_genes = "all"
  
  common_introns_path <- paste0(args$results_folder, "/common_introns_details_",
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters,collapse = "_"), ".rds")
  
  if( file.exists(common_introns_path) ) {
    
    message("Loading common introns...")
    common_introns_all_experiments <- readRDS(file = common_introns_path)
    
  } else {
    ## Load common introns
    message("Calculating common introns...")
    common_introns_all_experiments <- get_common_introns_across_experiments(RBPs = "all", required_clusters = required_clusters)
  }
  
  common_introns_all_experiments %>% distinct(RBP) %>% dplyr::count()
  common_introns_all_experiments %>% group_by(RBP) %>% distinct(ref_junID) %>% dplyr::count() %>% ungroup()
  common_introns_all_experiments_w_novel <- common_introns_all_experiments %>% 
    left_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_type, distance), by = "novel_junID")
  
  #################################################
  ## SUPPLEMENTARY FIGURE 7.a
  #################################################
  
  df_proportion_jxn <- map_df(all_projects, function(project_id) {
    # project_id <- all_projects[1]
    print(paste0(Sys.time(), " - ", project_id))
    common_introns_local_experiment <- common_introns_all_experiments_w_novel %>% filter(RBP == project_id)
    
    annotated_junc <- common_introns_local_experiment %>% distinct(ref_junID) %>% nrow()
    donor_junc <- common_introns_local_experiment %>% filter(novel_type == "novel_donor") %>% distinct(novel_junID) %>% nrow()
    acceptor_junc <- common_introns_local_experiment %>% filter(novel_type == "novel_acceptor") %>% distinct(novel_junID) %>% nrow()
    
    annotated_prop <- annotated_junc/(annotated_junc + donor_junc + acceptor_junc)
    donor_prop <- donor_junc/(annotated_junc + donor_junc + acceptor_junc)
    acceptor_prop <- acceptor_junc/(annotated_junc + donor_junc + acceptor_junc)
    
    return(data.frame(RBP = project_id, annotated_prop = annotated_prop, donor_prop = donor_prop, acceptor_prop = acceptor_prop))
  })
  
  
  df_proportion_jxn_tidy <- df_proportion_jxn %>%
    tidyr::gather(key = "type", value = "prop", -RBP) %>%
    mutate(type = str_remove_all(string = type, pattern = "_prop")) %>%
    filter(type != "annotated") %>%
    mutate(prop = prop * 100)
  
  df_proportion_jxn_tidy$type = factor(df_proportion_jxn_tidy$type, 
                                       levels = c("donor","acceptor"))
  
  unique_jxn_violin_plot <- ggplot(df_proportion_jxn_tidy, aes(type, prop, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point(size = 0.5, colour = "#333333") +
    geom_line( aes(group = RBP), colour = "#333333",linewidth = .2 )  +
    theme_light() +
    ylab("% unique junctions") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), text = element_text(colour = "black", size = 12), legend.position = "top") +
    scale_x_discrete(breaks = c("acceptor", "donor"), labels = c("Novel Acceptor", "Novel Donor")) +
    scale_fill_manual(values = c("#35B779FF", "#8d03b0"), breaks = c("donor", "acceptor"), labels = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) +
    custom_ggtheme + 
    theme(legend.position = "right")
  
  
  #################################################
  ## SUPPLEMENTARY FIGURE 7.b
  #################################################
  
  
  df_proportion_reads <- map_df(all_projects, function(project_id) {
    # project_id <- all_projects[1]
    print(paste0(Sys.time(), " - ", project_id))
    common_introns_local_experiment <- common_introns_all_experiments_w_novel %>% filter(RBP == project_id)
    
    annotated_reads <- common_introns_local_experiment %>% distinct(ref_junID, .keep_all = T) %>% pull(ref_sum_counts) %>% sum()
    donor_reads <- common_introns_local_experiment %>% filter(novel_type == "novel_donor") %>% distinct(novel_junID, .keep_all = T) %>% pull(novel_sum_counts) %>% sum()
    acceptor_reads <- common_introns_local_experiment %>% filter(novel_type == "novel_acceptor")%>% distinct(novel_junID, .keep_all = T) %>% pull(novel_sum_counts) %>% sum()
    
    annotated_prop <- annotated_reads/(annotated_reads + donor_reads + acceptor_reads)
    donor_prop <- donor_reads/(annotated_reads + donor_reads + acceptor_reads)
    acceptor_prop <- acceptor_reads/(annotated_reads + donor_reads + acceptor_reads)
    
    return(data.frame(RBP = project_id, annotated = annotated_prop, donor = donor_prop, acceptor = acceptor_prop))
  })
  
  df_proportion_reads_tidy <- df_proportion_reads %>% tidyr::gather(key = "type", value = "prop", -RBP) %>% filter(type != "annotated") %>% mutate(prop = prop * 100)
  df_proportion_reads_tidy$type = factor(df_proportion_reads_tidy$type, levels = c("donor", "acceptor"))
  
  reads_violin_plot <- ggplot(df_proportion_reads_tidy, aes(type, prop, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point(size = 0.5, colour = "#333333") +
    geom_line( aes(group = RBP), colour = "#333333",linewidth = .2  )  +
    theme_light() +
    ylab("% cumulative read counts") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), text = element_text(colour = "black", size = 12), legend.position = "top") +
    scale_x_discrete(breaks = c( "acceptor", "donor"), labels = c( "Novel Acceptor", "Novel Donor")) +
    scale_fill_manual(values = c( "#35B779FF", "#8d03b0"), breaks = c( "donor", "acceptor"),labels = c( "Novel Donor", "Novel Acceptor")) + 
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))+
    custom_ggtheme + 
    theme(legend.position = "right")
  
  ggpubr::ggarrange(unique_jxn_violin_plot, reads_violin_plot, common.legend = T, align = "h", ncol = 2, nrow = 1) +
    theme(plot.margin = margin(t = -5, r = -5, b = -10, l = -5), legend.box.margin=margin(b = -10))
  
  ggsave(file = paste0(args$figures_folder , "/supplementary_figure7ab.png"), width = 130, height = 50, dpi = 300, units = "mm")
  
  
  #################################################
  ## SUPPLEMENTARY FIGURE 7.c
  #################################################
  
  RBP_delta_MES <- common_introns_all_experiments_w_novel %>%
    left_join(y = master_introns %>% dplyr::select(ref_junID, ref_mes5ss, ref_mes3ss), by = "ref_junID")%>%
    left_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_mes5ss, novel_mes3ss), by = "novel_junID") %>% 
    distinct(novel_junID, .keep_all = T) %>%
    mutate(delta_ss5score = ref_mes5ss - novel_mes5ss, delta_ss3score = ref_mes3ss - novel_mes3ss) %>%
    dplyr::select(novel_junID, ref_junID, novel_type, delta_ss5score, delta_ss3score, RBP) %>%
    drop_na() %>%
    dplyr::select(delta_ss5score, delta_ss3score, novel_type) %>%
    gather(key = delta_type, value = deltaMES) %>%
    mutate(delta_type = ifelse(delta_type == "delta_ss5score", "MES Donor", "MES Acceptor")) %>%
    mutate(delta_type = delta_type %>% as.factor()) %>%
    mutate(deltaMES = deltaMES %>% as.double()) %>%
    filter(deltaMES != 0) 

  
  delta_mes <- ggplot(data = RBP_delta_MES )  +
    geom_density(aes(x = deltaMES, fill = delta_type), alpha = 0.9) +
    geom_vline(xintercept = 0) +
    ggforce::facet_col(~fct_rev(delta_type), strip.position = "right") +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("MES Donor", "MES Acceptor"), labels = c("Novel Donor", "Novel Acceptor")) +
    labs(x = "Delta MES", y = "Density") +
    theme_light() +
    custom_ggtheme + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), legend.box.margin = margin(b = -11) ) + 
    theme(legend.position = "none")
  
  delta_mes
  ggsave(file = paste0(args$figures_folder, "/supplementary_figure7c.png"), width = 45, height = 50, dpi = 300, units = "mm")
  
  
  
  #################################################
  ## SUPPLEMENTARY FIGURE 7.d
  #################################################
  
  limit_bp = 30
  
  common_introns_all_experiments_w_protein <- common_introns_all_experiments_w_novel %>% 
    left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding), by = "ref_junID")
  
  plot_distances_PC <- ggplot(common_introns_all_experiments_w_protein %>%
                                filter(protein_coding %in% c(0,100)) %>%
                                mutate(protein_coding_label = ifelse (protein_coding == 0, "Noncoding", "Protein-coding")) %>%
                                mutate(protein_coding_label = protein_coding_label %>% as.factor()) %>%
                                mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " ")) %>%
                                distinct(novel_junID, .keep_all = T) %>% 
                                filter(abs(distance) < limit_bp) ) + 
    geom_histogram(aes(x = distance, fill = fct_rev(novel_type)), bins = 60, binwidth = 1, position = "identity", alpha = 1, color = "black", linewidth = 0.1) +
    scale_x_continuous(breaks = seq(-limit_bp, limit_bp, length.out = 5)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("novel donor", "novel acceptor"), labels = c("Novel donor", "Novel acceptor")) +
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    facet_grid(fct_rev(protein_coding_label)~fct_rev(novel_type), scales = "free_y") +
    guides(fill = guide_legend(title = "Sample type: ", ncol = 1, nrow = 2 )) +
    theme_light() +
    custom_ggtheme + 
    theme(legend.position = "none", plot.margin = margin(t = 0, r = -5, b = 0, l = -5), legend.box.margin=margin(b = -10))
  plot_distances_PC
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 60), fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 33),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 30, ymax = 31), fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15, y = 48),  size = 3, label = "intron") +
    theme_void()
  
  plot_distances_PC / (distance_rectangle +  distance_rectangle) + patchwork::plot_layout(heights = c(8, 1))
  
  ggsave(file = paste0(args$figures_folder , "/supplementary_figure7d.png"), width = 125, height = 60, dpi = 300, units = "mm")
  
  #################################################
  ## SUPPLEMENTARY FIGURE 7.e
  #################################################
  
  df_modulo_experiments <- common_introns_all_experiments_w_novel %>%
    group_by(RBP) %>%
    distinct(novel_junID, .keep_all = T) %>%
    filter(abs(distance) <= 100) %>% 
    mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " ")) %>%
    mutate(modulo = abs(distance) %% 3) %>%
    ungroup
  
  df_modulo_experiments <- df_modulo_experiments %>% 
    group_by(RBP,modulo) %>%
    summarise(n = n()) %>%
    mutate(freq = (n / sum(n))*100) 
  
  df_modulo_experiments$modulo = factor(df_modulo_experiments$modulo, levels = c( "0", "1", "2"))
  
  ggplot(df_modulo_experiments, aes(x = freq, y = modulo)) +
    ggridges::geom_density_ridges_gradient() +
    ylab("Modulo3 of the distance") +
    xlab("% of novel junctions") +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(30,35,40), labels = c("30%","35%","40%")) +
    scale_y_discrete(expand = c(0,0.5,1,0))
  
  
  ggsave(file = paste0(args$figures_folder , "/supplementary_figure7e.png"), width = 50, height = 60, dpi = 300, units = "mm")
  
  
  
  #################################################
  ## SUPPLEMENTARY FIGURE 7.F
  #################################################
  
  df_all_introns_tidy <- common_introns_all_experiments_w_protein %>%
    dplyr::select(ref_junID, MSR_D, MSR_A, RBP, protein_coding) %>%
    gather(key = "MSR_type", value = "MSR", -ref_junID, -RBP, - protein_coding) %>%
    mutate(type_label = ifelse(MSR == 0, "Accurate Splicing", "Mis-splicing")) 
  
  df_all_introns_tidy <- df_all_introns_tidy %>%
    mutate(percentile_group = case_when(MSR == 0 ~ "0",
                                        MSR > 0 & MSR <= 0.2 ~ "(0,.2]",
                                        MSR > 0.2 & MSR <= 0.4 ~ "(.2,.4]",
                                        MSR > 0.4 & MSR <= 0.6 ~ "(.4,.6]",
                                        MSR > 0.6 & MSR <= 0.8 ~ "(.6,.8]",
                                        MSR > 0.8 & MSR <= 1 ~ "(.8,1]"))
  
  
  df_all_introns_tidy <- df_all_introns_tidy %>% group_by(RBP, MSR_type) %>% distinct(ref_junID, .keep_all = T) %>% ungroup()
  df_all_introns_tidy$MSR_type <- factor(df_all_introns_tidy$MSR_type, levels = c("MSR_D", "MSR_A"))
  
  plot1 <- ggplot(data = df_all_introns_tidy) + 
    geom_bar(aes(x = type_label, fill = MSR_type), position = "dodge") +
    ggplot2::labs(x = "", y = "Annotated introns") + theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("MSR_D","MSR_A"), labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) + 
    theme(plot.margin = margin(0,10,0,0), legend.box.margin=margin(b = -5))
  plot1
  
  ## 2. ZOOMED PLOT VERSION
  df_all_introns_tidy$percentile_group = factor(df_all_introns_tidy$percentile_group, levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  plot2 <- ggplot(data = df_all_introns_tidy %>% filter(percentile_group != "0")) +
    geom_bar(aes(x = percentile_group, fill = MSR_type),position = "dodge")+
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    ggforce::facet_zoom(ylim = c(0,22000), split = TRUE) +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("MSR_D","MSR_A"), labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  plot2
  
  ggpubr::ggarrange(plot1, plot2, labels = c("", ""), nrow = 2,  ncol = 1, heights = c(1,1.5), common.legend = T)
  ggsave(file = paste0(args$figures_folder,"/supplementary_figure7f.png"), width = 135, height = 70, dpi = 300, units = "mm")
  
  
  ########################################################
  ## SUPPLEMENTARY FIGURE 7.g
  ########################################################
  
  ## TIDY DATAFRAME
  df_biotype_result_tidy <- common_introns_all_experiments_w_protein %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(biotype = ifelse(protein_coding == 100, "PC", "non PC"))  %>%
    distinct(ref_junID, .keep_all = T) %>%
    group_by(ref_junID) %>%
    mutate(mean_coverage = (sum(ref_sum_counts)/sum(ref_n_individuals)) %>% log10()) %>%
    ungroup()
  
  df_biotype_result_tidy %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)
  
  ## PLOT
  plot_BS <- ggplot(data = df_biotype_result_tidy %>% distinct(ref_junID,.keep_all = T)) +
    geom_density(mapping = aes(x = mean_coverage, fill = biotype), alpha = 0.9) +
    ggtitle("Before subsampling") +
    xlab("log10 mean expression level") +
    theme_light() +
    custom_ggtheme +
    ggsci::scale_fill_npg(name = "Transcript biotype: ")
  
  plot_BS
  #### START SUBSAMPLIING
  
  ## lncRNA
  df_lncRNA_tidy <- df_biotype_result_tidy %>% 
    distinct(ref_junID, .keep_all = T) %>%
    dplyr::filter(biotype == "non PC")
  ## protein-coding
  df_protein_coding_tidy <- df_biotype_result_tidy %>% 
    distinct(ref_junID, .keep_all = T) %>%
    dplyr::filter(biotype  == "PC")
  
  data_combined <- rbind(df_protein_coding_tidy, df_lncRNA_tidy) %>% dplyr::select(-distance)
  
  ## QC
  if (!identical(data_combined$ref_junID %>% sort(),
                 df_biotype_result_tidy %>% distinct(ref_junID,.keep_all = T) %>% pull(ref_junID) %>% sort())) {
    stop("ERROR!")
  }
  
  
  if ( !file.exists(paste0(args$results_folder, "/MSR_subsample.rds")) ) {
    ## Subsampling introns to control by similarity in mean read coverage
    m.out <- MatchIt::matchit(biotype ~ mean_coverage, 
                              data = data_combined, 
                              distance = data_combined$mean_coverage,
                              method = "nearest", 
                              caliper = c(mean_coverage = 0.005), 
                              std.caliper = FALSE)
    
    subsample <- MatchIt::match.data(m.out)
    subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)
    
    saveRDS(object = subsample, file = paste0(args$results_folder, "/anchor_MSR_subsample.rds"))
  } else {
    subsample <- readRDS(file = paste0(args$results_folder, "/anchor_MSR_subsample.rds"))
  }
  
  
  plot_AS <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = biotype), alpha = 0.9) +
    ggtitle("After subsampling") +
    theme_light() +
    custom_ggtheme +
    xlab("log10 mean expresion level")+
    ggsci::scale_fill_npg(name = "Transcript biotype: ")
  
  ggpubr::ggarrange(plot_BS, plot_AS, labels = c("a", "b"), common.legend = T)
  
  #file_name <- paste0(args$figures_folder, "/anchor_MSR_subsampling")
  #ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  ## Analysing differences in MSR_D and MSR_A
  ## subset introns from lncRNA transcripts after subsampling
  subsample_lncRNA <- subsample %>% dplyr::filter(biotype=="non PC")
  ## subset introns from protein-coding transcripts after subsampling
  subsample_protein_coding <- subsample %>% dplyr::filter(biotype=="PC")
  ## QC
  if (intersect(subsample_protein_coding$ref_junID, subsample_lncRNA$ref_junID) %>% length() > 0 )  {
    stop("ERROR! some introns have been categorised as lncRNA and protein-coding.")
  }
  
  df_introns_biotype <- subsample %>% filter(protein_coding %in% c(0,100)) %>% mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC")) 
  
  any(df_introns_biotype %>% pull(ref_junID) %>% duplicated())
  
  df_introns_biotype$biotype = factor(df_introns_biotype$biotype, levels = c("non PC","PC"))
  df_introns_biotype <- df_introns_biotype %>% dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>% gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)
  df_introns_biotype$MSR_type = factor(df_introns_biotype$MSR_type, levels = c("MSR_A", "MSR_D"))

  print(paste0(Sys.time(), " - ", df_introns_biotype %>% nrow(), " - introns after tidying!"))
  
  df_introns_biotype %>% filter(biotype == "PC") %>% pull(MSR) %>% summary()
  df_introns_biotype %>% filter(biotype == "non PC") %>% pull(MSR) %>% summary()
  
  df_introns_biotype_tidy <- df_introns_biotype %>%
    group_by(biotype, MSR_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    mutate(percentile_group = case_when(MSR == 0 ~ "0",
                                        MSR > 0 & MSR <= 0.2 ~ "(0,.2]",
                                        MSR > 0.2 & MSR <= 0.4 ~ "(.2,.4]",
                                        MSR > 0.4 & MSR <= 0.6 ~ "(.4,.6]",
                                        MSR > 0.6 & MSR <= 0.8 ~ "(.6,.8]",
                                        MSR > 0.8 & MSR <= 1 ~ "(.8,1]"))
  
  df_introns_biotype_tidy$MSR_type = factor(df_introns_biotype_tidy$MSR_type, levels = c("MSR_D","MSR_A"))
  df_introns_biotype_tidy$percentile_group = factor(df_introns_biotype_tidy$percentile_group, levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  df_introns_biotype_tidy$biotype = factor(df_introns_biotype_tidy$biotype, levels = c("PC","non PC"))
  
  MSR_label <- c(MSR_D = "MSR Donor", MSR_A = "MSR Acceptor")
  
  ggplot(data = df_introns_biotype_tidy) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .3, color = "#333333")+
    #ggtitle("MSR Donor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    ggforce::facet_col(~MSR_type, strip.position = "top", labeller = labeller(MSR_type=MSR_label))+
    #ylim(c(0,23000))+
    scale_fill_manual(values = c("#333333","#cccccc"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding   ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = "",
                               ncol = 2, nrow = 1)) + 
    theme( plot.margin = margin(t = 0,b = 0,r = 0,l = 0),
           legend.box.margin=margin(b = -10) )+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
  
  file_name <- paste0(args$figures_folder, "/supplementary_figure7g")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 70, height = 60, units = "mm", dpi = 300)
  
  
  plotMSR_donor_zoomed <- ggplot(data = df_introns_biotype_tidy) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    ggforce::facet_col(~MSR_type, strip.position = "top", labeller = labeller(MSR_type=MSR_label))+
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    scale_y_continuous(limits =c(0,800), position = "right") +
    
    scale_fill_manual(values = c("#333333","#cccccc"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
    theme_light() +
    custom_ggtheme +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(title = NULL,  
                               ncol = 2, nrow = 1)) 
  plotMSR_donor_zoomed
  file_name <- paste0(args$figures_folder, "/supplementary_figure7g.png")
  ggplot2::ggsave(filename = file_name, width = 70, height = 60, units = "mm", dpi = 300)
  
}

supplementary_figure14 <- function () {
  
  
  #################################################
  ## LOAD COMMON INTRONS AND NOVEL JUNCTIONS
  ## FROM UPF1 and UPF2
  ## ADD BIOTYPE INFORMATION
  #################################################
  
  ## As we are only comparing UPF1 and UPF2, we get the common introns between the two genes
  
  target_genes = c("UPF1","UPF2")
  required_clusters <- c("case", "control")
  
  ## Get annotated common introns -------
  
  common_introns_path <- paste0(args$results_folder, "/common_introns_details_", 
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  NMD_intron <- if( file.exists(common_introns_path) ) {
     readRDS(file = common_introns_path)
  } else {
    get_common_introns_across_experiments(RBPs = target_genes,
                                          required_clusters = c("case", "control"))
  }
  
  ## Add information from the master INTRON table -----------------------------------
  ## Given that we are studying the NMD action, we only subset introns 
  ## from protein-coding transcripts
  
  
  NMD_intron_w_protein_coding <- NMD_intron %>%
    left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding,mean_phastCons17way5ss_100,mean_phastCons17way3ss_100,ref_mes5ss, ref_mes3ss),
              by = "ref_junID") %>%
    filter(protein_coding == 100)
  
  
  
  ## ONLY 37615 annotated introns from protein-coding transcripts overlap the two NMD datasets
  NMD_intron_w_protein_coding %>% dplyr::group_by(RBP, cluster) %>%  distinct(ref_junID) %>% dplyr::count() %>% ungroup() 
  
  
  #################################################
  ## GET COMMON INTRONS THAT ARE MIS-SPLICED &
  ## ACCURATELY SPLICED IN CONTROL AND CASE EXPERIMENTS
  #################################################
  
  ## 1. CONTROL EXPERIMENTS

  
  ## Get annotated introns with accurate splicing from CONTROL NMD experiments
  control_accurate_splicing <- aux_get_annotated_introns_accurate_splicing(data = NMD_intron_w_protein_coding,
                                                                           cluster_type = "control") 
  control_accurate_splicing %>% unique %>% length()
  
  ## Get introns with INACCURATE splicing across CONTROL NMD experiments
  control_mis_splicing <- aux_get_annotated_introns_inaccurate_splicing( data = NMD_intron_w_protein_coding,
                                                                         accurate_splicing_data = control_accurate_splicing, 
                                                                         cluster_type = "control", p_coding = 100) 
  control_mis_splicing %>% unique %>% length()
  
  
  
  
  ## 2. CASE EXPERIMENTS
  
  ## Get annotated introns with accurate splicing across CASE NMD experiments
  case_accurate_splicing <- aux_get_annotated_introns_accurate_splicing(data = NMD_intron_w_protein_coding,
                                                                        cluster_type = "case") 
  
  case_accurate_splicing%>% unique %>% length()
  
  ## Get annotated introns with INACCURATE splicing across CASE NMD experiments
  case_mis_splicing <- aux_get_annotated_introns_inaccurate_splicing( data = NMD_intron_w_protein_coding,
                                                                      accurate_splicing_data = case_accurate_splicing, 
                                                                      cluster_type = "case", p_coding = 100) 
  case_mis_splicing %>% unique %>% length()
  
  
  
  #################################################
  ## DEFINE NMD-RESISTANT AND NMD-SENSITIVE
  ## INTRONS
  #################################################
  
  ## 1. NMD RESISTANT INTRONS ---------------------------------------------------------------------------
  NMD_resistant_introns <- setdiff(control_accurate_splicing, case_mis_splicing) 
  
  
  ## Add master intron table info to calculate expression levels
  NMD_resistant <- NMD_intron_w_protein_coding %>%
    filter(ref_junID %in% NMD_resistant_introns) %>%
    mutate(NMD_type = "NMD Resistant",
           mean_expression = (ref_sum_counts/ref_n_individuals) %>% log10)  %>% 
    filter(cluster == "control")
  
  # NMD_resistant_introns %>% length()
  # NMD_resistant %>% dplyr::count(RBP, cluster)
  # NMD_resistant %>% dplyr::group_by(RBP, cluster) %>%  distinct(ref_junID) %>% dplyr::count() %>% ungroup() 
  
  
  
  ## 2. NMD SENSITIVE INTRONS ---------------------------------------------------------------------------
  NMD_sensitive_introns <- intersect(control_accurate_splicing, case_mis_splicing)
  
  
  ## Add master intron table info to calculate expression levels
  NMD_sensitive <- NMD_intron_w_protein_coding %>%
    filter(ref_junID %in% NMD_sensitive_introns)  %>%
    mutate(NMD_type = "NMD Sensitive",
           mean_expression = (ref_sum_counts/ref_n_individuals) %>% log10) %>% 
    filter(cluster == "control")
  
  # NMD_sensitive_introns %>% length()
  # NMD_sensitive %>% dplyr::count(RBP,cluster)
  # NMD_sensitive %>% dplyr::group_by(RBP, cluster) %>%  distinct(ref_junID) %>% dplyr::count() %>% ungroup() 
  
  
  
  
  #################################################
  ## SUBSAMPLE TO GET ANNOTATED INTRONS WITH 
  ## SIMILAR EXPRESSION LEVELS
  #################################################
  
  data_combined <- rbind(NMD_resistant, NMD_sensitive) %>%
    distinct(ref_junID, .keep_all = T)
  
  ## Subsampling introns to control by similarity in mean read coverage
  m.out <- MatchIt::matchit(NMD_type ~ mean_expression, 
                            data = data_combined, 
                            distance = data_combined$mean_expression,
                            method = "nearest", 
                            caliper = c(mean_expression = 0.005), 
                            std.caliper = FALSE)
  subsample <- MatchIt::match.data(m.out) %>% mutate(NMD_type = NMD_type %>% as.factor())
  
  # data_combined %>% dplyr::count(NMD_type)
  # subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(NMD_type)
  # subsample %>% dplyr::count(NMD_type)
  
  
  #################################################
  ## STATS FOR THE PAPER
  #################################################
  
  ## We found that mis-splicing events that were only evident in the context of NMD knockdown, were generated from annotated introns which had significantly lower phastCons17 and MES values at their 3’ss
  
  
  ## Statistical Tests
  run_wilcox_test <- function(data, col_test) {
    wilcox.test(x = data %>% filter(NMD_type == "NMD Sensitive") %>% pull(paste0(col_test)),
                y = data %>% filter(NMD_type == "NMD Resistant") %>% pull(paste0(col_test)),
                paired = T,
                alternative = "less")
  }
  
  
  ## Test 5'ss - phastCons17
  run_wilcox_test(data = subsample, "mean_phastCons17way5ss_100")
  
  
  ## Test 5'ss - MES
  run_wilcox_test(data = subsample, "ref_mes5ss")
  
  
  ## Test 3'ss - phastCons17
  run_wilcox_test(data = subsample, "mean_phastCons17way3ss_100")
  rstatix::wilcox_effsize(subsample, formula = mean_phastCons17way3ss_100 ~ NMD_type, paired = T)
  
  
  ## Test 3'ss - MES
  run_wilcox_test(data = subsample, "ref_mes3ss")
  rstatix::wilcox_effsize(subsample, formula = ref_mes3ss ~ NMD_type, paired = T)
  
  
  
  
  
  #################################################
  ## PLOTS
  #################################################
  
  plot_5ss <- ggplot(data = subsample, aes(x=mean_phastCons17way5ss_100, y=ref_mes5ss, colour=fct_rev(NMD_type))) +
    geom_point() +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "bottom") +
    xlab("phastCons17 5'ss")  +
    ylab("MaxEntScan 5'ss") +
    ylim(c(-40, 20)) 
  
  plot_3ss <- ggplot(data = subsample, aes(x=mean_phastCons17way3ss_100, y=ref_mes3ss, colour=fct_rev(NMD_type))) +
    geom_point() +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "bottom") +
    xlab("phastCons17 3'ss")  +
    ylab("MaxEntScan 3'ss") +
    ylim(c(-40, 20))
  
  
  ggpubr::ggarrange(ggExtra::ggMarginal(plot_5ss, groupColour = TRUE, groupFill = TRUE), 
                    ggExtra::ggMarginal(plot_3ss, groupColour = TRUE, groupFill = TRUE),
                    labels = c("b", "c"))
  
  
  
  ggsave(file = file.path(args$figures_folder, "/supplementary_figure14.png"), 
         width = 180, height = 90, dpi = 300, units = "mm")
  
  
  
}

supplementary_figure15_16 <- function () {
  
  ## LOAD COMMON INTRONS AND NOVEL JUNCTIONS
  
  target_genes <- c("UPF1", "UPF2")
  required_clusters <- c("case", "control")
  
  # Load common introns
  common_introns_path <- paste0(args$results_folder, "/common_introns_details_", 
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  common_introns_NMD <- if (file.exists(common_introns_path)) {
    message("Loading common introns...")
    readRDS(file = common_introns_path)
  } else {
    get_common_introns_across_experiments(RBPs = target_genes, required_clusters = required_clusters)
  }
  
  # Add protein-coding info
  common_introns_NMD <- common_introns_NMD %>%
    inner_join(master_introns %>% dplyr::select(ref_junID, protein_coding), by = "ref_junID")
  
  # Get novel junctions
  NMD_novel <- common_introns_NMD %>%
    filter(!is.na(novel_junID)) %>%
    left_join(master_novel_junctions %>% dplyr::select(novel_junID, novel_type, distance), by = "novel_junID") %>%
    mutate(type_colour = paste0(cluster, "_", novel_type),
           novel_type = ifelse(novel_type == "novel_donor", "Novel Donor", "Novel Acceptor"),
           cluster = ifelse(cluster == "case", "Knockdown (KD)", "Wild Type"))
  
  #################################################
  ## FUNCTION TO PLOT DISTANCES
  #################################################
  
  plot_nmd_distance <- function(data, gene_name, p_coding, limit_bp = 30) {
    
    # Histogram plot
    plot_data <- ggplot(data = data %>% #dplyr::select(RBP, distance, protein_coding) %>%
                          filter(RBP == gene_name,
                                 abs(distance) <= limit_bp,
                                 protein_coding == p_coding) ) + 
      geom_histogram(aes(x = distance, fill = type_colour), bins = limit_bp * 2, binwidth = 1, position = "stack") +
      facet_grid(fct_rev(novel_type) ~ fct_rev(cluster)) +
      ggtitle(paste0(gene_name, ifelse(p_coding == 100, " - Protein-coding", " - Noncoding"))) +
      xlab("Distance (in bp)") +
      ylab("Unique novel junctions") +
      theme_light() +
      scale_x_continuous(limits = c(-limit_bp, limit_bp), breaks = seq(-limit_bp, limit_bp, by = 3)) +
      scale_fill_manual(values = c("#9ce2c0","#35B779FF","#e99cfc", "#8c04ae"),
                        breaks = c("control_novel_donor", "case_novel_donor","control_novel_acceptor", "case_novel_acceptor"),
                        label = c("Donor - Wild Type", "Donor - KD", "Acceptor - Wild Type", "Acceptor - KD")) +
      guides(fill = guide_legend(title = NULL, ncol = 4, nrow = 1)) +
      custom_ggtheme
    
    if (p_coding == 0) {plot_data <- plot_data + ggpubr::rremove("legend")}
    
    # Add distance rectangle
    distance_rectangle <- ggplot() +
      geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100), fill = "grey", color = "black") +
      geom_text(aes(x = 15, y = 55), size = 3, label = "exon") +
      geom_rect(aes(xmin = -limit_bp, xmax = 0, ymin = 49, ymax = 51), fill = "grey", color = "black") +
      geom_text(aes(x = -15, y = 70), size = 3, label = "intron") +
      theme_void()
    
    # Combine plots
    combined_plot <- plot_data / (distance_rectangle + distance_rectangle) + 
      patchwork::plot_layout(heights = c(8, 1))
    
    return(combined_plot)
  }
  
  #################################################
  ## PLOT FOR EACH NMD GENE
  #################################################
  
  for (gene_name in target_genes) {
    
    # gene_name <- target_genes[1]
    
    # Plot for protein-coding
    PC_NMD <- plot_nmd_distance(data = NMD_novel, gene_name, p_coding = 100)
    
    # Plot for noncoding
    NPC_NMD <- plot_nmd_distance(NMD_novel, gene_name, p_coding = 0)
    
    # Combine and save
    combined_plots <- ggpubr::ggarrange(PC_NMD, NPC_NMD, ncol = 1, nrow = 2)
    
    file_name <- if (gene_name == "UPF1") {
      "supplementary_figure15.png"
    } else {
      "supplementary_figure16.png"
    }
    
    
    ggsave(file = file.path(args$figures_folder, file_name), 
           plot = combined_plots, width = 180, height = 180, dpi = 300, units = "mm")
  }
  
  
 
  
}

supplementary_figure17 <- function(RBP_novel) {
  
  RBP_novel <- RBP_novel %>%
    mutate(cluster = factor(cluster, levels = c("gene knockdown","control")))
  
  limit_bp = 200
  target_RBP <- "AQR"
  plot_aqr_long_distances <- ggplot(RBP_novel %>%
                                      filter(RBP == target_RBP) %>%
                                      mutate(novel_type = str_replace(string = novel_type,
                                                                      pattern = "_",
                                                                      replacement = " ")) %>% 
                                      mutate(novel_type = factor(novel_type, levels = c("novel donor", "novel acceptor"))) %>%
                                      filter(abs(distance) < limit_bp)) + 
    geom_histogram(aes(x = distance, fill = cluster),
                   bins = 60, 
                   binwidth = 1, 
                   position = "identity", 
                   alpha = 1, 
                   color = "#333333", 
                   linewidth = 0.02) +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp), 
                       breaks = seq(-limit_bp, limit_bp, length.out = 5)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(values =  c("#8d03b0", "#999999"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("AQR knockdown  ", "Control")) +
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    facet_grid(novel_type~RBP, 
               labeller = labeller(novel_type = c("novel donor" = "Novel donor", "novel acceptor" = "Novel acceptor"))) +
    theme_light() +
    custom_ggtheme
  
  
  distance_rectangle_longer <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 60), fill = "grey", color = "black") +
    geom_text(aes(x = 100, y = 30),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 30, ymax = 31), fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -100, y = 51),  size = 3, label = "intron") +
    theme_void()
  
  
  plot_aqr_long_distances / (distance_rectangle_longer) + patchwork::plot_layout(heights = c(8, 1))
  
  ggsave(file = paste0(args$figures_folder, "/", target_RBP, "_", limit_bp, "bp.png"), 
         width = 180, height = 90, dpi = 300, units = "mm")
  
}



##################################
## AUXILIARY FUNCTIONS
##################################

aux_get_annotated_introns_accurate_splicing <- function(data, cluster_type) {
  
  data %>% 
    filter(cluster == cluster_type) %>% 
    group_by(RBP) %>%
    filter(MSR_D == 0, MSR_A == 0) %>% 
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    group_by(ref_junID) %>%
    summarise(n = n()) %>%
    filter(n == 2) %>%
    ungroup() %>%
    pull(ref_junID) %>% 
    unique() %>%
    return()
  
}

aux_get_annotated_introns_inaccurate_splicing <- function(data, accurate_splicing_data, 
                                                      cluster_type, p_coding) {
  
  setdiff(data %>% 
            filter(cluster == cluster_type, protein_coding == p_coding) %>% 
            pull(ref_junID) %>% 
            unique, 
          accurate_splicing_data) %>%
    return()
  
}


GenerateMSRtests <- function(target_RBPs,
                             cluster_case = "case",
                             cluster_control = "control",
                             MSR_Table,
                             num_cores = 4,
                             file_output = "",
                             overwrite = F){
  
  
  if(!overwrite & file.exists(file_output)){
    MSR_tests <- readRDS(file_output)
    return(MSR_tests)
  }
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  MSR_tests <- foreach(i = seq(length(target_RBPs)), .packages = c("tidyverse")) %dopar%{
    # i <- 1
    project <- target_RBPs[i]
    
    MSR_project <- MSR_Table %>%
      dplyr::select(ref_junID,
                    case = paste0(project, "_", cluster_case),
                    control = paste0(project, "_", cluster_control)) %>%
      tidyr::drop_na() %>%
      tidyr::pivot_longer(cols = c("case", "control"), names_to = "group", values_to = "MSR") %>%
      dplyr::arrange(ref_junID)
    
    #MSR_case <- MSR_project %>% dplyr::filter(group == "case") %>% dplyr::pull(MSR)
    #MSR_control <- MSR_project %>% dplyr::filter(group == "control") %>% dplyr::pull(MSR)
    
    wilcox_test <- rstatix::wilcox_test(MSR_project, MSR ~ group, paired = TRUE, alternative = "greater")
    wilcox_effsize <- rstatix::wilcox_effsize(MSR_project, MSR ~ group, paired = TRUE, alternative = "greater")
    
    tibble(target_gene = project,
           statistical_test = NA,
           H0 = NA,
           H1 = NA,
           p.value = wilcox_test$p,
           effect_size = wilcox_effsize$effsize,
           magnitude = wilcox_effsize$magnitude)
  } %>% dplyr::bind_rows()
  parallel::stopCluster(cl)
  
  if(file_output != ""){
    MSR_tests %>% saveRDS(file_output)
  }
  
  return(MSR_tests)
}

##################################
## CALLS
##################################

# data_main_figure5_c()
