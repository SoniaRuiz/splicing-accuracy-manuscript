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

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/splicing_accuracy_manuscript_ENCODE_figures.R")


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
dependencies_folder <- paste0(base_folder, "/dependencies/")

results_folder <- file.path(base_folder, "results", project_name, gtf_version, "_paper_review/results")
dir.create(file.path(results_folder), recursive = TRUE, showWarnings = F)

figures_folder <- file.path(base_folder, "results", project_name, gtf_version, "_paper_review/figures")
dir.create(file.path(figures_folder), recursive = TRUE, showWarnings = F)


## QUERY MASTER TABLES 

query = paste0("SELECT * FROM 'metadata'")
master_metadata <- dbGetQuery(con, query) %>% as_tibble()
all_projects <- master_metadata$SRA_project %>% unique

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
  
  file_name <- paste0(figures_folder, "/avg_novel_jnx_per_tissue")
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
  
  file_name <- paste0(results_folder, "/common_introns_details_",
                      paste(RBPs,collapse = "_"), "_experiments_", 
                      paste(required_clusters,collapse = "_"), ".rds")
  saveRDS(object = common_introns_all_experiments, file = file_name)
  
  return(common_introns_all_experiments)
}

get_anchor_effect_control_knockdown <- function() {
  
  ## Load common introns
  
  required_clusters <- "control"
  target_genes = "all"
  
  common_introns_path <- paste0(results_folder, "/common_introns_details_",
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters,collapse = "_"), ".rds")
  
  if( file.exists(common_introns_path) ) {
    
    message("Loading common introns...")
    common_introns_all_experiments <- readRDS(file = common_introns_path)
    
  } else {
    ## Load common introns
    message("Calculating common introns...")
    common_introns_all_experiments <- get_common_introns_across_experiments(RBPs = "all",
                                                                            required_clusters = required_clusters)
  }

  common_introns_all_experiments %>%
    distinct(RBP) %>%
    dplyr::count()
  
  common_introns_all_experiments %>%
    group_by(RBP) %>%
    distinct(ref_junID) %>%
    dplyr::count() %>%
    ungroup()
  
  common_introns_all_experiments_w_novel <- common_introns_all_experiments %>%
    left_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_type, distance),
              by = "novel_junID")
  
  #################################################
  ## PROPORTION OF UNIQUE JUNCTIONS
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
    
    return(data.frame(RBP = project_id,
                      annotated_prop = annotated_prop,
                      donor_prop = donor_prop,
                      acceptor_prop = acceptor_prop))
    
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
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_x_discrete( breaks = c( "acceptor", "donor"),
                      labels = c( "Novel Acceptor", "Novel Donor")) +
    scale_fill_manual(values = c( "#35B779FF", "#8d03b0"),
                      breaks = c( "donor", "acceptor"),
                      labels = c( "Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1)) +
    custom_ggtheme + 
    theme(legend.position = "right")
  
  
  #################################################
  ## CUMMULATIVE NUMBER OF READS
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
    
    return(data.frame(RBP = project_id,
                      annotated = annotated_prop,
                      donor = donor_prop,
                      acceptor = acceptor_prop))
    
  })
  
  df_proportion_reads_tidy <- df_proportion_reads  %>%
    tidyr::gather(key = "type", value = "prop", -RBP) %>%
    filter(type != "annotated") %>%
    mutate(prop = prop * 100)
  
  df_proportion_reads_tidy$type = factor(df_proportion_reads_tidy$type, 
                                         levels = c("donor", "acceptor"))
  
  
  reads_violin_plot <- ggplot(df_proportion_reads_tidy, aes(type, prop, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point(size = 0.5, colour = "#333333") +
    geom_line( aes(group = RBP), colour = "#333333",linewidth = .2  )  +
    theme_light() +
    ylab("% cumulative read counts") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_x_discrete( breaks = c( "acceptor", "donor"),# "annotated_intron"),
                      labels = c( "Novel Acceptor", "Novel Donor")) +
    scale_fill_manual(values = c( "#35B779FF", "#8d03b0"),
                      breaks = c( "donor", "acceptor"),
                      labels = c( "Novel Donor", "Novel Acceptor")) + 
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))+
    custom_ggtheme + 
    theme(legend.position = "right")
  
  
  
  ggpubr::ggarrange(unique_jxn_violin_plot,
                    reads_violin_plot,
                    common.legend = T,
                    #labels = c("a", "b"),
                    align = "h",
                    ncol = 2,
                    nrow = 1) +
    theme( plot.margin = margin(t = -5,
                                r = -5,
                                b = -10,
                                l = -5),
           legend.box.margin=margin(b = -10))
  
  
  ggsave(file = paste0(figures_folder , "/anchor_RBP_unique_jxn_reads.png"),
         width = 130, height = 50, dpi = 300, units = "mm")
  
  #################################################
  ## MODULO
  #################################################
  
  df_modulo_experiments <- common_introns_all_experiments_w_novel %>%
    group_by(RBP) %>%
    distinct(novel_junID, .keep_all = T) %>%
    filter(abs(distance) <= 100) %>% 
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " ")) %>%
    mutate(modulo = abs(distance) %% 3) %>%
    ungroup
  
  df_modulo_experiments <- df_modulo_experiments %>% 
    group_by(RBP,modulo) %>%
    summarise(n = n()) %>%
    mutate(freq = (n / sum(n))*100) 
  
  df_modulo_experiments$modulo = factor(df_modulo_experiments$modulo, 
                                        levels = c( "0", "1", "2"))
  
  ggplot(df_modulo_experiments, aes(x = freq, y = modulo)) +
    ggridges::geom_density_ridges_gradient() +
    ylab("Modulo3 of the distance") +
    xlab("% of novel junctions") +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(30,35,40),
                       labels = c("30%","35%","40%")) +
    scale_y_discrete(expand = c(0,0.5,1,0))
  
  
  ggsave(file = paste0(figures_folder , "/anchor_RBP_modulo.png"),
         width = 50, height = 60, dpi = 300, units = "mm")
  
  
  
  
  #################################################
  ## DISTANCES
  #################################################
  
  
  limit_bp = 30
  
  common_introns_all_experiments_w_protein <- common_introns_all_experiments_w_novel %>%
    #filter(RBP == "U2AF2") %>%
    left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding),
              by = "ref_junID")
  
  plot_distances_PC <- ggplot(common_introns_all_experiments_w_protein %>%
                                filter(protein_coding %in% c(0,100)) %>%
                                mutate(protein_coding_label = ifelse (protein_coding == 0, "Noncoding", "Protein-coding")) %>%
                                mutate(protein_coding_label = protein_coding_label %>% as.factor()) %>%
                                mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " ")) %>%
                                distinct(novel_junID, .keep_all = T) %>% 
                                filter(abs(distance) < limit_bp) ) + 
    geom_histogram(aes(x = distance, fill = fct_rev(novel_type)),
                   bins = 60, 
                   binwidth = 1, 
                   position = "identity", 
                   alpha = 1, 
                   color = "black", 
                   linewidth = 0.1) +
    scale_x_continuous(breaks = seq(-limit_bp, limit_bp, length.out = 5)) + 
    # ggsci::scale_fill_npg()+
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("Novel donor", "Novel acceptor")) +
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    facet_grid(fct_rev(protein_coding_label)~fct_rev(novel_type),
               scales = "free_y") +
    
    guides(fill = guide_legend(title = "Sample type: ", ncol = 1, nrow = 2 )) +
    theme_light() +
    custom_ggtheme + 
    theme( legend.position = "none",
                            
                            plot.margin = margin(t = 0,
                                                 r = -5,
                                                 b = 0,
                                                 l = -5),
                            legend.box.margin=margin(b = -10))
  plot_distances_PC
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 60), fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 33),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 30, ymax = 31), fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15, y = 48),  size = 3, label = "intron") +
    theme_void()
  
  plot_distances_PC / (distance_rectangle +  distance_rectangle) + patchwork::plot_layout(heights = c(8, 1))
  
  
  ggsave(file = paste0(figures_folder , "/anchor_RBP_distances.png"),
         width = 125, height = 60, dpi = 300, units = "mm")
  
  
  
  
  #################################################
  ## DELTA MES
  #################################################
  
  
  ## DELTA MES VALUES 
  
  RBP_delta_MES <- common_introns_all_experiments_w_novel %>%
    filter(RBP == "U2AF2") %>%
    left_join(y = master_introns %>% dplyr::select(ref_junID, ref_mes5ss, ref_mes3ss),
              by = "ref_junID")%>%
    left_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_mes5ss, novel_mes3ss),
              by = "novel_junID") %>% 
    distinct(novel_junID, .keep_all = T) %>%
    mutate(delta_ss5score = ref_mes5ss - novel_mes5ss,
           delta_ss3score = ref_mes3ss - novel_mes3ss) %>%
    dplyr::select(novel_junID, ref_junID, novel_type,  
                  delta_ss5score, delta_ss3score, RBP) %>%
    drop_na() %>%
    dplyr::select(delta_ss5score, delta_ss3score, novel_type) %>%
    gather(key = delta_type, value = deltaMES) %>%
    mutate(delta_type = ifelse(delta_type == "delta_ss5score", "MES Donor", "MES Acceptor")) %>%
    mutate(delta_type = delta_type %>% as.factor()) %>%
    mutate(deltaMES = deltaMES %>% as.double()) %>%
    filter(deltaMES != 0) 
  
  
  
  ## PREPARE DATA FOR THE DENSITY PLOT
  
  delta_mes <- ggplot(data = RBP_delta_MES )  +
    geom_density(aes(x = deltaMES, fill = delta_type), alpha = 0.9) +
    geom_vline(xintercept = 0) +
    ggforce::facet_col(~fct_rev(delta_type), strip.position = "right") +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("MES Donor", "MES Acceptor"),
                      labels = c("Novel Donor", "Novel Acceptor")) +
    labs(x = "Delta MES", y = "Density") +
    theme_light() +
    custom_ggtheme + 
    theme( plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
           legend.box.margin = margin(b = -11) ) + 
    theme( legend.position = "none")
  
  delta_mes
  ggsave(file = paste0(figures_folder, "/anchor_RBP_deltaMES.png"),
         width = 45, height = 50, dpi = 300, units = "mm")
  
  
  
  
  #################################################
  ## MSR
  #################################################
  
  
  
  df_all_introns_tidy <- common_introns_all_experiments_w_protein %>%
    filter(RBP == "U2AF1") %>%
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
  
  df_all_introns_tidy$MSR_type = factor(df_all_introns_tidy$MSR_type, 
                                        levels = c("MSR_D", "MSR_A"))
  
  plot1 <- ggplot(data = df_all_introns_tidy) + 
    geom_bar(aes(x = type_label, fill = MSR_type), position = "dodge") +
    ggplot2::labs(x = "", y = "Annotated introns")+
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    #facet_grid(~target_gene) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1)) + 
    theme( plot.margin = margin(0,10,0,0),
           legend.box.margin=margin(b = -5) )
  
  
  plot1
  
  
  
  ## 2. ZOOMED PLOT VERSION
  
  df_all_introns_tidy$percentile_group = factor(df_all_introns_tidy$percentile_group,
                                                levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  plot2 <- ggplot(data = df_all_introns_tidy %>%
                    filter(percentile_group != "0")) +
    geom_bar(aes(x = percentile_group, fill = MSR_type),position = "dodge")+
    #ggtitle("All annotated introns") +
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    ggforce::facet_zoom(ylim=c(0,1000), split = TRUE) +
    #ylim(c(0,23000))+
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    #facet_grid(~target_gene) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2,
                               nrow = 1))
  
  
  ggpubr::ggarrange(plot1,
                    plot2,
                    labels = c("", ""),
                    nrow = 2, 
                    ncol = 1,
                    heights = c(1,1.5),
                    common.legend = T)
  
  ggsave(file = paste0(figures_folder,"/anchor_RBP_MSR.png"),
         width = 135, height = 70, dpi = 300, units = "mm")
  
  
  ########################################################
  ## MSR BIOTYPE
  ########################################################
  
  
  
  ## TIDY DATAFRAME
  df_biotype_result_tidy <- common_introns_all_experiments_w_protein %>%
    filter(RBP == "U2AF1") %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(biotype = ifelse(protein_coding == 100, "PC", "non PC"))  %>%
    distinct(ref_junID, .keep_all = T) %>%
    group_by(ref_junID) %>%
    mutate(mean_coverage = (sum(ref_sum_counts)/sum(ref_n_individuals)) %>% log10()) %>%
    ungroup()
  
  df_biotype_result_tidy %>%
    distinct(ref_junID, .keep_all = T) %>%
    dplyr::count(biotype)
  
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
  if ( !identical( data_combined$ref_junID %>% sort(), 
                   df_biotype_result_tidy %>% 
                   distinct(ref_junID,.keep_all = T) %>% 
                   pull(ref_junID) %>% 
                   sort()
  ) ) {
    print("ERROR!")
    break;
  }
  
  
  if ( !file.exists(paste0(results_folder, "/MSR_subsample.rds")) ) {
    ## Subsampling introns to control by similarity in mean read coverage
    m.out <- MatchIt::matchit(biotype ~ mean_coverage, 
                              data = data_combined, 
                              distance = data_combined$mean_coverage,
                              method = "nearest", 
                              caliper = c(mean_coverage = 0.005), 
                              std.caliper = FALSE)
    
    
    subsample <- MatchIt::match.data(m.out)
    subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)
    
    
    saveRDS(object = subsample,
            file = paste0(results_folder, "/anchor_MSR_subsample.rds"))
  } else {
    subsample <- readRDS(file = paste0(results_folder, "/anchor_MSR_subsample.rds"))
  }
  
  
  plot_AS <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = biotype), alpha = 0.9) +
    ggtitle("After subsampling") +
    theme_light() +
    custom_ggtheme +
    xlab("log10 mean expresion level")+
    ggsci::scale_fill_npg(name = "Transcript biotype: ")
  
  
  ggpubr::ggarrange(plot_BS,
                    plot_AS,
                    labels = c("a", "b"),
                    common.legend = T)
  
  file_name <- paste0(figures_folder, "/anchor_MSR_subsampling")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
  
  ## Analysing differences in MSR_D and MSR_A
  
  ## subset introns from lncRNA transcripts after subsampling
  subsample_lncRNA <- subsample %>%
    dplyr::filter(biotype=="non PC")
  ## subset introns from protein-coding transcripts after subsampling
  subsample_protein_coding <- subsample %>%
    dplyr::filter(biotype=="PC")
  ## QC
  if ( intersect(subsample_protein_coding$ref_junID, subsample_lncRNA$ref_junID) %>% length() > 0 )  {
    print("ERROR! some introns have been categorised as lncRNA and protein-coding.")
  }
  
  
  df_introns_biotype <- subsample %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC")) 
  
  any(df_introns_biotype %>%
        pull(ref_junID ) %>% 
        duplicated())
  
  
  df_introns_biotype$biotype = factor(df_introns_biotype$biotype, 
                                      levels = c("non PC","PC"))
  
  df_introns_biotype <- df_introns_biotype %>%
    dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)
  
  
  df_introns_biotype$MSR_type = factor(df_introns_biotype$MSR_type, 
                                       levels = c("MSR_A", "MSR_D"))
  
  
  print(paste0(Sys.time(), " - ", df_introns_biotype %>% nrow(), " - introns after tidying!"))
  
  
  df_introns_biotype %>% 
    filter(biotype == "PC") %>%
    pull(MSR) %>%
    summary()
  
  df_introns_biotype %>% 
    filter(biotype == "non PC") %>%
    pull(MSR) %>%
    summary()
  
  df_introns_biotype_tidy <- df_introns_biotype %>%
    group_by(biotype, MSR_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    mutate(percentile_group = case_when(MSR == 0 ~ "0",
                                        MSR > 0 & MSR <= 0.2 ~ "(0,.2]",
                                        MSR > 0.2 & MSR <= 0.4 ~ "(.2,.4]",
                                        MSR > 0.4 & MSR <= 0.6 ~ "(.4,.6]",
                                        MSR > 0.6 & MSR <= 0.8 ~ "(.6,.8]",
                                        MSR > 0.8 & MSR <= 1 ~ "(.8,1]"))
  
  df_introns_biotype_tidy$MSR_type = factor(df_introns_biotype_tidy$MSR_type, 
                                            levels = c("MSR_D","MSR_A"))
  df_introns_biotype_tidy$percentile_group = factor(df_introns_biotype_tidy$percentile_group, 
                                                    levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  df_introns_biotype_tidy$biotype = factor(df_introns_biotype_tidy$biotype, 
                                           levels = c("PC","non PC"))
  
  MSR_label <- c(
    MSR_D = "MSR Donor",
    MSR_A = "MSR Acceptor"
  )
  
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
  

  file_name <- paste0(figures_folder, "/anchor_MSR_biotype")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                  width = 70, height = 60, units = "mm", dpi = 300)
  
  
  plotMSR_donor_zoomed <- ggplot(data = df_introns_biotype_tidy) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    ggforce::facet_col(~MSR_type, strip.position = "top", labeller = labeller(MSR_type=MSR_label))+
    #ggtitle("MSR Donor") +
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
  file_name <- paste0(figures_folder, "/anchor_MSR_biotype_zoomed")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 70, height = 60, units = "mm", dpi = 300)
  
}


get_effect_size_data_all_RBPs <- function() {
  
  ###############################
  ## LOAD METADATA
  ###############################
  
  if (analysis_type=="shRNA") {
    
    metadata_RBPs <- master_metadata %>% 
      pivot_longer(c("Splicing regulation", "Spliceosome", "Exon Junction Complex", 
                     "NMD", "Novel RBP", "RNA modification"), names_to = "Category") %>%
      filter(value == 1) %>%
      dplyr::select(-value) %>%
      distinct(target_gene, sample_id, .keep_all = T)
    
    metadata_RBPs$Category %>% unique
  } else {
    metadata_RBPs <- master_metadata
  }
  
 
  
  
  ##########################
  ## GET THE EFFECT SIZES
  ##########################
  
  ## Load common introns ---------------------------------------
  
  target_genes = "all"
  required_clusters <- c("case", "control")
  
  common_introns_path <- paste0(results_folder, "/common_introns_details_", 
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  
  if ( !file.exists(common_introns_path)  ) {
    
    ## Get common introns
    message("Calculating common introns...")
    common_introns_all_experiments <- get_common_introns_across_experiments(RBPs = target_genes,
                                                                            required_clusters = c("case", "control"))
    
    common_introns_all_experiments %>% dplyr::count(RBP)
    
  } else {
    
    message("Loading common introns across '", target_genes, "' experiments...")
    common_introns_all_experiments <- readRDS(file = common_introns_path)
    
  }
  
  
  ## Run the wilcox test on MSR_D --------------------------
  
  if ( !file.exists(paste0(results_folder, "/ENCODE_effectsize_MSRD.rds")) ) {
    
    source(paste0(base_folder, "/knockdown_analysis/helper_functions.R"))
    
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
    
    MSR_D_tests <- generateMSRtests(target_RBPs = all_projects, 
                                    MSR_Table = MSR_D, 
                                    file_output = paste0(results_folder, "/ENCODE_effectsize_MSRD.rds"), 
                                    overwrite = T,
                                    num_cores = 8)
    
    ## Add the categories
    if (analysis_type=="shRNA") {
      MSR_D_tests <- MSR_D_tests %>%
        left_join(y = metadata_RBPs %>% dplyr::select(target_gene, Category),
                  by = "target_gene") 
    }
    
    ## Add bonferroni correction
    MSR_D_tests <- MSR_D_tests %>% mutate(FDR = NA, .after = p.value)
    MSR_D_tests$FDR <- p.adjust(MSR_D_tests$p.value, method = "fdr")
    
    ## Save results
    saveRDS(object = MSR_D_tests, file = paste0(results_folder, "/ENCODE_effectsize_MSRD.rds"))
    write_csv(x = MSR_D_tests %>%
                distinct(target_gene, .keep_all=T) %>%
                mutate(statistical_test = "Wilcoxon Rank text: rstatix::wilcox_test(data, formula, paired = TRUE, correct = TRUE, alternative = 'greater')",
                       H0 = "The MSR_D observations in case and control samples are symmetric about their median value.",
                       H1 = "The MSR_D observations in case samples are greater at their median value than the MSR_D observations in control samples."),
              file = paste0(results_folder, "/ENCODE_effectsize_MSRD.csv"), col_names = T)
    
  } else {
    
    message("Loading 'ENCODE_effectsize_MSRD.rds' file ...")
    MSR_D_tests <- readRDS(file = paste0(results_folder, "/ENCODE_effectsize_MSRD.rds"))
  }
  
  
  
  if ( !file.exists(paste0(results_folder, "/ENCODE_effectsize_MSRA.rds")) ) {
    
    source(paste0(base_folder, "/knockdown_analysis/helper_functions.R"))
    
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
    
    MSR_A_tests <- generateMSRtests(target_RBPs = all_projects, 
                                    MSR_Table = MSR_A, 
                                    file_output = paste0(results_folder, "/ENCODE_effectsize_MSRA.rds"),  
                                    overwrite = T, 
                                    num_cores = 8)
    
    ## Add the categories
    if (analysis_type=="shRNA") {
      MSR_A_tests <- MSR_A_tests %>%
        left_join(y = metadata_RBPs %>% dplyr::select(target_gene, Category),
                  by = "target_gene")
    }
    
    ## Add bonferroni correction
    MSR_A_tests <- MSR_A_tests %>% distinct(target_gene, .keep_all=T) %>% mutate(FDR = NA, .after = p.value)
    MSR_A_tests$FDR <- p.adjust(MSR_A_tests$p.value, method = "fdr")
    
    ## Save results
    saveRDS(object = MSR_A_tests, file = paste0(results_folder, "/ENCODE_effectsize_MSRA.rds"))
    write_csv(x = MSR_A_tests %>%
                distinct(target_gene, .keep_all=T) %>%
                mutate(statistical_test = "Wilcoxon Rank text: rstatix::wilcox_test(data, formula, paired = TRUE, correct = TRUE, alternative = 'greater')",
                       H0 = "The observations MSR_A in case samples vs control samples are symmetric about their median value.",
                       H1 = "The observations MSR_A in case samples are greater at their median value than the observations MSR_A in control samples."),
              file = paste0(results_folder, "/ENCODE_effectsize_MSRA.csv"))
  } else {
    message("Loading 'ENCODE_effectsize_MSRA.rds' file ...")
    MSR_A_tests <- readRDS(file = paste0(results_folder, "/ENCODE_effectsize_MSRA.rds"))
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
  
  MSR_D_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05) %>%
    distinct(target_gene, .keep_all=T) 
  
  MSR_D_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05) %>%
    distinct(target_gene, .keep_all=T) %>%
    pull(FDR) %>%
    summary
  
  
  
  ## MSR_A
  
  ## increase in mis-splicing rates in samples with gene knockdowns compared to untreated controls for 90% of the 54 genes considered
  
  (MSR_A_tests %>%
      filter(FDR <= 0.05) %>%
      distinct(target_gene) %>%
      nrow() * 100) / ( MSR_A_tests %>%
                          distinct(target_gene)%>%
                          nrow() )
  
  MSR_A_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05) %>%
    distinct(target_gene, .keep_all=T)
  
  MSR_A_tests %>%
    distinct(target_gene, .keep_all = T) %>%
    filter(FDR <= 0.05) %>%
    distinct(target_gene, .keep_all=T) %>%
    pull(FDR) %>%
    summary
  
 
  
  
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

  
plot_effect_size_data_all_RBPs <- function() {
  
  
  ## Run the wilcox test on MSR_A --------------------------
  
  
  message("Loading 'ENCODE_effectsize_MSRD.rds' file ...")
  MSR_D_tests <- readRDS(file = paste0(results_folder, "/ENCODE_effectsize_MSRD.rds"))
  
  message("Loading 'ENCODE_effectsize_MSRA.rds' file ...")
  MSR_A_tests <- readRDS(file = paste0(results_folder, "/ENCODE_effectsize_MSRA.rds"))
  
  
  ##########################
  ## TIDY EFFECT SIZES
  ##########################
  
  # Combine both MSR_A and MSR_D -------------------
  MSR_combined = rbind(MSR_A_tests %>% mutate(MSR_type = "MSR_A"), 
                       MSR_D_tests %>% mutate(MSR_type = "MSR_D")) %>% 
    filter(FDR <= 0.05 | target_gene == "SAFB2")## 'SAFB2' is a novel RBP that will be used as control
  MSR_combined$Category[is.na(MSR_combined$Category)] <- "Novel_RBP"
  
  max_genes <- 40
  
  
  # Filter the Splicing regulation category
  filter_splicing_regulation <- MSR_combined %>% 
    dplyr::select(-statistical_test, -H0, -H1) %>% 
    arrange(-effect_size) %>%
    #filter(Category == "Splicing_regulation") %>%
    distinct(target_gene, .keep_all = T) %>%
    head(max_genes) %>%
    pull(target_gene)
  
  MSR_graph_data <- MSR_combined %>%
    filter(target_gene %in% c(filter_splicing_regulation) ) %>% #| Category != c("Splicing regulation") ) %>%
    arrange(-effect_size) %>%
    mutate(target_gene = factor(target_gene, levels = .$target_gene %>% unique))
  
  if (analysis_type == "shRNA") {
    MSR_graph_data <- MSR_graph_data %>%
      mutate(# Use factors to sort the graph
        Category = factor(Category, levels = c("Splicing regulation", 
                                                "Spliceosome", 
                                                "NMD", 
                                                "Exon Junction Complex", 
                                                "Novel RBP")))  # Use factors to sort the graph)
    print(MSR_graph_data$Category %>% unique())
  }
  
  
 
  
  ##########################
  ## METADATA KEFF
  ##########################
  
  project_path <- "./ENCODE_SR/ENCODE_Metadata_Extraction/Metadata_results_experiments/"
  metadata_kEff_path <- paste0(project_path, "metadata_WB_kEff.tsv")
  metadata_kEff <- readr::read_delim(metadata_kEff_path, show_col_types = F) %>%
    mutate(kEff_text = ifelse(is.na(kEff_avg), kEff_avg, paste0(round(kEff_avg), "%"))) %>%
    mutate(kEff_text = kEff_text %>% as.factor())


  MSR_graph_data <- MSR_graph_data %>%
    left_join(y = metadata_kEff,
              by = "target_gene")
  MSR_graph_data %>% head()
  
  
  ##########################
  ## PLOT
  ##########################
  
  # Plot the graph
  plot_effectsize <- ggplot(MSR_graph_data, # %>%
                              #dplyr::filter(!(target_gene %in% c("PABPC1","RPS19","RPS3A", "EIF4G1"))), 
                            aes(x = target_gene, y = effect_size)) + 
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
    theme_light() 
  

  if (analysis_type == "shRNA") {
    plot_effectsize <- plot_effectsize +
      ggforce::facet_row(facets = vars(Category), 
                       scales = "free_x", space = "free",
                       #labeller = labeller(Category = category_labels),
                       drop = T,
                       shrink = T) 
  }
  plot_effectsize +
  custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.box = "horizontal") +
    guides(colour = guide_legend(override.aes = list(fill = '#999999'),
                                 title = NULL,
                                 label.position = "bottom",
                                 order = 1))
  
   plot_effectsize +
   ggnewscale::new_scale_fill() +
   geom_tile(stat = "identity",
             aes(y = 0.68, fill = kEff_avg, color = "No data\navailable"),
             linewidth = 0.5,
             width = 1,
             height = 0.045) +
   geom_text(aes(y = 0.76, label = kEff_text),
             color = "black",
             size = 2.5) +
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
     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
           legend.box = "horizontal") +
     guides(colour = guide_legend(override.aes = list(fill = '#999999'),
                                  title = NULL,
                                  label.position = "bottom",
                                  order = 1))
  
  # Save the graph
  ggsave(file = paste0(figures_folder, "/Effect_size_combined_top", max_genes, ".png"), 
         width = 180, height = 90, units = "mm", dpi = 300)
  
  
  
}


plot_AQR_U2AF2_acceptor <- function() {
  
  
  ###########################################
  ## LOAD COMMON INTRONS AND NOVEL FOR AQR 
  ## AND U2AF2
  ###########################################
  
  target_RBPs <- c("AQR", "U2AF2")
  
  required_clusters <- c("case", "control")
  
  common_introns_path <- paste0(results_folder, "/common_introns_details_", 
                                paste(target_RBPs, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  ## LOAD COMMON INTRONS
  
  if( file.exists(common_introns_path) ) {
    
    common_introns <- readRDS(file = common_introns_path)
    
  } else {
    ## Load common introns
    common_introns <- get_common_introns_across_experiments(RBPs = target_RBPs,
                                                            required_clusters = c("case", "control"))
  }
  
  RBP_intron <- common_introns  %>%
    mutate(cluster = ifelse(cluster == "case", "gene knockdown", "control")) %>%
    mutate(cluster = factor(cluster, levels = c("control", "gene knockdown"))) %>%
    as_tibble() %>%
    filter(RBP %in% target_RBPs) %>%
    inner_join(y = master_introns %>% dplyr::select(ref_junID, ref_mes5ss, ref_mes3ss),
               by = "ref_junID")
  
  
  RBP_intron %>%
    dplyr::group_by(RBP, cluster) %>%
    distinct(ref_junID) %>%
    dplyr::count() %>%
    ungroup()
  
  
  ## LOAD NOVEL JUNC FROM COMMON INTRONS
  
  RBP_novel <- RBP_intron %>%
    left_join(y = master_novel_junctions %>% 
                dplyr::select(novel_junID, novel_type, distance, novel_mes5ss, novel_mes3ss),
              by = "novel_junID")
  
  
  RBP_novel %>%
    dplyr::group_by(RBP, cluster) %>%
    distinct(novel_junID) %>%
    dplyr::count() %>%
    ungroup()
  
  
  
  RBP_novel_delta <- RBP_novel %>%
    dplyr::group_by(RBP, cluster) %>%
    distinct(novel_junID, .keep_all = T) %>%
    ungroup() %>%
    mutate(delta_ss5score = ref_mes5ss - novel_mes5ss) %>%
    mutate(delta_ss3score = ref_mes3ss - novel_mes3ss) 
  
  
  
  #################################################
  ## STATS FOR THE PAPER
  #################################################

  
  ## Delta 3 MES scores of novel acceptor junctions were significantly lower indicating that weaker acceptor splice sites were being used 
  
  ## AQR
  wilcox.test(x = RBP_novel_delta %>% filter(RBP == "AQR", cluster == "gene knockdown") %>% distinct(novel_junID, .keep_all = T) %>% pull(delta_ss3score),
              y = RBP_novel_delta %>% filter(RBP == "AQR", cluster == "control") %>% distinct(novel_junID, .keep_all = T) %>% pull(delta_ss3score),
              paired = F,
              alternative = "greater")
  
  rstatix::wilcox_effsize(data = RBP_novel_delta  %>% 
                            filter(RBP == "AQR") %>% 
                            distinct(novel_junID, .keep_all = T) %>% 
                            mutate(cluster = cluster %>% as.factor()),
                          formula = delta_ss3score ~ cluster)
  
  
  ## U2AF2
  wilcox.test(x = RBP_novel_delta %>% filter(RBP == "U2AF2", cluster == "gene knockdown") %>% pull(delta_ss3score),
              y = RBP_novel_delta %>% filter(RBP == "U2AF2", cluster == "control") %>% pull(delta_ss3score),
              paired = F,
              alternative = "greater")
  
  rstatix::wilcox_effsize(data = RBP_novel_delta  %>% 
                            filter(RBP == "U2AF2") %>%
                            distinct(novel_junID, .keep_all = T) %>% 
                            mutate(cluster = cluster %>% as.factor()),
                          formula = delta_ss3score ~ cluster,
                          paired = F)
  
  
  
  
  #################################################
  ## PLOT DELTA MES ONLY ACCEPTOR
  #################################################
  
  ## DISTANCES
  
  limit_bp = 30
  
  
  RBP_novel_acceptor_delta <- RBP_novel_delta %>%
    filter(novel_type == "novel_acceptor",
           RBP %in% c("AQR", "U2AF2")) %>%
    dplyr::select(delta_ss3score, cluster, RBP) %>%
    drop_na() %>%
    group_by(RBP, cluster) %>%
    mutate(medianMSR = delta_ss3score %>% median()) %>% 
    ungroup()  %>%
    group_by(RBP) %>%
    mutate(p.value = format(x = wilcox.test(delta_ss3score ~ cluster)$p.value,
                            digits = 2, scientific = T)) %>%
    mutate(p.value = ifelse((p.value %>% as.double()) < 0.001, 0.001, p.value))%>%
    ungroup() %>%
    mutate(RBP = factor(x = RBP, levels = target_RBPs)) %>%
    mutate(RBP = factor(x = RBP, levels = target_RBPs)) %>%
    mutate(p.value = ifelse(p.value == "0e+00", "2.2e-100",p.value)) %>%
    dplyr::select(delta_ss3score, cluster, medianMSR, RBP, p.value)
  
  
  delta_mes_acceptor <- ggplot(data = RBP_novel_acceptor_delta)  +
    geom_density(aes(x = delta_ss3score, fill = cluster), 
                 alpha = 0.8, linewidth = 0.3, 
                 color = "black") +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept=medianMSR, colour=cluster),
               linetype="dashed", linewidth=0.9) +
    facet_wrap(vars(RBP)) +
    geom_text(data = RBP_novel_acceptor_delta %>% group_by(RBP) %>% dplyr::select(p.value) %>% distinct(p.value), 
              aes(label = paste0("P<", p.value)), 
              x=25, y=0.08, family= "Arial", 
              size = 3.5, 
              colour = "#333333") +
    scale_fill_manual(values = c("#64037d", "#999999"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("KD", "Wild Type")) + 
    scale_colour_manual(values = c("#64037d", "#333333"),
                        breaks = c("gene knockdown", "control"),
                        labels = c("KD", "Wild Type")) + 
    labs(x = "Delta MES Acceptor", y = "Density") +
    theme_light() +
    guides(fill = guide_legend(title = "Sample type: ",
                               ncol = 2, nrow = 1 ),
           colour = guide_legend(title = "Median Delta MES: ",
                                 ncol = 2, nrow = 1 )) +
    custom_ggtheme 
  
  
  delta_mes_acceptor + 
    theme(legend.box = "horizontal") +
    theme( legend.margin = margin(-0.1,0.5,0,0, unit="cm"),
           legend.box.margin=margin(l = -5, b = -5))
  
  ggsave(file = paste0(figures_folder,"/", target_RBPs[1], "_", target_RBPs[2], "_acceptor_deltaMES.png"), 
         width = 180, height = 60, dpi = 300, units = "mm")
  
  
  
  
  
  
  ##############################################################
  ## DISTANCES 200 BP (AQR)
  ##############################################################
  
  
  
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
  
  ggsave(file = paste0(figures_folder, "/", target_RBP, "_", limit_bp, "bp.png"), 
         width = 180, height = 90, dpi = 300, units = "mm")
  
  
}


plot_UPF1_UPF2_distances <- function () {
  
  
  #################################################
  ## LOAD COMMON INTRONS AND NOVEL JUNCTIONS
  #################################################
  
  ## As we are only comparing UPF1 and UPF2, we get the common introns between the two genes
  
  target_genes = c("UPF1","UPF2")
  required_clusters <- c("case", "control")
  
  ## Get annotated common introns -------
  
  common_introns_path <- paste0(results_folder, "/common_introns_details_", 
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  if( file.exists(common_introns_path) ) {
    
    message("Loading common introns...")
    common_introns_NMD <- readRDS(file = common_introns_path)
    
  } else {
    ## Load common introns
    common_introns_NMD <- get_common_introns_across_experiments(RBPs = target_genes,
                                                                required_clusters = c("case", "control"))
  }
  
  ## Add protein-coding info
  common_introns_NMD <- common_introns_NMD %>%
    inner_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding),
              by = "ref_junID")
  
  
  ## Get novel junctions from common introns -------
  
  NMD_novel <- common_introns_NMD %>% 
    filter(!is.na(novel_junID)) %>%
    left_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_type, distance),
              by = "novel_junID")
  
  NMD_novel %>% 
    dplyr::count(RBP, cluster)
  
  
  NMD_novel <- NMD_novel %>%
    mutate(type_colour = paste0(cluster, "_", novel_type)) %>%
    mutate(novel_type = ifelse(novel_type == "novel_donor", "Novel Donor", "Novel Acceptor")) %>%
    mutate(cluster = ifelse(cluster == "case","Knockdown (KD)", "Wild Type"))
  
  #################################################
  ## STATS FOR THE PAPER
  #################################################
  
  
  ## Same number of annotated introns
  common_introns_NMD %>% 
    group_by(RBP, cluster) %>%
    distinct(ref_junID) %>%
    dplyr::count() %>%
    ungroup()
  
  ## Different number of novel junctions
  common_introns_NMD %>% dplyr::count(RBP, cluster)
  
  
  
  
  #################################################
  ## PLOT DISTANCES FOR "UPF1" 
  #################################################
  
  limit_bp <- 30
  
  for (gene_name in target_genes) {
    
    # gene_name <- target_genes[1]
    
    PC_NMD <- ggplot(data = NMD_novel %>%
                               filter(RBP == gene_name, 
                                      abs(distance) <= limit_bp,
                                      protein_coding == 100)) + 
      geom_histogram(aes(x = distance, fill = type_colour),
                     bins = limit_bp * 2,
                     binwidth = 1,
                     position = "stack"
      ) +
      ggplot2::facet_grid(fct_rev(novel_type)~fct_rev(cluster)) +
      ggtitle(paste0(gene_name, " - Protein-coding")) +
      xlab("Distance (in bp)") +
      ylab("Unique novel junctions") +
      theme_light() +
      scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                         breaks = seq(from = (limit_bp * -1), 
                                      to = limit_bp,
                                      by = 3)) +
      scale_fill_manual(values = c("#9ce2c0","#35B779FF","#e99cfc", "#8c04ae"),
                        breaks = c("control_novel_donor", "case_novel_donor","control_novel_acceptor", "case_novel_acceptor"),
                        label = c("Donor - Wild Type", "Donor - KD","Acceptor - Wild Type", "Acceptor - KD")) + 
      guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
                                 #override.aes = list(size = 3),
                                 ncol = 4, nrow = 1 ))+
      custom_ggtheme 
    
  
    
    NPC_NMD <- ggplot(data = NMD_novel %>%
                        filter(RBP == gene_name,
                               protein_coding == 0,
                               abs(distance) <= limit_bp)) + 
      geom_histogram(aes(x = distance, fill = type_colour),
                     bins = limit_bp * 2, binwidth = 1, position = "stack"
      ) +
      ggplot2::facet_grid(fct_rev(novel_type)~fct_rev(cluster)) +
      ggtitle(paste0(gene_name, " - Noncoding")) +
      xlab("Distance (in bp)") +
      ylab("Unique novel junctions") +
      theme_light() +
      scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                         breaks = seq(from = (limit_bp * -1), to = limit_bp,  by = 3)) +
      scale_fill_manual(values = c("#9ce2c0","#35B779FF","#e99cfc", "#8c04ae"),
                        breaks = c("control_novel_donor", "case_novel_donor","control_novel_acceptor", "case_novel_acceptor"),
                        label = c("Donor - Wild Type", "Donor - KD","Acceptor - Wild Type", "Acceptor - KD")) + 
      guides(fill = guide_legend(title = NULL, ncol = 4, nrow = 1 )) +
      custom_ggtheme+
      theme(legend.position = "none")


    
    distance_rectangle <- ggplot() +
      geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
                fill = "grey", color = "black") +
      geom_text(aes(x = 15, y = 55),  size = 3, label = "exon") +
      geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
                fill = "grey", alpha = 1, color = "black") +
      geom_text(aes(x = -15, y = 70),  size = 3, label = "intron") +
      theme_void()
    
    
    PC_NMD <- PC_NMD / (distance_rectangle + distance_rectangle ) +  patchwork::plot_layout(heights = c(8, 1))
    PC_NMD
    
    NPC_NMD <- (NPC_NMD + ggpubr::rremove("legend")) / (distance_rectangle + distance_rectangle ) +  patchwork::plot_layout(heights = c(8, 1))
    NPC_NMD   
    
    
    
    ggpubr::ggarrange(PC_NMD, 
                      NPC_NMD , 
                      vjust=.5, 
                      font.label = list(size = 12),
                      ncol = 1,
                      nrow = 2)
    
    ggsave(file = paste0(figures_folder, "/", gene_name, "_knockdown_distances.png"), 
           width = 180, height = 180, dpi = 300, units = "mm")
  }
 
  
}


plot_UPF1_UPF2_phastCons17_MES <- function () {
  
  
  #################################################
  ## LOAD COMMON INTRONS AND NOVEL JUNCTIONS
  ## FROM UPF1 and UPF2
  #################################################
  
  ## As we are only comparing UPF1 and UPF2, we get the common introns between the two genes
  
  target_genes = c("UPF1","UPF2")
  required_clusters <- c("case", "control")
  
  ## Get annotated common introns -------
  
  common_introns_path <- paste0(results_folder, "/common_introns_details_", 
                                paste(target_genes, collapse = "_"), "_experiments_", 
                                paste(required_clusters, collapse = "_"), ".rds")
  
  if( file.exists(common_introns_path) ) {
    
    NMD_intron <- readRDS(file = common_introns_path)
    
  } else {
    ## Load common introns
    NMD_intron <- get_common_introns_across_experiments(RBPs = target_genes,
                                                                required_clusters = c("case", "control"))
  }
  
  ## Add information from the master INTRON table -----------------------------------

  ## Given that we are studying the NMD action, we only subset introns 
  ## from protein-coding transcripts
  
  NMD_intron %>% dplyr::group_by(RBP, cluster) %>%  distinct(ref_junID) %>% dplyr::count() %>% ungroup() 
  
  NMD_intron_w_protein_coding <- NMD_intron %>%
    left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding,
                                                    mean_phastCons17way5ss_100,mean_phastCons17way3ss_100,
                                                    ref_mes5ss, ref_mes3ss),
               by = "ref_junID") %>%
    filter(protein_coding == 100)
  
  
  ## ONLY 37615 annotated introns from protein-coding transcripts overlap the two NMD datasets
  
  NMD_intron_w_protein_coding %>% dplyr::group_by(RBP, cluster) %>%  distinct(ref_junID) %>% dplyr::count() %>% ungroup() 

  
  
  #################################################
  ## GET COMMON INTRONS THAT ARE MIS-SPLICED &
  ## ACCURATELY SPLICED IN CONTROL AND CASE EXPERIMENTS
  #################################################
  
  ## 1. CONTROL EXPERIMENTS
  
  ## Get annotated introns with accurate splicing across CONTROL NMD experiments
  
  control_accurate_splicing <- NMD_intron_w_protein_coding %>% 
    filter(cluster == "control") %>% 
    group_by(RBP) %>%
    filter(MSR_D == 0, MSR_A == 0) %>% 
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    group_by(ref_junID) %>%
    summarise(n = n()) %>%
    filter(n == 2) %>%
    ungroup() %>%
    pull(ref_junID) %>% 
    unique()
  
  control_accurate_splicing %>% length()
  
  
  
  ## Get introns with INACCURATE splicing across CONTROL NMD experiments
  
  control_mis_splicing <- setdiff(NMD_intron_w_protein_coding %>% 
                                    filter(cluster == "control", protein_coding == 100) %>% 
                                    pull(ref_junID) %>% 
                                    unique, 
                                  control_accurate_splicing) 
  
  control_mis_splicing %>% unique %>% length()
  
  
  ## 2. CASE EXPERIMENTS
  
  ## Get annotated introns with accurate splicing across CASE NMD experiments
  
  case_accurate_splicing <- NMD_intron_w_protein_coding %>% 
    filter(cluster == "case")%>% 
    group_by(RBP) %>%
    filter(MSR_D == 0, MSR_A == 0) %>% 
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    group_by(ref_junID) %>%
    summarise(n = n()) %>%
    filter(n == 2) %>%
    ungroup() %>%
    pull(ref_junID) %>% 
    unique
  case_accurate_splicing %>% length()

  
  ## Get annotated introns with INACCURATE splicing across CASE NMD experiments
  case_mis_splicing <- setdiff(NMD_intron_w_protein_coding %>% 
                                 filter(cluster == "case") %>% 
                                 pull(ref_junID) %>% unique, 
                               case_accurate_splicing) %>% unique
  
  case_mis_splicing %>% unique %>% length()
  

  
  #################################################
  ## DEFINE NMD-RESISTANT AND NMD-SENSITIVE
  ## INTRONS
  #################################################
  
  
  ## 1. NMD RESISTANT INTRONS ---------------------------------------------------------------------------
  
  
  NMD_resistant_introns <- setdiff(control_accurate_splicing, case_mis_splicing) 
  NMD_resistant_introns %>% length()
  
  ## Add master intron table info to calculate expression levels
  NMD_resistant <- NMD_intron_w_protein_coding %>%
    filter(ref_junID %in% NMD_resistant_introns) %>%
    mutate(NMD_type = "NMD Resistant",
           mean_expression = (ref_sum_counts/ref_n_individuals) %>% log10)  %>% 
    filter(cluster == "control")
  
  NMD_resistant %>% dplyr::count(RBP, cluster)
  NMD_resistant %>% dplyr::group_by(RBP, cluster) %>%  distinct(ref_junID) %>% dplyr::count() %>% ungroup() 
  
  
  ## 2. NMD SENSITIVE INTRONS ---------------------------------------------------------------------------
  
  NMD_sensitive_introns <- intersect(control_accurate_splicing, case_mis_splicing)
  NMD_sensitive_introns %>% length()
  
  ## Add master intron table info to calculate expression levels
  NMD_sensitive <- NMD_intron_w_protein_coding %>%
    filter(ref_junID %in% NMD_sensitive_introns)  %>%
    mutate(NMD_type = "NMD Sensitive",
           mean_expression = (ref_sum_counts/ref_n_individuals) %>% log10) %>% 
    filter(cluster == "control")
  
  NMD_sensitive %>% dplyr::count(RBP,cluster)
  NMD_sensitive %>% dplyr::group_by(RBP, cluster) %>%  distinct(ref_junID) %>% dplyr::count() %>% ungroup() 
  
  

  
  #################################################
  ## SUBSAMPLE TO GET ANNOTATED INTRONS WITH 
  ## SIMILAR EXPRESSION LEVELS
  #################################################
  
  data_combined <- rbind(NMD_resistant, NMD_sensitive) %>%
    distinct(ref_junID, .keep_all = T)
  
  data_combined %>%
    dplyr::count(NMD_type)
  
  ## Subsampling introns to control by similarity in mean read coverage
  m.out <- MatchIt::matchit(NMD_type ~ mean_expression, 
                            data = data_combined, 
                            distance = data_combined$mean_expression,
                            method = "nearest", 
                            caliper = c(mean_expression = 0.005), 
                            std.caliper = FALSE)
  subsample <- MatchIt::match.data(m.out)
  subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(NMD_type)
  
  
  subsample %>% dplyr::count(NMD_type)
  
  
  #################################################
  ## STATS FOR THE PAPER
  #################################################
  
  ## We found that mis-splicing events that were only evident in the context of NMD knockdown, were generated from annotated introns which had significantly lower phastCons17 and MES values at their 3’ss
  
  ## Test 5'ss - phastCons17
  
  wilcox.test(x = subsample %>%
                filter(NMD_type == "NMD Sensitive") %>%
                pull(mean_phastCons17way5ss_100),
              y = subsample %>%
                filter(NMD_type == "NMD Resistant") %>%
                pull(mean_phastCons17way5ss_100),
              alternative = "less",
              paired = T)
  
  ## Test 5'ss - MES
  
  wilcox.test(x = subsample %>%
                filter(NMD_type == "NMD Sensitive") %>%
                pull(ref_mes5ss),
              y = subsample %>%
                filter(NMD_type == "NMD Resistant") %>%
                pull(ref_mes5ss),
              alternative = "less",
              paired = T)
  
  
  
  ## Test 3'ss - phastCons17
  
  wilcox.test(x = subsample %>%
                filter(NMD_type == "NMD Sensitive") %>%
                pull(mean_phastCons17way3ss_100),
              y = subsample %>%
                filter(NMD_type == "NMD Resistant") %>%
                pull(mean_phastCons17way3ss_100),
              
              alternative = "less",
              paired = T)
  rstatix::wilcox_effsize(data = subsample %>% 
                            mutate(NMD_type = NMD_type %>% as.factor()),
                          formula = mean_phastCons17way3ss_100 ~ NMD_type,
                          paired = T)
  
  
  ## Test 3'ss - MES
  
  wilcox.test(x = subsample %>%
                filter(NMD_type == "NMD Resistant") %>%
                pull(ref_mes3ss),
              y = subsample %>%
                filter(NMD_type == "NMD Sensitive") %>%
                pull(ref_mes3ss),
              alternative = "greater",
              paired = T)
  rstatix::wilcox_effsize(data = subsample %>% 
                            mutate(NMD_type = NMD_type %>% as.factor()),
                          formula = ref_mes3ss ~ NMD_type,
                          paired = T)
  
  
  
  
  #################################################
  ## PLOTS
  #################################################
  
  plot_5ss <- ggplot(data = subsample,
                     aes(x=mean_phastCons17way5ss_100, y=ref_mes5ss, colour=fct_rev(NMD_type))) +
    geom_point() +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "bottom") +
    xlab("phastCons17 5'ss")  +
    ylab("MaxEntScan 5'ss") +
    ylim(c(-40, 20)) 
  plot_5ss 
  
  plot_3ss <- ggplot(data = subsample,
                     aes(x=mean_phastCons17way3ss_100, y=ref_mes3ss, colour=fct_rev(NMD_type))) +
    geom_point() +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "bottom") +
    xlab("phastCons17 3'ss")  +
    ylab("MaxEntScan 3'ss") +
    ylim(c(-40, 20))
  plot_3ss 
  
  
  
  ggpubr::ggarrange(ggExtra::ggMarginal(plot_5ss, groupColour = TRUE, groupFill = TRUE), 
                    ggExtra::ggMarginal(plot_3ss, groupColour = TRUE, groupFill = TRUE),
                    labels = c("A", "B"))
  
  
  
  ggsave(file = paste0(figures_folder, "/NMD_sensitive_resistant_introns_scatter_plot.png"), 
         width = 180, height = 90, dpi = 300, units = "mm")
  
  
  
}


plot_CLIP_ENCODE_by_RBP <- function() {
  

  
  ################################
  ## LOAD METADATA & RBPs
  ################################
  
  # target_genes <- "TARDBP"
  all_projects <- master_metadata$SRA_project %>% unique()
  target_genes <- "all"
  target_genes <- "TARDBP"
  all_projects <- "TARDBP"
  
  required_clusters <- c("case", "control")
  
  file_name <- paste0(results_folder, 
                      "/common_introns_details_", target_genes, "_experiments_", 
                      paste(required_clusters,collapse = "_"), ".rds")
  
  if( file.exists(file_name) ) {
    
    message("Loading common introns across all experiments...")
    common_introns_all_experiments <- readRDS(file = file_name)
    
  } else {
    
    ## Calculate common introns
    common_introns_all_experiments <- get_common_introns_across_experiments(RBPs = target_genes,
                                                                            required_clusters = c("case", "control"))
  }
  
  ## Same number of introns across RBPs
  common_introns_all_experiments %>%
    group_by(RBP, cluster) %>%
    distinct(ref_junID) %>%
    dplyr::count() %>%
    ungroup()
  
  
  ## Different number of novel junctions across RBPs
  common_introns_all_experiments %>%
    group_by(RBP, cluster) %>%
    distinct(novel_junID) %>%
    dplyr::count() %>%
    ungroup()
  
  
  ####################################
  ## LOAD MSR VALUES
  ####################################
  
  # Load the MSR tables for current RBP
  MSR_RBPs <- common_introns_all_experiments %>%
    inner_join(y = master_introns %>% dplyr::select(ref_junID, seqnames, start, end, strand),
               by = "ref_junID") %>%
    dplyr::select(ref_junID,  
                  MSR_D, 
                  MSR_A,
                  RBP, 
                  cluster,
                  seqnames, start, end, strand) %>%
    group_by(RBP, cluster) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup()
  
  
  MSR_RBPs %>% head()

  ####################################
  ## MERGE ENCODE and ICLIP DATA
  ## Use introns with increasing levels
  ## of MSR
  ####################################
  
  if ( !file.exists(file = paste0(results_folder, "/iCLIP_ENCODE_all_introns_chisq.csv")) ) {
    
    iCLIP_binding_MSR_results <- map_df(all_projects, function(target_RBP) {
      
      # target_RBP <- "TARDBP"
      # target_RBP <- "HNRNPC"
      # target_RBP <- "UPF1"
      # target_RBP <- "UPF1"
      
      
      ## Only CLIP data from 'HepG2' and 'K562' cell lines
      iclip_RBP <- map_df(c("HepG2","K562"), function(celltype) {
        # celltype <- "HepG2"
        
        iCLIP_RBP_file_path <- file.path(dependencies_folder, "ENCORI", paste0("ENCORI_hg38_RBPTarget_", target_RBP, "_",celltype,".txt"))
        
        if ( file.exists(iCLIP_RBP_file_path) ) {
          readr::read_delim(file = iCLIP_RBP_file_path,
                            delim = "\t", col_names = T, skip = 3, show_col_types = FALSE) %>%
            return()
        } else {
          return(NULL)
        }
      })
      
      ## Only if there's iCLIP data available for the current RBP in HepG2|K562 cell lines
      
      if ( iclip_RBP %>% nrow() > 1 ) {
        
        #message(Sys.time(), " - ", target_RBP)
        
        map_df(c("MSR_D", "MSR_A"), function(MSR_type) {
          
          message(Sys.time(), " - ", target_RBP, " ", MSR_type)
          
          # MSR_type <- "MSR_D"
          # MSR_type <- "MSR_A"
          
          MSR_RBP <- MSR_RBPs %>%
            filter(RBP == target_RBP) %>%
            pivot_wider(id_cols = ref_junID,
                        names_from = c("cluster"),
                        values_from = all_of(MSR_type)) %>%
            mutate(type = MSR_type) %>% 
            inner_join(y = MSR_RBPs %>%
                         filter(RBP == target_RBP) %>%
                         dplyr::select(ref_junID, seqnames, start, end, strand),
                       by = "ref_junID") %>%
            distinct(ref_junID, type, .keep_all = T) %>%
            mutate(start = start - 100,
                   end = end + 100)
          
          MSR_RBP_inc <- MSR_RBP %>%
            filter(case > control)  %>%
            distinct(ref_junID, .keep_all = T) %>%
            GRanges()
          MSR_RBP_inc
          
          MSR_RBP_notinc <- MSR_RBP %>%
            filter(case <= control) %>%
            distinct(ref_junID, .keep_all = T) %>%
            GRanges()
          MSR_RBP_notinc
          
          if ( intersect(MSR_RBP_inc$ref_junID, 
                         MSR_RBP_notinc$ref_junID) %>% length() > 0) {
            print("ERROR! Some introns with increasing MSR values are also classified as introns with non-increasing MSR values!")
            break;
          }
          
          ####################################
          ## TIDY ICLIP DATA 
          ####################################
          
          
          
          iclip_RBP_gr <- iclip_RBP %>%
            dplyr::select(seqnames = "chromosome", start = "broadStart", end = "broadEnd", "strand" ) %>%
            GRanges()
          
          ####################################
          ## FIND OVERLAPS - MSR and iCLIP
          ####################################
          
          
          ## 1. Introns with increasing MSR values ----------------------------------------------------------
          
          MSR_RBP_inc_overlaps <- GenomicRanges::findOverlaps(query = MSR_RBP_inc,
                                                              subject = iclip_RBP_gr,
                                                              type = "any")
          ## w binding motifs
          MSR_RBP_inc_overlapp_iCLIP <- MSR_RBP_inc[S4Vectors::queryHits(MSR_RBP_inc_overlaps) %>% unique,]
          MSR_RBP_inc_overlapp_iCLIP
          
          ## without binding motifs
          MSR_RBP_inc_NOT_overlapp_iCLIP <- MSR_RBP_inc[-c(S4Vectors::queryHits(MSR_RBP_inc_overlaps) %>% unique),]
          MSR_RBP_inc_NOT_overlapp_iCLIP
          
          
          
          
          ## 2. Introns with NOT increasing MSR values --------------------------------------------------------
          
          MSR_RBP_notinc_overlaps <- GenomicRanges::findOverlaps(query = MSR_RBP_notinc,
                                                                 subject = iclip_RBP_gr,
                                                                 type = "any")
          ## w binding motifs
          MSR_RBP_notinc_overlapp_iCLIP <- MSR_RBP_notinc[S4Vectors::queryHits(MSR_RBP_notinc_overlaps) %>% unique,]
          MSR_RBP_notinc_overlapp_iCLIP 
          
          ## without binding motifs
          MSR_RBP_notinc_NOT_overlapp_iCLIP <- MSR_RBP_notinc[-c(S4Vectors::queryHits(MSR_RBP_notinc_overlaps) %>% unique),]
          MSR_RBP_notinc_NOT_overlapp_iCLIP 
          
          
          ## 3. Save data --------------------------------------------------------------------------------------
          
          ENCODE_ENCORI_MSR_RBP_motif <- 
            rbind(
              
              cbind(
                
                MSR_RBP_notinc[S4Vectors::queryHits(MSR_RBP_notinc_overlaps),] %>% 
                  as_tibble() %>%
                  dplyr::rename(MSR_case = case,
                                MSR_control = control),
                
                iclip_RBP_gr[S4Vectors::subjectHits(MSR_RBP_notinc_overlaps),] %>% 
                  as_tibble() %>% 
                  dplyr::rename(iCLIP_seqnames = seqnames,
                                iCLIP_start = start,
                                iCLIP_end = end,
                                iCLIP_width = width,
                                iCLIP_strand = strand)
              ) %>% 
                mutate(RBP_motif = T) %>%
                
                plyr::rbind.fill(MSR_RBP_notinc[-c(S4Vectors::queryHits(MSR_RBP_notinc_overlaps)),] %>% 
                                   as_tibble()%>%
                                   dplyr::rename(MSR_case = case,
                                                 MSR_control = control) %>%
                                   mutate(RBP_motif = F)) %>%
                as_tibble() %>%
                mutate(MSR_direction = "not_increasing") , 
              
              cbind(
                
                MSR_RBP_inc[S4Vectors::queryHits(MSR_RBP_inc_overlaps),] %>% 
                  as_tibble() %>%
                  dplyr::rename(MSR_case = case,
                                MSR_control = control),
                
                iclip_RBP_gr[S4Vectors::subjectHits(MSR_RBP_inc_overlaps),] %>% 
                  as_tibble() %>% 
                  dplyr::rename(iCLIP_seqnames = seqnames,
                                iCLIP_start = start,
                                iCLIP_end = end,
                                iCLIP_width = width,
                                iCLIP_strand = strand)) %>% #
                mutate(RBP_motif = T) %>%
                
                
                plyr::rbind.fill(MSR_RBP_inc[-c(S4Vectors::queryHits(MSR_RBP_inc_overlaps)),] %>% 
                                   as_tibble()%>%
                                   dplyr::rename(MSR_case = case,
                                                 MSR_control = control) %>%
                                   mutate(RBP_motif = F)) %>%
                
                
                as_tibble() %>%
                mutate(MSR_direction = "increasing") 
              
            ) %>%
            dplyr::rename(MSR_type = type) %>%
            mutate(RBP = target_RBP)
          
          
          ENCODE_ENCORI_MSR_RBP_motif
          
          ENCODE_ENCORI_MSR_RBP_motif %>%
            distinct(ref_junID, .keep_all =T) %>%
            mutate(start = start + 100,
                   end = end - 100) %>%
            dplyr::count(MSR_direction, RBP_motif )
          
          

          write.csv(x = ENCODE_ENCORI_MSR_RBP_motif %>%
                      inner_join(y = master_introns %>% dplyr::select(ref_junID, transcript_id),
                                 by = "ref_junID") %>%
                      inner_join(y = master_transcript %>% dplyr::select(id, gene_id),
                                 by = c("transcript_id" = "id")) %>%
                      inner_join(y = master_gene %>% dplyr::select(id, gene_name),
                                 by = c("gene_id" = "id")) %>%
                      dplyr::select(-c("transcript_id", "gene_id")) %>%
                      dplyr::relocate(gene_name, .after = "strand") %>%
                      mutate(start = start + 100,
                             end = end - 100),
                    file = paste0(results_folder, paste0("/ENCODE_ENCORI_",target_RBP,"_",MSR_type,".csv")),
                    row.names = F)
          
          ###################################################################################################
          
          RBP_results <- data.frame(RBP= target_RBP,
                                    MSR_tested = MSR_type,
                                    total_introns = MSR_RBP %>% distinct(ref_junID) %>% nrow(),
                                    n_introns =  c(MSR_RBP_inc_overlapp_iCLIP$ref_junID %>% unique(), 
                                                   MSR_RBP_inc_NOT_overlapp_iCLIP$ref_junID  %>% unique(), 
                                                   MSR_RBP_notinc_overlapp_iCLIP$ref_junID  %>% unique(), 
                                                   MSR_RBP_notinc_NOT_overlapp_iCLIP$ref_junID  %>% unique()),
                                    RBP_iCLIP_type = c(rep("RBP_motif",times=length(MSR_RBP_inc_overlapp_iCLIP)),
                                                       rep("RBP_no_motif",times=length(MSR_RBP_inc_NOT_overlapp_iCLIP)),
                                                       rep("RBP_motif",times=length(MSR_RBP_notinc_overlapp_iCLIP)),
                                                       rep("RBP_no_motif",times=length(MSR_RBP_notinc_NOT_overlapp_iCLIP))), 
                                    MSR_direction = c(rep("MSR_increasing",times=length(MSR_RBP_inc_overlapp_iCLIP)+length(MSR_RBP_inc_NOT_overlapp_iCLIP)),
                                                      rep("MSR_not_increasing",times=length(MSR_RBP_notinc_overlapp_iCLIP)+length(MSR_RBP_notinc_NOT_overlapp_iCLIP))))
          
          
          chisq <- chisq.test( table(RBP_results$MSR_direction, RBP_results$RBP_iCLIP_type) )
          

          MSR_RBP_inc_overlapp_iCLIP 
          
          saveRDS(object = RBP_results,
                  file = file.path(results_folder, paste0("/ENCORI_iCLIP_",target_RBP, "_", MSR_type, "_ChiSquare.rds")))
          
          return(  data.frame(RBP= target_RBP,
                              MSR_tested = MSR_type,
                              MSR_increasing_w_RBP_motif = MSR_RBP_inc_overlapp_iCLIP$ref_junID %>% length(), 
                              MSR_not_increasing_w_RBP_motif = MSR_RBP_inc_NOT_overlapp_iCLIP$ref_junID %>% length(), 
                              MSR_increasing_N_RBP_motif = MSR_RBP_notinc_overlapp_iCLIP$ref_junID %>% length(),  
                              MSR_not_increasing_N_RBP_motif = MSR_RBP_notinc_NOT_overlapp_iCLIP$ref_junID %>% length(), 
                              chisq_statistic = chisq$statistic,
                              chisq_method = chisq$method,
                              chisq_pvalue = chisq$p.value)
          )
          
        })
        
      } else {
        message("No iCLIP HepG2|K562 data for ", target_RBP)
      }

    })
    
    
    iCLIP_binding_MSR_results
    
    write.csv(x = iCLIP_binding_MSR_results,
              file = paste0(results_folder, "/iCLIP_ENCODE_all_introns_chisq.csv"),
              row.names = F)
    
  } else {
    iCLIP_binding_MSR_results <- read.table(file = paste0(results_folder, "/iCLIP_ENCODE_all_introns_chisq.csv"),
                                            header = T, sep = ",")
  }
  
  
  ##############################
  ## PAPER STATS
  ##############################
  
  iCLIP_binding_MSR_results %>%
    filter(RBP != "SAFB2") %>%
    pull(RBP) %>% unique %>% length()
  iCLIP_binding_MSR_results %>% filter(MSR_tested == "MSR_D",
                                       RBP != "SAFB2") %>% pull(chisq_pvalue) %>% summary
  
  iCLIP_binding_MSR_results %>% filter(MSR_tested == "MSR_A",
                                       RBP != "SAFB2") %>% pull(chisq_pvalue) %>% summary
  
  ##############################
  ## CALCULATE FOLD-CHANGES
  ##############################

  iCLIP_binding_MSR_results_prop <- iCLIP_binding_MSR_results %>% 
    rowwise() %>%
    
    mutate(MSR_increasing_w_motif = MSR_increasing_w_RBP_motif / (MSR_increasing_w_RBP_motif + MSR_increasing_N_RBP_motif),
           MSR_increasing_n_motif = MSR_increasing_N_RBP_motif / (MSR_increasing_w_RBP_motif + MSR_increasing_N_RBP_motif),
           
           MSR_not_increasing_w_motif = MSR_not_increasing_w_RBP_motif / (MSR_not_increasing_w_RBP_motif + MSR_not_increasing_N_RBP_motif),
           MSR_not_increasing_n_motif = MSR_not_increasing_N_RBP_motif / (MSR_not_increasing_w_RBP_motif + MSR_not_increasing_N_RBP_motif)) %>%
    
    mutate(fold_change_w_motif = log2(MSR_increasing_w_motif) - log2(MSR_not_increasing_w_motif)) %>% 
    mutate(fold_change_not_motif = log2(MSR_increasing_n_motif) - log2(MSR_not_increasing_n_motif)) 
  
  
  iCLIP_binding_MSR_results_prop <- iCLIP_binding_MSR_results_prop  %>% 
    dplyr::select(RBP, 
                  MSR_tested, 
                  "Introns with RBP binding sites from CLIP-seq data" = fold_change_w_motif,
                  "Introns without RBP binding sites" = fold_change_not_motif) %>%
    gather(type,  "2fold", -MSR_tested, -RBP) %>% 
    arrange( MSR_tested,  desc("2fold"))  %>%
    ungroup
  
  
  ########################################################
  ## ADD KNOCKDOWN EFFICIENCY
  ########################################################
  
  project_path <- "ENCODE_SR/ENCODE_Metadata_Extraction/Metadata_results_experiments/"
  metadata_kEff_path <- paste0(project_path, "metadata_WB_kEff.tsv")
  metadata_kEff <- readr::read_delim(metadata_kEff_path, show_col_types = F) %>%
    mutate(kEff_text = ifelse(is.na(kEff_avg), kEff_avg, paste0(round(kEff_avg), "%"))) %>%
    mutate(kEff_text = kEff_text %>% as.factor())
  
  
  iCLIP_binding_MSR_results_prop <- iCLIP_binding_MSR_results_prop %>%
    left_join(y = metadata_kEff,
              by = c("RBP" = "target_gene")) %>%
    filter(RBP != "SAFB2")
  iCLIP_binding_MSR_results_prop %>% head()
  
  
  
  
  ########################################################
  ## BAR PLOT
  ########################################################
  
  iCLIP_binding_MSR_results_prop$MSR_tested = factor( iCLIP_binding_MSR_results_prop$MSR_tested, 
                                                      levels = c( "MSR_D", 
                                                                  "MSR_A"),
                                                      labels = c( "MSR Donor", 
                                                                  "MSR Acceptor"))
  
  
  iCLIP_plot <- ggplot(data = iCLIP_binding_MSR_results_prop,
                       aes(x = RBP, y = `2fold`, group = factor(MSR_tested)) ) +
    geom_bar(
      stat = "identity",
      aes(x = tidytext::reorder_within(x = RBP, by = -`2fold`, within = type), y = `2fold`,
          fill = factor(MSR_tested)), position = position_dodge(width = 0.9) ) +
    ylab("Log2 Fold Change of introns\nwith increasing MSR in KD vs wild-type") +
    xlab("") +
    theme_light() +
    ggforce::facet_row(~ type, scales = "free_x") +
    tidytext::scale_x_reordered() +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2,  nrow = 1 ))  +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"), #8d03b0
                      breaks = c("MSR Donor","MSR Acceptor")) +
    theme( plot.margin = margin(t = 5, r = 5, b = -5, l = 5),
           legend.box.margin = margin(l = -10, b = -10, t = -5 ))
  
  iCLIP_plot +
    ggnewscale::new_scale_fill() +
    geom_tile(stat = "identity", 
              aes(y = 1.4, 
                  x = tidytext::reorder_within(x = RBP,
                                               by = -`2fold`,
                                               within = type),
                  fill = kEff_avg, color = "No data\navailable"), 
              linewidth = 0.5, width = 1, height = 0.1) + 
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
    guides(colour = guide_legend(override.aes = list(fill = '#999999'),
                                 title = NULL,
                                 label.position = "bottom",
                                 order = 1))
  
  
  ggsave(file = paste0(figures_folder, "/iCLIP_ENCODE_2fold_efficiency.png"), 
         width = 180, height = 75, units = "mm", dpi = 300)
  
  
}


##################################
## CALLS
##################################

plot_effect_size_data_all_RBPs()
