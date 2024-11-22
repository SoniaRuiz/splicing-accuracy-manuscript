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

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/manuscript_figures.R")

base_folder <- here::here()

## CONNECT TO THE DATABASE ------------------------------
supporting_reads <- 1
gtf_version <- 105
main_project <- paste0("splicing_",supporting_reads,"read")

database_path <- paste0(base_folder, "/database/", main_project, "/", gtf_version, "/", main_project, ".sqlite")

con <- dbConnect(RSQLite::SQLite(), database_path)
tables <- dbListTables(con)


## SET PATHS TO FOLDERS

args <-
  list(
    dependencies_folder = file.path(here::here(), "dependencies"),
    results_folder = file.path(here::here(), "results", main_project, gtf_version, "_paper_review", "results"),
    figures_folder = file.path(here::here(), "results", main_project, gtf_version, "_paper_review", "figures"),
    data_folder = file.path(here::here(), "results", main_project, gtf_version, "_paper_review", "data")
  )



dir.create(file.path(args$results_folder), recursive = TRUE, showWarnings = F)
dir.create(file.path(args$figures_folder), recursive = TRUE, showWarnings = F)
dir.create(file.path(args$data_folder), recursive = TRUE, showWarnings = F)


## QUERY MASTER TABLES 

query = paste0("SELECT * FROM 'metadata'")
df_metadata <- dbGetQuery(con, query) %>% distinct(sample_id, .keep_all = T) %>% as_tibble()
all_projects <- df_metadata$SRA_project %>% unique
df_metadata$cluster %>% unique

# df_metadata %>% dplyr::count(cluster) %>% print(n=100)

query <- paste0("SELECT * FROM 'intron'")
master_introns <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'novel'")
master_novel_junctions <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'transcript'")
master_transcripts <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'gene'")
master_genes <- dbGetQuery(con, query) %>% as_tibble()


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



########################################
## FUNCTIONS - produce figures for the
## paper
########################################


## SECTION 1 - GENERAL TESTS --------------------------------------


## 1. Novel donor and acceptor junctions are commonly detected and exceed the number of unique annotated introns by an average of 11-fold

stats_results_section_1 <- function() {
  
  ## We found that 268,988 (82.8%) annotated introns had at least a single associated novel donor or acceptor junction,
  ## with only 55,968 annotated introns appearing to be precisely spliced across all the samples and tissues studied.

  master_introns %>% distinct(ref_junID) %>% nrow()
  master_introns %>% dplyr::count(misspliced)
  
  
  ## Collectively, we detected 3,865,268 unique novel junctions, equating to 14 novel junctions per annotated intron. 
  master_novel_junctions %>% nrow() 
  master_novel_junctions %>% dplyr::count(novel_type)
  
  (master_novel_junctions %>% nrow()) / (master_introns %>% filter(misspliced==1) %>% nrow())
  
  ## originating from 32,026 genes and 201,541 transcripts, and covering 13,949 different human samples and 40 human tissues
  
  master_genes %>% distinct(gene_id) %>% nrow()
  master_transcripts %>% nrow()
  
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
                                                  filter(misspliced == 1) %>%
                                                  pull(n))) %>%
    round()
  
  
}




## SECTION 2 - MAIN FIGURES ---------------------------------------


## 2. Over 98% of novel donor and acceptor junctions are likely to be generated through splicing errors

main_figure2 <- function () {

  
  results_file <- file.path(args$data_folder, "figure2_a.csv")
  
  if (file.exists(results_file)) {
    
    df_reclassification_rates_tidy <- read_csv(file = results_file, col_names = T, show_col_types = F)
    
  } else {
   
    reclassification_file <- file.path(args$results_folder, "/reclassification_rates.rds")
    
    if (!file.exists(reclassification_file)) {
      
      
      df_reclassification_rates <- map_df(all_projects, function(project_id) {
        
        # project_id <- all_projects[1]
        
        print(paste0(Sys.time(), " - ", project_id))
        
        all_clusters <- df_metadata %>% filter(SRA_project == project_id) %>% distinct(cluster) %>% pull()
        
        versions <- c("97")
        dates <- c("26-May-2019")
  
        map_df(all_clusters, function(cluster) {
          
          # cluster <- all_clusters[1]
          print(paste0(cluster))
          
          map_df(versions, function(version) {
            
            # version <- versions[1]
            print(paste0("v", version))
            
            if ( file.exists(file.path(base_folder, "results", main_project, version, project_id, "base_data",
                                       paste0(project_id, "_", cluster, "_all_split_reads.rds"))) ) {
  
              
              
              ## ENSEMBL v97
              df_ensembl97 <-  readRDS(file = file.path(base_folder, "results", main_project, version, project_id, "base_data",
                                                        paste0(project_id, "_", cluster, "_all_split_reads.rds")) ) %>% as_tibble()
              df_ensembl97_introns <- df_ensembl97 %>% filter(type == "annotated")
              df_ensembl97_novel <- df_ensembl97 %>% filter(type %in% c("novel_donor", "novel_acceptor"))
              
              
              ## ENSEMBL v105
              df_ensembl105 <- readRDS(file = file.path(base_folder, "results", main_project, "105", project_id, "base_data",
                                                 paste0(project_id, "_", cluster, "_all_split_reads.rds")) ) %>% as_tibble()
              df_ensembl105_introns <- df_ensembl105 %>% filter(type == "annotated")
              df_ensembl105_novel <- df_ensembl105 %>% filter(type %in% c("novel_donor", "novel_acceptor"))
              
              
              rm(df_ensembl97)
              rm(df_ensembl105)
              
              
              ## NOVEL JUNCTIONS THAT KEPT ANNOTATION CATEGORY
              kept_annotation <- df_ensembl97_novel %>%
                filter(junID %in% df_ensembl105_novel$junID) %>% 
                distinct(junID, .keep_all = T)
              kept_annotation %>% nrow() %>% print()
              
              ## NOVEL JUNCTIONS THAT ENTERED ANNOTATION
              in_annotation <- df_ensembl97_novel %>%
                filter(junID %in% df_ensembl105_introns$junID) %>% 
                distinct(junID, .keep_all = T)
              in_annotation %>% nrow() %>% print()
              
              label <- paste0("Ensembl_v", version, " (", dates, ")")
            
              # (in_annotation %>% nrow()) / (in_annotation %>% nrow() +
              #                                       keep_annotation %>% nrow()) * 100
              
              return(data.frame(tissue = cluster,
                                reclassification_rates = label,
                                kept_annotation = kept_annotation %>% nrow() / df_ensembl97_novel %>% nrow(),
                                in_annotation = in_annotation %>% nrow() / df_ensembl97_novel %>% nrow() ))
            } else {
              return(NULL)
            }
            
          })
        })
      })
      
      saveRDS(object = df_reclassification_rates, file = reclassification_file)
      
      
    } else {
      df_reclassification_rates <- readRDS(file = reclassification_file)
    }
    
    
    ## We found that across all tissues, on average only 0.008 [0.005,0.012] of junctions defined as novel donor or acceptor 
    ## junctions using v97 were reclassified as annotated introns in v105, and thus part of a transcript structure (Fig.2a)
   
    ## Stats - all tissues
    df_reclassification_rates %>% pull(in_annotation) %>% min
    df_reclassification_rates %>% pull(in_annotation) %>% max
    df_reclassification_rates %>% pull(in_annotation) %>% mean
    
    
    ## Interestingly, we noted that the highest re-classification rates were observed amongst human brain tissues, 
    ## on average 0.009 [0.008,0.012]. 
    
    ## Stats - only brain
    df_reclassification_rates %>% filter(str_detect(tissue, pattern = "Brain")) %>% pull(in_annotation) %>% mean
    df_reclassification_rates %>% filter(str_detect(tissue, pattern = "Brain")) %>% pull(in_annotation) %>% max
    df_reclassification_rates %>% filter(str_detect(tissue, pattern = "Brain")) %>% pull(in_annotation) %>% min
    
    
    
    ## Prepare the dataframe before plotting it
    df_reclassification_rates_tidy <- df_reclassification_rates %>% 
      dplyr::select(kept_annotation, in_annotation, tissue) %>%
      tidyr::gather(key = "type", value = "prop", -tissue ) 
    
    df_reclassification_rates_tidy$type = factor(df_reclassification_rates_tidy$type, levels = c( "kept_annotation", "in_annotation" ) )
    
    df_reclassification_rates_tidy = df_reclassification_rates_tidy %>% 
      ungroup() %>% arrange(type , prop) %>% mutate(tissue = fct_inorder(tissue))
  

    write.csv(x = df_reclassification_rates_tidy, 
              file = file.path(args$data_folder, "figure2_a.csv"), row.names = F)
    
  } 
  
  
  ###################################
  ## RECLASSIFICATION RATES BAR PLOT
  ##################################
  
  # 
  # 
  # df_reclassification_rates_tidy <- df_reclassification_rates_tidy %>% 
  #   filter(!(tissue %in% c("Brain - Cortex", "Brain - Cerebellum")))  %>%
  #   arrange(type, desc(prop)) %>%
  #   mutate(tissue = as.factor(x = tissue))
  # 
  # # Create factor column with decreasing order TRUE
  # df_reclassification_rates_tidy$tissue <- factor(df_reclassification_rates_tidy$tissue, 
  #                                                 levels = df_reclassification_rates_tidy$tissue[order(df_reclassification_rates_tidy$prop, decreasing = TRUE)])
  # 
  # 
  

  
  df_reclassification_rates_tidy <- df_reclassification_rates_tidy %>%
    arrange(type, desc(prop)) %>%  # First by 'type', then by 'prop' descending
    mutate(tissue = factor(tissue, levels = unique(tissue)))  # Reorder tissue
  df_reclassification_rates_tidy$type = factor(df_reclassification_rates_tidy$type, levels = c("kept_annotation", "in_annotation"))

  ggplot(data = df_reclassification_rates_tidy, 
         aes(x = tissue, y = prop, fill = type))  +
    geom_bar(stat = "identity", position = "stack") + 
    ggforce::facet_zoom(ylim = c(0,0.013), zoom.size = 3) +
    ylab("Ratio of unique novel junctions") +
    xlab("") +
    scale_fill_manual(values = c( "#a6a6a6", "#1a1a1a" ),
                      breaks = c( "kept_annotation", "in_annotation" ),
                      labels = c( "Novel junctions mantaining category", "Novel junctions entering annotation" )) +
    theme_light() +
    custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) %>%
    return()
  
  
  ggplot2::ggsave(file.path(args$figures_folder, "main_figure2.png"), width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(file.path(args$figures_folder, "main_figure2.svg"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
  
}

main_figure2_b <- function() {
  
  ##################################################
  ## Proportion of unique junctions per GTEx tissue
  ##################################################

  results_file <- file.path(args$data_folder, "figure2_b.csv")
  
  if ( file.exists(results_file) ) {
    
    df_proportions_violin <- read_csv(file = results_file, col_names = T)
    
  } else {
    
    file_name <- file.path(args$results_folder, "/unique_donor_acceptor_jxn.rds")
    
    if ( !file.exists(file_name) )  {
      
      df_proportions <- map_df(all_projects, function(project_id) {
        
        # project_id <- all_projects[1]
        
        print(paste0(Sys.time(), " - ", project_id))
        
        all_clusters <- df_metadata %>%
          filter(SRA_project == project_id) %>%
          distinct(cluster) %>%
          pull()
        
        map_df(all_clusters, function(cluster_id) {
          
          # cluster_id <- all_clusters[1]
          
          ## Print the tissue
          print(paste0(cluster_id))
          
          samples <- df_metadata %>%
            dplyr::count(cluster) %>%
            filter(cluster == cluster_id) %>% 
            pull(n)
          
          
          ####################
          ## GET THE INTRONS
          ####################
          
          query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
          introns <- dbGetQuery(con, query) %>% as_tibble()
          query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
          introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
          
          
          ###########################
          ## GET THE NOVEL JUNCTIONS
          ###########################
          
          query <- paste0("SELECT DISTINCT novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
          novel_junctions <- dbGetQuery(con, query) %>% as_tibble() %>%
            inner_join(y = master_novel_junctions,
                       by = "novel_junID")
          
          
          introns %>% head()
          novel_junctions %>% head()
          
          
          ######################################
          ## Calculate proportion of unique jxn
          ######################################
          
          annotated_junc <- introns %>% distinct(ref_junID) %>% nrow()
          donor_junc <- novel_junctions %>% filter(novel_type == "novel_donor") %>% distinct(novel_junID) %>% nrow()
          acceptor_junc <- novel_junctions %>% filter(novel_type == "novel_acceptor") %>% distinct(novel_junID) %>% nrow()
          
          annotated_prop <- annotated_junc/(annotated_junc + donor_junc + acceptor_junc)
          donor_prop <- donor_junc/(annotated_junc + donor_junc + acceptor_junc)
          acceptor_prop <- acceptor_junc/(annotated_junc + donor_junc + acceptor_junc)
          
          
          
          ## Return the data.frame
          return(data.frame(tissue = cluster_id,
                            annotated_junc = annotated_junc ,
                            donor_junc = donor_junc,
                            acceptor_junc = acceptor_junc,
                            annotated_prop = annotated_prop,
                            donor_prop = donor_prop,
                            acceptor_prop = acceptor_prop,
                            samples = samples))
        })
      })
      
      ## Save results --------------------------------------------------------------------------
      
      saveRDS(object = df_proportions, file = file_name)
      write.csv(x = df_proportions, file = paste0(args$results_folder, "unique_donor_acceptor_jxn.csv"), row.names = F)
    } else {
      
      df_proportions <- readRDS(file = file.path(args$results_folder, "unique_donor_acceptor_jxn.rds"))
      
    }
    
    
    
    ###################
    ## STATS
    ###################
    
    ## We found that novel donor and acceptor junctions consistently accounted for 
    ## the majority of unique junctions detected (X%, range), 
    
    (df_proportions %>% 
       pull(donor_prop) %>% 
       mean() + df_proportions %>% 
       pull(acceptor_prop) %>% 
       mean()) * 100
    
    
    df_proportions %>%
      mutate(novel_prop = (donor_prop + acceptor_prop) * 100) %>%
      arrange(desc(novel_prop))
    
    ## While we detected an average of 241,044 unique annotated junctions across all tissues, 
    ## unique novel donor and acceptor junctions averaged 249,366 and 360,400 respectively. 
    
    df_proportions %>% 
      pull(annotated_junc) %>% 
      mean
    df_proportions %>% 
      pull(donor_junc) %>% 
      mean
    df_proportions %>% 
      pull(acceptor_junc) %>% 
      mean
    
    
    ## Prepare data to plot
    
    df_proportions <- df_proportions %>%
      dplyr::select(tissue, 
                    donor = donor_prop, 
                    acceptor = acceptor_prop, 
                    annotated_intron = annotated_prop) %>%
      tidyr::gather(key = "type", value = "prop", -tissue ) 
    
    df_proportions_violin <- df_proportions %>%
      filter(type != "annotated_intron") %>%
      mutate(prop = prop * 100)
    
    df_proportions_violin$type = factor(df_proportions_violin$type, 
                                        levels = c("donor","acceptor"))
    
    ## Save source data 
    write_csv(x = df_proportions_violin, file = file.path(args$data_folder, "figure2_b.csv"), col_names = T)
    
  }
  
  
   
  ######################
  ## VIOLIN PLOT
  ######################
  
  ## PLOT
  junx_violin_plot <- ggplot(df_proportions_violin, aes(type, prop, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point(size = 0.5, colour = "#333333") +
    geom_line( aes(group = tissue), colour = "#333333",linewidth = .2 )  +
    theme_light() +
    ylab("% unique junctions") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_x_discrete( breaks = c( "acceptor", "donor"),# "annotated_intron"),
                        labels = c( "Novel Acceptor", "Novel Donor")) +
    scale_fill_manual(values = c( "#35B779FF", "#8d03b0"),
                      breaks = c( "donor", "acceptor"),
                      labels = c( "Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1)) +
    custom_ggtheme + 
    theme(legend.position = "right")
    
  
  junx_violin_plot
  
  
  ## CALL THE FUNCTION 'get_unique_donor_acceptor_reads' TO OBTAIN THE CUMULATIVE NUMBER OF READS
  ## NEEDED TO CREATE A GGARANGE PLOT WITH A AND B GRAPHS
  
  reads_violin_plot <- main_figure2_c()
  
  
  ggpubr::ggarrange(junx_violin_plot, reads_violin_plot, common.legend = T, labels = c("b", "c"), align = "h", ncol = 2, nrow = 1)
  
  file_name <- paste0(args$figures_folder, "/main_figure2bc")
  ggplot2::ggsave(filename = paste0(file_name,".png"), width = 180, height = 90, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name,".svg"), width = 180, height = 90, units = "mm", dpi = 300)
  
}

main_figure2_c <- function() {
  
  results_file <- file.path(args$data_folder, "figure2_c.csv")
  
  if ( file.exists(results_file) ) {
    
    df_mean_counts_violin <- read_csv(file = results_file, col_names = T)
    
    
  } else {
    
  
    file_name <- file.path(args$results_folder, "unique_donor_acceptor_reads.rds")
    
    if ( !file.exists(file_name) )  {
       
      df_mean_counts <- map_df(all_projects, function(project_id) {
        
        print(paste0(Sys.time(), " - ", project_id))
  
        ## GET THE CLUSTERS
        all_clusters <- df_metadata %>%
          filter(SRA_project == project_id) %>%
          distinct(cluster) %>%
          pull()
        
        
        map_df(all_clusters, function(cluster_id) {
          
          # project_id <- "BRAIN"
          # cluster_id <- "Brain - Frontal Cortex (BA9)"
          
          print(cluster_id)
          
          samples <- df_metadata %>%
            dplyr::count(cluster) %>%
            filter(cluster == cluster_id) %>% 
            pull(n)
          
          ####################
          ## GET THE INTRONS
          ####################
          
          query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals FROM '", 
                          cluster_id, "_", project_id, "_nevermisspliced'")
          introns <- dbGetQuery(con, query) %>% as_tibble()
          query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals FROM '", 
                          cluster_id, "_", project_id, "_misspliced'")
          introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
          
          
          
          ###########################
          ## GET THE NOVEL JUNCTIONS
          ###########################
          
          query <- paste0("SELECT novel_junID, novel_sum_counts, novel_n_individuals 
                          FROM '", cluster_id, "_", project_id, "_misspliced'")
          novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
         
          novel_junctions <- novel_junctions %>%
            left_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_type),
                      by = "novel_junID") %>% 
            as_tibble()
          
          
          
          ###########################
          ## GET THE PROPORTIONS
          ###########################
          
          annotated <- introns %>%
            dplyr::distinct(ref_junID, .keep_all = T) %>%
            pull(ref_sum_counts) %>% 
            sum()
          
          acceptor <- novel_junctions %>%
            filter(novel_type == "novel_acceptor") %>%
            dplyr::distinct(novel_junID, .keep_all = T) %>%
            pull(novel_sum_counts) %>% 
            sum()
          
          donor <- novel_junctions %>%
            filter(novel_type == "novel_donor") %>%
            dplyr::distinct(novel_junID, .keep_all = T) %>%
            pull(novel_sum_counts) %>% 
            sum()
          
          annotated_p = annotated * 100 / (annotated + acceptor + donor)
          acceptor_p = acceptor * 100 / (annotated + acceptor + donor)
          donor_p = donor * 100 / (annotated + acceptor + donor)
          
          return(data.frame(tissue = cluster_id,
                            samples = samples,
                            type = c("annotated","acceptor", "donor"),
                            prop_counts = c(annotated_p, acceptor_p, donor_p),
                            sum_counts = c(annotated, acceptor, donor)))
          
        })
        
      })
      
      # ## Save data
      saveRDS(object = df_mean_counts, file = file_name)
      
    } else {
      
      df_mean_counts <- readRDS(file = file_name) %>% as_tibble()
    }
    
    
    ###################
    ## SOME STATS
    ###################
    
    ## we found that novel donor and acceptor junctions together accounted 
    ## for only 0.3 – 0.4% of all junction reads whereas annotated junctions 
    ## accounted for 99.4% - X% of junction reads across all tissues evaluated 
    
    df_mean_counts %>%
      filter(type == "annotated") %>%
      pull(prop_counts ) %>% 
      sort()
    
    df_mean_counts %>%
      filter(type != "annotated") %>%
      group_by(tissue) %>%
      mutate(sum_prop_counts = prop_counts %>% sum()) %>% 
      ungroup() %>%
      distinct(tissue, .keep_all = T) %>%
      pull(sum_prop_counts) %>% 
      sort()
    
  
    df_mean_counts %>%
      filter(tissue == "Brain - Frontal Cortex (BA9)")
    
    ########################
    ## PREPARE DATA TO PLOT
    ########################
    
    df_mean_counts_violin <- df_mean_counts %>%
      filter(type != "annotated") 
    
    df_mean_counts_violin$type = factor(df_mean_counts_violin$type, 
                                        levels = c("donor", "acceptor"))
    
    
    ## Save source data 
    write_csv(x = df_mean_counts_violin, file = file.path(args$data_folder, "figure2_c.csv"), col_names = T)
  }
  
  
  ######################
  ## VIOLIN PLOT
  ######################
  
  ## Plot
  reads_violin_plot <- ggplot(df_mean_counts_violin, aes(type, prop_counts, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point(size = 0.5, colour = "#333333") +
    geom_line( aes(group = tissue), colour = "#333333",linewidth = .2  )  +
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
    
  
  reads_violin_plot %>%
    return()
  
  
}




## 4. High sequence similarity between novel splice sites and their annotated pairs explains splicing errors

main_figure3_ab <- function() {
  


  file_name <- file.path(args$results_folder, "df_mes.rds")
  
  df_mes <-  if ( file.exists(file_name) )  {
    
    print("Loading MES results...")
    readRDS(file = file_name) %>% distinct(ref_junID, novel_junID, .keep_all = T) %>% as_tibble()
    
  } else {
    
    ###########################
    ## GET THE MES
    ###########################
    
    ## Get the MES scores corresponding to the mis-spliced MASTER introns

    introns <- master_introns %>%
      filter(misspliced == 1) %>%
      dplyr::select(ref_junID, ref_mes5ss, ref_mes3ss)
    
    ## Get the MASTER novel junctions
    novel_junctions <- master_novel_junctions %>% 
      dplyr::select(ref_junID, novel_junID, novel_type, novel_mes5ss, novel_mes3ss) %>% 
      as_tibble() 
    
    
    ###########################
    ## MERGE THE DATA AND RETURN
    ###########################
    
    df_mes <- novel_junctions %>% 
      left_join(y = introns, by = "ref_junID")  %>%
      mutate(diff_ss5score = ref_mes5ss - novel_mes5ss,
             diff_ss3score = ref_mes3ss - novel_mes3ss)
    
    saveRDS(object = df_mes %>% distinct(ref_junID, novel_junID, .keep_all = T) %>% as_tibble(),
            file = file_name)
      
  }
  
  
  ##########################
  ## DATA MAIN FIGURE
  ##########################
  
  df_delta_5ss <- df_mes %>%
    filter(novel_type == "novel_donor") %>% 
    dplyr::select(diff_ss5score) %>%
    tidyr::gather(key = "type", value = "MES") 
  
  df_delta_3ss <- df_mes %>%
    filter(novel_type == "novel_acceptor") %>% 
    dplyr::select(diff_ss3score) %>%
    tidyr::gather(key = "type", value = "MES") 
    
  # Save source data 
  write_csv(x = df_delta_5ss, file = file.path(args$data_folder, "figure3_a.csv"), col_names = T)
  write_csv(x = df_delta_3ss, file = file.path(args$data_folder, "figure3_b.csv"), col_names = T)

    
  ##################################
  ## PLOT MAIN FIGURE 3 (a and b)
  ##################################

  generate_plot <- function(df, breaks_label, label_text) {
    colour <- if (breaks_label == "diff_ss5score") { "#35B779FF" } else { "#8d03b0"}
    ggplot(data = df) +
      geom_density(mapping = aes(x = MES, fill = type)) +
      geom_vline(xintercept = 0) +
      xlab("Delta MES") +
      xlim(c(-40, 65)) +
      ylim(c(0, 0.11)) +
      theme_light() +
      scale_color_viridis_d() +
      scale_fill_manual(values = colour, breaks = breaks_label, labels = label_text) +
      custom_ggtheme +
      guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  }
 
  ## Plot
  deltaplot5ss <- generate_plot(df = df_delta_5ss, breaks_label = "diff_ss5score", label_text = "Delta MES Donor")
  deltaplot3ss <- generate_plot(df = df_delta_3ss, breaks_label = "diff_ss3score", label_text = "Delta MES Acceptor")
  ggpubr::ggarrange(deltaplot5ss, deltaplot3ss, labels = c("a", "b"), ncol = 2, nrow = 1)
  

  figure_name = file.path(args$figures_folder, "main_figure3_ab")
  ggplot2::ggsave(filename = paste0(figure_name, ".png"), width = 180, height = 50, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(figure_name, ".svg"), width = 180, height = 50, units = "mm", dpi = 300)
  
  
  ##########################
  ## STATS FOR THE PAPER
  ##########################
  
  ## "[...] As would be expected, we found that the majority of novel 5’ and 3’ splice sites were weaker than the corresponding annotated site with 82.5% of novel 5’..." 
 
  ((df_mes %>% 
      filter(novel_type == "novel_donor") %>%
      dplyr::select(intron_MES = ref_mes5ss, novel_donor_MES = novel_mes5ss) %>%
      mutate(MES_diff = intron_MES - novel_donor_MES) %>%
      pull(MES_diff) %>% 
      sign() %>% 
      table() %>% 
      as.data.frame() %>%
      filter(. != -1) %>%
      pull(Freq) %>%
      sum) * 100) / ( master_novel_junctions %>% 
                        dplyr::count(novel_type) %>%
                        filter(novel_type == "novel_donor") %>%
                        pull(n) )
 
  ## "[...] and 85.2% of novel 3’ sites having positive delta MES scores [...]"
  
  ((df_mes %>% 
      filter(novel_type == "novel_acceptor") %>%
      dplyr::select(intron_MES = ref_mes5ss, novel_acceptor_MES = novel_mes3ss) %>%
      mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
      pull(MES_diff) %>% 
      sign() %>% 
      table() %>% as.data.frame() %>%
      filter(. != -1) %>%
      pull(Freq) %>%
      sum) * 100) / ( master_novel_junctions %>% 
                        dplyr::count(novel_type) %>%
                        filter(novel_type == "novel_acceptor") %>%
                        pull(n) )
   
  
  ## "[...] Interestingly, this analysis demonstrated that novel 5’ splice sites had a modal delta value which was very close to zero "
  
  df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron_MES = ref_mes5ss, novel_donor_MES = novel_mes5ss) %>%
    mutate(MES_diff = intron_MES - novel_donor_MES) %>%
    pull(MES_diff) %>%
    get_mode()
  
  df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron_MES = ref_mes5ss, novel_donor_MES = novel_mes5ss) %>%
    mutate(MES_diff = intron_MES - novel_donor_MES) %>%
    pull(MES_diff) %>%
    summary()

  
  ## "[...] The modal delta value for novel 3’ splice sites was higher "
  
  df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron_MES = ref_mes3ss, novel_acceptor_MES = novel_mes3ss) %>%
    mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
    pull(MES_diff) %>%
    get_mode()
  
  
  df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron_MES = ref_mes3ss, novel_acceptor_MES = novel_mes3ss) %>%
    mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
    pull(MES_diff) %>%
    summary()
  
  
  df_delta_3ss$MES %>% min
  df_delta_3ss$MES %>% max
  
}

main_figure3_cd <- function() {
  
  limit_bp <- 30
  subsample = F
  
  results_file1 <- file.path(args$data_folder, "figure3_c.csv")
  results_file2 <- file.path(args$data_folder, "figure3_d.csv")
  
  
  if (file.exists(results_file1) && file.exists(results_file2)) {
    
    df_figure3_c <- read_csv(file = results_file1, col_names = T, show_col_types = F)
    df_figure3_d <- read_csv(file = results_file2, col_names = T, show_col_types = F)
    
  } else {
    
    ###############################
    ## GET DATA FOR FRONTAL CORTEX
    ###############################
    
    project_id <- "BRAIN"
    cluster_id <- "Brain - Frontal Cortex (BA9)"
    
    query <- paste0("SELECT tissue.novel_junID, tissue.ref_junID, tissue.ref_sum_counts, tissue.ref_n_individuals, novel.novel_type, novel.distance 
                    FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue 
                    INNER JOIN 'novel' ON novel.novel_junID = tissue.novel_junID")
    db_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
    
    
    
    ###############################
    ## PREPARE DATA BEFORE PLOT
    ###############################
  
    
    df_novel_tidy <- db_misspliced_introns %>%
      inner_join(master_introns %>% dplyr::select(ref_junID, protein_coding), by = "ref_junID") %>%
      filter(protein_coding %in% c(0,100)) %>%
      mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC"),
             novel_type = str_replace(string = novel_type, pattern = "_", replacement = " "),
             novel_type = str_to_title(novel_type)) 
    
    df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, levels = c("Novel Donor", "Novel Acceptor"))
    df_novel_tidy$type_PC = factor(df_novel_tidy$type_PC, levels = c("protein coding (PC)", "non PC"))
  
      
    if (subsample) {
      
      set.seed(100)
      
      df_novel_tidy_w_expression <- df_novel_tidy %>%
        group_by(type_PC) %>%
        distinct(ref_junID, .keep_all = T) %>%
        mutate(mean_expression = log10(ref_sum_counts/ref_n_individuals)) %>%
        ungroup()
      
      ## Subsampling introns to control by similarity in mean read coverage
      m.out <- MatchIt::matchit(type_PC ~ mean_expression, 
                                data = df_novel_tidy_w_expression %>% dplyr::select(-distance), 
                                distance = df_novel_tidy_w_expression$mean_expression,
                                method = "nearest", 
                                caliper = c(mean_expression = 0.005), 
                                std.caliper = FALSE)
      
      data_subsample <- MatchIt::match.data(m.out)

      
      ## Only consider matched annotated introns (i.e. introns with similar expression levels in protein-coding and non-coding transcripts)
      df_novel_tidy <- df_novel_tidy %>% dplyr::select(ref_junID, distance) %>%
        inner_join(y = data_subsample %>% dplyr::select(-c(distance, weights)), by = "ref_junID") 
      
    }
    
    df_novel_tidy %>% distinct(ref_junID, .keep_all = T) %>% filter(abs(distance) <= limit_bp) %>% dplyr::count(type_PC)
    df_novel_tidy %>% distinct(novel_junID, .keep_all = T) %>% filter(abs(distance) <= limit_bp) %>% dplyr::count(type_PC)
    
    df_figure3_c = rbind(df_novel_tidy %>% filter(type_PC == "protein coding (PC)", novel_type == "Novel Donor"),
                         df_novel_tidy %>% filter(type_PC == "protein coding (PC)", novel_type == "Novel Acceptor"))
    df_figure3_d = rbind(df_novel_tidy %>% filter(type_PC == "non PC", novel_type == "Novel Donor"),
                         df_novel_tidy %>% filter(type_PC == "non PC", novel_type == "Novel Acceptor"))
    
    # Save source data 
    write_csv(x = df_figure3_c, file = file.path(args$data_folder, "figure3_c.csv"), col_names = T)
    
    # Save source data 
    write_csv(x = df_figure3_d, file = file.path(args$data_folder, "figure3_d.csv"), col_names = T)
  
  } 
  

 
  #####################################
  ## FIGURES EXON AND INTRON DRAWS
  #####################################
  
  distance_rectangle_donor <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = 15, y = 80), size = 2.5, label = "intron") +
    geom_rect(aes(xmin = (limit_bp * -1), xmax = 0, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = -15, y = 55), size = 2.5, label = "exon") +
    theme_void()
  
  
  distance_rectangle_acceptor <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 2.5, label = "exon") +
    geom_rect(aes(xmin = (limit_bp * -1), xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15, y = 80),  size = 2.5, label = "intron") +
    theme_void()
  
  
  #####################################
  ## PROTEIN CODING -  DONOR & ACCEPTOR
  #####################################
  
  ## Donor
  plot_PC_donor <- ggplot(data = df_figure3_c %>% filter(novel_type == "Novel Donor")) + 
    geom_histogram(aes(x = distance, fill = novel_type), bins = limit_bp * 2, binwidth = 1, position = "stack") +
    ggtitle("Protein-coding transcripts") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c(limit_bp,(limit_bp * -1)), 
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6), trans = "reverse") +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  plot_PC_donor <- plot_PC_donor / distance_rectangle_donor +  patchwork::plot_layout(heights = c(8, 2))
  
  
  ## Acceptor
  plot_PC_acceptor <- ggplot(data = df_figure3_c %>% filter(novel_type == "Novel Acceptor")) + 
    geom_histogram(aes(x = distance, fill = novel_type), bins = limit_bp * 2, binwidth = 1, position = "stack") +
    ggtitle(" ") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp), breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  #plot_PC_acceptor <- plot_PC_acceptor + ggpubr::rremove("y.axis") + ggpubr::rremove("ylab") + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text") 
  plot_PC_acceptor <- plot_PC_acceptor / distance_rectangle_acceptor +  patchwork::plot_layout(heights = c(8, 2))
  
  ## Combine protein-coding donor and acceptor
  plot_PC <- ggpubr::ggarrange(plot_PC_donor, plot_PC_acceptor, ncol = 2, nrow = 1)
  
  
  #################################
  ## NON-CODING - DONOR & ACCEPTOR
  #################################
  
  ## Donor
  plot_NPC_donor <- ggplot(data = df_figure3_d %>% filter(novel_type == "Novel Donor")) + 
    geom_histogram(aes(x = distance, fill = novel_type), bins = limit_bp * 2, binwidth = 1, position = "stack") +
    ggtitle("Non-coding transcripts") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c(limit_bp,(limit_bp * -1)), breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6), trans = "reverse") +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  plot_NPC_donor <- plot_NPC_donor / distance_rectangle_donor +  patchwork::plot_layout(heights = c(8, 2))
  
  
  ## Acceptor
  plot_NPC_acceptor <- ggplot(data = df_figure3_d %>% filter(novel_type == "Novel Acceptor")) + 
    geom_histogram(aes(x = distance, fill = novel_type), bins = limit_bp * 2, binwidth = 1, position = "stack") +
    ggtitle(" ") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp), breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  #plot_NPC_acceptor <- plot_NPC_acceptor + ggpubr::rremove("y.axis") + ggpubr::rremove("ylab") + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text") 
  plot_NPC_acceptor <- plot_NPC_acceptor / distance_rectangle_acceptor +  patchwork::plot_layout(heights = c(8, 2))
  
  ## Combine non-coding donor and acceptor
  plot_NPC <- ggpubr::ggarrange(plot_NPC_donor, plot_NPC_acceptor,ncol = 2,nrow = 1)
  
  
  #####################################
  ## FIGURE LEGEND
  #####################################

  
  
  plot_legend <- ggpubr::get_legend(ggplot(data = df_figure3_c) + 
                                      geom_histogram(aes(x = distance, fill = novel_type))+
                                      scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("Novel Donor", "Novel Acceptor")) + 
                                      custom_ggtheme +
                                      theme(legend.spacing.x = unit(0.5, "lines")))
  
  
  #################################
  ## COMBINE ALL PLOTS 
  #################################
  
  combined_plot <- ggpubr::ggarrange(x = plot_PC, y = plot_NPC, 
                                     legend.grob = plot_legend, labels = c("c", "d"),
                                     ncol = 1, nrow = 2, 
                                     font.label = list(size = 14, color = "black", face = "bold", family = NULL))
  
  combined_plot 
  
  figure_name <- paste0(args$figures_folder, "/main_figure3cd")
  ggplot2::ggsave(filename = paste0(figure_name, ".png"), plot = combined_plot, width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(figure_name, ".svg"), plot = combined_plot, width = 180, height = 100, units = "mm", dpi = 300)

}

main_figure3_e <- function() {
  
  ############################
  ## LOAD MODULO 3 DATA 
  ############################

  file_name <- file.path(args$results_folder, "figure3_e.rds")
  
  if (!file.exists(file_name))  {
    
    df_modulo_tissues <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      print(paste0(Sys.time(), " - ", project_id))
      
      all_clusters <- df_metadata %>% filter(SRA_project == project_id) %>% distinct(cluster) %>% pull()
      
      map_df(all_clusters, function(cluster_id) {
        
        # cluster_id <- all_clusters[1]
        
        ## Print the tissue
        print(paste0(Sys.time(), " - ", cluster_id))
        
        #####################################
        ## Get the novel junctions from the current tissue
        #####################################
        
        query <- paste0("SELECT novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
        introns <- dbGetQuery(con, query) %>% as_tibble()
        
        ## Get the distance location
        query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                        paste(introns$novel_junID, collapse = ","),")")
        introns <- introns %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "novel_junID") %>% 
          as_tibble() 
        
        ## Add the transcript and MANE info
        query <- paste0("SELECT intron.ref_junID, intron.protein_coding, transcript.MANE
                        FROM 'intron' 
                        INNER JOIN transcript
                        ON intron.transcript_id = transcript.id
                        WHERE ref_junID IN (", paste(introns$ref_junID, collapse = ","),") 
                        AND transcript.MANE = 1")
        introns <- introns %>%
          inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "ref_junID") %>% 
          as_tibble() 
        
        df_novel_tidy <- introns %>%
          distinct(novel_junID, .keep_all = T) %>%
          filter(abs(distance) <= 100, MANE == 1) %>% 
          mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " ")) %>%
          mutate(type_p = ifelse(distance < 0, paste0(novel_type," intron"), paste0(novel_type," exon"))) %>% 
          mutate(modulo = abs(distance) %% 3)
        
        df_novel_tidy <- df_novel_tidy %>% 
          group_by(modulo) %>%
          summarise(n = n()) %>%
          mutate(freq = n / sum(n)) %>%
          mutate(tissue = cluster_id)
        
        return(df_novel_tidy)
      })
    })
    
    saveRDS(object = df_modulo_tissues, file = file_name)
    write_csv(x = df_modulo_tissues, file = file.path(args$data_folder, "figure3_e.csv"), col_names = T)
    
  } else {
    df_modulo_tissues <- readRDS(file = file_name) %>% as_tibble()
  }
  
  
  ###################
  ## DENSITY PLOT
  ###################
  
  df_modulo_tissues$modulo = factor(df_modulo_tissues$modulo, levels = c( "0", "1", "2"))
  
  ggplot(df_modulo_tissues, aes(x = freq, y = modulo)) +
    ggridges::geom_density_ridges_gradient() +
    ylab("Modulo3 of the distance") +
    xlab("% of novel junctions") +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(30,35,40), labels = c("30%","35%","40%")) +
    scale_y_discrete(expand = c(0,0.5,1,0))
  
  
  figure_name <- file.path(args$figures_folder, "main_figure3_e")
  ggplot2::ggsave(filename = paste0(figure_name, ".png"), width = 180, height = 50, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(figure_name, ".svg"), width = 180, height = 50, units = "mm", dpi = 300)
  
  
  
  ##############################
  ## STATS
  ##############################
  
  ## FCTX
  df_modulo_tissues %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)") %>%
    group_by(modulo) %>%
    mutate(mean_freq = freq %>% mean) %>%
    ungroup()
  
  df_modulo_tissues %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)") %>%
    group_by(modulo) %>%
    mutate(mean_freq = freq %>% mean) %>%
    ungroup() %>%
    mutate(modulo_type = ifelse(modulo == 0, "maintain_frame", "alter_frame")) %>%
    group_by(modulo_type) %>%
    distinct(mean_freq) %>%
    mutate(sum_freq = mean_freq %>% sum) 
  
  
  ## ALL TISSUES
  df_modulo_tissues %>%
    group_by(modulo) %>%
    mutate(mean_freq = freq %>% mean) %>%
    ungroup() %>%
    mutate(modulo_type = ifelse(modulo == 0, "maintain_frame", "alter_frame")) %>%
    group_by(modulo_type) %>%
    distinct(mean_freq) %>%
    mutate(sum_freq = mean_freq %>% sum) %>%
    pull(sum_freq)
  
  
}




## 6. Splicing error rates vary across introns and are likely to be underestimated in bulk RNA-seq data 

main_figure4_a <- function()  {
  
  
  results_file1 <- file.path(args$data_folder, "figure4_a.csv")
  
  if (file.exists(results_file1)) {
    
    df_all_introns_tidy <- read_csv(file = results_file1, col_names = T, show_col_types = F)
    
  } else {
    
    ## GET DATA FOR FRONTAL CORTEX
    con <- dbConnect(RSQLite::SQLite(), database_path)
    project_id <- "BRAIN"
    cluster_id <- "Brain - Frontal Cortex (BA9)"
    print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))
    
    query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", project_id, "_nevermisspliced' ")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", project_id, "_misspliced' ")
    introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
    
    ## Add information from the master table
    introns <- introns %>% left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding) %>% as_tibble(), by = "ref_junID") %>% as_tibble() 
    message(Sys.time(), " - ", introns %>% nrow(), " - introns!")
    
    ## 1. Compare mis-spliced vs never mis-spliced introns
    df_all_introns_tidy <- introns %>%
      mutate(mean_expression = ref_sum_counts/ref_n_individuals) %>%
      dplyr::select(ref_junID, MSR_D, MSR_A, ref_type, mean_expression,protein_coding) %>%
      gather(key = "MSR_type", value = "MSR", -ref_junID, -ref_type, -mean_expression,-protein_coding) %>%
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
    
    ## Save source data 
    write_csv(x = df_all_introns_tidy, file = file.path(args$data_folder, "figure4_a.csv"), col_names = T)
  }
  
  ###########################
  ## WILCOXON TEST
  ###########################
  
  
  df_all_introns_tidy_test <- df_all_introns_tidy %>%
    filter(ref_junID %in% intersect(df_all_introns_tidy %>% filter(MSR_type == "MSR_D") %>% pull(ref_junID),
                                    df_all_introns_tidy %>% filter(MSR_type == "MSR_A") %>% pull(ref_junID)) %>% unique()) %>%
    dplyr::select(-c(ref_type, type_label, percentile_group,mean_expression,protein_coding))%>%
    spread(MSR_type, value = MSR)
  
  
  # Focusing on frontal cortex brain tissue, we observed that while splicing errors were detected infrequently, 
  # with the MSRD and MSRA values highly skewed towards low values, there was considerable variation across introns 
  # (MSRD IQR=7.2e-04; MSRA IQR=1.9e-03). 
  IQR(x = df_all_introns_tidy_test$MSR_D)
  IQR(x = df_all_introns_tidy_test$MSR_A)
  
  
  ## Furthermore, consistent with the overall higher detection of novel acceptors as compared to novel donor junctions, 
  ## we observed a significant difference in the median of the two MSRD and MSRA distributions 
  ## (paired one-tailed Wilcoxon Rank-sum test, effect-size=0.09, P<0.001) 
  
  wilcox.test(x = df_all_introns_tidy_test$MSR_D,
              y = df_all_introns_tidy_test$MSR_A,
              alternative = "less")
  
  rstatix::wilcox_effsize(data = df_all_introns_tidy_test %>%
                            gather(key=MSR_type, value = MSR, -ref_junID) %>%
                            mutate(MSR_type = MSR_type %>% as.factor()),
                          formula = MSR ~ MSR_type,
                          alternative = "less")
  
  ###########################
  ## PLOT MSR
  ###########################
  
  plot1 <- ggplot(data = df_all_introns_tidy) + 
    geom_bar(aes(x = type_label, fill = MSR_type), position = "dodge") +
    ggplot2::labs(x = "", y = "Number of annotated introns")+
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("MSR_D","MSR_A"), labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  
  ## 2. Mis-splicing is more common at the acceptor
  df_all_introns_tidy$percentile_group = factor(df_all_introns_tidy$percentile_group, levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  plot2 <- ggplot(data = df_all_introns_tidy %>%
                    filter(percentile_group != "0")) + 
    geom_bar(aes(x = percentile_group, fill = MSR_type),position = "dodge")+
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    ggforce::facet_zoom(ylim=c(0,2500), split = TRUE) +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"), breaks = c("MSR_D","MSR_A"), labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  
  df_all_introns_tidy$ref_junID %>% unique %>% length()
  ggpubr::ggarrange(plot1, plot2, labels = c("", ""), nrow = 2, ncol = 1, heights = c(1,1.5), common.legend = T)
  
  figure_name <- file.path(args$figures_folder, "main_figure4_a")
  ggplot2::ggsave(filename = paste0(figure_name, ".png"), width = 180, height = 75, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(figure_name, ".svg"), width = 180, height = 75, units = "mm", dpi = 300)

}

main_figure4_bc <- function() {
  
  # Check if subsampling file exists
  subsample_file <- paste0(args$results_folder, "/supplementary_figure10.rds")
  
  subsample <- if (!file.exists(subsample_file)) {
    # Subsample using MatchIt
    supplementary_figure10()
  } else {
    readRDS(subsample_file)
  }
  

  subsample_MSR_tidy <- subsample %>%
    dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID) 
  subsample_MSR_tidy <- subsample %>% 
    dplyr::select(biotype, subclass, MSR_D, MSR_A, ref_junID) %>%
    arrange(subclass)
  
  
  # Get MSR Donor and Acceptor data
  subsample_MSR_D <- subsample_MSR_tidy %>%
    dplyr::select(biotype, subclass, MSR = MSR_D, ref_junID) %>%
    mutate(percentile_group = cut(MSR, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), include.lowest = F)) %>%
    mutate(percentile_group = percentile_group %>% as.character(),
           percentile_group = replace_na(percentile_group, "0"),
           percentile_group = factor(percentile_group, levels = c("0", "(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]", "(0.8,1]")))
  subsample_MSR_A <- subsample_MSR_tidy %>%
    dplyr::select(biotype, subclass, MSR = MSR_A, ref_junID) %>%
    mutate(percentile_group = cut(MSR, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), include.lowest = F)) %>%
    mutate(percentile_group = percentile_group %>% as.character(),
           percentile_group = replace_na(percentile_group, "0"),
           percentile_group = factor(percentile_group, levels = c("0", "(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]", "(0.8,1]")))
  
  
  # Plot MSR Donor and Acceptor
  plot_MSR <- function(df, title, zoomed = F) {
    plot <- ggplot(df %>% distinct(ref_junID, .keep_all = T)) +
      geom_bar(aes(x = percentile_group, fill = biotype), position = "dodge", linewidth = 0.5, color = "#333333") +
      ggtitle(title) + xlab("Mis-splicing ratio value group") + ylab("Number of annotated introns") +
      scale_fill_manual(values = c("#333333", "#999999"), 
                        breaks = c("PC", "non PC"),
                        labels = c("Protein-coding", "Non-coding")) +
      theme_light() + custom_ggtheme + guides(fill = guide_legend(ncol = 2, nrow = 1)) 
    if (zoomed) { plot + scale_y_continuous(limits =c(0,750), position = "right")} else {plot}
  }
  
  
  plot_MSR_donor <- plot_MSR(df = subsample_MSR_D, title = "MSR Donor")
  plot_MSR_acceptor <- plot_MSR(df = subsample_MSR_A, title = "MSR Acceptor")
  final_plot <- ggpubr::ggarrange(plot_MSR_donor, plot_MSR_acceptor, nrow = 1, ncol = 2, common.legend = TRUE)
  figure_name <- file.path(args$figures_folder, "/main_figure4_bc")
  ggsave(filename = paste0(figure_name, ".png"), final_plot, width = 180, height = 60, units = "mm", dpi = 300)
  ggsave(filename = paste0(figure_name, ".svg"), final_plot, width = 180, height = 60, units = "mm", dpi = 300)
  
  ## ZOMMED version for completion of figure4 bc 
  plot_MSR_donor_zoomed <- plot_MSR(df = subsample_MSR_D, title = "MSR Donor", zoomed = T)
  plot_MSR_acceptor_zoomed <- plot_MSR(df = subsample_MSR_A, title = "MSR Acceptor", zoomed = T)
  final_plot_zoomed <- ggpubr::ggarrange(plot_MSR_donor_zoomed, plot_MSR_acceptor_zoomed, nrow = 1, ncol = 2, common.legend = TRUE)
  figure_name_zoomed <- file.path(args$figures_folder, "/main_figure4_bc_zoomed")
  ggsave(filename = paste0(figure_name_zoomed, ".png"), final_plot_zoomed, width = 180, height = 60, units = "mm", dpi = 300)
  ggsave(filename = paste0(figure_name_zoomed, ".svg"), final_plot_zoomed, width = 180, height = 60, units = "mm", dpi = 300)
  
  
  #################################
  ## TEST
  #################################
  
  ## We found that splicing inaccuracies were more frequent amongst annotated introns from non-protein-coding transcripts 
  ## at both the 5’ss (paired one-tailed Wilcoxon Rank-sum test, effect-size=0.17, P<0.001) and 3’ss 
  ## (paired one-tailed Wilcoxon Rank-sum test, effect-size=0.19, P<0.001) 
  
  wilcox.test(x = subsample %>% filter(biotype =="PC") %>% pull(MSR_D),
              y = subsample %>% filter(biotype =="non PC") %>% pull(MSR_D),
              alternative = "less")
  rstatix::wilcox_effsize(data = subsample_stats %>% mutate(biotype = biotype %>% as.factor()),
                          formula = MSR ~ biotype,
                          alternative = "less")
  
  
  
  
  
  wilcox.test(x = subsample_stats %>% filter(biotype =="PC", MSR_type == "MSR_A") %>% pull(MSR),
              y = subsample_stats %>% filter(biotype =="non PC", MSR_type == "MSR_A") %>% pull(MSR),
              alternative = "less")
  rstatix::wilcox_effsize(data = subsample_stats %>% filter(MSR_type == "MSR_A") %>% mutate(biotype = biotype %>% as.factor()),
                          formula = MSR ~ biotype)
  
}

get_data_main_figure4_d <- function(common_introns = NULL, replace = F) {
  
  intron_size <- 100
  phastcons_type <- 17
  
  if (is.null(common_introns)) {
    common_introns <- F
  }
  
  folder_path <- file.path(args$results_folder, "ZIP/")
  dir.create(path = folder_path,showWarnings = F,recursive = T)

  ###############################
  ## QUERY THE DATABASE
  ## GET DATA FOR FRONTAL CORTEX (or all tissues if preferred)
  ###############################
  
  for (project_id in all_projects) {
    
    # project_id <- all_projects[1]
    ## GET THE CLUSTERS
    all_clusters <- df_metadata %>% filter(SRA_project == project_id) %>% distinct(cluster) %>% pull()
    
    for (cluster_id in all_clusters) { 
      
      # cluster_id <- all_clusters[1]
    
      if (!file.exists(paste0(folder_path, cluster_id, "_common_introns", common_introns, "_introns_phastcons", phastcons_type, "_intronsize", intron_size ,".rds")) || 
          replace) {
        
        message(project_id, " -> ", cluster_id, "...")
        
        ## GET MIS-SPLICED INTRONS AND JOIN WITH NOVEL MASTER DATA
        query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced'")
        introns_ms <- dbGetQuery(con, query) %>% as_tibble()
    
        ## JOIN WITH MASTER NOVEL JUNCTIONS
        introns_ms <- introns_ms %>% inner_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_type), by = c("novel_junID")) %>% as_tibble() 
        
        ## GET NEVER MIS-SPLICED INTRONS AND JOIN WITH MIS-SPLICED DATA
        query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced' " )
        introns <- dbGetQuery(con, query) %>% as_tibble()
        
        ## JOIN WITH MISSPLICED INTRONS
        introns <- plyr::rbind.fill(introns, introns_ms ) %>% as_tibble() 
        
        ###################################
        ## JOIN WITH MASTER INTRON TABLE
        ###################################
  
        introns <- introns %>% inner_join(y = master_introns %>% dplyr::select(-transcript_id), by = "ref_junID") %>% as_tibble()
        
        ## Discard introns with double the size of the sequence used for conservation and CDTS calculation
        introns %>% filter(ref_length < (intron_size * 2)) %>% distinct(ref_junID)
        introns <- introns %>% filter(ref_length >= (intron_size * 2)) 
        summary(introns)
        introns %>% distinct(ref_junID)
        
        #####################################
        ## JOIN WITH MASTER TRANSCRIPT TABLE
        #####################################
        
        query <- paste0("SELECT * FROM 'transcript' WHERE id IN (", paste(introns$transcript_id, collapse = ","),")")
        introns <- introns %>% inner_join(y = dbGetQuery(con, query), by = c("transcript_id" = "id")) %>%  as_tibble() 
        
        introns %>% nrow()
        introns %>% head()
        
        #####################################
        ## JOIN WITH MASTER GENE TABLE
        #####################################
        
        query <- paste0("SELECT * FROM 'gene' WHERE id IN (", base::paste(introns$gene_id, collapse = ","),")")
        introns <- introns %>% inner_join(y = dbGetQuery(con, query) %>% as_tibble(), by = c("gene_id" = "id")) %>% as_tibble() 
        
        introns %>% nrow()
        introns %>% head()
        
        introns$ref_junID %>% unique %>% length()
        
        
        #########################################################
        ## PREPARE THE DATA PRIOR MODELLING ---------------------
        #########################################################
        
        message("PhastCons type: ", phastcons_type)

        idb <- introns %>%
          distinct(ref_junID, novel_junID, .keep_all = T)  %>%
          dplyr::select(ref_junID,
                        intron_length = ref_length,
                        intron_5ss_score = ref_mes5ss,
                        intron_3ss_score = ref_mes3ss,
                        gene_num_transcripts = n_transcripts,
                        mean_CDTS5ss = paste0("mean_CDTS5ss_", intron_size),
                        mean_CDTS3ss = paste0("mean_CDTS3ss_", intron_size),
                        protein_coding,
                        mean_phastConsway5ss = paste0("mean_phastCons",phastcons_type,"way5ss_",intron_size),
                        mean_phastConsway3ss = paste0("mean_phastCons",phastcons_type,"way3ss_",intron_size),
                        MSR_D,
                        MSR_A) %>%
          mutate(#intron_length = intron_length %>% log10(),
                 intron_length = intron_length/100,
                 MSR_D = (MSR_D * 100000) %>% round(digits = 0),
                 MSR_A = (MSR_A * 100000) %>% round(digits = 0)) %>% 
          distinct(ref_junID,.keep_all = T) %>% 
          drop_na()
  
        if (common_introns) {
          common_annotated_introns <- readRDS(file = file.path(args$results_folder, "common_introns_all_tissues.rds"))
          idb <- idb %>% filter(ref_junID %in% common_annotated_introns$ref_junID)
        }
        
        ## SAVE DATA
        message("Saving 'idb' data from ", cluster_id, "...")
        saveRDS(object = idb,
                file = file.path(folder_path, paste0(cluster_id, "_common_introns", common_introns, "_introns_phastcons", phastcons_type, 
                                                     "_intronsize", intron_size ,".rds")))
      } else {
        
        message("Loading 'idb' data from ", cluster_id, "...")
        idb <- readRDS(file =  file.path(folder_path, 
                                         paste0(cluster_id, "_common_introns", common_introns,"_introns_phastcons", phastcons_type, "_intronsize", intron_size ,".rds"))) %>% drop_na()
        message("Data loaded!")
      }
      
      idb %>% summary
      
      if (!file.exists(file.path(folder_path, 
                                 paste0("ZINB_Donor_", cluster_id, "_common_introns", common_introns, "_introns_phastcons", phastcons_type, "_intronsize", intron_size ,".rds"))) || 
          replace) {
        
        ######################################################################################################
        ## ZERO-INFLATED NEGATIVE BINOMIAL POISSION MODELS ---------------------------------------------------
        ######################################################################################################
        
        message(Sys.time(), " - fitting a ZINB Donor model...")

        zeroinfl_model_D <- pscl::zeroinfl(MSR_D ~
                                             gene_num_transcripts +
                                             protein_coding +
                                             intron_length +
                                             intron_5ss_score +
                                             intron_3ss_score +
                                             mean_CDTS5ss +
                                             mean_CDTS3ss +
                                             mean_phastConsway5ss +
                                             mean_phastConsway3ss,
                                           data = idb)
        zeroinfl_model_D %>% summary %>% print()
        
        # cbind(zeroinfl_model_D$coefficients$count, confint(object =
        # zeroinfl_model_D,  level=0.95) ) %>% as.data.frame() %>%
        # tibble::rownames_to_column("covariate") %>% filter(str_detect(string =
        # covariate, pattern = "count")) %>% dplyr::rename(log.estimate = "V1",
        # log.conf.low = "2.5 %", log.conf.high = "97.5 %") %>% mutate(sign =
        # (sign(x =  (zeroinfl_model_D$coefficients$count) %>% unlist() %>%
        # unname()))) %>% mutate(log.estimate = (log.estimate %>% exp()) * sign,
        # log.conf.low = (log.conf.low %>% exp()) * sign, log.conf.high =
        # (log.conf.high %>% exp()) * sign)
        
        saveRDS(object = zeroinfl_model_D,
                file = file.path(folder_path, paste0("ZINB_Donor_", cluster_id, "_common_introns", common_introns, 
                                                     "_introns_phastcons", phastcons_type, 
                                                     "_intronsize", intron_size ,".rds")))
      }
      
      if (!file.exists(file.path(folder_path, paste0("/ZINB_Acceptor_", cluster_id, "_common_introns", common_introns, "_introns_phastcons", 
                                                     phastcons_type, "_intronsize", intron_size ,".rds"))) || replace) {
        
        message(Sys.time(), " - fitting a ZINB Acceptor model...")
        
        zeroinfl_model_A <- pscl::zeroinfl(MSR_A ~ 
                                             gene_num_transcripts +
                                             protein_coding + 
                                             intron_length +
                                             intron_5ss_score +
                                             intron_3ss_score +
                                             mean_CDTS5ss +
                                             mean_CDTS3ss +
                                             mean_phastConsway5ss +
                                             mean_phastConsway3ss,
                                           data = idb )
        
        zeroinfl_model_A %>% summary %>% print()
        
        saveRDS(object = zeroinfl_model_A,
                file = file.path(folder_path, paste0("/ZINB_Acceptor_", cluster_id, "_common_introns", common_introns, 
                                                     "_introns_phastcons", phastcons_type, "_intronsize", intron_size , ".rds")))
      }
    }
  }
}

main_figure4_d <- function() {
  
  project_id = "BRAIN"
  cluster_id = "Brain - Frontal Cortex (BA9)"
  intron_size = 100
  phastcons_type = 17
  common_introns = F
  
  ##########################
  ## LOAD THE DATA
  ##########################
  
  results_file1 <- file.path(args$data_folder, "figure4_d.csv")
  
  if (file.exists(results_file1)) {
    
    ZINB_plot_count_tidy <- read_csv(file = results_file1, col_names = T, show_col_types = F)
    
  } else {
    
    ZINB_tissue_MSR <- map_df(c("Donor", "Acceptor"), function(novel_type)  { 
      
      # novel_type <-  c("Donor", "Acceptor")[1]
      ## Load the ZIP model
      message(novel_type)
      ZINB_tissue <- readRDS(file = paste0(args$results_folder, "/ZIP/", 
                                           "/ZINB_", novel_type, "_", cluster_id, "_common_introns", common_introns, 
                                           "_introns_phastcons", phastcons_type, "_intronsize", intron_size ,".rds"))
     
      ## Use to obtain robust standard errors for the coefficients in the model
      ZINB_corrected_pvals <- lmtest::coeftest(ZINB_tissue, vcov = sandwich::sandwich)

      # Extracting coefficients
      ZINB_expCoef <- coef(ZINB_corrected_pvals)#ZINB_tissue$coefficients$count
      ZINB_expCoef <- as.data.frame(ZINB_expCoef, ncol = 2) %>%
        tibble::rownames_to_column("covariate") %>%
        mutate(sign = (sign(x =  ZINB_expCoef %>% unlist() %>% unname()))) %>%
        dplyr::select(-covariate)
      
      ## Extracting confidence intervals
      ZINB_confInt <- confint(object = ZINB_corrected_pvals,  level=0.95) %>% 
        as.data.frame() %>%
        tibble::rownames_to_column("covariate") 
      
      ## Join coefficients & confidence intervals
      ZINB_plot <- cbind(ZINB_expCoef, ZINB_confInt) %>%
        dplyr::rename(log.estimate = ZINB_expCoef, log.conf.low = "2.5 %", log.conf.high = "97.5 %") %>%
        mutate(pvalue = c((ZINB_corrected_pvals)[,4] %>% 
                            as.data.frame() %>% 
                            tibble::rownames_to_column("covariate") %>% 
                            pull(.))) %>%
        filter(str_detect(string = covariate, pattern = "count")) %>% 
        mutate(MSR = novel_type)
      
      ZINB_plot %>% return()
      
    })
    
    ###########################
    ## EXPONENTIATE THE LOGS
    ###########################
    
    ## 1. EXP COEFFICIENTS (given it is a log model)
    ZINB_tissue_MSR <- ZINB_tissue_MSR %>%
      mutate(log.estimate = (log.estimate %>% exp()) * sign,
             log.conf.low = (log.conf.low %>% exp()) * sign,
             log.conf.high = (log.conf.high %>% exp()) * sign)
    
    
    ##########################################
    ## TIDY THE COUNT MODEL
    ##########################################
    
    coef_names <- c("Num. transcripts", "Protein coding", "Intron Length", "MES 5'ss score", "MES 3'ss score", "CDTS 5'ss", "CDTS 3'ss", "PhastCons17 5'ss", "PhastCons17 3'ss")
    cov_names <- c("Gene-level","Gene-level", "Intron-level","Intron-level", "Intron-level","Intron-level", "Intron-level","Intron-level", "Intron-level")
    cov_order <- c("C","D","F","G","H","I","J","K","L")
    
    ## Filter the count model
    ZINB_plot_count <- ZINB_tissue_MSR %>%
      filter(str_detect(string = covariate, pattern = "tercept", negate = T),
             str_detect(string = covariate, pattern = "count")) %>%
      mutate(covariate_tidy = c(coef_names,coef_names))
    
    ## Order and tidy covariates
    ZINB_plot_count_tidy <- ZINB_plot_count %>%
      mutate(estimate_lab = paste0(round(x = log.estimate, digits = 2), " (", round(x = log.conf.low, digits = 2), ",", round(x = log.conf.high, digits = 2), ")")) %>%
      dplyr::mutate(p.value = signif(pvalue, digits=3)) %>%
      mutate(p.value = case_when(pvalue == 0 ~ "<2.2e-16", pvalue > 0 ~ formatC(pvalue, format='e',digits = 2) %>% as.character())) %>%
      mutate(covariate_group = c(cov_names, cov_names), covariate_order = c(cov_order, cov_order))
    
    ZINB_plot_count_tidy$covariate_order = factor(ZINB_plot_count_tidy$covariate_order, levels = cov_order)
    
    ## Save data
    write_csv(x = ZINB_plot_count_tidy, file = file.path(args$data_folder, "figure4_d.csv"), col_names = T)
    
  }
  
  
  ##########################################
  ## PLOT COUNT MODEL - by sections
  ##########################################
  
  ZINB_plot_count_tidy >- ZINB_plot_count_tidy %>%
    mutate(q = p.adjust(p = pvalue, method = "fdr")) %>%
    #dplyr::mutate(p.value = signif(q, digits=3)) %>%
    mutate(p.value = case_when(q == 0 ~ "<2.2e-16", q > 0 ~ formatC(q, format='e',digits = 2) %>% as.character()))
  
  ## Generate the center of the plot
  p_mid <- ggplot(data = ZINB_plot_count_tidy %>% filter(covariate_tidy != "Covariate"), 
                  mapping = aes(y = fct_rev(covariate_order), color = MSR)) + 
    theme_classic() +
    geom_point(aes(x=log.estimate, group = MSR), shape=15, size=3, position = position_dodge(width = .75)) +
    geom_linerange(aes(xmin = log.conf.low, xmax = log.conf.high, group = MSR), position = position_dodge(width = .75)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    labs(x="Exponentiated Coefficient", y="") +
    ggforce::facet_col(~covariate_group, scales = "free_y", space = 'free', strip.position = "left") + 
    theme(axis.line.y = element_blank(), axis.ticks.y= element_blank(),
          axis.text.y = element_blank(), axis.title.y= element_blank(),
          strip.text.y = element_blank(), strip.text = element_text(size = 3),
          text = element_text(size = 8),
          legend.position = "top")+
    scale_color_manual(values = c( "#35B779FF", "#8d03b0"), breaks = c( "Donor", "Acceptor"), labels = c( "Novel Donor", "Novel Acceptor"))  + 
    theme(plot.margin = margin(t = 0, r = -5, b = 0, l = -5), legend.box.margin=margin(b = -10, t = -5))
  p_mid 
  
  ## Generate the left section of the plot
  p_left <- ZINB_plot_count_tidy %>%
    mutate(covariate_tidy = ifelse(MSR == "Acceptor", "", covariate_tidy)) %>%
    ggplot(aes(y = fct_rev(covariate_order))) +
    geom_text(aes(x = 0, label = covariate_tidy, group= MSR), hjust = 0, size = 2.5, position = position_dodge(width = .75), fontface = "bold") +
    geom_text(aes(x = 1, label = estimate_lab, group= MSR), hjust = 0, size = 2.5, position = position_dodge(width = .75), fontface = ifelse(ZINB_plot_count_tidy$estimate_lab == "[95% CI]", "bold", "plain")) +
    theme_void() +
    ggforce::facet_col(~covariate_group,scales = "free_y", space = 'free',strip.position = "left") + 
    coord_cartesian(xlim = c(0,2)) +
    theme(title =  element_text(size = 8, colour = "black", face = "bold"), strip.text = element_text(size = 7))
  p_left
  
  
  ## Generate the right section of the plot
  p_right <- ZINB_plot_count_tidy %>%
    ggplot(aes(y = fct_rev(covariate_order))) +
    geom_text(aes(x = 0, y = fct_rev(covariate_order), label = p.value, group = MSR), hjust = 0, position = position_dodge(width = .75), size = 2.5,
               fontface = ifelse(ZINB_plot_count_tidy$p.value == "p-value", "bold", "plain") )  +
    ggforce::facet_col(~covariate_group, scales = "free_y", space = 'free',strip.position = "left") + 
    theme_void() +
    theme(axis.line.y = element_blank(), axis.ticks.y= element_blank(),
          axis.text.y= element_blank(), axis.title.y= element_blank(),
          strip.text.y = element_blank(), title =  element_text(size = 8, colour = "black", face = "bold"))
  p_right 
  
  
  library(patchwork)
  layout <- c(
    area(t = 0, l = 0, b = 30, r = 4), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
    area(t = 1, l = 5, b = 30, r = 9), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
    area(t = 0, l = 9.7, b = 30, r = 10)
  )
  # final plot arrangement
  p_left + ggtitle("Covariate") + 
    p_mid + ggtitle("[95% CI]") +
    p_right + ggtitle("p-value") + 
    patchwork::plot_layout(design = layout)  
  
  
  
  ## Save Plot (Main Figure 4.d)
  figure_name <- paste0(args$figures_folder, "/main_figure4d")
  ggplot2::ggsave(filename = paste0(figure_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(figure_name, ".svg"), width = 180, height = 90, units = "mm", dpi = 300)
  
}




## 8. Accuracy in splicing is affected by RNA-binding protein (RBP) expression changes

get_common_introns_across_tissues <- function () {
  
  all_projects <- df_metadata$SRA_project %>% unique()
  
  ################################
  ## CONNECT TO THE DATABASE
  ################################
  
  print(paste0(Sys.time(), " - getting unique and common junctions across clusters..."))
  
  
  ## Getting all introns that are common across all GTEx tissues -------------------------------------------
  
  all_introns <- list()
  
  for (project_id in all_projects) {
    
    message(project_id)
    
    ## GET THE CLUSTERS
    all_clusters <- df_metadata %>% filter(SRA_project == project_id) %>% distinct(cluster) %>% pull()

    for(cluster_id in all_clusters) { 
      
      query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      all_introns[[cluster_id]] <- introns$ref_junID
      
      message("Introns collected from '", cluster_id, "'")
    }
    
  }
  
  common_introns <- data.frame(ref_junID = Reduce(intersect,  all_introns))
  common_introns %>% head()
  common_introns %>% nrow()
  
  query <- paste0("SELECT * FROM 'intron' WHERE ref_junID IN (", paste(common_introns$ref_junID, collapse = ","), ")")
  introns <- dbGetQuery(con, query) %>% as_tibble()

  saveRDS(object = introns, file = paste0(args$results_folder, "/common_introns_all_tissues.rds"))
  
}

get_data_main_figure5_ab <- function() {
  
  intron_size = 100
  phastcons_type = 17
  common = F
  
  if (file.exists(file.path(args$results_folder, paste0("data_figure5_ab",common,".rds")))) {
    all_coefficient_tissues <- readRDS(file = file.path(args$results_folder, paste0("data_figure5_ab", common, ".rds")))
    
  } else {
    
    all_coefficient_tissues <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      message(project_id)
      
      ## GET THE CLUSTERS
      all_clusters <- df_metadata %>% filter(SRA_project == project_id) %>% distinct(cluster) %>% pull()
      
      map_df(all_clusters, function(cluster_id) { 
        
        # cluster_id <- all_clusters[1]
        
        map_df(c("Donor", "Acceptor"), function(novel_type) { 
          
          # novel_type <- "Donor"
          if (file.exists(paste0(args$results_folder, "/ZIP/ZINB_",novel_type,"_",cluster_id,
                                  "_common_introns", common, "_introns_phastcons", phastcons_type, "_intronsize", intron_size ,".rds")))  {
            
            message(novel_type, " - ", cluster_id)
            
            ZINB_tissue <- readRDS(file = paste0(args$results_folder, "/ZIP/ZINB_", novel_type,"_", cluster_id, 
                                                 "_common_introns", common, "_introns_phastcons", phastcons_type, "_intronsize", intron_size ,".rds"))

            message("Extracting coefficients...")
            ZINB_corrected_pvals <- lmtest::coeftest(ZINB_tissue, vcov = sandwich::sandwich)
            ZINB_expCoef <- ZINB_corrected_pvals$coefficients$count #stats::coef(object = ZINB_tissue$coefficients$count) 
            
            message("Returning dataframe...")
            as.data.frame(ZINB_expCoef, ncol = 2) %>%
              tibble::rownames_to_column("covariate") %>%
              mutate(cluster = cluster_id,
                     pvalue = (ZINB_tissue %>% summary)[1]$coefficients$count[,4],
                     MSR = novel_type) %>%
              return()
            
          } else {
            stop("Source data not found!")
          }
        })
      })
    })
    saveRDS(all_coefficient_tissues, file = file.path(args$results_folder, paste0("data_figure5_ab",common,".rds")))
  }
  
  
  ## Generate one plot per loop to produce Figure5 a & b
  
  for (MSR_type in c("Donor", "Acceptor")) {
    
    # MSR_type <- "Donor"
    # MSR_type <- "Acceptor"
    
    message(MSR_type)
    
    ## 1. Load the estimate variance across tissues
    all_coefficient_tissues_fdr <- all_coefficient_tissues %>%
      mutate(beta_sign = (sign(x = ZINB_expCoef))) %>%
      mutate(ZINB_expCoef = (ZINB_expCoef %>% exp()) * beta_sign) %>%
      filter(MSR == MSR_type) %>%
      mutate(ZINB_expCoef = ZINB_expCoef %>% as.double,
             q = p.adjust(pvalue, method = "fdr"))
    
    ## 2. only keep significant covariates
    all_coefficient_tissues_fdr <- all_coefficient_tissues_fdr %>%
      filter(q <= 0.05 | (is.nan(q))) %>%
      dplyr::select(-c(q, pvalue))
    
    
    ################################################################
    ## GET THE COUNT MODEL AND TIDY DATA FOR PLOTTING
    ################################################################
    
    all_coefficient_tissues_tidy <- all_coefficient_tissues_fdr %>%
      filter(str_detect(string = covariate, pattern = "Intercept", negate = T)) %>%
      dplyr::group_by(covariate) %>%
      spread(key = covariate, value = ZINB_expCoef) %>%
      ungroup()
    
    all_coefficient_tissues_tidy <- all_coefficient_tissues_tidy %>% 
      dplyr::rename("Gene num. transcripts" = "gene_num_transcripts",
                    "Intron Length" = "intron_length",
                    "Intron 5'ss MES score" = "intron_5ss_score",
                    "Intron 3'ss MES score" = "intron_3ss_score", 
                    "CDTS 5'ss" = "mean_CDTS5ss",
                    "CDTS 3'ss" = "mean_CDTS3ss", 
                    "PhastCons17 5'ss" = "mean_phastConsway5ss",
                    "PhastCons17 3'ss" = "mean_phastConsway3ss",
                    "Protein coding" = "protein_coding") %>%
      gather(feature, coefficient,-cluster) %>% as_tibble()
    
    label_faceting <- data.frame(
      label = c("Gene num. transcripts", "Protein coding", "Intron Length", "Intron 5'ss MES score", "Intron 3'ss MES score", "CDTS 5'ss", "CDTS 3'ss", "PhastCons17 5'ss", "PhastCons17 3'ss" ),
      subgroup = c("Gene\nLevel", "Gene\nLevel", "Intron Level", "Intron Level", "Intron Level", "Intron Level", "Intron Level", "Intron Level", "Intron Level")
      )
    
    all_coefficient_tissues_tidy <- all_coefficient_tissues_tidy %>%
      inner_join(y = label_faceting, by = c("feature" = "label")) %>% drop_na() %>%
      mutate(coefficient = coefficient %>% as.double())
    
    
    ################################################################
    ## Save source data 
    ################################################################
    
    file_name <- if (MSR_type == "Donor") { "figure5_a.csv" } else { "figure5_b.csv" }
    write_csv(x = all_coefficient_tissues_tidy, file = file.path(args$data_folder, file_name), col_names = T)
    
  }
  
}

main_figure5_ab <- function() {
  
  ## LOAD RESULTS
  results_file1 <- file.path(args$data_folder, "figure5_a.csv")
  results_file2 <- file.path(args$data_folder, "figure5_b.csv")
  
  if (!file.exists(results_file1) && !file.exists(results_file2)) { get_data_main_figure5_ab() } 
  
  if (file.exists(results_file1) && file.exists(results_file2)) {
    
    for (MSR_type in c("Donor", "Acceptor")) {
      
      if (MSR_type == "Donor") { 
        all_coefficient_tissues_tidy <- read_csv(file = results_file1, col_names = T, show_col_types = F) 
        color_boxplot <- "#35B779FF" 
        figure_name <- "main_figure5a"
      } else { 
        all_coefficient_tissues_tidy <- read_csv(file = results_file2, col_names = T, show_col_types = F)
        color_boxplot <- "#8d03b0"
        figure_name <- "main_figure5b"
      }
      
      ## PLOT
      all_coefficient_tissues_tidy$feature <- factor(all_coefficient_tissues_tidy$feature, 
                                                     levels = c("Gene num. transcripts", "Protein coding", "Intron Length", "Intron 5'ss MES score", "Intron 3'ss MES score", 
                                                                "CDTS 5'ss", "CDTS 3'ss", "PhastCons17 5'ss", "PhastCons17 3'ss") %>% rev())
      
      plotTissuesZINB <- ggplot(data = all_coefficient_tissues_tidy, aes(x = feature, y = coefficient, fill = feature)) + 
        geom_boxplot(fill = color_boxplot, 
                     aes(ymin = min(coefficient),
                         ymax = max(coefficient), 
                         lower = quantile(coefficient, .25), 
                         upper = quantile(coefficient, .75),
                         middle = median(coefficient))) +
        coord_flip() + 
        facet_grid(vars(subgroup), scales = "free", switch = "y", space = "free_y") +
        ylab("Distribution of the exponentiated coefficients across tissues (q<0.05)") +
        xlab(" ") +
        theme_light() +
        custom_ggtheme +
        scale_fill_manual(breaks = c("MSR_Donor","MSR_Acceptor"), labels = c("MSR_Donor","MSR_Acceptor")) +
        theme(axis.text.y = element_text(vjust = 0.5, hjust = 1)) +
        geom_hline(yintercept = 0,linetype='dotted')
      
      plotTissuesZINB
      
      dir.create(path = file.path(args$figures_folder, "ZIP"), recursive = T, showWarnings = F)
      file_name <- paste0(args$figures_folder, "/", figure_name)
      ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 60, units = "mm", dpi = 300)
      ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 60, units = "mm", dpi = 300)
    }
  }
}




## SECTION 2 - SUPPLEMENTARY FIGURES ---------------------------------------

supplementary_figure2_3 <- function () {
  
  
  df_jxn_sharing_all_tissues <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[9]
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(all_clusters, function(cluster_id) {
      
      # cluster_id <- all_clusters[1]
      ## Query the database
      
      if ( any(str_detect(tables, cluster_id)) ) {
        
        print(paste0(Sys.time(), " - ", cluster_id))
        
        query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue")
        db_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
        
        query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced' AS tissue")
        db_not_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
        
        plyr::rbind.fill(db_misspliced_introns, db_not_misspliced_introns) %>%
          mutate("project_id" = project_id, "cluster_id" = cluster_id) %>%
          return()
      }
    })
  }) %>%
    mutate(project_id = str_replace_all(string = project_id, pattern = "_", replacement = " "))
  
  ################################
  ## COMMON PLOT FUNCTION
  ################################
  
  generate_plot <- function(df, ylabel) {
    ggplot(data = df, aes(x = n_individuals)) +
      geom_histogram(boundary = 1, closed = "left")  +
      scale_x_continuous(breaks = c(1, 200, 400, 600, 800)) +
      facet_wrap(~project_id) +
      theme(
        legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
      ) +
      xlab("Number of individuals") +
      ylab(paste0(ylabel)) +
      theme_light() +
      custom_ggtheme
  }
  
  ################################
  ## Supplementary Figure 2
  ################################
  
  novel_data <- df_jxn_sharing_all_tissues %>%
    group_by(project_id) %>%
    distinct(novel_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::rename(n_individuals = novel_n_individuals) %>% 
    drop_na()
  
  novel_plot <- generate_plot(novel_data, ylabel = "Number of unique novel junctions")
  
  ggplot2::ggsave(plot = novel_plot,
                  filename = file.path(args$figures_folder, "supplementary_figure2.png"), width = 180, height = 180, units = "mm", dpi = 300)
  
  #################################
  ## Supplementary Figure 3
  #################################
  
  intron_data <- df_jxn_sharing_all_tissues %>%
    group_by(project_id) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup()  %>%
    dplyr::rename(n_individuals = ref_n_individuals) %>%
    drop_na()
  
  intron_plot <- generate_plot(intron_data, ylabel = "Number of unique annotated introns")
  
  ggplot2::ggsave(plot = intron_plot,
                  filename = file.path(args$figures_folder, "supplementary_figure3.png"), width = 180, height = 180, units = "mm", dpi = 300)
  
}

supplementary_figure4 <- function() {
  
  ## After accounting for sample number, we found that the highest numbers of unique 
  ## novel junctions were identified in X tissue with the lowest numbers in Y tissue, 
  
  db_metadata_tidy <- df_metadata %>%
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
    geom_bar(mapping = aes(x = cluster,  y = avg_n_novel_sample), stat = "identity") +
    ylab("Average number of unique novel\njunctions across samples") +
    xlab("GTEx tissue") +
    theme_light() +
    custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  
  ggplot2::ggsave(filename = file.path(args$figures_folder, "/supplementary_figure4.png"), width = 180, height = 100, units = "mm", dpi = 300)
}

supplementary_figure5 <- function () {
  

  ## CONNECT TO THE DATABASE
  
  project_id <- "BRAIN"
  cluster <- "Brain - Frontal Cortex (BA9)"
  
  reference_gtf_version <- 105
  gtf_versions <- c("76", "81", "90", "97", "104")
  dates <- c("2014", "2015", "2017", "2019", "2021")
  
  ## GET RECLASSIFICATION DETAILS
  
  df_reclassification_rates <- if ( file.exists( file.path(base_folder, "results", main_project, reference_gtf_version, project_id, "supplementary_figure5.rds")) ) {
    
    readRDS(file = paste0(getwd(), "/results/", project_id, "/v", reference_gtf_version, "/", main_project, "/supplementary_figure5.rds"))
    
  } else {
    
    map_df(gtf_versions, function(gtf_version) {
      
      # version <- gtf_versions[1]
      print(paste0("v", gtf_version))
      
      if ( file.exists(paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/",
                              project_id, "_", cluster, "_all_split_reads.rds")) ) {
        
        
        ## ENSEMBL old
        df_old <-  readRDS(file =  paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/",
                                          project_id, "_", cluster, "_all_split_reads.rds")) %>% as_tibble()
        db_introns_old <- df_old %>% filter(type == "annotated")
        db_novel_old <- df_old %>% filter(type %in% c("novel_donor", "novel_acceptor"))
        
        ## ENSEMBL v105
        df_new <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v105/", main_project, "/base_data/",
                                        project_id, "_", cluster, "_all_split_reads.rds") ) %>% as_tibble()
        
        db_introns_new <- df_new %>%
          filter(type == "annotated")
        db_novel_new <- df_new %>%
          filter(type %in% c("novel_donor", "novel_acceptor"))
        
        
        rm(df_old)
        rm(df_new)
        
        
        in_annotation <- db_novel_old %>%
          filter(junID %in% db_introns_new$junID) %>% 
          distinct(junID, .keep_all = T)
        in_annotation %>% nrow() %>% print()
        
        out_annotation <- db_introns_old %>%
          filter(!(junID %in% db_introns_new$junID)) %>% 
          distinct(junID, .keep_all = T) 
        out_annotation %>% nrow() %>% print()
        
        label <- NULL
        
        if (gtf_version == "76") {
          label <- paste0("v", gtf_version, "\n(", dates[1], ")")
        } else if (gtf_version == "81") {
          label <- paste0("v", gtf_version, "\n(", dates[2], ")")
        } else if (gtf_version == "90") {
          label <- paste0("v", gtf_version, "\n(", dates[3], ")")
        } else if (gtf_version == "97") {
          label <- paste0("v", gtf_version, "\n(", dates[4], ")")
        } else {
          label <- paste0("v", gtf_version, "\n(", dates[5], ")")
        }
        
        
        return(data.frame(tissue = cluster,
                          reclassification_rates = label,
                          in_annotation = (in_annotation %>% nrow() ) / db_novel_old %>% nrow(),
                          out_annotation = (out_annotation %>% nrow() ) / db_introns_old %>% nrow()))
      } else {
        return(NULL)
      }
      
      
    })
    
    file_name <- paste0(getwd(), "/results/", project_id, "/v", reference_gtf_version, "/", main_project,
                        "/supplementary_figure5.rds")
    saveRDS(object = df_reclassification_rates, file = file_name)
    
    
  } 
  
  ################################
  ## PREPARE DATA BEFORE PLOT
  ################################
  
  
  ## Tidy the dataframe prior display
  df_reclassification_rates$reclassification_rates = factor(df_reclassification_rates$reclassification_rates, 
                                                levels = c("v76\n(2014)", "v81\n(2015)", "v90\n(2017)",
                                                           "v97\n(2019)", "v104\n(2021)"))
  
  df_reclassification_rates <- df_reclassification_rates %>%
    mutate(in_annotation = in_annotation / 100,
           out_annotation = out_annotation / 100)
  
  ################################
  ## PLOT
  ################################
  
  
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data=df_reclassification_rates, aes(x = reclassification_rates, y = in_annotation, group=1)) +
    geom_line(linetype = "dashed", col = "red") +
    geom_point() +
    xlab(NULL) +
    ylab("re-classification rates\nas compared to Ensembl v105") +
    theme_light() +
    custom_ggtheme +
    guides(fill = "none")
  
  ggplot2::ggsave(filename = file.path(args$figures_folder, "supplementary_figure5.png"),  width = 180, height = 70, units = "mm", dpi = 300)
  
}

supplementary_figure8 <- function() {
  
  file_name <- file.path(args$results_folder, "df_mes.rds")
  
  df_mes <-  if ( file.exists(file_name) )  {
    print("Loading MES results...")
    readRDS(file = file_name) %>% distinct(ref_junID, novel_junID, .keep_all = T) %>% as_tibble()
  }else {
    stop("Please, run function 'main_figure3_ab' first")
  }
  
  df_5ss <- df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron = ref_mes5ss, novel_donor = novel_mes5ss) %>%
    gather(key = "junction_type", value = "ss5score")
  df_5ss %>% dplyr::count(junction_type)
  
  df_3ss <- df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron = ref_mes3ss, novel_acceptor = novel_mes3ss) %>%
    gather(key = "junction_type", value = "ss3score")
  df_3ss %>% dplyr::count(junction_type)
  
  ss5plot <- ggplot(df_5ss, aes(ss5score, fill = junction_type)) +
    geom_density(alpha = 0.9) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    xlab("MES Donor (5'ss) score") +
    scale_fill_manual(values = c("#999999", "#35B779FF"), breaks=c("intron", "novel_donor"), labels=c("Annotated Intron  ", "Novel Donor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = element_blank(), ncol = 2, nrow = 1))
  
  
  ss3plot <- ggplot(df_3ss, aes(ss3score, fill = junction_type)) +
    geom_density(alpha = 0.9) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    scale_fill_manual(values = c("#999999", "#64037d"), breaks = c("intron", "novel_acceptor"), labels = c("Annotated Intron  ", "Novel Acceptor")) +
    xlab("MES Acceptor (3'ss) score") +
    custom_ggtheme +
    guides(fill = guide_legend(title = element_blank(), ncol = 2, nrow = 1))
  
  
  ## COMBO
  ggpubr::ggarrange(ss5plot, ss3plot, labels = c("a", "b"), ncol = 2,  nrow = 1)
  
  file_name <- file.path(args$figures_folder, "/supplementary_fig3.png")
  ggplot2::ggsave(filename = file_name, width = 180, height = 80, units = "mm", dpi = 300)
  
}

supplementary_figure9 <- function() {
  
  
  #################################################
  ## LOAD AND PREPARE THE DATA
  #################################################
  
  file_name <- paste0(args$results_folder, "/df_mes.rds")
  df_mes <- readRDS(file = file_name) %>% as_tibble()
  
  df_mes_tidy <- df_mes %>%
    distinct(ref_junID, novel_junID, .keep_all = T) %>%
    left_join(y = master_novel_junctions %>% dplyr::select(novel_junID, distance), by = "novel_junID") %>%
    filter(abs(distance) <= 75) %>%
    mutate(mod3 = abs(distance)%%3)
  
  df_mes_tidy$novel_type = factor(df_mes_tidy$novel_type, levels = c("novel_donor", "novel_acceptor"))
  df_mes_tidy$novel_type = factor(df_mes_tidy$novel_type, levels = c("novel_donor", "novel_acceptor"))
  
  df_mes_tidy %>% dplyr::count(mod3)
  #################################################
  ## MODULO STATS
  #################################################
  
  df_mes_tidy_plot <- df_mes_tidy %>%
    dplyr::select( mod3, 
                   "Delta MES 5ss" = diff_ss5score, 
                   "Delta MES 3ss" = diff_ss3score, 
                   novel_type, 
                   distance) %>%
    mutate(mod3 = mod3 %>% as.factor()) %>%
    gather(type, score, -novel_type, - mod3, -distance) %>%
    filter((novel_type == "novel_donor" & type == "Delta MES 5ss") |
              (novel_type == "novel_acceptor" & type == "Delta MES 3ss")) %>%
    group_by(type,mod3) %>% 
    mutate(mean_score = mean(score)) %>% 
    ungroup()
  
  group_means  <- df_mes_tidy_plot %>%
    group_by(type,mod3) %>%
    summarise(mean = mean(score)) %>%
    ungroup()
  
  
  # p-values --------------------------------------------------------------------------------------------------
  stat.test.donor <- tibble::tribble(
    ~group1, ~group2,   ~p.adj,
    "0",     "1-2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 == 0) %>% pull(diff_ss5score),
                                 y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 > 0) %>% pull(diff_ss5score),
                                 alternative = "less",
                                 correct = T))$p.value %>% 
      formatC(format = "e", digits = 2)
  )
  stat.test.acceptor <- tibble::tribble(
    ~group1, ~group2,   ~p.adj,
    "0",     "1-2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 == 0) %>% pull(diff_ss3score),
                                 y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 > 0) %>% pull(diff_ss3score),
                                 alternative = "less",
                                 correct = T))$p.value  %>% formatC(format = "e", digits = 2)
  )
  
  
  # eff.sizes --------------------------------------------------------------------------------------------------
  
  eff.size.donor <- rstatix::wilcox_effsize(data = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor") %>%
                                              mutate(mod3 = mod3 %>% as.factor()),
                                            formula = diff_ss5score ~ mod3 )
  
  eff.size.acceptor <- rstatix::wilcox_effsize(data = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor") %>%
                                                 mutate(mod3 = mod3 %>% as.factor()),
                                               formula = diff_ss3score ~ mod3 )
  
  eff.size.acceptor
  
  
  ####################################################
  ## BOXPLOT
  ####################################################
  
  ggplot(df_mes_tidy_plot %>%  mutate(mod3 = ifelse(mod3 == "0", "0", "1-2"))) +
    geom_boxplot(aes(y = score, x = mod3, fill = mod3), alpha = 0.9)  +
    ggpubr::stat_pvalue_manual(
      rbind(stat.test.donor %>% mutate(type = "Delta MES 5ss"),
            stat.test.acceptor %>% mutate(type = "Delta MES 3ss")) %>%
        mutate(p.scient = format(p.adj, scientific = TRUE)), 
      y.position = 35, 
      size = 3,
      step.increase = 0.1,
      label = "p.scient"
    ) +
    facet_grid(~fct_rev(type)) +
    theme_light() +
    custom_ggtheme +
    ylab("Delta MES score") +
    xlab("Modulo3") +
    theme(legend.title = element_text(size = "7", colour = "black")) +
    guides(fill = guide_legend(title = "Modulo3")) +
    scale_fill_discrete(breaks = c("0", "1-2")) 
  
  
  ggplot2::ggsave(filename = file.path(args$figures_folder, "supplementary_figure9.png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
}

supplementary_figure10 <- function()  {
  
  ##############################################################################
  ## Supplementary Figure 10
  ## INTRON SUBSAMPLING BY MEAN EXPRESSION LEVELS
  ##############################################################################
  
  ## Load data from FCTX
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", project_id, "_nevermisspliced' ")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", project_id, "_misspliced' ")
  introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
  
  ## Add information from the master table
  introns <- introns %>% left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding) %>% as_tibble(), by = "ref_junID") %>% as_tibble() 
  message(Sys.time(), " - ", introns %>% nrow(), " - introns!")
  
  # Tidy dataframe
  df_biotype_result_tidy <- introns %>% filter(protein_coding %in% c(0, 100)) %>%
    mutate(biotype = ifelse(protein_coding == 100, "PC", "non PC"),
           mean_coverage = log10(ref_sum_counts/ref_n_individuals)) %>%
    distinct(ref_junID, .keep_all = TRUE)
  
  # Plot (a) Before subsampling
  plot_BS <- ggplot(df_biotype_result_tidy %>% distinct(ref_junID, .keep_all = TRUE)) +
    geom_density(aes(x = mean_coverage, fill = biotype), alpha = 0.9) +
    ggtitle("Before subsampling") + xlab("log10 mean expression level") +
    theme_light() + custom_ggtheme +
    ggsci::scale_fill_npg(name = "Transcript biotype: ")
  
  # Subsampling
  df_lncRNA_tidy <- df_biotype_result_tidy %>% filter(biotype == "non PC")
  df_protein_coding_tidy <- df_biotype_result_tidy %>% filter(biotype == "PC")
  data_combined <- rbind(df_protein_coding_tidy, df_lncRNA_tidy)
  
  data_combined %>% dplyr::count(biotype)
  
  # Check if subsampling file exists
  subsample_file <- paste0(args$results_folder, "/supplementary_figure10.rds")
  
  if (!file.exists(subsample_file)) {
    # Subsample using MatchIt
    m.out <- MatchIt::matchit(biotype ~ mean_coverage, data = data_combined, 
                              distance = data_combined$mean_coverage,
                              method = "nearest", caliper = c(mean_coverage = 0.005))
    subsample <- MatchIt::match.data(m.out) # %>% distinct(ref_junID, .keep_all = TRUE)
    saveRDS(subsample, subsample_file)
  } else {
    subsample <- readRDS(subsample_file)
    
  }
  subsample %>% dplyr::count(biotype)
  
  # Plot (b) After subsampling
  plot_AS <- ggplot(subsample) +
    geom_density(aes(x = mean_coverage, fill = biotype), alpha = 0.9) +
    ggtitle("After subsampling") + xlab("log10 mean expression level") +
    theme_light() + custom_ggtheme +
    ggsci::scale_fill_npg(name = "Transcript biotype: ")
  
  
  # Combine and save plots
  combined_plot <- ggpubr::ggarrange(plot_BS, plot_AS, labels = c("a", "b"), common.legend = TRUE)
  ggsave(filename = file.path(args$figures_folder, "supplementary_figure10.png"), combined_plot, width = 180, height = 90, units = "mm", dpi = 300)
  
  return(subsample)
  # ############################################################################################
  # 
  # 
  # ## Analysing differences in MSR_D and MSR_A
  # 
  # ## subset introns from lncRNA transcripts after subsampling
  # subsample_lncRNA <- subsample %>%
  #   dplyr::filter(biotype=="non PC")
  # ## subset introns from protein-coding transcripts after subsampling
  # subsample_protein_coding <- subsample %>%
  #   dplyr::filter(biotype=="PC")
  # ## QC
  # if ( intersect(subsample_protein_coding$ref_junID, subsample_lncRNA$ref_junID) %>% length() > 0 )  {
  #   print("ERROR! some introns have been categorised as lncRNA and protein-coding.")
  # }
  
  
  # ##############################################################################
  # ## Figure 4 c & d
  # ## CORRECT BY PER-INTRON MEAN EXPRESSION LEVELS
  # ##############################################################################


  ## PREPARE DATA TO PLOT AND GET STATS

  df_introns_biotype <- subsample #%>%
    #filter(protein_coding %in% c(0,100)) %>%
    #mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC"))

  any(df_introns_biotype %>% pull(ref_junID ) %>% duplicated())
  df_introns_biotype$biotype = factor(df_introns_biotype$biotype, levels = c("non PC","PC"))


  df_introns_biotype <- df_introns_biotype %>%
    dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)


  df_introns_biotype$MSR_type = factor(df_introns_biotype$MSR_type, levels = c("MSR_A", "MSR_D"))
  print(paste0(Sys.time(), " - ", df_introns_biotype %>% nrow(), " - introns after tidying!"))


  df_introns_biotype %>% filter(biotype == "PC") %>% pull(MSR) %>% summary()
  df_introns_biotype %>% filter(biotype == "non PC") %>% pull(MSR) %>% summary()


  ###########################
  ## PLOT BY MSR GROUP
  ###########################

  df_introns_biotype_tidy <- df_introns_biotype %>%
    group_by(biotype, MSR_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    mutate(percentile_group = case_when(MSR == 0 ~ "0",
                                        MSR > 0 & MSR <= 0.2 ~ "(0,.2]",
                                        MSR > 0.2 & MSR <= 0.4 ~ "(.2,.4]",
                                        MSR > 0.4 & MSR <= 0.6 ~ "(.4,.6]",
                                        MSR > 0.6 & MSR <= 0.8 ~ "(.6,.8]",
                                        MSR > 0.8 & MSR <= 1 ~ "(.8,1]"))

  df_introns_biotype_tidy$MSR_type = factor(df_introns_biotype_tidy$MSR_type, levels = c("MSR_D","MSR_A"))
  df_introns_biotype_tidy$percentile_group = factor(df_introns_biotype_tidy$percentile_group, levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  df_introns_biotype_tidy$biotype = factor(df_introns_biotype_tidy$biotype, levels = c("PC","non PC"))
  
  
  
  # #----------------------------------------------------------------------------
  # ## Save SOURCE DONOR data
  # write_csv(x = df_introns_biotype_tidy %>% filter(MSR_type == "MSR_D"),
  #           file = file.path(args$data_folder, "figure4_b.csv"), col_names = T)
  # ## Save SOURCE ACCEPTOR data
  # write_csv(x = df_introns_biotype_tidy %>% filter(MSR_type == "MSR_A"),
  #           file = file.path(args$data_folder, "figure4_c.csv"), col_names = T)
  # #----------------------------------------------------------------------------
  
  
  

  plotMSR_donor <- ggplot(data = df_introns_biotype_tidy %>% filter(MSR_type == "MSR_D")) +
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    ggtitle("MSR Donor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    #ylim(c(0,23000))+
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding   ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = "",
                               ncol = 2, nrow = 1))
  plotMSR_donor


  ## Non-protein-coding
  plotMSR_acceptor <- ggplot(data = df_introns_biotype_tidy %>% filter(MSR_type == "MSR_A")) +
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = 0.5, color = "#333333")+
    ggtitle("MSR Acceptor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"), labels = c("Protein-coding ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(ncol = 2, nrow = 1))


  ggpubr::ggarrange(plotMSR_donor, plotMSR_acceptor,
                    nrow = 1, ncol = 2, common.legend = T)#,
  #labels = c("b", "c"))
  # 
  # 
  # 
  # 
  # print(paste0(Sys.time(), " - saving plot..."))
  # ggplot2::ggsave(file.path(args$figures_folder, "/panel4bc.png"), width = 180, height = 90, units = "mm", dpi = 300)
  # 
  # 
  # 
  # ## PRINT THE ZOOMED VERSION OF FIGURE 4 B & C -------------------------------------
  # 
  # plotMSR_donor_zoomed <- ggplot(data = df_introns_biotype_tidy %>% 
  #                                  filter(MSR_type == "MSR_D")) + 
  #   geom_bar(aes(x = percentile_group, fill = biotype),
  #            position = "dodge", linewidth = .5, color = "#333333")+
  #   ggtitle("MSR Donor") +
  #   xlab("Mis-splicing ratio value group") +
  #   ylab("") +
  #   scale_y_continuous(limits =c(0,800), position = "right") +
  #   
  #   scale_fill_manual(values = c("#333333","#999999"),
  #                     breaks = c("PC","non PC"),
  #                     labels = c("Protein-coding","Non-protein-coding")) +
  #   theme_light() +
  #   custom_ggtheme +
  #   guides(fill = guide_legend(title = NULL,  
  #                              ncol = 2, nrow = 1)) 
  # plotMSR_donor_zoomed
  # 
  # 
  # ## Non-protein-coding
  # plotMSR_acceptor_zoomed <- ggplot(data = df_introns_biotype_tidy %>% 
  #                                     filter(MSR_type == "MSR_A")) + 
  #   geom_bar(aes(x = percentile_group, fill = biotype),
  #            position = "dodge", linewidth = .5, color = "#333333")+
  #   ggtitle("MSR Acceptor") +
  #   xlab("Mis-splicing ratio value group") +
  #   ylab("") +
  #   scale_y_continuous(limits =c(0,800), position = "right") +
  #   
  #   scale_fill_manual(values = c("#333333","#999999"),
  #                     breaks = c("PC","non PC"),
  #                     labels = c("Protein-coding","Non-protein-coding")) +
  #   theme_light() +
  #   custom_ggtheme +
  #   guides(fill = guide_legend(title = NULL, 
  #                              ncol = 2, nrow = 1))
  # plotMSR_acceptor_zoomed
  # 
  # 
  # ggpubr::ggarrange(plotMSR_donor_zoomed,
  #                   plotMSR_acceptor_zoomed,
  #                   nrow = 1,
  #                   ncol = 2,
  #                   common.legend = T)#,
  # 
  # 
  # print(paste0(Sys.time(), " - saving plot..."))
  # file_name <- paste0(args$figures_folder, "/panel4bc-zoomed")
  # ggplot2::ggsave(paste0(file_name, ".png"), 
  #                 width = 180, height = 90, units = "mm", dpi = 300)
  
  
  ####################################
  ## STATS
  ####################################
  
  
  # ## 1. We started by focusing on frontal cortex and observed that while splicing errors are detected 
  # ## infrequently with the MSRD and MSRA values highly skewed towards low values, there was considerable 
  # ## variation across introns (MSRA interquartile range = XXXXX; MSRD interquartile range = XXXX)
  # 
  # introns %>%
  #   pull(MSR_D) %>% 
  #   IQR()
  # 
  # 
  # introns %>%
  #   pull(MSR_A) %>% 
  #   IQR()
  # 
  # 
  # 
  # ## 3. Given that NMD activity would be expected to reduce the detection of splicing errors amongst mRNA transcripts, 
  # ## we separately assessed mis-splicing of annotated introns solely used in protein-coding transcripts (N=35,726 based 
  # ## on Ensembl v105) and those that were exclusively used in non-coding transcripts (N=8,945 based on Ensembl v105).
  # 
  # introns %>% filter(protein_coding == 100) %>% distinct(ref_junID) %>% nrow()
  # introns %>% filter(protein_coding == 0) %>% distinct(ref_junID) %>% nrow()
  # 
  # df_introns_biotype %>% filter(biotype == "PC") %>% distinct(ref_junID) %>% nrow()
  # df_introns_biotype %>% filter(biotype == "non PC") %>% distinct(ref_junID) %>% nrow()
  # 
  # 
  # ## We found that mis-splicing was more frequent amongst annotated introns from non-protein-coding transcripts at both the 5’ss
  # ## 2. Splicing noise is more frequent in introns from non-protein coding transcripts 
  # ## than in introns from protein-coding transcripts even after correcting by mean read coverage
  # 
  # ## splicing noise at the donor is more frequent in introns from NPC vs PC at the 5'ss
  # wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_D") %>% pull(MSR),
  #             y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_D") %>% pull(MSR),
  #             alternative = "less",
  #             paired = T,
  #             correct = T)
  # rstatix::wilcox_effsize(data = df_introns_biotype %>% 
  #                           filter(MSR_type == "MSR_D")%>%
  #                           mutate(biotype = biotype %>% as.factor()),
  #                         formula = MSR ~biotype,
  #                         paired = T)
  # 
  # ## splicing noise at the donor is more frequent in introns from NPC vs PC at the 3'ss
  # wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_A") %>% pull(MSR),
  #             y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_A") %>% pull(MSR),
  #             alternative = "less",
  #             paired = T,
  #             correct = T)
  # rstatix::wilcox_effsize(data = df_introns_biotype %>% 
  #                           filter(MSR_type == "MSR_A")%>%
  #                           mutate(biotype = biotype %>% as.factor()),
  #                         formula = MSR ~ biotype,
  #                         paired = T)
  
}

supplementary_figure11 <- function() {
  
  project_id1 <- "SKIN"
  cluster_id1 <- "Skin - Sun Exposed (Lower leg)"
  project_id2 <- "SKIN"
  cluster_id2 <- "Skin - Not Sun Exposed (Suprapubic)"
  stats=F
  
  tables <- c(paste0(cluster_id1, "_", project_id1), paste0(cluster_id2, "_", project_id2))
  supplementary_figure11_path <- file.path(args$results_folder, "supplementary_figure11_after_subsampling.rds")
  
  if (!file.exists(supplementary_figure11_path)) {
    
    database_introns <- map_df(tables, function(table) {
      
      message(Sys.time(), " - ", table)
      
      # Fetch data from database tables
      df_introns <- run_sql_query(con, table, columns = "ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals, transcript_id", "_nevermisspliced")
      df_introns <- bind_rows(df_introns, run_sql_query(con, table, columns = "ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals, transcript_id", "_misspliced"))
      
      # Join with introns, transcripts, and genes
      ref_junID_str <- paste(df_introns$ref_junID, collapse = ",")
      df_introns <- df_introns %>%
        left_join(dbGetQuery(con, paste0("SELECT DISTINCT ref_junID, ref_length, ref_mes5ss, ref_mes3ss, protein_coding FROM 'intron' WHERE ref_junID IN (", ref_junID_str, ")")) %>% as_tibble(), by = "ref_junID") %>%
        left_join(dbGetQuery(con, paste0("SELECT * FROM 'transcript' WHERE id IN (", paste(df_introns$transcript_id, collapse = ","), ")")) %>% as_tibble(), by = c("transcript_id" = "id")) %>%
        left_join(dbGetQuery(con, paste0("SELECT * FROM 'gene' WHERE id IN (", paste(df_introns$gene_id, collapse = ","), ")")) %>% as_tibble(), by = c("gene_id" = "id"))
      
      df_introns %>% mutate(tissue = table)
      
    })
    
    # Tidy and subsample
    df_database_introns_tidy <- database_introns %>%
      group_by(tissue) %>%
      distinct(ref_junID, .keep_all = TRUE) %>%
      ungroup() %>%
      mutate(mean_coverage = log10(ref_sum_counts / ref_n_individuals))
    
    saveRDS(df_database_introns_tidy, file = file.path(args$results_folder, "supplementary_figure11_before_subsampling.rds"))
    
    # Subsampling
    data_combined <- df_database_introns_tidy %>%
      group_by(ref_junID) %>%
      filter(n() == 2) %>%
      ungroup() %>%
      mutate(tissue_int = ifelse(tissue == "Skin - Sun Exposed (Lower leg)_SKIN", 0, 1))
    
    message("Start data subsampling...")
    
    set.seed(100)
    m.out <- MatchIt::matchit(tissue_int~mean_coverage, data = data_combined)
    subsample <- MatchIt::match.data(m.out) 
    
    message("Data subsampling finished!")
    
    saveRDS(object = subsample, file = file.path(args$results_folder, "supplementary_figure11_after_subsampling.rds"))
    
    df_database_introns <- df_database_introns_tidy
    subsample_introns <- subsample
    rm(subsample, df_database_introns_tidy)
    gc()
    
  } else {
    df_database_introns <- readRDS(file.path(args$results_folder, "supplementary_figure11_before_subsampling.rds"))
    subsample_introns <- readRDS(file.path(args$results_folder, "supplementary_figure11_after_subsampling.rds"))
  }
  
  # Visualization
  plot_coverage <- function(df, title) {
    ggplot(data = df %>% distinct(ref_junID, .keep_all = TRUE)) +
      geom_density(aes(x = mean_coverage, fill = tissue), alpha = 0.8) +
      ggtitle(title) +
      xlab("log10 mean expression level") +
      ggsci::scale_fill_npg() +
      theme_light() +
      custom_ggtheme +
      theme(legend.position = 'top', legend.spacing.x = unit(0.3, 'cm'))
  }
  
  
  df_database_introns_tidy <- df_database_introns %>%
    mutate(tissue = recode(tissue, "Skin - Sun Exposed (Lower leg)" = "Skin - Sun Exposed", "Skin - Not Sun Exposed (Suprapubic)" = "Skin - Not Sun Exposed"))
  
  subsample_introns_tidy <- subsample_introns %>%
    mutate(tissue = recode(tissue, 
                           "Skin - Sun Exposed (Lower leg)" = "Skin - Sun Exposed", 
                           "Skin - Not Sun Exposed (Suprapubic)" = "Skin - Not Sun Exposed"))
  
  df_database_introns_tidy %>% dplyr::count(tissue)
  subsample_introns_tidy %>% dplyr::count(tissue)
  
  plot_coverage_bf <- plot_coverage(df_database_introns_tidy, "Before subsampling")
  plot_coverage_af <- plot_coverage(subsample_introns_tidy, "After subsampling")
  
  ggpubr::ggarrange(plot_coverage_bf, plot_coverage_af, labels = c("a", "b"), common.legend = TRUE)
  
  ggplot2::ggsave(file.path(args$figures_folder, "supplementary_figure11.png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  # Statistics if needed
  if (stats) {
    
    subsample_introns_tidy_paired_MSRD <- subsample_introns_tidy %>% 
      dplyr::select(MSR_D, subclass, tissue) %>%
      spread(key = tissue, MSR_D)
      
    
    donor_test <- wilcox.test(subsample_introns_tidy_paired_MSRD$`Skin - Sun Exposed (Lower leg)_SKIN`,
                              subsample_introns_tidy_paired_MSRD$`Skin - Not Sun Exposed (Suprapubic)_SKIN`,
                              #alternative = "greater",
                              paired = TRUE, 
                              correct = TRUE)
    
    
    subsample_introns_tidy_paired_MSRA <- subsample_introns_tidy %>% 
      dplyr::select(MSR_A, subclass, tissue) %>%
      spread(key = tissue, MSR_A)
    
    acceptor_test <- wilcox.test(subsample_introns_tidy_paired_MSRA$`Skin - Sun Exposed (Lower leg)_SKIN`,
                                 subsample_introns_tidy_paired_MSRA$`Skin - Not Sun Exposed (Suprapubic)_SKIN`,
                                 #alternative = "greater",
                                 paired = TRUE, 
                                 correct = TRUE)
    
    list(donor_test = donor_test, acceptor_test = acceptor_test)
  }
}

## We Found that the expression levels of 107 RBPs (FDR<0.04) and 5 essential NMD genes67 (FDR<0.04) decreased with age in multiple tissues (Supplementary Figure 18a,b) (Supplementary Table 7). 
## Focusing on brain tissue alone, 40% of the 115 RBPs studied had decreased expression levels with age (FDR<0.04) (Supplementary Figure 18b).

supplementary_figure18 <- function() {
  
  
  source("~/PROJECTS/recount3-database-project/scripts/27_age_effect_uncorrected_TPM_lm.R")
  source("~/PROJECTS/splicing-accuracy-manuscript/scripts/99_utils.R")
  
  
  ###################################
  ## LOAD THE GENE LIST
  ###################################
  
  gene_type <- "NMD"
  
  if (gene_type == "NMD") {
    gene_list <- data.frame(id = c("ENSG00000005007", "ENSG00000151461", "ENSG00000169062", "ENSG00000157106", "ENSG00000198952", "ENSG00000070366", "ENSG00000116698"),
                            name = c("UPF1", "UPF2", "UPF3", "SMG1", "SMG5", "SMG6", "SMG7"))
  } else {
    gene_list <- all_RBPs <- xlsx::read.xlsx(file = file.path(args$dependencies_folder, '/RBPs_subgroups.xlsx'), header = TRUE, sheetIndex = 1) %>% as_tibble() %>% distinct(name, .keep_all = T)
  }
  
  
  ###################################
  ## LINEAR REGRESSION TO TEST IF THE
  ## COVARIATE AGE AFFECTS THE LOG10 TPM LEVELS
  ###################################
  
  if ( !file.exists(file.path(args$results_folder, paste0(gene_type, "_genes_age_lm.rds"))) ) {
    
    gene_age_lm <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      
      message(Sys.time(), " - ", project_id)
      local_results_folder <- file.path(base_folder, "results", main_project, "/", gtf_version, project_id)
      
      if ( file.exists(paste0(local_results_folder, "/", gene_type, "_expression/", project_id, "_tpm_ensembl", gtf_version,".rds")) &&
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
        recount_tpm <- readRDS(file = file.path(base_folder, "results/tpm/", paste0(project_id, "_tpm.rds"))) %>%
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
                file = paste0(local_results_folder, "/", gene_type, "_expression/", 
                              project_id, "_tpm_ensembl", gtf_version, ".rds"))
        
        
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
                                                 results.folder = paste0(local_results_folder,  "/", gene_type, "_expression/"))
      
      return(lm_output %>%
               mutate(project = project_id))
      
      
    })
    
    saveRDS(object = gene_age_lm %>% as_tibble(),
            file = file.path(args$results_folder,"/", paste0(gene_type, "_genes_age_lm.rds")) ) 
    
    write_csv(x = gene_age_lm,
              file = paste0(args$results_folder,"/", paste0(gene_type, "_genes_age_lm.csv")))
    
  } else {
    message("Loading '", paste0(gene_type, "_genes_age_lm.rds"), "' file ...")
    gene_age_lm <- readRDS(file = file.path(args$results_folder, paste0(gene_type, "_genes_age_lm.rds")) ) 
  }
  
  tissues_data_tidy <- gene_age_lm %>%  filter(covariate == "gtex.age") 
  
  
  ##################################
  ## GET STATS FOR THE PAPER
  ##################################
  
  ## We formally assessed this in the GTEx dataset, and found that the expression levels of 107 RBPs (FDR<0.04) and 5 essential NMD genes67 (FDR<0.04) decreased with age in multiple tissues
  
  tissues_data_tidy %>%
    filter(q <= 0.05, Estimate < 0) %>%
    pull(name) %>%
    unique() %>%
    length()
  
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
  
  brain_data_tidy <- tissues_data_tidy %>%
    filter(covariate == "gtex.age",
           project == "BRAIN") %>%
    distinct(name, .keep_all =T)
  
  
  ((brain_data_tidy %>%
      filter(q <= 0.05, Estimate < 0) %>%
      nrow) * 100) / (brain_data_tidy %>% nrow)
  
  brain_data_tidy$name
  
  brain_data_tidy %>%
    filter(q <= 0.05, Estimate < 0) %>%
    pull(q) %>%
    summary()
  
  ##################################
  ## PLOT
  ##################################
  
  ## Check the NMD factors with TPM values affected by age
  gene_age_lm_tidy <- tissues_data_tidy %>%
    filter(covariate == "gtex.age") %>%
    mutate(q = ifelse(Estimate < 0, q, NA)) %>%
    mutate(q = ifelse(q > 0.05, NA, q)) %>%
    mutate(`log10(q)` = q %>% log10()) %>%
    group_by(name) %>%
    mutate(name = factor(name, levels=(gene_age_lm$name)[order(gene_age_lm$name %>% unique %>% dplyr::desc())] %>% unique))  %>%
    ungroup()
  
  #gene_age_lm_tidy$name <- factor(gene_age_lm_tidy$name, levels=(gene_age_lm_tidy$name)[order(gene_age_lm_tidy$Estimate)])
  
  ggplot(data = gene_age_lm_tidy %>%
           mutate(project = str_replace(project, pattern = "_", replacement = " ")) ) + 
    geom_tile(mapping = aes(x = fct_rev(name), 
                            y = fct_rev(project), 
                            fill = `log10(q)`, 
                            colour = "q>=0.05")) +
    scale_fill_gradient(low = "red", 
                        high = "white", 
                        na.value = '#cccccc') +
    scale_colour_manual(values = c( "q>=0.05" = "#cccccc")) +
    xlab("") + 
    ylab("") + 
    theme_light()  +
    custom_ggtheme +
    guides(colour = guide_legend(override.aes = list(fill = '#cccccc'),
                                 title = "",
                                 label.position = "bottom",
                                 order = 2)) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 5),
          axis.text.y = element_text(size = 5, colour = "black"),
          legend.position = "top",
          legend.text = element_text(size = 5, colour = "black"),
          legend.title = element_text(size = 5, colour = "black"),
          legend.box.margin = margin(l = -10, r = -10, b = -10, t = -5
          ))
  
  
  
  file_name <- if (gene_type == "NMD") {
    file.path(args$figures_folder, "supplementary_figure18_a.png")
  } else {
    file.path(args$figures_folder, "supplementary_figure18_b.png")
  }
  ggplot2::ggsave(filename = file_name, width = 180, height = 90, units = "mm", dpi = 300)
  
  
}

supplementary_figure28 <- function() {
  
  
  #############################################
  ## PREPARE DATA BEFORE PLOT
  #############################################
  
  ## Join both datasets of junction lengths
  df_all_lengths_tidy <- rbind(master_introns %>%
                                 dplyr::select(length = ref_length) %>%
                                 mutate(type = "annotated intron"),
                               master_novel_junctions %>%
                                 dplyr::select(length = novel_length) %>%
                                 mutate(type = "novel junction"))
  df_all_lengths_tidy %>% head()
  
  
  mode_length_annotated <- get_mode(df_all_lengths_tidy %>%
                                      filter(type == "annotated intron") %>%
                                      pull(length))
  mode_length_novel <- get_mode(df_all_lengths_tidy %>%
                                  filter(type == "novel junction") %>%
                                  pull(length))
  
  #############################################
  ## PLOT IMPLIED INTRON LENGTHS
  #############################################
  
  df_all_lengths_tidy$type = factor( df_all_lengths_tidy$type, 
                                     levels = c( "novel junction", "annotated intron" ) )
  
  
  ggplot(data = df_all_lengths_tidy %>% filter(length <= 200)) + 
    geom_density(aes(x = length, colour = type), 
                 alpha = 0.7, stat = "count", position = "identity") +
    geom_vline(xintercept = c(mode_length_annotated, mode_length_novel), 
               colour = c("#661100", "#009E73"),
               linetype="dotted", 
               linewidth = 1) +
    annotate(x = mode_length_annotated, y = 1000, label = mode_length_annotated %>% as.character(), vjust = 2, geom = "label") +
    annotate(x = mode_length_novel, y = 1000, label = mode_length_novel %>% as.character(), vjust = 2, geom = "label") +
    xlab("Implied intron length (in bp)") +
    ylab("Number of unique junctions") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12")) +
    guides(colour = guide_legend(title = "Split Reads Category: ",
                                 ncol = 2,
                                 nrow = 1)) +
    scale_colour_manual(values = c("#661100", "#009E73"), 
                        name = "Category",
                        breaks = c("annotated intron", "novel junction"),
                        labels = c("Annotated intron", "Novel junction")) +
    custom_ggtheme
  
  ## Save plot
  ggplot2::ggsave(filename = file.path(args$figures_folder, "supplementary_figure28.png"), width = 180, height = 90, units = "mm", dpi = 300)
  
}


##################################
## SUPPLEMENTARY TABLES
##################################

supplementary_table1 <- function() {

  
  SRA_projects <- df_metadata$SRA_project %>% unique()
  
  df_read_count <- map_df(SRA_projects, function(project_id) {
    
    # project_id <- SRA_projects[1]
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    
    map_df(all_clusters, function(cluster_id) {
      
      # cluster_id <- all_clusters[1]
      
      message(Sys.time(), " - ", cluster_id)
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts 
                      FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      db_never <- dbGetQuery(con, query) %>% as_tibble()
      db_never$ref_sum_counts
      
      query <- paste0("SELECT ref_junID, novel_junID, ref_sum_counts, novel_sum_counts 
                      FROM '", cluster_id, "_", project_id, "_misspliced'")
      db_misspliced <- dbGetQuery(con, query) %>% as_tibble()
      db_misspliced$ref_sum_counts
      
      novel_donor <- db_misspliced  %>%
        filter(novel_junID %in% (master_novel_junctions %>% 
                                   filter(novel_type == "novel_donor") %>%
                                   pull(novel_junID)))
      
      novel_acceptor <- db_misspliced  %>%
        filter(novel_junID %in% (master_novel_junctions %>% 
                                   filter(novel_type == "novel_acceptor") %>%
                                   pull(novel_junID)))
      
      
      return( data.frame(annotated_unique_junctions = c(db_never$ref_junID,
                                                        db_misspliced$ref_junID) %>% unique() %>% length(),
                         annotated_median_read_count = c(db_never$ref_sum_counts,
                                                         db_misspliced$ref_sum_counts) %>% median(),
                         annotated_max_read_count = c(db_never$ref_sum_counts,
                                                      db_misspliced$ref_sum_counts) %>% max(),
                         annotated_min_read_count = c(db_never$ref_sum_counts,
                                                      db_misspliced$ref_sum_counts) %>% min(),
                         novel_donor_unique_junctions = novel_donor$novel_junID %>% unique() %>% length(),
                         novel_donor_median_read_count = novel_donor$novel_sum_counts %>% median(),
                         novel_donor_max_read_count = novel_donor$novel_sum_counts %>% max(),
                         novel_donor_min_read_count = novel_donor$novel_sum_counts %>% min(),
                         novel_acceptor_unique_junctions = novel_acceptor$novel_junID %>% unique() %>% length(),
                         novel_acceptor_median_read_count = novel_acceptor$novel_sum_counts %>% median(),
                         novel_acceptor_max_read_count = novel_acceptor$novel_sum_counts %>% max(),
                         novel_acceptor_min_read_count = novel_acceptor$novel_sum_counts %>% min(),
                         tissue = cluster_id) )
      
    })
    
  })
  
  write.csv(x = df_read_count %>% relocate(tissue), file = file.path(args$results_folder, "/supplementary_table1.csv"), row.names = F)  
  
  

  ## Focusing on frontal cortex, we found that this equated to a median read count 
  ## of 2,700 for annotated junctions with novel donor or acceptor events having a 
  ## median read count of only 2 in both cases. 
  
  df_read_count %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)") %>%
    as.data.frame()
  
}

supplementary_table2 <- function() {
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  limit_bp <- 30
  
  df_distances_all_tissues <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    
    map_df(all_clusters, function(cluster_id) {
      
      # cluster_id <- all_clusters[1]
      
      query <- paste0("SELECT tissue.novel_junID, tissue.ref_junID, novel.novel_type, novel.distance 
                      FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue
                      INNER JOIN 'novel' ON novel.novel_junID = tissue.novel_junID
                      WHERE distance >= ", (limit_bp*-1), " AND distance <= ", limit_bp)
      db_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
      
      ## Stats - Donor - Negative
      mode_negative_distance_donor <- db_misspliced_introns %>% 
        filter(distance < 0, novel_type == "novel_donor") %>%
        pull(distance) %>%
        get_mode()
      min_negative_distance_donor <- db_misspliced_introns %>% 
        filter(distance < 0, novel_type == "novel_donor") %>%
        pull(distance) %>%
        min()
      max_negative_distance_donor <- db_misspliced_introns %>% 
        filter(distance < 0, novel_type == "novel_donor") %>%
        pull(distance) %>%
        max()
      
      ## Stats - Donor - Positive
      mode_positive_distance_donor <- db_misspliced_introns %>% 
        filter(distance > 0, novel_type == "novel_donor") %>%
        pull(distance) %>%
        get_mode()
      min_positive_distance_donor <- db_misspliced_introns %>% 
        filter(distance > 0, novel_type == "novel_donor") %>%
        pull(distance) %>%
        min()
      max_positive_distance_donor <- db_misspliced_introns %>% 
        filter(distance > 0, novel_type == "novel_donor") %>%
        pull(distance) %>%
        max()
      
      
      ## Stats - Acceptor - Negative
      mode_negative_distance_acceptor <- db_misspliced_introns %>% 
        filter(distance < 0, novel_type == "novel_acceptor") %>%
        pull(distance) %>%
        get_mode()
      min_negative_distance_acceptor <- db_misspliced_introns %>% 
        filter(distance < 0, novel_type == "novel_acceptor") %>%
        pull(distance) %>%
        min()
      max_negative_distance_acceptor <- db_misspliced_introns %>% 
        filter(distance < 0, novel_type == "novel_acceptor") %>%
        pull(distance) %>%
        max()
      
      ## Stats - Acceptor - Negative
      mode_positive_distance_acceptor <- db_misspliced_introns %>% 
        filter(distance > 0, novel_type == "novel_acceptor") %>%
        pull(distance) %>%
        get_mode()
      min_positive_distance_acceptor <- db_misspliced_introns %>% 
        filter(distance > 0, novel_type == "novel_acceptor") %>%
        pull(distance) %>%
        min()
      max_positive_distance_acceptor <- db_misspliced_introns %>% 
        filter(distance > 0, novel_type == "novel_acceptor") %>%
        pull(distance) %>%
        max()
      
      data.frame(tissue = cluster_id,
                 mode_negative_distance_donor = mode_negative_distance_donor,
                 min_negative_distance_donor = min_negative_distance_donor,
                 max_negative_distance_donor = max_negative_distance_donor,
                 
                 mode_positive_distance_donor = mode_positive_distance_donor,
                 min_positive_distance_donor = min_positive_distance_donor,
                 max_positive_distance_donor = max_positive_distance_donor,
                 
                 mode_negative_distance_acceptor = mode_negative_distance_acceptor,
                 min_negative_distance_acceptor = min_negative_distance_acceptor,
                 max_negative_distance_acceptor = max_negative_distance_acceptor,
                 
                 mode_positive_distance_acceptor = mode_positive_distance_acceptor,
                 min_positive_distance_acceptor = min_positive_distance_acceptor,
                 max_positive_distance_acceptor = max_positive_distance_acceptor) %>%
        return()
      
      
    })
    
  })
    
  write.csv(x = df_distances_all_tissues, file = file.path(args$results_folder, "supplementary_table2.csv"), row.names = F)
}

supplementary_table3 <- function()  {
  
  ############################
  ## LOAD/GENERATE DATA
  ############################
  
  file_name <- file.path(args$results_folder, "supplementary_table3.csv")
  
  print(paste0(Sys.time(), " - getting MSR across tissues..."))
  
  
  MSR_all_tissues <- map_df(all_projects, function(project_id) {
    
    #############################
    ## GET THE CLUSTERS
    #############################
    
    print(paste0(Sys.time(), " - ", project_id))
    
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    
    map_df(all_clusters, function(cluster_id) {
      
      print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))
      
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                      FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                      FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      introns <- introns %>%
        left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding) %>% as_tibble(),
                  by = "ref_junID") %>% 
        as_tibble() 
      
      print(paste0(Sys.time(), " - ", introns %>% nrow(), " - introns!"))
      
      ##############################################################################
      ## CORRECT BY MEAN COVERAGE
      ##############################################################################
      
      
      ## TIDY DATAFRAME
      df_biotype_result_tidy <- introns %>%
        filter(protein_coding %in% c(0,100)) %>%
        mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC"))  %>%
        dplyr::rename(biotype = type_PC) %>%
        distinct(ref_junID, .keep_all = T) %>%
        group_by(ref_junID) %>%
        mutate(mean_expression = (ref_sum_counts/ref_n_individuals) %>% log10()) %>%
        ungroup()
      
      
      #### START SUBSAMPLIING
      
      print(paste0(Sys.time(), " - starting subsampling process..."))
      set.seed(100)
      ## Subsampling introns to control by similarity in mean read coverage
      m.out <- MatchIt::matchit(biotype ~ mean_expression, 
                                data = df_biotype_result_tidy, 
                                distance = df_biotype_result_tidy$mean_expression,
                                method = "nearest", 
                                caliper = c(mean_expression = 0.005), 
                                std.caliper = FALSE)
      subsample2 <- MatchIt::match.data(m.out)
      df_introns_biotype <- subsample2
      
      # subsample2 %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)
      
      
      any(df_introns_biotype %>% pull(ref_junID) %>% duplicated())
      
      df_introns_biotype$biotype = factor(df_introns_biotype$biotype, levels = c("PC", "non PC"))
      df_introns_biotype <- df_introns_biotype %>%
        dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
        gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)
      
      
      ###################################################
      ## TEST 1 - DONOR IS LESS MISSPLICED THAN ACCEPTOR
      ###################################################
      
      print(paste0(Sys.time(), " - starting wilcoxon test #1..."))
      
      ## Introns are less miss-spliced at the donor than at the acceptor
      wilcox_MSR_test <- wilcox.test(x = df_introns_biotype %>% filter(MSR_type == "MSR_D") %>% pull(MSR),
                                     y = df_introns_biotype %>% filter(MSR_type == "MSR_A") %>% pull(MSR),
                                     alternative = "less",
                                     paired = T,
                                     correct = T)
      effsize_MSR_test <- rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                                                    filter(biotype == "PC") %>%
                                                    mutate(MSR_type = MSR_type %>% as.factor()),
                                                  formula = MSR ~ MSR_type,
                                                  paired = T)
      
      
      ##############################################
      ## TEST 2 - PC IS LESS MISSPLICED THAN NON-PC
      ##############################################
      
      print(paste0(Sys.time(), " - starting wilcoxon test #2..."))
      
      ## Introns from Non-PC transcripts are more mis-spliced at the donor than protein-coding introns
      wilcox_D <- wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_D") %>% pull(MSR),
                              y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_D") %>% pull(MSR),
                              alternative = "less",
                              paired = T,
                              correct = T)
      effsize_D <- rstatix::wilcox_effsize(data = df_introns_biotype %>%  
                                             filter(MSR_type == "MSR_D") %>%
                                             mutate(biotype = biotype %>% as.factor()),
                                           formula = MSR ~ biotype,
                                           paired = T)
      
      
      
      print(paste0(Sys.time(), " - starting wilcoxon test #3..."))
      
      ## Introns from Non-PC transcripts are more mis-spliced at the acceptor than protein-coding introns
      wilcox_A <- wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_A") %>% pull(MSR),
                              y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_A") %>% pull(MSR),
                              alternative = "less",
                              paired = T,
                              correct = T)
      effsize_A <- rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                                             filter(MSR_type == "MSR_A")%>%
                                             mutate(biotype = biotype %>% as.factor()),
                                           formula = MSR ~ biotype,
                                           paired = T)
      
      
      ##############################################
      ## SAVE RESULTS
      ##############################################
      
      print(paste0(Sys.time(), " - saving results..."))
      
      data.frame(test = c("one-sided paired Wilcoxon signed rank test with continuity correction #1", 
                          "one-sided paired Wilcoxon signed rank test with continuity correction #2",
                          "one-sided paired Wilcoxon signed rank test with continuity correction #3"),
                 H1 = c("Introns are overall less miss-spliced at the donor than at the acceptor splice site.", 
                        "Introns from Non-PC transcripts are overall more mis-spliced at the donor than introns from PC transcripts.",
                        "Introns from Non-PC transcripts are overall more mis-spliced at the acceptor than introns from PC transcripts."),
                 pvalue = c(wilcox_MSR_test$p.value,
                            wilcox_D$p.value,
                            wilcox_A$p.value),
                 effsize = c(effsize_MSR_test$effsize %>% unname,
                             effsize_D$effsize %>% unname,
                             effsize_A$effsize %>% unname))  %>%
        mutate(tissue = cluster_id) %>%
        return()
      
    })
    
  })
  
  write.csv(x = MSR_all_tissues, file = file_name, row.names = F)
  
  
  
}



##################################
## AUX FUNCTIONS
##################################

# Helper function for querying and merging
run_sql_query <- function(con, table, columns, table_suffix) {
  query <- paste0("SELECT DISTINCT ", columns, " FROM '", table, table_suffix, "'")
  dbGetQuery(con, query) %>% as_tibble()
}

get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}

##################################
## CALLS
##################################

#supplementary_figure10()
get_data_main_figure4_d(common_introns=F, replace=T)