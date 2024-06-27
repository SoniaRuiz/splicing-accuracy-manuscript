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

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/splicing_accuracy_manuscript_figures.R")

base_folder <- here::here()

## CONNECT TO THE DATABASE ------------------------------
supporting_reads <- 1
gtf_version <- 111
main_project <- paste0("GTEX_",supporting_reads,"read_subsampleFALSE")

database_path <- paste0(base_folder, "/database/", main_project, "/", gtf_version, "/", main_project, ".sqlite")

con <- dbConnect(RSQLite::SQLite(), database_path)
tables <- dbListTables(con)


## SET PATHS TO FOLDERS

dependencies_folder <- paste0(base_folder, "/dependencies/")

results_folder <- file.path(base_folder, "results", main_project, gtf_version, "_paper_review/results")
dir.create(file.path(results_folder), recursive = TRUE, showWarnings = F)

figures_folder <- file.path(base_folder, "results", main_project, gtf_version, "_paper_review/figures")
dir.create(file.path(figures_folder), recursive = TRUE, showWarnings = F)


## QUERY MASTER TABLES 

query = paste0("SELECT * FROM 'metadata'")
df_metadata <- dbGetQuery(con, query) %>% distinct(sample_id, .keep_all = T) %>% as_tibble()
all_projects <- df_metadata$SRA_project %>% unique

query <- paste0("SELECT * FROM 'intron'")
master_introns <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'novel'")
master_novel_junctions <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'transcript'")
master_transcripts <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'gene'")
master_genes <- dbGetQuery(con, query) %>% as_tibble()


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

#query <- paste0("SELECT * from 'Adipose - Subcutaneous_ADIPOSE_TISSUE_nevermisspliced'")
#dbGetQuery(con, query) %>% as_tibble()

########################################
## FUNCTIONS - produce figures for the
## paper
########################################


## SECTION 1 ---------------------------------------------


## 1. Novel donor and acceptor junctions are commonly detected and exceed the number of unique annotated introns by an average of 11-fold

get_database_stats <- function() {
  
  ## We found that 268,988 (82.8%) annotated introns had at least a single associated novel donor or acceptor junction,
  ## with only 55,968 annotated introns appearing to be precisely spliced across all the samples and tissues studied.

  master_introns %>% head()
  master_introns %>% distinct(ref_junID) %>% nrow()
  master_introns %>% dplyr::count(misspliced)
  
  
  ## Collectively, we detected 3,865,268 unique novel junctions, equating to 14 novel junctions per annotated intron. 
  master_novel_junctions %>% head
  master_novel_junctions %>% nrow() 
  master_novel_junctions %>% dplyr::count(novel_type)
  
  (master_novel_junctions %>% nrow()) / (master_introns %>% filter(misspliced=="Yes") %>% nrow())
  
  ## Collectively, we detected 31,811 genes and 199,551 transcripts
  master_transcripts %>% nrow()
  
  ## JOIN WITH GENE DATA
  master_genes %>% distinct(gene_id) %>% nrow()
  
  ## Novel junctions exceed in X fold to annotated introns
  (master_novel_junctions %>% distinct(novel_junID) %>% nrow()) / (master_introns %>% distinct(ref_junID) %>% nrow())
  
  
  ## Percentage of mis-spliced introns
  (((master_introns %>%
      dplyr::count(misspliced) %>%
      filter(misspliced == "Yes") %>%
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
    inner_join(y = df_novel_jxn_count,
               by = "cluster") %>%
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
  
  
  
  tables <- dbListTables(con)
  
  ## Methods: number of samples by RIN number
  df_metadata %>% nrow() 
  
  df_metadata %>%
    filter(rin >= 8) %>% nrow()
  df_metadata %>%
    filter(rin >= 7) %>% nrow()
  df_metadata %>%
    filter(rin >= 6) %>% nrow()
  
  if ( df_metadata %>% filter(rin < 6) %>% nrow() > 1 ) {
    print("ERROR! Some of the samples considered present a RIN number lower than 6!")
    break;
  }
  
  
  df_metadata %>%
    dplyr::count(cluster) %>%
    print(n = 50)
  
}

get_junc_length <- function() {
  
  
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
  

  
  
  
  ggplot(data = df_all_lengths_tidy %>%
           filter(length <= 200)) + 
    geom_density(aes(x = length, 
                   colour = type), 
               alpha = 0.7,
               stat = "count",
               position = "identity") +
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
  file_name <- paste0(figures_folder, "/junction_length")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
}

get_intron_data_summary <- function() {
  
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  #######################################
  ## GET DATA FROM THE SELECTED TISSUE
  
  query <- paste0("SELECT * 
                  FROM '", cluster_id, "_", project_id, "_nevermisspliced' " )
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT * 
                  FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
  introns %>% nrow()
  
  introns %>% head()
  
  
  #######################################
  ## ADD DATA FROM MASTER TABLES
  
  fctx_introns <- introns %>%
    inner_join(y = master_introns %>% dplyr::select(-transcript_id),
               by = "ref_junID") %>% 
    as_tibble() 
  fctx_introns %>% nrow()
  
  ## JOIN WITH TRANSCRIPT DATA
  query <- paste0("SELECT * FROM 'transcript'")
  fctx_introns <- fctx_introns %>%
    inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
               by = c("transcript_id" = "id")) %>% 
    as_tibble() 
  fctx_introns %>% nrow()
  
  
  ## JOIN WITH GENE DATA
  query <- paste0("SELECT * FROM 'gene'")
  introns <- introns %>%
    inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
               by = c("gene_id" = "id")) %>% 
    as_tibble() 
  introns %>% nrow()
  
  introns <- introns %>%
    distinct(ref_junID, .keep_all = T)
  
  introns %>% nrow()
  
  ###########################
  
  #########################
  ## TIDY DATA
  #########################
  
  idb <- introns %>%
    dplyr::select(ref_junID,
                  intron_length = ref_length,
                  intron_5ss_score = ref_mes5ss,
                  intron_3ss_score = ref_mes3ss,
                  gene_length = gene_width,
                  gene_tpm = gene_tpm,
                  gene_num_transcripts = n_transcripts,
                  CDTS_5ss = mean_CDTS5ss_35,
                  CDTS_3ss = mean_CDTS3ss_35,
                  protein_coding = protein_coding,
                  
                  mean_phastCons4way_5ss_35 = mean_phastCons4way5ss_35,
                  mean_phastCons4way_3ss_35 = mean_phastCons4way3ss_35,
                  mean_phastCons4way_5ss_70 = mean_phastCons4way5ss_70,
                  mean_phastCons4way_3ss_70 = mean_phastCons4way3ss_70,
                  mean_phastCons4way_5ss_100 = mean_phastCons4way5ss_100,
                  mean_phastCons4way_3ss_100 = mean_phastCons4way3ss_100,
                  
                  mean_phastCons20way_5ss_35 = mean_phastCons20way5ss_35,
                  mean_phastCons20way_3ss_35 = mean_phastCons20way3ss_35,
                  mean_phastCons20way_5ss_70 = mean_phastCons20way5ss_70,
                  mean_phastCons20way_3ss_70 = mean_phastCons20way3ss_70,
                  mean_phastCons20way_5ss_100 = mean_phastCons20way5ss_100,
                  mean_phastCons20way_3ss_100 = mean_phastCons20way3ss_100,
                  
                  MSR_D,
                  MSR_A) 
  
  #########################################################
  
  
  idb$gene_length %>% summary()
  idb$gene_tpm %>% summary()
  idb$protein_coding %>% summary()
  idb$gene_num_transcripts %>% summary()

  idb$intron_length %>% summary()
  idb$intron_5ss_score %>% summary()
  idb$intron_3ss_score %>% summary()
  idb$CDTS_5ss %>% summary()
  idb$CDTS_3ss %>% summary()
  
  
  ## Conservation - phastCons4
  
  idb$mean_phastCons4way_5ss_35 %>% summary()
  idb$mean_phastCons4way_3ss_35 %>% summary()
  
  idb$mean_phastCons4way_5ss_70 %>% summary()
  idb$mean_phastCons4way_3ss_70 %>% summary()
  
  idb$mean_phastCons4way_5ss_100 %>% summary()
  idb$mean_phastCons4way_3ss_100 %>% summary()

  
  ## Conservation - phastCons20
  
  idb$mean_phastCons20way_5ss_35 %>% summary()
  idb$mean_phastCons20way_3ss_35 %>% summary()
  
  idb$mean_phastCons20way_5ss_70 %>% summary()
  idb$mean_phastCons20way_3ss_70 %>% summary()
  
  idb$mean_phastCons20way_5ss_100 %>% summary()
  idb$mean_phastCons20way_3ss_100 %>% summary()

  
  
  ## MSR
  
  idb$MSR_D %>% summary()
  idb$MSR_A %>% summary()

  idb %>%
    filter(MSR_D == 0) %>% nrow
  idb %>%
    filter(MSR_D > 0) %>% nrow
  idb %>%  nrow

  idb %>% nrow()


  plot(idb$MSR_A,
       idb$CDTS_5ss, pch = 16, cex = 1.3, col = "blue")
  abline(glm(MSR_A ~ CDTS_5ss, data = idb))
  #########################################################
  
}

plot_junction_sharing <- function () {
  
 
  df_jxn_sharing_all_tissues <- map_df(all_projects, function(project_id) {
    
    
    # project_id <- all_projects[9]
    print(paste0(Sys.time(), " - ", project_id))
    
    #############################
    ## GET THE CLUSTERS
    #############################
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(all_clusters, function(cluster_id) {
      
      # cluster_id <- all_clusters[1]
      ## Query the database
      
      if ( any(str_detect(tables, cluster_id)) ) {
        
        print(paste0(Sys.time(), " - ", cluster_id))
        
        query <- paste0("SELECT *
                        FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue")
        db_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
        
        query <- paste0("SELECT *
                        FROM '", cluster_id, "_", project_id, "_nevermisspliced' AS tissue")
        db_not_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
        
        db_all <-  plyr::rbind.fill(db_misspliced_introns, 
                                    db_not_misspliced_introns) %>%
          mutate("project_id" = project_id,
                 "cluster_id" = cluster_id)
        
        return(db_all)
      }
      
      
    })
  })
  
  
  
  
  ################################
  ## PLOT NOVEL JUNCTION SHARING
  ################################
  
  
  df_jxn_sharing_all_tissues_novel <- df_jxn_sharing_all_tissues %>%
    group_by(project_id) %>%
    distinct(novel_junID, .keep_all = T) %>%
    ungroup() %>%
    drop_na()
  
  
  ## Explore the data
  df_jxn_sharing_all_tissues_novel %>% 
    filter(project_id == "BRAIN") %>%
    dplyr::count(novel_n_individuals) %>%
    as.data.frame()
  df_jxn_sharing_all_tissues_novel$novel_n_individuals %>% unique %>% summary
  
  
  ## Plot
  ggplot(df_jxn_sharing_all_tissues_novel %>%
           mutate(project_id = str_replace_all(string = project_id,
                                               pattern = "_",
                                               replacement = " ")), 
         aes(x = novel_n_individuals)) +
    geom_histogram(boundary = 1, closed = "left")  +
    scale_x_continuous(breaks = c(1, 200, 400, 600, 800)) +
    facet_wrap(~project_id) +
    #labs(title = 'Novel junction sharing across individuals per GTEx tissue',
    #     caption = "Each novel junctions has at least 2 supporting split reads in at least 2 samples from different tissues") +
    theme(
      legend.position = "none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Number of individuals") +
    ylab("Number of unique novel junctions") +
    theme_light() +
    custom_ggtheme
  
  
  file_name <- paste0(figures_folder, "/jxn_sharing_all_tissues_novel")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 180, units = "mm", dpi = 300)
  
  
  #################################
  ## PLOT ANNOTATED INTRON SHARING
  #################################
  
  df_jxn_sharing_all_tissues_intron <- df_jxn_sharing_all_tissues %>%
    group_by(project_id) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    drop_na()

  
  ggplot(df_jxn_sharing_all_tissues_intron %>%
           mutate(project_id = str_replace_all(string = project_id,
                                               pattern = "_",
                                               replacement = " ")), 
         aes(x = ref_n_individuals)) +
    geom_histogram(boundary = 1, closed = "left")  +
    scale_x_continuous(breaks = c(1, 200, 400, 600, 800)) +
    facet_wrap(~project_id) +
    #labs(title = 'Annotated intron sharing across individuals per GTEx tissue',
    #     caption = "Each represented anotated intron has at least 2 supporting split reads in at least 2 samples from different tissues") +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8) ) +
    xlab("Number of individuals") +
    ylab("Number of unique annotated introns") +
    theme_light() +
    custom_ggtheme
  
  file_name <- paste0(figures_folder, "/jxn_sharing_all_tissues_annotated")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 180, units = "mm", dpi = 300)
  
  
}





## 2. Over 98% of novel donor and acceptor junctions are likely to be generated through splicing errors

get_reclassification_rates_all_tissues <- function () {

  
  #############################
  ## CONNECT TO THE DATABASE
  #############################
 
  contamination_file <- paste0(results_folder, "/contamination_rates.rds")
  
  if ( !file.exists(contamination_file) ) {
    
    
    df_contamination <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      # project_id <- all_projects[12]
      # project_id <-"KIDNEY"
      
      print(paste0(Sys.time(), " - ", project_id))
      
      all_clusters <- df_metadata %>%
        filter(SRA_project == project_id) %>%
        distinct(cluster) %>%
        pull()
      
      versions <- c("97")
      dates <- c("26-May-2019")

      map_df(all_clusters, function(cluster) {
        
        # cluster <- all_clusters[1]
        print(paste0(cluster))
        
        map_df(versions, function(version) {
          
          # version <- versions[1]
          print(paste0("v", version))
          
          if (file.exists(paste0(getwd(), "/results/", project_id, "/v", version, "/", main_project, "/base_data/",
                                 project_id, "_", cluster, "_all_split_reads.rds"))) {
            
            data_old <- paste0(getwd(), "/results/", project_id, "/v", version, "/", main_project, "/base_data/",
                               project_id, "_", cluster, "_all_split_reads.rds")
            
            ## ENSEMBL v97
            df_old <-  readRDS(file = data_old) %>% as_tibble()
            
            db_introns_old <- df_old %>% 
              filter(type == "annotated")
            
            db_novel_old <- df_old %>%
              filter(type %in% c("novel_donor", "novel_acceptor"))
            
            
            ## ENSEMBL v105
            df_new <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v105/", main_project, "/base_data/",
                                            project_id, "_", cluster, "_all_split_reads.rds")) %>% as_tibble()
            db_introns_new <- df_new %>%
              filter(type == "annotated")
            
            db_novel_new <- df_new %>%
              filter(type %in% c("novel_donor", "novel_acceptor"))
            
            
            rm(df_old)
            rm(df_new)
            
            
            ## NOVEL JUNCTIONS THAT KEPT ANNOTATION CATEGORY
            kept_annotation <- db_novel_old %>%
              filter(junID %in% db_novel_new$junID) %>% 
              distinct(junID, .keep_all = T)
            kept_annotation %>% nrow() %>% print()
            
            ## NOVEL JUNCTIONS THAT ENTERED ANNOTATION
            in_annotation <- db_novel_old %>%
              filter(junID %in% db_introns_new$junID) %>% 
              distinct(junID, .keep_all = T)
            in_annotation %>% nrow() %>% print()
            
            ## INTRONS THAT EXITED ANNOTATION
            # out_annotation <- db_introns_old %>%
            #   filter( !(junID %in% db_introns_new$junID) ) %>% 
            #   distinct(junID, .keep_all = T) 
            # out_annotation %>% nrow() %>% print()
            
            # keep_annotation <- db_novel_new %>%
            #   filter(junID %in% db_novel_old$junID) %>% 
            #   distinct(junID, .keep_all = T) 
            # keep_annotation %>% nrow() %>% print()
        
            
            label <- paste0("Ensembl_v", version, " (", dates, ")")
            
            
            # (in_annotation %>% nrow()) / (in_annotation %>% nrow() +
            #                                       keep_annotation %>% nrow()) * 100
            
            return(data.frame(tissue = cluster,
                              contamination_rates = label,
                              kept_annotation = kept_annotation %>% nrow() / db_novel_old %>% nrow(),
                              in_annotation = in_annotation %>% nrow() / db_novel_old %>% nrow() ))
          } else {
            return(NULL)
          }
          
        })
      })
    })
    
    
    contamination_path <- paste0(getwd(), "/results/_paper/results/")
    dir.create(file.path(contamination_path), recursive = TRUE, showWarnings = T)
    saveRDS(object = df_contamination, file = contamination_file)
    
    
  } else {
    df_contamination <- readRDS(file = contamination_file)
  }
  
  df_contamination <- df_contamination %>%
    mutate(in_annotation = in_annotation / 100,
           kept_annotation = kept_annotation / 100)
  
  
  ## all tissues
  df_contamination %>%
    pull(in_annotation) %>% 
    min
  df_contamination %>%
    pull(in_annotation) %>% 
    max
  
  ## Only brain
  df_contamination %>%
    filter(str_detect(tissue, pattern = "Brain")) %>%
    pull(in_annotation) %>% 
    mean
  df_contamination %>%
    filter(str_detect(tissue, pattern = "Brain")) %>%
    pull(in_annotation) %>% 
    max
  df_contamination %>%
    filter(str_detect(tissue, pattern = "Brain")) %>%
    pull(in_annotation) %>% 
    min
  
  ## Prepare the dataframe before plotting it
  df_contamination_tidy <- df_contamination %>% 
    dplyr::select(kept_annotation, in_annotation, tissue) %>%
    tidyr::gather(key = "type", value = "prop", -tissue ) 
  
  df_contamination_tidy$type = factor( df_contamination_tidy$type, 
                                       levels = c( "kept_annotation", "in_annotation" ) )
  
  df_contamination_tidy = df_contamination_tidy %>% 
    ungroup() %>%
    arrange(type , prop) %>%
    mutate(tissue = fct_inorder(tissue))
  
  
  # colours <- ifelse(str_detect(string = as.factor(df_contamination_tidy$tissue), pattern = "Brain"), "red", "black")
  
  ###############################
  ## CONTAMINATION RATES 
  ## BAR PLOT
  ###############################
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data = df_contamination_tidy %>%
           filter(!(tissue %in% c("Brain - Cortex", "Brain - Cerebellum")))) +
    geom_bar(mapping = aes(x = tissue, y = prop, fill = type), stat = "identity", position = "identity") + 
    #coord_flip() +
    ggforce::facet_zoom(ylim = c(0,0.013), zoom.size = 3) +
    ylab("re-classification rates") +
    xlab("Tissue") +
    scale_fill_manual(values = c( "#a6a6a6", "#1a1a1a" ),
                      breaks = c( "kept_annotation", "in_annotation" ),
                      labels = c( "novel junctions mantaining category", "novel junctions entering annotation" )) +
    #scale_fill_viridis_d(option = "A") +
    theme_light() +
    custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) %>%
    return()
  
  figures_path <- paste0(getwd(), "/results/_paper/figures/")
  dir.create(file.path(figures_path), recursive = TRUE, showWarnings = T)
  
  
  file_name <- paste0(figures_path, "/contamination_rates_all_tissues")
  
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  file_name <- paste0(figures_path, "/contamination_rates_all_tissues_large")
  
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 130, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 130, units = "mm", dpi = 300)
  
  
  
  
}

get_reclassification_rates_all_tissues_stats <- function() {
  
  
  contamination_file <- paste0(getwd(), "/results/_paper/results/contamination_rates.rds")
  df_contamination <- readRDS(file = contamination_file)
  
  df_contamination_tidy <- df_contamination %>% 
    dplyr::select(kept_annotation, in_annotation, tissue) %>%
    tidyr::gather(key = "type", value = "prop", -tissue ) 
  
  df_contamination_tidy$type = factor( df_contamination_tidy$type, 
                                       levels = c( "kept_annotation", "in_annotation" ) )
  
  #########
  ## Stats
  #########
  
  
  # We found that across all tissues, on average only 0.8% (1.22-0.54 %) 
  ## of junctions defined as novel donor or acceptor junctions using Ensembl v97 
  ## were reclassified as annotated junctions in v105
  
  df_contamination_tidy %>%
    filter(type == "in_annotation") %>% 
    mutate(average_prop = prop %>% mean()) %>%
    arrange(desc(prop))
  
  
  ## Furthermore, these findings strongly suggest that the vast majority of novel donor 
  ## and acceptor junctions are generated through mis-splicing events, with on average only 
  ## < X% being explained by contamination with junctions originating from stable transcripts.
  
  df_contamination_tidy %>%
    filter(type == "in_annotation") %>% 
    mutate(average_prop = prop %>% mean())
  
  
  rm(df_contamination)
  gc()
  
}

get_reclassification_rates_FCTX <- function () {
  
  #############################
  ## CONNECT TO THE DATABASE
  #############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  
  project_id <- "BRAIN"
  cluster <- "Brain - Frontal Cortex (BA9)"
  reference_gtf_version <- 105
  gtf_versions <- c("76", "81", "90", "97", "104")
  dates <- c("2014", "2015", "2017", "2019", "2021")
  
  ## Load saved file
  if ( !file.exists( paste0(getwd(), "/results/", project_id, "/v", reference_gtf_version, "/", main_project,
                            "/", cluster, "_contamination_rates.rds") ) ) {
    
    df_contamination <- map_df(gtf_versions, function(gtf_version) {
          
      # version <- gtf_versions[1]
      print(paste0("v", gtf_version))
      
      if ( file.exists(paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/",
                              project_id, "_", cluster, "_all_split_reads.rds")) ) {
        
        data_old <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", main_project, "/base_data/",
                           project_id, "_", cluster, "_all_split_reads.rds")
        
        ## ENSEMBL old
        df_old <-  readRDS(file = data_old) %>%
          as_tibble()
        db_introns_old <- df_old %>% 
          filter(type == "annotated")
        db_novel_old <- df_old %>%
          filter(type %in% c("novel_donor", "novel_acceptor"))
        
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
                          contamination_rates = label,
                          in_annotation = (in_annotation %>% nrow() ) / db_novel_old %>% nrow(),
                          out_annotation = (out_annotation %>% nrow() ) / db_introns_old %>% nrow()))
      } else {
        return(NULL)
      }
          

    })
    
    
    file_name <- paste0(getwd(), "/results/", project_id, "/v", reference_gtf_version, "/", main_project,
                        "/", cluster, "_contamination_rates.rds")
    saveRDS(object = df_contamination, file = file_name)
    
    
  } else {
    df_contamination <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", reference_gtf_version, "/", main_project,
                                              "/", cluster, "_contamination_rates.rds"))
  }
  
  
  ## Tidy the dataframe prior display
  df_contamination$contamination_rates = factor(df_contamination$contamination_rates, 
                                                levels = c("v76\n(2014)",
                                                           "v81\n(2015)",
                                                           "v90\n(2017)",
                                                           "v97\n(2019)",
                                                           "v104\n(2021)"))


  ################################
  ## LINES PLOT
  ################################
  df_contamination <- df_contamination %>%
    mutate(in_annotation = in_annotation / 100,
           out_annotation = out_annotation / 100)
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data=df_contamination, aes(x = contamination_rates, 
                                    y = in_annotation,
                                    group=1)) +
    geom_line(linetype = "dashed", col = "red") +
    geom_point() +
    xlab(NULL) +
    ylab("re-classification rates\nas compared to Ensembl v105") +
    theme_light() +
    custom_ggtheme +
    guides(fill = "none")
  
  
  file_name <- paste0(getwd(), "/results/_paper/figures/contamination_rates_FCTX_lineplot")
  ggplot2::ggsave(paste0(file_name, ".svg"),  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"),  width = 180, height = 70, units = "mm", dpi = 300)

}







## 3. Mis-splicing is more common at acceptor than donor splice sites

get_unique_donor_acceptor_jxn <- function() {
  
  ##################################################
  ## Proportion of unique junctions per GTEx tissue
  ##################################################
 
  
  file_name <- file.path(results_folder, "/unique_donor_acceptor_jxn.rds")

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
    write.csv(x = df_proportions, file = paste0(results_folder, "unique_donor_acceptor_jxn.csv"), row.names = F)
  } else {
    df_proportions <- read.csv(file = paste0(results_folder, "unique_donor_acceptor_jxn.csv"))

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
  
  
  ######################
  ## VIOLIN PLOT
  ######################
  
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
  
  reads_violin_plot <- get_unique_donor_acceptor_reads()
  
  
  ggpubr::ggarrange(junx_violin_plot,
                    reads_violin_plot,
                    common.legend = T,
                    labels = c("b", "c"),
                    align = "h",
                    ncol = 2,
                    nrow = 1)
  
  file_name <- paste0(figures_folder, "/panel2bc")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
  # ## Save the figure 3
  # file_name <- paste0(getwd(), "/results/_paper/figures/percent_unique_donor_acceptor_violin")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), 
  #                 width = 183, height = 130, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), 
  #                 width = 90, height = 90, units = "mm", dpi = 300)
}

get_unique_donor_acceptor_reads <- function() {
  
  file_name <- paste0(results_folder, "/unique_donor_acceptor_reads.rds")
  
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
    
    df_mean_counts <- readRDS(file = file_name) %>%
      as_tibble()
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
  
  ######################
  ## VIOLIN PLOT
  ######################
  
  df_mean_counts_violin <- df_mean_counts %>%
    filter(type != "annotated") 
  
  df_mean_counts_violin$type = factor(df_mean_counts_violin$type, 
                                      levels = c("donor", "acceptor"))
  
  
  reads_violin_plot <- ggplot(, aes(type, prop_counts, fill = type)) + 
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

get_mean_read_count_all_tissues <- function() {
  
  #tables <- c("Brain - Frontal Cortex (BA9)_BRAIN_misspliced",
  #            "Brain - Frontal Cortex (BA9)_BRAIN_nevermisspliced")
  
  SRA_projects <- df_metadata$SRA_project %>% unique()
  
 
  if (!file.exists(paste0(results_folder, "/05_supplementary_table1.csv"))) {
    
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
    write.csv(x = df_read_count %>% relocate(tissue),
              file = paste0(results_folder, "/05_supplementary_data_read_count.csv"),
              row.names = F)  
  } else { 
    df_read_count <- read_delim(file = paste0(results_folder, "/05_supplementary_data_read_count.csv"))
  }
  
  
  
  
  

 
  
  ## Focusing on frontal cortex, we found that this equated to a median read count 
  ## of 2,700 for annotated junctions with novel donor or acceptor events having a 
  ## median read count of only 2 in both cases. 
  
  df_read_count %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)") %>%
    as.data.frame()
  
}





## 4. High sequence similarity between novel splice sites and their annotated pairs explains splicing errors

get_maxentscan_score <- function() {
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  file_name <- paste0(results_folder, "/df_mes.rds")
  
  if ( !file.exists(file_name) )  {
    
    ###########################
    ## GET THE MES
    ###########################
    
    ## Get the MES scores corresponding to the mis-spliced MASTER introns

    introns <- master_introns %>%
      filter(misspliced == "Yes") %>%
      dplyr::select(ref_junID, ref_mes5ss, ref_mes3ss)
    
    
    ## Get the MASTER novel junctions
    
    novel_junctions <- master_novel_junctions %>% 
      dplyr::select(ref_junID, novel_junID, novel_type, novel_mes5ss, novel_mes3ss) %>% 
      as_tibble() 
    
    
    ###########################
    ## MERGE THE DATA AND RETURN
    ###########################
    
    df_mes <- novel_junctions %>% 
      left_join(y = introns,
                by = "ref_junID")  %>%
      mutate(diff_ss5score = ref_mes5ss - novel_mes5ss,
             diff_ss3score = ref_mes3ss - novel_mes3ss)
    
    saveRDS(object = df_mes %>% as_tibble(), file = file_name)
    
  } else {
    print("Loading MES results...")
    df_mes <- readRDS(file = file_name) %>% as_tibble()
  }
  
  
  df_mes <- df_mes %>%
    distinct(ref_junID, novel_junID, .keep_all = T)
  
  
  
  #############################
  ## COMBINED FIGURES 
  ## MAXENTSCAN score & DISTANCES
  #############################
  
  
  ## ss5score ----------------------------------------------------------
  
    
  df_5ss <- df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron = ref_mes5ss, novel_donor = novel_mes5ss) %>%
    gather(key = "junction_type", value = "ss5score")
  
  ss5plot <- ggplot(df_5ss, aes(ss5score, fill = junction_type)) +
    geom_density(alpha = 0.9) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    xlab("MES Donor (5'ss) score") +
    scale_fill_manual(values = c("#999999", "#35B779FF"), 
                      breaks=c("intron", "novel_donor"),
                      labels=c("Annotated Intron  ", "Novel Donor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = element_blank(),
                               ncol = 2, nrow = 1))
  
  # ss5plot
  # file_name <- paste0(getwd(), "/results/_paper/figures/MES_ss5plot")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 90, height = 90, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 90, height = 90, units = "mm", dpi = 300)
  # 
  
  ## ss3score -----------------------------------------------------------
  
  df_3ss <- df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron = ref_mes3ss, novel_acceptor = novel_mes3ss) %>%
    gather(key = "junction_type", value = "ss3score")
  
  ss3plot <- ggplot(df_3ss, aes(ss3score, fill = junction_type)) +
    geom_density(alpha = 0.9) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    scale_fill_manual(values = c("#999999", "#64037d"), 
                      breaks = c("intron", "novel_acceptor"),
                      labels = c("Annotated Intron  ", "Novel Acceptor")) +
    xlab("MES Acceptor (3'ss) score") +
    custom_ggtheme +
    guides(fill = guide_legend(title = element_blank(), ncol = 2, nrow = 1))
  
  
  # ss3plot
  # file_name <- paste0(getwd(), "/results/_paper/figures/MES_ss3plot")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 90, height = 90, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 90, height = 90, units = "mm", dpi = 300)
  # 
  
  ## COMBO
  ggpubr::ggarrange( ss5plot, 
                     ss3plot, 
                     labels = c("a", "b"),
                     ncol = 2, 
                     nrow = 1)
  
  file_name <- paste0(figures_folder, "/supplementary_fig3")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 80, units = "mm", dpi = 300)
  
  
  ## Delta MES -----------------------------------------------------------
  
  df_delta_5ss <- df_mes %>%
    filter(novel_type == "novel_donor") %>% 
    dplyr::select(diff_ss5score) %>%
    tidyr::gather(key = "type", value = "MES") 
  
  ## Plot
  deltaplot5ss <- ggplot(data = df_delta_5ss) +
    geom_density(mapping = aes(x = MES, fill = type)) +
    geom_vline(xintercept = 0) +
    xlab("Delta MES") +
    xlim(c(-40, 65)) +
    ylim(c(0, 0.11)) +
    theme_light() +
    scale_color_viridis_d() +
    scale_fill_manual(values =  c("#35B779FF"),
                      breaks = c("diff_ss5score"),
                      labels = c("Delta MES Donor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  # deltaplot5ss
  # file_name <- paste0(getwd(), "/results/_paper/figures/MES_delta_5ss")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 90, height = 90, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 90, height = 90, units = "mm", dpi = 300)
  # 
 
  df_delta_3ss <- df_mes %>%
    filter(novel_type == "novel_acceptor") %>% 
    dplyr::select(diff_ss3score) %>%
    tidyr::gather(key = "type", value = "MES") 
  
  deltaplot3ss <- ggplot(data = df_delta_3ss) +
    geom_density(mapping = aes(x = MES, fill = type)) +
    geom_vline(xintercept = 0) +
    xlab("Delta MES") +
    xlim(c(-40, 65)) +
    ylim(c(0, 0.11)) +
    theme_light() +
    scale_color_viridis_d() +
    scale_fill_manual(values =  c("#64037d"),
                      breaks = c("diff_ss3score"),
                      labels = c("Delta MES Acceptor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  # deltaplot3ss
  
  # file_name <- paste0(getwd(), "/results/_paper/figures/MES_delta_3ss")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 90, height = 90, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 90, height = 90, units = "mm", dpi = 300)
  # 
  ##############
  ## COMBO
  ##############
  
  ggpubr::ggarrange( deltaplot5ss, 
                     deltaplot3ss, 
                     labels = c("a", "b"),
                     ncol = 2, 
                     nrow = 1)
  
  file_name <- paste0(figures_folder, "/panel3ab")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 50, units = "mm", dpi = 300)
  file_name <- paste0(figures_folder, "/panel3ab_large")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
  ##########################
  ## STATS FOR THE PAPER
  ##########################
  
  ## As would be expected, we found that the majority of novel 5’ and 3’ splice sites were 
  ## weaker than the corresponding annotated site with 82.5% of novel 5’ 
 
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
 


  
  ## and 85.2% of novel 3’ sites having positive delta MES scores
  
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
  
  
  
   
  
  ## Interestingly, this analysis demonstrated that novel 5’ splice sites had a modal delta value which was very close to zero
  
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
  
  # ((df_mes %>% 
  #     filter(novel_type == "novel_donor") %>%
  #     dplyr::select(intron_MES = ref_mes5ss, novel_donor_MES = novel_mes5ss) %>%
  #     mutate(MES_diff = intron_MES - novel_donor_MES) %>%
  #     pull(MES_diff) %>% 
  #     sign() %>% 
  #     table() %>% 
  #     as.data.frame() %>%
  #     filter(. == 0) %>%
  #     pull(Freq)) * 100) / ( master_novel_junctions %>% 
  #                       dplyr::count(novel_type) %>%
  #                       filter(novel_type == "novel_donor") %>%
  #                       pull(n) )
  
  
  ## The modal delta value for novel 3’ splice sites was higher 
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
  
  
  # ((df_mes %>% 
  #     filter(novel_type == "novel_acceptor") %>%
  #     dplyr::select(intron_MES = ref_ss5score, novel_acceptor_MES = novel_ss3score) %>%
  #     mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
  #     pull(MES_diff) %>% 
  #     sign() %>% 
  #     table() %>% 
  #     as.data.frame() %>%
  #     filter(. == 0) %>%
  #     pull(Freq)) * 100) / ( db_novel %>% 
  #                              dplyr::count(novel_type) %>%
  #                              filter(novel_type == "novel_acceptor") %>%
  #                              pull(n) )
  
  df_delta_3ss$MES %>% min
  df_delta_3ss$MES %>% max
  
}





## 5. Novel junctions associated with protein-coding transcripts are predicted to be deleterious in 63.5% of cases

get_distances <- function() {
  
  subsample <- F
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  limit_bp <- 30
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  query <- paste0("SELECT tissue.novel_junID, tissue.ref_junID, tissue.ref_sum_counts, tissue.ref_n_individuals, novel.novel_type, novel.distance 
                  FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue 
                  INNER JOIN 'novel' ON novel.novel_junID = tissue.novel_junID")
  db_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  
  #############################
  ## PLOT THE DISTANCES GRAPH 
  #############################
  
  
  df_novel_tidy <- db_misspliced_introns %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  

  plot_all <- ggplot(data = df_novel_tidy) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggforce::facet_col(vars(novel_type)) +
    #ggtitle(paste0("All biotypes\n")) +
    #ylim(y_axes) +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("novel donor  ", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
                               override.aes = list(size = 3),
                               ncol = 2, nrow = 1 )) +
    custom_ggtheme
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15, y = 70),  size = 3, label = "intron") +
    theme_void()
  
  
  plot_all <- plot_all / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  plot_all
  
  
  ## Save plot
  file_name <- paste0(figures_folder, "/FCTX_distances_all")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                  width = 180, height = 100, units = "mm", dpi = 300)
  

  #########################################
  ## GENERATE PLOT FOR ALL ################
  #########################################
  
  ggplot(data = df_novel_tidy) + 
    geom_density(aes(x = distance, fill = novel_type)) +
    ggplot2::facet_grid(vars(novel_type)) +
    #ggtitle(paste0("All biotypes\n")) +
    #ylim(y_axes) +
    xlab("Distance (in bp)") +
    #ylab("Unique novel junctions") +
    theme_light() +
    #scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
    #                   breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("Novel Donor  ", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
                               override.aes = list(size = 3),
                               ncol = 2, nrow = 1 )) +
    custom_ggtheme
  
  file_name <- paste0(figures_folder, "/FCTX_distances_large_all")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
  ###############################
  ## PREPARE DATA BEFORE PLOT
  ###############################

  df_master_introns <- master_introns %>%
    dplyr::select(ref_junID, protein_coding) %>% as_tibble()
  
  df_novel_tidy <- db_misspliced_introns %>%
    inner_join(df_master_introns,
               by = "ref_junID") %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC"),
           novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "),
           novel_type = str_to_title(novel_type)) %>%
    filter(abs(distance) <= limit_bp)
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("Novel Donor", "Novel Acceptor"))
  
  df_novel_tidy$type_PC = factor(df_novel_tidy$type_PC, 
                                 levels = c("protein coding (PC)", "non PC"))
  
  df_novel_tidy %>%
    distinct(novel_junID, .keep_all = T) %>%
    dplyr::count(type_PC)
    
  if ( subsample ) {
    
    df_biotype_result_tidy <- df_novel_tidy %>%
      distinct(ref_junID, .keep_all = T) %>%
      group_by(ref_junID) %>%
      mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
      mutate(mean_coverage = mean_coverage %>% log10()) %>%
      ungroup()
    
    ## Subsampling introns to control by similarity in mean read coverage
    m.out <- MatchIt::matchit(type_PC ~ mean_coverage, 
                              data = df_biotype_result_tidy %>% dplyr::select(-distance), 
                              distance = df_biotype_result_tidy$mean_coverage,
                              method = "nearest", 
                              caliper = c(mean_coverage = 0.005), 
                              std.caliper = FALSE)
    data_subsample <- MatchIt::match.data(m.out)
    data_subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(type_PC)
    
    df_novel_tidy <- df_novel_tidy %>%
      filter(ref_junID %in% data_subsample$ref_junID )
    
  }
  
  
  df_novel_tidy %>%
    dplyr::count(type_PC)
  

  
  
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
  
  plot_legend <- ggpubr::get_legend(ggplot(data = df_novel_tidy %>% 
                                             filter(type_PC == "protein coding (PC)") ) + 
                                      geom_histogram(aes(x = distance, fill = novel_type))+
                                      scale_fill_manual(values = c("#35B779FF","#64037d"),
                                                        breaks = c("Novel Donor", "Novel Acceptor")) + 
                                      custom_ggtheme +
                                      theme(legend.margin = margin(t = 0, r = 5, b = 0, l = 5, unit="cm"), 
                                            legend.box.margin = margin(l = 0, b = -5)))
  
  ###############################
  ## PROTEIN CODING & DONOR
  ###############################
  
  plot_PC_donor <- ggplot(data = df_novel_tidy %>% 
                              filter(type_PC == "protein coding (PC)",
                                     novel_type == "Novel Donor")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack") +
    #ylim(c(0,250))+
    ggtitle("Protein-coding transcripts") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c(limit_bp,(limit_bp * -1)),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6),
                       trans = "reverse") +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  plot_PC_donor
  

  
  plot_PC_donor <- plot_PC_donor / distance_rectangle_donor +  patchwork::plot_layout(heights = c(8, 2))
  plot_PC_donor
  
  ###############################
  ## PROTEIN CODING & ACCEPTOR
  ###############################
  
  plot_PC_acceptor <- ggplot(data = df_novel_tidy %>% 
                               filter(type_PC == "protein coding (PC)",
                                      novel_type == "Novel Acceptor")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack") +
    ggtitle(" ") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    #ylim(c(0,250))+
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
    
  
  plot_PC_acceptor <- plot_PC_acceptor + 
    ggpubr::rremove("y.axis") + 
    ggpubr::rremove("ylab")+ 
    ggpubr::rremove("y.ticks") + 
    ggpubr::rremove("y.text") 
  
  
  plot_PC_acceptor <- plot_PC_acceptor / distance_rectangle_acceptor +  patchwork::plot_layout(heights = c(8, 2))
  plot_PC_acceptor
  
  
  ########################################
  ## PROTEIN CODING - DONOR & ACCEPTOR
  ########################################
  
  plot_PC <- ggpubr::ggarrange(plot_PC_donor,
                               plot_PC_acceptor,
                               ncol = 2,
                               nrow = 1)
  
  
  plot_PC
  
  file_name <- paste0(figures_folder, "/FCTX_distancesPC_subsampled", subsample)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                  width = 180, height = 90, 
                  units = "mm", dpi = 300)
  
  
  
  ###############################
  ## NONPROTEIN CODING & DONOR
  ###############################
  
  plot_NPC_donor <- ggplot(data = df_novel_tidy %>% 
                            filter(type_PC == "non PC",
                                   novel_type == "Novel Donor")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack") +
    #ylim(c(0,250))+
    ggtitle("Protein-coding transcripts") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c(limit_bp,(limit_bp * -1)),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6),
                       trans = "reverse") +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  plot_NPC_donor <- plot_NPC_donor / distance_rectangle_donor +  patchwork::plot_layout(heights = c(8, 2))
  plot_NPC_donor
  
  ###############################
  ## NONPROTEIN CODING & ACCEPTOR
  ###############################
  
  plot_NPC_acceptor <- ggplot(data = df_novel_tidy %>% 
                               filter(type_PC == "non PC",
                                      novel_type == "Novel Acceptor")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack") +
    ggtitle(" ") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    #ylim(c(0,250))+
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  
  plot_NPC_acceptor <- plot_NPC_acceptor + 
    ggpubr::rremove("y.axis") + 
    ggpubr::rremove("ylab")+ 
    ggpubr::rremove("y.ticks") + 
    ggpubr::rremove("y.text") 
  
  
  plot_NPC_acceptor <- plot_NPC_acceptor / distance_rectangle_acceptor +  patchwork::plot_layout(heights = c(8, 2))
  plot_NPC_acceptor
  
  
  ########################################
  ## PROTEIN CODING - DONOR & ACCEPTOR
  ########################################
  
  plot_NPC <- ggpubr::ggarrange(plot_NPC_donor,
                                plot_NPC_acceptor,
                               ncol = 2,
                               nrow = 1)
  
  
  plot_NPC
  
  file_name <- paste0(figures_folder, "/FCTX_distancesNPC_subsampled", subsample)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                  width = 180, height = 90, 
                  units = "mm", dpi = 300)
  
  ###############################
  ## PROTEIN CODING & NONCODING 
  ###############################
  
  ## GENERATE LEGEND
  
  plot <- ggpubr::ggarrange(plot_PC, 
                            plot_NPC,
                            common.legend = T,
                            legend = "top",
                            legend.grob = plot_legend,
                            labels = c("c", "d"),
                            #align = "v",
                            ncol = 1, 
                            nrow = 2,
                            font.label = list(size = 14, color = "black", face = "bold", family = NULL))
  
  
  plot
  
  
  file_name <- paste0(figures_folder, "/FCTX_distances_biotype_subsampled", subsample)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)

  
}

get_distances_all_tissues <- function() {
  
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  limit_bp <- 30
  
  
  if (!file.exists(paste0(results_folder, "/05_supplementary_table2.csv"))) {
    
    df_distances_all_tissues <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      
      ## GET THE CLUSTERS FOR THE CURRENT TISSUE
      
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
    
    
    write.csv(x = df_distances_all_tissues,
              file = paste0(results_folder, "/05_supplementary_data_distances_all_tissues.csv"),
              row.names = F)
  } else {
    df_distances_all_tissues <- read_delim(file = paste0(results_folder, "/05_supplementary_data_distances_all_tissues.csv"))
  }
  
  df_distances_all_tissues$mode_negative_distance_donor
  df_distances_all_tissues$mode_positive_distance_donor
  
  df_distances_all_tissues$mode_negative_distance_acceptor
  df_distances_all_tissues$mode_positive_distance_acceptor
  
  df_distances_all_tissues %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)")
  
}

get_modulo <- function() {
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  file_path <- paste0(results_folder, "/df_modulo_basic_tissues_100bpfilter.rds")
  
  if ( !file.exists(file_path) )  {
    
    df_modulo_tissues <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      print(paste0(Sys.time(), " - ", project_id))
      
      all_clusters <- df_metadata %>%
        filter(SRA_project == project_id) %>%
        distinct(cluster) %>%
        pull()
      
      map_df(all_clusters, function(cluster_id) {
        
        # cluster_id <- all_clusters[1]
        
        ## Print the tissue
        print(paste0(Sys.time(), " - ", cluster_id))
        
        #####################################
        ## Get the novel junctions from the current tissue
        #####################################
        
        query <- paste0("SELECT novel_junID 
                        FROM '", cluster_id, "_", project_id, "_misspliced'")
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
          mutate(novel_type = str_replace(string = novel_type,
                                          pattern = "_",
                                          replacement = " ")) %>%
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
    
    saveRDS(object = df_modulo_tissues, file = file_path)
    
  } else {
    
    df_modulo_tissues <- readRDS(file = file_path) %>%
      as_tibble()
    
  }
  
  df_modulo_tissues$modulo = factor(df_modulo_tissues$modulo, 
                                    levels = c( "0", "1", "2"))
  
  
  ################
  ## DENSITY PLOT
  ################
  
  # df_modulo_tissues <- df_modulo_tissues %>%
  #   mutate(freq = freq * 100)
  
  ggplot(df_modulo_tissues, aes(x = freq, y = modulo)) +
    ggridges::geom_density_ridges_gradient() +
    ylab("Modulo3 of the distance") +
    xlab("% of novel junctions") +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(30,35,40),
                       labels = c("30%","35%","40%")) +
    scale_y_discrete(expand = c(0,0.5,1,0))
  
  
  file_name <- paste0(figures_folder, "/panel3e")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 50, units = "mm", dpi = 300)
  
  
  file_name <- paste0(figures_folder, "/panel3e_large")
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 70, units = "mm", dpi = 300)
  
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

get_MES_similarity_cryptic_splice_sites <- function() {
  
  
  ## Get the delta MES

  file_name <- paste0(results_folder, "/df_mes.rds")
  
  if ( !file.exists(file_name) )  {
    
    ###########################
    ## GET THE MES
    ###########################
    
    ## Get the MES scores corresponding to the MASTER mis-spliced introns

    introns <- master_introns %>% 
      filter(misspliced == "Yes") %>%
      dplyr::select(ref_junID, ref_mes5ss, ref_mes3ss) %>% 
      as_tibble()
    
    
    ## Get the scores corresponding to the MASTER novel junctions

    novel_junctions <- master_novel_junctions %>%
      dplyr::select(ref_junID, novel_junID, novel_type, novel_mes5ss, novel_mes3ss, distance ) %>% 
      as_tibble() 
    
    ###########################
    ## MERGE THE DATA AND RETURN
    ###########################
    
    df_mes <- novel_junctions %>% 
      left_join(y = introns,
                by = "ref_junID")  %>%
      mutate(diff_ss5score = ref_mes5ss - novel_mes5ss,
             diff_ss3score = ref_mes3ss - novel_mes3ss)
    
    ## Save Data
    saveRDS(object = df_mes %>% as_tibble(), file = file_name)
    
  } else {
    df_mes <- readRDS(file = file_name) %>% 
      as_tibble()
  }
  
  
  df_mes_tidy <- df_mes %>%
    distinct(ref_junID, novel_junID, .keep_all = T) %>%
    filter(abs(distance) <= 75) %>%
    mutate(mod3 = abs(distance)%%3)
  
  df_mes_tidy$novel_type = factor(df_mes_tidy$novel_type, 
                                  levels = c("novel_donor", "novel_acceptor"))
  
  df_mes_tidy$novel_type = factor(df_mes_tidy$novel_type, 
                                  levels = c("novel_donor", "novel_acceptor"))
  
  #################################################
  ## PLOT DATA
  #################################################
  
  df_mes_tidy_plot <- df_mes_tidy %>%
    dplyr::select( mod3, 
                   "Delta MES 5ss" = diff_ss5score, 
                   "Delta MES 3ss" = diff_ss3score, 
                   novel_type, 
                   distance) %>%
    mutate(mod3 = mod3 %>% as.factor()) %>%
    gather( type, score, -novel_type, - mod3, -distance) %>%
    filter( (novel_type == "novel_donor" & type == "Delta MES 5ss") |
              (novel_type == "novel_acceptor" & type == "Delta MES 3ss")) %>%
    #mutate(type = type %>% factor(levels = c("Delta MES 5ss","Delta MES 3ss"))) %>%
    #mutate(mod3 = mod3 %>% factor(levels = c("2","1","0"))) %>%
    group_by(type,mod3) %>% 
    mutate(mean_score = mean(score)) %>% 
    ungroup()
  
  group_means  <- df_mes_tidy_plot %>%
    group_by(type,mod3) %>%
    summarise(mean = mean(score))
  
  
  # p-values --------------------------------------------------------------------------------------------------
  stat.test.donor <- tibble::tribble(
    ~group1, ~group2,   ~p.adj,
    "0",     "1-2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 == 0) %>% pull(diff_ss5score),
                                 y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 > 0) %>% pull(diff_ss5score),
                                 alternative = "less",
                                 correct = T))$p.value %>% 
      formatC(format = "e", digits = 2)#,
    # "0",     "2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 == 0) %>% pull(diff_ss5score),
    #                              y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 == 2) %>% pull(diff_ss5score),
    #                              correct = T))$p.value  %>% formatC(format = "e", digits = 2),
    # "1",     "2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 == 1) %>% pull(diff_ss5score),
    #                            y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_donor", mod3 == 2) %>% pull(diff_ss5score),
    #                            correct = T))$p.value %>% formatC(format = "e", digits = 2),
  )
  stat.test.acceptor <- tibble::tribble(
    ~group1, ~group2,   ~p.adj,
    "0",     "1-2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 == 0) %>% pull(diff_ss3score),
                                 y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 > 0) %>% pull(diff_ss3score),
                                 alternative = "less",
                                 correct = T))$p.value  %>% formatC(format = "e", digits = 2)#,
    # "0",     "2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 == 0) %>% pull(diff_ss3score),
    #                            y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 == 2) %>% pull(diff_ss3score),
    #                            correct = T))$p.value %>% formatC(format = "e", digits = 2),
    # "1",     "2", (wilcox.test(x = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 == 1) %>% pull(diff_ss3score),
    #                            y = df_mes_tidy %>% dplyr::filter(novel_type == "novel_acceptor", mod3 == 2) %>% pull(diff_ss3score),
    #                            correct = T))$p.value %>% formatC(format = "e", digits = 2)
  )
  
  
  # eff.sizes --------------------------------------------------------------------------------------------------
  
  eff.size.donor <- rstatix::wilcox_effsize(data = df_mes_tidy %>% 
                                              dplyr::filter(novel_type == "novel_donor") %>%
                                              mutate(mod3 = mod3 %>% as.factor()),
                                            formula = diff_ss5score ~ mod3 )
  
  eff.size.acceptor <- rstatix::wilcox_effsize(data = df_mes_tidy %>% 
                                                 dplyr::filter(novel_type == "novel_acceptor") %>%
                                                 mutate(mod3 = mod3 %>% as.factor()),
                                               formula = diff_ss3score ~ mod3 )
  
  eff.size.acceptor
  
  
  ####################################################
  ## BOXPLOT
  ####################################################
  
  ggplot(df_mes_tidy_plot %>%
           mutate(mod3 = ifelse(mod3 == "0", "0", "1-2"))) +
    
    geom_boxplot(aes(y = score, x = mod3, fill = mod3), alpha = 0.9)  +
    ggpubr::stat_pvalue_manual(
      rbind(stat.test.donor %>%
              mutate(type = "Delta MES 5ss"),
            stat.test.acceptor  %>%
              mutate(type = "Delta MES 3ss")) %>%
        mutate(p.scient = format(p.adj, scientific = TRUE)), 
      y.position = 35, 
      size = 3,
      step.increase = 0.1,
      label = "p.scient"
      
    ) +
    facet_grid(~fct_rev(type)) +
    theme_light() +
    custom_ggtheme +
    #ylim(c(-40,95)) +
    ylab("Delta MES score") +
    xlab("Modulo3") +
    theme(legend.title = element_text(size = "7", colour = "black")) +
    guides(fill = guide_legend(title = "Modulo3")) +
    scale_fill_discrete(breaks = c("0", "1-2")) 
  
  file_name <- paste0(figures_folder, "/deltames_modulo_100bp_exons_boxplot_one_tailed")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
}





## 6. Splicing error rates vary across introns and are likely to be underestimated in bulk RNA-seq data 

get_MSR_FCTX <- function()  {
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                  FROM '", cluster_id, "_", project_id, "_nevermisspliced' ")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                  FROM '", cluster_id, "_", project_id, "_misspliced' ")
  introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
  
  introns <- introns %>%
    left_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding) %>% as_tibble(),
              by = "ref_junID") %>% 
    as_tibble() 
  
  introns
  
  print(paste0(Sys.time(), " - ", introns %>% nrow(), " - introns!"))
  
  ###########################
  ## PLOT ALL INTRONS
  ###########################
  
  ## 1. Compare mis-spliced vs never mis-spliced introns
  
  df_all_introns_tidy <- introns %>%
    dplyr::select(ref_junID, MSR_D, MSR_A, ref_type) %>%
    gather(key = "MSR_type", value = "MSR", -ref_junID, -ref_type) %>%
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
    ggplot2::labs(x = "", y = "Number of annotated introns")+
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  ## 2. Mis-splicing is more common at the acceptor
  
  df_all_introns_tidy$percentile_group = factor(df_all_introns_tidy$percentile_group, 
                                                levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  plot2 <- ggplot(data = df_all_introns_tidy %>%
                    filter(percentile_group != "0")) + 
    geom_bar(aes(x = percentile_group, fill = MSR_type),position = "dodge")+
    #ggtitle("All annotated introns") +
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    ggforce::facet_zoom(ylim=c(0,2500), split = TRUE) +
    #ylim(c(0,23000))+
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR Donor","MSR Acceptor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  

  
  ggpubr::ggarrange(plot1,
                    plot2,
                    labels = c("", ""),
                    nrow = 2, ncol = 1,
                    heights = c(1,1.5),
                    common.legend = T)
  

  
  file_name <- paste0(figures_folder, "/panel4a")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 75, units = "mm", dpi = 300)
  
  
  file_name <- paste0(figures_folder, "/panel4a_large")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 110, units = "mm", dpi = 300)
  
  
  ##############################################################################
  ## CORRECT BY PER-INTRON MEAN EXPRESSION LEVELS
  ##############################################################################
  
  
  ## TIDY DATAFRAME
  df_biotype_result_tidy <- introns %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC"))  %>%
    dplyr::rename(biotype = type_PC) %>%
    distinct(ref_junID, .keep_all = T) %>%
    group_by(ref_junID) %>%
    mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
    mutate(mean_coverage = mean_coverage %>% log10()) %>%
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
  
  data_combined <- rbind(df_protein_coding_tidy, df_lncRNA_tidy) 
  
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
            file = paste0(results_folder, "/MSR_subsample.rds"))
  } else {
    subsample <- readRDS(file = paste0(results_folder, "/MSR_subsample.rds"))
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
  
  file_name <- paste0(figures_folder, "/MSR_FCTX_subsampling")
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
  
  
  ################################
  ## TIDY DATA BEFORE PLOTTING
  ################################
  
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
  
  
  ###########################
  ## PLOT BY MSR GROUP
  ###########################
  
  # scales::show_col(scales::viridis_pal(option = "B")(40))
  # scales::show_col(scales::viridis_pal(option = "C")(40))
  # scales::show_col(scales::viridis_pal(option = "D")(40))
  # scales::show_col(scales::viridis_pal(option = "E")(40))
  # scales::show_col(scales::viridis_pal(option = "F")(40))
  # scales::show_col(scales::viridis_pal(option = "G")(40))
  # scales::show_col(scales::viridis_pal(option = "H")(40))
  
  
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
  
  plotMSR_donor <- ggplot(data = df_introns_biotype_tidy %>% 
                            filter(MSR_type == "MSR_D")) + 
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
  plotMSR_acceptor <- ggplot(data = df_introns_biotype_tidy %>% 
                               filter(MSR_type == "MSR_A")) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = 0.5, color = "#333333")+
    ggtitle("MSR Acceptor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    #ylim(c(0,23000))+
    
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(ncol = 2, nrow = 1))
  
  
  ggpubr::ggarrange(plotMSR_donor,
                    plotMSR_acceptor,
                    nrow = 1,
                    ncol = 2,
                    common.legend = T)#,
                    #labels = c("b", "c"))
  
  
  
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0(figures_folder, "/panel4bc")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
  
  ####################################
  
  
  
  plotMSR_donor_zoomed <- ggplot(data = df_introns_biotype_tidy %>% 
                                   filter(MSR_type == "MSR_D")) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    ggtitle("MSR Donor") +
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    scale_y_continuous(limits =c(0,800), position = "right") +
    
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL,  
                               ncol = 2, nrow = 1)) 
  plotMSR_donor_zoomed
  
  
  ## Non-protein-coding
  plotMSR_acceptor_zoomed <- ggplot(data = df_introns_biotype_tidy %>% 
                                   filter(MSR_type == "MSR_A")) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    ggtitle("MSR Acceptor") +
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    scale_y_continuous(limits =c(0,800), position = "right") +
    
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, 
                               ncol = 2, nrow = 1))
  plotMSR_acceptor_zoomed
  
  
  ggpubr::ggarrange(plotMSR_donor_zoomed,
                    plotMSR_acceptor_zoomed,
                    nrow = 1,
                    ncol = 2,
                    common.legend = T)#,
                    #labels = c("b", "c"))
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0(figures_folder, "/panel4bc-zoomed")
  ggplot2::ggsave(paste0(file_name, ".png"), 
                  width = 180, height = 90, units = "mm", dpi = 300)
  
  
  ####################################
  ## STATS
  ####################################
  
 
  ## 1. We started by focusing on frontal cortex and observed that while splicing errors are detected 
  ## infrequently with the MSRD and MSRA values highly skewed towards low values, there was considerable 
  ## variation across introns (MSRA interquartile range = XXXXX; MSRD interquartile range = XXXX)
  
  introns %>%
    pull(MSR_D) %>% 
    IQR()
  
  
  introns %>%
    pull(MSR_A) %>% 
    IQR()
  
  
  
  ## 3. Given that NMD activity would be expected to reduce the detection of splicing errors amongst mRNA transcripts, 
  ## we separately assessed mis-splicing of annotated introns solely used in protein-coding transcripts (N=35,726 based 
  ## on Ensembl v105) and those that were exclusively used in non-coding transcripts (N=8,945 based on Ensembl v105).
  
  introns %>% filter(protein_coding == 100) %>% distinct(ref_junID) %>% nrow()
  introns %>% filter(protein_coding == 0) %>% distinct(ref_junID) %>% nrow()
  
  df_introns_biotype %>% filter(biotype == "PC") %>% distinct(ref_junID) %>% nrow()
  df_introns_biotype %>% filter(biotype == "non PC") %>% distinct(ref_junID) %>% nrow()
  
  
  ## We found that mis-splicing was more frequent amongst annotated introns from non-protein-coding transcripts at both the 5’ss
  ## 2. Splicing noise is more frequent in introns from non-protein coding transcripts 
  ## than in introns from protein-coding transcripts even after correcting by mean read coverage
  
  ## splicing noise at the donor is more frequent in introns from NPC vs PC at the 5'ss
  wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_D") %>% pull(MSR),
              y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_D") %>% pull(MSR),
              alternative = "less",
              paired = T,
              correct = T)
  rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                            filter(MSR_type == "MSR_D")%>%
                            mutate(biotype = biotype %>% as.factor()),
                          formula = MSR ~biotype,
                          paired = T)
  
  ## splicing noise at the donor is more frequent in introns from NPC vs PC at the 3'ss
  wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_A") %>% pull(MSR),
              y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_A") %>% pull(MSR),
              alternative = "less",
              paired = T,
              correct = T)
  rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                            filter(MSR_type == "MSR_A")%>%
                            mutate(biotype = biotype %>% as.factor()),
                          formula = MSR ~ biotype,
                          paired = T)
  
}

get_MSR_tissues <- function()  {
 
  ############################
  ## LOAD/GENERATE DATA
  ############################
  
  file_name <- paste0(results_folder, "/05_supplementary_table3.csv")
  
  print(paste0(Sys.time(), " - getting MSR across tissues..."))
  
  if ( !file.exists(file_name) ) {
    
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
    
  } else {
    MSR_all_tissues <- read.csv(file = file_name)
  }
  
}





## 7. Local sequence conservation is the most important predictor of mis-splicing

get_ZIP_single_tissue <- function(project_id = "BRAIN",
                                  cluster_id = "Brain - Frontal Cortex (BA9)",
                                  intron_size = 100,
                                  phastcons_type = 17,
                                  common_introns = F) {
  
  
  
  ###############################
  ## QUERY THE DATABASE
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  if ( !file.exists(paste0(results_folder,
                           "/", cluster_id,
                           "_common_introns", common_introns,
                           "_introns_phastcons", phastcons_type,
                           "_intronsize", intron_size ,".rds")) ) {
    
    message("Intron size: ", intron_size, "...")
    
    
    ## GET MIS-SPLICED INTRONS AND JOIN WITH NOVEL MASTER DATA
    query <- paste0("SELECT *
                    FROM '", cluster_id, "_", project_id, "_misspliced'")
    introns_ms <- dbGetQuery(con, query) %>% as_tibble()

    
    ## JOIN WITH MASTER NOVEL TABLE
    introns_ms <- introns_ms %>%
      inner_join(y = master_novel_junctions %>% 
                   dplyr::select(novel_junID, novel_type) %>% as_tibble(),
                 by = c("novel_junID")) %>% 
      as_tibble() 
    
    
    ## GET NEVER MIS-SPLICED INTRONS AND JOIN WITH MIS-SPLICED DATA
    query <- paste0("SELECT * 
                    FROM '", cluster_id, "_", project_id, "_nevermisspliced' " )
    introns <- dbGetQuery(con, query) %>%
      as_tibble()
    
    

    ## JOIN WITH MISSPLICED INTRONS
    introns <- plyr::rbind.fill(introns, introns_ms ) %>% as_tibble() 
    
    
    ###################################
    ## JOIN WITH MASTER INTRON TABLE
    ###################################
    
     
    # query <- paste0("SELECT seqnames, start, end, strand, ref_junID, ref_coordinates, ref_length, mean_phastCons", 
    #                 phastcons_type, "way5ss_", intron_size,", mean_phastCons", phastcons_type, "way3ss_", 
    #                 intron_size, ", mean_CDTS5ss_", intron_size,", mean_CDTS3ss_", intron_size,
    #                 ", ref_mes5ss, ref_mes3ss, protein_coding
    #                 FROM 'intron' " )
    
    
    introns <- introns %>%
      inner_join(y = master_introns %>%
                   filter(ref_length >= (intron_size * 2)) %>%
                   dplyr::select(-transcript_id),
                 by = "ref_junID") %>% 
      as_tibble()
    
    
    # if ( common_introns ) {
    #   all_common_introns <- readRDS(file = paste0(results_folder, "/common_introns/common_introns_all_tissues.rds"))
    #   introns <- introns %>%
    #     filter(ref_junID %in% all_common_introns$ref_junID)
    # } 
    
    
    summary(introns)
    
    #####################################
    ## JOIN WITH MASTER TRANSCRIPT TABLE
    #####################################
    
    query <- paste0("SELECT * 
                    FROM 'transcript' 
                    WHERE id IN (", paste(introns$transcript_id, collapse = ","),")")
    introns <- introns %>%
      inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                 by = c("transcript_id" = "id")) %>% 
      as_tibble() 
    
    introns %>% nrow()
    introns %>% head()
    
    
    #####################################
    ## JOIN WITH MASTER GENE TABLE
    #####################################
    
    query <- paste0("SELECT *
                    FROM 'gene' 
                    WHERE id IN (", base::paste(introns$gene_id, collapse = ","),")")
    introns <- introns %>%
      inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                 by = c("gene_id" = "id")) %>% 
      as_tibble() 
    
    introns %>% nrow()
    introns %>% head()
    
    introns$ref_junID %>% unique %>% length()
    
    #########################################################
    ## PREPARE THE DATA PRIOR MODELLING ---------------------
    #########################################################
    
    
    # introns <- introns %>% dplyr::mutate(novel_type = replace_na(novel_type, 0)) 
    
    message("PhastCons type: ", phastcons_type)
    
 
    
    idb <- introns %>%
      distinct(ref_junID, novel_junID, .keep_all = T)  %>%
      dplyr::select(ref_junID,
                    novel_junID,
                    novel_type, 
                    ref_type,
                    intron_length = ref_length,
                    intron_5ss_score = ref_mes5ss,
                    intron_3ss_score = ref_mes3ss,
                    gene_length = gene_width,
                    gene_tpm = gene_tpm,
                    gene_num_transcripts = n_transcripts,
                    mean_CDTS5ss = paste0("mean_CDTS5ss_", intron_size),
                    mean_CDTS3ss = paste0("mean_CDTS3ss_", intron_size),
                    protein_coding,
                    mean_phastConsway5ss = paste0("mean_phastCons",phastcons_type,"way5ss_",intron_size),
                    mean_phastConsway3ss = paste0("mean_phastCons",phastcons_type,"way3ss_",intron_size),
                    
                    MSR_D,
                    MSR_A,
                    
                    ref_sum_counts,
                    novel_sum_counts) 
    
    
    idb_tidy <- idb %>%
      mutate(MSR_D = MSR_D * 100000,
             MSR_A = MSR_A * 100000) %>%
      mutate(MSR_D = MSR_D %>% round(digits = 0),
             MSR_A = MSR_A %>% round(digits = 0)) %>% 
      dplyr::select(ref_junID,
                    gene_length,
                    gene_tpm,
                    gene_num_transcripts,
                    protein_coding,
                    intron_length, 
                    intron_5ss_score, 
                    intron_3ss_score,
                    mean_CDTS5ss,
                    mean_CDTS3ss,
                    mean_phastConsway5ss, 
                    mean_phastConsway3ss,
                    
                    MSR_D,
                    MSR_A)  %>%
      distinct(ref_junID,.keep_all = T) %>%
      drop_na() %>%
      filter(gene_tpm > 0)
    
    
    
    summary(idb_tidy)
    
    ## SAVE DATA
    message("Saving 'idb' data from ", cluster_id, "...")
    saveRDS(object = idb_tidy,
            file = paste0(results_folder,"/", 
                          cluster_id,
                          "_common_introns", common_introns, 
                          "_introns_phastcons", phastcons_type, 
                          "_intronsize", intron_size ,".rds"))
  } 
  else {
    
    message("Loading 'idb' data from ", cluster_id, "...")
    idb_tidy <- readRDS(file =  paste0(results_folder,"/", cluster_id,
                                       "_common_introns", common_introns,
                                       "_introns_phastcons", phastcons_type,
                                       "_intronsize", intron_size ,".rds")) %>%
      drop_na()
    message("Data loaded!")
  }
  
  
  
  
  
  ######################################################################################################
  ## ZERO-INFLATED NEGATIVE BINOMIAL POISSION MODELS ---------------------------------------------------
  ######################################################################################################
  
  message(Sys.time(), " - fitting a ZINB Donor model...")
  
  
  zeroinfl_model_D <- pscl::zeroinfl(MSR_D ~ 
                                       #gene_tpm +
                                       #gene_length +
                                       gene_num_transcripts +
                                       protein_coding + 
                                       intron_length +
                                       intron_5ss_score +
                                       intron_3ss_score +
                                       mean_CDTS5ss +
                                       mean_CDTS3ss +
                                       mean_phastConsway5ss +
                                       mean_phastConsway3ss,
                                     data = idb_tidy %>%
                                       distinct(ref_junID,.keep_all = T) %>%
                                       mutate(intron_length = intron_length %>% log(),
                                              gene_length = gene_length %>% log(),
                                              gene_tpm = gene_tpm %>% log()))
  zeroinfl_model_D %>% summary %>% print()
  
  saveRDS(object = zeroinfl_model_D,
          file = paste0(results_folder,
                        "/ZINB_Donor_", cluster_id, 
                        "_common_introns", common_introns, 
                        "_introns_phastcons", phastcons_type, 
                        "_intronsize", intron_size ,".rds") )
  
  message(Sys.time(), " - fitting a ZINB Acceptor model...")
  
  zeroinfl_model_A <- pscl::zeroinfl(MSR_A ~ 
                                       # gene_tpm +
                                       # gene_length +
                                       gene_num_transcripts +
                                       protein_coding + 
                                       intron_length +
                                       intron_5ss_score +
                                       intron_3ss_score +
                                       mean_CDTS5ss +
                                       mean_CDTS3ss +
                                       mean_phastConsway5ss +
                                       mean_phastConsway3ss,
                                     data = idb_tidy %>%
                                       distinct(ref_junID,.keep_all = T) %>%
                                       mutate(intron_length = intron_length %>% log(),
                                              gene_tpm = gene_tpm %>% log()))
  
  zeroinfl_model_A %>% summary %>% print()
  
  
  saveRDS(object = zeroinfl_model_A,
          file = paste0(results_folder,
                        "/ZINB_Acceptor_", cluster_id, 
                        "_common_introns", common_introns, 
                        "_introns_phastcons", phastcons_type, 
                        "_intronsize", intron_size , ".rds"))
  
}

plot_ZIP_single_tissue <- function() {
  
  
  project_id = "BRAIN"
  cluster_id = "Brain - Frontal Cortex (BA9)"
  intron_size = 100
  phastcons_type = 17
  common_introns = F
  
  
  ##########################
  ## LOAD THE DATA
  ##########################
  
  
  ZINB_tissue_MSR <- map_df(c("Donor", "Acceptor"), function(novel_type)  { 
    
    
    # novel_type <-  "Donor"
    # novel_type <-  "Acceptor"
    
    message(Sys.time(), " - ", novel_type)
    ZINB_tissue <- readRDS(file = paste0(results_folder, "/",
                                         "/ZINB_", novel_type, "_", cluster_id, 
                                         "_common_introns", common_introns, 
                                         "_introns_phastcons", phastcons_type, 
                                         "_intronsize", intron_size ,".rds"))
    
    ZINB_corrected_pvals <- lmtest::coeftest(ZINB_tissue, vcov = sandwich::sandwich)
    ZINB_corrected_pvals
    
    
    ##############################################################################
    ## EXP COEFFICIENTS
    ##############################################################################
    
    ## Exponentiated log coefficients
    ZINB_expCoef <- (coef(ZINB_corrected_pvals))
    ZINB_expCoef <- as.data.frame(ZINB_expCoef, ncol = 2) %>%
      tibble::rownames_to_column("covariate") %>%
      mutate(sign = (sign(x = coef(ZINB_corrected_pvals) %>% unlist() %>% unname())))
    
    
    ## Get confidence intervals
    ZINB_confInt <- confint(object = ZINB_corrected_pvals,  level=0.95) %>%
      as.data.frame(ZINB_expCoef, ncol = 3) #%>%
    #mutate(`2.5 %` = `2.5 %` * (sign(confint(object = ZINB_corrected_pvals,  level=0.95)[,1] %>% unlist() %>% unname())) ) %>%
    #mutate(`97.5 %` = `97.5 %` * (sign(confint(object = ZINB_corrected_pvals,  level=0.95)[,2] %>% unlist() %>% unname())))
    
    
    ## Join data
    ZINB_plot <- cbind(ZINB_expCoef, ZINB_confInt) %>%
      dplyr::rename(log.estimate = ZINB_expCoef,
                    log.conf.low = "2.5 %",
                    log.conf.high = "97.5 %") %>%
      mutate(pvalue = c( ZINB_corrected_pvals[,4] %>% 
                           as.data.frame() %>%
                           tibble::rownames_to_column("covariate") %>%
                           filter(covariate != "Log(theta)") %>%
                           pull(.)))
    #(ZINB_table$coefficients$count %>% as.data.frame()) %>%
    #    ZINB_corrected_pvals %>% as.data.frame()
    #   tibble::rownames_to_column("covariate") %>%
    #   filter(covariate != "Log(theta)") %>%
    #   pull("Pr(>|z|)"),
    # (ZINB_table$coefficients$zero %>% as.data.frame())[,4]))
    
    ZINB_plot %>%
      mutate(MSR = novel_type)
  })
  
  
  ###########################
  ## TIDY RESULTS
  ###########################
  
  ZINB_tissue_MSR <- ZINB_tissue_MSR %>%
    mutate(log.estimate = ifelse(pvalue > 0.05, 0, log.estimate ),
           log.conf.low = ifelse(pvalue > 0.05, 0, log.conf.low ),
           log.conf.high = ifelse(pvalue > 0.05, 0, log.conf.high )) %>%
    mutate(log.estimate = ifelse(covariate == "count_intron_length", log.estimate/100, log.estimate ),
           log.conf.low = ifelse(covariate == "count_intron_length", log.conf.low/100, log.conf.low ),
           log.conf.high = ifelse(covariate == "count_intron_length", log.conf.high/100, log.conf.high ))  %>%
    mutate(log.estimate = log.estimate %>% exp(),
           log.conf.low = log.conf.low %>% exp(),
           log.conf.high = log.conf.high %>% exp() ) %>%
    mutate(log.estimate = log.estimate  * sign,
           log.conf.low = log.conf.low  * sign,
           log.conf.high = log.conf.high   * sign ) 
  
  ##########################################
  ## COUNT MODEL
  ##########################################
  
  coef_names <- c(#"Gene TPM", 
    "Num. transcripts", "Protein coding",
    "Intron Length", "MES 5'ss score", "MES 3'ss score",
    "CDTS 5'ss", "CDTS 3'ss", "PhastCons17 5'ss", "PhastCons17 3'ss")
  
  cov_names <- c(#"Gene-level",
    "Gene-level","Gene-level",
    "Intron-level","Intron-level","Intron-level","Intron-level", "Intron-level","Intron-level", "Intron-level")
  
  cov_order <- c(#"B",
    "C","D",
    "F","G","H","I","J","K","L")
  
  ZINB_plot_count <- ZINB_tissue_MSR %>%
    filter(str_detect(string = covariate, pattern = "tercept", negate = T),
           str_detect(string = covariate, pattern = "count")) %>%
    mutate(covariate_tidy = c(coef_names,coef_names))
  
  ZINB_plot_count_tidy <- ZINB_plot_count %>%
    mutate(estimate_lab = paste0(round(x = log.estimate,
                                       digits = 2), " (", 
                                 round(x = log.conf.low,
                                       digits = 2), ",", 
                                 round(x = log.conf.high,
                                       digits = 2), ")")) %>%
    dplyr::mutate(p.value = signif(pvalue, digits=3)) %>%
    mutate(p.value = case_when(pvalue == 0 ~ "<2.2e-16", pvalue > 0 ~ formatC(pvalue, format='e',digits = 2) %>% as.character())) %>%
    mutate(covariate_group = c(cov_names, cov_names),
           covariate_order = c(cov_order, cov_order) )
  
  ZINB_plot_count_tidy$covariate_order = factor(ZINB_plot_count_tidy$covariate_order, 
                                                levels = cov_order)
  
  
  ##########################################
  ## PLOT COUNT MODEL
  ##########################################
  
  
  p_mid <- ggplot(data = ZINB_plot_count_tidy %>%
                    filter(covariate_tidy != "Covariate"),
                  mapping = aes(y = fct_rev(covariate_order), color = MSR)) + 
    theme_classic() +
    geom_point(aes(x=log.estimate, group = MSR), shape=15, size=3, position = position_dodge(width = .75)) +
    geom_linerange(aes(xmin=log.conf.low, xmax=log.conf.high, group = MSR), position = position_dodge(width = .75)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    labs(x="Exp. Estimate", y="") +
    #coord_cartesian(ylim=c(1,12), xlim=c(-.5, 2)) +
    #annotate("text", x = -.3, y = 12, label = "Decrease in mis-splicing activity") +
    #annotate("text", x = 0.5, y = 12, label = "Increase in mis-splicing activity") + 
    ggforce::facet_col(~covariate_group, scales = "free_y", space = 'free',
                       strip.position = "left") + 
    theme(axis.line.y = element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y= element_blank(),
          axis.title.y= element_blank(),
          strip.text.y = element_blank(),
          strip.text = element_text(size = 3),
          legend.position = "top")+
    scale_color_manual(values = c( "#35B779FF", "#8d03b0"),
                       breaks = c( "Donor", "Acceptor"),
                       labels = c( "Novel Donor", "Novel Acceptor"))  + 
    theme( plot.margin = margin(t = 0,
                                r = -5,
                                b = 0,
                                l = -5),
           legend.box.margin=margin(b = -10, t = -5))#+
  p_mid 
  
  
  p_left <- ZINB_plot_count_tidy %>%
    mutate(covariate_tidy = ifelse(MSR == "Acceptor", "", covariate_tidy)) %>%
    ggplot(aes(y = fct_rev(covariate_order))) +
    geom_text(aes(x = 0, label = covariate_tidy, group= MSR),
              hjust = 0, 
              size = 2.5, position = position_dodge(width = .75),
              fontface = "bold") +
    geom_text(aes(x = 1, label = estimate_lab, group= MSR), 
              hjust = 0,
              size = 2.5, position = position_dodge(width = .75),
              fontface = ifelse(ZINB_plot_count_tidy$estimate_lab == "[95% CI]", "bold", "plain")) +
    theme_void() +
    ggforce::facet_col(~covariate_group,scales = "free_y",
                       space = 'free',strip.position = "left") + 
    coord_cartesian(xlim = c(0,2)) +
    theme(title =  element_text(size = 8, colour = "black", face = "bold"),
          strip.text = element_text(size = 7))
  p_left
  
  
  p_right <- ZINB_plot_count_tidy %>%
    ggplot(aes(y = fct_rev(covariate_order))) +
    geom_text( aes(x = 0, 
                   y = fct_rev(covariate_order), 
                   label = p.value,
                   group = MSR),
               hjust = 0, position = position_dodge(width = .75), 
               size = 2.5,
               fontface = ifelse(ZINB_plot_count_tidy$p.value == "p-value", "bold", "plain") )  +
    ggforce::facet_col(~covariate_group, scales = "free_y",
                       space = 'free',strip.position = "left") + 
    theme_void() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y= element_blank(),
          axis.title.y= element_blank(),
          strip.text.y = element_blank(),
          title =  element_text(size = 8, 
                                colour = "black",
                                face = "bold"))
  p_right 
  
  
  library(patchwork)
  layout <- c(
    area(t = 0, l = 0, b = 30, r = 4), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
    area(t = 1, l = 5, b = 30, r = 9), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
    area(t = 0, l = 9.7, b = 30, r = 10)
  )
  # final plot arrangement
  p_left +
    ggtitle("  Covariate            [95% CI]") +
    p_mid +
    p_right +
    ggtitle("               p-value") + 
    patchwork::plot_layout(design = layout)  
  
  
  file_name <- paste0(figures_folder, "/ZINB_COUNT_",
                      project_id, "_", cluster_id, 
                      "_common_introns", common_introns, 
                      "_ZINB_phastcons", phastcons_type, 
                      "_intronsize", intron_size)
  ggplot2::ggsave(paste0(file_name, "_noOutliers.png"), 
                  width = 180, height = 90, units = "mm", dpi = 300)
  
  
  
  ##########################################
  ## ZERO MODEL
  ##########################################
  
  # ZINB_plot_zero <- ZINB_plot %>%
  #   filter(str_detect(string = covariate, pattern = "tercept", negate = T),
  #          str_detect(string = covariate, pattern = "zero")) %>%
  #   mutate(covariate_tidy = coef_names)
  # 
  # ZINB_plot_zero_tidy <- ZINB_plot_zero %>%
  #   mutate(estimate_lab = paste0(round(x = log.estimate,
  #                                      digits = 2), " (", 
  #                                round(x = log.conf.low,
  #                                      digits = 2), ",", 
  #                                round(x = log.conf.high,
  #                                      digits = 2), ")")) %>%
  #   # round p-values to two decimal places, except in cases where p < .001
  #   mutate(p.value = case_when(pvalue < .001 ~ "<.001",
  #                              round(pvalue, 2) == .05 ~ as.character(round(pvalue,3)),
  #                              pvalue < .01 ~ str_pad( # if less than .01, go one more decimal place
  #                                as.character(round(pvalue, 3)),
  #                                width = 4,
  #                                pad = "0",
  #                                side = "right"
  #                              ),
  #                              TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
  #                                as.character(round(pvalue, 2)),
  #                                width = 4,
  #                                pad = "0",
  #                                side = "right"
  #                              ))) %>%
  #   mutate(covariate_group = c(#"Gene-level",
  #                              "Gene-level","Gene-level",
  #                              "Intron-level","Intron-level","Intron-level","Intron-level", "Intron-level","Intron-level", "Intron-level"),
  #          covariate_order = c(#"B",
  #                              "C","D",
  #                              "F","G","H","I","J","K","L")) 
  # 
  # ZINB_plot_zero_tidy$covariate_order = factor(ZINB_plot_zero_tidy$covariate_order, 
  #                                              levels = c(#"A",
  #                                                #"B",
  #                                                "C","D",
  #                                                "F","G","H","I","J","K","L"))
  # 
  # 
  # 
  # # ZINB_plot_count_tidy <- rbind(plot_labels <- c("","","","","","Covariate","(95% CI)", "p-value"),
  # #                               ZINB_plot_count_tidy) 
  # # 
  # 
  # 
  # 
  # ##########################################
  # ## PLOT COUNT MODEL
  # ##########################################
  # 
  # 
  # p_mid <- ggplot(data = ZINB_plot_zero_tidy %>%
  #                   filter(covariate_tidy != "Covariate"),
  #                 mapping = aes(y = fct_rev(covariate_order))) + 
  #   theme_classic() +
  #   geom_point(aes(x=log.estimate), shape=15, size=3) +
  #   geom_linerange(aes(xmin=log.conf.low, xmax=log.conf.high)) +
  #   geom_vline(xintercept = 0, linetype="dashed") +
  #   labs(x="Exp. Estimate", y="") +
  #   #coord_cartesian(ylim=c(1,12), xlim=c(-.5, 2)) +
  #   #annotate("text", x = -.3, y = 12, label = "Decrease in mis-splicing activity") +
  #   #annotate("text", x = 0.5, y = 12, label = "Increase in mis-splicing activity") + 
  #   ggforce::facet_col(~covariate_group, scales = "free_y", strip.position = "left") + 
  #   theme(axis.line.y = element_blank(),
  #         axis.ticks.y= element_blank(),
  #         axis.text.y= element_blank(),
  #         axis.title.y= element_blank(),
  #         strip.text.y = element_blank())
  # p_mid
  # 
  # 
  # p_left <- ZINB_plot_zero_tidy %>%
  #   ggplot(aes(y = fct_rev(covariate_order))) +
  #   geom_text(aes(x = 0, label = covariate_tidy),
  #             hjust = 0, 
  #             size = 2.5,
  #             fontface = "bold") +
  #   geom_text(aes(x = 1, label = estimate_lab), 
  #             hjust = 0,
  #             size = 2.5,
  #             fontface = ifelse(ZINB_plot_zero_tidy$estimate_lab == "[95% CI]", "bold", "plain")) +
  #   theme_void() +
  #   ggforce::facet_col(~covariate_group,scales = "free_y",strip.position = "left") + 
  #   coord_cartesian(xlim = c(0,2)) +
  #   theme(title =  element_text(size = 8, 
  #                               colour = "black",
  #                               face = "bold"))
  # p_left
  # 
  # 
  # p_right <- ZINB_plot_zero_tidy %>%
  #   ggplot(aes(y = fct_rev(covariate_order))) +
  #   geom_text( aes(x = 0, y = fct_rev(covariate_order), 
  #                  label = p.value),
  #              hjust = 0, 
  #              size = 2.5,
  #              fontface = ifelse(ZINB_plot_zero_tidy$p.value == "p-value", "bold", "plain") )  +
  #   ggforce::facet_col(~covariate_group, scales = "free_y",strip.position = "left") + 
  #   theme_void() +
  #   theme(axis.line.y = element_blank(),
  #         axis.ticks.y= element_blank(),
  #         axis.text.y= element_blank(),
  #         axis.title.y= element_blank(),
  #         strip.text.y = element_blank(),
  #         title =  element_text(size = 8, 
  #                               colour = "black",
  #                               face = "bold"))
  # p_right 
  # 
  # 
  # library(patchwork)
  # layout <- c(
  #   area(t = 0, l = 0, b = 30, r = 4), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  #   area(t = 1, l = 5, b = 30, r = 9), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  #   area(t = 0, l = 10, b = 30, r = 10)
  # )
  # # final plot arrangement
  # p_left +
  #   ggtitle("Covariate            [95% CI]") +
  #   p_mid +
  #   p_right +
  #   ggtitle("  p-value") + 
  #   patchwork::plot_layout(design = layout)  
  # 
  # 
  # file_name <- paste0(figures_folder, "/ZINB_ZERO_",novel_type,"_", project_id, "_", cluster_id, 
  #                     "_common_introns", common_introns, 
  #                     "_phastcons", phastcons_type, 
  #                     "_intronsize", intron_size)
  # ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 70, units = "mm", dpi = 300)
  # }
  #   
  # }
  
  
  
  
  ##############################################################################
  
  # 
  # coef_names <- c(#"Gene Length" = "gene_length",
  #   "Gene TPM" = "gene_tpm",
  #   "Gene num. transcripts" = "gene_num_transcripts",
  #   "Protein coding" = "protein_coding",
  #   "Intron Length" = "intron_length",
  #   "Intron 5'ss MES score" = "intron_5ss_score",
  #   "Intron 3'ss MES score" = "intron_3ss_score",
  #   "CDTS 5'ss" = "CDTS_5ss",
  #   "CDTS 3'ss" = "CDTS_3ss",
  #   "PhastCons17 5'ss" = "mean_phastConsway5ss",
  #   "PhastCons17 3'ss" = "mean_phastConsway3ss")
  # 
  
  
  
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
    
    print(paste0(Sys.time(), " - ", project_id))
    
    ## GET THE CLUSTERS
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
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
      
      
      all_introns[[cluster_id]] <- introns$ref_junID
      
      print(paste0(Sys.time(), " - intron IDs collected from '", cluster_id, "'"))
    }
    
  }
  
  common_introns <- data.frame(ref_junID = Reduce(intersect,  all_introns))
  common_introns %>% head()
  common_introns %>% nrow()
  
  query <- paste0("SELECT * 
                  FROM 'intron' 
                  WHERE ref_junID IN (", paste(common_introns$ref_junID, collapse = ","), ")")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  
  file_name <- paste0(results_folder, "/common_introns_all_tissues.rds")
  saveRDS(object = introns, file = file_name)
  
}

get_ZIP_across_tissues <- function() {
  
  
  all_projects <- df_metadata$SRA_project %>% unique()
  
  ## LOAD COMMON INTRONS ACROSS TISSUES
  
  #results_folder <- file.path(base_folder, "/natcomms_review/results")
  #
  
  # project_id <- "BRAIN"
  # cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  phastcons_type <- 17
  intron_size <- 100
  
  for (project_id in all_projects) {
    
    # project_id <- all_projects[1]
    # project_id <- "PANCREAS"
    # project_id <- "BLOOD"
    
    message(Sys.time(), " - ", project_id)
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    for (cluster_id in all_clusters) {
      
      # cluster_id <- all_clusters[1]
      # cluster_id <- all_clusters[2]
      
      message(Sys.time(), " - ", cluster_id)
      
      get_ZIP_single_tissue(project_id,
                            cluster_id,
                            intron_size,
                            phastcons_type,
                            common_introns = F) 
    }
  }
}

plot_ZIP_variance_across_tissues <- function() {
  
  
  #############################
  ## LOAD RESULTS
  #############################
  
  phastcons_type <- 17
  intron_size <- 100
  common <- FALSE
  
  all_coefficient_tissues <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    # project_id <- "BRAIN"
    
    print(paste0(Sys.time(), " - ", project_id))
    
    ## GET THE CLUSTERS
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(all_clusters, function(cluster_id) { 
      
      # cluster_id <- all_clusters[1]
      # cluster_id <- "Esophagus - Mucosa"
      
      map_df(c("Donor", "Acceptor"), function(novel_type) { 
        
        # novel_type <- "Donor"
        if ( file.exists(paste0(results_folder, "/ZINB_",novel_type,"_",
                                cluster_id,
                                "_common_introns", common,
                                "_introns_phastcons", phastcons_type, 
                                "_intronsize", intron_size ,".rds")) )  {
          
          message(novel_type, " - ", cluster_id)
          
          ZINB_tissue <- readRDS(file = paste0(results_folder,
                                               "/ZINB_", novel_type,"_", cluster_id, 
                                               "_common_introns", common,
                                               "_introns_phastcons", phastcons_type, 
                                               "_intronsize", intron_size ,".rds"))
          ZINB_tissue %>% summary()
          message("'ZINB_tissue' loaded!")
          
          ZINB_expCoef <- (coef((lmtest::coeftest(ZINB_tissue, vcov = sandwich::sandwich)))) #* 
          #(sign(x = coef((lmtest::coeftest(ZINB_tissue, vcov = sandwich::sandwich))) %>% unlist() %>% unname()))
          
          ZINB_expCoef <- as.data.frame(ZINB_expCoef, ncol = 2) %>%
            tibble::rownames_to_column("covariate") %>%
            mutate(cluster = cluster_id) %>%
            mutate(pvalue = lmtest::coeftest(ZINB_tissue, vcov = sandwich::sandwich)[,4])
          # c( (ZINB_table$coefficients$count %>% as.data.frame()) %>%
          #                        tibble::rownames_to_column("covariate") %>%
          #                        filter(covariate != "Log(theta)") %>%
          #                        pull("Pr(>|z|)"),
          #                      (ZINB_table$coefficients$zero %>% as.data.frame())[,4]))
          
          ZINB_expCoef %>% 
            mutate(MSR = novel_type) %>%
            return()
        }
      })
    })
  })
  
  
  saveRDS(all_coefficient_tissues,
          file = paste0(results_folder, "/all_coefficient_tissues.rds"))
  
  for (MSR_type in c("Donor", "Acceptor")) {
    
    # MSR_type <- "Donor"
    # MSR_type <- "Acceptor"
    
    ## 1. Load the estimate variance across tissues
    all_coefficient_tissues_fdr <- all_coefficient_tissues %>%
      mutate(beta_sign = (sign(x = ZINB_expCoef))) %>%
      ## Intron length was log transformed prior zero-inflation modelling, so it needs to be divided by 100
      mutate(ZINB_expCoef = ifelse(covariate == "count_intron_length", ZINB_expCoef/100, ZINB_expCoef) ) %>%
      mutate(ZINB_expCoef = ZINB_expCoef %>% exp() ) %>%
      mutate(ZINB_expCoef = ZINB_expCoef *  beta_sign) %>%
      filter(MSR == MSR_type) %>%
      mutate(q = p.adjust(pvalue, method = "fdr"),
             ZINB_expCoef = ZINB_expCoef %>% as.double) %>%
      as_tibble()
    
    ## 2. only keep significant covariates
    all_coefficient_tissues_fdr <- all_coefficient_tissues_fdr %>%
      filter(q <= 0.05 | (is.nan(q))) %>%
      dplyr::select(-c(q, pvalue))
    
    
    ################################################################
    ## COEFFICIENT OF VARIATION
    ################################################################
    
    # df_glm_estimates %>%
    #   filter(q <= 0.05) %>%
    #   dplyr::count(feature)
    
    
    all_coefficient_tissues_tidy <- all_coefficient_tissues_fdr %>%
      #filter(q < 0.05) %>%
      #mutate(estimate = ifelse(q > 0.05, NA, estimate))%>%
      filter(str_detect(string = covariate, pattern = "count_"),
             str_detect(string = covariate, pattern = "Intercept",negate = T)) %>%
      dplyr::group_by(covariate) %>%
      spread(key = covariate, value = ZINB_expCoef) %>%
      ungroup()
    
    all_coefficient_tissues_tidy <- all_coefficient_tissues_tidy %>% 
      dplyr::rename("Gene num. transcripts" = "count_gene_num_transcripts",
                    #"Gene TPM" = "count_gene_tpm", 
                    "Intron Length" = "count_intron_length",
                    "Intron 5'ss MES score" = "count_intron_5ss_score",
                    "Intron 3'ss MES score" = "count_intron_3ss_score", 
                    "CDTS 5'ss" = "count_mean_CDTS5ss",
                    "CDTS 3'ss" = "count_mean_CDTS3ss", 
                    "PhastCons17 5'ss" = "count_mean_phastConsway5ss",
                    "PhastCons17 3'ss" = "count_mean_phastConsway3ss",
                    "Protein coding" = "count_protein_coding")  %>%
      gather(feature, coefficient,-cluster) %>% as_tibble()
    
    
    
    label_faceting <- data.frame(label = c(#"Gene TPM",
      "Gene num. transcripts",
      "Protein coding",
      "Intron Length",
      "Intron 5'ss MES score",
      "Intron 3'ss MES score", 
      "CDTS 5'ss",
      "CDTS 3'ss",
      "PhastCons17 5'ss",
      "PhastCons17 3'ss" ),
      subgroup = c(# "Gene Level",
        "Gene\nLevel",
        "Gene\nLevel",
        "Intron Level",
        "Intron Level",
        "Intron Level", 
        "Intron Level",
        "Intron Level",
        "Intron Level",
        "Intron Level"))
    
    all_coefficient_tissues_tidy <- all_coefficient_tissues_tidy %>%
      inner_join(y = label_faceting,
                 by = c("feature" = "label")) %>%
      drop_na()
    
    ## TIDY VALUES
    all_coefficient_tissues_tidy$feature <- factor( all_coefficient_tissues_tidy$feature, 
                                                    levels = c(# "Gene TPM",
                                                      "Gene num. transcripts",
                                                      "Protein coding",
                                                      "Intron Length",
                                                      "Intron 5'ss MES score",
                                                      "Intron 3'ss MES score", 
                                                      "CDTS 5'ss",
                                                      "CDTS 3'ss",
                                                      "PhastCons17 5'ss",
                                                      "PhastCons17 3'ss") %>% rev() )
  
    if (MSR_type == "Donor") {
      color_boxplot <- "#35B779FF"
    } else {
      color_boxplot <- "#8d03b0"
    }
    
    ## PLOT
    plotTissuesZINB <- ggplot(data = all_coefficient_tissues_tidy %>%
                                mutate(coefficient = coefficient %>% as.double()), 
                              aes(x = feature, y = coefficient, fill = feature) ) + 
      geom_boxplot(fill = color_boxplot) +
      coord_flip() +
      #ggforce::facet_col(vars(type)) + 
      facet_grid(vars(subgroup), scales = "free", switch = "y", space = "free_y")  +
      #ggtitle(graph_title) +
      ylab("Distribution of the significant estimate beta values (q<0.05)") +
      xlab(" ") +
      theme_light() +
      custom_ggtheme +
      scale_fill_manual(breaks = c("MSR_Donor","MSR_Acceptor"),
                        labels = c("MSR_Donor","MSR_Acceptor")) +
      theme(axis.text.y = element_text(vjust = 0.5, hjust = 1)) +
      geom_hline(yintercept = 0,linetype='dotted')
    
    
    plotTissuesZINB
    file_name <- paste0(figures_folder, 
                        "/ZINB_",MSR_type,
                        "_alltissues_phastcons", phastcons_type, 
                        "_common_introns", common,
                        "_intronsize", intron_size)
    ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                    width = 180, height = 60, units = "mm", dpi = 300)
  }

}

compare_tissues_somatic_mutations <- function(project_id1 = "SKIN",
                                              cluster_id1 = "Skin - Sun Exposed (Lower leg)",
                                              project_id2 = "SKIN",
                                              cluster_id2 = "Skin - Not Sun Exposed (Suprapubic)",
                                              stats = F) {
  project_id1 = "SKIN"
  cluster_id1 = "Skin - Sun Exposed (Lower leg)"
  project_id2 = "SKIN"
  cluster_id2 = "Skin - Not Sun Exposed (Suprapubic)"
  # 
  # project_id1 = "BRAIN"
  # cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)"
  # project_id2 = "PANCREAS"
  # cluster_id2 = "Pancreas"
  #
  # project_id1 = "BRAIN"
  # cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)"
  # project_id2 = "BLOOD"
  # cluster_id2 = "Whole Blood"
  
  tables <- c(paste0(cluster_id1, "_", project_id1), 
              paste0(cluster_id2, "_", project_id2)) #c("Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)")
  
  
  if ( !file.exists(paste0(results_folder, "/somatic_mutations_subsampled_", 
                           project_id1, "_", project_id2, ".rds")) ) {
    
    database_introns <- map_df(tables, function(table) {
      
      message(Sys.time(), " - ", table)
      # table <- tables[1]
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals, transcript_id 
                      FROM '", table, "_nevermisspliced'")
      df_introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals, transcript_id 
                      FROM '", table, "_misspliced'")
      df_introns <- rbind(df_introns, dbGetQuery(con, query) %>% as_tibble())
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_length, ref_mes5ss, ref_mes3ss, 
                  protein_coding FROM 'intron' WHERE ref_junID IN (",
                      paste(df_introns$ref_junID, collapse = ","),")")
      df_introns <- df_introns %>%
        left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = "ref_junID") %>% 
        as_tibble() 
      
      ## MERGE WITH TRANSCRIPTS
      query <- paste0("SELECT * 
                      FROM 'transcript' WHERE id IN (",
                      paste(df_introns$transcript_id, collapse = ","),")")
      df_introns <- df_introns %>%
        left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = c("transcript_id" = "id")) %>% 
        as_tibble() 
      
      
      ## MERGE WITH GENES
      query <- paste0("SELECT * 
                      FROM 'gene' WHERE id IN (",
                      paste(df_introns$gene_id, collapse = ","),")")
      df_introns <- df_introns %>%
        left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = c("gene_id" = "id")) %>% 
        as_tibble() 
      
      
      ## Add cluster_info and return
      df_introns %>%
        mutate(tissue = table) %>%
        return()
      
    } )
    
    
    ## TIDY DATAFRAME
    df_database_introns_tidy <- database_introns %>%
      mutate(tissue = tissue %>% as.factor()) %>%
      group_by(tissue) %>%
      distinct(ref_junID, .keep_all = T) %>%
      ungroup() %>%
      group_by(tissue, ref_junID) %>%
      mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
      mutate(mean_coverage = mean_coverage %>% log10()) %>%
      ungroup()  
    
    saveRDS(object = df_database_introns_tidy,
            file = paste0(results_folder, "/introns_", project_id1, "_", project_id2, ".rds"))
    
    ###########################
    ## SUBSAMPLING
    ###########################
    
    ## Choose only common intron across the two tissues
    data_combined <- df_database_introns_tidy %>%
      group_by(tissue) %>%
      distinct(ref_junID, .keep_all = T) %>%
      ungroup() %>%
      group_by(ref_junID) %>%
      mutate(n = n()) %>%
      ungroup() %>%
      #dplyr::select(tissue, mean_coverage, ref_junID, n) %>%
      filter(n == 2) %>%
      dplyr::select(-n)
    
    data_combined$tissue %>% unique
    data_combined$mean_coverage %>% summary
    data_combined %>% dplyr::count(ref_junID) %>% filter(n!=2)
    data_combined%>% dplyr::count(tissue)
    
    print(paste0(Sys.time(), " - subsampling..."))
    
    ## SUBSAMPLE BY MEAN COVERAGE
    ## Subsampling introns to control by similarity in mean read coverage
    m.out <- MatchIt::matchit(tissue ~ mean_coverage, 
                              data = data_combined, 
                              distance = data_combined$mean_coverage,
                              method = "nearest", 
                              caliper = .005, 
                              std.caliper = FALSE)
    subsample <- MatchIt::match.data(m.out)
    subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(tissue)
    
    saveRDS(object = subsample,
            file = paste0(results_folder, "/somatic_mutations_subsampled_", project_id1, "_", project_id2, ".rds"))
    
    
    df_database_introns <- df_database_introns_tidy
    subsample_introns <- subsample
    
    rm(subsample)
    rm(df_database_introns_tidy)
    gc()
    
    print(paste0(Sys.time(), " - subsampling finished!"))
    
  } else {
    
    df_database_introns <- readRDS(file = paste0(results_folder, "/introns_", project_id1, "_", project_id2, ".rds"))
    subsample_introns <- readRDS(paste0(results_folder, "/somatic_mutations_subsampled_", project_id1, "_", project_id2, ".rds"))
    
  }
  
  ##########################
  ## PLOT - COVERAGE
  ##########################
  subsample_introns %>%
    dplyr::count(tissue)
  
  df_database_introns_tidy <- df_database_introns %>%
    mutate(tissue = str_remove(string = tissue, pattern = "_SKIN"))%>%
    mutate(tissue = ifelse(tissue == "Skin - Sun Exposed (Lower leg)", "Skin - Sun Exposed", "Skin - Not Sun Exposed"))
  
  subsample_introns_tidy <- subsample_introns %>%
    mutate(tissue = str_remove(string = tissue, pattern = "_SKIN"))%>%
    mutate(tissue = ifelse(tissue == "Skin - Sun Exposed (Lower leg)", "Skin - Sun Exposed", "Skin - Not Sun Exposed"))
  
  ## VISUALISE MEAN READ COVERAGE BEFORE SUBSAMPLING
  plot_coverage_bf <- ggplot(data = df_database_introns_tidy %>% distinct(ref_junID,.keep_all = T)) +
    geom_density(mapping = aes(x = mean_coverage, fill = tissue), alpha = 0.8) +
    ggtitle("Before subsampling") +
    #scale_fill_discrete(name = "") +
    xlab("log10 mean expression level") +
    ggsci::scale_fill_npg() +
    theme_light() +
    custom_ggtheme  +
    theme(legend.position = 'top', 
          legend.spacing.x  = unit(0.3, 'cm'))
  plot_coverage_bf
  
  ## VISUALISE MEAN READ COVERAGE AFTER SUBSAMPLING
  plot_coverage_af <- ggplot(data = subsample_introns_tidy) +
    geom_density(mapping = aes(x = mean_coverage, fill = tissue), alpha = 0.8) +
    ggtitle("After subsampling") +
    #ggtitle("Mean read coverage per annotated intron across all samples\nfrom 54 GTEx v8 tissues - Subsampling performed.") +
    xlab("log10 mean expression level") +
    ggsci::scale_fill_npg() +
    theme_light() +
    custom_ggtheme
  plot_coverage_af
  
  ggpubr::ggarrange(plot_coverage_bf,
                    plot_coverage_af,
                    labels = c("a", "b"),
                    common.legend = T)
  
  file_name <- paste0(figures_folder, "/", project_id1, "_", project_id2, "_somatic_mutations_coverage")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  ################################
  ## ONLY EVALUATE COMMON INTRONS
  ################################
  
  ## TIDY THE DATA
  ## ONLY CHECK COMMON INTRONS
  
  common_introns <- subsample_introns %>%
    group_by(tissue) %>%
    distinct(ref_junID) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == 2) 
  
  df_coverage <- subsample_introns %>%
    filter(ref_junID %in% common_introns$ref_junID)
  
  df_coverage$tissue = factor(df_coverage$tissue, 
                              levels = c(paste0(cluster_id1, "_", project_id1),
                                         paste0(cluster_id2, "_", project_id2)))
  
  
  df_coverage %>% dplyr::count(tissue)
  

  ####################################
  ## STATS
  ####################################

  if (stats) {
    
    df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_D") %>% pull(MSR) %>% summary
    df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_D") %>% pull(MSR) %>% summary
    
    ## 1. Splicing noise is more frequent in introns from non-protein coding transcripts 
    ## than in introns from protein-coding transcripts even after correcting by mean read coverage

    
    ## splicing noise at the donor is more frequent in introns from  vs "Skin - Sun Exposed (Lower leg)" 
    wilcox.test(x = df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_D") %>% pull(MSR),
                y = df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_D") %>% pull(MSR),
                #alternative = "greater",
                paired = T,
                correct = T)
    
    
    ## splicing noise at the acceptor 
    
    df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_A") %>% pull(MSR) %>% summary
    df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_A") %>% pull(MSR) %>% summary
    
    wilcox.test(x = df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_A") %>% pull(MSR),
                y = df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_A") %>% pull(MSR),
                #alternative = "greater",
                paired = T,
                correct = T)
    
    

  }
  
}





## We Found that the expression levels of 107 RBPs (FDR<0.04) and 5 essential NMD genes67 (FDR<0.04) decreased with age in multiple tissues (Supplementary Figure 18a,b) (Supplementary Table 7). 
## Focusing on brain tissue alone, 40% of the 115 RBPs studied had decreased expression levels with age (FDR<0.04) (Supplementary Figure 18b).

RBP_NMD_expression_across_tissues <- function() {
  
  
  source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/scripts/27_age_effect_uncorrected_TPM_lm.R")
  source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/scripts/99_utils.R")
  
  
  ###################################
  ## LOAD THE GENE LIST
  ###################################
  
  gene_type <- "RBP"
  
  if (gene_type == "NMD") {
    gene_list <- data.frame(id = c("ENSG00000005007", "ENSG00000151461", "ENSG00000169062", "ENSG00000157106", "ENSG00000198952", "ENSG00000070366", "ENSG00000116698"),
                            name = c("UPF1", "UPF2", "UPF3", "SMG1", "SMG5", "SMG6", "SMG7"))
  } else {
    gene_list <- all_RBPs <- xlsx::read.xlsx(file = file.path(here::here(dependencies_folder, '/RBPs_subgroups.xlsx')), 
                                             header = TRUE, sheetIndex = 1) %>% as_tibble() %>% distinct(name, .keep_all = T)
  }
  
  
  ###################################
  ## LINEAR REGRESSION TO TEST IF THE
  ## COVARIATE AGE AFFECTS THE LOG10 TPM LEVELS
  ###################################
  
  if ( !file.exists(file.path(results_folder, paste0(gene_type, "_genes_age_lm.rds"))) ) {
    
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
            file = file.path(results_folder,"/", paste0(gene_type, "_genes_age_lm.rds")) ) 
    
    write_csv(x = gene_age_lm,
              file = paste0(results_folder,"/", paste0(gene_type, "_genes_age_lm.csv")))
    
  } else {
    message("Loading '", paste0(gene_type, "_genes_age_lm.rds"), "' file ...")
    gene_age_lm <- readRDS(file = file.path(results_folder, paste0(gene_type, "_genes_age_lm.rds")) ) 
  }
  
  
  tissues_data_tidy <- gene_age_lm %>%
    #mutate(Estimate = Estimate / 100) %>%
    filter(covariate == "gtex.age") 
  
  
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
  
  # write.csv(x = brain_data_tidy %>%
  #             dplyr::select(name, Estimate, q, pval) %>%
  #             as.data.frame,
  #           file = paste0(results_folder, "/_paper_review/results/brain_age_RPBs_lm.csv"))
  # 
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
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 5),
          axis.text.y = element_text(size = 5, colour = "black"),
          legend.position = "top",
          legend.text = element_text(size = 5, colour = "black"),
          legend.title = element_text(size = 5, colour = "black"),
          legend.box.margin = margin(l = -10, r = -10, b = -10, t = -5
          ))
  
  
  
  file_name <- paste0(figures_folder, "/age_all_effect_", gene_type, ".png")
  ggplot2::ggsave(filename = file_name, width = 180, height = 90, units = "mm", dpi = 300)
  
  
}


##################################
## CALLS
##################################

get_ZIP_across_tissues()
  