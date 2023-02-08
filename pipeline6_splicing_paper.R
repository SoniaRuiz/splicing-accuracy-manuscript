library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(DBI)
library(ggforce)
library(doParallel)

####################################################
## CONNECT TO THE SPLICING DATABASE ################
####################################################

# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline6_splicing_paper.R")

## CONNECT TO THE DATABASE ------------------------------

setwd("~/PROJECTS/splicing-project-recount3/")
gtf_version <- 105
main_project <- "splicing"
# database_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v",
#                         gtf_version, "/", main_project, "/", main_project, ".sqlite")
getwd()
database_path <- paste0(getwd(), "/database/v", gtf_version, "/", main_project,
                          "/", main_project, ".sqlite")


con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)

## GET FROM MASTER TABLE
query = paste0("SELECT * FROM 'master'")
df_metadata <- dbGetQuery(con, query) %>% as_tibble()

all_projects <- df_metadata$SRA_project %>% unique
all_projects %>% length() %>% print()


get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}


########################################
## FUNCTIONS - produce figures for the
## paper
########################################


## SECTION 1 ---------------------------------------------



get_database_stats <- function() {
  
  
  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'master'")
  db_metadata <- dbGetQuery(con, query) %>% as_tibble()
  
  
  ## Methods: number of samples by RIN number
  db_metadata %>% nrow() 
  
  db_metadata %>%
    filter(rin >= 8) %>% nrow()
  db_metadata %>%
    filter(rin >= 7) %>% nrow()
  db_metadata %>%
    filter(rin >= 6) %>% nrow()
  
  if ( db_metadata %>% filter(rin < 6) %>% nrow() > 1 ) {
    print("ERROR! Some of the samples considered present a RIN number lower than 6!")
    break;
  }
  
  
  db_metadata %>%
    dplyr::count(cluster) %>%
    print(n = 50)
  
  
  query <- paste0("SELECT * from 'intron'")
  db_introns <- dbGetQuery(con, query) 
  db_introns %>% distinct(ref_junID) %>% nrow()
  db_introns %>%
    dplyr::count(misspliced)
 
  query <- paste0("SELECT * from 'novel'")
  db_novel <- dbGetQuery(con, query) 
  db_novel %>% nrow() 
  
  db_novel %>% 
    dplyr::count(novel_type)
  
  ## Novel junctions exceed in X fold to annotated introns
  (db_novel %>% distinct(novel_junID) %>% nrow()) / (db_introns %>% filter(misspliced == 1) %>% distinct(ref_junID) %>% nrow())
  
  ## Percentage of mis-spliced introns
  (((db_introns %>%
      dplyr::count(misspliced) %>%
      filter(misspliced == 1) %>%
      pull(n)) * 100 )) /
  ((db_introns %>% distinct(ref_junID) %>% nrow()) %>%
    round(digits = 1))
  
  
  ## Collectively, we detected X novel junctions
  db_novel %>% distinct(novel_junID) %>% nrow()
  
  
  ## equating to 12 novel junctions per an annotated junction.
  ((novel %>% distinct(novel_junID) %>% nrow()) / (db_introns %>%
                                                  dplyr::count(misspliced) %>%
                                                  filter(misspliced == 1) %>%
                                                  pull(n))) %>%
    round()
  
  
  ## After accounting for sample number, we found that the highest numbers of unique 
  ## novel junctions were identified in X tissue with the lowest numbers in Y tissue, 
  
  
  db_metadata_tidy <- db_metadata %>%
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
      query <- paste0("SELECT COUNT(DISTINCT novel_junID) FROM '", cluster_id, "_", project_id, "_misspliced'")
      novel_junctions <- dbGetQuery(con, query) %>% as.integer()
      
      return(data.frame(cluster = cluster_id,
                        n_novel = novel_junctions))
    })
    
  })
    
  
  db_metadata_tidy %>%
    inner_join(y = df_novel_jxn_count,
               by = "cluster")
    
  
  db_metadata_tidy %>%
    inner_join(y = df_novel_jxn_count,
               by = "cluster") %>% 
    mutate(n_novel_sample = n_novel/n) %>% 
    arrange(n_novel_sample)
}

get_contamination_rates_all_tissues <- function () {

  
  #############################
  ## CONNECT TO THE DATABASE
  #############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  contamination_file <- paste0(getwd(), "/results/_paper/results/contamination_rates.rds")
  
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
                              kept_annotation = (kept_annotation %>% nrow() * 100) / db_novel_old %>% nrow(),
                              in_annotation = (in_annotation %>% nrow() * 100) / db_novel_old %>% nrow()))
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
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data = df_contamination_tidy) +
    geom_bar(mapping = aes(x = tissue, y = prop, fill = type), 
             stat = "identity", position = "identity") + 
    #coord_flip() +
    ggforce::facet_zoom(ylim = c(0,1.5), 
                        zoom.size = 3) +
    ylab("% unique novel junctions") +
    xlab("Tissue") +
    scale_fill_manual(values = c( "#666666", "#F79044FF" ),
                      breaks = c( "kept_annotation", "in_annotation" ),
                      labels = c( "novel junctions mantaining category", "novel junctions entring annotation" )) +
    #scale_fill_viridis_d(option = "A") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          #axis.text.x = element_text(colour = "black", size = "12"),
          legend.position = "top",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) %>%
    return()
  
  figures_path <- paste0(getwd(), "/results/_paper/figures/")
  dir.create(file.path(figures_path), recursive = TRUE, showWarnings = T)
  file_name <- paste0(figures_path, "/contamination_rates_all_tissues")
  
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 250, height = 140, units = "mm", dpi = 300)
  
  
  
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

get_contamination_rates_FCTX <- function () {
  
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
                          in_annotation = (in_annotation %>% nrow() * 100) / db_novel_old %>% nrow(),
                          out_annotation = (out_annotation %>% nrow() * 100) / db_introns_old %>% nrow()))
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
  
  
  df_contamination$contamination_rates = factor(df_contamination$contamination_rates, 
                                                levels = c("v76\n(2014)",
                                                           "v81\n(2015)",
                                                           "v90\n(2017)",
                                                           "v97\n(2019)",
                                                           "v104\n(2021)"))


  ################################
  ## LINES PLOT
  ################################
  
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data=df_contamination, aes(x = contamination_rates, 
                                    y = in_annotation,
                                    group=1)) +
    geom_line(linetype = "dashed", col = "red") +
    geom_point() +
    xlab(NULL) +
    ylab("contamination rates (%)\nas compared to Ensembl v105") +
    #scale_fill_viridis_d(option = "E")  +
    #ggtitle(paste0("Brain - Frontal Cortex (BA9)\nContamination rates compared with Ensembl v105")) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          legend.position = NULL,
          axis.text.y = element_text(vjust = 0.5,hjust = 1)) +
    guides(fill = "none")
  
  
  file_name <- paste0(getwd(), "/results/_paper/figures/contamination_rates_FCTX_lineplot")
  ggplot2::ggsave(paste0(file_name, ".svg"),  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"),  width = 70, height = 70, units = "mm", dpi = 300)
  
  
  ################################
  ## BAR PLOT
  ################################
  
  ggplot(data = df_contamination) +
    geom_bar(aes(x = contamination_rates, 
                 y = in_annotation ), 
             stat = "identity") + 
    coord_flip() +
    xlab(NULL) +
    ylab("contamination rates (%) as compared to Ensembl v105") +
    #scale_fill_viridis_d(option = "E")  +
    #ggtitle(paste0("Brain - Frontal Cortex (BA9)\nContamination rates compared with Ensembl v105")) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          legend.position = NULL,
          axis.text.y = element_text(#angle = 50, 
            vjust = 0.5,
            hjust = 1)) +
    guides(fill = "none")
  
  file_name <- paste0(getwd(), "/results/_paper/figures/contamination_rates_FCTX")
  ggplot2::ggsave(paste0(file_name, ".svg"),  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"),  width = 183, height = 183, units = "mm", dpi = 300)
  
}

get_unique_donor_acceptor_jxn <- function() {
  
  ##################################################################################################
  ## Proportion of each type of unique junctions per GTEx tissue.
  ##################################################################################################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  file_name <- paste0(getwd(), "/results/_paper/results/unique_donor_acceptor_jxn.rds")

 
  if ( !file.exists(file_name) )  {
    
    df_proportions <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[5]
      
      print(paste0(Sys.time(), " - ", project_id))
      
      all_clusters <- df_metadata %>%
        filter(SRA_project == project_id) %>%
        distinct(cluster) %>%
        pull()
      
      print(all_clusters)
      
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
        novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
        query <- paste0("SELECT DISTINCT novel_junID, novel_type FROM 'novel' WHERE novel_junID IN (", 
                        paste(novel_junctions$novel_junID, collapse = ","), ")")
        novel_junctions <- novel_junctions %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "novel_junID")
        
        
        introns %>% head()
        novel_junctions %>% head()
        
        
        ###############################
        ## Calculate the proportions
        ###############################
        
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
    figures_path <- paste0(getwd(), "/results/_paper/results/")
    dir.create(file.path(figures_path), recursive = TRUE, showWarnings = T)
    
    saveRDS(object = df_proportions, file = file_name)
    write.csv(x = df_proportions,
              file = paste0(getwd(), "/results/_paper/results/unique_donor_acceptor_jxn.csv"),
              row.names = F)
  } else {
    df_proportions <- read.csv(file = paste0(getwd(), "/results/_paper/results/unique_donor_acceptor_jxn.csv"))
    
    df_proportions %>%
      arrange(desc(annotated_prop))
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
    geom_point() +
    geom_line( aes(group = tissue) )  +
    theme_light() +
    ylab("% unique junctions") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_x_discrete( breaks = c( "acceptor", "donor"),# "annotated_intron"),
                        labels = c( "novel acceptor", "novel donor")) +
    scale_fill_manual(values = c( "#35B779FF", "#8d03b0"),
                      breaks = c( "donor", "acceptor"),
                      labels = c( "novel donor", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))+
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.text.x = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          legend.text = element_text(colour = "black", size = "12"),
          plot.caption = element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") 
  
  junx_violin_plot
  
  ## Save the figure 3
  file_name <- paste0(getwd(), "/results/_paper/figures/percent_unique_donor_acceptor_violin")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), 
                  width = 183, height = 130, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                  width = 100, height = 100, units = "mm", dpi = 300)
  
  ######################
  ## BAR PLOT
  ######################
  
  df_proportions <- df_proportions %>%
    dplyr::select(tissue, 
                  donor = donor_prop, 
                  acceptor = acceptor_prop, 
                  annotated_intron = annotated_prop) %>%
    tidyr::gather(key = "type", value = "prop", -tissue ) 
  
  df_proportions %>% head()

  
  ## First, order by junction type
  df_proportions$type <- factor(df_proportions$type, 
                             levels = c( "acceptor", 
                                         "annotated_intron",
                                         "donor"))
  
  df_proportions_final <- df_proportions %>% 
    ungroup() %>%
    arrange(type , prop) %>%
    mutate(tissue = fct_inorder(tissue))
  
  
  # colours <- ifelse(str_detect(string = as.factor(df_proportions_final$tissue), 
  #                              pattern = "Brain"), "red", "black")
  
  
  ggplot(df_proportions_final %>% filter(type %in% c("donor", "acceptor")), 
         aes(x = tissue, y = prop * 100, 
             group = type, fill = type)) +
    geom_col(position = position_dodge()) +
    #geom_text(aes(label = round(x = prop, digits = 2)), 
    #          position = position_stack(vjust = 0.5, reverse = TRUE), colour = "white", size = 1.5) +
    coord_flip() +
    theme_light() +
    ylab("unique novel junctions across samples (%)") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          #axis.text.x = element_text(color = colours),
          axis.text.y = element_text(#color = colours,
                                     #angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          axis.title.x = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "9"),
          plot.title = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    scale_fill_manual(values = c( "#0D0887FF", "#EF7F4FFF"),
                      breaks = c( "acceptor", "donor"),# "annotated_intron"),
                      labels = c( "novel acceptor", "novel donor")) + #, "annotated intron")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  ## Save the figure 3
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/percent_unique_donor_acceptor")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), 
                  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                  width = 183, height = 183, units = "mm", dpi = 300)
}

get_unique_donor_acceptor_reads <- function() {
  
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  file_name <- paste0(getwd(), "/results/_paper/results/unique_donor_acceptor_reads.rds")
  
  if ( !file.exists(file_name) )  {
     
    df_mean_counts <- map_df(all_projects, function(project_id) {
      
      #############################
      ## GET THE CLUSTERS
      #############################
      print(paste0(Sys.time(), " - ", project_id))
      
      
      all_clusters <- df_metadata %>%
        filter(SRA_project == project_id) %>%
        distinct(cluster) %>%
        pull()
      
      print(all_clusters)
      
      
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
        
        query <- paste0("SELECT ref_junID, ref_sum_counts, ref_n_individuals FROM '", 
                        cluster_id, "_", project_id, "_nevermisspliced'")
        introns <- dbGetQuery(con, query) %>% as_tibble()
        query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals FROM '", 
                        cluster_id, "_", project_id, "_misspliced'")
        introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
        
        
        
        ###########################
        ## GET THE NOVEL JUNCTIONS
        ###########################
        
        query <- paste0("SELECT novel_junID, novel_sum_counts, novel_n_individuals FROM '", cluster_id, "_", project_id, "_misspliced'")
        novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
        query <- paste0("SELECT novel_junID, novel_type FROM 'novel' WHERE novel_junID IN (", 
                        paste(novel_junctions$novel_junID, collapse = ","), ")")
        novel_junctions <- novel_junctions %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
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
  

  ######################
  ## VIOLIN PLOT
  ######################
  
  df_mean_counts_violin <- df_mean_counts %>%
    filter(type != "annotated") 
  
  df_mean_counts_violin$type = factor(df_mean_counts_violin$type, 
                                      levels = c("donor", "acceptor"))
  
  
  reads_violin_plot <- ggplot(df_mean_counts_violin, aes(type, prop_counts, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point() +
    geom_line( aes(group = tissue) )  +
    theme_light() +
    ylab("% read counts") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_x_discrete( breaks = c( "acceptor", "donor"),# "annotated_intron"),
                      labels = c( "novel acceptor", "novel donor")) +
    scale_fill_manual(values = c( "#35B779FF", "#8d03b0"),
                      breaks = c( "donor", "acceptor"),
                      labels = c( "novel donor", "novel acceptor")) + 
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))+
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.text.x = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          legend.text = element_text(colour = "black", size = "12"),
          plot.caption = element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") 
  
  reads_violin_plot
  file_name <- paste0(getwd(), "/results/_paper/figures/unique_donor_acceptor_reads_violin")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 150, height = 150, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 100, height = 100, units = "mm", dpi = 300)
  
  
  gga
  junx_violin_plot
  
  ggpubr::ggarrange(junx_violin_plot,
                    reads_violin_plot,
                    common.legend = T,
                    labels = c("c", "d"),
                    align = "v",
                    ncol = 2,
                    nrow = 1)
  
  file_name <- paste0(getwd(), "/results/_paper/figures/donor_acceptor_junx_reads_violin")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 150, height = 150, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 200, height = 100, units = "mm", dpi = 300)
  
  #####################
  ## BAR PLOT
  #####################
  
  ## PREPARE DATA FOR PLOTTING
  ## Spread
  df_mean_counts <- df_mean_counts %>% 
    dplyr::select( -sum_counts) %>%
    #filter(tissue %in% all_clusters) %>%
    tidyr::spread(key = "type", value = "prop_counts")

  
  
  ## Gather
  df_mean_counts_tidy <- df_mean_counts %>% 
    #dplyr::select(annotated = annotated_p, acceptor = acceptor_p, donor = donor_p, tissue) %>%
    tidyr::gather(key = "type", value = "mean", -tissue, - samples )
  
  
  df_mean_counts_tidy$type = factor(df_mean_counts_tidy$type, 
                                    levels = c("annotated", "acceptor", "donor"))
  
  df_mean_counts_tidy = df_mean_counts_tidy %>% 
    ungroup() %>%
    arrange(type , desc(mean)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  #colours <- ifelse(str_detect(string = as.factor(df_mean_counts_tidy$tissue), pattern = "Brain"), "red", "black")
  
  
  ggplot(data = df_mean_counts_tidy %>% filter(type %in% c("acceptor", "donor"))) +
    geom_col(mapping = aes(x = tissue, y = mean, fill = type), 
             position = position_dodge(), 
             alpha = 0.9) +
    #coord_cartesian( ylim=c(0,1) ) +
    coord_flip() +
    xlab("") +
    ylab("read counts across samples (%)") +
    theme_light() +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = c(0,0.5,0.75,1)) +
    
    scale_fill_manual(values =  c("#A41F9AFF","#0D0887FF","#EF7F4FFF"),#c( "#440154FF", "#35B779FF", "#39558CFF"),
                      breaks = c("annotated", "acceptor", "donor"),
                      labels = c("annotated introns", "novel acceptor", "novel donor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          axis.text.y = element_text(#color = colours,
                                     #angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) 
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/unique_donor_acceptor_reads")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  
}

get_mean_read_count_FTCX <- function() {
  
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  tables <- dbListTables(con)
  
  #tables <- c("Brain - Frontal Cortex (BA9)_BRAIN_misspliced",
  #            "Brain - Frontal Cortex (BA9)_BRAIN_nevermisspliced")
  
  SRA_projects <- df_metadata$SRA_project %>% unique()
  

  
  query <- paste0("SELECT novel_junID, novel_type from 'novel'")
  db_novel <- dbGetQuery(con, query) %>% as_tibble()
  
  
  df_read_count <- map_df(SRA_projects, function(project_id) {
    
    # project_id <- SRA_projects[1]
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    
    map_df(all_clusters, function(cluster_id) {
      
      # cluster_id <- all_clusters[1]
      
      query <- paste0("SELECT ref_sum_counts FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      db_never <- dbGetQuery(con, query) %>% as_tibble()
      db_never$ref_sum_counts
      
      query <- paste0("SELECT ref_sum_counts, novel_sum_counts, novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
      db_misspliced <- dbGetQuery(con, query) %>% as_tibble()
      db_misspliced$ref_sum_counts
      
      novel_donor_sum_count <- db_misspliced  %>%
        filter(novel_junID %in% (db_novel %>% 
                                   filter(novel_type == "novel_donor") %>%
                                   pull(novel_junID))) %>%
        pull(novel_sum_counts)
      
      novel_acceptor_sum_count <- db_misspliced  %>%
        filter(novel_junID %in% (db_novel %>% 
                                   filter(novel_type == "novel_acceptor") %>%
                                   pull(novel_junID))) %>%
        pull(novel_sum_counts)
      
      
      return( data.frame(annotated_median_count = c(db_never$ref_sum_counts,
                                                    db_misspliced$ref_sum_counts) %>% median(),
                         annotated_max_count = c(db_never$ref_sum_counts,
                                                 db_misspliced$ref_sum_counts) %>% max(),
                         annotated_min_count = c(db_never$ref_sum_counts,
                                                 db_misspliced$ref_sum_counts) %>% min(),
                         novel_donor_median_count = novel_donor_sum_count %>% median(),
                         novel_donor_max_count = novel_donor_sum_count %>% max(),
                         novel_donor_min_count = novel_donor_sum_count %>% min(),
                         novel_acceptor_median_count = novel_acceptor_sum_count %>% median(),
                         novel_acceptor_max_count = novel_acceptor_sum_count %>% max(),
                         novel_acceptor_min_count = novel_acceptor_sum_count %>% min(),
                         tissue = cluster_id) )
      
    })
    
   
  })
  
  
  
  write.csv(x = df_read_count %>%
              relocate(tissue),
            file = paste0(getwd(), "/results/_paper/05_supplementary_data_read_count.csv"),
            row.names = F)

 
  
  ## Focusing on frontal cortex, we found that this equated to a median read count 
  ## of 2,700 for annotated junctions with novel donor or acceptor events having a 
  ## median read count of only 2 in both cases. 
  
  df_read_count %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)")
  
}

get_maxentscan_score <- function() {
  
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  
  file_name <- paste0(getwd(), "/results/_paper/results/df_mes.rds")
  
  if ( !file.exists(file_name) )  {
    
    df_mes <- map_df(all_projects, function(project_id) {
      
      #############################
      ## GET THE CLUSTERS
      #############################
      
      print(paste0(Sys.time(), " - ", project_id))
      
      # clusters_used <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
      #                                        project_id, "/raw_data/all_clusters_used.rds"))
      
      all_clusters <- df_metadata %>%
        filter(SRA_project == project_id) %>%
        distinct(cluster) %>%
        pull()
      
      print(all_clusters)
      
      map_df(all_clusters, function(cluster_id) {
        
        ## Print the tissue
        print(paste0(Sys.time(), " - ", cluster_id))
        
        samples <- df_metadata %>%
          dplyr::count(cluster) %>%
          filter(cluster == cluster_id) %>% 
          pull(n)
        
        
        ####################
        ## GET THE INTRONS
        ####################
        
        #query <- paste0("SELECT ref_junID, ref_type FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
        #introns <- dbGetQuery(con, query) %>% as_tibble()
        
        query <- paste0("SELECT DISTINCT ref_junID, ref_type FROM '", cluster_id, "_", project_id, "_misspliced'")
        introns <- dbGetQuery(con, query) %>% as_tibble() #rbind(introns, dbGetQuery(con, query) %>% as_tibble())
        
        query <- paste0("SELECT ref_junID, ref_ss5score, ref_ss3score FROM 'intron' WHERE ref_junID IN (",
                        paste(introns$ref_junID, collapse = ","),")")
        introns <- introns %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "ref_junID") %>% 
          as_tibble() 
        
        
        ###########################
        ## GET THE NOVEL JUNCTIONS
        ###########################
        
        query <- paste0("SELECT ref_junID, novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
        novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
        
        query <- paste0("SELECT novel_junID, novel_type, novel_ss5score, novel_ss3score FROM 'novel' WHERE novel_junID IN (", 
                        paste(novel_junctions$novel_junID, collapse = ","), ")")
        novel_junctions <- novel_junctions %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "novel_junID") %>% 
          as_tibble() 
        
        
        ###########################
        ## MERGE THE DATA AND RETURN
        ###########################
        
        df_merged <- introns %>% 
          dplyr::select(ref_junID, ref_ss5score, ref_ss3score, ref_type) %>%
          left_join(y = novel_junctions %>% 
                      dplyr::select(ref_junID, novel_junID, 
                                    novel_ss5score, novel_ss3score, 
                                    novel_type),
                    by = "ref_junID")  %>%
          mutate(diff_ss5score = ref_ss5score - novel_ss5score,
                 diff_ss3score = ref_ss3score - novel_ss3score,
                 tissue = cluster_id)
        
        df_merged %>% return()
        
      })
    })
    
    saveRDS(object = df_mes %>%
              filter(ref_type != "never") %>%
              distinct(novel_junID, .keep_all = T) %>% 
              as_tibble(), file = file_name)
    
  } else {
    df_mes <- readRDS(file = file_name) %>% as_tibble()
  }
  
  
  df_mes <- df_mes %>%
    filter(ref_type != "never") %>%
    distinct(novel_junID, .keep_all = T)
  
  
  
  #############################
  ## COMBINED FIGURES 
  ## MAXENTSCAN score & DISTANCES
  #############################
  
  
  ## ss5score ----------------------------------------------------------
  
  df_5ss <- df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron = ref_ss5score, novel_donor = novel_ss5score) %>%
    gather(key = "junction_type", value = "ss5score")
  
  ss5plot <- ggplot(df_5ss, aes(ss5score, fill = junction_type)) +
    geom_density(alpha = 0.8) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    scale_fill_manual(values = c("#A41F9AFF", "#EF7F4FFF"), 
                      breaks=c("intron", "novel_donor"),
                      labels=c("intron  ", "novel donor")) +
    theme(axis.line = element_line(colour = "black" ),
          axis.text = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          legend.position = "top") +
    xlab("5' splice site MES score") +
    guides(fill = guide_legend(title = element_blank(),
                               ncol = 2, nrow = 1))
  
  ss5plot
  file_name <- paste0(getwd(), "/results/_paper/figures/MES_ss5plot")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ## ss3score -----------------------------------------------------------
  
  df_3ss <- df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron = ref_ss3score, novel_acceptor = novel_ss3score) %>%
    gather(key = "junction_type", value = "ss3score")
  
  ss3plot <- ggplot(df_3ss, 
                    aes(ss3score, fill = junction_type)) +
    geom_density(alpha = 0.8) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    scale_fill_manual(values = c("#A41F9AFF", "#0D0887FF"), 
                      breaks = c("intron", "novel_acceptor"),
                      labels = c("intron  ", "novel acceptor")) +
    theme(axis.line = element_line(colour = "black" ),
          axis.text = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          legend.position = "top") +
    xlab("3' splice site MES score") +
    guides(fill = guide_legend(title = element_blank(), ncol = 2, nrow = 1))
  
  
  ss3plot
  file_name <- paste0(getwd(), "/results/_paper/figures/MES_ss3plot")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ## COMBO
  
  ggpubr::ggarrange( ss5plot, 
                     ss3plot, 
                     labels = c("a", "b"),
                     ncol = 2, nrow = 1)
  file_name <- paste0(getwd(), "/results/_paper/figures/MES_combo_supplementary")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
  ## Delta MES -----------------------------------------------------------
  
  
  df_delta_5ss <- df_mes %>%
    filter(novel_type == "novel_donor") %>% 
    dplyr::select(diff_ss5score) %>%
    tidyr::gather(key = "type", value = "MES") 
  
  
  ## Plot
  deltaplot5ss <- ggplot(data = df_delta_5ss) +
    geom_density(mapping = aes(x = MES,
                               fill = type)) +
    geom_vline(xintercept = 0) +
    xlab("Delta MaxEntScan 5'ss score") +
    xlim(c(-40, 65)) +
    ylim(c(0, 0.11)) +
    theme_light() +
    scale_color_viridis_d() +
    scale_fill_manual(values =  c("#35B779FF"),
                      breaks = c("diff_ss5score"),
                      labels = c("Delta MES 5'ss")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "11"),
          axis.title = element_text(colour = "black", size = "11"),
          legend.text = element_text(size = "11"),
          legend.title = element_text(size = "11"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  deltaplot5ss
  file_name <- paste0(getwd(), "/results/_paper/figures/MES_delta_5ss")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 100, height = 100, units = "mm", dpi = 300)
  
 
  df_delta_3ss <- df_mes %>%
    filter(novel_type == "novel_acceptor") %>% 
    dplyr::select(diff_ss3score) %>%
    tidyr::gather(key = "type", value = "MES") 
  
  deltaplot3ss <- ggplot(data = df_delta_3ss) +
    geom_density(mapping = aes(x = MES, fill = type)) +
    geom_vline(xintercept = 0) +
    xlab("Delta MaxEntScan 3'ss score") +
    xlim(c(-40, 65)) +
    ylim(c(0, 0.11)) +
    theme_light() +
    scale_color_viridis_d() +
    scale_fill_manual(values =  c("#64037d"),
                      breaks = c("diff_ss3score"),
                      labels = c("Delta MES 3'ss")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "11"),
          axis.title = element_text(colour = "black", size = "11"),
          legend.text = element_text(size = "11"),
          legend.title = element_text(size = "11"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  deltaplot3ss
  
  file_name <- paste0(getwd(), "/results/_paper/figures/MES_delta_3ss")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 100, height = 100, units = "mm", dpi = 300)
  
  ##############
  ## COMBO
  ##############
  
  ggpubr::ggarrange( deltaplot5ss, 
                     deltaplot3ss, 
                     labels = c("a", "b"),
                     ncol = 2, 
                     nrow = 1)
  file_name <- paste0(getwd(), "/results/_paper/figures/MES_delta_combo")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
  ##############
  ## STATS
  ##############
  
  ## As would be expected, we found that the majority of novel 5’ and 3’ splice sites were 
  ## weaker than the corresponding annotated site with 82.5% of novel 5’ 
 
  
  query <- paste0("SELECT * from 'novel'")
  db_novel <- dbGetQuery(con, query) 
  
 
  ((df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron_MES = ref_ss5score, novel_donor_MES = novel_ss5score) %>%
    mutate(MES_diff = intron_MES - novel_donor_MES) %>%
    pull(MES_diff) %>% 
    sign() %>% 
    table() %>% as.data.frame() %>%
    filter(. != -1) %>%
    pull(Freq) %>%
    sum) * 100) / ( db_novel %>% 
                      dplyr::count(novel_type) %>%
                      filter(novel_type == "novel_donor") %>%
                      pull(n) )
 


  
  ## and 85.2% of novel 3’ sites having positive delta MES scores
  
  ((df_mes %>% 
      filter(novel_type == "novel_acceptor") %>%
      dplyr::select(intron_MES = ref_ss5score, novel_acceptor_MES = novel_ss3score) %>%
      mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
      pull(MES_diff) %>% 
      sign() %>% 
      table() %>% as.data.frame() %>%
      filter(. != -1) %>%
      pull(Freq) %>%
      sum) * 100) / ( db_novel %>% 
                        dplyr::count(novel_type) %>%
                        filter(novel_type == "novel_acceptor") %>%
                        pull(n) )
  
  
  
  ## 
  
  ## Interestingly, this analysis demonstrated that novel 5’ splice sites had a modal delta value which was very close to zero
  
  df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron_MES = ref_ss5score, novel_donor_MES = novel_ss5score) %>%
    mutate(MES_diff = intron_MES - novel_donor_MES) %>%
    pull(MES_diff) %>%
    get_mode()
  
  df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron_MES = ref_ss5score, novel_donor_MES = novel_ss5score) %>%
    mutate(MES_diff = intron_MES - novel_donor_MES) %>%
    pull(MES_diff) %>%
    summary()
  
  ((df_mes %>% 
      filter(novel_type == "novel_donor") %>%
      dplyr::select(intron_MES = ref_ss5score, novel_donor_MES = novel_ss5score) %>%
      mutate(MES_diff = intron_MES - novel_donor_MES) %>%
      pull(MES_diff) %>% 
      sign() %>% 
      table() %>% 
      as.data.frame() %>%
      filter(. == 0) %>%
      pull(Freq)) * 100) / ( db_novel %>% 
                        dplyr::count(novel_type) %>%
                        filter(novel_type == "novel_donor") %>%
                        pull(n) )
  
  
  ## The modal delta value for novel 3’ splice sites was higher 
  df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron_MES = ref_ss5score, novel_acceptor_MES = novel_ss3score) %>%
    mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
    pull(MES_diff) %>%
    get_mode()
  
  df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron_MES = ref_ss5score, novel_acceptor_MES = novel_ss3score) %>%
    mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
    pull(MES_diff) %>%
    summary()
  
  
  ((df_mes %>% 
      filter(novel_type == "novel_acceptor") %>%
      dplyr::select(intron_MES = ref_ss5score, novel_acceptor_MES = novel_ss3score) %>%
      mutate(MES_diff = intron_MES - novel_acceptor_MES) %>%
      pull(MES_diff) %>% 
      sign() %>% 
      table() %>% 
      as.data.frame() %>%
      filter(. == 0) %>%
      pull(Freq)) * 100) / ( db_novel %>% 
                               dplyr::count(novel_type) %>%
                               filter(novel_type == "novel_acceptor") %>%
                               pull(n) )
  
  df_delta_3ss$MES %>% min
  df_delta_3ss$MES %>% max
  
}

get_modulo_fctx_PC <- function() {
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  reference_gtf_version <- 105
  

  
  query <- paste0("SELECT novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (", paste(introns$novel_junID, collapse = ","),")")
  introns <- introns %>%
    inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
               by = "novel_junID") %>% as_tibble() 
  
  
  
  df_novel_tidy <- introns %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  
  
  ###############################
  ## PROTEIN CODING
  ###############################
  
  
  query <- paste0("SELECT ref_junID, protein_coding, lncRNA FROM 'intron'")
  df_master_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  df_novel_tidy <- introns %>%
    inner_join(df_master_introns,
               by = "ref_junID") %>%
    filter(protein_coding == 0) %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))

  
  ###############################
  ## GET THE MODULO
  ###############################
  
  df_novel_tidy %>%
    mutate(modulo3 = abs(distance %% 3)) %>%
    dplyr::count(modulo3) %>%
    mutate(n_prop = n / (df_novel_tidy$novel_junID %>% unique() %>% length()))
  
  
  
  
}

get_modulo_basic_multiple_tissue <- function() {
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  file_path <- paste0(getwd(), "/results/_paper/results/df_modulo_basic_tissues_100bpfilter.rds")
  
  if ( !file.exists(file_path) )  {
    
    df_modulo_tissues <- map_df(all_projects, function(project_id) {
      
      print(paste0(Sys.time(), " - ", project_id))
      
      all_clusters <- df_metadata %>%
        filter(SRA_project == project_id) %>%
        distinct(cluster) %>%
        pull()
      
      map_df(all_clusters, function(cluster_id) {
        
        # cluster_id <- all_clusters[1]
        
        ## Print the tissue
        print(paste0(Sys.time(), " - ", cluster_id))
        
        
        ## Get the novel junctions from the current tissue
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
                        AND MANE = 1")
        introns <- introns %>%
          inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "ref_junID") %>% 
          as_tibble() 
        introns$MANE %>% unique
        
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
  
  library(ggridges)
  
  df_modulo_tissues <- df_modulo_tissues %>%
    mutate(freq = freq * 100)
  
  ggplot(df_modulo_tissues, aes(x = freq, y = modulo)) +
    geom_density_ridges_gradient() +
    ylab("Modulo3 of the distance") +
    xlab("% of novel junctions") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          axis.text.y = element_text(size = "10", vjust = 1, hjust = 1),
          legend.text = element_text(colour = "black",size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          plot.caption = element_text(colour = "black",size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "none") +
    scale_x_continuous(breaks = c(30,35,40),
                       labels = c(30,35,40)) +
    scale_y_discrete(expand = c(0,0,1,0))
  
  
  file_name <- paste0(getwd(), "/results/_paper/figures/modulo3_alltissues_density")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 80, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 80, units = "mm", dpi = 300)
  
  
  ################
  ## BAR PLOT
  ################

  df_modulo_tissues <- df_modulo_tissues %>% 
    ungroup() %>%
    arrange(modulo, freq) %>%
    mutate(tissue = fct_inorder(tissue))
  
  # colours <- ifelse(str_detect(string = as.factor(df_modulo_tissues$tissue), 
  #                              pattern = "Brain"), "red", "black")
  
  plot_all_tissues <- ggplot(data = df_modulo_tissues, 
                             aes(x = factor(tissue), y = freq*100, fill = factor(modulo))) + 
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip()+
    scale_fill_viridis_d() +
    ylab("% of novel junctions") +
    xlab("Tissue") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          axis.text.y = element_text(#colour = colours,
                                     size = "10",
                                     #angle = 75, 
                                     vjust = 1,
                                     hjust = 1),
          legend.text = element_text(colour = "black",size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          strip.text.y = element_blank(),
          plot.caption = element_text(colour = "black",size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo 3: ",
                               ncol = 3, 
                               nrow = 1))
  
  
  plot_all_tissues
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/modulo3_alltissues")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
 
  
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

get_modulo_complete_multiple_tissue <- function() {
  
  #############################
  ## Supplementary Figure 4.1 
  ## MODULO
  #############################
  
  if ( !file.exists( file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                   main_project,"/results/df_modulo_tissues_100bpfilter.rds") ) ) {
    

    query <- paste0("SELECT * FROM 'master'")
    df_metadata <- dbGetQuery(con, query) 
    
    
    all_projects <- df_metadata$SRA_project %>% unique()
    
    
    df_modulo_tissues <- map_df(all_projects, function(project_id) {
      
      print(paste0(Sys.time(), " - ", project_id))
      
      all_clusters <- df_metadata %>%
        filter(SRA_project == project_id) %>%
        distinct(cluster) %>%
        pull()
      
      map_df(all_clusters, function(cluster_id) {
        
        # cluster_id <- all_clusters[1]
        print(cluster_id)
        
        ## Load IDB
        
        ## Print the tissue
        print(paste0(Sys.time(), " - ", cluster_id))
        
        samples <- df_metadata %>%
          dplyr::count(cluster) %>%
          filter(cluster == cluster_id) %>% 
          pull(n)
        
        
        query <- paste0("SELECT novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
        introns <- dbGetQuery(con, query) %>% as_tibble()
        query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                        paste(introns$novel_junID, collapse = ","),")")
        introns <- merge(x = introns,
                         y = dbGetQuery(con, query) %>% as_tibble(),
                         by = "novel_junID",
                         all.x = T) %>% 
          as_tibble() 
        query <- paste0("SELECT ref_junID, protein_coding FROM 'intron' WHERE ref_junID IN (",
                        paste(introns$ref_junID, collapse = ","),") AND MANE = 1")
        introns <- merge(x = introns,
                         y = dbGetQuery(con, query) %>% as_tibble(),
                         by = "ref_junID",
                         all.x = T) %>% 
          as_tibble() 
        
        
        df_novel_tidy <- introns %>%
          filter(abs(distance) <= 100) %>%
          mutate(novel_type = str_replace(string = novel_type,
                                          pattern = "_",
                                          replacement = " ")) %>%
          mutate(type_p = ifelse(distance < 0, paste0(novel_type," intron"), paste0(novel_type," exon"))) %>% 
          mutate(modulo = abs(distance) %% 3)
        
        df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                          levels = c("novel donor", "novel acceptor"))
        df_novel_tidy$type_p = factor(df_novel_tidy$type_p, 
                                      levels = c("novel acceptor intron", 
                                                 "novel acceptor exon",
                                                 "novel donor intron", 
                                                 "novel donor exon"))
        
        
        ## Number of modulo 0 unique introns
        modulo0_donor_exon <- (df_novel_tidy %>%
                                 filter(modulo == 0, type_p == "novel donor exon") %>% nrow()) * 100 / 
          (df_novel_tidy %>%
             filter(type_p == "novel donor exon") %>% 
             nrow())
        modulo0_donor_intron <- (df_novel_tidy %>%
                                   filter(modulo == 0, type_p == "novel donor intron") %>% nrow()) * 100 / 
          (df_novel_tidy %>%
             filter(type_p == "novel donor intron") %>% 
             nrow())
        modulo0_aceptor_exon <- (df_novel_tidy %>%
                                   filter(modulo == 0, type_p == "novel acceptor exon") %>% nrow()) * 100 /
          (df_novel_tidy %>%
             filter(type_p == "novel acceptor exon") %>% 
             nrow())
        modulo0_aceptor_intron <- (df_novel_tidy %>%
                                     filter(modulo == 0, type_p == "novel acceptor intron") %>% nrow()) * 100 / 
          (df_novel_tidy %>%
             filter(type_p == "novel acceptor intron") %>% 
             nrow())
        
        ## Number of modulo 1
        modulo1_donor_exon <- (df_novel_tidy %>%
                                 filter(modulo == 1, type_p == "novel donor exon") %>% nrow()) * 100 / 
          (df_novel_tidy %>% 
             filter(type_p == "novel donor exon") %>% 
             nrow())
        modulo1_donor_intron <- (df_novel_tidy %>%
                                   filter(modulo == 1, type_p == "novel donor intron") %>% nrow()) * 100 / 
          (df_novel_tidy %>% 
             filter(type_p == "novel donor intron") %>% 
             nrow())
        modulo1_aceptor_exon <- (df_novel_tidy %>%
                                   filter(modulo == 1, type_p == "novel acceptor exon") %>% nrow()) * 100 / 
          (df_novel_tidy %>%
             filter(type_p == "novel acceptor exon") %>% 
             nrow())
        modulo1_aceptor_intron <- (df_novel_tidy %>%
                                     filter(modulo == 1, type_p == "novel acceptor intron") %>% nrow()) * 100 / 
          (df_novel_tidy %>%
             filter(type_p == "novel acceptor intron") %>% 
             nrow())
        
        ## Number of modulo 2
        modulo2_donor_exon <- (df_novel_tidy %>%
                                 filter(modulo == 2, type_p == "novel donor exon") %>% nrow()) * 100 / 
          (df_novel_tidy %>% 
             filter(type_p == "novel donor exon") %>% 
             nrow())
        modulo2_donor_intron <- (df_novel_tidy %>%
                                   filter(modulo == 2, type_p == "novel donor intron") %>% nrow()) * 100 / 
          (df_novel_tidy %>% 
             filter(type_p == "novel donor intron") %>% 
             nrow())
        modulo2_aceptor_exon <- (df_novel_tidy %>%
                                   filter(modulo == 2, type_p == "novel acceptor exon") %>% nrow()) * 100 / 
          (df_novel_tidy %>% 
             filter(type_p == "novel acceptor exon") %>% 
             nrow())
        modulo2_aceptor_intron <- (df_novel_tidy %>%
                                     filter(modulo == 2, type_p == "novel acceptor intron") %>% nrow()) * 100 / 
          (df_novel_tidy %>% 
             filter(type_p == "novel acceptor intron") %>% 
             nrow())
        
        return(data.frame(tissue = cluster_id,
                          modulo = c(modulo0_donor_exon, 
                                     modulo0_donor_intron, 
                                     modulo0_aceptor_exon, 
                                     modulo0_aceptor_intron, 
                                     modulo1_donor_exon, 
                                     modulo1_donor_intron, 
                                     modulo1_aceptor_exon, 
                                     modulo1_aceptor_intron, 
                                     modulo2_donor_exon, 
                                     modulo2_donor_intron, 
                                     modulo2_aceptor_exon, 
                                     modulo2_aceptor_intron),
                          novel_type = c("donor exon", 
                                         "donor intron", 
                                         "aceptor exon", 
                                         "aceptor intron", 
                                         "donor exon", 
                                         "donor intron", 
                                         "aceptor exon", 
                                         "aceptor intron", 
                                         "donor exon", 
                                         "donor intron", 
                                         "aceptor exon", 
                                         "aceptor intron"),
                          type = c("modulo0","modulo0","modulo0","modulo0",
                                   "modulo1","modulo1","modulo1","modulo1",
                                   "modulo2","modulo2","modulo2","modulo2")))
        
        
      })
    })
    saveRDS(df_modulo_tissues, 
            file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                          main_project,"/results/df_modulo_tissues_100bpfilter.rds"))
    
  } else {
    df_modulo_tissues <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                               main_project,"/results/df_modulo_tissues_100bpfilter.rds")) %>%
      as_tibble()
  }
  
  #######################
  ## DENSITY PLOT
  #######################
  
  ggplot(df_modulo_tissues, aes(x = modulo , y = type , fill=type )) +
    geom_density_ridges_gradient() +
    facet_grid(vars(novel_type)) +
    scale_fill_viridis_d(name = "Freq.Mod3", option = "D") +
    #labs(title = 'Modulo3 of the distances to annotated splice sites')+
    ylab("Modulo3") +
    xlab("Percentage") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          axis.text.y = element_text(#colour = colours,
            size = "10",
            #angle = 75, 
            vjust = 1,
            hjust = 1),
          legend.text = element_text(colour = "black",size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          #strip.text.y = element_blank(),
          plot.caption = element_text(colour = "black",size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo3 type: "))
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project,"/figures/modulo_alltissues_exonintron_density")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  #######################
  ## BAR PLOT
  #######################
  
  df_modulo_tissues <- df_modulo_tissues %>%
    mutate(type = str_remove(type, pattern = "modulo")) %>%
    mutate(novel_t = ifelse(str_detect(novel_type, "donor"), "novel donor", "novel acceptor")) %>%
    mutate(locus_t = ifelse(str_detect(novel_type, "intron"), "intron", "exon"))
  
  #df_modulo_tissues$novel_type = factor(df_modulo_tissues$novel_type, 
  #                                      levels = c( "donor exon", "donor intron", "aceptor exon", "aceptor intron"))
  
  df_modulo_tissues <- df_modulo_tissues %>% 
    ungroup() %>%
    arrange(novel_t, locus_t, desc(modulo)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  # colours <- ifelse(str_detect(string = as.factor(df_modulo_tissues$tissue), 
  #                              pattern = "Brain"), "red", "black")
  # 
  
  # scales::show_col(scales::viridis_pal()(40))
  df_modulo_tissues %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)")
  
  plot_all <- ggplot(data = df_modulo_tissues) +
    geom_bar(mapping = aes(x = tissue, 
                           y = modulo, 
                           fill = factor(type)), 
             stat = "identity", 
             position = position_dodge()) +
    facet_grid(novel_t~locus_t) +
    coord_flip() +
    xlab("Tissue") +
    ylab("Percentage of novel junctions (%)") +
    theme_light() +
    #ggforce::facet_col(~novel_type) +
    scale_fill_viridis_d()  +
    #scale_fill_manual(values = c("#35B779FF","#440154FF","#31688EFF"),#"#FDE725FF", "#35B779FF", "#31688EFF","#440154FF"
    #                  breaks = c("modulo0","modulo1","modulo2"),
    #                  labels = c("modulo0","modulo1","modulo2")) +
    #scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          axis.text.y = element_blank(),#text(color = colours, vjust = 1, 
          #hjust = 1, size = "7"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo 3:", ncol = 3, nrow = 1))
  
  plot_all
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project,"/figures/modulo_alltissues_exonintron")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  #####################################
  ## STATS FOR THE PAPER
  #####################################
  
  df_modulo_tissues %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)") %>%
    dplyr::group_by(novel_type, type) %>%
    mutate(modulo_mean = modulo %>% mean) %>%
    distinct(modulo_mean) %>%
    arrange(novel_type, desc(modulo_mean))
    
  
}



get_junc_length <- function() {
  
  
  tables <- dbListTables(con)
  tables
  
  query <- paste0("SELECT * from 'intron'")
  db_introns <- dbGetQuery(con, query) %>% as_tibble()
  db_introns %>% distinct(ref_junID) %>% nrow()
  db_introns %>%
    dplyr::count(misspliced)
  
  
  query <- paste0("SELECT * from 'novel'")
  db_novel <- dbGetQuery(con, query) %>% as_tibble()
  db_novel %>% distinct(novel_junID) %>% nrow()
  db_novel %>% dplyr::count(novel_type)
  
  
  
  
  ## HISTOGRAM PLOT
  df_all_lengths_tidy <- rbind(db_introns %>%
                                 dplyr::select(length = ref_length) %>%
                                 mutate(type = "annotated intron"),
                               db_novel %>%
                                 dplyr::select(length = novel_length) %>%
                                 mutate(type = "novel junction"))
  
  
  get_mode( df_all_lengths_tidy %>%
              filter(type == "annotated intron") %>%
              pull(length) )
  df_all_lengths_tidy %>%
    filter(type == "annotated intron") %>%
    pull(length) %>% mean()
  
  
  df_all_lengths_tidy$type = factor( df_all_lengths_tidy$type, 
                                     levels = c( "novel junction", "annotated intron" ) )
  
  
  ggplot(data = df_all_lengths_tidy %>%
           drop_na() %>%
           filter(length <= 200)) + 
    geom_count(aes(x = length, 
                   colour = type), 
               #dplyr::mutate(df_all_lengths, z = FALSE), 
               alpha = 0.7,
               stat = "count",
               position = "identity") +
    #geom_histogram(aes(x = length, fill = type), binwidth = 5, 
    #               dplyr::mutate(df_all_lengths, z = TRUE), 
    #               alpha = 0.7, position = "identity") +
    # ggforce::facet_zoom(xlim = c(0, 200),
    #                    ylim = c(0, 710000),
    #                    zoom.data = z, horizontal = FALSE) +
    xlim(c(0,200))+
    ## ggtitle(paste0("Length of the annotated, partially unannotated and completely unannotated junctions.\nAll samples used - ", tissue)) +
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
                        labels = c("annotated intron", "novel junctions")) +
    theme(legend.position = "top",
          legend.text = element_text(size = 11)) %>%
    return()
  
  ## Save plot
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project,"/figures/junction_length")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
}

get_distances <- function() {
  
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  limit_bp <- 30
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  query <- paste0("SELECT tissue.novel_junID, tissue.ref_junID, novel.novel_type, novel.distance 
                  FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue
                  INNER JOIN 'novel' ON novel.novel_junID = tissue.novel_junID")
  db_misspliced_introns <- dbGetQuery(con, query) %>% as_tibble()
  # query <- paste0("SELECT * 
  #                 FROM 'novel' WHERE novel_junID IN (",
  #                 paste(introns$novel_junID, collapse = ","),")")
  # introns <- merge(x = introns,
  #                  y = dbGetQuery(con, query) %>% as_tibble(),
  #                  by = "novel_junID",
  #                  all.x = T) %>% 
  #   as_tibble() 
  
  
 
  #############################
  ## PLOT THE DISTANCES GRAPH 
  #############################
  
  
  df_novel_tidy <- db_misspliced_introns %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  
  #title <- paste0("Distances - ", cluster)
  
  
  #################################
  ## GENERATE PLOT ################
  #################################
  
  plot_all <- ggplot(data = df_novel_tidy) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggplot2::facet_grid(vars(novel_type)) +
    ggtitle(paste0("All biotypes\n")) +
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
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "10"),
          axis.text.x = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          strip.text = element_text(colour = "black", size = "10"), 
          legend.text = element_text(colour = "black", size = "10"),
          plot.caption = element_text(colour = "black", size = "10"),
          plot.title = element_text(colour = "black", size = "10"),
          legend.title = element_text(colour = "black", size = "10"),
          legend.position = "top") 
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 5, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size =5, label = "intron") +
    theme_void()


  plot_all <- plot_all / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  plot_all
  
  file_name <- paste0(getwd(), "/results/_paper/figures/distances_FCTX")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  

  ###############################
  ## PROTEIN CODING
  ###############################
 
  query <- paste0("SELECT ref_junID, protein_coding, lncRNA FROM 'intron'")
  df_master_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  df_novel_tidy <- db_misspliced_introns %>%
    inner_join(df_master_introns,
               by = "ref_junID") %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  
  df_novel_tidy$type_PC = factor(df_novel_tidy$type_PC, 
                                 levels = c("protein coding (PC)", "non PC"))
  
 
  
  #################################
  ## GENERATE PLOT ################
  #################################
  
  plot_PC <- ggplot(data = df_novel_tidy %>% 
                      filter(type_PC == "protein coding (PC)")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack") +
    ggplot2::facet_grid(vars(novel_type)) +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("novel donor", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "11"),
          axis.text.x = element_text(colour = "black", size = "11"),
          axis.title = element_text(colour = "black", size = "11"),
          strip.text.y = element_blank(), 
          legend.text = element_text(colour = "black", size = "11"),
          plot.caption = element_text(colour = "black", size = "11"),
          plot.title = element_text(colour = "black", size = "11"),
          legend.title = element_text(colour = "black", size = "11"),
          legend.position = "top") 
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 4, label = "exon") +
    geom_rect(aes(xmin = (limit_bp * -1), xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 80),  size = 4, label = "intron") +
    theme_void()
  
  
  plot_PC <- plot_PC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  plot_PC
  
  file_name <- paste0(getwd(), "/results/_paper/figures/FCTX_distancesPC")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 120, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 100, height = 100, units = "mm", dpi = 300)
  
  
  
  df_novel_tidy %>% filter(lncRNA == 100) %>% distinct(ref_junID) %>% nrow()
  
  plot_NPC <- ggplot(data = df_novel_tidy %>% 
                       filter(type_PC == "non PC")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggplot2::facet_grid(vars(novel_type)) +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("novel donor", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "11"),
          axis.text.x = element_text(colour = "black", size = "11"),
          axis.title = element_text(colour = "black", size = "11"),
          strip.text.y = element_blank(), 
          legend.text = element_text(colour = "black", size = "11"),
          plot.caption = element_text(colour = "black", size = "11"),
          plot.title = element_text(colour = "black", size = "11"),
          legend.title = element_text(colour = "black", size = "11"),
          legend.position = "top") 
  
  plot_NPC <- plot_NPC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  plot_NPC
  
  file_name <- paste0(getwd(), "/results/_paper/figures/FCTX_distancesNPC_lncRNA")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 100, height = 100, units = "mm", dpi = 300)
  
  
  #######################################
  ## COMBINE ALL 3 PLOTS
  #######################################
  
  plot <- ggpubr::ggarrange(plot_PC,
                            plot_NPC,
                            common.legend = T,
                            legend = "top",
                            labels = c("c", "d"),
                            align = "v",
                            ncol = 2, nrow = 1)


  
  
  plot

  
  file_name <- paste0(getwd(), "/results/_paper/figures/distances_FCTX_all")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 369, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 100, units = "mm", dpi = 300)

  
  #########################################
  ## STATS
  #########################################
  
  df_novel_tidy %>% 
    distinct(ref_junID, .keep_all = T) %>%
    dplyr::count(type_PC)

}




#################################
## MSR AND LINEAR MODELS
#################################

get_MSR_FCTX <- function()  {
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  print(paste0(Sys.time(), " - ", project_id, " - ", cluster_id))
  
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", 
                  cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", 
                  cluster_id, "_", project_id, "_misspliced'")
  introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())

  query <- paste0("SELECT DISTINCT ref_junID, protein_coding, lncRNA FROM 'intron' WHERE ref_junID IN (",
                  paste(introns$ref_junID, collapse = ","),")")
  introns <- introns %>%
    left_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = "ref_junID") %>% 
    as_tibble() 
  
  introns
  
  print(paste0(Sys.time(), " - ", introns %>% nrow(), " - introns!"))
  
  
  
  ###########################
  ## PLOT ALL INTRONS
  ###########################
  
  df_all_introns_tidy <- introns %>%
    dplyr::select(ref_junID, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -ref_junID)
  
  
  df_all_introns_tidy$MSR_type = factor(df_all_introns_tidy$MSR_type, 
                                        levels = c("MSR_D", "MSR_A"))
  
  df_all_introns_tidy$MSR_type %>% unique
  
  df_all_introns_tidy <- df_all_introns_tidy %>%
    mutate(percentile_group = case_when(MSR == 0 ~ "0",
                                        MSR > 0 & MSR <= 0.2 ~ "(0,.2]",
                                        MSR > 0.2 & MSR <= 0.4 ~ "(.2,.4]",
                                        MSR > 0.4 & MSR <= 0.6 ~ "(.4,.6]",
                                        MSR > 0.6 & MSR <= 0.8 ~ "(.6,.8]",
                                        MSR > 0.8 & MSR <= 1 ~ "(.8,1]"))
  
  df_all_introns_tidy$percentile_group = factor(df_all_introns_tidy$percentile_group, 
                                                levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  plotMSR <- ggplot(data = df_all_introns_tidy) + 
    geom_bar(aes(x = percentile_group,
                 fill = MSR_type),
             position = "dodge")+
    #ggtitle("All annotated introns") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    ggforce::facet_zoom(ylim=c(0,2500), split = TRUE) +
    #ylim(c(0,23000))+
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR Donor","MSR Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "8"),
          
          axis.title = element_text(colour = "black", size = "11"),
          strip.text =  element_text(colour = "black", size = "11"),
          plot.title = element_text(colour = "black", size = "11"),
          legend.text = element_text(size = "11"),
          legend.title = element_text(size = "11"),
          legend.position = "top",
          axis.text.x = element_text(hjust = 0.5, vjust = 1)) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  plotMSR
  
  file_name <- paste0(getwd(), "/results/_paper/figures/MSR_FCTX_all")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 160, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 80, units = "mm", dpi = 300)
  
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
    mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
    mutate(mean_coverage = mean_coverage %>% log10()) %>%
    ungroup()
  
  df_biotype_result_tidy %>%
    distinct(ref_junID, .keep_all = T) %>%
    dplyr::count(biotype)
  
  ## PLOT
  plot_BS <- ggplot(data = df_biotype_result_tidy %>% distinct(ref_junID,.keep_all = T)) +
    geom_density(mapping = aes(x = mean_coverage, fill = biotype), alpha = 0.8) +
    ggtitle("Before subsampling") +
    theme(legend.position = "top",
          legend.text = element_text(size = 12),
          text = element_text(size = 12)) +
    scale_fill_discrete(name = "Transcript biotype: ") +
    xlab("log10 mean read coverage")
  
  
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
  }
  
  ## Subsampling introns to control by similarity in mean read coverage
  m.out <- MatchIt::matchit(biotype ~ mean_coverage, 
                            data = data_combined, 
                            distance = data_combined$mean_coverage,
                            method = "nearest", 
                            caliper = c(mean_coverage = 0.005), 
                            std.caliper = FALSE)
  subsample <- MatchIt::match.data(m.out)
  
  subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)
  
  plot_AS <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = biotype), alpha = 0.8) +
    ggtitle("After subsampling") +
    #ggtitle("Mean read coverage per annotated intron across all samples\nfrom 54 GTEx v8 tissues - Subsampling performed.") +
    scale_fill_discrete(name = "Transcript biotype: ") +
    theme(legend.position = "top",
          legend.text = element_text(size = 12),
          text = element_text(size = 12)) +
    xlab("log10 mean read coverage")
  
  
  ggpubr::ggarrange(plot_BS,
                    plot_AS,
                    labels = c("a", "b"),
                    common.legend = T)
  
  file_name <- paste0(getwd(), "/results/_paper/figures/MSR_FCTX_subsampling")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
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
  
  df_introns_biotype <- subsample # %>%
    #filter(protein_coding %in% c(0,100)) %>%
    #mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC")) 
  
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
  
  scales::show_col(scales::viridis_pal(option = "B")(40))
  scales::show_col(scales::viridis_pal(option = "C")(40))
  scales::show_col(scales::viridis_pal(option = "D")(40))
  scales::show_col(scales::viridis_pal(option = "E")(40))
  scales::show_col(scales::viridis_pal(option = "F")(40))
  scales::show_col(scales::viridis_pal(option = "G")(40))
  scales::show_col(scales::viridis_pal(option = "H")(40))
  
  
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
    ylim(c(0,23000))+
    theme_light() +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          strip.text =  element_text(colour = "black", size = "10"),
          plot.title = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          legend.position = "top",
          axis.text.x = element_text(#angle = 20, 
                                     hjust = 0.5,
                                     vjust = 1)) +
    guides(fill = guide_legend(title = NULL,
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
    ylim(c(0,23000))+
    theme_light() +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          strip.text =  element_text(colour = "black", size = "10"),
          plot.title = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          legend.position = "top",
          axis.text.x = element_text(#angle = 20, 
                                     hjust = 0.5,
                                     vjust = 1)) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, nrow = 1))
  
  
  ggpubr::ggarrange(plotMSR_donor,
                    plotMSR_acceptor,
                    nrow = 1,
                    ncol = 2,
                    common.legend = T,
                    labels = c("b", "c"))
  
   
  
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0(getwd(), "/results/_paper/figures/MSR_FCTX_final")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 160, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 90, units = "mm", dpi = 300)
  

  
  
  
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
  
  
  ## 2. Furthermore, consistent with the overall higher detection of novel acceptor as compared to novel donor junctions,
  ## we observed a significant difference in the MSRD and MSRA density plots (paired Wilcoxon rank sum test with continuity 
  ## correction, pval = 4.76e-09). 
  
  wilcox.test(x = introns %>% pull(MSR_D),
              y = introns %>% pull(MSR_A),
              alternative = "less",
              paired = T,
              correct = T)
  
  wilcox.test(x = df_introns_biotype %>% filter(MSR_type == "MSR_D") %>% pull(MSR),
              y = df_introns_biotype %>% filter(MSR_type == "MSR_A") %>% pull(MSR),
              alternative = "less",
              paired = T,
              correct = T)
  
  
  ## 3. Given that NMD activity would be expected to reduce the detection of splicing errors amongst mRNA transcripts, 
  ## we separately assessed mis-splicing of annotated introns solely used in protein-coding transcripts (N=35,726 based 
  ## on Ensembl v105) and those that were exclusively used in non-coding transcripts (N=8,945 based on Ensembl v105).
  
  introns %>% filter(protein_coding == 100) %>% distinct(ref_junID) %>% nrow()
  introns %>% filter(protein_coding == 0) %>% distinct(ref_junID) %>% nrow()
  
  df_introns_biotype %>% filter(biotype == "PC") %>% distinct(ref_junID) %>% nrow()
  df_introns_biotype %>% filter(biotype == "non PC") %>% distinct(ref_junID) %>% nrow()
  
  
  ## 1. Splicing noise occurs more commonly at acceptor than at donor splice sites 
  ## regardless of whether the intron is exclusively used in protein-coding 
  ## or in non-protein-coding transcripts 

  
  ## Introns from PC transcripts are less miss-spliced at the donor than at the acceptor - protein coding
  wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_D") %>% pull(MSR),
              y = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_A") %>% pull(MSR),
              alternative = "less",
              paired = T,
              correct = T)
  rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                            filter(biotype == "PC") %>%
                            mutate(MSR_type = MSR_type %>% as.factor()),
                          formula = MSR ~ MSR_type,
                          paired = T)
  
  
  
  ## Introns from Non-PC transcripts are less miss-spliced at the donor than at the acceptor - non-protein-coding
  wilcox.test(x = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_D") %>% pull(MSR),
              y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_A") %>% pull(MSR),
              alternative = "less",
              paired = T,
              correct = T)
  rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                            filter(biotype == "non PC")%>%
                            mutate(MSR_type = MSR_type %>% as.factor()),
                          formula = MSR ~ MSR_type,
                          paired = T)
  
  
  
  ## 2. Splicing noise is more frequent in introns from non-protein coding transcripts 
  ## than in introns from protein-coding transcripts even after correcting by mean read coverage
  
  ## splicing noise at the donor is more frequent in introns from NPC vs PC 
  wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_D") %>% pull(MSR),
              y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_D") %>% pull(MSR),
              alternative = "less",
              paired = T,
              correct = T)
  rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                            filter(MSR_type == "MSR_D")%>%
                            mutate(biotype = biotype %>% as.factor()),
                          formula = MSR ~ biotype,
                          paired = T)
  
  
  ## splicing noise at the acceptor is more frequent in introns from NPC vs PC 
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
  ## CONNECT TO THE DATABASE
  ############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  ############################
  ## LOAD/GENERATE DATA
  ############################
  
  file_name <- paste0(paste0(getwd(), "/results/_paper/results/MSR_all_tissues.rds"))
  
  
  if ( !file.exists(file_name) ) {
    
    df_MSR_biotype <- map_df(all_projects, function(project_id) {
      
      #############################
      ## GET THE CLUSTERS
      #############################
      
      print(paste0(Sys.time(), " - ", project_id))
      
      # clusters_used <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
      #                                        project_id, "/raw_data/all_clusters_used.rds"))
      
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
        
        query <- paste0("SELECT DISTINCT ref_junID, protein_coding, lncRNA FROM 'intron' WHERE ref_junID IN (",
                        paste(introns$ref_junID, collapse = ","),")")
        introns <- introns %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
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
          mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
          mutate(mean_coverage = mean_coverage %>% log10()) %>%
          ungroup()
        
        df_biotype_result_tidy %>%
          distinct(ref_junID, .keep_all = T) %>%
          dplyr::count(biotype)
      
        
        #### START SUBSAMPLIING
        
        ## Subsampling introns to control by similarity in mean read coverage
        m.out <- MatchIt::matchit(biotype ~ mean_coverage, 
                                  data = df_biotype_result_tidy, 
                                  distance = df_biotype_result_tidy$mean_coverage,
                                  method = "nearest", 
                                  caliper = c(mean_coverage = 0.005), 
                                  std.caliper = FALSE)
        subsample <- MatchIt::match.data(m.out)
        
        subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(biotype)
        
        ##########################################
        ## DONOR IS LESS MISSPLICED THAN ACCEPTOR
        ##########################################
        
        df_introns_biotype <- subsample
        
        any(df_introns_biotype %>%
              pull(ref_junID) %>% duplicated())
        
        
        df_introns_biotype$biotype = factor(df_introns_biotype$biotype, 
                                            levels = c("PC", "non PC"))
        
        df_introns_biotype <- df_introns_biotype %>%
          dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
          gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)
        
        
        ## Introns from PC transcripts are less miss-spliced at the donor than at the acceptor
        wilcox_PC <- wilcox.test(x = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_D") %>% pull(MSR),
                                 y = df_introns_biotype %>% filter(biotype == "PC", MSR_type == "MSR_A") %>% pull(MSR),
                                 alternative = "less",
                                 paired = T,
                                 correct = T)
        effsize_PC <- rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                                                filter(biotype == "PC") %>%
                                                mutate(MSR_type = MSR_type %>% as.factor()),
                                              formula = MSR ~ MSR_type,
                                              paired = T)
        
        
        
        ## Introns from Non-PC transcripts are less miss-spliced at the donor than at the acceptor
        wilcox_NPC <- wilcox.test(x = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_D") %>% pull(MSR),
                                  y = df_introns_biotype %>% filter(biotype == "non PC", MSR_type == "MSR_A") %>% pull(MSR),
                                  alternative = "less",
                                  paired = T,
                                  correct = T)
        effsize_NPC <- rstatix::wilcox_effsize(data = df_introns_biotype %>% 
                                                 filter(biotype == "non PC")%>%
                                                 mutate(MSR_type = MSR_type %>% as.factor()),
                                               formula = MSR ~ MSR_type,
                                               paired = T)
        
        
        
        ##########################################
        ## PC IS LESS MISSPLICED THAN NON-PC
        ##########################################
        
        
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
        
        data.frame(test = c("PC_MSRD_vs_MSRA", 
                            "NPC_MSRD_vs_MSRA",
                            "MSRD_PC_vs_NPC",
                            "MSRA_PC_vs_NPC"),
                   hypothesis = c("Introns from PC transcripts are less miss-spliced at the donor than at the acceptor", 
                                  "Introns from Non-PC transcripts are less miss-spliced at the donor than at the acceptor", 
                                  "Introns from Non-PC transcripts are more mis-spliced at the donor than protein-coding introns",
                                  "Introns from Non-PC transcripts are more mis-spliced at the acceptor than protein-coding introns"),
                   pvalue = c(wilcox_PC$p.value,
                              wilcox_NPC$p.value,
                              wilcox_D$p.value,
                              wilcox_A$p.value),
                   effsize = c(effsize_PC$effsize %>% unname,
                               effsize_NPC$effsize %>% unname,
                               effsize_D$effsize %>% unname,
                               effsize_A$effsize %>% unname))  %>%
          mutate(tissue = cluster_id) %>%
          return()
        
      
        
        })
      
    })
    saveRDS(object = df_MSR_biotype, file = file_name)
    
  } else {
    df_MSR_biotype <- readRDS(file = file_name)
  }
  
  
  ###########################
  ## PLOT EFFECT SIZE 
  ## ACROSS TISSUES
  ###########################
  
  #scales::show_col(colours = viridis::viridis_pal(option = "A")(n = 20))
  
  # df_MSR_biotype <- df_MSR_biotype %>%
  #   arrange(tissue, effsize)
  # 
  # df_MSR_biotype <- df_MSR_biotype %>% 
  #   ungroup() %>%
  #   arrange(desc(type), effsize) %>%
  #   mutate(tissue = fct_inorder(tissue))
  # 
  # ggplot(data = df_MSR_biotype) +
  #   geom_bar(aes(x = tissue, 
  #                y = effsize, 
  #                fill = type),
  #            stat = "identity", 
  #            position = "dodge") +
  #   coord_flip()+
  #   xlab("") +
  #   ylab("Effect size") +
  #   theme_minimal() + 
  #   theme(axis.line = element_line(colour = "black"), 
  #         axis.text.x = element_text(colour = "black", size = "12"),
  #         axis.text.y = element_text(colour = "black", size = "12"),
  #         axis.title = element_text(colour = "black", size = "12"),
  #         legend.text = element_text(colour = "black", size = "12"),
  #         legend.title = element_text(colour = "black", size = "12"),
  #         legend.position = "top",
  #         panel.grid.major.x = element_blank(),
  #         panel.grid.major.y = element_blank(),
  #         axis.ticks = element_line(colour = "black", size = 2)) +  
  #   guides(fill = guide_legend(title = "Transcript Biotype:",
  #                                ncol = 2, 
  #                                nrow = 1))    + 
  #   scale_fill_manual(values = c("#56147DFF","#FD9A6AFF"),
  #                     breaks = c("PC", "NPC"),
  #                     labels = c("PC", "NPC")) 
  # 
  # 
  # file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
  #                     main_project, "/figures/MSR_effsize_alltissues")
  # ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  # ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
}

get_lm_single_tissue <- function() {
  
  ###############################
  ## QUERY THE DATABASE
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  
  ###################################
  ## GET TPM
  
  base_folder <- paste0(getwd(),"/results/", project_id, "/v", gtf_version, "/", main_project, "/")
  samples_used <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version, 
                                        "/", main_project, "/base_data/", project_id, "_", cluster_id, "_samples_used.rds"))
  
  
  tpm <- readRDS(file = paste0(base_folder, "results/tpm/",
                               "/", project_id, "_", cluster_id, "_tpm.rds")) %>% 
    dplyr::select(gene_id = gene, all_of(samples_used))
  
  
  tpm <- tpm  %>%
    dplyr::mutate(tpm_median = matrixStats::rowMedians(x = as.matrix(.[2:(ncol(tpm))]))) %>%
    dplyr::select(gene_id, tpm_median) 
  
  #######################################
  
  
  query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
  
  query <- paste0("SELECT ref_junID, ref_length, ref_ss5score, ref_ss3score, 
                  ref_cons5score, ref_cons3score, ref_CDTS5score, ref_CDTS3score, 
                  protein_coding
                  FROM 'intron' 
                  WHERE ref_junID IN (", paste(introns$ref_junID, collapse = ","),")")
  introns <- introns %>%
    left_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = "ref_junID") %>% 
    as_tibble() 
  
  ## JOIN WITH TRANSCRIPT DATA
  query <- paste0("SELECT * FROM 'transcript' WHERE id IN (", paste(introns$transcript_id, collapse = ","),")")
  introns <- introns %>%
    left_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = c("transcript_id" = "id")) %>% 
    as_tibble() 
  
  
  ## JOIN WITH GENE DATA
  query <- paste0("SELECT *
                  FROM 'gene' WHERE id IN (",
                  base::paste(introns$gene_id, collapse = ","),")")
  introns <- introns %>%
    left_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = c("gene_id" = "id")) %>% 
    as_tibble() 
  
  
  #########################
  ## TIDY DATA
  #########################
  
  idb <- introns %>%
    left_join(y = tpm,
              by = c("gene_id.y"="gene_id")) %>% 
    dplyr::select(-gene_id) %>%
    dplyr::rename(gene_tpm = tpm_median) %>%
    dplyr::select(ref_junID,
                  intron_length = ref_length,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score,
                  gene_length = gene_width,
                  gene_tpm = gene_tpm,
                  gene_num_transcripts = n_transcripts,
                  CDTS_5ss = ref_CDTS5score,
                  CDTS_3ss = ref_CDTS3score,
                  protein_coding = protein_coding,
                  mean_phastCons20way_5ss = ref_cons5score,
                  mean_phastCons20way_3ss = ref_cons3score,
                  MSR_D,
                  MSR_A)
  
  # idb <- introns %>%
  #   as.data.frame() %>%
  #   dplyr::distinct(ref_junID, .keep_all = T) %>%
  #   dplyr::rename(intron_length = ref_length,
  #                 intron_5ss_score = ref_ss5score,
  #                 intron_3ss_score = ref_ss3score,
  #                 gene_length = gene_width,
  #                 gene_tpm = gene_tpm,
  #                 gene_num_transcripts = n_transcripts,
  #                 CDTS_5ss = ref_CDTS5score,
  #                 CDTS_3ss = ref_CDTS3score,
  #                 protein_coding = protein_coding,
  #                 mean_phastCons20way_5ss = ref_cons5score,
  #                 mean_phastCons20way_3ss = ref_cons3score)
  
  
  #########################
  ## LINEAR MODELS
  #########################
  
  
  fit_donor <- lm(MSR_D ~ 
                    gene_length +
                    gene_tpm +
                    gene_num_transcripts +
                    protein_coding +
                    intron_length +
                    intron_5ss_score + 
                    intron_3ss_score +
                    CDTS_5ss +
                    CDTS_3ss +
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss, 
                  data = idb)
  fit_donor %>% summary()
  
  # fit_donor$coefficients[names(fit_donor$coefficients)[summary(fit_donor)$coefficients[,4] < 0.05]]  
    
  fit_acceptor <- lm(MSR_A ~ 
                       gene_length +
                       gene_tpm +
                       gene_num_transcripts +
                       protein_coding +
                       intron_length + 
                       intron_5ss_score + 
                       intron_3ss_score +
                       CDTS_5ss +
                       CDTS_3ss +
                       mean_phastCons20way_5ss +
                       mean_phastCons20way_3ss, 
                     data = idb)
  fit_acceptor %>% summary()
  
  model_names <- c("MSR Donor", "MSR Acceptor")
  
  coef_names <- c("Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  "Protein coding" = "protein_coding",
                  "Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score",
                  "CDTS 5'ss" = "CDTS_5ss",
                  "CDTS 3'ss" = "CDTS_3ss",
                  "PhastCons20 5'ss" = "mean_phastCons20way_5ss",
                  "PhastCons20 3'ss" = "mean_phastCons20way_3ss")
  
  
 
  
  ##############################################################################
  
  plotLM <- jtools::plot_summs(fit_donor, 
                               fit_acceptor,
                               #scale = TRUE, 
                               #robust = T,
                               #n.sd = 2,
                               #pvals = TRUE,
                               legend.title = "",
                               #plot.distributions = TRUE,
                               ci_level = 0.95,
                               coefs = coef_names,
                               colors = c("#35B779FF","#64037d"),
                               model.names = model_names,
                               plot.distributions = F,
                               rescale.distributions = T,
                               point.shape = T) +
    guides(colour = guide_legend(ncol = 2, nrow = 1)) +
    theme_minimal() + 
    theme(axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = "black", size = "10"),
          axis.text.y = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top",
          #panel.grid.major.x = element_blank(),
          #panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 2)) +  
    ylab("Covariates") +
    geom_hline(yintercept = seq(from = 0,
                                to = length((fit_donor$coefficients %>% names)[-1]) + .5,
                                by = 1),
               colour = "#999999") +
    guides(colour = guide_legend(ncol = 2, nrow = 1))     +
    scale_x_continuous(breaks = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025),
                       labels = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025))
  
  
  
  tiles_data <- rbind(summary(fit_donor)$coefficients[,4] %>% 
                        as_tibble() %>%
                        mutate(value = ifelse(value == 0, max(value), value)) %>%
                        mutate(q = p.adjust(value,method = "bonferroni")) %>%
                        mutate(q = value %>% log10())%>%
                        mutate(names = names(summary(fit_donor)$coefficients[,4] ),
                               type = "MSR Donor") %>%
                        dplyr::rename(pvalue = value) %>%
                        filter(names != "(Intercept)") %>%
                        rev(),
                      summary(fit_acceptor)$coefficients[,4] %>% 
                        as_tibble() %>%
                        mutate(value = ifelse(value == 0, max(value), value)) %>%
                        mutate(q = p.adjust(value,method = "bonferroni")) %>%
                        mutate(q = q %>% log10())%>%
                        mutate(names = names(summary(fit_acceptor)$coefficients[,4] ),
                               type = "MSR Acceptor") %>%
                        dplyr::rename(pvalue = value) %>%
                        filter(names != "(Intercept)") %>%
                        rev())
  
  tiles_data <- tiles_data %>%
    inner_join(y = data.frame(col_name = c("gene_length","gene_tpm","gene_num_transcripts","protein_coding","intron_length","intron_5ss_score","intron_3ss_score",
                                              "CDTS_5ss","CDTS_3ss","mean_phastCons20way_5ss","mean_phastCons20way_3ss"),
                                 col_label = c("Gene Length","Gene TPM","Gene num. transcripts","Protein coding",
                                               "Intron Length","Intron 5'ss MES score","Intron 3'ss MES score","CDTS 5'ss","CDTS 3'ss",
                                               "PhastCons20 5'ss","PhastCons20 3'ss")) %>% as_tibble(),
               by = c("names" = "col_name"))
  
  tiles_data$col_label <- factor(tiles_data$col_label, levels= c("Gene Length","Gene TPM","Gene num. transcripts","Protein coding",
                                                         "Intron Length","Intron 5'ss MES score","Intron 3'ss MES score","CDTS 5'ss","CDTS 3'ss",
                                                         "PhastCons20 5'ss","PhastCons20 3'ss") %>% rev())
  tiles_data$type <- factor(tiles_data$type, levels= c("MSR Donor",
                                                         "MSR Acceptor"))
  
  ggpubr::ggarrange(plotLM +
                      theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 0.9)),
                    ggplot() + 
                      geom_tile(data = tiles_data, mapping = aes(x = type, y = col_label, fill = q))+
                      scale_fill_gradient(low = "red", high = "white") +
                      theme(legend.position = "top") +
                      xlab(" ") + 
                      ylab("") + 
                      theme_minimal() + 
                      theme(axis.line = element_line(colour = "black"), 
                            axis.text.x = element_text(colour = "black", size = "10"),
                            axis.text.y = element_text(colour = "black", size = "12"),
                            axis.title = element_text(colour = "black", size = "12"),
                            legend.text = element_text(colour = "black", size = "12"),
                            legend.title = element_text(colour = "black", size = "12"),
                            legend.position = "top",
                            #panel.grid.major.x = element_blank(),
                            #panel.grid.major.y = element_blank(),
                            axis.ticks = element_line(colour = "black", linewidth = 2))+
                      theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust = 0.9)),
                    ncol = 2,
                    nrow = 1)
    
    # geom_point(data = summary(fit_acceptor)$coefficients[,4] %>% 
    #              as_tibble() %>%
    #              mutate(names = names(summary(fit_acceptor)$coefficients[,4] ))  %>%
    #              dplyr::rename(pvalue = value) %>%
    #              filter(names != "(Intercept)") %>%
    #              rev(),
    #            aes(x = 0.003,
    #                y = c(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1), 
    #                size = pvalue),
    #            colour = "#64037d") +
    # scale_size(range=c(8, 1), breaks=c(0, 1e-174, 1e-31, 1e-12, 1e-1))  +
    # theme(legend.position="top", legend.box="vertical", legend.margin=margin())
  
  
  
  file_name <- paste0(getwd(), "/results/_paper/figures/lm_FTCX")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 200, height = 90, units = "mm", dpi = 300)
  
  
  ############################################
  ## STATS
  ############################################
  
  fit_donor %>% summary()
  fit_acceptor %>% summary()
  
  
  
}

get_common_introns_across_tissues <- function () {
  
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
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
      
      query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      
      all_introns[[cluster_id]] <- introns$ref_junID
      
      print(paste0(Sys.time(), " - intron IDs collected from '", cluster_id, "'"))
    }
    
  }
  
  common_introns <- data.frame(ref_junID = Reduce(intersect,  all_introns))
  common_introns %>% head()
  common_introns %>% nrow()
  
  query <- paste0("SELECT * FROM 'intron' WHERE ref_junID IN (", paste(common_introns$ref_junID, collapse = ","), ")")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  
  file_name <- paste0(getwd(), "/results/_paper/results/common_introns.rds")
  saveRDS(object = introns, file = file_name)
  
}

get_estimate_variance_across_tissues <- function() {
  
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  ## LOAD COMMON INTRONS
  common_introns <- readRDS(file = paste0(getwd(), "/results/_paper/results/common_introns.rds"))

  ################################
  ## QUERY THE DATABASE PER TISSUE
  ################################

  df_estimates <- map_df(all_projects, function(project_id) {
    
    
    # project_id <- all_projects[1]
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
      print(cluster_id)
      
      
      ###################################
      ## GET TPM
      
      base_folder <- paste0(getwd(),"/results/", project_id, "/v", gtf_version, "/", main_project, "/")
      samples_used <- readRDS(file = paste0(getwd(), "/results/", project_id, "/v", gtf_version, 
                                            "/", main_project, "/base_data/", project_id, "_", cluster_id, "_samples_used.rds"))
    
      tpm <- readRDS(file = paste0(base_folder, "results/tpm/",
                                   "/", project_id, "_", cluster_id, "_tpm.rds")) %>% 
        dplyr::select(gene_id = gene, all_of(samples_used))
      
      tpm <- tpm  %>%
        dplyr::mutate(tpm_median = matrixStats::rowMedians(x = as.matrix(.[2:(ncol(tpm))]))) %>%
        dplyr::select(gene_id, gene_tpm = tpm_median) 
      
      ##################################
      
      query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      
      introns <- introns %>%
        filter(ref_junID %in% common_introns$ref_junID)
      
      ## Add MASTER INTRON INFO
      query <- paste0("SELECT ref_junID, ref_length AS intron_length, ref_ss5score AS intron_5ss_score, 
                      ref_ss3score AS intron_3ss_score, 
                      ref_cons5score AS mean_phastCons20way_5ss, 
                      ref_cons3score AS mean_phastCons20way_3ss, 
                      ref_CDTS5score AS CDTS_5ss, 
                      ref_CDTS3score AS CDTS_3ss, 
                      protein_coding, transcript_id 
                      FROM 'intron' WHERE ref_junID IN (",
                      paste(introns$ref_junID, collapse = ","),")")
      introns <- merge(x = introns,
                       y = dbGetQuery(con, query) %>% as_tibble(),
                       by = "ref_junID",
                       all.x = T) %>% 
        as_tibble() 
      
      
      ## JOIN WITH TRANSCRIPT DATA
      query <- paste0("SELECT * FROM 'transcript' WHERE id IN (", paste(introns$transcript_id, collapse = ","),")")
      introns <- introns %>%
        left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = c("transcript_id" = "id")) %>% 
        as_tibble() 
      
      
      ## JOIN WITH GENE DATA
      query <- paste0("SELECT *
                  FROM 'gene' WHERE id IN (",
                      base::paste(introns$gene_id, collapse = ","),")")
      introns <- introns %>%
        left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = c("gene_id" = "id")) %>% 
        as_tibble()  %>%
        dplyr::rename(gene_length = gene_width)
      
      # query <- paste0("SELECT id, n_transcripts AS gene_num_transcripts, 
      #                 gene_width AS gene_length, gene_name FROM 'gene' WHERE id IN (",
      #                 paste(introns$gene_id, collapse = ","),")")
      # introns <- merge(x = introns,
      #                  y = dbGetQuery(con, query) %>% as_tibble(),
      #                  by.x = "gene_id",
      #                  by.y = "id",
      #                  all.x = T) %>% 
      #   as_tibble() 
      
      
      ## JOIN WITH TPM DATA
      introns <- introns %>%
        left_join(y = tpm,
                  by = c("gene_id.y" = "gene_id"))
      
      introns <- introns %>%
        dplyr::rename(gene_num_transcripts = n_transcripts )
      
      
      introns %>% distinct(ref_junID) %>% nrow() %>% print()
      
      introns
      
      model <- lm(MSR_D ~
                    intron_length +
                    intron_5ss_score +
                    intron_3ss_score +
                    gene_tpm +
                    gene_length +
                    gene_num_transcripts +
                    protein_coding +
                    CDTS_5ss +
                    CDTS_3ss +
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss,
                  data = introns)

      #MSR_Donor_list[[tissue]] <- model

      MSR_Donor <- data.frame(feature = cluster_id,
                              intron_length = model$coefficients["intron_length"] %>% unname(),
                              intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                              intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                              gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                              gene_length = model$coefficients["gene_length"] %>% unname(),
                              #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),

                              gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                              protein_coding = model$coefficients["protein_coding"] %>% unname(),
                              CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                              CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                              mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                              mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname())



      #model %>% summary() %>% print()
      MSR_Donor %>% mutate(type = "MSR_Donor")

      MSR_Donor <- cbind(MSR_Donor %>%
              gather(key = feature, value = estimate),
            pval = ((model %>% summary())$coefficients  %>% as.data.frame())[-1,4]) %>%
        mutate(tissue = cluster_id)

      ## Acceptor
      model <- lm(MSR_A ~
                    intron_length +
                    intron_5ss_score +
                    intron_3ss_score +
                    gene_tpm +
                    gene_length +
                    gene_num_transcripts +
                    # u2_intron +
                    # clinvar +
                    protein_coding +
                    CDTS_5ss +
                    CDTS_3ss +
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss,
                  data = introns)
      #
      #MSR_Acceptor_list[[tissue]] <- model

      #ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
      #model$coefficients[-ind_sign] <- 0

      MSR_Acceptor <- data.frame(feature = cluster_id,
                                 intron_length = model$coefficients["intron_length"] %>% unname(),
                                 intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                                 intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                                 gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                                 gene_length = model$coefficients["gene_length"] %>% unname(),
                                 #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),

                                 gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                                 protein_coding = model$coefficients["protein_coding"] %>% unname(),
                                 CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                                 CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                                 mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                                 mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname())

      #model %>% summary() %>% print()
      model$coefficients

      MSR_Acceptor <- cbind(MSR_Acceptor %>%
                           gather(key = feature, value = estimate),
                         pval = ((model %>% summary())$coefficients  %>% as.data.frame())[-1,4]) %>%
        mutate(tissue = cluster_id)

      return( rbind(MSR_Donor %>% mutate(type = "MSR_Donor"),
                    MSR_Acceptor %>% mutate(type = "MSR_Acceptor")) )
      
    })
  })
  
  
  #############################
  ## SAVE RESULTS
  #############################
  
  file_name <- paste0(getwd(), "/results/_paper/results/variance_estimate_tissues.rds")
  saveRDS(object = df_estimates, file = file_name)
  
}

plot_estimate_variance_across_tissues <- function() {
  

  df_estimate <- readRDS(file = paste0(getwd(), "/results/_paper/results/variance_estimate_tissues.rds")) %>%
    group_by(type) %>%
    mutate(q = p.adjust(pval, method = "fdr")) %>%
    ungroup() %>%
    drop_na()
  
  df_estimate %>%
    # group_by(tissue) %>%
    dplyr::count(feature)
  
  ## Interestingly, this analysis also demonstrated that while the effect of sequence conservation on splicing 
  ## noise was always significant, its effect size was variable (XXXXXrange of betas, Figure). 
  
  ## Donor - pval
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_5ss",
           type == "MSR_Donor") %>%
    ungroup() %>%
    pull(pval_corrected) %>% summary()
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_3ss",
           type == "MSR_Donor") %>%
    ungroup() %>%
    pull(pval_corrected) %>% summary()
  
  
  ## Acceptor - pval
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_5ss",
           type == "MSR_Acceptor") %>%
    ungroup() %>%
    pull(pval_corrected) %>% summary()

  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_3ss",
           type == "MSR_Acceptor") %>%
    ungroup() %>%
    pull(pval_corrected) %>% summary()
  
  
  
  
  
  ## Donor - effect size
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_5ss",
           type == "MSR_Donor") %>%
    ungroup() %>%
    pull(estimate) %>% summary()
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_3ss",
           type == "MSR_Donor") %>%
    ungroup()  %>%
    pull(estimate) %>% summary()
  
  
  ## Acceptor - effect size
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_5ss",
           type == "MSR_Acceptor") %>%
    ungroup()%>%
    pull(estimate) %>% summary()
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_3ss",
           type == "MSR_Acceptor") %>%
    ungroup() %>%
    pull(estimate) %>% summary()
 

  
  
  
  
  
  df_estimate <- df_estimate %>%
    #group_by(type) %>%
    #filter(pval_corrected <= 0.05) %>%
    #ungroup() %>%
    dplyr::select(-c(pval, pval_corrected)) %>%
    drop_na()
  
  MSR_Donor <- df_estimate %>%
    filter(type == "MSR_Donor") %>%
    mutate(type = "MSR Donor") %>%
    spread(key=feature, value= estimate)
  MSR_Acceptor <- df_estimate %>%
    filter(type == "MSR_Acceptor") %>%
    mutate(type = "MSR Acceptor") %>%
    spread(key=feature, value= estimate)
  
  #############################
  ## PLOT - ESTIMATE VARIANCE
  #############################
  
  ## DONOR
  
  MSR_Donor_tidy <- MSR_Donor %>% 
    dplyr::rename("Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  "Protein coding" = "protein_coding",
                  "Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  "CDTS 5'ss" = "CDTS_5ss",
                  "CDTS 3'ss" = "CDTS_3ss",
                  "PhastCons20 5'ss" = "mean_phastCons20way_5ss",
                  "PhastCons20 3'ss" = "mean_phastCons20way_3ss") %>%
    gather(tissue, feature, -type, -tissue) %>% as_tibble() %>%
    mutate(type = "MSR Donor")
  

  
  MSR_Acceptor_tidy <- MSR_Acceptor %>%
    dplyr::rename("Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  "Protein coding" = "protein_coding",
                  "Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  "CDTS 5'ss" = "CDTS_5ss",
                  "CDTS 3'ss" = "CDTS_3ss",
                  "PhastCons20 5'ss" = "mean_phastCons20way_5ss",
                  "PhastCons20 3'ss" = "mean_phastCons20way_3ss") %>%
    gather(tissue, feature, -type, -tissue) %>%
    mutate(type = "MSR Acceptor")
  
  
  label_faceting <- data.frame(label = c( "Gene Length",
                                           "Gene TPM",
                                           "Gene num. transcripts",
                                           "Protein coding",
                                           "Intron Length",
                                           "Intron 5'ss MES score",
                                           "Intron 3'ss MES score", 
                                           "CDTS 5'ss",
                                           "CDTS 3'ss",
                                           "PhastCons20 5'ss",
                                           "PhastCons20 3'ss" ),
                               subgroup = c( "gene level",
                                             "gene level",
                                             "gene level",
                                             "gene level",
                                             "intron level",
                                             "intron level",
                                             "intron level", 
                                             "intron level",
                                             "intron level",
                                             "intron level",
                                             "intron level" ))
  
  MSR_Donor_tidy <- MSR_Donor_tidy %>%
    inner_join(y = label_faceting,
               by = c("tissue" = "label"))
  MSR_Acceptor_tidy <- MSR_Acceptor_tidy %>%
    inner_join(y = label_faceting,
               by = c("tissue" = "label"))
  
    ## TIDY VALUES
  MSR_Donor_tidy$tissue <- factor( MSR_Donor_tidy$tissue, 
                                   levels = c( "Gene Length",
                                               "Gene TPM",
                                               "Gene num. transcripts",
                                               "Protein coding",
                                               "Intron Length",
                                               "Intron 5'ss MES score",
                                               "Intron 3'ss MES score", 
                                               "CDTS 5'ss",
                                               "CDTS 3'ss",
                                               "PhastCons20 5'ss",
                                               "PhastCons20 3'ss") %>% rev() )
  MSR_Acceptor_tidy$tissue <- factor( MSR_Acceptor_tidy$tissue, 
                                      levels = c( "Gene Length",
                                                  "Gene TPM",
                                                  "Gene num. transcripts",
                                                  "Protein coding",
                                                  "Intron Length",
                                                  "Intron 5'ss MES score",
                                                  "Intron 3'ss MES score", 
                                                  "CDTS 5'ss",
                                                  "CDTS 3'ss",
                                                  "PhastCons20 5'ss",
                                                  "PhastCons20 3'ss" ) %>% rev() )
  ## PLOT
  plotTissues5ssLM <- ggplot(data = MSR_Donor_tidy, aes(tissue, feature, fill=feature)) + 
    geom_boxplot(fill = "#35B779FF") +
    coord_flip() +
    #ggforce::facet_col(vars(type)) + 
    facet_grid(vars(subgroup), scales = "free", switch = "y", space = "free_y")  +
    #ggtitle(graph_title) +
    ylab("Distribution of the significant beta values (q<0.05)") +
    xlab(" ") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    scale_fill_manual(breaks = c("MSR_Donor","MSR_Acceptor"),
                       labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.text.y = element_text(#angle = 70, 
                                     vjust = 0.5,
                                     hjust = 1)) +
    geom_hline(yintercept = 0,linetype='dotted')
  
  
  plotTissues5ssLM
  file_name <- paste0(getwd(), "/results/_paper/figures/lm_donor_alltissues")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 120, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 80, units = "mm", dpi = 300)
  
  
  
  plotTissues3ssLM <- ggplot(data = MSR_Donor_tidy, aes(tissue, feature, fill=feature)) + 
    geom_boxplot(fill = "#8d03b0") +
    coord_flip() +
    #ggforce::facet_col(vars(type)) + 
    facet_grid(vars(subgroup), scales = "free", switch = "y", space = "free_y")  +
    #ggtitle(graph_title) +
    ylab("Distribution of the significant beta values (q<0.05)") +
    xlab(" ") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    scale_fill_manual(breaks = c("MSR_Donor","MSR_Acceptor"),
                      labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.text.y = element_text(#angle = 70, 
      vjust = 0.5,
      hjust = 1)) +
    geom_hline(yintercept = 0,linetype='dotted')
  
  
  plotTissues3ssLM
  file_name <- paste0(getwd(), "/results/_paper/figures/lm_accepor_alltissues")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 120, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 80, units = "mm", dpi = 300)
  
  
  #######################################
  ## TISSUES WITH THE MOST EXTREME CONSERVATION
  #######################################
  
  ## Donor - cons 5'ss (desc) and (asc)
  df_estimate %>%
    filter(type == "MSR_Donor",
           feature == "mean_phastCons20way_5ss") %>%
    arrange(desc(estimate)) %>%
    distinct(tissue, .keep_all = T)
  df_estimate %>%
    filter(type == "MSR_Donor",
           feature == "mean_phastCons20way_5ss") %>%
    arrange(estimate) %>%
    distinct(tissue, .keep_all = T)
  
  
  ## Donor - cons 3'ss (desc) and (asc)
  df_estimate %>%
    filter(type == "MSR_Donor",
           feature == "mean_phastCons20way_3ss") %>%
    arrange(desc(estimate)) %>%
    distinct(tissue, .keep_all = T)
  df_estimate %>%
    filter(type == "MSR_Donor",
           feature == "mean_phastCons20way_3ss") %>%
    arrange(estimate) %>%
    distinct(tissue, .keep_all = T)
  
  
  ## Acceptor - cons 5'ss (desc) and (asc)
  df_estimate %>%
    filter(type == "MSR_Acceptor",
           feature == "mean_phastCons20way_5ss") %>%
    arrange(desc(estimate)) %>%
    distinct(tissue, .keep_all = T)
  df_estimate %>%
    filter(type == "MSR_Acceptor",
           feature == "mean_phastCons20way_5ss") %>%
    arrange(estimate) %>%
    distinct(tissue, .keep_all = T)
  
  
  ## Donor - cons 3'ss (desc) and (asc)
  df_estimate %>%
    filter(type == "MSR_Acceptor",
           feature == "mean_phastCons20way_3ss") %>%
    arrange(desc(estimate)) %>%
    distinct(tissue, .keep_all = T)
  df_estimate %>%
    filter(type == "MSR_Acceptor",
           feature == "mean_phastCons20way_3ss") %>%
    arrange(estimate) %>%
    distinct(tissue, .keep_all = T)
  
}


#################################
## SOMATIC MUTATIONS - SKIN
#################################


# project_id1 = "SKIN"
# cluster_id1 = "Skin - Sun Exposed (Lower leg)"
# project_id2 = "SKIN"
# cluster_id2 = "Skin - Not Sun Exposed (Suprapubic)"
# stats = F
# 
# 
# project_id1 = "BRAIN"
# cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)"
# project_id2 = "PANCREAS"
# cluster_id2 = "Pancreas"

project_id1 = "BRAIN"
cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)"
project_id2 = "BLOOD"
cluster_id2 = "Whole Blood"

compare_tissues_somatic_mutations <- function(project_id1 = "SKIN",
                                              cluster_id1 = "Skin - Sun Exposed (Lower leg)",
                                              project_id2 = "SKIN",
                                              cluster_id2 = "Skin - Not Sun Exposed (Suprapubic)",
                                              stats = F) {
  
  tables <- c(paste0(cluster_id1, "_", project_id1), 
              paste0(cluster_id2, "_", project_id2)) #c("Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)")
  
  
  if ( !file.exists(paste0(getwd(), "/results/_paper/results/somatic_mutations_subsampled_", 
                           project_id1, "_", project_id2, ".rds")) ) {
    
    database_introns <- map_df(tables, function(table) {
      
      # table <- tables[1]
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals, transcript_id FROM '", 
                      table, "_nevermisspliced'")
      df_introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals, transcript_id FROM '", 
                      table, "_misspliced'")
      df_introns <- rbind(df_introns, dbGetQuery(con, query) %>% as_tibble())
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_length, ref_ss5score, ref_ss3score, 
                  ref_cons5score, ref_cons3score, ref_CDTS5score, ref_CDTS3score, 
                  protein_coding, lncRNA FROM 'intron' WHERE ref_junID IN (",
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
            file = paste0(getwd(), "/results/_paper/results/introns_", project_id1, "_", project_id2, ".rds"))
    
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
            file = paste0(getwd(), "/results/_paper/results/somatic_mutations_subsampled_", project_id1, "_", project_id2, ".rds"))
    
    print(paste0(Sys.time(), " - subsampling finished!"))
  } else {
    
    df_database_introns_tidy <- readRDS(file = paste0(getwd(), "/results/_paper/results/introns_", project_id1, "_", project_id2, ".rds"))
    subsample <- readRDS(paste0(getwd(), "/results/_paper/results/somatic_mutations_subsampled_", project_id1, "_", project_id2, ".rds"))
    
  }
  
  ##########################
  ## PLOT - COVERAGE
  ##########################
  
  ## VISUALISE MEAN READ COVERAGE BEFORE SUBSAMPLING
  plot_coverage_bf <- ggplot(data = df_database_introns_tidy %>% distinct(ref_junID,.keep_all = T)) +
    geom_density(mapping = aes(x = mean_coverage, fill = tissue), alpha = 0.8) +
    ggtitle("Before subsampling") +
    theme(legend.position = "top",
          legend.text = element_text(size = 12),
          text = element_text(size = 12)) +
    scale_fill_discrete(name = "Tissue: ") +
    xlab("log10 mean read coverage")
  plot_coverage_bf
  
  ## VISUALISE MEAN READ COVERAGE AFTER SUBSAMPLING
  plot_coverage_af <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = tissue), alpha = 0.8) +
    ggtitle("After subsampling") +
    #ggtitle("Mean read coverage per annotated intron across all samples\nfrom 54 GTEx v8 tissues - Subsampling performed.") +
    scale_fill_discrete(name = "Tissue: ") +
    theme(legend.position = "top",
          legend.text = element_text(size = 12),
          text = element_text(size = 12)) +
    xlab("log10 mean read coverage")
  plot_coverage_af
  
  ggpubr::ggarrange(plot_coverage_bf,
                    plot_coverage_af,
                    labels = c("a", "b"),
                    common.legend = T)
  file_name <- paste0(getwd(), "/results/_paper/figures/", project_id1, "_", project_id2, "_somatic_mutations_coverage")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  ################################
  ## ONLY EVALUATE COMMON INTRONS
  ################################
  
  ## TIDY THE DATA
  ## ONLY CHECK COMMON INTRONS
  
  common_introns <- subsample %>%
    group_by(tissue) %>%
    distinct(ref_junID) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == 2) 
  
  df_coverage <- subsample %>%
    filter(ref_junID %in% common_introns$ref_junID)
  
  df_coverage$tissue = factor(df_coverage$tissue, 
                              levels = c(paste0(cluster_id1, "_", project_id1),
                                         paste0(cluster_id2, "_", project_id2)))
  
  
  df_coverage %>% dplyr::count(tissue)
  
  
  df_coverage <- df_coverage %>%
    dplyr::select(ref_junID, tissue, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -tissue, -ref_junID)
  
  df_coverage$MSR_type = factor(df_coverage$MSR_type, 
                                       levels = c("MSR_A", "MSR_D"))
  print(paste0(Sys.time(), " - ", df_coverage %>% nrow(), " - introns after tidying!"))
  

  
  
  ################################
  ## PLOT - MSR
  ################################
  
  plotSKIN_sunexposed <- ggplot(data = df_coverage %>% 
                         filter(tissue == paste0(cluster_id1, "_", project_id1))) + 
    geom_bar(aes(x = MSR,  fill = MSR_type), 
             stat = "bin", 
             bins = 30,
             position = "dodge")  +
    ggtitle(paste0(cluster_id1)) +
    xlab("Mis-splicing ratio") +
    ylab("Intron count") +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text =  element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  plotSKIN_sunexposed
  
  
  
  ## SKIN Non-Sunexposed
  plotSKIN_nonexposed <- ggplot(data = df_coverage %>% 
                          filter(tissue == paste0(cluster_id2, "_", project_id2))) + 
    #geom_density(aes(x = MSR, fill = MSR_type), alpha = 0.8) +
    geom_bar(aes(x = MSR,  fill = MSR_type), 
             stat = "bin", 
             bins = 30,
             position = "dodge") +
    ggtitle(cluster_id2) +
    xlab("Mis-splicing ratio") +
    ylab("Intron count") +
    ylim(c(0,25000)) +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text =  element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  plotSKIN_nonexposed
  
  
  ggpubr::ggarrange(plotSKIN_sunexposed,
                    plotSKIN_nonexposed,
                    #labels = c("a", "d"),
                    common.legend = T)
  
  file_name <- paste0(getwd(), "/results/_paper/figures/", project_id1, "_", project_id2, "_somatic_mutations_MSR")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  ####################################
  ## STATS
  ####################################

  if (stats) {
    
    df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_D") %>% pull(MSR) %>% summary
    df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_D") %>% pull(MSR) %>% summary
    
    ## 1. Splicing noise is more frequent in introns from non-protein coding transcripts 
    ## than in introns from protein-coding transcripts even after correcting by mean read coverage
    
    ## splicing noise at the donor is more frequent in introns from NPC vs "Skin - Sun Exposed (Lower leg)" 
    wilcox.test(x = df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_D") %>% pull(MSR),
                y = df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_D") %>% pull(MSR),
                #alternative = "greater",
                paired = T,
                correct = T)
    rstatix::wilcox_effsize(data = df_coverage %>% 
                              filter(MSR_type == "MSR_D")%>%
                              mutate(tissue = tissue %>% as.factor()),
                            formula = MSR ~ tissue,
                            paired = T)
    
    
    df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_A") %>% pull(MSR) %>% summary
    df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_A") %>% pull(MSR) %>% summary
    
    ## splicing noise at the acceptor is more frequent in introns from NPC vs PC 
    wilcox.test(x = df_coverage %>% filter(tissue == paste0(cluster_id1, "_", project_id1), MSR_type == "MSR_A") %>% pull(MSR),
                y = df_coverage %>% filter(tissue == paste0(cluster_id2, "_", project_id2), MSR_type == "MSR_A") %>% pull(MSR),
                #alternative = "greater",
                paired = T,
                correct = T)
    rstatix::wilcox_effsize(data = df_coverage %>% 
                              filter(MSR_type == "MSR_A")%>%
                              mutate(tissue = tissue %>% as.factor()),
                            formula = MSR ~ tissue,
                            paired = T)
    
    
    
    # ####################################
    # ## LM
    # ####################################
    # 
    # df_coverage %>%
    #   inner_join(y = skin_introns %>%
    #                dplyr::select(ref_junID,
    #                              gene_id, ref_length, ref_ss5score, ref_ss3score,
    #                              ref_cons5score, ref_cons3score, ref_CDTS5score,
    #                              ref_CDTS3score,
    #                              protein_coding, n_transcripts, gene_width, gene_name),
    #              by = "ref_junID")
    # 
    # 
    # 
    # #########################
    # ## TIDY DATA
    # #########################
    # 
    # idb <- introns %>%
    #   as.data.frame() %>%
    #   dplyr::distinct(ref_junID, .keep_all = T) %>%
    #   dplyr::rename(intron_length = ref_length,
    #                 intron_5ss_score = ref_ss5score,
    #                 intron_3ss_score = ref_ss3score,
    #                 gene_length = gene_width,
    #                 gene_tpm = gene_tpm,
    #                 gene_num_transcripts = n_transcripts,
    #                 CDTS_5ss = ref_CDTS5score,
    #                 CDTS_3ss = ref_CDTS3score,
    #                 protein_coding = protein_coding,
    #                 mean_phastCons20way_5ss = ref_cons5score,
    #                 mean_phastCons20way_3ss = ref_cons3score)
    # 
    # 
    # #########################
    # ## LINEAR MODELS
    # #########################
    # 
    # 
    # fit_donor <- lm(MSR_D ~ 
    #                   intron_length +
    #                   intron_5ss_score + 
    #                   intron_3ss_score +
    #                   gene_length +
    #                   gene_tpm +
    #                   gene_num_transcripts +
    #                   protein_coding +
    #                   CDTS_5ss +
    #                   CDTS_3ss +
    #                   mean_phastCons20way_5ss +
    #                   mean_phastCons20way_3ss, 
    #                 data = idb)
    # 
    # 
    # # fit_donor$coefficients[names(fit_donor$coefficients)[summary(fit_donor)$coefficients[,4] < 0.05]]  
    # 
    # fit_acceptor <- lm(MSR_A ~ 
    #                      intron_length + 
    #                      intron_5ss_score + 
    #                      intron_3ss_score +
    #                      gene_length +
    #                      gene_tpm +
    #                      gene_num_transcripts +
    #                      protein_coding +
    #                      CDTS_5ss +
    #                      CDTS_3ss +
    #                      mean_phastCons20way_5ss +
    #                      mean_phastCons20way_3ss, 
    #                    data = idb)
    # fit_acceptor %>% summary()
    # 
    # model_names <- c("MSR_Donor", "MSR_Acceptor")
    # 
    # 
    # coef_names <- c("Intron Length" = "intron_length",
    #                 "Intron 5'ss MES score" = "intron_5ss_score",
    #                 "Intron 3'ss MES score" = "intron_3ss_score", 
    #                 #"Inter. 5'ss & 3'ss MES score" = "intron_5ss_score:intron_3ss_score", 
    #                 #"Intron Type u2" = "u2_intronTRUE",
    #                 #"Intron ClinVar mutation" = "clinvarTRUE",
    #                 "Gene Length" = "gene_length",
    #                 "Gene TPM" = "gene_tpm",
    #                 "Gene num. transcripts" = "gene_num_transcripts",
    #                 "CDTS 5'ss" = "CDTS_5ss",
    #                 "CDTS 3'ss" = "CDTS_3ss",
    #                 "Cons 5'ss" = "mean_phastCons20way_5ss",
    #                 "Cons 3'ss" = "mean_phastCons20way_3ss",
    #                 "Protein coding" = "protein_coding")
    # 
    # plotLM <- jtools::plot_summs(fit_donor, 
    #                              fit_acceptor,
    #                              #scale = TRUE, 
    #                              #robust = T,
    #                              #n.sd = 2,
    #                              #pvals = TRUE,
    #                              legend.title = "",
    #                              #plot.distributions = TRUE,
    #                              ci_level = 0.95,
    #                              coefs = coef_names,
    #                              colors = c("#35B779FF","#440154FF"),
    #                              model.names = model_names,
    #                              plot.distributions = T,
    #                              rescale.distributions = T,
    #                              
    #                              point.shape = T
    # ) + 
    #   theme_minimal() + 
    #   theme(axis.line = element_line(colour = "black"), 
    #         axis.text.x = element_text(colour = "black", size = "12"),
    #         axis.text.y = element_text(colour = "black", size = "12"),
    #         axis.title = element_text(colour = "black", size = "12"),
    #         legend.text = element_text(colour = "black", size = "12"),
    #         legend.title = element_text(colour = "black", size = "12"),
    #         legend.position = "top",
    #         panel.grid.major.x = element_blank(),
    #         panel.grid.major.y = element_blank(),
    #         axis.ticks = element_line(colour = "black", size = 2)) +  
    #   ylab("Covariates") +
    #   geom_hline(yintercept = seq(from = 0,
    #                               to = length((fit_donor$coefficients %>% names)[-1]) + .5,
    #                               by = 1)) +
    #   guides(colour = guide_legend(#title = NULL,
    #     ncol = 2, 
    #     nrow = 1))     +
    #   scale_x_continuous(breaks = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025),
    #                      labels = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025))
    # 
    # plotLM 
    # plotLM +
    #   annotate("text",
    #            y = c(1.65, 2.65, 3.65, 4.65, 5.65, 6.65, 7.65, 8.65, 9.65, 10.65, 11.65),
    #            x = -0.0025, 
    #            label = c( paste("pval",(summary(fit_donor)$coefficients[,4] %>% 
    #                                       as_tibble() %>%
    #                                       mutate(value = format(value, format = "e", digits = 2)))$value[-1] %>% rev(), sep = ":")),
    #            family = "", 
    #            fontface = 1, 
    #            size=3,
    #            colour = c("#35B779FF") )  +
    #   annotate("text",
    #            y = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 7.25, 8.25, 9.25, 10.25, 11.25),
    #            x = -0.0025, 
    #            label = c( paste("pval",(summary(fit_acceptor)$coefficients[,4] %>% 
    #                                       as_tibble() %>%
    #                                       mutate(value = format(value, format = "e", digits = 2)))$value[-1] %>% rev(), sep = ":")),
    #            family = "", 
    #            fontface = 1, 
    #            size=3,
    #            colour = c("#440154FF") )
    # 
    # file_name <- paste0(getwd(), "/results/_paper/figures/lm_FTCX")
    # ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
    # ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  }
}


#################################
## ENCODE KNOCKDOWN
#################################

## EFFECT-SIZE GRAPH

plot_effect_size_data <- function() {
  
  source(paste0(getwd(), "/code/helper_functions.R"))
  
  ## LOAD RBPS
  
  metadata_filtered <- readr::read_delim("/home/grocamora/RytenLab-Research/Additional_files/ENCODE_files/metadata_filtered.tsv", show_col_types = F)
  
  target_RBPs <- metadata_filtered %>%
    dplyr::filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
    dplyr::pull(target_gene) %>%
    unique()
  
  metadata_RBPs <- metadata_filtered %>% 
    filter(target_gene %in% target_RBPs) %>%
    pivot_longer(c("Splicing_regulation", "Spliceosome", "Exon_junction_complex", "NMD"), names_to = "Category") %>%
    filter(value == 1) %>%
    dplyr::select(-value) %>%
    distinct(target_gene, sample_id, .keep_all = T)
  
  required_clusters <- metadata_RBPs %>% pull(experiment_type) %>% unique

  
  ## LOAD COMMON INTRONS
  
  overwrite = F
  num_cores = 10 
  global_common_introns_path = paste0(getwd(), "/code/variables/global_common_introns_filtered.rds")
  global_common_novel_path = paste0(getwd(), "/code/variables/global_common_novel_filtered.rds")
  
  RBPs_path <- path.expand("/home/grocamora/RytenLab-Research/03-ENCODE_RBP_Automatization/RBPs")
  
  # Load global common introns
  if(!file.exists(global_common_introns_path)){
    global_common_introns <- generateCommonIntronsParallel(target_RBPs, RBPs_path, metadata_RBPs, 
                                                           required_clusters, global_common_introns_path, 
                                                           overwrite, num_cores)
  }else{
    global_common_introns <- readRDS(global_common_introns_path)
  }
  
  # Load global common novel
  if(!file.exists(global_common_novel_path)){
    common_ref_coordinates <- global_common_introns %>% pull(ref_coordinates) %>% unique
    global_common_novel <- generateCommonNovelParallel(target_RBPs, common_ref_coordinates, RBPs_path, 
                                                       metadata_RBPs, required_clusters, global_common_novel_path, 
                                                       overwrite, num_cores)
  }else{
    global_common_novel <- readRDS(global_common_novel_path)
  }
  
  
  ## GET THE EFFECT SIZE
  
  # Variable paths
  overwrite = F
  num_cores = 16 
  MSR_A_tests_path <- paste0(getwd(), "/code/variables/MSR_A_tests.rds")
  MSR_D_tests_path <- paste0(getwd(), "/code/variables/MSR_D_tests.rds")
  
  # Load the MSR tables
  MSR_D <- global_common_introns %>%
    dplyr::select(ref_coordinates, MSR_A = ref_missplicing_ratio_tissue_NA, MSR_D = ref_missplicing_ratio_tissue_ND, target_gene, cluster) %>%
    pivot_wider(id_cols = ref_coordinates, 
                names_from = c("target_gene", "cluster"), 
                values_from = c("MSR_D"))
  
  MSR_A <- global_common_introns %>%
    dplyr::select(ref_coordinates, MSR_A = ref_missplicing_ratio_tissue_NA, MSR_D = ref_missplicing_ratio_tissue_ND, target_gene, cluster) %>%
    pivot_wider(id_cols = ref_coordinates, 
                names_from = c("target_gene", "cluster"), 
                values_from = c("MSR_A"))
  
  ## Execute the wilcox test on MSR_A and MSR_D
  if (!file.exists(paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.rds"))) {
    
    MSR_A_tests <- generateMSRtests(target_RBPs, MSR = MSR_A, file_output = "", overwrite = T, num_cores = num_cores)
    MSR_D_tests <- generateMSRtests(target_RBPs = target_RBPs, 
                                    MSR = MSR_D, 
                                    file_output = "", 
                                    overwrite = T,
                                    num_cores = num_cores)
    
    ## Add the categories
    MSR_A_tests <- addMSRcategories(MSR_A_tests, metadata_RBPs)
    MSR_D_tests <- addMSRcategories(MSR_D_tests, metadata_RBPs)
    
    ## Add bonferroni correction
    MSR_A_tests <- addBonferroniCorrection(MSR_A_tests)
    MSR_D_tests <- addBonferroniCorrection(MSR_D_tests)
    
    
    saveRDS(object = MSR_D_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.rds"))
    xlsx::write.xlsx2(x = MSR_D_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.xlsx"), 
                      sheetName = "MSRD_test", row.names = F, append = T)
    
    saveRDS(object = MSR_A_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.rds"))
    xlsx::write.xlsx2(x = MSR_A_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.xlsx"), 
                     sheetName = "MSRA_test", row.names = F, append = T)
  } else {
    MSR_D_tests <- readRDS(file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.rds"))
    MSR_A_tests <- readRDS(file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.rds"))
  }
  
  # Combine both MSR_A and MSR_D
  MSR_combined = rbind(MSR_A_tests %>% mutate(MSR_type = "MSR_A"), 
                       MSR_D_tests %>% mutate(MSR_type = "MSR_D")) %>%
    filter(p.value.bonferroni <= 0.05)
  
  # Filter the Splicing regulation category
  filter_splicing_regulation <- MSR_combined %>% 
    select(-statistical_test, -H0, -H1) %>% 
    arrange(-effect_size) %>%
    filter(Category == "Splicing_regulation") %>%
    distinct(target_gene, .keep_all = T) %>%
    head(splicing_max) %>%
    pull(target_gene)
  
  MSR_graph_data <- MSR_combined %>%
    filter(target_gene %in% filter_splicing_regulation | Category != "Splicing_regulation") %>%
    arrange(-effect_size) %>%
    mutate(target_gene = factor(target_gene, levels = .$target_gene %>% unique),                                           # Use factors to sort the graph
           Category = factor(Category, levels = c("Splicing_regulation", "Spliceosome", "NMD", "Exon_junction_complex")))  # Use factors to sort the graph
  
  # Plot the graph
  ggplot(MSR_graph_data) + 
    geom_bar(aes(x = target_gene, y = effect_size, fill = MSR_type), 
             stat = "identity", color = "black", linewidth = 0.25, width = 0.80, position = "dodge") + 
    geom_hline(aes(yintercept = 0.1), linewidth = 0.25) + geom_hline(aes(yintercept = 0.3), linewidth = 0.25) + geom_hline(aes(yintercept = 0.5), linewidth = 0.25) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), 
                       breaks = seq(0, 0.6, 0.1),
                       labels = c("0", "Small", "0.2", "Moderate", "0.4", "Large", "0.6")) +
    scale_x_discrete(expand = expansion(add = c(0.7, 0.7))) +
    coord_flip() +
    viridis::scale_fill_viridis(option="viridis", discrete = T, begin = 0.20, end = 0.75, 
                                name = "Splice site:",
                                labels = c("MSR_A" = "Acceptor", "MSR_D" = "Donor"),
                                guide = guide_legend(reverse = T)) +
    labs(x = "Target shRNA knockdown gene", y = "Effect size") + 
    ggforce::facet_col(vars(Category), scales = "free_y", space = "free",
                       labeller = labeller(Category = category_labels)) +
    custom_gg_theme
  
  # Save the graph
  ggsave(file = paste0("images/Effect_size_combined_top", splicing_max, ".png"), width = 183, height = (183/23)*(splicing_max+13), units = "mm", dpi = 300)
}


## AQR AND U2AF1

plot_data_AQR_U2AF1 <- function() {
  
  
  
  ## All plots only contain common introns across RBP projects 
  
  overwrite = F
  num_cores = 10 
  global_common_introns_path = paste0(getwd(), "/code/variables/global_common_introns_filtered.rds")
  global_common_novel_path = paste0(getwd(), "/code/variables/global_common_novel_filtered.rds")
  
  RBPs_path <- path.expand("/home/grocamora/RytenLab-Research/03-ENCODE_RBP_Automatization/RBPs")
  
  
  ################################################
  ## LOAD COMMON INTRONS
  ################################################
  
  # Load global common introns
  if(!file.exists(global_common_introns_path)){
    global_common_introns <- generateCommonIntronsParallel(studied_RBPs = target_RBPs, 
                                                           RBP_path = RBPs_path, 
                                                           metadata = metadata_RBPs, 
                                                           required_clusters = required_clusters, 
                                                           file_output = global_common_introns_path, 
                                                           overwrite = overwrite, 
                                                           num_cores = num_cores)
  }else{
    global_common_introns <- readRDS(global_common_introns_path)
  }
  
  
  # Load global common novel
  if(!file.exists(global_common_novel_path)){
    common_ref_coordinates <- global_common_introns %>% pull(ref_coordinates) %>% unique
    global_common_novel <- generateCommonNovelParallel(target_RBPs, common_ref_coordinates, RBPs_path, 
                                                       metadata_RBPs, required_clusters, global_common_novel_path, 
                                                       overwrite, num_cores)
  }else{
    global_common_novel <- readRDS(global_common_novel_path)
  }
  
  scales::show_col(colours = viridis::viridis_pal(option = "A")(n = 20))
  scales::show_col(colours = viridis::viridis_pal(option = "B")(n = 20))
  scales::show_col(colours = viridis::viridis_pal(option = "C")(n = 20))
  scales::show_col(colours = viridis::viridis_pal(option = "D")(n = 20))
  scales::show_col(colours = viridis::viridis_pal(option = "E")(n = 20))
  scales::show_col(colours = viridis::viridis_pal(option = "F")(n = 20))
  
  #################################################
  ## PLOT AQR AND U2AF1
  #################################################
  
  
  target_RBPs = c("AQR", "U2AF1")
  
  RBP_intron <- global_common_introns %>% filter(target_gene %in% target_RBPs)
  RBP_novel <- global_common_novel %>% filter(target_gene %in% target_RBPs)
  
  
  RBP_intron <- RBP_intron %>%
    mutate(cluster = ifelse(cluster == "case", "gene knockdown", "control"))%>%
    mutate(cluster = factor(cluster, levels = c("control", "gene knockdown")))
  RBP_novel <- RBP_novel %>%
    mutate(cluster = ifelse(cluster == "case", "gene knockdown", "control")) %>%
    mutate(cluster = factor(cluster, levels = c("gene knockdown","control" )))
  
  ## DISTANCES
  
  limit_bp = 30
  
  plot_aqr_u2af1 <- ggplot(RBP_novel %>%
                             filter(novel_type == "novel_acceptor") %>%
                             mutate(novel_type = str_replace(string = novel_type,
                                                             pattern = "_",
                                                             replacement = " ")) %>%
                             filter(abs(distance) < limit_bp)) + 
    geom_histogram(aes(x = distance, fill = cluster),
                   bins = 60, 
                   binwidth = 1, 
                   position = "identity", 
                   alpha = 1, 
                   color = "black", 
                   linewidth = 0.2) +
    #facet_grid(vars(target_gene))+
    scale_x_continuous(expand = expansion(mult = c(0, 0)), 
                       limits = c((limit_bp * -1), limit_bp), 
                       breaks = seq(-limit_bp, limit_bp, length.out = 5)) + 
    
    scale_fill_manual(values = c("#D6445CFF","#333333"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("gene knockdown  ", "control")) +
    
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    # ggforce::facet_col(vars(novel_type), 
    #                    labeller = labeller(novel_type = c("novel_acceptor" = "Novel acceptor"))) +
    facet_wrap(vars(target_gene)) +
    guides(fill = guide_legend(title = "Sample type: ", #title = "Junction category & Strand",
                               #override.aes = list(size = 3),
                               ncol = 2, nrow = 1 )) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "11"),
          axis.text.x = element_text(colour = "black", size = "11",
                                     hjust = 0.1),
          axis.title = element_text(colour = "black", size = "11"),
          strip.text = element_text(colour = "black", size = "11"), 
          legend.text = element_text(colour = "black", size = "11"),
          plot.caption = element_text(colour = "black", size = "11"),
          plot.title = element_text(colour = "black", size = "11"),
          legend.title = element_text(colour = "black", size = "11"),
          legend.position = "none") 
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 60), fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 30),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 30, ymax = 31), fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15, y = 50),  size = 3, label = "intron") +
    theme_void()
  
  
  distances_acceptor <- plot_aqr_u2af1 / (distance_rectangle +  distance_rectangle) + patchwork::plot_layout(heights = c(8, 1))
  distances_acceptor
  
  ggsave(file = paste0(getwd(), "/results/_paper/figures/distance_stacked_AQR_U2AF1.png"), 
         width = 183, height = 110, dpi = 300, units = "mm")
  
  
  
  ## DISTANCES LONGER ------------------------
  
  limit_bp = 200
  target_RBP <- "AQR"
  plot_aqr_u2af1 <- ggplot(RBP_novel %>%
                             filter(target_gene == target_RBP) %>%
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
                   color = "black", 
                   linewidth = 0.1) +
    #facet_grid(vars(target_gene))+
    scale_x_continuous(expand = expansion(mult = c(0, 0)), 
                       limits = c((limit_bp * -1), limit_bp), 
                       breaks = seq(-limit_bp, limit_bp, length.out = 5)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(values = c("#D6445CFF","#333333"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("gene knockdown  ", "control")) +
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    #ggtitle(paste0(target_RBP)) + 
    facet_grid(vars(novel_type), 
               labeller = labeller(novel_type = c("novel donor" = "Novel donor", "novel acceptor" = "Novel acceptor"))) +
    facet_grid(vars(novel_type)) +
    guides(fill = guide_legend(title = "Sample type: ", #title = "Junction category & Strand",
                               override.aes = list(size = 3),
                               ncol = 2, nrow = 1 )) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.text.x = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          legend.text = element_text(colour = "black", size = "12"),
          plot.caption = element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "none") 
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 60), fill = "grey", color = "black") +
    geom_text(aes(x = 100, y = 30),  size = 5, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 30, ymax = 31), fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -100, y = 51),  size = 5, label = "intron") +
    theme_void()
  
  
  plot_aqr_u2af1 / (distance_rectangle) + patchwork::plot_layout(heights = c(8, 1))
  
  ggsave(file = paste0(getwd(), "/results/_paper/figures/distance_stacked_AQR_200bp.png"), 
         width = 200, height = 100, dpi = 300, units = "mm")
  
  
  ## DELTA MES -----------------------------------------
  
  MES_gene <- "AQR"
  
  RBP_merged1 <- RBP_novel %>% 
    filter(target_gene == MES_gene) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(novel_coordinates, ref_junID, cluster, novel_ss3score, target_gene) %>%
    left_join(RBP_intron %>% 
                filter(target_gene == MES_gene) %>%
                dplyr::select(ref_junID, ref_ss5score, ref_ss3score) %>% distinct(),
              by = "ref_junID") %>%
    mutate(delta_ss3score = ref_ss3score - novel_ss3score) %>%
    ungroup()
  
  RBP_delta1 <- RBP_merged1 %>%
    dplyr::select(target_gene, delta_ss3score,cluster ) %>%
    gather(key = "delta_type", value = "delta", delta_ss3score, target_gene) %>%
    mutate(delta = delta %>% as.double()) %>%
    drop_na() %>%
    mutate(delta_type = ifelse(delta_type == "delta_ss5score", "delta MES 5'ss", "delta MES 3'ss")) %>%
    mutate(delta_type = delta_type %>% as.factor())  %>%
    mutate(RBP = MES_gene)
  
  
  MES_gene <- "U2AF1"
  
  RBP_merged2 <- RBP_novel %>% 
    filter(target_gene == MES_gene) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(novel_coordinates, ref_junID, cluster, novel_ss3score, target_gene) %>%
    left_join(RBP_intron %>% 
                filter(target_gene == MES_gene) %>%
                dplyr::select(ref_junID, ref_ss5score, ref_ss3score) %>% distinct(),
              by = "ref_junID") %>%
    mutate(delta_ss3score = ref_ss3score - novel_ss3score) %>%
    ungroup()
  
  RBP_delta2 <- RBP_merged %>%
    dplyr::select(target_gene, delta_ss3score,cluster ) %>%
    gather(key = "delta_type", value = "delta", delta_ss3score, target_gene) %>%
    mutate(delta = delta %>% as.double()) %>%
    drop_na() %>%
    mutate(delta_type = ifelse(delta_type == "delta_ss5score", "delta MES 5'ss", "delta MES 3'ss")) %>%
    mutate(delta_type = delta_type %>% as.factor()) %>%
    mutate(RBP = MES_gene)
  
  
  
  delta_mes <- ggplot(rbind(RBP_delta1,RBP_delta2)) +
    geom_density(aes(x = delta, fill = cluster), 
                 alpha = 0.8, linewidth = 0.3, color = "black", adjust = 1) +
    geom_vline(xintercept = 0) +
    facet_wrap(vars(RBP)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
    
    scale_fill_manual(values = c("#D6445CFF","#333333"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("gene knockdown  ", "control")) + 
    labs(x = "Delta 3' MaxEntScan score", y = "Density") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "11"),
          axis.text.x = element_text(colour = "black", size = "11",
                                     hjust = 0.1),
          axis.title = element_text(colour = "black", size = "11"),
          strip.text = element_text(colour = "black", size = "11"), 
          legend.text = element_text(colour = "black", size = "11"),
          plot.caption = element_text(colour = "black", size = "11"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "11"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Sample type: ",
                               ncol = 2, nrow = 1 )) 
  
  ggpubr::ggarrange(delta_mes,
                    distances_acceptor,
                    ncol = 1,
                    nrow = 2,
                    common.legend = T)
  #custom_gg_theme + theme(legend.key = element_rect(color = "black", linewidth = .05))
  
  ggsave(file = paste0(getwd(),"/results/_paper/figures/delta_MaxEntScan_acceptor.png"), 
         width = 180, height = 170, dpi = 300, units = "mm")
  
  
  
  
  ## This analysis identified a significant reduction (pval) in the strength of the novel 3’ss compared 
  ## to paired annotated sites in case versus control samples, suggesting that the spliceosome was no 
  ## longer able to distinguish between splicing signals accurately at acceptor sites. 
  
  
  ## the alternative has to be greater as a weaker novel splice site will produce higher delta MES values
  wilcox.test(x = RBP_delta %>%
                filter(cluster == "gene knockdown",
                       delta_type == "delta MES 3'ss") %>%
                pull(delta),
              y = RBP_delta %>%
                filter(cluster == "control",
                       delta_type == "delta MES 3'ss") %>%
                pull(delta),
              alternative = "greater")
  
  
}

#################################
## GO ENRICHMENT
#################################

## TODO: not only obtain common introns across age groups, but also obtain common introns across age groups with similar coverage.
## Then analyse changes in splicing noise across age


plot_GO_Enrichment <- function() {
  
  
  project_id <- "BRAIN"
  
  ## Read data.frame containing only common introns across age groups in BRAIN
  folder_results <- paste0(getwd(), "/database/v", gtf_version, "/age_subsampled/")
  
  ## MSR_D
  file_name <- paste0(folder_results, "/results/df_MSRD_common_introns_", project_id,".rds")
  df_MSRD <- readRDS(file = file_name)
  df_MSRD_increasing <- df_MSRD %>%
    filter((`20-39` < `40-59` &
              `40-59` < `60-79`)) 

  ## MSR_A
  file_name <- paste0(folder_results, "/results/df_MSRA_common_introns_", project_id,".rds")
  df_MSRA <- readRDS(file = file_name)
  df_MSRA_increasing <- df_MSRA %>%
    filter((`20-39` < `40-59` &
              `40-59` < `60-79`)) 
  
  
  ## All genes background
  bg_genes <- c(df_MSRD$gene_name, df_MSRA$gene_name) %>% unique()
  
  
  #############################################
  ## MSR_D
  #############################################
  
  ego_MSRD <- clusterProfiler::enrichGO(
                  gene          = df_MSRD_increasing$gene_name %>% unique(),
                  universe      = bg_genes,
                  keyType       = "SYMBOL",
                  OrgDb         = "org.Hs.eg.db", ##Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
                  ont           = "ALL",
                  pAdjustMethod = "bonferroni",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  
  ## PLOT 1
  MSR_D_GO_plot <- barplot(ego_MSRD, 
                           showCategory = 20) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
    xlab("Gene count") +
    facet_grid(ONTOLOGY~., scale = "free") +
    theme_minimal() + 
    theme(text = element_text(colour = "black",size = 12),
          axis.line = element_line(colour = "black"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.text.y = element_text(colour = "black",
                                     vjust = 0.3,
                                     hjust = 1)) +
    guides(fill = guide_legend(title = "pval",
                               label.position = "bottom") ) 
  MSR_D_GO_plot
  
  folder_figures <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/", 
                           main_project, "/figures/")
  dir.create(file.path(folder_figures), recursive = TRUE, showWarnings = T)
  file_name <- paste0(folder_figures, "/GOEnrichment_MSR_D_increasing")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  ## PLOT 2
  
  edox <- clusterProfiler::setReadable(ego_MSRD, OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL')
  edox2 <- enrichplot::pairwise_termsim(edox,  showCategory = 30)
  p1 <- enrichplot::treeplot(edox2, fontsize = 2.5) +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = "10"),
          legend.position = "top") +
  guides(colour = guide_legend(title = "pval",
                               label.position = "bottom") ) 
  p1$layers[[7]]$aes_params$size <- 2.5
  p1
  
  
  file_name <- paste0(folder_figures, "/GOEnrichment_MSR_D_increasing2")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  #############################################
  ## MSR_A
  #############################################
  
  ego_MSRA <- clusterProfiler::enrichGO(
                                    gene          = df_MSRA_increasing$gene_name %>% unique(), 
                                    universe      = bg_genes,
                                    keyType       = "SYMBOL",
                                    OrgDb         = "org.Hs.eg.db", ##Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
                                    ont           = "ALL",
                                    pAdjustMethod = "bonferroni",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)
  
  
  ## PLOT 1
  MSR_A_GO_plot <- barplot(ego_MSRA, showCategory = 20)  +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
    xlab("Gene count") +
    facet_grid(ONTOLOGY~., scale = "free") +
    theme(text = element_text(colour = "black",size = 12), 
          axis.line = element_line(colour = "black"), 
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.text.y = element_text(colour = "black", 
                                     vjust = 0.3,
                                     hjust = 1),
          axis.text.x = element_text(colour = "black")) +
    guides(fill = guide_legend(title = "pval",
                               label.position = "bottom") ) 

  MSR_A_GO_plot
  file_name <- paste0(folder_figures, "/GOEnrichment_MSR_A_increasing")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  

  
  ## PLOT 2
  edox <- clusterProfiler::setReadable(ego_MSRA, 'org.Hs.eg.db', 'SYMBOL')
  edox2 <- enrichplot::pairwise_termsim(edox)
  p1 <- enrichplot::treeplot(edox2, fontsize = 2.5) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title = "pval",
                                 label.position = "bottom") ) 
  p1$layers[[7]]$aes_params$size <- 2.5
  p1
  file_name <- paste0(folder_figures, "/GOEnrichment_MSR_A_increasing2")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
}


##################################
## OTHER
##################################


prepare_gtex_samples <- function() {
  
  
  if (file.exists("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects.rds")) {
    
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects.rds")
    
    all_projects_used <- NULL
    all_clusters_used <- NULL
    all_samples_used <- NULL
    
    for (project_id in all_projects) {
      
      
      clusters <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                                        project_id, "/raw_data/all_clusters.rds"))
      
      print("All clusters saved!")
      
      
      ## SAVE CLUSTERS SP ------------------------------------------------------
      
      clusters_used <- NULL
      
      for (cluster in clusters) {
        
        if (!(cluster %in% c("Brain - Cortex", "Brain - Cerebellum")) &&
            !(project_id %in% c("BREAST", "CERVIX_UTERI", "FALLOPIAN_TUBE", "OVARY",
                                "PROSTATE", "TESTIS", "UTERUS", "VAGINA"))) {
          
          ## Load samples
          samples <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                                           project_id, "/results/base_data/", cluster, "/", 
                                           project_id, "_", cluster, "_samples.rds"))
          
          if (samples %>% length() >= 70) {
            clusters_used <- c(clusters_used, cluster)
            all_projects_used <- c(all_projects_used, project_id) %>% unique()
            all_clusters_used <- c(all_clusters_used, cluster) %>% unique()
            all_samples_used <- c(all_samples_used, samples) %>% unique()
          }
        }
      }
      
      clusters_used <- clusters_used %>% unique()
      print(paste0(Sys.time(), " - ", project_id, " finished!"))
      saveRDS(object = clusters_used,
              file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                            project_id, "/raw_data/all_clusters_used.rds"))
    }
    
    all_projects_used %>% length()
    all_clusters_used %>% length()
    all_samples_used %>% length()
    
    saveRDS(object = all_projects_used %>% unique(),
            file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects_used.rds")
    
  } else {
    all_projects <- c( "ADIPOSE_TISSUE", "ADRENAL_GLAND", "BLADDER", "BLOOD", "BLOOD_VESSEL", "BONE_MARROW", "BRAIN", "BREAST",
                       "CERVIX_UTERI", "COLON", "ESOPHAGUS", "FALLOPIAN_TUBE", "HEART", "KIDNEY", "LIVER", "LUNG",
                       "MUSCLE","NERVE", "OVARY", "PANCREAS", "PITUITARY", "PROSTATE", "SALIVARY_GLAND", "SKIN",
                       "SMALL_INTESTINE", "SPLEEN", "STOMACH", "TESTIS", "THYROID", "UTERUS", "VAGINA") %>% sort()
    saveRDS(object = all_projects,
            file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects.rds")
  }
  saveRDS(object = data.frame(project = c( "ADIPOSE_TISSUE", "ADRENAL_GLAND", "BLADDER", "BLOOD", 
                                           "BLOOD_VESSEL", "BONE_MARROW", "BRAIN", "BREAST",
                                           "CERVIX_UTERI", "COLON", "ESOPHAGUS", "FALLOPIAN_TUBE", 
                                           "HEART", "KIDNEY", "LIVER", "LUNG",
                                           "MUSCLE","NERVE", "OVARY", "PANCREAS", 
                                           "PITUITARY", "PROSTATE", "SALIVARY_GLAND", "SKIN",
                                           "SMALL_INTESTINE", "SPLEEN", "STOMACH", "TESTIS", 
                                           "THYROID", "UTERUS", "VAGINA"),
                              samples = c(1293, 274, 21, 1048,
                                          1398, 204, 2931, 482,
                                          19, 822, 1577, 9,
                                          942, 98, 251, 655,
                                          881, 659, 195, 360,
                                          301, 263, 178, 1940,
                                          193, 255, 384, 410,
                                          706, 159, 173),
                              raw_jxn = c(5768983, 2333266, 822305, 4690652,
                                          5691345, 2507266, 9484210, 3604618,
                                          807104, 4478932, 6106677, 617634,
                                          3991264, 1571657, 2057473, 4660024,
                                          3726485, 4242128, 2279666, 2329118,
                                          3134613, 2706136, 2240651, 7186905,
                                          2415045, 2602164, 2851632, 9220923,
                                          4636674, 2001206, 2140753)),
          file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects_summary.rds")
}


 

##################################
## CALLS
##################################

get_modulo_basic_multiple_tissue()

# compare_tissues_somatic_mutations(project_id1 = "BRAIN",
#                                   cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)",
#                                   project_id2 = "BLOOD",
#                                   cluster_id2 = "Whole Blood")
  