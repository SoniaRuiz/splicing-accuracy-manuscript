library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(DBI)

####################################################
## GENERATE SPLICING DATABASE ######################
####################################################

# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline6_splicing_paper.R")

## CONNECT TO THE DATABASE ------------------------------

gtf_version <- 105
main_project <- "splicing"
database_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v",
                        gtf_version, "/", main_project, "/", main_project, ".sqlite")


getwd()
con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)

## GET FROM MASTER TABLE
query = paste0("SELECT * FROM 'master'")
df_metadata <- dbGetQuery(con, query) 

all_projects <- df_metadata$SRA_project %>% unique
all_projects %>% length() %>% print()




########################################
## FUNCTIONS
########################################


## SECTION 1 ---------------------------------------------

get_database_stats <- function() {
  
  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'master'")
  db_metadata <- dbGetQuery(con, query) %>% as_tibble()
  
  db_metadata %>% nrow() 
  db_metadata %>%
    filter(rin >= 7) %>% nrow()
  db_metadata %>%
    dplyr::count(tissue) %>%
    print(n = 50)
  
  
  query <- paste0("SELECT * from 'intron'")
  db_introns <- dbGetQuery(con, query) 
  db_introns %>% distinct(ref_junID) %>% nrow()
  db_introns %>%
    dplyr::count(misspliced)
  
  db_introns %>%
    dplyr::count(gene_id)
  
  
  query <- paste0("SELECT * from 'novel'")
  novel <- dbGetQuery(con, query) 
  novel %>% distinct(novel_junID) %>% nrow()
  novel %>% dplyr::count(novel_type)
  
  
}

get_contamination_rates_all_tissues <- function () {

  
  #############################
  ## CONNECT TO THE DATABASE
  #############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  if ( !file.exists( paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                            main_project, "/results/contamination_rates.rds") )) {
    
    
    df_contamination <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      # project_id <- all_projects[12]
      # project_id <-"KIDNEY"
      
      print(paste0(Sys.time(), " - ", project_id))
      
      
      if (all_tissues) {
        all_clusters <- df_metadata %>%
          filter(SRA_project == project_id) %>%
          distinct(cluster) %>%
          pull()
        versions <- c("97")
        dates <- c("26-May-2019")
      } else {
        all_clusters <- "Brain - Frontal Cortex (BA9)"
        versions <- c("76", "81", "90", "97", "104")
        dates <- c("18-Jul-2014", "07-Jul-2015", "28-Jul-2017", "26-May-2019", "19-Mar-2021")
      }
      
      map_df(all_clusters, function(cluster) {
        
        # cluster <- all_clusters[1]
        print(paste0(cluster))
        
        map_df(versions, function(version) {
          
          # version <- versions[1]
          print(paste0("v", version))
          
          if (file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                                 project_id, "/v", version, "/", main_project, "_project/raw_data/",
                                 project_id, "_", cluster, "_all_split_reads_sample_tidy.rds"))) {
            
            data_old <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                               project_id, "/v", version, "/", main_project, "_project/raw_data/",
                               project_id, "_", cluster, "_all_split_reads_sample_tidy.rds")
            
            ## ENSEMBL v97
            df_old <-  readRDS(file = data_old) %>%
              as_tibble()
            db_introns_old <- df_old %>% 
              filter(type == "annotated")
            db_novel_old <- df_old %>%
              filter(type %in% c("novel_donor", "novel_acceptor"))
            
            ## ENSEMBL v105
            df_new <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                                            project_id, "/v105/", main_project, "_project/raw_data/",
                                            project_id, "_", cluster, "_all_split_reads_sample_tidy.rds")) %>%
              as_tibble()
            db_introns_new <- df_new %>%
              filter(type == "annotated")
            db_novel_new <- df_new %>%
              filter(type %in% c("novel_donor", "novel_acceptor"))
            
            
            rm(df_old)
            rm(df_new)
            
            
            in_annotation <- db_introns_new %>%
              filter(junID %in% db_novel_old$junID) %>% 
              distinct(junID, .keep_all = T)
            in_annotation %>% nrow() %>% print()
            
            out_annotation <- db_novel_new %>%
              filter(junID %in% db_introns_old$junID) %>% 
              distinct(junID, .keep_all = T) 
            out_annotation %>% nrow() %>% print()
            
            label <- NULL
            
            if (!all_tissues) {
              if (version == "76") {
                label <- paste0("Ensembl_v", version, " (", dates[1], ")")
              } else if (version == "81") {
                label <- paste0("Ensembl_v", version, " (", dates[2], ")")
              } else if (version == "90") {
                label <- paste0("Ensembl_v", version, " (", dates[3], ")")
              } else {
                label <- paste0("Ensembl_v", version, " (", dates[4], ")")
              } 
            } else {
              label <- paste0("Ensembl_v", version, " (", dates, ")")
            }
            
            
            return(data.frame(tissue = cluster,
                              contamination_rates = label,
                              in_annotation = (in_annotation %>% nrow() * 100) / db_novel_old %>% nrow(),
                              out_annotation = (out_annotation %>% nrow() * 100) / db_introns_old %>% nrow()))
          } else {
            return(NULL)
          }
          
        })
      })
    })
    
    
    contamination_path <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                 main_project, "/results/contamination_rates.rds")
    dir.create(file.path(contamination_path), recursive = TRUE, showWarnings = T)
    saveRDS(object = df_contamination,
            file = paste0(contamination_path, "/contamination_rates.rds"))
    
    
  } else {
    df_contamination <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                              main_project, "/results/contamination_rates.rds"))
  }
  
  
  df_contamination_tidy <- df_contamination %>% 
    dplyr::select(in_annotation, out_annotation, tissue) %>%
    tidyr::gather(key = "type", value = "prop", -tissue ) 
  
  df_contamination_tidy$type = factor( df_contamination_tidy$type, 
                                       levels = c( "in_annotation", "out_annotation" ) )
  
  df_contamination_tidy = df_contamination_tidy %>% 
    ungroup() %>%
    arrange(type , prop) %>%
    mutate(tissue = fct_inorder(tissue))
  
  # colours <- ifelse(str_detect(string = as.factor(df_contamination_tidy$tissue), pattern = "Brain"), "red", "black")
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data = df_contamination_tidy) +
    #geom_bar(data = novel_104, mapping = aes(x = novel_p_individuals, y = ..count..), stat = "count", color = "black") +
    geom_bar(mapping = aes(x = tissue, y = prop, fill = type), 
             stat = "identity", position = "identity") + 
    coord_flip() +
    ylab("contamination rates (%)") +
    xlab(" ") +
    #ggtitle(paste0("Contamination rates in Ensembl v97 vs v105")) +
    #scale_y_continuous(breaks = c(0,0.3,0.6,0.9)) +
    #scale_x_binned() +
    #scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100))+
    scale_fill_manual(values = c("#440154FF","#FDE725FF"),
                      breaks = c("in_annotation", "out_annotation"),
                      labels = c("introns entring annotation", "introns exiting annotation")) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          axis.text.x = element_text(colour = "black", size = "10"),
          legend.position = "top",
          axis.text.y = element_text(#color = colours,
            #                          angle = 70, 
            vjust = 0.5,
            hjust = 1)
    ) +
    guides(fill = guide_legend(title = NULL, ncol = 1, nrow = 2)) %>%
    return()
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/contamination_rates_all_tissues")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  ## Stats
  df_contamination_tidy %>%
    group_by(type) %>%
    mutate(average_prop = prop %>% mean()) %>%
    ungroup() %>%
    arrange(type, desc(prop)) %>%
    print(n=90)
  
  

  
  rm(df_contamination)
  gc()
}

get_contamination_rates_FCTX <- function (all_tissues = F) {
  
  #############################
  ## CONNECT TO THE DATABASE
  #############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  
  project_id <- "BRAIN"
  cluster <- "Brain - Frontal Cortex (BA9)"
  versions <- c("76", "81", "90", "97", "104")
  dates <- c("Jul-2014", "Jul-2015", "Jul-2017", "May-2019", "Mar-2021")
  
  if ( !file.exists( paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                            project_id, "/v105/", main_project, 
                            "_project/results/contamination/contamination_rates.rds") ) ) {
    
    df_contamination <- map_df(versions, function(version) {
          
      # version <- versions[1]
      print(paste0("v", version))
      
      if (file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                             project_id, "/v", version, "/", main_project, "_project/raw_data/",
                             project_id, "_", cluster, "_all_split_reads_sample_tidy.rds"))) {
        
        data_old <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                           project_id, "/v", version, "/", main_project, "_project/raw_data/",
                           project_id, "_", cluster, "_all_split_reads_sample_tidy.rds")
        
        ## ENSEMBL v97
        df_old <-  readRDS(file = data_old) %>%
          as_tibble()
        db_introns_old <- df_old %>% 
          filter(type == "annotated")
        db_novel_old <- df_old %>%
          filter(type %in% c("novel_donor", "novel_acceptor"))
        
        ## ENSEMBL v105
        df_new <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                                        project_id, "/v105/", main_project, "_project/raw_data/",
                                        project_id, "_", cluster, "_all_split_reads_sample_tidy.rds")) %>%
          as_tibble()
        db_introns_new <- df_new %>%
          filter(type == "annotated")
        db_novel_new <- df_new %>%
          filter(type %in% c("novel_donor", "novel_acceptor"))
        
        
        rm(df_old)
        rm(df_new)
        
        
        in_annotation <- db_introns_new %>%
          filter(junID %in% db_novel_old$junID) %>% 
          distinct(junID, .keep_all = T)
        in_annotation %>% nrow() %>% print()
        
        out_annotation <- db_novel_new %>%
          filter(junID %in% db_introns_old$junID) %>% 
          distinct(junID, .keep_all = T) 
        out_annotation %>% nrow() %>% print()
        
        label <- NULL
        
     
        if (version == "76") {
          label <- paste0("Ensembl v", version, "\n(", dates[1], ")")
        } else if (version == "81") {
          label <- paste0("Ensembl v", version, "\n(", dates[2], ")")
        } else if (version == "90") {
          label <- paste0("Ensembl v", version, "\n(", dates[3], ")")
        } else if (version == "97") {
          label <- paste0("Ensembl v", version, "\n(", dates[4], ")")
        } else {
          label <- paste0("Ensembl v", version, "\n(", dates[5], ")")
        }
        
        
        return(data.frame(tissue = cluster,
                          contamination_rates = label,
                          in_annotation = (in_annotation %>% nrow() * 100) / db_novel_old %>% nrow(),
                          out_annotation = (out_annotation %>% nrow() * 100) / db_introns_old %>% nrow()))
      } else {
        return(NULL)
      }
          

    })
    
    
    contamination_path <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                                 project_id, "/v105/", main_project, "_project/results/contamination")
    dir.create(file.path(contamination_path), recursive = TRUE, showWarnings = T)
    saveRDS(object = df_contamination,
            file = paste0(contamination_path, "/contamination_rates.rds"))
    
    
  } else {
    df_contamination <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                                              project_id, "/v105/", main_project, 
                                              "_project/results/contamination/contamination_rates.rds"))
  }
  
  
  df_contamination$contamination_rates = factor(df_contamination$contamination_rates, 
                                                levels = c("Ensembl v76\n(Jul-2014)",
                                                           "Ensembl v81\n(Jul-2015)",
                                                           "Ensembl v90\n(Jul-2017)",
                                                           "Ensembl v97\n(May-2019)",
                                                           "Ensembl v104\n(Mar-2021)"))


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
          axis.text.y = element_text(#angle = 50, 
            vjust = 0.5,
            hjust = 1)) +
    guides(fill = "none")
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/contamination_rates_FCTX_lineplot")
  ggplot2::ggsave(paste0(file_name, ".svg"),  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"),  width = 183, height = 183, units = "mm", dpi = 300)
  
  
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
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/contamination_rates_FCTX")
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

  if ( exists("df_proportions") )  {
    saveRDS(object = df_proportions,
            file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                          main_project, "/results/unique_donor_acceptor_jxn.rds"))
  } else {
    df_proportions <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                            main_project, "/results/unique_donor_acceptor_jxn.rds"))
    
    df_proportions %>%
      arrange(desc(annotated_prop))
  }
  
  
  
  ###################
  ## STATS
  ###################
  
  
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
  
  df_proportions_violin <- df_proportions %>%
    filter(type != "annotated_intron") %>%
    mutate(prop = prop * 100)
  
  df_proportions_violin$type = factor(df_proportions_violin$type, 
                                      levels = c("acceptor", "donor"))
   

  ggplot(df_proportions_violin, aes(type, prop, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point() +
    geom_line( aes(group = tissue) )  +
    theme_light() +
    ylab("% unique novel junctions") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_x_discrete( breaks = c( "acceptor", "donor"),# "annotated_intron"),
                        labels = c( "novel acceptor", "novel donor")) +
    scale_fill_manual(values = c( "#4393c3", "#f4a582"),
                      breaks = c( "acceptor", "donor"),# "annotated_intron"),
                      labels = c( "novel acceptor", "novel donor")) + #, "annotated intron")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  
  
  ## Save the figure 3
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/percent_unique_donor_acceptor_violin")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), 
                  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
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
  
  if ( exists("df_mean_counts") )  {
     
    saveRDS(object = df_mean_counts,
            file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                          main_project, "/results/unique_donor_acceptor_reads.rds"))
    
  } else {
    
    df_mean_counts <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                            main_project, "/results/unique_donor_acceptor_reads.rds")) %>%
      as_tibble()
  }
  
  ###################
  ## SOME STATS
  ###################
  
  df_mean_counts %>%
    filter(type == "annotated") %>%
    pull(prop_counts ) %>% 
    mean()
  df_mean_counts %>%
    filter(type == "acceptor") %>%
    pull(prop_counts ) %>% 
    mean()
  df_mean_counts %>%
    filter(type == "donor") %>%
    pull(prop_counts ) %>% 
    mean()
  
  
  ######################
  ## VIOLIN PLOT
  ######################
  
  df_mean_counts_violin <- df_mean_counts %>%
    filter(type != "annotated") 
  
  df_mean_counts_violin$type = factor(df_mean_counts_violin$type, 
                                      levels = c("acceptor", "donor"))
  
  
  ggplot(df_mean_counts_violin,# %>%
           #filter(tissue %in% c("Muscle - Skeletal", "Whole Blood","Brain - Cerebellum")), 
         aes(type, prop_counts, fill = type)) + 
    geom_violin(trim = FALSE) +
    geom_point() +
    geom_line( aes(group = tissue) )  +
    theme_light() +
    ylab("% read counts across samples") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          text = element_text(colour = "black", size = 12),
          legend.position = "top") +
    scale_x_discrete( breaks = c( "acceptor", "donor"),# "annotated_intron"),
                      labels = c( "novel acceptor", "novel donor")) +
    scale_fill_manual(values = c( "#4393c3", "#f4a582"),
                      breaks = c( "acceptor", "donor"),# "annotated_intron"),
                      labels = c( "novel acceptor", "novel donor")) + #, "annotated intron")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/unique_donor_acceptor_reads_violin")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
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
  
  tables <- c("Brain - Frontal Cortex (BA9)_BRAIN_misspliced",
              "Brain - Frontal Cortex (BA9)_BRAIN_nevermisspliced")
  
  df_metadata %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)") %>%
    as_tibble()
  
  annotated_sum_count <- NULL
  novel_sum_count <- NULL
  novel_donor_sum_count <- NULL
  novel_acceptor_sum_count <- NULL
  
  query <- paste0("SELECT novel_junID, novel_type from 'novel'")
  db_novel <- dbGetQuery(con, query) %>% as_tibble()
  
  for (table in tables) {
    
    if ( !(table %in% c("intron", "novel", "gene", "mane", "master")) ) {
      
      print(paste0("Table: '", table, "'!"))
      
      query <- paste0("SELECT * from '", table,"'")
      db <- dbGetQuery(con, query) %>% as_tibble()
      
      annotated_sum_count <- c(annotated_sum_count, 
                               db$ref_sum_counts)
      novel_sum_count <- c(novel_sum_count, 
                           db$novel_sum_counts)
      
      if ( table == "Brain - Frontal Cortex (BA9)_BRAIN_misspliced" ) {
        novel_donor_sum_count <- c(novel_donor_sum_count, db %>%
                                     filter(novel_junID %in% (db_novel %>% 
                                                                filter(novel_type == "novel_donor") %>%
                                                                pull(novel_junID))) %>%
                                     pull(novel_sum_counts))
        
        novel_acceptor_sum_count <- c(novel_acceptor_sum_count, db %>%
                                        filter(novel_junID %in% (db_novel %>% 
                                                                   filter(novel_type == "novel_acceptor") %>%
                                                                   pull(novel_junID))) %>%
                                        pull(novel_sum_counts))
      }

      
    }
  }
  annotated_sum_count %>% summary()
  novel_donor_sum_count %>% summary()
  novel_acceptor_sum_count %>% summary()
  
 
  
}

get_modulo_basic_multiple_tissue <- function() {
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  
  
  file_path <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/results/df_modulo_basic_tissues_100bpfilter.rds")
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
        
        samples <- df_metadata %>%
          dplyr::count(cluster) %>%
          filter(cluster == cluster_id) %>% 
          pull(n)
        
        
        
        
        query <- paste0("SELECT novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
        introns <- dbGetQuery(con, query) %>% as_tibble()
        query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                        paste(introns$novel_junID, collapse = ","),")")
        introns <- introns %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "novel_junID") %>% 
          as_tibble() 
        query <- paste0("SELECT ref_junID, protein_coding FROM 'intron' WHERE ref_junID IN (",
                        paste(introns$ref_junID, collapse = ","),")  AND MANE = 1")
        introns <- introns %>%
          left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                    by = "ref_junID") %>% 
          as_tibble() 
        
        
        df_novel_tidy <- introns %>%
          filter(abs(distance) <= 100) %>%
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
  ## VIOLIN PLOT
  ################
  library(ggridges)
  
  df_modulo_tissues <- df_modulo_tissues %>%
    mutate(freq = freq * 100)
  
  ggplot(df_modulo_tissues, aes(x = freq, y = modulo, fill = modulo)) +
    geom_density_ridges_gradient() +
    scale_fill_viridis_d(name = "Percentage Mod3:", option = "D") +
    #labs(title = 'Modulo3 of the distances to annotated splice sites')+
    ylab("Modulo3") +
    xlab("Percentage of novel junctions") +
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
         
          plot.caption = element_text(colour = "black",size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo3 type: "))
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/modulo3_alltissues_density")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
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
    dplyr::group_by(novel_type, type) %>%
    mutate(modulo_mean = modulo %>% mean) %>%
    distinct(modulo_mean) %>%
    arrange(novel_type, desc(modulo_mean))
    
  
}

get_maxentscan_score <- function() {
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
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
      
      query <- paste0("SELECT ref_junID, ref_type FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_type FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
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
  
  
  if ( exists("df_mes") )  {
    saveRDS(object = df_mes %>% as_tibble(),
            file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                          main_project, "/results/df_mes.rds"))
  } else {
    df_mes <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                            main_project, "/results/df_mes.rds")) %>%
      as_tibble()
  }
  
  
  
  
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
                               ncol = 1, nrow = 2))
  
  ss5plot
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/panel2_ss5plot")
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
    guides(fill = guide_legend(title = element_blank(), ncol = 1, nrow = 2))
  
  
  ss3plot
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/panel2_ss3plot")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ## Delta MES -----------------------------------------------------------
  
  df_delta <- df_mes %>% 
    dplyr::select(diff_ss5score, diff_ss3score) %>%
    tidyr::gather(key = "type", value = "mean") %>%
    dplyr::filter(mean != 0) %>%
    mutate(type = type %>% as.factor()) %>%
    mutate(type = relevel(type, ref = "diff_ss5score"))
  
  ## Plot
  deltaplot <- ggplot(data = df_delta) +
    geom_density(mapping = aes(x = mean, fill = type), 
                 alpha = 0.8) +
    geom_vline(xintercept = 0) +
    xlab("Delta MaxEntScan score") +
    theme_light() +
    scale_color_viridis_d() +
    scale_fill_manual(values =  c("#EF7F4FFF","#0D0887FF"),
                      breaks = c("diff_ss5score", "diff_ss3score"),
                      labels = c("Delta 5'ss  ", "Delta 3'ss")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "10"),
          axis.title = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "10"),
          legend.title = element_text(size = "10"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 1, nrow = 2))
  
  deltaplot
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/panel2_deltaplot")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ##############
  ## COMBO
  ##############
  
  ggpubr::ggarrange( ss5plot, ss3plot, deltaplot, 
                     labels = c("c", "d", "e"),
                     ncol = 3, nrow = 1)
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/MES_combo")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 360, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 360, height = 183, units = "mm", dpi = 300)
 
}

get_distances <- function() {
  
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  limit_bp <- 30
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  query <- paste0("SELECT novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                  paste(introns$novel_junID, collapse = ","),")")
  introns <- merge(x = introns,
                   y = dbGetQuery(con, query) %>% as_tibble(),
                   by = "novel_junID",
                   all.x = T) %>% 
    as_tibble() 
  
  
 
  #############################
  ## PLOT THE DISTANCES GRAPH 
  #############################
  
  
  df_novel_tidy <- introns %>%
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
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
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
    geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
    theme_void()


  plot_all <- plot_all / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/distances_FCTX")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  

  ###############################
  ## PROTEIN CODING
  ###############################
 
  query <- paste0("SELECT ref_junID, protein_coding, lncRNA FROM 'intron'")
  df_master_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  df_novel_tidy <- introns %>%
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
  
  df_novel_tidy %>% filter(type_PC == "protein coding (PC)") %>% distinct(ref_junID) %>% nrow()
  
  
  plot_PC <- ggplot(data = df_novel_tidy %>% 
                      filter(type_PC == "protein coding (PC)")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggplot2::facet_grid(vars(novel_type)) +
    ggtitle(paste0("Protein coding (PC)")) +
    #ylim(y_axes) +
    
    xlab("Distances (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
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
    geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
    geom_rect(aes(xmin = (limit_bp * -1), xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
    theme_void()
  
  
  plot_PC <- plot_PC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/FCTX_distancesPC")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  df_novel_tidy %>% filter(lncRNA == 100) %>% distinct(ref_junID) %>% nrow()
  
  plot_NPC <- ggplot(data = df_novel_tidy %>% 
                       filter(type_PC == "non PC")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggplot2::facet_grid(vars(novel_type)) +
    ggtitle(paste0("Non-protein coding (NPC)\n")) +
    #ylim(y_axes) +
    xlab("Distances (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
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
  
  
  plot_NPC <- plot_NPC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/FCTX_distancesNPC_lncRNA")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  #######################################
  ## COMBINE ALL 3 PLOTS
  #######################################
  
  plot <- ggpubr::ggarrange(plot_all,
                            plot_PC,
                            plot_NPC,
                            common.legend = T,
                            labels = c("f", "g", "h"),
                            align = "v",
                            ncol = 3,
                            nrow = 1)


  plot

  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/distances_FCTX_all")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 369, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 369, height = 183, units = "mm", dpi = 300)


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
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/MSR_FCTX_subsampling")
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
                                      levels = c("PC", "non PC"))
  
  df_introns_biotype <- df_introns_biotype %>%
    dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)
  
  
  df_introns_biotype$MSR_type = factor(df_introns_biotype$MSR_type, 
                                      levels = c("MSR_A", "MSR_D"))
  
  
  print(paste0(Sys.time(), " - ", df_introns_biotype %>% nrow(), " - introns after tidying!"))
  
  ###########################
  ## PLOT MIS-SPLICING RATIO 
  ###########################
  
  plotMSR_PC <- ggplot(data = df_introns_biotype %>% 
                         filter(biotype == "PC")) + 
    #geom_density(aes(x = MSR, fill = MSR_type), alpha = 0.8) +
    #ggforce::facet_zoom(xlim=c(0,0.02)) +
    geom_bar(aes(x = MSR,  fill = MSR_type), 
             stat = "bin", 
             bins = 30,
             position = "dodge")  +
    ggtitle("Protein-coding") +
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
  
  
  
  ## Non-protein-coding
  plotMSR_NPC <- ggplot(data = df_introns_biotype %>% 
                          filter(biotype == "non PC")) + 
    #geom_density(aes(x = MSR, fill = MSR_type), alpha = 0.8) +
    geom_bar(aes(x = MSR,  fill = MSR_type), 
             stat = "bin", 
             bins = 30,
             position = "dodge")  +
    #ggforce::facet_zoom(xlim=c(0,0.02)) +
    #facet_wrap(vars(type_PC)) +
    ggtitle("Non protein-coding") +
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
  plotMSR_NPC
  
  
  ggpubr::ggarrange(plotMSR_PC,
                    plotMSR_NPC,
                    #labels = c("a", "d"),
                    common.legend = T)
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/MSR_FCTX")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  ####################################
  ## STATS
  ####################################
  
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
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/results/MSR_all_tissues.rds")
  
  
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
        
        query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
        introns <- dbGetQuery(con, query) %>% as_tibble()
        query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", project_id, "_misspliced'")
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
  
  
  query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
  
  query <- paste0("SELECT ref_junID, ref_length, ref_ss5score, ref_ss3score, 
                  ref_cons5score, ref_cons3score, ref_CDTS5score, ref_CDTS3score, 
                  protein_coding
                  FROM 'intron' 
                  WHERE ref_junID IN (",
                  paste(introns$ref_junID, collapse = ","),")")
  introns <- introns %>%
    left_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = "ref_junID") %>% 
    as_tibble() 
  
  
  query <- paste0("SELECT id, n_transcripts, gene_width, gene_name 
                  FROM 'gene' WHERE id IN (",
                  paste(introns$gene_id, collapse = ","),")")
  introns <- introns %>%
    left_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = c("gene_id" = "id")) %>% 
    as_tibble() 
  
  
  #########################
  ## TIDY DATA
  #########################
  
  idb <- introns %>%
    as.data.frame() %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    dplyr::rename(intron_length = ref_length,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score,
                  gene_length = gene_width,
                  gene_tpm = gene_tpm,
                  gene_num_transcripts = n_transcripts,
                  CDTS_5ss = ref_CDTS5score,
                  CDTS_3ss = ref_CDTS3score,
                  protein_coding = protein_coding,
                  mean_phastCons20way_5ss = ref_cons5score,
                  mean_phastCons20way_3ss = ref_cons3score)
  
  
  #########################
  ## LINEAR MODELS
  #########################
  
  
  fit_donor <- lm(MSR_D ~ 
                    intron_length +
                    intron_5ss_score + 
                    intron_3ss_score +
                    gene_length +
                    gene_tpm +
                    gene_num_transcripts +
                    protein_coding +
                    CDTS_5ss +
                    CDTS_3ss +
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss, 
                  data = idb)
 
  
  # fit_donor$coefficients[names(fit_donor$coefficients)[summary(fit_donor)$coefficients[,4] < 0.05]]  
    
  fit_acceptor <- lm(MSR_A ~ 
                       intron_length + 
                       intron_5ss_score + 
                       intron_3ss_score +
                       gene_length +
                       gene_tpm +
                       gene_num_transcripts +
                       protein_coding +
                       CDTS_5ss +
                       CDTS_3ss +
                       mean_phastCons20way_5ss +
                       mean_phastCons20way_3ss, 
                     data = idb)
  fit_acceptor %>% summary()
  
  model_names <- c("MSR_Donor", "MSR_Acceptor")
  
  
  coef_names <- c("Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  #"Inter. 5'ss & 3'ss MES score" = "intron_5ss_score:intron_3ss_score", 
                  #"Intron Type u2" = "u2_intronTRUE",
                  #"Intron ClinVar mutation" = "clinvarTRUE",
                  "Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  "CDTS 5'ss" = "CDTS_5ss",
                  "CDTS 3'ss" = "CDTS_3ss",
                  "Cons 5'ss" = "mean_phastCons20way_5ss",
                  "Cons 3'ss" = "mean_phastCons20way_3ss",
                  "Protein coding" = "protein_coding")
  
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
                             colors = c("#35B779FF","#440154FF"),
                             model.names = model_names,
                             plot.distributions = T,
                             rescale.distributions = T,
                   
                             point.shape = T
                             ) + 
    theme_minimal() + 
    theme(axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = "black", size = "12"),
          axis.text.y = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2)) +  
    ylab("Covariates") +
    geom_hline(yintercept = seq(from = 0,
                                to = length((fit_donor$coefficients %>% names)[-1]) + .5,
                                by = 1)) +
    guides(colour = guide_legend(#title = NULL,
                                 ncol = 2, 
                                 nrow = 1))     +
    scale_x_continuous(breaks = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025),
                       labels = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025))

  plotLM 
  plotLM +
  annotate("text",
           y = c(1.65, 2.65, 3.65, 4.65, 5.65, 6.65, 7.65, 8.65, 9.65, 10.65, 11.65),
           x = -0.0025, 
           label = c( paste("pval",(summary(fit_donor)$coefficients[,4] %>% 
                                      as_tibble() %>%
                                      mutate(value = format(value, format = "e", digits = 2)))$value[-1] %>% rev(), sep = ":")),
           family = "", 
           fontface = 1, 
           size=3,
           colour = c("#35B779FF") )  +
    annotate("text",
             y = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 7.25, 8.25, 9.25, 10.25, 11.25),
             x = -0.0025, 
             label = c( paste("pval",(summary(fit_acceptor)$coefficients[,4] %>% 
                                        as_tibble() %>%
                                        mutate(value = format(value, format = "e", digits = 2)))$value[-1] %>% rev(), sep = ":")),
             family = "", 
             fontface = 1, 
             size=3,
             colour = c("#440154FF") )
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/lm_FTCX")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  # if (save_results) {
  #   plotLM
  #   filename <- paste0(folder_name, "/images/lm_MSR.png")
  #   ggplot2::ggsave(filename = filename, 
  #                   width = 183, height = 183, units = "mm", dpi = 300)
  # } else {
  #   plot %>% return()
  # }
  
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
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/results/common_introns.rds")
  saveRDS(object = introns, file = file_name)
  
}

get_estimate_variance_across_tissues <- function() {
  
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  ## LOAD COMMON INTRONS
  common_introns <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                          main_project, "/results/common_introns.rds"))

  ################################
  ## QUERY THE DATABASE PER TISSUE
  ################################

  df_estimates <- map_df(all_projects, function(project_id) {
    
    print(paste0(Sys.time(), " - ", project_id))
    
    #############################
    ## GET THE CLUSTERS
    #############################
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(all_clusters, function(cluster_id) {
 
      print(cluster_id)
      
      query <- paste0("SELECT ref_junID, MSR_D, MSR_A, gene_tpm FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT ref_junID, MSR_D, MSR_A, gene_tpm FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      query <- paste0("SELECT ref_junID, ref_length AS intron_length, ref_ss5score AS intron_5ss_score, 
                      ref_ss3score AS intron_3ss_score, 
                      ref_cons5score AS mean_phastCons20way_5ss, 
                      ref_cons3score AS mean_phastCons20way_3ss, 
                      ref_CDTS5score AS CDTS_5ss, 
                      ref_CDTS3score AS CDTS_3ss, 
                      protein_coding, gene_id 
                      FROM 'intron' WHERE ref_junID IN (",
                      paste(introns$ref_junID, collapse = ","),")")
      introns <- merge(x = introns,
                       y = dbGetQuery(con, query) %>% as_tibble(),
                       by = "ref_junID",
                       all.x = T) %>% 
        as_tibble() 
      
      
      query <- paste0("SELECT id, n_transcripts AS gene_num_transcripts, 
                      gene_width AS gene_length, gene_name FROM 'gene' WHERE id IN (",
                      paste(introns$gene_id, collapse = ","),")")
      introns <- merge(x = introns,
                       y = dbGetQuery(con, query) %>% as_tibble(),
                       by.x = "gene_id",
                       by.y = "id",
                       all.x = T) %>% 
        as_tibble() 
      
      
      model <- lm(MSR_D ~ 
                    intron_length + 
                    intron_5ss_score *
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
      
      
      #model %>% summary() %>% print()
      #MSR_Donor_list[[tissue]] <- model
      
      
      #ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
      #model$coefficients[-ind_sign] <- 0
      
      MSR_Donor <- data.frame(tissue = cluster_id,
                              intron_length = model$coefficients["intron_length"] %>% unname(),
                              intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                              intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                              gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                              gene_length = model$coefficients["gene_length"] %>% unname(),
                              #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),
                              protein_coding = model$coefficients["protein_coding"] %>% unname(),
                              gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                              CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                              CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                              mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                              mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname(),
                              pval = ((model %>% summary())$coefficients  %>% as.data.frame())[,4])
      
      
      ## Acceptor
      model <- lm(MSR_A ~
                    intron_length + 
                    intron_5ss_score *
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
      #model %>% summary() %>% print()
      #MSR_Acceptor_list[[tissue]] <- model
      
      #ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
      #model$coefficients[-ind_sign] <- 0
      
      MSR_Acceptor <- data.frame(tissue = cluster_id,
                                 intron_length = model$coefficients["intron_length"] %>% unname(),
                                 intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                                 intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                                 gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                                 gene_length = model$coefficients["gene_length"] %>% unname(),
                                 #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),
                                 protein_coding = model$coefficients["protein_coding"] %>% unname(),
                                 gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                                 CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                                 CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                                 mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                                 mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname(),
                                 pval = ((model %>% summary())$coefficients  %>% as.data.frame())[,4])
      
      return(rbind(MSR_Donor %>% mutate(type = "MSR_Donor"), 
                   MSR_Acceptor %>% mutate(type = "MSR_Acceptor")))
      
    })
  })
  
  
  #############################
  ## SAVE RESULTS
  #############################
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/results/variance_estimate_tissues.rds")
  saveRDS(object = df_estimates, file = file_name)
  
}

plot_estimate_variance_across_tissues <- function() {
  

  df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                                       main_project, "/results/variance_estimate_tissues.rds")) %>%
    group_by(type) %>%
    mutate(pval_corrected =p.adjust(pval)) %>%
    filter(pval_corrected < 0.05) %>%
    dplyr::select(-c(pval, pval_corrected))
 
  MSR_Donor <- df_estimate %>%
    filter(type == "MSR_Donor") 
  MSR_Acceptor <- df_estimate %>%
    filter(type == "MSR_Acceptor")
  
  #############################
  ## PLOT - ESTIMATE VARIANCE
  #############################
  
  ## DONOR
  
  MSR_Donor_tidy <- MSR_Donor %>% 
    dplyr::rename("Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  "Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  "CDTS 5'ss" = "CDTS_5ss",
                  "CDTS 3'ss" = "CDTS_3ss",
                  "Conservation 5'ss" = "mean_phastCons20way_5ss",
                  "Conservation 3'ss" = "mean_phastCons20way_3ss", 
                  "Protein coding" = "protein_coding") %>%
    gather(tissue, feature, -type, -tissue) %>% as_tibble() %>%
    mutate(type = "MSR_Donor")
  
  MSR_Donor_tidy %>%
    group_by()
  
  MSR_Acceptor_tidy <- MSR_Acceptor %>%
    dplyr::rename("Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  "Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  "CDTS 5'ss" = "CDTS_5ss",
                  "CDTS 3'ss" = "CDTS_3ss",
                  "Conservation 5'ss" = "mean_phastCons20way_5ss",
                  "Conservation 3'ss" = "mean_phastCons20way_3ss", 
                  "Protein coding" = "protein_coding") %>%
    gather(tissue, feature, -type, -tissue) %>%
    mutate(type = "MSR_Acceptor")
  
  MSR_tidy <- rbind(MSR_Donor_tidy, MSR_Acceptor_tidy) %>%
    mutate(type = factor(type, 
                         levels = c("MSR_Donor", "MSR_Acceptor")))
  
  
  # MSR_tidy <- MSR_tidy %>%
  #   filter(!(tissue %in% c("CDTS_5ss",
  #                          "CDTS_3ss",
  #                          "mean_phastCons20way_5ss",
  #                          "mean_phastCons20way_3ss")))
  
  plotTissuesLM <- ggplot(MSR_tidy, aes(tissue, feature)) + 
    geom_boxplot() +
    coord_flip() +
    facet_grid(vars(type)) +
    #ggtitle(graph_title) +
    ylab("Distribution of the significant beta values (p.adjust<0.05)") +
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
  
  
  # ggpubr::ggarrange(ggpubr::ggarrange(plotMSR,           
  #                                     plotLM,
  #                                     ncol = 2,
  #                                     labels = c("a", "b")),
  #                   plotTissuesLM, 
  #                   common.legend = F,
  #                   labels = c("", "c"),
  #                   #ncol = 2, 
  #                   nrow = 2,
  #                   #widths = c(1,1,2,2),
  #                   align = "v")
  
  
  plotTissuesLM
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/lm_alltissues")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 143, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 143, units = "mm", dpi = 300)
}



#################################
## SOMATIC MUTATIONS - SKIN
#################################

compare_skin_tissues_somatic_mutations <- function() {
  
  project_id <- "SKIN"
  tables <- c("Skin - Sun Exposed (Lower leg)",
              "Skin - Not Sun Exposed (Suprapubic)")
  
  
  if ( !file.exists("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/SKIN/v105/age_subsampled_project/results/somatic_mutations_comparison.rds") ) {
    
    skin_introns <- map_df(tables, function(cluster_id) {
      
      # cluster_id <- tables[1]
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals, gene_id FROM '", 
                      cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT DISTINCT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals, gene_id FROM '", 
                      cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_length, ref_ss5score, ref_ss3score, 
                  ref_cons5score, ref_cons3score, ref_CDTS5score, ref_CDTS3score, 
                  protein_coding, lncRNA FROM 'intron' WHERE ref_junID IN (",
                      paste(introns$ref_junID, collapse = ","),")")
      introns <- introns %>%
        left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = "ref_junID") %>% 
        as_tibble() 
      
      
      query <- paste0("SELECT id, n_transcripts, gene_width, gene_name 
                      FROM 'gene' WHERE id IN (",
                      paste(introns$gene_id, collapse = ","),")")
      introns <- introns %>%
        left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = c("gene_id" = "id")) %>% 
        as_tibble() 
      
      introns %>%
        mutate(tissue = cluster_id)
      
    } )
    
    
    query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced'")
    introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
    
    query <- paste0("SELECT ref_junID, ref_length, ref_ss5score, ref_ss3score, 
                  ref_cons5score, ref_cons3score, ref_CDTS5score, ref_CDTS3score, 
                  protein_coding
                  FROM 'intron' 
                  WHERE ref_junID IN (",
                    paste(introns$ref_junID, collapse = ","),")")
    introns <- introns %>%
      left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = "ref_junID") %>% 
      as_tibble() 
    
    
    query <- paste0("SELECT id, n_transcripts, gene_width, gene_name 
                  FROM 'gene' WHERE id IN (",
                    paste(introns$gene_id, collapse = ","),")")
    introns <- introns %>%
      left_join(y = dbGetQuery(con, query) %>% as_tibble(),
                by = c("gene_id" = "id")) %>% 
      as_tibble() 
    
    ## TIDY DATAFRAME
    df_skin_tidy <- skin_introns %>%
      mutate(tissue = tissue %>% as.factor()) %>%
      group_by(tissue) %>%
      distinct(ref_junID, .keep_all = T) %>%
      ungroup() %>%
      group_by(tissue, ref_junID) %>%
      mutate(mean_coverage = ref_sum_counts/ref_n_individuals) %>%
      mutate(mean_coverage = mean_coverage %>% log10()) %>%
      ungroup()
    
    saveRDS(object = df_skin_tidy,
            file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/SKIN/v105/age_subsampled_project/results/skin_introns.rds")
    
    ###########################
    ## SUBSAMPLING
    ###########################
    
    ## SUBSAMPLE BY MEAN COVERAGE
    ## Subsampling introns to control by similarity in mean read coverage
    m.out <- MatchIt::matchit(tissue ~ mean_coverage, 
                              data = df_skin_tidy, 
                              distance = df_skin_tidy$mean_coverage,
                              method = "nearest", 
                              caliper = c(mean_coverage = 0.005), 
                              std.caliper = FALSE)
    subsample <- MatchIt::match.data(m.out)
    subsample %>% distinct(ref_junID, .keep_all = T) %>% dplyr::count(tissue)
    
    saveRDS(object = subsample,
            file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/SKIN/v105/age_subsampled_project/results/somatic_mutations_comparison.rds")
  } else {
    df_skin_tidy <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/SKIN/v105/age_subsampled_project/results/skin_introns.rds")
    subsample <- readRDS("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/SKIN/v105/age_subsampled_project/results/somatic_mutations_comparison.rds")
    
    subsample %>%
      dplyr::select(-distance, -weights, -subclass) %>%
      inner_join(y = skin_introns %>%
                   dplyr::select(ref_junID,
                                 gene_id, ref_length, ref_ss5score, ref_ss3score,
                                 ref_cons5score, ref_cons3score, ref_CDTS5score,
                                 ref_CDTS3score,
                                 n_transcripts, gene_width),
                 by = "ref_junID")
  }
  
  ##########################
  ## PLOT - COVERAGE
  ##########################
  
  ## VISUALISE MEAN READ COVERAGE BEFORE SUBSAMPLING
  plot_coverage_bf <- ggplot(data = df_skin_tidy %>% distinct(ref_junID,.keep_all = T)) +
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
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/skin_somatic_mutations_coverage")
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
                              levels = c("Skin - Sun Exposed (Lower leg)",
                                          "Skin - Not Sun Exposed (Suprapubic)"))
  
  df_coverage %>%
    dplyr::count(tissue)
  
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
                         filter(tissue == "Skin - Sun Exposed (Lower leg)")) + 
    geom_bar(aes(x = MSR,  fill = MSR_type), 
             stat = "bin", 
             bins = 30,
             position = "dodge")  +
    ggtitle("Skin - Sun Exposed (Lower leg)") +
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
                          filter(tissue == "Skin - Not Sun Exposed (Suprapubic)")) + 
    #geom_density(aes(x = MSR, fill = MSR_type), alpha = 0.8) +
    geom_bar(aes(x = MSR,  fill = MSR_type), 
             stat = "bin", 
             bins = 30,
             position = "dodge") +
    ggtitle("Skin - Not Sun Exposed (Suprapubic)") +
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
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/skin_somatic_mutations_MSR")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  ####################################
  ## STATS
  ####################################

  df_coverage %>% filter(tissue == "Skin - Sun Exposed (Lower leg)", MSR_type == "MSR_D") %>% pull(MSR) %>% summary
  df_coverage %>% filter(tissue == "Skin - Not Sun Exposed (Suprapubic)", MSR_type == "MSR_D") %>% pull(MSR) %>% summary
  
  ## 1. Splicing noise is more frequent in introns from non-protein coding transcripts 
  ## than in introns from protein-coding transcripts even after correcting by mean read coverage
  
  ## splicing noise at the donor is more frequent in introns from NPC vs "Skin - Sun Exposed (Lower leg)" 
  wilcox.test(x = df_coverage %>% filter(tissue == "Skin - Sun Exposed (Lower leg)", MSR_type == "MSR_D") %>% pull(MSR),
              y = df_coverage %>% filter(tissue == "Skin - Not Sun Exposed (Suprapubic)", MSR_type == "MSR_D") %>% pull(MSR),
              #alternative = "greater",
              paired = T,
              correct = T)
  rstatix::wilcox_effsize(data = df_coverage %>% 
                            filter(MSR_type == "MSR_D")%>%
                            mutate(tissue = tissue %>% as.factor()),
                          formula = MSR ~ tissue,
                          paired = T)
  
  
  df_coverage %>% filter(tissue == "Skin - Sun Exposed (Lower leg)", MSR_type == "MSR_A") %>% pull(MSR) %>% summary
  df_coverage %>% filter(tissue == "Skin - Not Sun Exposed (Suprapubic)", MSR_type == "MSR_A") %>% pull(MSR) %>% summary
  
  ## splicing noise at the acceptor is more frequent in introns from NPC vs PC 
  wilcox.test(x = df_coverage %>% filter(tissue == "Skin - Sun Exposed (Lower leg)", MSR_type == "MSR_A") %>% pull(MSR),
              y = df_coverage %>% filter(tissue == "Skin - Not Sun Exposed (Suprapubic)", MSR_type == "MSR_A") %>% pull(MSR),
              #alternative = "greater",
              paired = T,
              correct = T)
  rstatix::wilcox_effsize(data = df_coverage %>% 
                            filter(MSR_type == "MSR_A")%>%
                            mutate(tissue = tissue %>% as.factor()),
                          formula = MSR ~ tissue,
                          paired = T)
  
  
  
  ####################################
  ## LM
  ####################################
  
  df_coverage %>%
    inner_join(y = skin_introns %>%
                 dplyr::select(ref_junID,
                               gene_id, ref_length, ref_ss5score, ref_ss3score,
                               ref_cons5score, ref_cons3score, ref_CDTS5score,
                               ref_CDTS3score,
                               protein_coding, n_transcripts, gene_width, gene_name),
               by = "ref_junID")
  
  
  
  #########################
  ## TIDY DATA
  #########################
  
  idb <- introns %>%
    as.data.frame() %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    dplyr::rename(intron_length = ref_length,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score,
                  gene_length = gene_width,
                  gene_tpm = gene_tpm,
                  gene_num_transcripts = n_transcripts,
                  CDTS_5ss = ref_CDTS5score,
                  CDTS_3ss = ref_CDTS3score,
                  protein_coding = protein_coding,
                  mean_phastCons20way_5ss = ref_cons5score,
                  mean_phastCons20way_3ss = ref_cons3score)
  
  
  #########################
  ## LINEAR MODELS
  #########################
  
  
  fit_donor <- lm(MSR_D ~ 
                    intron_length +
                    intron_5ss_score + 
                    intron_3ss_score +
                    gene_length +
                    gene_tpm +
                    gene_num_transcripts +
                    protein_coding +
                    CDTS_5ss +
                    CDTS_3ss +
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss, 
                  data = idb)
  
  
  # fit_donor$coefficients[names(fit_donor$coefficients)[summary(fit_donor)$coefficients[,4] < 0.05]]  
  
  fit_acceptor <- lm(MSR_A ~ 
                       intron_length + 
                       intron_5ss_score + 
                       intron_3ss_score +
                       gene_length +
                       gene_tpm +
                       gene_num_transcripts +
                       protein_coding +
                       CDTS_5ss +
                       CDTS_3ss +
                       mean_phastCons20way_5ss +
                       mean_phastCons20way_3ss, 
                     data = idb)
  fit_acceptor %>% summary()
  
  model_names <- c("MSR_Donor", "MSR_Acceptor")
  
  
  coef_names <- c("Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  #"Inter. 5'ss & 3'ss MES score" = "intron_5ss_score:intron_3ss_score", 
                  #"Intron Type u2" = "u2_intronTRUE",
                  #"Intron ClinVar mutation" = "clinvarTRUE",
                  "Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  "CDTS 5'ss" = "CDTS_5ss",
                  "CDTS 3'ss" = "CDTS_3ss",
                  "Cons 5'ss" = "mean_phastCons20way_5ss",
                  "Cons 3'ss" = "mean_phastCons20way_3ss",
                  "Protein coding" = "protein_coding")
  
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
                               colors = c("#35B779FF","#440154FF"),
                               model.names = model_names,
                               plot.distributions = T,
                               rescale.distributions = T,
                               
                               point.shape = T
  ) + 
    theme_minimal() + 
    theme(axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = "black", size = "12"),
          axis.text.y = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2)) +  
    ylab("Covariates") +
    geom_hline(yintercept = seq(from = 0,
                                to = length((fit_donor$coefficients %>% names)[-1]) + .5,
                                by = 1)) +
    guides(colour = guide_legend(#title = NULL,
      ncol = 2, 
      nrow = 1))     +
    scale_x_continuous(breaks = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025),
                       labels = c(-0.01,-0.0075,-0.005,-0.0025, 0, 0.0025))
  
  plotLM 
  plotLM +
    annotate("text",
             y = c(1.65, 2.65, 3.65, 4.65, 5.65, 6.65, 7.65, 8.65, 9.65, 10.65, 11.65),
             x = -0.0025, 
             label = c( paste("pval",(summary(fit_donor)$coefficients[,4] %>% 
                                        as_tibble() %>%
                                        mutate(value = format(value, format = "e", digits = 2)))$value[-1] %>% rev(), sep = ":")),
             family = "", 
             fontface = 1, 
             size=3,
             colour = c("#35B779FF") )  +
    annotate("text",
             y = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 7.25, 8.25, 9.25, 10.25, 11.25),
             x = -0.0025, 
             label = c( paste("pval",(summary(fit_acceptor)$coefficients[,4] %>% 
                                        as_tibble() %>%
                                        mutate(value = format(value, format = "e", digits = 2)))$value[-1] %>% rev(), sep = ":")),
             family = "", 
             fontface = 1, 
             size=3,
             colour = c("#440154FF") )
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/paper/",
                      main_project, "/figures/lm_FTCX")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
}

#################################
## GO ENRICHMENT
#################################

## TODO: not only obtain common introns across age groups, but also obtain common introns across age groups with similar coverage.
## Then analyse changes in splicing noise across age


plot_GO_Enrichment <- function() {
  
  
  project_id <- "BRAIN"
  
  ## Read data.frame containing only common introns across age groups in BRAIN
  folder_results <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", 
                           gtf_version, "/age_subsampled/")
  
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
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
}

supplementary_age_stratification <- function() {
  
  # project_id <- "MUSCLE"
  # project_id <- "BLOOD"
  
  source(paste0("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline4-2_age_stratification.R"))
  
  ref <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  
  age_samples_clusters_tidy <- age_stratification_init_data(project_id)
  age_supergroups <- age_samples_clusters_tidy$age_group %>% unique()
  
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", project_id,
                        "/results/pipeline3/missplicing-ratio/age/")
  
  ## Load the IDBs
  df_age_groups_intron <- map_df(age_supergroups, function(age_group) {
    
    # age_group <- age_supergroups[1]
    readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_introns.rds")) %>%
      #mutate(ref_junID = ref_junID %>% as.integer()) %>%
      mutate(sample_type = age_group) %>%
      return()
    
  })
  
  df_age_groups_novel <- map_df(age_supergroups, function(age_group) {
    
    #print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
    readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_novel.rds")) %>%
      #mutate(ref_junID = ref_junID %>% as.integer()) %>%
      #mutate(novel_junID = novel_junID %>% as.integer()) %>%
      mutate(sample_type = age_group) %>%
      return()
    
  }) 
  
  ## Get the junctions that are common across age supergroups
  
  common_introns_misspliced <- df_age_groups_novel %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_supergroups %>% length()) %>%
    pull(ref_junID)
  
  common_introns <- df_age_groups_intron %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_supergroups %>% length()) %>%
    pull(ref_junID)
  
  # setdiff(common_introns_misspliced,common_introns)
  # setdiff(common_introns,common_introns_misspliced)
  # common_introns_misspliced %>% length()
  # common_introns %>% length()
  # intersect(common_introns_misspliced,common_introns) %>% length()
  
  ## Filter by the common junctions
  df_age_groups_intron_tidy <- merge(x = df_age_groups_intron %>% data.table::as.data.table(),
                                     y = data.table::data.table(ref_junID = common_introns_misspliced),
                                     by = "ref_junID",
                                     all.y = T)
  
  saveRDS(object = df_age_groups_intron_tidy,
          file = paste0(folder_root, "/common_introns_age.rds"))
  
  df_age_groups_novel_tidy <- merge(x = df_age_groups_novel %>% data.table::as.data.table(),
                                    y = data.table::data.table(ref_junID = common_introns_misspliced),
                                    by = "ref_junID",
                                    all.y = T)
  
  
  
  bg_genes <- df_age_groups_novel_tidy %>% distinct(gene_name) %>% pull() %>% unlist %>% unique()
  
  
  
  ##########################################################
  ## DISTANCES
  ##########################################################
  
  distance_plot <- age_stratification_plot_distances(df = df_age_groups_novel_tidy, 
                                                     age_levels = c("60-79","40-59","20-39" ), 
                                                     distance_limit = 30)
  
  muscle_distance_plot <- distance_plot
  blood_distance_plot <- distance_plot
  ##########################################################
  ## MSR PLOTS
  ##########################################################
  
  MSR_plot <- age_stratification_plot_MSR(df = df_age_groups_intron_tidy,
                                          age_levels = c("20-39","40-59","60-79" )) 
  
  muscle_MSR_plot <- MSR_plot
  blood_MSR_plot <- MSR_plot
  
  
  ## ARRANGE -----------------------------------------------
  
  ggpubr::ggarrange(muscle_distance_plot,
                    blood_distance_plot,
                    muscle_MSR_plot,
                    blood_MSR_plot,
                    nrow = 2,
                    ncol = 2,
                    labels = c("a", "b", "c", "d"))
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/supplementary_age_stratification_panel.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ##########################################################
  ## GO ENRICHMENT
  ##########################################################
  
  ## MSR_D -------------------------------------------------
  
  df_MSRD <- df_age_groups_intron_tidy %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_D = ref_missplicing_ratio_tissue_ND,
                  gene_name) %>%
    mutate(MSR_D = MSR_D %>% round(digits = 4)) %>%
    spread(sample_type, MSR_D)
  
  
  genes_MSRD_increasing <- df_MSRD %>%
    filter(`20-39` < `40-59`,
           `40-59` < `60-79`) %>%
    mutate_if(is.list, simplify_all) %>%   
    mutate(`20-39` = replace(`20-39`, `20-39` == 0, 0.0000001)) %>%
    mutate(`60-79` = replace(`60-79`, `60-79` == 0, 0.0000001)) %>%
    mutate(fold_increase = gtools::foldchange(`20-39`, `60-79`)) %>%
    mutate(MSR_type = "MSR_Donor") %>%
    dplyr::select(gene_name, fold_increase, MSR_type) %>%
    arrange(desc(fold_increase))
  
  
  ego_D <- clusterProfiler::enrichGO(
    gene          = genes_MSRD_increasing$gene_name %>% unlist() %>% unique(), #match(x = (genes_MSRD_increasing$gene_name %>% unlist() %>% unique()), bg_genes),
    universe      = bg_genes,
    keyType       = "SYMBOL",
    OrgDb         = "org.Hs.eg.db", ##Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)
  # head(ego)
  
  
  MSR_D_GO_plot <- barplot(ego_D, showCategory = 20)
  MSR_D_GO_plot <- MSR_D_GO_plot  +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
    xlab("Gene count") +
    #coord_flip() +
    facet_grid(ONTOLOGY~., scale = "free") +
    #theme_minimal() + 
    theme(axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = "black", size = "9"),
          #axis.text.y = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(colour = "black", size = "9"),
          legend.title = element_text(colour = "black", size = "9"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.text.y = element_text(colour = "black", size = "9",#angle =10, 
                                     vjust = 0.3,
                                     hjust = 1))
  MSR_D_GO_plot
  
  edox <- clusterProfiler::setReadable(ego_D, 'org.Hs.eg.db', 'SYMBOL')
  edox2 <- enrichplot::pairwise_termsim(edox)
  
  
  enrichplot::treeplot(edox2)+
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
    #theme_minimal() + 
    theme(#axis.line = element_line(colour = "black"), 
      #axis.text.x = element_text(colour = "black", size = "9"),
      #axis.text.y = element_text(colour = "black", size = "9"),
      axis.title = element_text(colour = "black", size = "9"),
      legend.text = element_text(colour = "black", size = "9"),
      legend.title = element_text(colour = "black", size = "9"),
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      # panel.grid.major.y = element_blank(),
      axis.ticks = element_line(colour = "black", size = 2))#,
  #axis.text.y = element_text(colour = "black", size = "9",#angle =10, 
  #                           vjust = 0.3,
  #                           hjust = 1))
  
  ## MSR_A -------------------------------------------------
  
  
  df_MSRA <- df_age_groups_intron_tidy %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_A = ref_missplicing_ratio_tissue_NA,
                  gene_name) %>%
    mutate(MSR_A = MSR_A %>% round(digits = 4)) %>%
    spread(sample_type, MSR_A)
  
  
  genes_MSRA_increasing <- df_MSRA %>%
    filter(`20-39` < `40-59`,
           `40-59` < `60-79`) %>%
    mutate_if(is.list, simplify_all) %>%
    mutate(`20-39` = replace(`20-39`, `20-39` == 0, 0.0000001)) %>%
    mutate(`60-79` = replace(`60-79`, `60-79` == 0, 0.0000001)) %>%
    mutate(fold_increase = gtools::foldchange(`20-39`, `60-79`)) %>%
    mutate(MSR_type = "MSR_Acceptor") %>%
    dplyr::select(gene_name, fold_increase, MSR_type) %>%
    arrange(desc(fold_increase))
  
  
  
  ego_A <- clusterProfiler::enrichGO(
    gene          = genes_MSRA_increasing$gene_name %>% unlist() %>% unique(), #match(x = (genes_MSRD_increasing$gene_name %>% unlist() %>% unique()), bg_genes),
    universe      = bg_genes,
    keyType       = "SYMBOL",
    OrgDb         = "org.Hs.eg.db", ##Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)
  #head(ego)
  MSR_A_GO_plot <- barplot(ego_A, showCategory = 10)  +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
    xlab("Gene count") +
    #coord_flip() +
    facet_grid(ONTOLOGY~., scale = "free") +
    #theme_minimal() + 
    theme(axis.line = element_line(colour = "black"), 
          #axis.text.x = element_text(colour = "black", size = "9"),
          #axis.text.y = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(colour = "black", size = "9"),
          legend.title = element_text(colour = "black", size = "9"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.text.y = element_text(colour = "black", 
                                     size = "9",#angle =10, 
                                     vjust = 0.3,
                                     hjust = 1),
          axis.text.x = element_text(colour = "black", 
                                     size = "9",
                                     angle = 40, 
                                     vjust = 1,
                                     hjust = 1)) +
    guides(fill = guide_legend(title = "pval")) 
  
  MSR_A_GO_plot
  
  
  
  # # ## MSR_D AND MSR_A TOGETHER --------------------------------
  # # 
  # # ggpubr::ggarrange(MSR_D_GO_plot, 
  # #                   MSR_A_GO_plot,
  # #                   labels = c("a", "b"),
  # #                   ncol = 2, 
  # #                   nrow = 1)
  # # 
  # # 
  # # 
  # # MSR_D_GO_plot
  # # MSR_A_GO_plot
  # # 
  # # #########
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # genes_MSRD_dcreasing <- df_MSRD %>%
  #   filter(`20-39` > `40-59`,
  #          `40-59` > `60-79`) %>%
  #   mutate_if(is.list, simplify_all) %>%  
  #   mutate(`20-39` = replace(`20-39`, `20-39` == 0, 0.0000001)) %>%
  #   mutate(`60-79` = replace(`60-79`, `60-79` == 0, 0.0000001)) %>%
  #   mutate(fold_decrease = gtools::foldchange(`20-39`, `60-79`)) %>%
  #   arrange(desc(fold_decrease))
  # 
  # intersect(setdiff(genes_MSRD_increasing$gene_name, 
  #                   genes_MSRD_dcreasing$gene_name), bg_genes)
  # 
  # 
  # 
  # 
  # gene_enrichment <- gprofiler2::gost(query = setdiff(genes_MSRD_increasing$gene_name, 
  #                                                     genes_MSRD_dcreasing$gene_name),
  #                                     custom_bg = bg_genes,
  #                                     organism = "hsapiens",
  #                                     ordered_query = F,
  #                                     correction_method = "bonferroni",
  #                                     significant = T)
  # 
  # GO_enrichment <- gene_enrichment$result %>% filter(str_detect(source, pattern = "GO")) %>%
  #   mutate(go_type = str_sub(source, start = 4, end = 5))
  # 
  # GO_enrich_reduc <- rutils::go_reduce(pathway_df = data.frame(go_type = GO_enrichment$go_type,
  #                                                              go_id = GO_enrichment$term_id,
  #                                                              go_term = GO_enrichment$term_name),
  #                                      orgdb = "org.Hs.eg.db",
  #                                      threshold = 0.7,
  #                                      scores = NULL,
  #                                      measure = "Wang"
  # )
  # 
  # df_result <- merge(x = gene_enrichment$result %>% filter(source %in% c("GO:BP", "GO:MF", "GO:CC", "KEGG")),
  #                    y = GO_enrich_reduc,
  #                    by.x = "term_id",
  #                    by.y = "go_id",
  #                    all.x = T)
  # 
  # ## JOIN with KEGG results and add pvalues
  # 
  # rbind(df_result %>%
  #         as_tibble() %>%
  #         dplyr::select(parent_term, p_value, query_size, intersection_size, source) %>%
  #         group_by(parent_term, source) %>%
  #         mutate(p_value = p_value %>% max) %>%
  #         ungroup() %>%
  #         distinct(parent_term, .keep_all = T),
  #       df_result %>% 
  #         filter(source == "KEGG")%>%
  #         dplyr::select(p_value, query_size, intersection_size, source, term_name) %>%
  #         dplyr::rename(parent_term = term_name)) %>%
  #   arrange(p_value) %>%
  #   mutate(p_value = p_value %>% formatC(format = "e", digits = 2)) %>%
  #   #filter(source == "KEGG") %>%
  #   DT::datatable()
  # 
  
  #############################################################################
  ## RBPs
  #############################################################################
  
  projects_id <- c("BRAIN", "MUSCLE", "BLOOD")
  
  df_tissues <- map_df(projects_id, function(project_id) {
    
    # project_id <- projects_id[1]
    
    df_analysis_lm_uncorrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", project_id,
                                                         "/results/pipeline3/rbp/tpm_lm_all.csv")) %>%
      mutate(tissue = project_id) %>%
      drop_na() %>%
      distinct(name, .keep_all = T)
    
    return(df_analysis_lm_uncorrected)
    
  })
  
  RBP_tissue_comparison <- ggplot(df_tissues) + 
    geom_density(aes(x=Estimate, fill=tissue), alpha = 0.7) +
    facet_grid(vars(type))+ 
    #ggtitle("Age effect in the predicction of the uncorrected TPM values of\n208 spliceosomal RBPs across samples of each body region") +
    xlab("Age Effect (i.e. estimate)") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          strip.text = element_text(colour = "black", size = "9"),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 3, 
                               nrow = 1))
  
  
  ## TOGETHER ----------------------------------------------------
  
  ggpubr::ggarrange(distance_plot,
                    MSR_plot,
                    MSR_A_GO_plot,
                    RBP_tissue_comparison,
                    nrow = 2,
                    ncol = 2,
                    labels = c("c", "d", "e", "f"))
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/age_stratification_panel.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
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

compare_skin_tissues_somatic_mutations()

  