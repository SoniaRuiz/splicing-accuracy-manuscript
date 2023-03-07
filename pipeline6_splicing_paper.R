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
# all_projects %>% length() %>% print()


get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}


custom_ggtheme <-  theme(text = element_text(size = 7, family="Arial", colour = "black"),
                         legend.text = element_text(size = "7", family="Arial", colour = "black"),
                         axis.ticks = element_line(colour = "black", linewidth = 2),
                         axis.text = element_text(size = 7, family="Arial", colour = "black"),
                         axis.line = element_line(colour = "black"),
                         axis.title = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.y = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.x = element_text(size = 7, family="Arial", colour = "black", 
                                                    hjust = 0.5, vjust = 0.5),
                         strip.text = element_text(size = 7, family="Arial", colour = "black"),
                         legend.position = "top",
                         legend.box = "vertical")

########################################
## FUNCTIONS - produce figures for the
## paper
########################################


## SECTION 1 ---------------------------------------------


## 1. Novel donor and acceptor junctions are commonly detected and exceed the number of unique annotated introns by an average of 11-fold

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
  
  ## We found that 268,988 (82.8%) annotated introns had at least a single associated novel donor or acceptor junction,
  ## with only 55,968 annotated introns appearing to be precisely spliced across all the samples and tissues studied.
  query <- paste0("SELECT * from 'intron'")
  db_introns <- dbGetQuery(con, query) 
  db_introns %>% distinct(ref_junID) %>% nrow()
  db_introns %>%
    dplyr::count(misspliced)
  
  
  ## Collectively, we detected 3,865,268 unique novel junctions, equating to 14 novel junctions per annotated intron. 
   query <- paste0("SELECT * from 'novel'")
  db_novel <- dbGetQuery(con, query) 
  db_novel %>% nrow() 
  
  db_novel %>% 
    dplyr::count(novel_type)
  
  (db_novel %>% nrow()) / (db_introns %>% filter(misspliced==1) %>% nrow())
  
  
  
  ## Novel junctions exceed in X fold to annotated introns
  (db_novel %>% distinct(novel_junID) %>% nrow()) / (db_introns %>% distinct(ref_junID) %>% nrow())
  
  ## Percentage of mis-spliced introns
  (((db_introns %>%
      dplyr::count(misspliced) %>%
      filter(misspliced == 1) %>%
      pull(n)) * 100 )) /
  ((db_introns %>% distinct(ref_junID) %>% nrow()) %>%
    round(digits = 1))
  
  
  ## Collectively, we detected X novel junctions
  db_novel %>% distinct(novel_junID) %>% nrow()
  
  
  ## equating to 14 novel junctions per an annotated junction.
  ((db_novel %>% distinct(novel_junID) %>% nrow()) / (db_introns %>%
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
  
  db_metadata_tidy %>%
    inner_join(y = df_novel_jxn_count,
               by = "cluster") %>% 
    mutate(n_novel_sample = n_novel/n) %>% 
    arrange(desc(n_novel_sample))
}

get_junc_length <- function() {
  
  
  tables <- dbListTables(con)
  tables
  
  ## Get all annotated introns
  query <- paste0("SELECT DISTINCT ref_junID, ref_length from 'intron'")
  db_introns <- dbGetQuery(con, query) %>% as_tibble()
  db_introns %>% distinct(ref_junID) %>% nrow()
  db_introns %>% head()
  
  ## Get all novel junctions
  query <- paste0("SELECT DISTINCT novel_junID, novel_length from 'novel'")
  db_novel <- dbGetQuery(con, query) %>% as_tibble()
  db_novel %>% distinct(novel_junID) %>% nrow()
  
  ## Join both datasets of junction lengths
  df_all_lengths_tidy <- rbind(db_introns %>%
                                 dplyr::select(length = ref_length) %>%
                                 mutate(type = "annotated intron"),
                               db_novel %>%
                                 dplyr::select(length = novel_length) %>%
                                 mutate(type = "novel junction"))
  df_all_lengths_tidy %>% head()
  
  
  ## Get the mode length of the annotated introns
  get_mode( df_all_lengths_tidy %>%
              filter(type == "annotated intron") %>%
              pull(length) )
  
  ## Get the mean length of the annotated introns
  df_all_lengths_tidy %>%
    filter(type == "annotated intron") %>%
    pull(length) %>% mean()
  
  
  ## Get the mode length of the novel junctions
  get_mode( df_all_lengths_tidy %>%
              filter(type == "novel junction") %>%
              pull(length) )
  
  ## Get the mean length of the novel junctions
  df_all_lengths_tidy %>%
    filter(type == "novel junction") %>%
    pull(length) %>% mean()
  
  
  
  #############################################
  ## PLOT IMPLIED INTRON LENGTHS
  #############################################
  
  df_all_lengths_tidy$type = factor( df_all_lengths_tidy$type, 
                                     levels = c( "novel junction", "annotated intron" ) )
  
  
  ggplot(data = df_all_lengths_tidy %>%
           drop_na() %>%
           filter(length <= 200)) + 
    geom_count(aes(x = length, 
                   colour = type), 
               alpha = 0.7,
               stat = "count",
               position = "identity") +
    xlim(c(0,200))+
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
    custom_ggtheme %>%
    return()
  
  ## Save plot
  file_name <- paste0(getwd(), "/results/_paper/figures/junction_length")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 150, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 150, height = 100, units = "mm", dpi = 300)
  
}




## 2. Over 98% of novel donor and acceptor junctions are likely to be generated through splicing errors

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
  ggplot(data = df_contamination_tidy) +
    geom_bar(mapping = aes(x = tissue, y = prop, fill = type), stat = "identity", position = "identity") + 
    #coord_flip() +
    ggforce::facet_zoom(ylim = c(0,1.5), zoom.size = 3) +
    ylab("% unique novel junctions") +
    xlab("Tissue") +
    scale_fill_manual(values = c( "#a6a6a6", "#1a1a1a" ),
                      breaks = c( "kept_annotation", "in_annotation" ),
                      labels = c( "novel junctions mantaining category", "novel junctions entring annotation" )) +
    #scale_fill_viridis_d(option = "A") +
    theme_light() +
    custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1).
          ax) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) %>%
    return()
  
  figures_path <- paste0(getwd(), "/results/_paper/figures/")
  dir.create(file.path(figures_path), recursive = TRUE, showWarnings = T)
  file_name <- paste0(figures_path, "/contamination_rates_all_tissues")
  
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 100, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 100, units = "mm", dpi = 300)
  
  
  
  
}

get_contamination_rates_all_tissues_stats <- function() {
  
  
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
  
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data=df_contamination, aes(x = contamination_rates, 
                                    y = in_annotation,
                                    group=1)) +
    geom_line(linetype = "dashed", col = "red") +
    geom_point() +
    xlab(NULL) +
    ylab("contamination rates (%)\nas compared to Ensembl v105") +
    theme_light() +
    custom_ggtheme +
    guides(fill = "none")
  
  
  file_name <- paste0(getwd(), "/results/_paper/figures/contamination_rates_FCTX_lineplot")
  ggplot2::ggsave(paste0(file_name, ".svg"),  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"),  width = 180, height = 90, units = "mm", dpi = 300)

}




## 3. Mis-splicing is more common at acceptor than donor splice sites

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
  
  reads_violin_plot <- get_unique_donor_acceptor_reads()
  
  # file_name <- paste0(getwd(), "/results/_paper/figures/unique_donor_acceptor_reads_violin")
  # ggplot2::ggsave(paste0(file_name, ".svg"), width = 150, height = 150, units = "mm", dpi = 300)
  # ggplot2::ggsave(paste0(file_name, ".png"), width = 90, height = 90, units = "mm", dpi = 300)
  # 
  
  
  ggpubr::ggarrange(junx_violin_plot,
                    reads_violin_plot,
                    common.legend = T,
                    labels = c("b", "c"),
                    align = "h",
                    ncol = 2,
                    nrow = 1)
  
  file_name <- paste0(getwd(), "/results/_paper/figures/panel2bc")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 90, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
  # ## Save the figure 3
  # file_name <- paste0(getwd(), "/results/_paper/figures/percent_unique_donor_acceptor_violin")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), 
  #                 width = 183, height = 130, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), 
  #                 width = 90, height = 90, units = "mm", dpi = 300)
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
    
  
  reads_violin_plot
  
  
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
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      db_never <- dbGetQuery(con, query) %>% as_tibble()
      db_never$ref_sum_counts
      
      query <- paste0("SELECT ref_junID, novel_junID, ref_sum_counts, novel_sum_counts FROM '", cluster_id, "_", project_id, "_misspliced'")
      db_misspliced <- dbGetQuery(con, query) %>% as_tibble()
      db_misspliced$ref_sum_counts
      
      novel_donor <- db_misspliced  %>%
        filter(novel_junID %in% (db_novel %>% 
                                   filter(novel_type == "novel_donor") %>%
                                   pull(novel_junID)))
      
      novel_acceptor <- db_misspliced  %>%
        filter(novel_junID %in% (db_novel %>% 
                                   filter(novel_type == "novel_acceptor") %>%
                                   pull(novel_junID)))
      
      
      return( data.frame(annotated_unique_junctions = c(db_never$ref_junID,
                                                        db_misspliced$ref_junID) %>% unique() %>% length(),
                         annotated_median_count = c(db_never$ref_sum_counts,
                                                    db_misspliced$ref_sum_counts) %>% median(),
                         annotated_max_count = c(db_never$ref_sum_counts,
                                                 db_misspliced$ref_sum_counts) %>% max(),
                         annotated_min_count = c(db_never$ref_sum_counts,
                                                 db_misspliced$ref_sum_counts) %>% min(),
                         novel_donor_unique_junctions = novel_donor$novel_junID %>% unique() %>% length(),
                         novel_donor_median_count = novel_donor$novel_sum_counts %>% median(),
                         novel_donor_max_count = novel_donor$novel_sum_counts %>% max(),
                         novel_donor_min_count = novel_donor$novel_sum_counts %>% min(),
                         novel_acceptor_unique_junctions = novel_acceptor$novel_junID %>% unique() %>% length(),
                         novel_acceptor_median_count = novel_acceptor$novel_sum_counts %>% median(),
                         novel_acceptor_max_count = novel_acceptor$novel_sum_counts %>% max(),
                         novel_acceptor_min_count = novel_acceptor$novel_sum_counts %>% min(),
                         tissue = cluster_id) )
      
    })
    
   
  })
  
  
  
  write.csv(x = df_read_count %>%
              relocate(tissue),
            file = paste0(getwd(), "/results/_paper/results/05_supplementary_data_read_count.csv"),
            row.names = F)

 
  
  ## Focusing on frontal cortex, we found that this equated to a median read count 
  ## of 2,700 for annotated junctions with novel donor or acceptor events having a 
  ## median read count of only 2 in both cases. 
  
  df_read_count %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)")
  
}




## 4. High sequence similarity between novel splice sites and their annotated pairs explains splicing errors

get_maxentscan_score <- function() {
  
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  
  file_name <- paste0(getwd(), "/results/_paper/results/df_mes.rds")
  
  if ( !file.exists(file_name) )  {
    
  
    ###########################
    ## GET THE MES
    ###########################
    
    
    ## Get the MES scores corresponding to the mis-spliced introns
    query <- paste0("SELECT ref_junID, ref_ss5score, ref_ss3score FROM 'intron' WHERE  misspliced = 1 ")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    
    ## Get the novel junctions
    query <- paste0("SELECT ref_junID, novel_junID, novel_type, novel_ss5score, novel_ss3score FROM 'novel'")
    novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
    
    ###########################
    ## MERGE THE DATA AND RETURN
    ###########################
    
    df_mes <- novel_junctions %>% 
      left_join(y = introns,
                by = "ref_junID")  %>%
      mutate(diff_ss5score = ref_ss5score - novel_ss5score,
             diff_ss3score = ref_ss3score - novel_ss3score)
    
    saveRDS(object = df_mes %>% as_tibble(), file = file_name)
    
  } else {
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
    dplyr::select(intron = ref_ss5score, novel_donor = novel_ss5score) %>%
    gather(key = "junction_type", value = "ss5score")
  
  ss5plot <- ggplot(df_5ss, aes(ss5score, fill = junction_type)) +
    geom_density(alpha = 0.8) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    xlab("MES Donor (5'ss) score") +
    scale_fill_manual(values = c("#A41F9AFF", "#EF7F4FFF"), 
                      breaks=c("intron", "novel_donor"),
                      labels=c("Annotated Intron  ", "Novel Donor")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = element_blank(),
                               ncol = 2, nrow = 1))
  
  ss5plot
  # file_name <- paste0(getwd(), "/results/_paper/figures/MES_ss5plot")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 90, height = 90, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 90, height = 90, units = "mm", dpi = 300)
  # 
  
  ## ss3score -----------------------------------------------------------
  
  df_3ss <- df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron = ref_ss3score, novel_acceptor = novel_ss3score) %>%
    gather(key = "junction_type", value = "ss3score")
  
  ss3plot <- ggplot(df_3ss, aes(ss3score, fill = junction_type)) +
    geom_density(alpha = 0.8) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    scale_fill_manual(values = c("#A41F9AFF", "#0D0887FF"), 
                      breaks = c("intron", "novel_acceptor"),
                      labels = c("Annotated Intron  ", "Novel Acceptor")) +
    xlab("MES Acceptor (3'ss) score") +
    custom_ggtheme +
    guides(fill = guide_legend(title = element_blank(), ncol = 2, nrow = 1))
  
  
  ss3plot
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
  file_name <- paste0(getwd(), "/results/_paper/figures/supplementary_fig3")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 90, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
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
  
  deltaplot5ss
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
  
  deltaplot3ss
  
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
  
  file_name <- paste0(getwd(), "/results/_paper/figures/panel3ab")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 175, height = 55, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 175, height = 55, units = "mm", dpi = 300)
  
  
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
    table() %>% 
      as.data.frame() %>%
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




## 5. Novel junctions associated with protein-coding transcripts are predicted to be deleterious in 63.5% of cases

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
  
  # file_name <- paste0(getwd(), "/results/_paper/figures/distances_FCTX")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ###############################
  ## PROTEIN CODING
  ###############################
  
  query <- paste0("SELECT ref_junID, protein_coding, lncRNA FROM 'intron'")
  df_master_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  df_novel_tidy <- db_misspliced_introns %>%
    inner_join(df_master_introns,
               by = "ref_junID") %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC"),
           novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "),
           novel_type = str_to_title(novel_type))
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("Novel Donor", "Novel Acceptor"))
  
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
    ggtitle("Protein-coding transcripts") +
    xlab("Distance (in bp)") +
    ylab("Unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme +
    theme(legend.position = "none") +
    
    #########
    geom_segment(data = data.frame(x = 5.8, xend = 9.1, linewidth = 1,
                                   y = 980, yend = 980,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 5.9, xend = 5.9, linewidth = 1, 
                        y = 960, yend = 980,
                        colour = "#333333",
                        novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 9.0, xend = 9.0, linewidth = 1,
                                   y = 960, yend = 980,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend) ) +
    #########
  
    geom_segment(data = data.frame(x = 8.8, xend = 12.1, linewidth =1,
                                   y = 780, yend = 780,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 8.9, xend = 8.9, linewidth =1, 
                                   y = 760, yend = 780,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 12.0, xend = 12.0, linewidth =1,
                                   y = 760, yend = 780,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend) ) +
    #########
    geom_segment(data = data.frame(x = 11.8, xend = 15.1, linewidth =1,
                                 y = 550, yend = 550,
                                 colour = "#333333",
                                 novel_type="Novel Acceptor"),
               aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 11.9, xend = 11.9, linewidth =1, 
                                   y = 530, yend = 550,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 15.0, xend = 15.0, linewidth =1,
                                   y = 530, yend = 550,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend) ) +
    

    #########
    geom_segment(data = data.frame(x = 14.8, xend = 18.1, linewidth =1,
                                 y = 440, yend = 440,
                                 colour = "#333333",
                                 novel_type="Novel Acceptor"),
               aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 14.9, xend = 14.9, linewidth =1, 
                                   y = 420, yend = 440,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 18.0, xend = 18.0, linewidth =1,
                                   y = 420, yend = 440,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend) ) +
    
    #########
    geom_segment(data = data.frame(x = 17.8, xend = 21.1, linewidth =1,
                                 y = 350, yend = 350,
                                 colour = "#333333",
                                 novel_type="Novel Acceptor"),
               aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 17.9, xend = 17.9, linewidth =1, 
                                   y = 330, yend = 350,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 21.0, xend = 21.0, linewidth =1,
                                   y = 330, yend = 350,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend) ) +
    
    
    #########
  geom_segment(data = data.frame(x = 20.8, xend = 24.1, linewidth =1,
                                 y = 290, yend = 290,
                                 colour = "#333333",
                                 novel_type="Novel Acceptor"),
               aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 20.9, xend = 20.9, linewidth =1, 
                                   y = 270, yend = 290,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend)) +
    geom_segment(data = data.frame(x = 24.0, xend = 24.0, linewidth =1,
                                   y = 270, yend = 290,
                                   colour = "#333333",
                                   novel_type="Novel Acceptor"),
                 aes(x=x,y=y,yend=yend,xend=xend) ) +
    
    ggplot2::facet_grid(vars(factor(novel_type, levels=c('Novel Donor','Novel Acceptor'))))
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp * -1), xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 80),  size = 3, label = "intron") +
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
    ylab("") +
    theme_light() +
    ggtitle("Non-coding transcripts")+
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      breaks = c("Novel Donor", "Novel Acceptor")) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1 )) +
    custom_ggtheme 
  
  plot_legend <- ggpubr::get_legend(plot_NPC)
  
  plot_NPC <- plot_NPC+ ggpubr::rremove("legend")
  
  plot_NPC <- plot_NPC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  plot_NPC
  
  
  # file_name <- paste0(getwd(), "/results/_paper/figures/FCTX_distancesNPC_lncRNA")
  # ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  # ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 100, height = 100, units = "mm", dpi = 300)
  
  
  #######################################
  ## COMBINE ALL 2 PLOTS
  #######################################
  
  plot <- ggpubr::ggarrange(plot_PC,
                            plot_NPC + ggpubr::rremove("ylab")+ ggpubr::rremove("legend"),
                            common.legend = T,
                            legend = "top",
                            
                            labels = c("c", "d"),
                            align = "v",
                            ncol = 2, 
                            nrow = 1,
                            legend.grob = plot_legend)
  
  
  
  
  plot
  
  
  file_name <- paste0(getwd(), "/results/_paper/figures/panel3cd")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 90, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
  
  #########################################
  ## STATS
  #########################################
  
  df_novel_tidy %>% 
    distinct(ref_junID, .keep_all = T) %>%
    dplyr::count(type_PC)
  
}

get_distances_all_tissues <- function() {
  
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  limit_bp <- 30
  
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  df_distances_all_tissues <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    ## GET THE CLUSTERS FOR THE CURRENT TISSUE
    
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    print(all_clusters)
    
    map_df(all_clusters, function(cluster_id) {
      
      # cluster_id <- all_clusters[1]
      
      query <- paste0("SELECT tissue.novel_junID, tissue.ref_junID, novel.novel_type, novel.distance 
                      FROM '", cluster_id, "_", project_id, "_misspliced' AS tissue
                      INNER JOIN 'novel' ON novel.novel_junID = tissue.novel_junID
                      WHERE distance >= -30 AND distance <= 30")
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
  
  df_distances_all_tissues %>%
    filter(tissue == "Brain - Frontal Cortex (BA9)")
  
  write.csv(x = df_distances_all_tissues,
            file = paste0(getwd(), "/results/_paper/results/05_supplementary_data_distances_all_tissues.csv"),
            row.names = F)
}

get_modulo <- function() {
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  
  file_path <- paste0(getwd(), "/results/_paper/results/df_modulo_basic_tissues_100bpfilter.rds")
  
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
        query <- paste0("SELECT intron.ref_junID, intron.protein_coding, intron.MANE, transcript.MANE AS tMANE
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
  
  library(ggridges)
  
  df_modulo_tissues <- df_modulo_tissues %>%
    mutate(freq = freq * 100)
  
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
  
  
  file_name <- paste0(getwd(), "/results/_paper/figures/panel3e")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 55, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 55, units = "mm", dpi = 300)
  
  
  
  
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
  

  file_name <- paste0(getwd(), "/results/_paper/figures/panel4a")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 75, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 75, units = "mm", dpi = 300)
  
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
    custom_ggtheme +
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
    custom_ggtheme +
    xlab("log10 mean read coverage")
  
  
  ggpubr::ggarrange(plot_BS,
                    plot_AS,
                    labels = c("a", "b"),
                    common.legend = T)
  
  file_name <- paste0(getwd(), "/results/_paper/figures/MSR_FCTX_subsampling")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 90, units = "mm", dpi = 300)
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
    ylim(c(0,23000))+
    theme_light() +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
    custom_ggtheme +
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
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, nrow = 1))
  
  
  ggpubr::ggarrange(plotMSR_donor,
                    plotMSR_acceptor,
                    nrow = 1,
                    ncol = 2,
                    common.legend = T,
                    labels = c("b", "c"))
  
  
  
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0(getwd(), "/results/_paper/figures/panel4bc")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 55, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 55, units = "mm", dpi = 300)
  
  ####################################
  
  
  
  plotMSR_donor_zoomed <- ggplot(data = df_introns_biotype_tidy %>% 
                                   filter(MSR_type == "MSR_D")) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    ggtitle("MSR Donor") +
    xlab("Mis-splicing ratio value group") +
    ylab("") +
    scale_y_continuous(limits =c(0,750), position = "right") +
    theme_light() +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
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
    scale_y_continuous(limits =c(0,750), position = "right") +
    theme_light() +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding","Non-protein-coding")) +
    custom_ggtheme +
    guides(fill = guide_legend(title = NULL, 
                               ncol = 2, nrow = 1))
  plotMSR_acceptor_zoomed
  
  
  ggpubr::ggarrange(plotMSR_donor_zoomed,
                    plotMSR_acceptor_zoomed,
                    nrow = 1,
                    ncol = 2,
                    common.legend = T,
                    labels = c("b", "c"))
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0(getwd(), "/results/_paper/figures/panel4bc-zoomed")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 55, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 55, units = "mm", dpi = 300)
  
  
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
  
  file_name <- paste0(paste0(getwd(), "/results/_paper/results/05_supplementary_data_MSR_all_tissues.csv"))
  
  
  if ( !file.exists(file_name) ) {
    
    MSR_all_tissues <- map_df(all_projects, function(project_id) {
      
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
                   H1 = c("Introns from PC transcripts are less miss-spliced at the donor than at the acceptor", 
                                  "Introns from Non-PC transcripts are less miss-spliced at the donor than at the acceptor", 
                                  "Introns from Non-PC transcripts are more mis-spliced at the donor than introns from PC transcripts",
                                  "Introns from Non-PC transcripts are more mis-spliced at the acceptor than introns from PC transcripts"),
                   V = c(wilcox_PC$statistic %>% unname(),
                         wilcox_NPC$statistic %>% unname(),
                         wilcox_D$statistic %>% unname(),
                         wilcox_A$statistic %>% unname()),
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
    
    write.csv(x = MSR_all_tissues, file = file_name, row.names = F)
    
  } else {
    MSR_all_tissues <- read.csv(file = file_name)
  }
  
}




## 7. Local sequence conservation is the most important predictor of mis-splicing

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
  ## GET DATA FROM THE DATABASE
  
  query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
  introns %>% nrow()
  
  query <- paste0("SELECT ref_junID, ref_length, ref_ss5score, ref_ss3score, 
                  ref_cons5score, ref_cons3score, ref_CDTS5score, ref_CDTS3score, protein_coding
                  FROM 'intron' 
                  WHERE ref_junID IN (", paste(introns$ref_junID, collapse = ","),")")
  introns <- introns %>%
    inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = "ref_junID") %>% 
    as_tibble() 
  introns %>% nrow()
  
  ## JOIN WITH TRANSCRIPT DATA
  query <- paste0("SELECT * FROM 'transcript' WHERE id IN (", paste(introns$transcript_id, collapse = ","),")")
  introns <- introns %>%
    inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
              by = c("transcript_id" = "id")) %>% 
    as_tibble() 
  introns %>% nrow()
  
  
  ## JOIN WITH GENE DATA
  query <- paste0("SELECT *
                  FROM 'gene' WHERE id IN (",
                  base::paste(introns$gene_id, collapse = ","),")")
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
                  #mean_phastCons20way_5ss = ref_cons5score,
                  #mean_phastCons20way_3ss = ref_cons3score,
                  MSR_D,
                  MSR_A)
  
  
  idb %>% nrow()
  
  idb <-idb %>% 
    left_join(y = db_introns_tidy %>% 
                as_tibble() %>%
                distinct(ref_junID, .keep_all = T) %>%
                dplyr::select(ref_junID, 
                              mean_phastCons20way_5ss = phastCons20way_3ss_mean,
                              mean_phastCons20way_3ss = phastCons20way_5ss_mean),
              by = "ref_junID")
  
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
                               robust = T,
                               #n.sd = 2,
                               #pvals = TRUE,
                               legend.title = "",
                               #plot.distributions = TRUE,
                               ci_level = 0.95,
                               coefs = coef_names,
                               colors = c("#35B779FF","#64037d"),
                               groups =  list("Gene Level" = c("Gene Length",
                                                               "Gene TPM",
                                                               "Gene num. transcripts",
                                                               "Protein coding"),
                                              "Intron Level" = c("Intron Length",
                                                                 "Intron 5'ss MES score",
                                                                 "Intron 3'ss MES score",
                                                                 "CDTS 5'ss",
                                                                 "CDTS 3'ss",
                                                                 "PhastCons20 5'ss",
                                                                 "PhastCons20 3'ss")),
                               #facet.cols = 2,
                               facet.label.pos = "left",
                               model.names = model_names,
                               plot.distributions = F,
                               rescale.distributions = T,
                               point.shape = T) +
    guides(colour = guide_legend(ncol = 2, nrow = 1)) +
    theme_light()  +
    facet_col(~group, 
              space = "free", 
              scales = "free_y",
              strip.position = "left") + 
    custom_ggtheme +
    ylab("Covariates") +
    geom_hline(yintercept = seq(from = 0,
                                to = length((fit_donor$coefficients %>% names)[-1]) + .5,
                                by = 1),
               linewidth = 0.2,
               colour = "#999999" )
  
  plotLM
  
  tiles_data <- rbind(summary(fit_donor)$coefficients[,4] %>% 
                        as_tibble() %>%
                        mutate(q = p.adjust(value,method = "fdr")) %>%
                        mutate(value = ifelse(q > 0.05, NA, q)) %>%
                        mutate(value = ifelse(value == 0, 2.2e-16, value)) %>%
                        mutate(q = value %>% log10())%>%
                        mutate(names = names(summary(fit_donor)$coefficients[,4] ),
                               type = "MSR Donor") %>%
                        dplyr::rename(pvalue = value) %>%
                        filter(names != "(Intercept)") %>%
                        rev(),
                      
                      summary(fit_acceptor)$coefficients[,4] %>% 
                        as_tibble() %>%
                        mutate(value = ifelse(value == 0, 2.2e-16, value)) %>%
                        mutate(q = p.adjust(value,method = "fdr")) %>%
                        mutate(value = ifelse(value > 0.05, NA, value)) %>%
                        mutate(q = q %>% log10())%>%
                        mutate(names = names(summary(fit_acceptor)$coefficients[,4] ),
                               type = "MSR Acceptor") %>%
                        dplyr::rename(pvalue = value) %>%
                        filter(names != "(Intercept)") %>%
                        rev()) 
  
  tiles_data <- tiles_data %>%
    inner_join(y = data.frame(col_name = c("gene_length","gene_tpm","gene_num_transcripts","protein_coding",
                                           "intron_length","intron_5ss_score","intron_3ss_score",
                                           "CDTS_5ss","CDTS_3ss","mean_phastCons20way_5ss","mean_phastCons20way_3ss"),
                              col_label = c("Gene Length","Gene TPM","Gene num. transcripts","Protein coding",
                                            "Intron Length","Intron 5'ss MES score","Intron 3'ss MES score","CDTS 5'ss","CDTS 3'ss",
                                            "PhastCons20 5'ss","PhastCons20 3'ss"),
                              group_label = c("Gene Level","Gene Level","Gene Level","Gene Level",
                                              "Intron Level","Intron Level","Intron Level","Intron Level","Intron Level",
                                              "Intron Level","Intron Level")) %>% as_tibble(),
               by = c("names" = "col_name"))
  
  tiles_data$col_label <- factor(tiles_data$col_label, levels= c("Gene Length","Gene TPM","Gene num. transcripts","Protein coding",
                                                         "Intron Length","Intron 5'ss MES score","Intron 3'ss MES score","CDTS 5'ss","CDTS 3'ss",
                                                         "PhastCons20 5'ss","PhastCons20 3'ss") %>% rev())
  tiles_data$type <- factor(tiles_data$type, levels= c("MSR Donor",
                                                         "MSR Acceptor"))

  
  ggpubr::ggarrange(plotLM ,
                    ggplot() + 
                      geom_tile(data = tiles_data %>%
                                  dplyr::rename("log10(q)"=q), 
                                mapping = aes(x = type, y = col_label, fill = `log10(q)`, colour = "q>0.05")) +
                      scale_fill_gradient(low = "red",high = "white", na.value = '#cccccc') +
                      scale_colour_manual(values = c( "q>0.05" = "#cccccc")) +
                      xlab(" ") + 
                      ylab("") + 
                      facet_col(~group_label,
                                space = "free",
                                scales = "free_y",
                                strip.position = "left") +
                      theme_light() +
                      custom_ggtheme  + 
                      theme(legend.box = "horizontal",
                            legend.box.margin=margin(b = -11,l = -20))+
                      guides(colour = guide_legend(override.aes = list(fill = '#cccccc'),
                                                   title = NULL,
                                                   label.position = "bottom",
                                                   order = 2)),
                    ncol = 2,
                    nrow = 1)

  
  file_name <- paste0(getwd(), "/results/_paper/figures/panel4c")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 180, height = 70, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 70, units = "mm", dpi = 300)
  
  

}




## 8. Accuracy in splicing is affected by RNA-binding protein (RBP) expression changes

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
  
  file_name <- paste0(getwd(), "/results/_paper/results/common_introns_all_tissues.rds")
  saveRDS(object = introns, file = file_name)
  
}

get_estimate_variance_across_tissues <- function() {
  
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  all_projects <- df_metadata$SRA_project %>% unique()
  
  ## LOAD COMMON INTRONS
  common_introns <- readRDS(file = paste0(getwd(), "/results/_paper/results/common_introns_all_tissues.rds"))

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
      ## GET Median TPM value across the tissues
      
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
      ## Get splicing data for the current tissue corresponding to the set of common introns
      
      query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      ## Only common introns across tissues
      introns <- introns %>%
        filter(ref_junID %in% common_introns$ref_junID) %>% 
        distinct(ref_junID, .keep_all = T)
      
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
      introns <- introns %>%
        inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                   by = "ref_junID") %>% 
        as_tibble() 
      
      
      ## JOIN WITH TRANSCRIPT DATA
      query <- paste0("SELECT * FROM 'transcript' WHERE id IN (", paste(introns$transcript_id, collapse = ","),")")
      introns <- introns %>%
        inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                  by = c("transcript_id" = "id")) %>% 
        as_tibble() 
      
      
      ## JOIN WITH GENE DATA
      query <- paste0("SELECT *
                  FROM 'gene' WHERE id IN (",
                      base::paste(introns$gene_id, collapse = ","),")")
      introns <- introns %>%
        inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                   by = c("gene_id" = "id")) %>% 
        as_tibble()  %>%
        dplyr::rename(gene_length = gene_width)
      
      
      ## JOIN WITH TPM DATA
      introns <- introns %>%
        left_join(y = tpm,
                  by = c("gene_id.y" = "gene_id"))
      
      introns <- introns %>%
        dplyr::rename(gene_num_transcripts = n_transcripts )
      
      
      ####################################
      ## LINEAR MODELS
      
      ## 1. Model the MSR at the Donor
      
      donor_model <- lm(MSR_D ~
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
      donor_model %>% summary()

      #MSR_Donor_list[[tissue]] <- donor_model

      MSR_Donor <- data.frame(feature = donor_model$coefficients %>% names(),
                              estimate = donor_model$coefficients %>% unname()) %>%
        inner_join(y = ((donor_model %>% summary())$coefficients) %>% 
                     as_tibble(rownames = "feature") %>%
                     dplyr::select(feature,`Pr(>|t|)`),
                   by = "feature") %>%
        mutate(tissue = cluster_id,
               type = "MSR_Donor") %>%
        dplyr::rename(pval = `Pr(>|t|)`)

      
      
      ## Acceptor
      
      ## 2. Model the MSR at the Acceptor splice site
      
      acceptor_model <- lm(MSR_A ~
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
      acceptor_model %>% summary()

      MSR_Acceptor <- data.frame(feature = acceptor_model$coefficients %>% names(),
                              estimate = acceptor_model$coefficients %>% unname()) %>%
        inner_join(y = ((acceptor_model %>% summary())$coefficients) %>% 
                     as_tibble(rownames = "feature") %>%
                     dplyr::select(feature,`Pr(>|t|)`),
                   by = "feature") %>%
        mutate(tissue = cluster_id,
               type = "MSR_Acceptor") %>%
        dplyr::rename(pval = `Pr(>|t|)`)


      return( rbind(MSR_Donor, MSR_Acceptor) )
      
    })
  })
  
  
  #############################
  ## SAVE RESULTS
  #############################
  
  file_name <- paste0(getwd(), "/results/_paper/results/variance_estimate_lm_all_tissues.rds")
  saveRDS(object = df_estimates, file = file_name)
  
}

plot_estimate_variance_across_tissues <- function() {
  

  ## 1. Load the estimate variance across tissues
  df_estimate <- readRDS(file = paste0(getwd(), "/results/_paper/results/variance_estimate_lm_all_tissues.rds")) %>%
    group_by(type) %>%
    mutate(q = p.adjust(pval, method = "fdr")) %>%
    ungroup()
  
  ## Check the number of non-significant covariates
  df_estimate %>%
    filter(q <= 0.05)
  df_estimate %>%
    filter(q > 0.05) 
  
  ## We only keep significant covariates
  
  df_estimate <- df_estimate %>%
    filter(q <= 0.05) %>%
    dplyr::select(-c(q, pval))
  
  
  ## Separate the estimate values by MSR
  
  MSR_Donor <- df_estimate %>%
    filter(type == "MSR_Donor") %>%
    mutate(type = "MSR Donor") %>%
    spread(key = feature, value = estimate)
  
  MSR_Acceptor <- df_estimate %>%
    filter(type == "MSR_Acceptor") %>%
    mutate(type = "MSR Acceptor") %>%
    spread(key = feature, value = estimate)
  
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
                               subgroup = c( "Gene Level",
                                             "Gene Level",
                                             "Gene Level",
                                             "Gene Level",
                                             "Intron Level",
                                             "Intron Level",
                                             "Intron Level", 
                                             "Intron Level",
                                             "Intron Level",
                                             "Intron Level",
                                             "Intron Level"))
  
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
  plotTissuesMSRDonor <- ggplot(data = MSR_Donor_tidy, aes(tissue, feature, fill = feature)) + 
    geom_boxplot(fill = "#35B779FF") +
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
    theme(axis.text.y = element_text(#angle = 70, 
                                     vjust = 0.5,
                                     hjust = 1)) +
    geom_hline(yintercept = 0,linetype='dotted')
  
  
  plotTissuesMSRDonor
  file_name <- paste0(getwd(), "/results/_paper/figures/panel5a")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 57, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 57, units = "mm", dpi = 300)
  
  
  
  plotTissues3ssLM <- ggplot(data = MSR_Acceptor_tidy, aes(tissue, feature, fill = feature)) + 
    geom_boxplot(fill = "#8d03b0") +
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
    theme(axis.text.y = element_text(#angle = 70, 
      vjust = 0.5,
      hjust = 1)) +
    geom_hline(yintercept = 0,linetype='dotted')
  
  
  plotTissues3ssLM
  file_name <- paste0(getwd(), "/results/_paper/figures/panel5b")
  ggplot2::ggsave(filename = paste0(file_name, ".svg"), width = 180, height = 57, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(file_name, ".png"), width = 180, height = 57, units = "mm", dpi = 300)
  
  
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
  
  
  #######################################
  ## STATS
  #######################################
  
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
    pull(q) %>% 
    summary()
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_3ss",
           type == "MSR_Donor") %>%
    ungroup() %>%
    pull(q) %>% summary()
  
  
  ## Acceptor - pval
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_5ss",
           type == "MSR_Acceptor") %>%
    ungroup() %>%
    pull(q) %>% 
    summary()
  
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_3ss",
           type == "MSR_Acceptor") %>%
    ungroup() %>%
    pull(q) %>% summary()
  
  
  
  
  
  ## Donor - effect size
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_5ss",
           type == "MSR_Donor") %>%
    ungroup() %>%
    pull(estimate) %>% 
    summary()
  df_estimate %>%
    group_by(type) %>%
    filter(feature == "mean_phastCons20way_3ss",
           type == "MSR_Donor") %>%
    ungroup()  %>%
    pull(estimate) %>% 
    summary()
  
  
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
  
}

compare_tissues_somatic_mutations <- function(project_id1 = "SKIN",
                                              cluster_id1 = "Skin - Sun Exposed (Lower leg)",
                                              project_id2 = "SKIN",
                                              cluster_id2 = "Skin - Not Sun Exposed (Suprapubic)",
                                              stats = F) {
  # project_id1 = "SKIN"
  # cluster_id1 = "Skin - Sun Exposed (Lower leg)"
  # project_id2 = "SKIN"
  # cluster_id2 = "Skin - Not Sun Exposed (Suprapubic)"
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
  
  df_database_introns_tidy<-  df_database_introns_tidy %>%
    mutate(tissue = str_remove(string = tissue,
                               pattern = "_SKIN"))
  
  ## VISUALISE MEAN READ COVERAGE BEFORE SUBSAMPLING
  plot_coverage_bf <- ggplot(data = df_database_introns_tidy %>% distinct(ref_junID,.keep_all = T)) +
    geom_density(mapping = aes(x = mean_coverage, fill = tissue), alpha = 0.8) +
    ggtitle("Before subsampling") +
    theme(legend.position = "top",
          legend.text = element_text(size = 11),
          text = element_text(size = 12)) +
    scale_fill_discrete(name = "") +
    xlab("log10 mean read coverage")
  plot_coverage_bf
  
  ## VISUALISE MEAN READ COVERAGE AFTER SUBSAMPLING
  plot_coverage_af <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = tissue), alpha = 0.8) +
    ggtitle("After subsampling") +
    #ggtitle("Mean read coverage per annotated intron across all samples\nfrom 54 GTEx v8 tissues - Subsampling performed.") +
    scale_fill_discrete(name = "") +
    theme(legend.position = "top",
          legend.text = element_text(size = 11),
          text = element_text(size = 12)) +
    xlab("log10 mean read coverage")
  plot_coverage_af
  
  ggpubr::ggarrange(plot_coverage_bf,
                    plot_coverage_af,
                    labels = c("a", "b"),
                    common.legend = T)
  
  file_name <- paste0(getwd(), "/results/_paper/figures/", project_id1, "_", project_id2, "_somatic_mutations_coverage")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 183, height = 150, units = "mm", dpi = 300)
  
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
    
    ## splicing noise at the donor is more frequent in introns from  vs "Skin - Sun Exposed (Lower leg)" 
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
    
    ## splicing noise at the acceptor 
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

plot_effect_size_data <- function() {
  
  source(paste0(getwd(), "/code/helper_functions.R"))
  
  category_labels <- c("Splicing_regulation" = "Splicing regulation", 
                       "Spliceosome"="Spliceosome", 
                       "Exon_junction_complex" = "EJC", 
                       "NMD" = "NMD")
  max_genes <- 7
  
  ## LOAD RBPS ---------------
  
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

  
  ## LOAD COMMON INTRONS ----------------
  
  overwrite = F
  num_cores = 10 
  global_common_introns_path = paste0(getwd(), "/code/variables/global_common_introns_filtered.rds")
  global_common_novel_path = paste0(getwd(), "/code/variables/global_common_novel_filtered.rds")
  
  RBPs_path <- path.expand("/home/grocamora/RytenLab-Research/03-ENCODE_RBP_Automatization/RBPs")
  
  # Load global common introns
  if( !file.exists(global_common_introns_path) ) {
    global_common_introns <- generateCommonIntronsParallel(target_RBPs, RBPs_path, metadata_RBPs, 
                                                           required_clusters, global_common_introns_path, 
                                                           overwrite, num_cores)
  }else{
    global_common_introns <- readRDS(global_common_introns_path)
  }
  
  # Load global common novel
  if ( !file.exists(global_common_novel_path) ) {
    common_ref_coordinates <- global_common_introns %>% pull(ref_coordinates) %>% unique
    global_common_novel <- generateCommonNovelParallel(target_RBPs, common_ref_coordinates, RBPs_path, 
                                                       metadata_RBPs, required_clusters, global_common_novel_path, 
                                                       overwrite, num_cores)
  }else{
    global_common_novel <- readRDS(global_common_novel_path)
  }
  
  
  ## GET THE EFFECT SIZE ---------------------
  
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
  
  ## Execute the wilcox test on MSR_A and MSR_D --------------------------
  ## Execute the wilcox test on MSR_A and MSR_D
  if ( !file.exists(paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.rds")) ) {
    
    
    MSR_D_tests <- generateMSRtests(target_RBPs = target_RBPs, 
                                    MSR = MSR_D, 
                                    file_output = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.rds"), 
                                    overwrite = T,
                                    num_cores = num_cores)
    
    ## Add the categories
    MSR_D_tests <- addMSRcategories(MSR_D_tests, metadata_RBPs)
    
    ## Add bonferroni correction
    MSR_D_tests <- addBonferroniCorrection(MSR_D_tests)
    
    
    saveRDS(object = MSR_D_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.rds"))
    xlsx::write.xlsx2(x = MSR_D_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.xlsx"), 
                      sheetName = "MSRD_test", row.names = F, append = T)
  } else {
    MSR_D_tests <- readRDS(file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRD.rds"))
  }
  if ( !file.exists(paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.rds")) ) {
    
    MSR_A_tests <- generateMSRtests(target_RBPs, MSR = MSR_A, file_output = "", overwrite = T, num_cores = num_cores)
    
    ## Add the categories
    MSR_A_tests <- addMSRcategories(MSR_A_tests, metadata_RBPs)
    
    ## Add bonferroni correction
    MSR_A_tests <- addBonferroniCorrection(MSR_A_tests)
    
    saveRDS(object = MSR_A_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.rds"))
    xlsx::write.xlsx2(x = MSR_A_tests, file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.xlsx"), 
                     sheetName = "MSRA_test", row.names = F, append = T)
  } else {
    MSR_A_tests <- readRDS(file = paste0(getwd(), "/results/_paper/results/ENCODE_effectsize_MSRA.rds"))
  }
  
  # Combine both MSR_A and MSR_D -------------------
  MSR_combined = rbind(MSR_A_tests %>% mutate(MSR_type = "MSR_A"), 
                       MSR_D_tests %>% mutate(MSR_type = "MSR_D")) %>%
    filter(p.value.bonferroni <= 0.05)
  
  # Filter the Splicing regulation category
  filter_splicing_regulation <- MSR_combined %>% 
    dplyr::select(-statistical_test, -H0, -H1) %>% 
    arrange(-effect_size) %>%
    filter(Category == "Splicing_regulation") %>%
    distinct(target_gene, .keep_all = T) %>%
    head(max_genes) %>%
    pull(target_gene)
  
  MSR_graph_data <- MSR_combined %>%
    filter(target_gene %in% filter_splicing_regulation | Category != "Splicing_regulation") %>%
    arrange(-effect_size) %>%
    mutate(target_gene = factor(target_gene, levels = .$target_gene %>% unique),                                           # Use factors to sort the graph
           Category = factor(Category, levels = c("Splicing_regulation", "Spliceosome", "NMD", "Exon_junction_complex")))  # Use factors to sort the graph
  
  
  ######
  ## METADATA KEFF
  ######
  
  project_path <- "/home/grocamora/RytenLab-Research/09-ENCODE_Gene_expression/"
  metadata_kEff_path <- paste0(project_path, "variables/metadata_kEff.tsv")
  metadata_kEff <- readr::read_delim(metadata_kEff_path, show_col_types = F) %>%
    mutate(kEff_text = ifelse(is.na(kEff_avg), kEff_avg, paste0(round(kEff_avg), "%"))) %>%
    mutate(kEff_text = kEff_text %>% as.factor())
  
  
  MSR_graph_data <- MSR_graph_data %>%
    inner_join(y = metadata_kEff,
               by = "target_gene")
  MSR_graph_data %>% head()
  
  ######
  ## PLOT
  ######
  
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
    labs(x = "Target gene shRNA knockdown", y = "Probability of superior MSR in\ngene knockdown vs. untreated samples") + 
    facet_row(facets = vars(Category), 
               scales = "free_x", space = "free",
               labeller = labeller(Category = category_labels),
               drop = T,
              shrink = T) +
    scale_fill_manual(values = c("#35B779FF","#64037d"),
                      labels = c("MSR_A" = "MSR Acceptor", "MSR_D" = "MSR Donor"),
                      breaks = c("MSR_D", "MSR_A"))  +
    guides(fill = guide_legend(title = NULL, order = 2, ncol = 2,  nrow = 1 )) +
    theme_light() 
    
  plot_effectsize
  
  plot_effectsize +
    ggnewscale::new_scale_fill() +
    #geom_hline(data = data.frame(hline = c(0.1, 0.3, 0.5)), aes(yintercept = hline), linewidth = 0.25) +
    geom_tile(stat = "identity", aes(y = 0.68, fill = kEff_avg, color = "No data\navailable"), 
              linewidth = 0.5, width = 1, height = 0.045) + 
    geom_text(aes(y = 0.76, label = kEff_text),
              color = "black", 
              size = 3) +
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
    theme(axis.text.x = element_text(angle = 90),
          legend.box = "horizontal") +
    guides(colour = guide_legend(override.aes = list(fill = '#999999'),
                                 title = NULL,
                                 label.position = "bottom",
                                 order = 1))

  
  # Save the graph
  ggsave(file = paste0(getwd(), "/results/_paper/figures/Effect_size_combined_top", max_genes, ".png"), 
         width = 180, height = 90, units = "mm", dpi = 300)
  ggsave(file = paste0(getwd(), "/results/_paper/figures/Effect_size_combined_top", max_genes, ".svg"), 
         width = 180, height = 90, units = "mm", dpi = 300)
}

plot_data_AQR_U2AF2 <- function() {
  
  
  
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
  } else {
    global_common_introns <- readRDS(global_common_introns_path)
    global_common_introns %>% distinct(ref_coordinates) %>% nrow()
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
  
  # scales::show_col(colours = viridis::viridis_pal(option = "A")(n = 20))
  # scales::show_col(colours = viridis::viridis_pal(option = "B")(n = 20))
  # scales::show_col(colours = viridis::viridis_pal(option = "C")(n = 20))
  # scales::show_col(colours = viridis::viridis_pal(option = "D")(n = 20))
  # scales::show_col(colours = viridis::viridis_pal(option = "E")(n = 20))
  # scales::show_col(colours = viridis::viridis_pal(option = "F")(n = 20))
  
  #################################################
  ## PLOT AQR AND U2AF1
  #################################################
  
  
  target_RBPs = c("U2AF2","AQR")
  
  RBP_intron <- global_common_introns %>% filter(target_gene %in% target_RBPs)
  RBP_novel <- global_common_novel %>% filter(target_gene %in% target_RBPs)
  
  
  RBP_novel %>%
    dplyr::count(target_gene,cluster, novel_type)
  
  RBP_intron <- RBP_intron %>%
    mutate(cluster = ifelse(cluster == "case", "gene knockdown", "control"))%>%
    mutate(cluster = factor(cluster, levels = c("control", "gene knockdown")))
  RBP_novel <- RBP_novel %>%
    mutate(cluster = ifelse(cluster == "case", "gene knockdown", "control")) %>%
    mutate(cluster = factor(cluster, levels = c("gene knockdown","control" )))
  
  ## DISTANCES
  
  limit_bp = 30
  
  # RBP_novel %>%
  #   filter(novel_type == "novel_acceptor") %>%
  #   mutate(novel_type = str_replace(string = novel_type,
  #                                   pattern = "_",
  #                                   replacement = " ")) %>%
  #   filter(abs(distance) < limit_bp,
  #          abs(distance) == 1)
  
  setdiff(x = RBP_intron %>%
            filter(target_gene == "AQR") %>%
            pull(ref_coordinates) %>%
            unique,
          y = RBP_intron %>%
            filter(target_gene == "U2AF2") %>%
            pull(ref_coordinates) %>%
            unique)
  
  
  setdiff(x = RBP_novel %>%
            filter(target_gene == "AQR") %>%
            pull(ref_coordinates) %>%
            unique,
          y = RBP_intron %>%
            filter(target_gene == "AQR") %>%
            pull(ref_coordinates) %>%
            unique)
  setdiff(x = RBP_novel %>%
            filter(target_gene == "U2AF2") %>%
            pull(ref_coordinates) %>%
            unique,
          y = RBP_intron %>%
            filter(target_gene == "U2AF2") %>%
            pull(ref_coordinates) %>%
            unique)
  

  
  plot_aqr_u2af2 <- ggplot(RBP_novel %>%
                             filter(novel_type == "novel_acceptor") %>%
                             mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " ")) %>%
                             filter(abs(distance) < limit_bp) %>%
                             mutate(target_gene = factor(x = target_gene, levels = target_RBPs))) + 
    geom_histogram(aes(x = distance, fill = cluster),
                   bins = 60, 
                   binwidth = 1, 
                   position = "identity", 
                   alpha = 1, 
                   color = "black", 
                   linewidth = 0.1) +
    scale_x_continuous(breaks = seq(-limit_bp, limit_bp, length.out = 5)) + 
    # ggsci::scale_fill_npg()+
    scale_fill_manual(values = c("#91D1C2B2", "#999999"),
                     breaks = c("gene knockdown", "control"),
                     labels = c("gene knockdown  ", "control")) +
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    facet_wrap(vars(target_gene)) +
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
  
  
  distances_acceptor <- plot_aqr_u2af2 / (distance_rectangle +  distance_rectangle) + patchwork::plot_layout(heights = c(8, 1))
  distances_acceptor
  
  #mypal = ggsci::pal_npg("nrc", alpha = 0.7)(10)
  #scales::show_col(mypal)
  
  # ggsave(file = paste0(getwd(), "/results/_paper/figures/distance_stacked_",target_RBPs[1],"_",target_RBPs[2],".png"), 
  #        width = 183, height = 110, dpi = 300, units = "mm")
  
  
  
  
  
  ## DELTA MES -----------------------------------------
  
  RBP_novel <- RBP_novel %>%
    mutate(cluster = factor(cluster, levels = c("control","gene knockdown")))
  

  ## DELTA VALUES FOR THE FIRST RBP
  RBP_merged1 <- RBP_novel %>% 
    filter(target_gene == target_RBPs[1]) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(novel_coordinates, ref_junID, cluster, novel_ss3score, target_gene) %>%
    left_join(RBP_intron %>% 
                filter(target_gene == target_RBPs[1]) %>%
                dplyr::select(ref_junID, ref_ss3score) %>% distinct(),
              by = "ref_junID") %>%
    mutate(delta_ss3score = ref_ss3score - novel_ss3score) %>%
    ungroup()
  RBP_delta1 <- RBP_merged1 %>%
    dplyr::select(target_gene, delta_ss3score, cluster ) %>%
    mutate(delta_ss3score = delta_ss3score %>% as.double()) %>%
    drop_na() %>%
    mutate(delta_type = "Delta MES Acceptor") %>%
    mutate(delta_type = delta_type %>% as.factor())  %>%
    mutate(RBP = target_RBPs[1])
  
  
  ## Delta 3’ MES scores of novel acceptor junctions were significantly lower (pval<9.4e-04) indicating that weaker acceptor splice sites were being used 
  
  wilcox.test(x = RBP_delta1 %>% filter(target_gene == "U2AF2", cluster == "gene knockdown") %>% pull(delta_ss3score),
              y = RBP_delta1 %>% filter(target_gene == "U2AF2", cluster == "control") %>% pull(delta_ss3score),
              paired = F,
              alternative = "greater")

  ## DELTA VALUES FOR THE SECOND RBP
  RBP_merged2 <- RBP_novel %>% 
    filter(target_gene == target_RBPs[2]) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(novel_coordinates, ref_junID, cluster, novel_ss3score, target_gene) %>%
    left_join(RBP_intron %>% 
                filter(target_gene == target_RBPs[2]) %>%
                dplyr::select(ref_junID, ref_ss3score) %>% distinct(),
              by = "ref_junID") %>%
    mutate(delta_ss3score = ref_ss3score - novel_ss3score) %>%
    ungroup()
  RBP_delta2 <- RBP_merged2 %>%
    dplyr::select(target_gene, delta_ss3score,cluster ) %>%
    mutate(delta_ss3score = delta_ss3score %>% as.double()) %>%
    drop_na() %>%
    mutate(delta_type = "Delta MES Acceptors") %>%
    mutate(delta_type = delta_type %>% as.factor()) %>%
    mutate(RBP = target_RBPs[2])
  
  wilcox.test(x = RBP_delta2 %>% filter(target_gene == "AQR", cluster == "gene knockdown") %>% pull(delta_ss3score),
              y = RBP_delta2 %>% filter(target_gene == "AQR", cluster == "control") %>% pull(delta_ss3score),
              paired = F,
              alternative = "greater")
  
  
  ## MERGE DELTA VALUES FOR THE TWO RBPS
  delta_mes_acceptor <- rbind(RBP_delta1, RBP_delta2) %>%
    group_by(RBP, cluster) %>%
    mutate(medianMSR = delta_ss3score %>% median()) %>% 
    ungroup() %>%
    group_by(RBP) %>%
    mutate(p.value = format(x = wilcox.test(delta_ss3score ~ cluster)$p.value,
                            digits = 2, scientific = T)) %>%
    ungroup() %>%
    mutate(target_gene = factor(x = target_gene,levels = target_RBPs)) %>%
    mutate(RBP = factor(x = RBP, levels = target_RBPs)) %>%
    mutate(p.value = ifelse(p.value == "0e+00", "2.2e-100",p.value))
  


  
  delta_mes <- ggplot(data = delta_mes_acceptor)  +
    geom_density(aes(x = delta_ss3score, fill = cluster), 
                 alpha = 0.8, linewidth = 0.3, 
                 color = "black") +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept=medianMSR, colour=cluster),
               linetype="dashed", linewidth=0.9) +
    facet_wrap(vars(RBP)) +
    geom_text(data = delta_mes_acceptor, 
              aes(label = paste0("p = ", p.value)), 
              x=25, y=0.08, family= "Arial", size = 3, 
              colour = "#333333") +
    scale_fill_manual(values = c("#91D1C2B2", "#999999"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("Gene knockdown  ", "Control")) + 
    scale_colour_manual(values = c("#91D1C2B2", "#999999"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("Gene knockdown  ", "Control")) + 
    labs(x = "Delta MES Acceptor", y = "Density") +
    theme_light() +
    guides(fill = guide_legend(title = "Sample type: ",
                               ncol = 2, nrow = 1 ),
           colour = guide_legend(title = "Median Delta MES: ",
                               ncol = 2, nrow = 1 )) +
    custom_ggtheme 
  delta_mes + theme(legend.box = "horizontal")
  
  
  ggpubr::ggarrange(x = distances_acceptor,
                    y = delta_mes + theme(legend.box = "horizontal"),
                    labels = c("a", "b"),
                    ncol = 1,
                    nrow = 2,
                    common.legend = T)
  
  
  
  ggsave(file = paste0(getwd(),"/results/_paper/figures/panel6ab.png"), 
         width = 180, height = 140, dpi = 300, units = "mm")
  ggsave(file = paste0(getwd(),"/results/_paper/figures/panel6ab.svg"), 
         width = 180, height = 140, dpi = 300, units = "mm")
  
 

  ## DISTANCES LONGER ------------------------
  
  RBP_novel <- RBP_novel %>%
    mutate(cluster = factor(cluster, levels = c("gene knockdown","control")))
  
  limit_bp = 200
  target_RBP <- "AQR"
  plot_aqr_long_distances <- ggplot(RBP_novel %>%
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
                   color = "#333333", 
                   linewidth = 0.02) +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp), 
                       breaks = seq(-limit_bp, limit_bp, length.out = 5)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(values =  c("#91D1C2B2", "#999999"),
                      breaks = c("gene knockdown", "control"),
                      labels = c("Gene knockdown  ", "Control")) +
    labs(x = "Distance (bp)", y = "Number of unique novel junctions") + 
    facet_grid(novel_type~target_gene, 
               labeller = labeller(novel_type = c("novel donor" = "Novel donor", "novel acceptor" = "Novel acceptor"))) +
    theme_light() +
    custom_ggtheme +
    theme(legend.position = "none") 
  
  
  distance_rectangle_longer <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 60), fill = "grey", color = "black") +
    geom_text(aes(x = 100, y = 30),  size = 3, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 30, ymax = 31), fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -100, y = 51),  size = 3, label = "intron") +
    theme_void()
  
  
  plot_aqr_long_distances / (distance_rectangle_longer) + patchwork::plot_layout(heights = c(8, 1))
  
  ggsave(file = paste0(getwd(), "/results/_paper/figures/panel6c.png"), 
         width = 180, height = 65, dpi = 300, units = "mm")
  ggsave(file = paste0(getwd(), "/results/_paper/figures/panel6c.svg"), 
         width = 180, height = 65, dpi = 300, units = "mm")
  
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

get_distances_all_tissues()

# compare_tissues_somatic_mutations(project_id1 = "BRAIN",
#                                   cluster_id1 = "Brain - Nucleus accumbens (basal ganglia)",
#                                   project_id2 = "BLOOD",
#                                   cluster_id2 = "Whole Blood")
  