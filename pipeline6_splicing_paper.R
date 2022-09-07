library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(DBI)

####################################################
## GENERATE SPLICING DATABASE ######################
####################################################

source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_SQL_generation.R")

## CONNECT TO THE DATABASE ------------------------------
database_path <- "./database/introverse_qc.sqlite"
con <- dbConnect(RSQLite::SQLite(), "./database/splicing.sqlite")
dbListTables(con)

age <- F

hg38 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
hg_MANE <- rtracklayer::import(con = "/data/references/MANE/MANE.GRCh38.v1.0.ensembl_genomic.gtf")

SRA_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
SRA_projects %>% length() %>% print()



## Remove the tables -----------------------------------------------------------

remove_tables(all = F)
dbListTables(con)


## Create the tables -----------------------------------------------------------

create_master_table(age, QC = T)

create_mane_table()

create_gene_table()

create_intron_table()

create_novel_table()

create_cluster_tables()


prepare_gtex_samples <- function() {
  
  
  if (file.exists("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")) {
    
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
    
    all_projects_used <- NULL
    all_clusters_used <- NULL
    all_samples_used <- NULL
    
    for (project_id in all_projects) {
    
      
      clusters <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                        project_id, "/raw_data/all_clusters.rds"))
      
      print("All clusters saved!")
      
      
      ## SAVE CLUSTERS SP ------------------------------------------------------
      
      clusters_used <- NULL
      
      for (cluster in clusters) {
        
        if (!(cluster %in% c("Brain - Cortex", "Brain - Cerebellum")) &&
            !(project_id %in% c("BREAST", "CERVIX_UTERI", "FALLOPIAN_TUBE", "OVARY",
                                "PROSTATE", "TESTIS", "UTERUS", "VAGINA"))) {
          
          ## Load samples
          samples <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
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
              file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/raw_data/all_clusters_used.rds"))
    }
    
    all_projects_used %>% length()
    all_clusters_used %>% length()
    all_samples_used %>% length()
    
    saveRDS(object = all_projects_used %>% unique(),
            file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
    
  } else {
    all_projects <- c( "ADIPOSE_TISSUE", "ADRENAL_GLAND", "BLADDER", "BLOOD", "BLOOD_VESSEL", "BONE_MARROW", "BRAIN", "BREAST",
                       "CERVIX_UTERI", "COLON", "ESOPHAGUS", "FALLOPIAN_TUBE", "HEART", "KIDNEY", "LIVER", "LUNG",
                       "MUSCLE","NERVE", "OVARY", "PANCREAS", "PITUITARY", "PROSTATE", "SALIVARY_GLAND", "SKIN",
                       "SMALL_INTESTINE", "SPLEEN", "STOMACH", "TESTIS", "THYROID", "UTERUS", "VAGINA") %>% sort()
    saveRDS(object = all_projects,
            file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects.rds")
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
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_summary.rds")
}


########################################
## DATABASE STATS
########################################

tables <- dbListTables(con)

query <- paste0("SELECT * from 'intron'")
dbGetQuery(con, query) %>% distinct(ref_junID) %>% nrow()

query <- paste0("SELECT * from 'novel'")
dbGetQuery(con, query) %>% distinct(ref_junID) %>% nrow()
dbGetQuery(con, query) %>% distinct(novel_junID) %>% nrow()

dbGetQuery(con, query) %>% dplyr::count(novel_type)


query <- paste0("SELECT * from 'master'")
db_metadata <- dbGetQuery(con, query) %>% as_tibble()
db_metadata %>%
  mutate(rin = rin %>% as.double()) %>%
  filter(rin >= 6)

########################################
## FUNCTIONS
########################################


## SECTION 1 ---------------------------------------------

get_contamination_rates <- function(all_tissues = F) {
  
  # all_tissues <- F
  
  if (all_tissues) {
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  } else {
    all_projects <- "BRAIN"
  }
  
  df_contamination <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    
    if (all_tissues) {
      versions <- c("97")
      dates <- c("26-May-2019")
    } else {
      all_clusters <- all_clusters[5]
      versions <- c("76", "81", "90", "97", "104")
      dates <- c("18-Jul-2014", "07-Jul-2015", "28-Jul-2017", "26-May-2019", "19-Mar-2021")
    }
    
    map_df(all_clusters, function(cluster) {
      
      # cluster <- all_clusters[1]
      print(paste0(cluster))
      
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/base_data/", cluster)
      
      map_df(versions, function(version) {
        
        # version <- versions[1]
        print(paste0("v", version))
        
        #db_introns_old <- readRDS(file = paste0(folder_name, "/v", version, "/", cluster, "_db_introns.rds"))
        #db_novel_old <- readRDS(file = paste0(folder_name, "/v", version, "/", cluster, "_db_novel.rds"))
        
        # db_introns_new <- readRDS(file = paste0(folder_name, "/v105/", cluster, "_db_introns.rds"))
        # db_novel_new <- readRDS(file = paste0(folder_name, "/v105/", cluster, "_db_novel.rds"))
        
        ## ENSEMBL v97
        df_old <-  readRDS(file = paste0(folder_name, "/", cluster, "_annotated_SR_details_length_", version,".rds")) %>%
          as_tibble()
        db_introns_old <- df_old %>% 
          filter(type == "annotated")
        db_novel_old <- df_old %>%
          filter(type %in% c("novel_donor", "novel_acceptor"))
        
        ## ENSEMBL v105
        df_new <- readRDS(file = paste0(folder_name, "/", cluster, "_annotated_SR_details_length_105.rds")) %>%
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
        in_annotation %>% nrow()
        
        
        
        out_annotation <- db_novel_new %>%
          filter(junID %in% db_introns_old$junID) %>% 
          distinct(junID, .keep_all = T) 
        out_annotation %>% nrow()
        
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
        
      })
    })
  })
      
      
      
  if (all_tissues) {
    
    if (exists("df_contamination")) {
      saveRDS(object = df_contamination,
              file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/contamination_rates_raw.rds"))
    } else {
      df_contamination <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/contamination_rates_raw.rds"))
    }
    
    
    df_contamination_tidy <- df_contamination %>% 
      dplyr::select(in_annotation, out_annotation, tissue) %>%
      tidyr::gather(key = "type", value = "prop", -tissue ) 
    
    df_contamination_tidy$type = factor(df_contamination_tidy$type, 
                                        levels = c( "in_annotation","out_annotation"))
    
    df_contamination_tidy = df_contamination_tidy %>% 
      ungroup() %>%
      arrange(type , prop) %>%
      mutate(tissue = fct_inorder(tissue))
    
    colours <- ifelse(str_detect(string = as.factor(df_contamination_tidy$tissue), pattern = "Brain"), "red", "black")
    
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
            axis.text.y = element_text(color = colours,
            #                          angle = 70, 
                                        vjust = 0.5,
                                        hjust = 1)
            ) +
      guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) %>%
      return()
    
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/contamination_rates_tissues.svg")
    ggplot2::ggsave(file_name, width = 183, height = 183, 
                    units = "mm", dpi = 300)
    
    ## Stats
    df_contamination_tidy %>% filter(type == "in_annotation") %>% pull(prop) %>% mean()
    rm(df_contamination)
    
  } else {
    
    if (exists("df_contamination")) {
      saveRDS(object = df_contamination,
              file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                            project_id, "/results/pipeline3/contamination_rates.rds"))
    } else {
      df_contamination <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                                project_id, "/results/pipeline3/contamination_rates.rds"))
    }
    
    df_contamination$contamination_rates = factor(df_contamination$contamination_rates, 
                                                  levels = c( 
                                                              
                                                               
                                                              
                                                              "Ensembl_v104 (26-May-2019)",
                                                              "Ensembl_v97 (26-May-2019)",
                                                              "Ensembl_v90 (28-Jul-2017)",
                                                              "Ensembl_v81 (07-Jul-2015)", 
                                                              "Ensembl_v76 (18-Jul-2014)"))
    
    ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
    ggplot(data = df_contamination) +
      geom_bar(aes(x = contamination_rates, 
                   y = in_annotation, fill =contamination_rates ), 
               stat = "identity") + 
      coord_flip() +
      xlab(NULL) +
      ylab("contamination rates (%)") +
      scale_fill_viridis_d(option = "E")  +
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
    
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/contamination_rates_FCTX.svg")
    ggplot2::ggsave(file_name,  width = 183, height = 183,
                    units = "mm", dpi = 300)
    
    rm(df_contamination)
  }
  
  
  
  
}

get_unique_donor_acceptor_jxn <- function() {
  
  ##################################################################################################
  ## Proportion of each type of unique junctions per GTEx tissue.
  ##################################################################################################
  
  df_proportions <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[5]
    
    print(paste0(Sys.time(), " - ", project_id))
    
    clusters_used <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id,
             cluster %in% clusters_used) %>%
      distinct(cluster) %>%
      pull()
    
    print(all_clusters)
    
    if (!identical(clusters_used %>% sort(), all_clusters %>% sort())) {
      print("Error! Both set of clusters are not identical.")
      break;
    }
    
    map_df(all_clusters, function(cluster_id) {
    
      # cluster_id <- all_clusters[1]
      
      ## Print the tissue
      print(paste0(cluster_id))
      
      samples <- df_metadata %>%
        dplyr::count(cluster) %>%
        filter(cluster == cluster_id) %>% 
        pull(n)
      
      #folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id)
      #samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster,  "_samples.rds"))
      #
      
      ####################
      ## GET THE INTRONS
      ####################
      
      query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      
      # introns <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
      #                                  cluster, "/v105/", cluster, "_db_introns.rds"))
      
      ###########################
      ## GET THE NOVEL JUNCTIONS
      ###########################
      
      query <- paste0("SELECT DISTINCT novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
      novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
      query <- paste0("SELECT DISTINCT novel_junID, novel_type FROM 'novel' WHERE novel_junID IN (", paste(novel_junctions$novel_junID, collapse = ","), ")")
      novel_junctions <- merge(x = novel_junctions,
                               y = dbGetQuery(con, query) %>% as_tibble(),
                               by = "novel_junID",
                               all.x = T)
      
      # novel_junctions <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
      #                                          cluster, "/v105/", cluster, "_db_novel.rds"))
      
      
      introns %>% head()
      novel_junctions %>% head()
      
      
      ###############################
      ## Calculate the proportions
      ###############################
      
      annotated_junc <- (introns %>% distinct(ref_junID) %>% nrow()) 
      donor_junc <- (novel_junctions %>% filter(novel_type == "novel_donor") %>% distinct(novel_junID) %>% nrow())
      acceptor_junc <- (novel_junctions %>% filter(novel_type == "novel_acceptor") %>% distinct(novel_junID) %>% nrow())
      
      annotated_prop <- annotated_junc/(annotated_junc + donor_junc + acceptor_junc)
      donor_prop <- donor_junc/(annotated_junc + donor_junc + acceptor_junc)
      acceptor_prop <- acceptor_junc/(annotated_junc + donor_junc + acceptor_junc)
      
      
      
      
      
      # ## Generate the data.frane
      # df <- data.frame(tissue = tissue,
      #                  prop = annotated_prop,
      #                  type = "intron")
      # 
      # df <- rbind(df,
      #             data.frame(tissue = tissue,
      #                        prop = donor_prop,
      #                        type = "donor"))
      # 
      # df <- rbind(df,
      #             data.frame(tissue = tissue,
      #                        prop = acceptor_prop,
      #                        type = "acceptor"))
      
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

  if (exists("df_proportions"))  {
    saveRDS(object = df_proportions,
            file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/paper/unique_donor_acceptor_jxn.rds")
  } else {
    df_proportions <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/paper/unique_donor_acceptor_jxn.rds")
    
    df_proportions %>%
      arrange(desc(annotated_prop))
  }
  
  ###################
  ## STATS
  ###################
  
  df_proportions %>% 
    pull(annotated_junc) %>% mean
  df_proportions %>% 
    pull(donor_junc) %>% mean
  df_proportions %>% 
    pull(acceptor_junc) %>% mean
  
  #df_proportions <- df_proportions %>% 
  #cor.test(x = df_proportions$acceptor_prop,
  #         y = df_proportions$samples)
  
  df_proportions <- df_proportions %>%
    dplyr::select(tissue, 
                  donor = donor_prop, 
                  acceptor = acceptor_prop, 
                  annotated_intron = annotated_prop) %>%
    tidyr::gather(key = "type", value = "prop", -tissue ) 
  
  df_proportions %>% head()
  
  # ## Stats
  # df_proportions %>%
  #   mutate(annotated_raw = annotated_junc * samples) %>% pull(annotated_raw) %>% mean
  #   arrange(desc(annotated_raw))
  # df_proportions %>%
  #   mutate(donor_raw = donor_junc * samples) %>% pull(donor_raw) %>% mean
  #   arrange(desc(donor_raw))
  # df_proportions %>%
  #   mutate(acceptor_raw = acceptor_junc * samples) %>% pull(acceptor_raw) %>% mean
  #   arrange(desc(acceptor_raw))
    
    
    
  ##############################
  ## PLOT THE DATA
  ##############################
  
  ## First, order by junction type
  df_proportions$type <- factor(df_proportions$type, 
                             levels = c( "acceptor", 
                                         "annotated_intron",
                                         "donor"))
  
  df_proportions_final <- df_proportions %>% 
    ungroup() %>%
    arrange(type , prop) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df_proportions_final$tissue), 
                               pattern = "Brain"), "red", "black")
  
  
  ggplot(df_proportions_final %>% filter(type %in% c("donor", "acceptor")), 
         aes(x = tissue, y = prop * 100, 
             group = type, fill = type)) +
    geom_col(position = position_dodge()) +
    #geom_text(aes(label = round(x = prop, digits = 2)), 
    #          position = position_stack(vjust = 0.5, reverse = TRUE), colour = "white", size = 1.5) +
    coord_flip() +
    theme_light() +
    ylab("unique num. junctions across samples (%)") +
    xlab("") +
    #scale_fill_viridis_d()  +
    #ggtitle("Percentage of unique introns and unique novel junctions per tissue") +
    #scale_x_discrete("tissue", breaks = breaks) +
    
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          #axis.text.x = element_text(color = colours),
          axis.text.y = element_text(color = colours,
                                     #angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          axis.title.x = element_text(colour = "black", size = "10"),
          legend.text = element_text(size = "9"),
          plot.title = element_text(size = "9"),
          legend.title = element_text(size = "9")) +
    
    #scale_fill_manual(values =  c("#0D0887FF","#EF7F4FFF","#A41F9AFF"
    scale_fill_manual(values = c( "#0D0887FF", "#EF7F4FFF"),
                      breaks = c( "acceptor", "donor"),# "annotated_intron"),
                      labels = c( "novel acceptor", "novel donor")) + #, "annotated intron")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 1, 
                               nrow = 3))
  
  
  ## Save the figure 3
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/percent_unique_donor_acceptor.svg"
  ggplot2::ggsave(filename = file_name, 
                  width = 183, height = 183, units = "mm", dpi = 300)
}

get_unique_donor_acceptor_reads <- function() {
  
  df_mean_counts <- map_df(all_projects, function(project_id) {
    
    #############################
    ## GET THE CLUSTERS
    #############################
    
    print(paste0(Sys.time(), " - ", project_id))
    
    clusters_used <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id,
             cluster %in% clusters_used) %>%
      distinct(cluster) %>%
      pull()
    
    if (!identical(clusters_used %>% sort(), all_clusters %>% sort())) {
      print("Error! Both set of clusters are not identical.")
    }
    
    
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
      query <- paste0("SELECT novel_junID, novel_type FROM 'novel' WHERE novel_junID IN (", paste(novel_junctions$novel_junID, collapse = ","), ")")
      novel_junctions <- merge(x = novel_junctions,
                               y = dbGetQuery(con, query) %>% as_tibble(),
                               by = "novel_junID",
                               all.x = T) %>% as_tibble()
      
    
      ###########################
      ## GET THE PROPORTIONS
      ###########################

      
      annotated <- introns %>%
        dplyr::distinct(ref_junID, .keep_all = T) %>%
        mutate(ref_mean_counts = ref_sum_counts / ref_n_individuals ) %>%
        pull(ref_sum_counts) %>% 
        mean()
      
      ###########################################################
      
      introns %>%
        dplyr::distinct(ref_junID, .keep_all = T) %>%
        pull(ref_sum_counts) %>%
        summary()
      novel_junctions %>%
        filter(novel_type == "novel_acceptor") %>%
        dplyr::distinct(novel_junID, .keep_all = T) %>%
        pull(novel_sum_counts) %>% 
        summary()
      novel_junctions %>%
        filter(novel_type == "novel_donor") %>%
        dplyr::distinct(novel_junID, .keep_all = T) %>%
        pull(novel_sum_counts) %>% 
        summary()
      
      ###########################################################
      
      acceptor <- novel_junctions %>%
        filter(novel_type == "novel_acceptor") %>%
        dplyr::distinct(novel_junID, .keep_all = T) %>%
        mutate(novel_mean_counts = novel_sum_counts / novel_n_individuals ) %>%
        pull(novel_sum_counts) %>% 
        sum()
      
      
      
      
      donor <- novel_junctions %>%
        filter(novel_type == "novel_donor") %>%
        dplyr::distinct(novel_junID, .keep_all = T) %>%
        mutate(novel_mean_counts = novel_sum_counts / novel_n_individuals ) %>%
        pull(novel_sum_counts) %>% 
        sum()
      
      annotated_p = annotated * 100 / (annotated + acceptor + donor)
      acceptor_p = acceptor * 100 / (annotated + acceptor + donor)
      donor_p = donor * 100 / (annotated + acceptor + donor)
      
      return(data.frame(tissue = cluster_id,
                        samples = samples,
                        type = c("annotated","acceptor", "donor"),
                        prop_counts = c(annotated_p, acceptor_p, donor_p),
                        mean_counts = c(annotated, acceptor, donor)))
    
    })
    
  })
  
  # ## Save data
  
  if (exists("df_mean_counts"))  {
    saveRDS(object = df_mean_counts,
            file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/paper/unique_donor_acceptor_reads.rds")
  } else {
    df_mean_counts <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/paper/unique_donor_acceptor_reads.rds") %>%
      as_tibble()
    
    df_mean_counts %>%
      filter(type == "annotated") %>%
      pull(mean_counts) %>% mean
    df_mean_counts %>%
      filter(type == "acceptor") %>%
      pull(mean_counts) %>% mean
    
    
    df_mean_counts %>%
      filter(type == "annotated") %>%
      pull(prop_counts ) %>% mean
    df_mean_counts %>%
      filter(type == "acceptor") %>%
      pull(prop_counts ) %>% mean
    df_mean_counts %>%
      filter(type == "donor") %>%
      pull(prop_counts ) %>% mean
  }
  
  # saveRDS(df_mean_counts, file = "/home/sruiz/PROJECTS/splicing-project/results/paper/figure1_data.rds")
  
  # df_mean_counts <- readRDS("/home/sruiz/PROJECTS/splicing-project/results/paper/figure1_data.rds")
  #all_clusters <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_clusters_used.rds")
  
  ## Spread
  df_mean_counts <- df_mean_counts %>% 
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
  
  colours <- ifelse(str_detect(string = as.factor(df_mean_counts_tidy$tissue), pattern = "Brain"), "red", "black")
  
  
  ggplot(data = df_mean_counts_tidy %>% filter(type %in% c("acceptor", "donor"))) +
    geom_col(mapping = aes(x = tissue, y = mean, fill = type), 
             position = position_dodge(), 
             alpha = 0.9) +
    coord_cartesian(ylim=c(0,1)) +
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
          axis.text.y = element_text(color = colours,
                                     #angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9")) +
    guides(fill = guide_legend(title = NULL, ncol = 1, nrow = 3)) #+
  #facet_zoom(ylim=c(0,10))
  
  
  # library(scales)
  # show_col(viridis_pal()(20)) 
  # viridis_pal()(20)
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/unique_donor_acceptor_reads.svg")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  write.csv(x = df_mean_counts,
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/unique_donor_acceptor_reads.csv", row.names = FALSE)
  
}

get_mean_read_count_FTCX <- function() {
  
  tables <- c("Brain - Frontal Cortex (BA9)_BRAIN_misspliced")#,
              #"Brain - Frontal Cortex (BA9)_BRAIN_nevermisspliced")
  
  annotated_sum_count <- NULL
  novel_sum_count <- NULL
  novel_donor_sum_count <- NULL
  novel_acceptor_sum_count <- NULL
  
  query <- paste0("SELECT novel_junID, novel_type from 'novel'")
  db_novel <- dbGetQuery(con, query) %>% as_tibble()
  
  for (table in tables) {
    
    if (!(table %in% c("intron", "novel", "gene", "mane", "master"))) {
      
      print(paste0("Table: '", table, "'!"))
      
      query <- paste0("SELECT * from '", table,"'")
      db <- dbGetQuery(con, query) %>% as_tibble()
      
      annotated_sum_count <- c(annotated_sum_count, db$ref_sum_counts)
      novel_sum_count <- c(novel_sum_count, db$novel_sum_counts)
      
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
  annotated_sum_count %>% length()
  novel_donor_sum_count %>% length()
  novel_acceptor_sum_count %>% length()
  
  (annotated_sum_count %>% sum()) / (annotated_sum_count %>% length() + 
                                       novel_donor_sum_count %>% length() +
                                       novel_acceptor_sum_count %>% length())
  (novel_donor_sum_count %>% sum()) / (annotated_sum_count %>% length() + 
                                       novel_donor_sum_count %>% length() +
                                       novel_acceptor_sum_count %>% length())
  (novel_acceptor_sum_count %>% sum()) / (annotated_sum_count %>% length() + 
                                       novel_donor_sum_count %>% length() +
                                       novel_acceptor_sum_count %>% length())
  novel_sum_count %>% mean
  novel_donor_sum_count %>% mean 
  novel_acceptor_sum_count %>% mean 
  novel_donor_sum_count %>% sum() 
  novel_acceptor_sum_count %>% sum() 
  (novel_donor_sum_count %>% mean +
      novel_acceptor_sum_count %>% mean) /2
  
}

get_maxentscan_score <- function() {
  
  df_mes <- map_df(all_projects, function(project_id) {
    
    #############################
    ## GET THE CLUSTERS
    #############################
    
    print(paste0(Sys.time(), " - ", project_id))
    
    # clusters_used <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
    #                                        project_id, "/raw_data/all_clusters_used.rds"))
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    # if (!identical(clusters_used, all_clusters)) {
    #   print("Error! Both set of clusters are not identical.")
    # }
    
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
      
      #folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id)
      #samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster,  "_samples.rds"))
      #
      
      ####################
      ## GET THE INTRONS
      ####################
      
      query <- paste0("SELECT ref_junID, ref_type FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT DISTINCT ref_junID, ref_type FROM '", cluster_id, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      query <- paste0("SELECT ref_junID, ref_ss5score, ref_ss3score FROM 'intron' WHERE ref_junID IN (",
                      paste(introns$ref_junID, collapse = ","),")")
      introns <- merge(x = introns,
                       y = dbGetQuery(con, query) %>% as_tibble(),
                       by = "ref_junID",
                       all.x = T) %>% as_tibble() 
      
      
      # introns <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
      #                                  cluster, "/v105/", cluster, "_db_introns.rds"))
      
      ###########################
      ## GET THE NOVEL JUNCTIONS
      ###########################
      
      query <- paste0("SELECT ref_junID, novel_junID FROM '", cluster_id, "_", project_id, "_misspliced'")
      novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
      query <- paste0("SELECT novel_junID, novel_type, novel_ss5score, novel_ss3score FROM 'novel' WHERE novel_junID IN (", 
                      paste(novel_junctions$novel_junID, collapse = ","), ")")
      novel_junctions <- merge(x = novel_junctions,
                               y = dbGetQuery(con, query) %>% as_tibble(),
                               by = "novel_junID",
                               all.x = T) %>% as_tibble() 
      
      # folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id)
      # samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster,  "_samples.rds"))
      # 
      # db_introns <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
      #                                    cluster, "/v105/", cluster, "_db_introns.rds"))
      # db_novel <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
      #                                   cluster, "/v105/", cluster, "_db_novel.rds"))
      
      df_merged <- merge(x = introns %>% 
                           dplyr::select(ref_junID, ref_ss5score,
                                         ref_ss3score, ref_type),
                         y = novel_junctions %>% 
                           dplyr::select(ref_junID, novel_junID, 
                                         novel_ss5score, novel_ss3score, novel_type),
                         by = "ref_junID")  %>%
        mutate(diff_ss5score = ref_ss5score - novel_ss5score,
               diff_ss3score = ref_ss3score - novel_ss3score,
               tissue = cluster_id)
      
      df_merged %>% return()
    })
  })
  
  # all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_clusters_used.rds"))
  # 
  # df_mes <- df_mes %>%
  #   filter(tissue %in% all_clusters)
  
  saveRDS(object = df_mes %>% as_tibble(), 
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/paper/df_mes.rds")
  
  
  # df_mes <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_mes.rds")
  
  
  #############################
  ## COMBINED FIGURES 
  ## MAXENTSCAN score & DISTANCES
  #############################
  
  
  ## ss5score ----------------------------------------------------------
  
  
  df_5ss <- df_mes %>% 
    filter(novel_type == "novel_donor") %>%
    dplyr::select(intron = ref_ss5score, novel_donor = novel_ss5score) %>%
    gather(key = "junction_type", value = "ss5score")
  
  ss5plot <- ggplot(df_5ss, aes(ss5score, 
                                #group = fct_reorder(junction_type, ss5score, .fun = median, .desc = TRUE), 
                                fill = junction_type)) +
    geom_density(alpha = 0.8) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    scale_fill_manual(values = c("#A41F9AFF", "#EF7F4FFF"), 
                      breaks=c("intron", "novel_donor"),
                      labels=c("intron  ", "novel donor")) +
    theme(axis.line = element_line(colour = "black" ),
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          legend.text = element_text(size = "14"),
          legend.title = element_text(size = "14"),
          legend.position = "top") +
    xlab("5' splice site MES score") +
    guides(fill = guide_legend(title = element_blank(),
                               ncol = 2, nrow = 1))
  
  ss5plot
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel2_ss5plot.svg"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
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
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          legend.text = element_text(size = "14"),
          legend.title = element_text(size = "14"),
          legend.position = "top") +
    xlab("3' splice site MES score") +
    guides(fill = guide_legend(title = element_blank(), ncol = 2, nrow = 1))
  
  
  ss3plot
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel2_ss3plot.svg"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
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
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          legend.text = element_text(size = "14"),
          legend.title = element_text(size = "14"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  deltaplot
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel2_deltaplot.svg"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
}


##################################
## MODULO, DISTANCE Y MSR
##################################

## SECTION 2 ---------------------------------------------

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
    ggtitle("All transcripts") +
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
                               ncol = 3 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.text.x = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          strip.text = element_text(colour = "black", size = "14"), 
          legend.text = element_text(colour = "black", size = "14"),
          plot.caption = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "14"),
          legend.title = element_text(colour = "black", size = "14"),
          legend.position = "top") 
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
    geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
    theme_void()


  plot_all / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel2_distances.svg"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)

  
  
}

get_distances_PC <- function() {
  
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
  
  
  plot_PC <- ggplot(data = df_novel_tidy %>% filter(type_PC == "protein coding (PC)")) + 
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
                               ncol = 3 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.text.x = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          strip.text = element_text(colour = "black", size = "14"), 
          legend.text = element_text(colour = "black", size = "14"),
          plot.caption = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "14"),
          legend.title = element_text(colour = "black", size = "14"),
          legend.position = "top") 
  
  
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
    geom_rect(aes(xmin = (limit_bp * -1), xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
    theme_void()
  
  
  plot_PC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel2_distancesPC.svg"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  df_novel_tidy %>% filter(type_PC == "non PC") %>% distinct(ref_junID) %>% nrow()
  
  plot_NPC <- ggplot(data = df_novel_tidy %>% filter(type_PC == "non PC")) + 
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
                               ncol = 3 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.text.x = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          strip.text = element_text(colour = "black", size = "14"), 
          legend.text = element_text(colour = "black", size = "14"),
          plot.caption = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "14"),
          legend.title = element_text(colour = "black", size = "14"),
          legend.position = "top") 
  
  
  plot_NPC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel2_distancesNPC.svg"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  #######################################
  ## COMBINE ALL PLOTS
  #######################################
  
  # plot <- ggpubr::ggarrange(ss5plot,
  #                           ss3plot,
  #                           deltaplot,
  #                           plot_all / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1)), 
  #                           plot_PC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1)), 
  #                           plot_NPC / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1)),
  #                           #common.legend = T,
  #                           labels = c("a", "b", "c", 
  #                                      "d", "f", "g"),
  #                           align = "v",
  #                           ncol = 2, 
  #                           nrow = 3)
  # 
  # 
  # 
  # 
  # 
  # plot
  # 
  # file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel_2.svg"
  # ggplot2::ggsave(file_name, width = 183, height = 143, units = "mm", dpi = 300)


}

##IGNORE
get_modulo_basic_single_tissue <- function() {
  
  
  con <- dbConnect(RSQLite::SQLite(), "./database/splicing.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
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
  query <- paste0("SELECT ref_junID, protein_coding FROM 'intron' WHERE ref_junID IN (",
                  paste(introns$ref_junID, collapse = ","),")")
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
    mutate(module = abs(distance) %% 3)
  
  df_novel_tidy <- df_novel_tidy %>% 
    group_by(module) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  
  plot_single_tissue <- ggplot(data = df_novel_tidy, 
         aes(x = factor(module), y = freq*100, fill = factor(module))) + 
    geom_bar(stat = "identity", position = "dodge") +
    #ggtitle(paste0(title)) +
    geom_text(aes(label=round(freq*100, digits = 1)), size = 3, 
              position=position_dodge(width=0.9), vjust=-0.2, colour = "red") +
    #scale_x_continuous(breaks = c(0, 1, 2)) +
    #scale_y_continuous(labels = scales::percent) +
    scale_fill_viridis_d() +
    ylab("% of novel junctions") +
    xlab("Modulo 3") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          #axis.text.x = element_text(angle = 50, 
          #                           vjust = 1,
          #                           hjust = 1),
          legend.text = element_text(colour = "black",size = "9"),
          strip.text = element_text(colour = "black", size = "9"), 
          strip.text.y = element_blank(),
          plot.caption = element_text(colour = "black",size = "9"),
          legend.title = element_text(colour = "black", size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo 3: ",
                               ncol = 3, 
                               nrow = 1))
  
  ### PLOT ALL TISSUES
  
  
  df_novel_tidy <- df_novel %>%
    filter(abs(distance) <= 100) %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " ")) %>%
    mutate(type_p = ifelse(distance < 0, paste0(novel_type," intron"), paste0(novel_type," exon"))) %>% 
    mutate(module = abs(distance) %% 3)
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  df_novel_tidy$type_p = factor(df_novel_tidy$type_p, 
                                levels = c("novel acceptor intron", 
                                           "novel acceptor exon",
                                           "novel donor intron", 
                                           "novel donor exon"))
  title <- paste0("Modulo 3 - ", cluster)
  
  df_novel_tidy <- df_novel_tidy %>%
    group_by(module, type_p) %>% 
    summarise(count = n()) %>% 
    mutate(perc = count/sum(count))
  
  plot_fctx <- ggplot(data = df_novel_tidy, 
                 aes(x = factor(type_p), y = perc*100, fill = factor(module))) + 
    geom_bar(stat = "identity", position = "dodge") +
    #ggtitle(paste0(title)) +
    ggplot2::geom_text(aes(label=round(perc*100, digits = 1)), 
              size = 2,
              position=position_dodge(width=0.9), vjust=-0.2, colour = "red") +
    #scale_x_continuous(breaks = c(0, 1, 2)) +
    #scale_y_continuous(labels = scales::percent) +
    scale_fill_viridis_d() +
    ylab("% of novel junctions") +
    xlab("novel splice site location") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
    
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          axis.text.x = element_text(angle = 50, 
                                     vjust = 1,
                                     hjust = 1),
          legend.text = element_text(colour = "black",size = "9"),
          strip.text = element_text(colour = "black", size = "9"), 
          strip.text.y = element_blank(),
          plot.caption = element_text(colour = "black",size = "9"),
          legend.title = element_text(colour = "black", size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo 3: ",
                               ncol = 3, 
                               nrow = 1))
  
  plot_fctx
  
  
  ggpubr::ggarrange(plot_single_tissue,            
                    plot_all_tissues,
                    plot_fctx, 
                    plot_all,
                    common.legend = T,
                    labels = c("a", "b", "c", "d"),
                    ncol = 2, 
                    nrow = 2,
                    #widths = c(1,1,2,2),
                    align = "v")
  plot_single_tissue
  plot_all_tissues
  plot_fctx
  plot_all
  legend("a","b","c","d")
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel2.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
}

## NOT IGNORE
get_modulo_basic_multiple_tissue <- function() {
  
  ############################
  ## CONNECT TO THE DATABASE
  ############################
  
  con <- dbConnect(RSQLite::SQLite(), "database/splicing.sqlite")
  dbListTables(con)
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
                      paste(introns$ref_junID, collapse = ","),")  AND MANE = 1")
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
      
      df_novel_tidy <- df_novel_tidy %>% 
        group_by(modulo) %>%
        summarise(n = n()) %>%
        mutate(freq = n / sum(n)) %>%
        mutate(tissue = cluster_id)
      
      return(df_novel_tidy)
    })
  })
  
  saveRDS(df_modulo_tissues, file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_basic_tissues_100bpfilter.rds")
  # df_modulo_tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_basic_tissues_100bpfilter.rds")
  
  #############
  
  df_modulo_tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_basic_tissues_100bpfilter.rds")
  df_modulo_tissues$modulo = factor(df_modulo_tissues$modulo, 
                                        levels = c( "0", "1", "2"))
  df_modulo_tissues <- df_modulo_tissues %>% 
    ungroup() %>%
    arrange(modulo, freq) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df_modulo_tissues$tissue), 
                               pattern = "Brain"), "red", "black")
  
  plot_all_tissues <- ggplot(data = df_modulo_tissues, 
         aes(x = factor(tissue), y = freq*100, fill = factor(modulo))) + 
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip()+
    #ggtitle(paste0(title)) +
    #geom_text(aes(label=round(freq*100, digits = 1)), 
    #          position=position_dodge(width=0.9), vjust=-0.2, colour = "red") +
    #scale_x_continuous(breaks = c(0, 1, 2)) +
    #scale_y_continuous(labels = scales::percent) +
    scale_fill_viridis_d() +
    ylab("% of novel junctions") +
    xlab("Tissue") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          axis.text.y = element_text(colour = colours,
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
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel3_mod3alltissues.svg"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  ## (df_modulo_tissues %>% filter(modulo == 2) %>% pull(freq) %>% mean + df_modulo_tissues %>% filter(modulo == 2) %>% pull(freq) %>% mean) * 100
}

## IGNORE
get_modulo_tissues_v1 <- function() {
  
  #############################
  ## Supplementary Figure 4.1 
  ## MODULO
  #############################
  
  if (!file.exists("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues.rds")) {
    
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
    
    df_modulo_tissues <- map_df(all_projects, function(project_id) {
      
      # project_id <- all_projects[1]
      
      print(paste0(Sys.time(), " - ", project_id))
      
      all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                             project_id, "/raw_data/all_clusters_used.rds"))
      
      map_df(all_clusters, function(cluster) {
        
        # cluster <- all_clusters[1]
        print(cluster)
        
        ## Load IDB
        db_novel <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                          "/results/pipeline3/missplicing-ratio/", cluster, "/v105/", cluster, "_db_novel.rds"))
        ## Load novel junction database for current tissue
        db_novel %>% head()
        
        
        df_novel_tidy <- db_novel %>%
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
                                 filter(modulo == 0, type_p == "novel donor exon") %>% nrow()) * 100 / (df_novel_tidy %>%
                                                                                                          filter(type_p == "novel donor exon") %>% 
                                                                                                          nrow())
        modulo0_donor_intron <- (df_novel_tidy %>%
                                   filter(modulo == 0, type_p == "novel donor intron") %>% nrow()) * 100 / (df_novel_tidy %>%
                                                                                                              filter(type_p == "novel donor intron") %>% 
                                                                                                              nrow())
        modulo0_aceptor_exon <- (df_novel_tidy %>%
                                   filter(modulo == 0, type_p == "novel acceptor exon") %>% nrow()) * 100 / (df_novel_tidy %>%
                                                                                                                   filter(type_p == "novel acceptor exon") %>% 
                                                                                                                   nrow())
        modulo0_aceptor_intron <- (df_novel_tidy %>%
                                   filter(modulo == 0, type_p == "novel acceptor intron") %>% nrow()) * 100 / (df_novel_tidy %>%
                                                                                                               filter(type_p == "novel acceptor intron") %>% 
                                                                                                               nrow())
                                                                                                         
        ## Number of modulo 1
        modulo1_donor_exon <- (df_novel_tidy %>%
                                 filter(modulo == 1, type_p == "novel donor exon") %>% nrow()) * 100 / (df_novel_tidy %>% 
                                                                                                              filter(type_p == "novel donor exon") %>% 
                                                                                                              nrow())
        modulo1_donor_intron <- (df_novel_tidy %>%
                                 filter(modulo == 1, type_p == "novel donor intron") %>% nrow()) * 100 / (df_novel_tidy %>% 
                                                                                                          filter(type_p == "novel donor intron") %>% 
                                                                                                          nrow())
        modulo1_aceptor_exon <- (df_novel_tidy %>%
                                   filter(modulo == 1, type_p == "novel acceptor exon") %>% nrow()) * 100 / (df_novel_tidy %>%
                                                                                                               filter(type_p == "novel acceptor exon") %>% 
                                                                                                               nrow())
        modulo1_aceptor_intron <- (df_novel_tidy %>%
                                   filter(modulo == 1, type_p == "novel acceptor intron") %>% nrow()) * 100 / (df_novel_tidy %>%
                                                                                                               filter(type_p == "novel acceptor intron") %>% 
                                                                                                               nrow())
        
        ## Number of modulo 2
        modulo2_donor_exon <- (df_novel_tidy %>%
                                 filter(modulo == 2, type_p == "novel donor exon") %>% nrow()) * 100 / (df_novel_tidy %>% 
                                                                                                          filter(type_p == "novel donor exon") %>% 
                                                                                                          nrow())
        modulo2_donor_intron <- (df_novel_tidy %>%
                                 filter(modulo == 2, type_p == "novel donor intron") %>% nrow()) * 100 / (df_novel_tidy %>% 
                                                                                                          filter(type_p == "novel donor intron") %>% 
                                                                                                          nrow())
        modulo2_aceptor_exon <- (df_novel_tidy %>%
                              filter(modulo == 2, type_p == "novel acceptor exon") %>% nrow()) * 100 / (df_novel_tidy %>% 
                                                                                                          filter(type_p == "novel acceptor exon") %>% 
                                                                                                         nrow())
        modulo2_aceptor_intron <- (df_novel_tidy %>%
                                   filter(modulo == 2, type_p == "novel acceptor intron") %>% nrow()) * 100 / (df_novel_tidy %>% 
                                                                                                               filter(type_p == "novel acceptor intron") %>% 
                                                                                                               nrow())
        
        return(data.frame(tissue = cluster,
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
    saveRDS(df_modulo_tissues, file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues.rds")
    
  } else {
    df_modulo_tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues.rds")  
  }
  
  
  
  df_modulo_tissues$novel_type = factor(df_modulo_tissues$novel_type, 
                                        levels = c( "donor exon", "donor intron", "aceptor exon", "aceptor intron"))
  
  df_modulo_tissues <- df_modulo_tissues %>% 
    ungroup() %>%
    arrange(type , desc(modulo)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df_modulo_tissues$tissue), 
                               pattern = "Brain"), "red", "black")

  # library(scales)
  # show_col(viridis_pal()(40))
  
  plot <- ggplot(data = df_modulo_tissues) +
    geom_bar(mapping = aes(x = tissue, 
                           y = modulo, 
                           fill = factor(type)), 
             stat = "identity", 
             position = position_dodge()) +
    xlab("Tissue") +
    ylab("Percentage of junctions (%)") +
    theme_light() +
    facet_grid(~novel_type) +
    # scale_fill_viridis_d(option = "cividis")  +
    scale_fill_manual(values = c("#35B779FF","#440154FF","#440154FF"),
                      breaks = c("modulo0","modulo1","modulo2"),
                      labels = c("modulo0","modulo1","modulo2")) +
    #scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "14"),
          strip.text = element_text(colour = "black", size = "14"),
          axis.text.x = element_blank(),
          legend.text = element_text(size = "14"),
          legend.title = element_text(size = "14"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1))
  
  plot
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/modulo_all_tissues.png"
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  return(plot)
  
  

}

## NOT IGNORE
get_modulo_tissues_v2 <- function() {
  
  #############################
  ## Supplementary Figure 4.1 
  ## MODULO
  #############################
  
  if (!file.exists("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues.rds")) {
    
    con <- dbConnect(RSQLite::SQLite(), "database/splicing.sqlite")
    dbListTables(con)
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
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues_100bpfilter.rds")
    
  } else {
    df_modulo_tissues <- readRDS(file = 
                                   "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues_100bpfilter.rds")  
  }
  
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
  
  colours <- ifelse(str_detect(string = as.factor(df_modulo_tissues$tissue), 
                               pattern = "Brain"), "red", "black")
  
  
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
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel3_mod3alltissues100bp.svg"
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
}

get_missplicing_ratio_PC <- function()  {
  
  
  ###############################
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  query <- paste0("SELECT ref_junID, MSR_D,  MSR_A, ref_type FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT ref_junID, MSR_D,  MSR_A, ref_type FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
  query <- paste0("SELECT ref_junID, protein_coding FROM 'intron' WHERE ref_junID IN (",
                  paste(introns$ref_junID, collapse = ","),")")
  introns <- merge(x = introns,
                   y = dbGetQuery(con, query) %>% as_tibble(),
                   by = "ref_junID",
                   all.x = T) %>% 
    as_tibble() 
  
 
  
  introns %>% head()
  introns %>% nrow()
  
  # folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
  #                       cluster_id, "/v105/")
  # 
  # 
  # if (file.exists(paste0(folder_root, "/", cluster_id, "_db_introns.rds"))) {
  #   df <- readRDS(file = paste0(folder_root, "/", cluster_id, "_db_introns.rds"))
  # } else {
  #   print(paste0(paste0(folder_root, "/", cluster_id, "_db_introns.rds"), " file doesn't exist."))
  # }
  
  
  
  ################################
  ## TIDY DATA FOR PROTEIN-CODING
  ################################
  
  df_tidy <- introns %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC")) 
  
  df_tidy$type_PC = factor(df_tidy$type_PC, 
                                    levels = c("PC", "non PC"))
  
  
  df_tidy <- df_tidy %>%
    dplyr::select(type_PC, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -type_PC)
  
  ###########################
  ## PLOT MIS-SPLICING RATIO 
  ###########################
  
  plotMSR_PC <- ggplot(data = df_tidy %>% filter(type_PC == "PC")) + 
    geom_density(aes(x = MSR, fill = MSR_type), alpha = 0.8) +
    ggforce::facet_zoom(xlim=c(0,0.05)) +
    ggtitle("Protein-coding") +
    xlab("Mis-splicing ratio") +
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
  plotMSR_NPC <- ggplot(data = df_tidy %>% filter(type_PC == "non PC")) + 
    geom_density(aes(x = MSR, fill = MSR_type), alpha = 0.8) +
    ggforce::facet_zoom(xlim=c(0,0.05)) +
    #facet_wrap(vars(type_PC)) +
    ggtitle("Non protein-coding") +
    xlab("Mis-splicing ratio") +
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
  
  ggpubr::ggarrange(plotMSR_PC,
                    plotMSR_NPC,
                    labels = c("c", "d"), common.legend = T)
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel3_MSR.svg"
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)

  
  wilcox.test(x = df_tidy %>% filter(type_PC == "PC", MSR_type == "MSR_D") %>% pull(MSR),
              y = df_tidy %>% filter(type_PC == "PC", MSR_type == "MSR_A") %>% pull(MSR),
              alternative = "less")
  wilcox.test(x = df_tidy %>% filter(type_PC == "non PC", MSR_type == "MSR_D") %>% pull(MSR),
              y = df_tidy %>% filter(type_PC == "non PC", MSR_type == "MSR_A") %>% pull(MSR),
              alternative = "less")
}


#################################
## LINEAR MODELS
#################################

## SECTION 3 ---------------------------------------------

get_lm_single_tissue <- function() {
  
  
  ###############################
  ## QUERY THE DATABASE
  ## GET DATA FOR FRONTAL CORTEX
  ###############################
  
  project_id <- "BRAIN"
  cluster_id <- "Brain - Frontal Cortex (BA9)"
  
  
  query <- paste0("SELECT * 
                  FROM '", cluster_id, "_", project_id, "_nevermisspliced'")
  introns <- dbGetQuery(con, query) %>% as_tibble()
  query <- paste0("SELECT * 
                  FROM '", cluster_id, "_", project_id, "_misspliced'")
  introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
  
  query <- paste0("SELECT ref_junID, ref_length, ref_ss5score, ref_ss3score, 
                  ref_cons5score, ref_cons3score, ref_CDTS5score, ref_CDTS3score, 
                  protein_coding
                  FROM 'intron' 
                  WHERE ref_junID IN (",
                  paste(introns$ref_junID, collapse = ","),")")
  introns <- merge(x = introns,
                   y = dbGetQuery(con, query) %>% as_tibble(),
                   by = "ref_junID",
                   all.x = T) %>% 
    as_tibble() 
  
  
  query <- paste0("SELECT id, n_transcripts, gene_width, gene_name 
                  FROM 'gene' WHERE id IN (",
                  paste(introns$gene_id, collapse = ","),")")
  introns <- merge(x = introns,
                   y = dbGetQuery(con, query) %>% as_tibble(),
                   by.x = "gene_id",
                   by.y = "id",
                   all.x = T) %>% 
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
  fit_donor %>% summary()
  
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
                             robust = T,
                             #inner_ci_level = .75,
                             #n.sd = 2,
                             pvals = TRUE,
                             legend.title = "",
                             #plot.distributions = TRUE,
                             ci_level = 0.95,
                             coefs = coef_names,
                             colors = c("#35B779FF","#440154FF"),
                             model.names = model_names) + 
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
    geom_hline(yintercept = seq(from = 0.5,
                                to = length((fit_donor$coefficients %>% names)[-1]) + .5,
                                by = 1)) +
    guides(colour = guide_legend(#title = NULL,
                                 ncol = 2, 
                                 nrow = 1))
  
  plotLM
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel3_LM.svg"
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
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
  
  
  ################################
  ## CONNECT TO THE DATABASE
  ################################
  
  print(paste0(Sys.time(), " - getting unique and common junctions across clusters..."))
  
  
  ## Getting all introns that are common across all GTEx tissues -------------------------------------------
  
  all_introns <- list()
  
  for (project_id in all_projects) {
    
    print(paste0(Sys.time(), " - ", project_id))
    
    #############################
    ## GET THE CLUSTERS
    #############################
    
    print(paste0(Sys.time(), " - ", project_id))
    
    clusters_used <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    print(clusters_used)
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id,
             cluster %in% clusters_used) %>%
      distinct(cluster) %>%
      pull()
    
    if (!identical(clusters_used, all_clusters)) {
      print("Error! Both set of clusters are not identical.")
    }
    
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
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/common_introns.rds")
  saveRDS(object = introns, file = file_name)
  
}

get_estimate_variance_across_tissues <- function() {
  
  
  ## LOAD COMMON INTRONS
  common_introns <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/common_introns.rds")

  
  ################################
  ## QUERY THE DATABASE PER TISSUE
  ################################

  df_estimates <- map_df(all_projects, function(project_id) {
    
    print(paste0(Sys.time(), " - ", project_id))
    
    #############################
    ## GET THE CLUSTERS
    #############################
    
    print(paste0(Sys.time(), " - ", project_id))
    
    clusters_used <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    all_clusters <- df_metadata %>%
      filter(SRA_project == project_id,
             cluster %in% clusters_used) %>%
      distinct(cluster) %>%
      pull()
    
    if (!identical(clusters_used %>% sort(), all_clusters %>% sort())) {
      print("Error! Both set of clusters are not identical.")
    }
    
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
                      ref_CDTS3score AS CDTS_3ss, protein_coding, gene_id FROM 'intron' WHERE ref_junID IN (",
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
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/variance_estimate.rds")
  saveRDS(object = df_estimates, file = file_name)
  
}

plot_estimate_variance_across_tissues <- function() {
  

  df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/variance_estimate.rds")) %>%
    group_by(type) %>%
    mutate(pval_corrected =p.adjust(pval)) %>%
    filter(pval_corrected <= 0.05) %>%
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
    ylab("Distribution of the significant beta values (pval<=0.05)") +
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
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel3_tissuesLM.svg")
  ggplot2::ggsave(filename = file_name, width = 183, height = 143, units = "mm", dpi = 300)
  
  # if (save_results) {
  #   plotTissuesLM
  #   
  #   ## SAVE THE RESULTS
  #   if (brain_tissues) {
  #     file_name <- "/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_brain_tissues.png"
  #   } else {
  #     file_name <- "/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_all_tissues.png"
  #   }
  #   
  #   ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  #   
  # } else {
  #   plot %>% return()
  # }
}



#################################
## AGE STRATIFICATION
#################################

## SECTION 4 ---------------------------------------------

get_common_age_stratification <- function(project_id = "BRAIN") {
  
  # project_id <- "BRAIN"
  
  source(paste0("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline4-2_age_stratification.R"))
  
  #ref <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  
  age_samples_clusters_tidy <- age_stratification_init_data(project_id = project_id)
  age_supergroups <- age_samples_clusters_tidy$age_group %>% unique()
  
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                        "/results/pipeline3/missplicing-ratio/age/")
  
  ## Load the IDBs (the intron and the novel tables)
  df_age_groups_intron <- map_df(age_supergroups, function(age_group) {
    readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_introns.rds")) %>%
      mutate(sample_type = age_group) %>%
      return()
  })
  df_age_groups_novel <- map_df(age_supergroups, function(age_group) {
    readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_novel.rds")) %>%
      mutate(sample_type = age_group) %>%
      return()
  })
  
  ## Get the introns that are commonly mis-spliced across age supergroups
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
  

  ## Filter the INTRONS table by the common mis-spliced introns
  df_age_groups_intron_tidy <- merge(x = df_age_groups_intron %>% data.table::as.data.table(),
                                     y = data.table::data.table(ref_junID = common_introns),
                                     by = "ref_junID",
                                     all.y = T)
  
  if ((df_age_groups_intron_tidy$ref_junID %>% length()) / 
      (age_supergroups %>% length())) {
    print("CORRECT!")
  } else {
    print("ERROR!")
  }

  ## Filter the NOVEL table by the common mis-spliced introns
  df_age_groups_novel_tidy <- merge(x = df_age_groups_novel %>% data.table::as.data.table(),
                                    y = data.table::data.table(ref_junID = common_introns_misspliced),
                                    by = "ref_junID",
                                    all.y = T)
  if ((df_age_groups_novel_tidy$ref_junID %>% length()) / 
      (age_supergroups %>% length())) {
    print("CORRECT!")
  } else {
    print("ERROR!")
  }
  
  
  bg_genes <- df_age_groups_novel_tidy %>% distinct(gene_name) %>% pull() %>% unlist %>% unique()
  #bg_genes %>% length()
  #df_age_groups_intron_tidy%>% distinct(gene_name) %>% pull() %>% unlist %>% unique() %>% length()
  
  
  ##########################################################
  ## DISTANCES
  ##########################################################
  
  distance_plot <- age_stratification_plot_distances(df = df_age_groups_novel_tidy, 
                                                     age_levels = c("60-79","40-59","20-39" ), 
                                                     distance_limit = 30)
  
  ##########################################################
  ## MSR PLOTS
  ##########################################################
  
  MSR_plot <- age_stratification_plot_MSR(df = df_age_groups_intron_tidy,
                                          age_levels = c("20-39","40-59","60-79" )) 
  
  
  
  ##########################################################
  ## GO ENRICHMENT
  ##########################################################
  
  ## MSR_D -------------------------------------------------

  df_MSRD <- df_age_groups_tidy %>%
    dplyr::select(ref_junID,
           sample_type,
           MSR_D,
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
    arrange(fold_increase)
  
  genes_MSRD_discard <- df_MSRD %>%
    filter(`20-39` > `40-59` |
             `20-39` > `60-79`) %>%
    distinct(gene_name) %>% 
    pull() %>% unlist()
  
  genes_incrD <- setdiff(genes_MSRD_increasing$gene_name %>% unlist() %>% unique, genes_MSRD_discard)
  saveRDS(object = genes_incrD,
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRD.rds")
  
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
  
  
  enrichplot::treeplot(edox2) +
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
    #filter(ref_junID %in% df_age_groups_novel_tidy$ref_junID) %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_A = ref_missplicing_ratio_tissue_NA,
                  gene_name) %>%
    mutate(MSR_A = MSR_A %>% round(digits = 4)) %>%
    spread(sample_type, MSR_A)

  genes_MSRA_discard <- df_MSRA %>%
    filter(`20-39` > `40-59` |
             `20-39` > `60-79`) %>%
    distinct(gene_name) %>% 
    pull() %>% unlist()
  
  genes_MSRA_increasing <- df_MSRA %>%
    filter(`20-39` < `40-59`,
           `40-59` < `60-79`) %>%
    mutate_if(is.list, simplify_all) %>%
    mutate(`20-39` = replace(`20-39`, `20-39` == 0, 0.0000001)) %>%
    mutate(`60-79` = replace(`60-79`, `60-79` == 0, 0.0000001)) %>%
    mutate(fold_increase = gtools::foldchange(`20-39`, `60-79`)) %>%
    mutate(MSR_type = "MSR_Acceptor") %>%
    dplyr::select(gene_name, fold_increase, MSR_type) %>%
    arrange(fold_increase)
  
  ########################################################
  
  genes_incrA <- setdiff(genes_MSRA_increasing$gene_name %>% unlist() %>% unique, genes_MSRA_discard)
  saveRDS(object = genes_incrA,
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRA.rds")
  
 
  
  ########################################################

  query_genes <- genes_MSRA_increasing$gene_name %>% unlist() %>% unique()
  #query_genes <- setdiff(genes_MSRA_increasing$gene_name %>% unlist() %>% unique(), genes_discard)
  
  ego_A <- clusterProfiler::enrichGO(
                                    gene          = query_genes, 
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
    guides(fill = guide_legend(title = "pval",
                               label.position = "bottom") 
           ) 

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
    
    df_analysis_lm_uncorrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
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

supplementary_age_stratification <- function() {
  
  # project_id <- "MUSCLE"
  # project_id <- "BLOOD"
  
  source(paste0("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline4-2_age_stratification.R"))
  
  ref <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  
  age_samples_clusters_tidy <- age_stratification_init_data(project_id)
  age_supergroups <- age_samples_clusters_tidy$age_group %>% unique()
  
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
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
    
    df_analysis_lm_uncorrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
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
## ENCORI STUFF
##################################






 





  