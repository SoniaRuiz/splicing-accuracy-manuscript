library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)

####################################################
## PAPER ###########################################
####################################################


get_contamination_rates <- function(all_tissues = T) {
  
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
      arrange(type , desc(prop)) %>%
      mutate(tissue = fct_inorder(tissue))
    
    colours <- ifelse(str_detect(string = as.factor(df_contamination_tidy$tissue), pattern = "Brain"), "red", "black")
    
    ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
    ggplot(data = df_contamination_tidy) +
      #geom_bar(data = novel_104, mapping = aes(x = novel_p_individuals, y = ..count..), stat = "count", color = "black") +
      geom_bar(mapping = aes(x = tissue, y = prop, fill = type), 
               stat = "identity", position = "identity") + 
      xlab(" ") +
      ylab("contamination rates (%)") +
      #ggtitle(paste0("Contamination rates in Ensembl v97 vs v105")) +
      #scale_y_continuous(breaks = c(0,0.3,0.6,0.9)) +
      #scale_x_binned() +
      #scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100))+
      scale_fill_manual(values = c("#440154FF","#FDE725FF"),
                        breaks = c("in_annotation", "out_annotation"),
                        labels = c("introns entring annotation", "introns exiting annotation")) +
      theme_light() +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12"),
            legend.text = element_text(size = "12"),
            legend.title = element_text(size = "12"),
            legend.position = "top",
            axis.text.x = element_text(color = colours,
                                       angle = 70, 
                                       vjust = 1,
                                       hjust = 1)) +
      guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) %>%
      return()
    
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/contamination_rates_tissues.png")
    ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
    
    ## Stats
    df_contamination_tidy %>% filter(type == "in_annotation") %>% pull(prop) %>% mean()
    
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
                                                  levels = c( "Ensembl_v76 (18-Jul-2014)", 
                                                              "Ensembl_v81 (07-Jul-2015)", 
                                                              "Ensembl_v90 (28-Jul-2017)", 
                                                              "Ensembl_v97 (26-May-2019)", 
                                                              "Ensembl_v104 (26-May-2019)"))
    
    ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
    ggplot(data = df_contamination) +
      geom_bar(aes(x = contamination_rates, y = in_annotation), 
               stat = "identity") + 
      xlab(NULL) +
      ylab("contamination rates (%)") +
      #ggtitle(paste0("Brain - Frontal Cortex (BA9)\nContamination rates compared with Ensembl v105")) +
      theme_light() +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12"),
            legend.text = element_text(size = "12"),
            legend.title = element_text(size = "12"),
            legend.position = NULL,
            axis.text.x = element_text(angle = 50, 
                                       vjust = 1,
                                       hjust = 1)) +
      guides(fill = "none")
    
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/contamination_rates_FCTX.png")
    ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
    
  }
  
  
  
  
}

get_unique_donor_acceptor_jxn <- function() {
  
  ##################################################################################################
  ## Proportion of each type of unique junctions per GTEx tissue.
  ##################################################################################################
  
  
  all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  
  
  df_proportions <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
  
    
    map_df(all_clusters, function(cluster) {
    
      # cluster <- all_clusters[1]
      
      ## Print the tissue
      print(paste0(Sys.time(), " - ", cluster))
      
      folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id)
      samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster,  "_samples.rds"))
      
      
      
      introns <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
                                       cluster, "/v105/", cluster, "_db_introns.rds"))
      novel_junctions <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
                                               cluster, "/v105/", cluster, "_db_novel.rds"))
      
      
      introns %>% head()
      novel_junctions %>% head()
      
      
      
      ## Calculate the proportions
      annotated_junc <- (introns %>% distinct(ref_junID) %>% nrow()) / (samples %>% length())
      donor_junc <- (novel_junctions %>% filter(novel_type == "novel_donor") %>% distinct(novel_junID) %>% nrow()) / (samples %>% length())
      acceptor_junc <- (novel_junctions %>% filter(novel_type == "novel_acceptor") %>% distinct(novel_junID) %>% nrow()) / (samples %>% length())
      
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
      return(data.frame(tissue = cluster,
                        annotated_junc = annotated_junc ,
                        donor_junc = donor_junc,
                        acceptor_junc = acceptor_junc,
                        annotated_prop = annotated_prop,
                        donor_prop = donor_prop,
                        acceptor_prop = acceptor_prop,
                        samples = samples %>% length()))
    })
  })
  
  ## Save results --------------------------------------------------------------------------

  if (exists("df_proportions"))  {
    saveRDS(object = df_proportions,
            file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/unique_donor_acceptor_jxn.rds")
  } else {
    df_proportions <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/unique_donor_acceptor_jxn.rds")
  }
  

  df_proportions %>% head()
  
  df_prop_tidy <- df_proportions %>% 
    dplyr::select(tissue, donor = donor_junc, acceptor = acceptor_junc) %>%
    tidyr::gather(key = "type", value = "prop", -tissue ) 
  
  
  ## Stats
  df_proportions %>%
    mutate(annotated_raw = annotated_junc * samples) %>% pull(annotated_raw) %>% mean
    arrange(desc(annotated_raw))
  df_proportions %>%
    mutate(donor_raw = donor_junc * samples) %>% pull(donor_raw) %>% mean
    arrange(desc(donor_raw))
  df_proportions %>%
    mutate(acceptor_raw = acceptor_junc * samples) %>% pull(acceptor_raw) %>% mean
    arrange(desc(acceptor_raw))
    
    
    
  ##############################
  ## Plot the proportions 
  ##############################
  
    
  df_proportions %>% head()
  df_prop_tidy %>% head()
  
  
  ## First, order by junction type
  df_prop_tidy$type = factor(df_prop_tidy$type, 
                             levels = c( "acceptor", "donor"))
  df2 = df_prop_tidy %>% 
    ungroup() %>%
    arrange(type , desc(prop)) %>%
    mutate(tissue = fct_inorder(tissue))
  colours <- ifelse(str_detect(string = as.factor(df2$tissue), pattern = "Brain"), "red", "black")
  
  
  ggplot(df2, aes(x = tissue, 
                  y = prop, 
                  group = type, 
                  fill = type)) +
    geom_col(position = position_dodge()) +
    #geom_text(aes(label = round(x = prop, digits = 2)), 
    #          position = position_stack(vjust = 0.5, reverse = TRUE), colour = "white", size = 1.5) +
    #coord_flip() +
    theme_light() +
    ylab("unique novel junctions") +
    xlab("") +
    #scale_fill_viridis_d()  +
    #ggtitle("Percentage of unique introns and unique novel junctions per tissue") +
    #scale_x_discrete("tissue", breaks = breaks) +
    
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          #axis.text.y = element_text(color = colours),
          axis.text.x = element_text(color = colours,
                                     angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          plot.title = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    
    #scale_fill_manual(values =  c("#0D0887FF","#EF7F4FFF","#A41F9AFF"
    scale_fill_manual(values = c( "#0D0887FF", "#EF7F4FFF"),
                      breaks = c( "acceptor", "donor" ),
                      labels = c( "novel acceptor", "novel donor")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 4, 
                               nrow = 1))
  
  
  ## Save the figure 3
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/percent_unique_donor_acceptor.png"
  ggplot2::ggsave(filename = file_name, 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  # ############
  # 
  # gtex_tpm <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_gene_median_tpm_tidy.rds")
  # acceptor_genes <- c("ENSG00000161547") #"ENSG00000115524", "ENSG00000161547", "ENSG00000160201","ENSG00000169249"
  # df_acceptor_tpm <- gtex_tpm %>%
  #   filter(gene_id %in% acceptor_genes) %>%
  #   gather(key = tissue, value = tpm, -gene_id)  
  # df_acceptor_tpm <- cbind(df_acceptor_tpm,
  #                          data.frame(tissue_name = tissues_tidy))
  # 
  # df_acceptor_tpm
  # 
  # df <- df_prop_tidy %>% 
  #   filter(type == "acceptor_prop")
  # 
  # df_merged <- merge(df,
  #                    df_acceptor_tpm,
  #                    by = "tissue_name")
  # 
  # cor.test(x = df_merged$prop,
  #          y = df_merged$tpm)
}

get_unique_donor_acceptor_reads <- function() {
  
  ##############################################################################################
  ## Distribution of mean counts of novel donor, novel acceptor and annotated junctions 
  ## across all tissues
  ##############################################################################################
  
  
  all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  
  
  df_mean_counts <- map_df(all_projects, function(project_id) {
    
    # project_id <- all_projects[1]
    
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    map_df(all_clusters, function(cluster) {
    
      # tissue <- tissues[11]
      print(cluster)
      
      ## Load IDB
      
      folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id)
      samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster,  "_samples.rds"))
      
      db_intron <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
                                       cluster, "/v105/", cluster, "_db_introns.rds"))
      db_novel <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
                                               cluster, "/v105/", cluster, "_db_novel.rds"))
      
      
      annotated <- db_intron %>%
        dplyr::distinct(ref_junID, .keep_all = T) %>%
        dplyr::mutate(ref_mean_counts = ref_sum_counts / ref_n_individuals) %>%
        pull(ref_mean_counts) %>% 
        mean()
      
      acceptor <- db_novel %>%
        filter(novel_type == "novel_acceptor") %>%
        dplyr::distinct(novel_junID, .keep_all = T) %>%
        mutate(novel_mean_counts = novel_sum_counts / novel_n_individuals ) %>%
        pull(novel_mean_counts) %>% 
        mean()
      
      donor <- db_novel %>%
        filter(novel_type == "novel_donor") %>%
        dplyr::distinct(novel_junID, .keep_all = T) %>%
        mutate(novel_mean_counts = novel_sum_counts / novel_n_individuals ) %>%
        pull(novel_mean_counts) %>% 
        mean()

      
      return(data.frame(tissue = cluster,
                        samples = samples %>% length(),
                        type = c("annotated","acceptor", "donor"),
                        counts = c(annotated, acceptor, donor)))
    
    })
    
  })
  
  # ## Save data
  
  if (exists("df_mean_counts"))  {
    saveRDS(object = df_mean_counts,
            file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/unique_donor_acceptor_reads.rds")
  } else {
    df_mean_counts <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/unique_donor_acceptor_reads.rds")
  }
  
  # saveRDS(df_mean_counts, file = "/home/sruiz/PROJECTS/splicing-project/results/paper/figure1_data.rds")
  
  # df_mean_counts <- readRDS("/home/sruiz/PROJECTS/splicing-project/results/paper/figure1_data.rds")
  
  
  ## Spread
  df_mean_counts <- df_mean_counts %>% 
    tidyr::spread(key = "type", value = "counts")

  
  ## Convert to percentage
  df_mean_counts <- df_mean_counts %>% 
    mutate(annotated_p = annotated * 100 / (annotated + acceptor + donor),
           acceptor_p = acceptor * 100 / (annotated + acceptor + donor),
           donor_p = donor * 100 / (annotated + acceptor + donor))
  
  ## Gather
  df_mean_counts_tidy <- df_mean_counts %>% 
    dplyr::select(annotated = annotated_p, acceptor = acceptor_p, donor = donor_p, tissue) %>%
    tidyr::gather(key = "type", value = "mean", -tissue )
  
  
  df_mean_counts_tidy$type = factor(df_mean_counts_tidy$type, 
                                    levels = c("annotated", "acceptor", "donor"))
  
  df_mean_counts_tidy = df_mean_counts_tidy %>% 
    ungroup() %>%
    arrange(type , mean) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df_mean_counts_tidy$tissue), pattern = "Brain"), "red", "black")
  
  
  ggplot(data = df_mean_counts_tidy) +
    geom_col(mapping = aes(x = tissue, y = mean, fill = type), 
             #position = position_dodge(width = 0.6), 
             alpha = 0.9) +
    xlab("") +
    ylab("percent of avg. reads across samples (%)") +
    theme_light() +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = c(0,10,25,50,75,100)) +
    #coord_cartesian(ylim=c(1.1,50)) +
    scale_fill_manual(values =  c("#0D0887FF","#EF7F4FFF","#A41F9AFF"),#c( "#440154FF", "#35B779FF", "#39558CFF"),
                      breaks = c("annotated", "acceptor", "donor"),
                      labels = c("annotated", "novel acceptor", "novel donor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          axis.text.x = element_text(color = colours,
                                     angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) #+
  #facet_zoom(ylim=c(0,10))
  
  
  # library(scales)
  # show_col(viridis_pal()(20)) 
  # viridis_pal()(20)
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/unique_donor_acceptor_reads.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  write.csv(x = df_mean_counts,
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/unique_donor_acceptor_reads.csv", row.names = FALSE)
  
  
  # ############
  # 
  # gtex_tpm <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_gene_median_tpm_tidy.rds")
  # NMD_genes <- c("ENSG00000005007", "ENSG00000151461", "ENSG00000169062", "ENSG00000125351")
  # acceptor_genes <- c("ENSG00000169249") #"ENSG00000115524", "ENSG00000161547", "ENSG00000160201","ENSG00000169249"
  # df_acceptor_tpm <- gtex_tpm %>%
  #   filter(gene_id %in% acceptor_genes) %>%
  #   gather(key = tissue, value = tpm, -gene_id)  
  # df_acceptor_tpm <- cbind(df_acceptor_tpm,
  #                          data.frame(tissue_name = tissues_tidy))
  # 
  # df_acceptor_tpm
  # 
  # df <- df_mean_counts_tidy %>% 
  #   filter(type == "acceptor")
  # 
  # df_merged <- merge(df,
  #                    df_acceptor_tpm,
  #                    by = "tissue_name")
  # 
  # df_merged  %>% head
  # 
  # cor.test(x = df_merged$mean,
  #          y = df_merged$tpm)
  # 
  # ## Donor
  # df <- df_mean_counts_tidy %>% 
  #   filter(type == "donor")
  # 
  # df_merged <- merge(df,
  #                    df_acceptor_tpm,
  #                    by = "tissue_name")
  # 
  # df_merged  %>% head
  # cor.test(x = df_merged$mean,
  #          y = df_merged$tpm)
  # 
  # 
  # ## Annotated
  # df <- df_mean_counts_tidy %>% 
  #   filter(type == "annotated")
  # 
  # df_merged <- merge(df,
  #                    df_acceptor_tpm,
  #                    by = "tissue_name")
  # 
  # df_merged  %>% head
  # cor.test(x = df_merged$mean,
  #          y = df_merged$tpm)
  
  
}

get_maxentscan_score <- function() {
  
  #############################
  ## Figure 5 ##
  ## MAXENTSCAN score
  #############################
  
  all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  
  
  df_mes <- map_df(all_projects, function(project_id) {
    
    
    print(paste0(Sys.time(), " - ", project_id))
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    map_df(all_clusters, function(cluster) {
      
      # tissue <- tissues[11]
      print(cluster)
      
      ## Load IDB
      
      folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id)
      samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster,  "_samples.rds"))
      
      db_introns <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
                                         cluster, "/v105/", cluster, "_db_introns.rds"))
      db_novel <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
                                        cluster, "/v105/", cluster, "_db_novel.rds"))
      
      df_merged <- merge(x = db_introns %>% dplyr::select(ref_junID, 
                                                          ref_ss5score, ref_ss3score, ref_type),
                         y = db_novel %>% dplyr::select(ref_junID, 
                                                        novel_junID, 
                                                        novel_ss5score, novel_ss3score, novel_type),
                         by = "ref_junID")  %>%
        mutate(diff_ss5score = ref_ss5score - novel_ss5score,
               diff_ss3score = ref_ss3score - novel_ss3score,
               tissue = cluster)
      
      return(df_merged)
    })
  })
  
  saveRDS(object = df_mes, 
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_mes.rds")
  
  
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
    scale_fill_manual(values = c("#0D0887FF", "#A41F9AFF"), 
                      breaks=c("intron", "novel_donor"),
                      labels=c("intron  ", "novel donor")) +
    theme(axis.line = element_line(colour = "black" ),
          axis.text = element_text(colour = "black", size = "8"),
          axis.title = element_text(colour = "black", size = "8"),
          legend.text = element_text(size = "8"),
          legend.title = element_text(size = "8"),
          legend.position = "top") +
    xlab("5' splice site MES score") +
    #ggtitle("Distributions of the 5'ss ME scores obtained from\n7484 samples and 40 GTEx tissues") +
    guides(fill = guide_legend(title = element_blank(),
                               ncol = 2, nrow = 1))# + 
    #facet_grid(~tissue)
  
  
  
  ## ss3score -----------------------------------------------------------
  
  
  df_3ss <- df_mes %>% 
    filter(novel_type == "novel_acceptor") %>%
    dplyr::select(intron = ref_ss3score, novel_acceptor = novel_ss3score) %>%
    gather(key = "junction_type", value = "ss3score")
  
  ss3plot <- ggplot(df_3ss, aes(ss3score, 
                                #group = fct_reorder(junction_type, ss3score, .desc = TRUE), 
                                fill = junction_type)) +
    geom_density(alpha = 0.8) +
    ylim(c(0, 0.3)) +
    xlim(c(-40, 20)) +
    theme_light() +
    scale_fill_manual(values = c("#0D0887FF", "#EF7F4FFF"), 
                      breaks = c("intron", "novel_acceptor"),
                      labels = c("intron  ", "novel acceptor")) +
    theme(axis.line = element_line(colour = "black" ),
          axis.text = element_text(colour = "black", size = "8"),
          axis.title = element_text(colour = "black", size = "8"),
          legend.text = element_text(size = "8"),
          legend.title = element_text(size = "8"),
          legend.position = "top") +
    xlab("3' splice site MES score") +
    guides(fill = guide_legend(title = element_blank(), ncol = 2, nrow = 1))  #+ 
    #facet_grid(~tissue)
  
  
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
    scale_fill_manual(values =  c("#35B779FF","#440154FF"),
                      breaks = c("diff_ss5score", "diff_ss3score"),
                      labels = c("Delta 5'ss  ", "Delta 3'ss")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "8"),
          axis.title = element_text(colour = "black", size = "8"),
          #axis.text.x = element_text(vjust = 1,
          #                           hjust = 1),
          legend.text = element_text(size = "8"),
          legend.title = element_text(size = "8"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  
  #############################
  ## DISTANCES
  #############################
  
  ## Combine plots -------------------------------------------------
  
  plot <- ggpubr::ggarrange(ss5plot, 
                            ss3plot,
                            deltaplot,
                            #common.legend = T,
                             
                            plot_all, 
                            plot_PC, 
                            plot_NPC,
                            #common.legend = T,
                            labels = c("a", "b", "c", "d", "e", "f"),
                            ncol = 2, 
                            nrow = 3,
                            align = "v")
  

  
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
}


get_distances <- function() {
  
  limit_bp <- 30
  cluster <- "Brain - Frontal Cortex (BA9)"
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  
  if (file.exists(paste0(folder_root, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_root, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_root, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
  }
  
  #############################
  ## PLOT THE DISTANCES GRAPH 
  #############################
  
  
  df_novel_tidy <- df_novel %>%
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
    ggtitle("All transcripts biotype") +
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
          axis.text = element_text(colour = "black", size = "8"),
          axis.text.x = element_text(colour = "black", size = "8"),
          axis.title = element_text(colour = "black", size = "8"),
          strip.text = element_text(colour = "black", size = "7"), 
          legend.text = element_text(colour = "black", size = "8"),
          plot.caption = element_text(colour = "black", size = "8"),
          plot.title = element_text(colour = "black", size = "8"),
          legend.title = element_text(colour = "black", size = "8"),
          legend.position = "top") 
  

  plot_all
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/distances.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)

  
  
}

get_distances_PC <- function() {
  
  limit_bp <- 30
  
  cluster <- "Brain - Frontal Cortex (BA9)"
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  
  if (file.exists(paste0(folder_root, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_root, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_root, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
  }
  
  
  
  
  #############################
  ## PLOT THE DISTANCES GRAPH 
  #############################

  
  df_novel_tidy <- df_novel %>%
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
          axis.text = element_text(colour = "black", size = "8"),
          axis.text.x = element_text(colour = "black", size = "8"),
          axis.title = element_text(colour = "black", size = "8"),
          strip.text = element_text(colour = "black", size = "7"), 
          legend.text = element_text(colour = "black", size = "8"),
          plot.caption = element_text(colour = "black", size = "8"),
          plot.title = element_text(colour = "black", size = "8"),
          legend.title = element_text(colour = "black", size = "8"),
          legend.position = "top") 
  
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
          axis.text = element_text(colour = "black", size = "8"),
          axis.text.x = element_text(colour = "black", size = "8"),
          axis.title = element_text(colour = "black", size = "8"),
          strip.text = element_text(colour = "black", size = "7"), 
          legend.text = element_text(colour = "black", size = "8"),
          plot.caption = element_text(colour = "black", size = "8"),
          plot.title = element_text(colour = "black", size = "8"),
          legend.title = element_text(colour = "black", size = "8"),
          legend.position = "top") 
  
  
  plot <- ggpubr::ggarrange(plot_all, 
                            plot_PC, 
                            plot_NPC,
                            common.legend = T,
                            labels = c("a", "b", "c"),
                            align = "v",
                            ncol = 2, 
                            nrow = 2)
  

  plot
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/distances_all.png"
  ggplot2::ggsave(file_name, width = 183, height = 143, units = "mm", dpi = 300)

  
  
  ## RELEASE SOME MEMORY
  rm(df_tidy)
  rm(title)
  rm(file_name)
  rm(y_axes)
}

get_modulo_basic_single_tissue <- function() {
  
  cluster <- "Brain - Frontal Cortex (BA9)"
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  
  if (file.exists(paste0(folder_root, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_root, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_root, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
  }
  ## Load novel junction database for current tissue
  df_novel %>% head()
  
  
  df_novel_tidy <- df_novel %>%
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

get_modulo_basic_multiple_tissue <- function() {
  
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
      db_intron <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                         "/results/pipeline3/missplicing-ratio/", cluster, "/v105/", cluster, "_db_introns.rds"))
      db_intron <- db_intron %>%
        filter(MANE == T)
      db_novel <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                        "/results/pipeline3/missplicing-ratio/", cluster, "/v105/", cluster, "_db_novel.rds"))
      ## Load novel junction database for current tissue
      db_novel %>% head()
      
      
      df_novel_tidy <- db_novel %>%
        filter(abs(distance) <= 100,
               ref_junID %in% db_intron$ref_junID) %>%
        mutate(novel_type = str_replace(string = novel_type,
                                        pattern = "_",
                                        replacement = " ")) %>%
        mutate(type_p = ifelse(distance < 0, paste0(novel_type," intron"), paste0(novel_type," exon"))) %>% 
        mutate(modulo = abs(distance) %% 3)
      
      df_novel_tidy <- df_novel_tidy %>% 
        group_by(modulo) %>%
        summarise(n = n()) %>%
        mutate(freq = n / sum(n)) %>%
        mutate(tissue = cluster)
      
      return(df_novel_tidy)
    })
  })
  
  saveRDS(df_modulo_tissues, file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_basic_tissues_100bpfilter.rds")
  
  #############
  
  df_modulo_tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_basic_tissues_100bpfilter.rds")
  df_modulo_tissues$modulo = factor(df_modulo_tissues$modulo, 
                                        levels = c( "0", "1", "2"))
  df_modulo_tissues <- df_modulo_tissues %>% 
    ungroup() %>%
    arrange(modulo, desc(freq)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df_modulo_tissues$tissue), 
                               pattern = "Brain"), "red", "black")
  
  plot_all_tissues <- ggplot(data = df_modulo_tissues, 
         aes(x = factor(tissue), y = freq*100, fill = factor(modulo))) + 
    geom_bar(stat = "identity", position = "dodge") +
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
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          axis.text.x = element_blank(),#text(colour = colours,
                                     #angle = 75, 
                                     #vjust = 1,
                                     #hjust = 1),
          legend.text = element_text(colour = "black",size = "9"),
          strip.text = element_text(colour = "black", size = "9"), 
          strip.text.y = element_blank(),
          plot.caption = element_text(colour = "black",size = "9"),
          legend.title = element_text(colour = "black", size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo 3: ",
                               ncol = 3, 
                               nrow = 1))
    
  plot_all_tissues
  ## (df_modulo_tissues %>% filter(modulo == 2) %>% pull(freq) %>% mean + df_modulo_tissues %>% filter(modulo == 2) %>% pull(freq) %>% mean) * 100
}

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

get_modulo_tissues_v2 <- function() {
  
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
        db_intron <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                          "/results/pipeline3/missplicing-ratio/", cluster, "/v105/", cluster, "_db_introns.rds"))
        db_intron <- db_intron %>%
          filter(MANE == T)
        db_novel <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                          "/results/pipeline3/missplicing-ratio/", cluster, "/v105/", cluster, "_db_novel.rds"))
        ## Load novel junction database for current tissue
        db_novel %>% head()
        
        
        df_novel_tidy <- db_novel %>%
          filter(abs(distance) <= 100,
                 ref_junID %in% db_intron$ref_junID) %>%
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
    saveRDS(df_modulo_tissues, file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues_100bpfilter.rds")
    
  } else {
    df_modulo_tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/df_modulo_tissues.rds")  
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
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          plot.title = element_text(colour = "black", size = "9"),
          strip.text = element_text(colour = "black", size = "9"),
          axis.text.x = element_blank(), #text(color = colours, angle = 60, vjust = 1, hjust = 1,size = "11"),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Modulo 3:", ncol = 3, nrow = 1))
  
  plot_all
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/modulo_all_tissues100bp.png"
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  return(plot)
  
  
  
}

get_missplicing_ratio_PC <- function()  {
  
  
  cluster <- "Brain - Frontal Cortex (BA9)"
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  
  
  if (file.exists(paste0(folder_root, "/", cluster, "_db_introns.rds"))) {
    df <- readRDS(file = paste0(folder_root, "/", cluster, "_db_introns.rds"))
  } else {
    print(paste0(paste0(folder_root, "/", cluster, "_db_introns.rds"), " file doesn't exist."))
  }
  
  df %>% head()
  df %>% nrow()
  
  
  
  df_tidy <- df %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) 
  
  ##########################
  ## PLOT MIS-SPLICING RATIO 
  ##########################
  
  y_axes_values <- density(x = df_tidy %>% filter(protein_coding == 100) %>% pull(ref_missplicing_ratio_tissue_ND))$y
  y_axes <- c(0, y_axes_values %>% max() + 
                y_axes_values %>% max() / 4)
  
  plot_PC <- ggplot(data = df_tidy %>% filter(protein_coding == 100)) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"),
                 alpha = 0.8) +
    geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
                 alpha = 0.8) +
    ggtitle(paste0("Protein coding (PC)\n")) +
    xlab("Mis-splicing ratio") +
    #ylim(y_axes) +
    ggforce::facet_zoom(xlim = c(0,0.15)) +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("MSR_D","MSR_A")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          strip.text.x =  element_text(colour = "black", size = "9"),
          plot.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  plot_NPC <- ggplot(data = df_tidy %>% filter(protein_coding == 0)) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"),
                 alpha = 0.8) +
    geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
                 alpha = 0.8) +
    ggtitle(paste0("Non-protein coding (NPC)\n")) +
    xlab("Mis-splicing ratio") +
    ggforce::facet_zoom(xlim = c(0,0.15)) +
    #ylim(y_axes) +
    ylab("") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("MSR_D","MSR_A")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          plot.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  plotMSR <- ggpubr::ggarrange(plot_PC, 
                            plot_NPC,
                            common.legend = T,
                            #labels = c("a.1","a.2"),
                            align = "hv",
                            ncol = 2, 
                            nrow = 1)

  
  
  
  
  ggpubr::ggarrange(plotMSR,            
                    plotLM,
                    plotTissuesLM, 
                    common.legend = F,
                    labels = c("a", "b", "c"),
                    ncol = 2, 
                    nrow = 2,
                    #widths = c(1,1,2,2),
                    align = "v")
  
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel3.png")
  ggplot2::ggsave(filename = file_name, width = 183, height = 143, units = "mm", dpi = 300)
  
}

get_lm_single_tissue <- function() {
  
  #########################
  ## LOAD AND TIDY DATA
  #########################
  cluster <- "Brain - Frontal Cortex (BA9)"
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  

  idb <- readRDS(file = paste0(folder_root, "/", cluster, "_db_introns.rds"))
  
  # idb <- idb %>%
  #   as.data.frame() %>%
  #   dplyr::distinct(ref_junID, .keep_all = T) %>%
  #   dplyr::rename(intron_length = width,
  #                 intron_5ss_score = ref_ss5score,
  #                 intron_3ss_score = ref_ss3score,
  #                 gene_length = gene_width,
  #                 gene_tpm = tpm_median_rct3,
  #                 gene_num_transcripts = n_transcripts,
  #                 CDTS_5ss = CDTS_5ss_mean,
  #                 CDTS_3ss = CDTS_3ss_mean,
  #                 mean_phastCons20way_5ss = phastCons20way_5ss_mean,
  #                 mean_phastCons20way_3ss = phastCons20way_3ss_mean)
  idb <- idb %>%
    as.data.frame() %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    dplyr::rename(intron_length = width,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score,
                  gene_length = gene_width,
                  gene_tpm = tpm_median_rct3,
                  gene_num_transcripts = n_transcripts,
                  CDTS_5ss = CDTS_5ss_mean,
                  CDTS_3ss = CDTS_3ss_mean,
                  protein_coding = protein_coding,
                  mean_phastCons20way_5ss = phastCons20way_5ss_mean,
                  mean_phastCons20way_3ss = phastCons20way_3ss_mean)
  
  #########################
  ## LINEAR MODELS
  #########################
  
  
  fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ 
                    intron_length +
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
  
  fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ 
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
  
  model_names <- c("MSR_D", "MSR_A")
  
  
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
          axis.text.x = element_text(colour = "black", size = "9"),
          axis.text.y = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(colour = "black", size = "9"),
          legend.title = element_text(colour = "black", size = "9"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2)) +  
    ylab("Predictors") +
    geom_hline(yintercept = seq(from = 0.5,
                                to = length((fit_donor$coefficients %>% names)[-1]) + .5,
                                by = 1)) +
    guides(colour = guide_legend(#title = NULL,
                                 ncol = 2, 
                                 nrow = 1))
  
  
  if (save_results) {
    plotLM
    filename <- paste0(folder_name, "/images/lm_MSR.png")
    ggplot2::ggsave(filename = filename, 
                    width = 183, height = 183, units = "mm", dpi = 300)
  } else {
    plot %>% return()
  }
  
}

get_estimate_variance <- function() {
  

  df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/variance_estimate.rds"))
  graph_title <- "Distribution of the estimate values across 54 GTEx tissues"
  
  
  
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
    gather(tissue, feature, -type, -tissue) %>%
    mutate(type = "MSR_D")
  
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
    mutate(type = "MSR_A")
  
  MSR_tidy <- rbind(MSR_Donor_tidy, MSR_Acceptor_tidy) %>%
    mutate(type = factor(type, 
                         levels = c("MSR_D", "MSR_A")))
  
  
  # MSR_tidy <- MSR_tidy %>%
  #   filter(!(tissue %in% c("CDTS_5ss",
  #                          "CDTS_3ss",
  #                          "mean_phastCons20way_5ss",
  #                          "mean_phastCons20way_3ss")))
  
  plotTissuesLM <- ggplot(MSR_tidy, aes(tissue, feature)) + 
    geom_boxplot() +
    facet_grid(vars(type)) +
    #ggtitle(graph_title) +
    ylab("Distribution of the estimate") +
    xlab("Predictor") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          strip.text = element_text(colour = "black", size = "8"),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    scale_fill_manual(breaks = c("MSR_Donor","MSR_Acceptor"),
                      labels = c("MSR_D","MSR_A")) +
    theme(axis.text.x = element_text(angle = 70, 
                                     vjust = 1,
                                     hjust = 1)) +
    geom_hline(yintercept = 0,linetype='dotted')
  
  
  if (save_results) {
    plotTissuesLM
    
    ## SAVE THE RESULTS
    if (brain_tissues) {
      file_name <- "/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_brain_tissues.png"
    } else {
      file_name <- "/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_all_tissues.png"
    }
    
    ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
    
  } else {
    plot %>% return()
  }
}



####################################################
## RUNNING ANALYSES USING THE INTRON DATABASE ######
####################################################
  
  

## DISTANCES ----------------------------------------

summarise_distances_tissues <- function(project_id = "BRAIN",
                                        gtf_version = 105,
                                        all_tissues = F) {
  
  if (all_tissues) {
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  } else {
    all_projects <- project_id
  }
  
  all_data <- map_df(all_projects, function(project_id) {
  
    print(paste0(project_id))
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    all_data <- map_df(all_clusters, function(cluster) {
      
      print(cluster)
      
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/pipeline3/missplicing-ratio/",
                            cluster, "/v", gtf_version, "/")
      
      df <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds")) %>%
        as_tibble()
      
      df_tidy_donor <- df %>%
        filter(novel_type == "novel_donor") %>%
        pull(distance) %>%
        summary() %>% as.array %>%as.data.frame() %>%
        mutate(tissue = cluster) %>%
        mutate("NovelType" = "Novel Donor")
      
      df_tidy_donor <- rbind(df_tidy_donor,
                             data.frame(Var1 = "Mode (+)",
                                        Freq = get_mode(df %>%
                                                          filter(novel_type == "novel_donor", distance > 0) %>%
                                                          pull(distance)),
                                        tissue = cluster,
                                        "NovelType" = "Novel Donor"))
      
      df_tidy_donor <- rbind(df_tidy_donor,
                             data.frame(Var1 = "Mode (-)",
                                        Freq = get_mode(df %>%
                                                          filter(novel_type == "novel_donor", distance < 0) %>%
                                                          pull(distance)),
                                        tissue = cluster,
                                        "NovelType" = "Novel Donor"))
      ##############
      
      df_tidy_acceptor <- df %>%
        filter(novel_type == "novel_acceptor") %>%
        pull(distance) %>%
        summary() %>% 
        as.array %>%
        as.data.frame() %>%
        mutate(tissue = cluster) %>%
        mutate("NovelType" = "Novel Acceptor")
      
      df_tidy_acceptor <- rbind(df_tidy_acceptor,
                                data.frame(Var1 = "Mode (+)",
                                           Freq = get_mode(df %>%
                                                             filter(novel_type == "novel_acceptor", distance > 0) %>%
                                                             pull(distance)),
                                           tissue = cluster,
                                           "NovelType" = "Novel Acceptor"))
      
      df_tidy_acceptor <- rbind(df_tidy_acceptor,
                                data.frame(Var1 = "Mode (-)",
                                           Freq = get_mode(df %>%
                                                             filter(novel_type == "novel_acceptor", distance < 0) %>%
                                                             pull(distance)),
                                           tissue = cluster,
                                           "NovelType" = "Novel Acceptor"))
      
      
      
      return(rbind(df_tidy_donor, df_tidy_acceptor))
    })
    
    return(all_data)
    
  })
  
  saveRDS(object = all_data,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/summary_distances.rds"))
}

## MSR ----------------------------------------

summarise_MSR_tissues <- function(project_id = "BRAIN",
                                  gtf_version = 105,
                                  all_tissues = F) {
  
  if (all_tissues) {
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  } else {
    all_projects <- project_id
  }
  
  all_data <- map_df(all_projects, function(project_id) {
    
    print(paste0(project_id))
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    all_data <- map_df(all_clusters, function(cluster) {
      
      print(cluster)
      
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/pipeline3/missplicing-ratio/",
                            cluster, "/v", gtf_version, "/")
      
      
      df <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>% as_tibble()
      df %>% head()
      df %>% nrow()
      
      
      df_tidy_donor <- df %>%
        pull(ref_missplicing_ratio_tissue_ND) %>%
        summary() %>% 
        as.array %>%
        as.data.frame() %>%
        mutate(tissue = cluster) %>%
        mutate("MSR_type" = "MSR_Donor")
      
      
      ##############
      
      df_tidy_acceptor <- df %>%
        pull(ref_missplicing_ratio_tissue_NA) %>%
        summary() %>% 
        as.array() %>%
        as.data.frame() %>%
        mutate(tissue = cluster) %>%
        mutate("MSR_type" = "MSR_Acceptor")
      
      
      return(rbind(df_tidy_donor, df_tidy_acceptor))
    })
    
  })
  
  saveRDS(object = all_data,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/summary_MSRs.rds"))
}

## LINEAR MODELS ------------------------------------

summarise_lm_tissues <- function(gtf_version = 105,
                                 brain = F) {
  
  if (!brain) {
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  } else {
    all_projects <- project_id
  }
  
  all_data <- map_df(all_projects[-6], function(project_id) {
    
    print(paste0(project_id))
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    all_data <- map_df(all_clusters, function(cluster) {
    
      # cluster <- gtex_tissues[11]
      print(cluster)
      
      ## Load the IDB
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/pipeline3/missplicing-ratio/",
                            cluster, "/v", gtf_version, "/")
      
      idb <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
      
      idb <- idb %>%
        as.data.frame() %>%
        dplyr::distinct(ref_junID, .keep_all = T) %>%
        dplyr::rename(intron_length = width,
                      intron_5ss_score = ref_ss5score,
                      intron_3ss_score = ref_ss3score,
                      gene_length = gene_width,
                      gene_tpm = tpm_median_rct3,
                      gene_num_transcripts = n_transcripts,
                      protein_coding = protein_coding,
                      CDTS_5ss = CDTS_5ss_mean,
                      CDTS_3ss = CDTS_3ss_mean,
                      mean_phastCons20way_5ss = phastCons20way_5ss_mean,
                      mean_phastCons20way_3ss = phastCons20way_3ss_mean)
      
      ## Donor
      model <- lm(ref_missplicing_ratio_tissue_ND ~ 
                    intron_length + 
                    intron_5ss_score *
                    intron_3ss_score +
                    gene_tpm +
                    gene_length +
                    # u2_intron +
                    # clinvar +
                    protein_coding +
                    CDTS_5ss + 
                    CDTS_3ss + 
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss +
                    gene_num_transcripts,  
                  data = idb)
      
      
      
      
      MSR_Donor <- data.frame(tissue = cluster,
                              "MSR_Type" = "Donor",
                              intron_length = model$coefficients["intron_length"] %>% unname(),
                              intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                              intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                              gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                              gene_length = model$coefficients["gene_length"] %>% unname(),
                              #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),
                              protein_coding = model$coefficients["protein_coding"] %>% unname(),
                              gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                              intron_5ss_scoreintron_3ss_score = model$coefficients["intron_5ss_score:intron_3ss_score"] %>% unname(),
                              CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                              CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                              mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                              mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname())
      
      
      ## Acceptor
      model <- lm(ref_missplicing_ratio_tissue_NA ~
                    intron_length + 
                    intron_5ss_score *
                    intron_3ss_score +
                    gene_tpm +
                    gene_length +
                    # u2_intron +
                    # clinvar +
                    protein_coding +
                    CDTS_5ss + 
                    CDTS_3ss + 
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss +
                    gene_num_transcripts,  
                  data = idb)
      
      MSR_Acceptor <- data.frame(tissue = cluster,
                                 "MSR_Type" = "Acceptor",
                                 intron_length = model$coefficients["intron_length"] %>% unname(),
                                 intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                                 intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                                 gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                                 gene_length = model$coefficients["gene_length"] %>% unname(),
                                 #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),
                                 protein_coding = model$coefficients["protein_coding"] %>% unname(),
                                 gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                                 intron_5ss_scoreintron_3ss_score = model$coefficients["intron_5ss_score:intron_3ss_score"] %>% unname(),
                                 CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                                 CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                                 mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                                 mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname())
      
      return(rbind(MSR_Donor,MSR_Acceptor))
    })
  })
 
  
  saveRDS(object = all_data, 
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/summary_lm.rds"))
}

get_common_introns <- function (gtf_version = 105,
                                brain = F) {
  
  if (brain) {
    all_projects <- "BRAIN"
  } else {
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
  }
  
  print(paste0(Sys.time(), " - getting unique and common junctions across clusters..."))
  
  
  ## Getting all introns that are common across all GTEx tissues -------------------------------------------
  
  all_introns <- list()
  
  for (project_id in all_projects) {
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    for(cluster in all_clusters) { 
      
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/pipeline3/missplicing-ratio/",
                            cluster, "/v", gtf_version, "/")
      
      all_introns[[cluster]] <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
        distinct(ref_junID) %>% 
        pull()
      
      print(paste0(Sys.time(), " - intron IDs collected from '", cluster, "'"))
      
    }
  }
  
  common_introns <- data.frame(ref_junID = Reduce(intersect,  all_introns))
  common_introns %>% head()
  common_introns %>% nrow()
  
  all_common_introns_data <- map_df(all_projects, function(project_id) {
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    map_df(all_clusters, function(cluster) { 
      
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/pipeline3/missplicing-ratio/",
                            cluster, "/v", gtf_version, "/")
      
      readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
        filter(ref_junID %in% common_introns$ref_junID) %>%
        unnest(gene_name) %>%
        dplyr::select(ref_junID, ref_type, gene_name) %>%
        distinct(ref_junID, .keep_all = T) %>% 
        return()
      
    })
    
  })
  
  
  if (brain) {
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                        project_id, "/results/pipeline3/common_introns.rds")
  } else{
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/common_introns.rds")
  }
  
  saveRDS(object = all_common_introns_data, file = file_name)
  
}

get_estimate_variance <- function(gtf_version = 105,
                                  brain = F) {
  
  ## LOAD COMMON INTRONS
  
  if (brain) {
    all_projects <- "BRAIN"
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                        project_id, "/results/pipeline3/common_introns.rds")
  } else{
    all_projects <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds")
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/common_introns.rds")
    
  }
  
  common_introns <- readRDS(file = file_name)
  
  
  ## LOOP PER PROJECT
  
  df_estimates <- map_df(all_projects[-6], function(project_id) {
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    map_df(all_clusters, function(cluster) {
    
      # tissue <- clusters[1]
      print(cluster)
      
      ## Load the IDB
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/pipeline3/missplicing-ratio/", 
                            cluster, "/v", gtf_version, "/")
      
      idb <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) #%>%
      
      idb <- idb  %>%
        filter(ref_junID %in% common_introns$ref_junID) %>%
        dplyr::distinct(ref_junID, .keep_all = T) %>%
        filter(u2_intron == T) %>%
        dplyr::rename(intron_length = width,
                      intron_5ss_score = ref_ss5score,
                      intron_3ss_score = ref_ss3score,
                      gene_length = gene_width,
                      gene_tpm = tpm_median_rct3,
                      protein_coding  = protein_coding,
                      gene_num_transcripts = n_transcripts,
                      CDTS_5ss = CDTS_5ss_mean,
                      CDTS_3ss = CDTS_3ss_mean,
                      mean_phastCons20way_5ss = phastCons20way_5ss_mean,
                      mean_phastCons20way_3ss = phastCons20way_3ss_mean)  %>%
        filter(gene_tpm > 0)  %>% 
        filter(intron_length < gene_length)
      
      # idb %>% nrow() %>% print()
      
      
      ## Donor
      model <- lm(ref_missplicing_ratio_tissue_ND ~ 
                    intron_length + 
                    intron_5ss_score *
                    intron_3ss_score +
                    gene_tpm +
                    gene_length +
                    protein_coding +
                    CDTS_5ss + 
                    CDTS_3ss + 
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss +
                    gene_num_transcripts,  
                  data = idb)
      
      
      #model %>% summary() %>% print()
      #MSR_Donor_list[[tissue]] <- model
      
      
      #ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
      #model$coefficients[-ind_sign] <- 0
      
      MSR_Donor <- data.frame(tissue = cluster,
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
                                    mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname())
      
      
      ## Acceptor
      model <- lm(ref_missplicing_ratio_tissue_NA ~
                    intron_length + 
                    intron_5ss_score *
                    intron_3ss_score +
                    gene_tpm +
                    gene_length +
                    # u2_intron +
                    # clinvar +
                    protein_coding +
                    CDTS_5ss + 
                    CDTS_3ss + 
                    mean_phastCons20way_5ss +
                    mean_phastCons20way_3ss +
                    gene_num_transcripts,  
                  data = idb)
      #model %>% summary() %>% print()
      #MSR_Acceptor_list[[tissue]] <- model
      
      #ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
      #model$coefficients[-ind_sign] <- 0
      
      MSR_Acceptor <- data.frame(tissue = cluster,
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
                                       mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname())
      
      return(rbind(MSR_Donor %>% mutate(type = "MSR_Donor"), 
                   MSR_Acceptor %>% mutate(type = "MSR_Acceptor")))
      
    })
  })
  
  
  #############################
  ## SAVE RESULTS
  #############################
    
  if (brain) {
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                        project_id, "/results/pipeline3/variance_estimate")
  } else{
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/variance_estimate")
  }
  
  saveRDS(object = df_estimates, file = paste0(file_name, ".rds"))
  
}




## CHANGES IN TMP VALUES ----------------------------

get_changes_tpm <- function() {
  
  gtex_tissues <-  readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")
  
  clusters <- gtex_tissues[c(1:5,11,12,17)]
  
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/")
  
  ## Get common introns across two main brain tissues
  db_list <- NULL
  for (cluster in clusters) {
    
    print(cluster)  
    # cluster <- clusters[1]
    file_stats <- paste0(folder_root, "/", cluster, "/v104/", cluster, "_db_introns.rds")
    df_introns <- readRDS(file = file_stats)
    
    db_list[[cluster]] <- df_introns$ref_junID
    #select(ref_junID, ref_missplicing_ratio_tissue_ND, ref_missplicing_ratio_tissue_NA, tpm)
  }
  common_introns <- data.frame(ref_junID = Reduce(intersect,  db_list))
  
  
  db_common <- data.frame(ref_junID = as.character(),
                          gene =  as.character(),
                          tpm = as.double(),
                          MSR_D = as.double(),
                          MSR_A = as.double(),
                          tissue = as.double())
  
  for (cluster in clusters) {
    
    print(cluster)  
    
    # cluster <- clusters[1]
    file_stats <- paste0(folder_root, "/", cluster, "/v104/", cluster, "_db_introns.rds")
    df_introns <- readRDS(file = file_stats)
    
    df_introns <- df_introns %>%
      filter(ref_junID %in% common_introns$ref_junID) %>% 
      filter(u2_intron == T | u12_intron == T) %>%
      dplyr::select(ref_junID, tpm, ref_missplicing_ratio_tissue_ND, ref_missplicing_ratio_tissue_NA) %>%
      distinct(ref_junID, .keep_all = T)
    
    df_introns %>% nrow %>% print
    db_common <- rbind(db_common,
                       df_introns %>% 
                         mutate(tissue = cluster) %>%
                         dplyr::select(ref_junID, tpm, MSR_D = ref_missplicing_ratio_tissue_ND, MSR_A=  ref_missplicing_ratio_tissue_NA, tissue))
    
  }
  
  
  db_common %>% nrow()
  db_common %>% head()
  
  db_common <- split(db_common, db_common$tissue)
  
  i <- 1
  target <- 4
  
  for (i in c(1:length(db_common))) {
    
    if (target != i){
      db_common_merged <- merge(x = db_common[[target]],
                                y = db_common[[i]],
                                by = "ref_junID")
      
      db_common_merged %>% head()
      db_common_merged <- db_common_merged %>%
        mutate(delta_tpm = tpm.x - tpm.y) %>%
        mutate(delta_MSR_D = MSR_D.x - MSR_D.y) %>%
        mutate(delta_MSR_A = MSR_A.x - MSR_A.y) 
      
      
      db_common_merged %>% head()
      
      print((db_common %>% names())[target])
      print((db_common %>% names())[i])
      
      cor.test(x = db_common_merged$delta_tpm,
               y = db_common_merged$delta_MSR_D) %>% print()
      
      cor.test(x = db_common_merged$delta_tpm,
               y = db_common_merged$delta_MSR_A) %>% print()
      
    }
  }
  
  
  
  
  
  
  
  
  
}


## GENE ENRICHMENT ----------------------------------

get_GO_enrichment <- function(cluster,
                              all_split_reads_details_104,
                              folder_name) {
  
  library(gprofiler2)
  #library(rutils)
  
  
  #########################################################
  ## LOAD THE MIS-SPLICING RATIO DATA AND REPLACE BY ZEROES
  #########################################################
  
  df <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/all_introns_tidy.rds"))
  df %>% head()
  df %>% nrow()
  
  
  
  ## ENRICHMENT AT NEVER -----------------------------------------------------------
  
  
  df_never <- df %>%
    filter(final_type == "never") %>%
    filter(final_type != "acceptor") %>%
    unnest(gene_id) 
  
  df_never %>% head() %>% dplyr::select(gene_id)
  
  
  ranked_genes <- df_never$gene_id %>% unique()
  ranked_genes %>% length() %>% print()
  
  gene_enrichment_never_desc <- gprofiler2::gost(query = ranked_genes,
                                                 organism = "hsapiens",
                                                 ordered_query = F,
                                                 correction_method = "bonferroni",
                                                 sources = c("GO"))
  
  
  gene_enrichment_never_desc$result
  
  
  
  
  ## ENRICHMENT AT DONOR -----------------------------------------------------------
  
  
  df_donor <- df %>%
    filter(final_type != "never") %>%
    filter(final_type != "acceptor") %>%
    unnest(gene_id) %>%
    group_by(gene_id) %>%
    mutate(missplicing_donor = ref_missplicingratio_ND_tissues %>% min()) %>%
    ungroup() %>%
    arrange(ref_missplicingratio_ND_tissues)
  
  df_donor %>% head() %>% dplyr::select(gene_id)
  
  
  ranked_genes <- df_donor$gene_id %>% unique()
  ranked_genes %>% length() %>% print()
  
  gene_enrichment_donor_desc <- gprofiler2::gost(query = ranked_genes,
                                                 organism = "hsapiens",
                                                 ordered_query = T,
                                                 correction_method = "bonferroni",
                                                 sources = c("GO"))
  
  gene_enrichment_donor_desc$result
  
  
  
  ## ENRICHMENT AT ACCEPTOR -----------------------------------------------------------
  
  
  df_acceptor <- df %>%
    filter(final_type != "never") %>%
    filter(final_type != "donor") %>%
    unnest(gene_id) %>%
    group_by(gene_id) %>%
    mutate(missplicing_acceptor = ref_missplicingratio_NA_tissues %>% min()) %>%
    ungroup() %>%
    arrange(ref_missplicingratio_NA_tissues)
  
  df_acceptor %>% head() %>% dplyr::select(gene_id)
  
  
  ranked_genes <- df_acceptor$gene_id %>% unique()
  ranked_genes %>% length() %>% print()
  
  gene_enrichment_acceptor_desc <- gprofiler2::gost(query = ranked_genes,
                                                    organism = "hsapiens",
                                                    ordered_query = T,
                                                    correction_method = "bonferroni",
                                                    sources = c("GO"))
  
  gene_enrichment_acceptor_desc$result
  
}


## PROPORTION OF INTRONS ----------------------------

get_proportion_introns_clusters <- function(clusters,
                                            folder_root,
                                            folder_results) {
  
  
  
  df_missplicing_all <- map_df(clusters, function(cluster) {  # cluster <- clusters[11]
    
    
    ## Load the intron database
    
    print(paste0(Sys.time(), " - loading ", cluster, " data..."))
    folder_name <- paste0(folder_root, "/", cluster, "/")
    
    df <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
    
    
    ## Return a dataframe with the intron proportions
    
    return(data.frame(tissue = cluster,
                      
                      ref_junc = df %>%
                        dplyr::distinct(ref_junID) %>%
                        nrow(),
                      
                      genes = df %>%
                        dplyr::distinct(gene_id) %>%
                        nrow(),
                      
                      both_junc = df %>%
                        filter(ref_type == "both") %>% nrow(),
                      
                      donor_junc = df %>%
                        filter(ref_type == "donor") %>% nrow(),
                      
                      acceptor_junc =  df %>%
                        filter(ref_type == "acceptor") %>% nrow(),
                      
                      none_junc = df %>%
                        filter(ref_type == "never") %>% nrow(),
                      
                      
                      prop_both = (df %>%
                                     filter(ref_type == "both") %>% nrow()) / df %>%
                        dplyr::distinct(ref_junID) %>%
                        nrow(),
                      
                      prop_donor = (df %>%
                                      filter(ref_type == "donor") %>% nrow()) / df %>%
                        dplyr::distinct(ref_junID) %>%
                        nrow(),
                      
                      prop_acceptor =  (df %>%
                                          filter(ref_type == "acceptor") %>% nrow()) / df %>%
                        dplyr::distinct(ref_junID) %>%
                        nrow(),
                      
                      prop_none = (df %>%
                                     filter(ref_type == "never") %>% nrow()) / df %>%
                        dplyr::distinct(ref_junID) %>%
                        nrow(),
                      
                      both_genes = df %>%
                        filter(ref_type == "both") %>% 
                        dplyr::distinct(gene_id) %>%
                        unname() %>% unlist() %>% length(),
                      
                      donor_genes = df %>%
                        filter(ref_type == "donor") %>% 
                        dplyr::distinct(gene_id) %>%
                        unname() %>% unlist() %>% length(),
                      
                      acceptor_genes = df %>%
                        filter(ref_type == "acceptor") %>% 
                        dplyr::distinct(gene_id) %>%
                        unname() %>% unlist() %>% length(),
                      
                      none_genes = df %>%
                        filter(ref_type == "never") %>% 
                        dplyr::distinct(gene_id) %>%
                        unname() %>% unlist() %>% length(),
                      
                      #both_unique_genes = setdiff(genes_both, c(genes_donor, genes_acceptor, genes_none)) %>% length(),
                      
                      #donor_unique_genes = setdiff(genes_donor, c(genes_both, genes_acceptor, genes_none)) %>% length(),
                      
                      #acceptor_unique_genes = setdiff(genes_acceptor, c(genes_both, genes_donor, genes_none)) %>% length(),
                      
                      #none_unique_genes = setdiff(genes_none, c(genes_both, genes_donor, genes_acceptor)) %>% length(),
                      
                      donor_min = df %>% dplyr::select(ref_missplicing_ratio_tissue_ND) %>% min(), 
                      donor_mode = df %>% pull(ref_missplicing_ratio_tissue_ND) %>% get_mode(), 
                      donor_median = df %>% pull(ref_missplicing_ratio_tissue_ND) %>% median(), 
                      donor_mean = df %>% pull(ref_missplicing_ratio_tissue_ND) %>% mean(),
                      donor_max = df %>% dplyr::select(ref_missplicing_ratio_tissue_ND) %>% max(), 
                      
                      acceptor_min = df %>% dplyr::select(ref_missplicing_ratio_tissue_NA) %>% min(), 
                      acceptor_mode = df %>% pull(ref_missplicing_ratio_tissue_NA) %>% get_mode(), 
                      acceptor_median = df %>% pull(ref_missplicing_ratio_tissue_NA) %>% median(), 
                      acceptor_mean = df %>% pull(ref_missplicing_ratio_tissue_NA) %>% mean(),
                      acceptor_max = df %>% dplyr::select(ref_missplicing_ratio_tissue_NA) %>% max()
    ))
    
  })
  
  df_missplicing_all %>%  
    print()
  
  saveRDS(object = df_missplicing_all,
          file = paste0(folder_results, "/db_intron_prop_all_clusters.rds"))
  
}

plot_proportion_introns_clusters <- function(clusters,
                                             folder_results,
                                             save_result = F) {
  
  
  library(dplyr)
  library(forcats)
  
  
  if (any(clusters == "PD")) {
    label = round(df$data * 100, digits = 2)
    xlabel = "sample type"
  } else {
    label = ""
    xlabel = "tissue"
  }
  
  
  ## Load data -----------------------------------------------------------------------
  
  df_missplicing_all <- readRDS(file = paste0(folder_results, "/db_intron_prop_all_clusters.rds"))
  
  
  
  ## Plot results --------------------------------------------------------------------
  
  df <- data.frame(data = df_missplicing_all$prop_both,
                   tissue = df_missplicing_all$tissue,
                   type = "both")
  df <- rbind(df, 
              data.frame(data = df_missplicing_all$prop_acceptor,
                         tissue = df_missplicing_all$tissue,
                         type = "acceptor"))
  
  df <- rbind(df, 
              data.frame(data = df_missplicing_all$prop_donor,
                         tissue = df_missplicing_all$tissue,
                         type = "donor"))
  
  
  df <- rbind(df, 
              data.frame(data = df_missplicing_all$prop_none,
                         tissue = df_missplicing_all$tissue,
                         type = "never"))
  
  
  
  df$type = factor(df$type, 
                   levels = c( "never", "donor", "acceptor", "both" ))
  
  
  df2 = df %>% 
    ungroup() %>%
    arrange(type , desc(data)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df2$tissue), pattern = "Brain"), "red", "black")
  
  
  
  ggplot(df2, aes(x = tissue, 
                  y = data, 
                  group = type, 
                  fill = type)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    #coord_flip() +
    theme_light() +
    ylab("proportion of introns") +
    xlab(xlabel) +
    scale_fill_viridis_d()  +
    ggtitle("Proportion of intron type per GTEx tissue.") +
    #scale_x_discrete("tissue", breaks = breaks) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          #axis.text.y = element_text(color = colours),
          axis.text.x = element_text(color = colours,
                                     angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    scale_fill_manual(values = c("#FDE725FF", "#35B779FF", "#31688EFF","#440154FF"),
                      breaks = c("never", "acceptor", "donor", "both"),
                      labels = c("never", "acceptor", "donor", "both")) +
    guides(fill = guide_legend(title = "Intron type: ",
                               ncol = 4, 
                               nrow = 1)) 
  
  
  
  
  
  if (save_result) {
    file_name <- paste0(folder_results, "/proportion_junction_tissue_never_ordered.png")
    ggplot2::ggsave(filename = file_name,
                    width = 183, height = 183, units = "mm", dpi = 300)
  }
  
  
}



## COMMON JUNCTIONS ----------------------------------




QC_common_junctions <- function(clusters = gtex_tissues[c(7:17)],
                                folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/",
                                folder_save_name) {
  
  all_introns_tidy <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_brain_tidy.rds")
  
  intron_position <- sample(x = c(1:(all_introns_tidy %>% nrow())),size = 1)
  
  ## Filter the common intron randomly selected
  intron <- all_introns_tidy[intron_position,]
  
  ## Check that each single tissue from the IDB contains that intron
  all_introns <- map_df(clusters, function(cluster) {
    folder_name <- paste0(folder_root, "/", cluster)
    df <- readRDS(file = paste0(folder_name, "/v104/", cluster, "_db_introns.rds")) %>%
      as.data.frame() %>%
      dplyr::filter(ref_junID == intron$ref_junID) %>%
      mutate(tissue = cluster)
    
    return(df)
  })
  
  ## Count the number of tissues
  if (all_introns$tissue %>% unique() %>% length() != clusters %>% length()) {
    print("Error! The selected intron has not been found in all tissues!")
  }
  
  all_introns$ref_type %>% unique()
  intron$final_type %>% unique()
  intron$ref_missplicing_ratio_tissue_ND
  intron$ref_missplicing_ratio_tissue_NA
  intron$ref_missplicingratio_ND_tissues
  intron$ref_missplicingratio_NA_tissues
  
}

GO_enrichment_common_junctions <- function() {
  
  if (!exists("all_introns_tidy")) {
    all_introns_tidy <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/all_introns_brain_tidy.rds")
    all_introns_tidy %>% head()
  }
  
  
  both_misspliced_junc <- all_introns_tidy %>%
    filter(final_type == "both") 
  both_misspliced_junc %>% distinct(ref_junID, .keep_all = T) %>% nrow()
  both_genes <- both_misspliced_junc %>% unnest(gene_name) %>% distinct(gene_name) %>% pull()
  
  
  donor_misspliced_junc <- all_introns_tidy %>%
    filter(final_type == "donor")
  donor_misspliced_junc %>% distinct(ref_junID, .keep_all = T) %>% nrow()
  donor_genes <- donor_misspliced_junc %>% unnest(gene_name) %>% distinct(gene_name) %>% pull()
  
  
  acceptor_misspliced_junc <- all_introns_tidy %>%
    filter(final_type == "acceptor")
  acceptor_misspliced_junc %>% distinct(ref_junID, .keep_all = T) %>% nrow()
  acceptor_genes <- acceptor_misspliced_junc %>% unnest(gene_name) %>% distinct(gene_name) %>% pull()
  
  
  never_misspliced_junctions <- all_introns_tidy %>%
    filter(final_type == "never")
  never_misspliced_junctions %>% distinct(ref_junID, .keep_all = T) %>% nrow()
  never_genes <- never_misspliced_junctions %>% unnest(gene_name) %>% distinct(gene_name) %>% pull()
  
  
  ## UNIQUE GENES
  both_genes_unique <- setdiff(both_genes, c(never_genes, donor_genes, acceptor_genes))
  donor_genes_unique <- setdiff(donor_genes, c(never_genes, both_genes, acceptor_genes))
  acceptor_genes_unique <- setdiff(acceptor_genes, c(never_genes, donor_genes, both_genes))
  never_genes_unique <- setdiff(never_genes, c(both_genes, donor_genes, acceptor_genes))
  
  
  
  
  #######################################
  ## GO ENRICHMENT
  #######################################
  
  all_genes <- all_introns_tidy %>% unnest(gene_name) %>% distinct(gene_name) %>% pull() 
  
  ## 'Both' type genes
  gene_enrichment_both <- gprofiler2::gost(query = both_genes_unique,
                                           custom_bg = all_genes,
                                           organism = "hsapiens",
                                           ordered_query = F,
                                           correction_method = "bonferroni",
                                           sources = c("GO"))
  gene_enrichment_both$result %>% dplyr::select(p_value, term_name)
  
  ## 'Donor' type genes
  gene_enrichment_donor <- gprofiler2::gost(query = donor_genes_unique,
                                            custom_bg = all_genes,
                                            organism = "hsapiens",
                                            ordered_query = F,
                                            correction_method = "bonferroni",
                                            sources = c("GO"))
  gene_enrichment_donor$result %>% dplyr::select(p_value, term_name)
  
  
  ## 'Acceptor' type genes
  gene_enrichment_acceptor <- gprofiler2::gost(query = acceptor_genes_unique,
                                               custom_bg = all_genes,
                                               organism = "hsapiens",
                                               ordered_query = F,
                                               correction_method = "bonferroni",
                                               sources = c("GO"))
  gene_enrichment_acceptor$result %>% select(p_value, term_name)
  
  
  
  ## 'Never' type genes
  gene_enrichment_never <- gprofiler2::gost(query = never_genes_unique,
                                            custom_bg = all_genes,
                                            organism = "hsapiens",
                                            ordered_query = F,
                                            correction_method = "bonferroni",
                                            sources = c("GO"))
  gene_enrichment_never$result %>% dplyr::select(p_value, term_name)
  
  #never_genes_unique_name <- gene_ID_to_gene_name(never_genes_unique)
}

cdts_cons_scores_common_junctions <- function() {
  
  if (!exists("all_introns_tidy")) {
    all_introns_tidy <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/all_introns_brain_tidy.rds")
    all_introns_tidy %>% head()
  }
  
  
  both_misspliced_junc <- all_introns_tidy %>%
    filter(final_type == "both") 
  
  
  donor_misspliced_junc <- all_introns_tidy %>%
    filter(final_type == "donor")
  
  acceptor_misspliced_junc <- all_introns_tidy %>%
    filter(final_type == "acceptor")
  
  
  never_misspliced_junctions <- all_introns_tidy %>%
    filter(final_type == "never")
  
  
  #################################
  ## PLOT CONSERVATION
  #################################
  
  df_plot_conservation <- rbind(both_misspliced_junc,
                                never_misspliced_junctions) %>%
    dplyr::select(mean_phastCons20way_5ss, 
           mean_phastCons20way_3ss,
           final_type) %>%
    gather(key = cons_type, value = cons_value, -final_type) %>%
    mutate(cons_type = cons_type %>% as.factor())
  
  df_plot_conservation %>% head()
  
  
  
  ggplot(df_plot_conservation) +
    geom_density(aes(x = cons_value, 
                     fill = final_type), alpha = 0.6) +
    facet_grid(vars(cons_type))+ 
    theme_light() +
    #scale_fill_manual(values = c("#35B779FF","#440154FF"),
    #                  breaks = c("#35B779FF","#440154FF"),
    #                  labels = c("novel donor","novel acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  
  #################################
  ## PLOT CONSTRAINT
  #################################
  
  df_plot_CDTS <- rbind(both_misspliced_junc,
                        never_misspliced_junctions) %>%
    dplyr::select(CDTS_5ss, CDTS_3ss, final_type) %>%
    gather(key = CDTS_type, value = CDTS_value, -final_type)
  
  df_plot_CDTS %>% head()
  
  
  
  ggplot(df_plot_CDTS) +
    geom_density(aes(x = CDTS_value, 
                     fill = final_type), alpha = 0.6) +
    facet_grid(vars(CDTS_type))+ 
    theme_light() +
    #scale_fill_manual(values = c("#35B779FF","#440154FF"),
    #                  breaks = c("#35B779FF","#440154FF"),
    #                  labels = c("novel donor","novel acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
}



## RBP HeatMap ----------------------------------------


RBP_expression_heatmap <- function() {
  
  gtex_tpm <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_gene_median_tpm_tidy.rds")
  
  ## Filter the TPM file to obtain only the expression levels for the NMD genes -----------------------------
  NMD_genes <- data.frame(id = c("ENSG00000005007", "ENSG00000151461", "ENSG00000169062", "ENSG00000125351"),
                          name = c("UPF1", "UPF2", "UPF3A", "UPF3B"))
  
  RBP_genes <- data.frame(id = c("ENSG00000115524", "ENSG00000160201", "ENSG00000063244", "ENSG00000136450"),
                          name = c("SF3B1", "U2AF1", "U2AF2", "SRSF1"))
  
  genes <- rbind(NMD_genes, RBP_genes)
  
  df_NMD_tpm <- gtex_tpm %>%
    filter(gene_id %in% genes$id) %>% 
    gather(key = tissue, value = tpm, -gene_id)
  
  df_NMD_tpm <- merge(x = df_NMD_tpm,
                      y = genes,
                      by.x = "gene_id",
                      by.y = "id") %>%
    dplyr::rename(gene_name = name)
  
  
  ggplot(df_NMD_tpm, aes(x=tissue, y=gene_name)) +
    geom_tile(aes(fill = tpm)) +
    theme_light() +
    ylab("gene name")+
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "14"),
          axis.text.x = element_text(angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          legend.text = element_text(size = "14"),
          legend.title = element_text(size = "14")) 
}



##########################
## NOVEL EXON DISCOVERY ##
##########################

find_all_original_introns <- function(cluster,
                                      samples,
                                      split_read_counts,
                                      all_split_reads_details_104,
                                      folder_name) {
  
  ## Get number of introns (annotated junctions) in ensembl v104
  all_split_reads_details_104 %>%
    filter(type == "annotated") %>%
    nrow()
  
  ## Get number of introns (annotated junctions) in my intron database
  db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  db_introns %>% nrow()
  
  ## Make sure all introns from my intron database can be found on the ensembl v104 object
  if ((db_introns %>% nrow()) != (all_split_reads_details_104 %>%
                                  filter(type == "annotated") %>%
                                  filter(junID %in% db_introns$ref_junID) %>%
                                  nrow())) {
    print("Error: some introns from the intron database cannot be found on the original ensembl v104 data.")
  }
  
  ## Get number of novel junctions from my novel junction database
  db_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  db_novel %>% nrow()
  
  
  ## Make sure all introns from my novel database can be found on the intron database
  if ((db_novel$ref_junID %>% unique() %>% length()) != (db_introns %>%
                                                         filter(ref_type != "never") %>%
                                                         distinct(ref_junID) %>%
                                                         nrow())) {
    print("Error: some introns from the novel database cannot be found on the intron database.")
  }
  
  
  ## Get stats from the number of each type of intron type
  
  donor <- db_introns %>%
    filter(ref_type == "donor")
  donor %>% nrow()
  
  acceptor <- db_introns %>%
    filter(ref_type == "acceptor")
  acceptor %>% nrow()
  
  both <- db_introns %>%
    filter(ref_type == "both")
  both %>% nrow()
  
  never <- db_introns %>%
    filter(ref_type == "never")
  never %>% nrow()
  
  (never %>% nrow()) + (donor %>% nrow()) + (acceptor %>% nrow()) + (both %>% nrow())
  
  
  
  ## Get introns from the original ensembl v104 data that have not been included on the intron database
  diff_introns_IDs <- setdiff(all_split_reads_details_104 %>% 
                                filter(type == "annotated") %>% 
                                dplyr::distinct(junID) %>%
                                pull(junID), 
                              db_introns$ref_junID)
  diff_introns_IDs %>% length()
  
  
  
  
  
  ## Get introns mis-spliced but not attached to any junction
  db_notpaired <- readRDS(file = paste0("~/PROJECTS/splicing-project/results/pipeline3/distances/", 
                                        cluster, "/not-misspliced/", cluster, "_all_misspliced_not_paired.rds"))
  db_notpaired %>% unique() %>% length()
  
  
  
  
  ## From them, get introns from the original ensembl v104 data that present different genes at both ends
  introns_gene_ambiguous <- all_split_reads_details_104 %>%
    filter(junID %in% diff_introns_IDs) %>%
    filter(!(junID %in% db_notpaired)) %>%
    filter(type == "annotated") %>%
    filter(gene_id_end %>% as.character() != gene_id_start %>% as.character()) %>%
    filter(gene_id_end %>% as.character() != gene_id_start %>% as.character()) %>%
    distinct(junID)
  introns_gene_ambiguous %>% nrow()
  
  
  intersect(introns_gene_ambiguous$junID,db_notpaired)
  
  ## Get the rest of introns 
  introns_other <- all_split_reads_details_104 %>%
    filter(type == "annotated") %>%
    filter(junID %in% diff_introns_IDs) %>%
    filter(!(junID %in% db_notpaired)) %>%
    filter(!(junID %in% introns_gene_ambiguous$junID))
  introns_other %>% nrow()
  
  
  
  ### NOVEL DATA
  
  
  
  
  
  novel_donor <- db_novel %>%
    filter(ref_junID %in% donor$ref_junID) 
  novel_donor$novel_junID %>% unique %>% length()
  novel_donor$ref_junID %>% unique %>% length()
  
  
  
  
  
  
  novel_acceptor <- db_novel %>%
    filter(ref_junID %in% acceptor$ref_junID)
  novel_acceptor$novel_junID %>% unique %>% length()
  novel_acceptor$ref_junID %>% unique %>% length()
  
  
  
  both_novel <- db_novel %>%
    filter(ref_junID %in% both$ref_junID) 
  both_novel$novel_junID %>% unique %>% length()
  both_novel$ref_junID %>% unique %>% length()
  
  
  (both_novel$novel_junID %>% unique() %>% length()) + 
    (novel_acceptor$novel_junID %>% unique() %>% length()) + 
    (novel_donor$novel_junID %>% unique() %>% length())
  
  
  ## Get novel junctions from the original ensembl v104 data that present different genes at both ends
  
  diff_novel_junc <- all_split_reads_details_104 %>%
    filter(type %in% c("novel_donor", "novel_acceptor")) %>%
    filter(!(junID %in% (db_novel$novel_junID %>% unique())))
  diff_novel_junc %>% nrow()
  
  novel_gene_ambiguous <- diff_novel_junc %>%
    filter(gene_id_end %>% as.character() != gene_id_start %>% as.character()) %>%
    distinct(junID)
  novel_gene_ambiguous %>% nrow()
}

get_novel_annotation_incorporation <- function(cluster) {
  
  folder_new <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", cluster, "/")
  folder_old <- paste0("/home/sruiz/PROJECTS/splicing-project/results/base_data_old_version/", cluster, "/missplicing-ratio/")
  
  
  
  ## LOAD ANNOTATED INTRONS FROM NEWEST ANNOTATION
  file_name <- paste0(folder_new, "/", cluster, "_db_introns.rds")
  df_introns_new <- readRDS(file = file_name)
  df_introns_new %>% distinct(ref_junID) %>% nrow()
  
  file_name <- paste0(folder_new, "/", cluster, "_db_novel.rds")
  df_novel_new <- readRDS(file = file_name)
  df_novel_new  %>% distinct(novel_junID) %>% nrow()
  
  
  
  ## LOAD NOVEL JUNCTIONS DATABASE FROM OLDEST ANNOTATION
  file_name <- paste0(folder_old, "/", cluster, "_db_introns.rds")
  df_introns_old <- readRDS(file = file_name)
  df_introns_old %>% distinct(ref_junID) %>% nrow()
  
  file_name <- paste0(folder_old, "/", cluster, "_db_novel.rds")
  df_novel_old <- readRDS(file = file_name)
  df_novel_old %>% distinct(novel_junID) %>% nrow()
  
  
  
  
  df_novel_old %>%
    filter(novel_junID %in% df_introns_new$ref_junID) %>%
    mutate(modulo_distance = abs(distance) %% 3) %>%
    dplyr::select(seqnames,
           start,
           end, strand, width, novel_type, distance, modulo_distance,
           novel_n_individuals, novel_mean_counts, gene_name) %>%
    arrange(gene_name %>% as.character())
  
  df_introns_new %>%
    filter(gene_name == "SNCA") 
  
  
  ## PLOT NEW NOVEL JUNCTIONS EXPRESSED IN MORE THAN 96 OF INDIV. WITH MORE THAN 5 READS
  
  df_novel_new <- df_novel_new %>%
    mutate(color = ifelse(novel_n_individuals >= 96 & novel_mean_counts >= 5, "red", "black"))
  
  ggplot(data = df_novel_new) +
    geom_point(mapping = aes(x = novel_n_individuals,
                             y = novel_mean_counts, 
                             colour = color)) + 
    scale_color_manual(values = c("red" = "red", 
                                  "black" = "black"),
                       labels = c(paste0("individuals < 96 & mean_counts < 5 (",df_novel_new %>% filter(color == "black") %>% nrow,")"),  
                                  paste0("individuals >= 96 & mean_counts >= 5 (", df_novel_new %>% filter(color == "red") %>% nrow,")"))) +
    theme_light() +
    ylab("mean number of reads") +
    xlab("number of individuals") +
    ggtitle(paste0("Novel junctions expression - FCTX\nNovel junctions represented: ", 
                   df_novel_new %>% distinct(novel_junID) %>% nrow(), 
                   "\nReference transcriptome version: March 2021")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = NULL,
                                ncol = 1, 
                                nrow = 2))
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/Brain-FrontalCortex_BA9/images/novel_junctions_expression.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  ## SAVE NOVEL JUNCTIONS EXPRESSED IN MORE THAN 96 OF INDIV WITH MORE THAN 5 READS IN A BED FILE
  
  df_bed <- df_novel_new %>%
    filter(novel_n_individuals >= 96, novel_mean_counts >= 5) %>%
    dplyr::select(chrom = seqnames, chromStart = start, chromEnd = end, strand, name = novel_junID)
  write.table(df_bed, "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/new_introns.bed")
  
  ###############################
  
  # all_split_reads_details_97_old <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data_old_version/Brain-FrontalCortex_BA9/Brain-FrontalCortex_BA9_annotated_SR_details_length.rds")
  # 
  # 
  # all_old_novel <- all_split_reads_details_97_old %>%
  #   filter(type %in% c("novel_donor", "novel_acceptor"))
  # 
  # 
  # all_old_novel %>%
  #   filter(junID %in% df_introns_new$ref_junID) %>%
  #   mutate(modulo_distance = abs(distance) %% 3) %>%
  #   select(seqnames,
  #          start,
  #          end, strand, width, novel_type, distance, modulo_distance,
  #          novel_n_individuals, novel_mean_counts, gene_name) %>%
  #   arrange(gene_name %>% as.character())
  
}


get_novel_exon_discovery <- function(cluster,
                                     folder_name) {
  
  ## LOAD INTRONS DATABASE DATA
  file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
  df_introns <- readRDS(file = file_name)
  
  
  if (any(df_introns$seqnames %in% c("X", "Y", "MT"))) {
    print("Error: some of the reference introns have been found within the sex/MT chromosomes!")
  }
  
  df_introns <- df_introns %>%
    dplyr::filter(ref_type == "both")#,
  #gene_biotype == "protein_coding")
  
  df_introns %>% head()
  df_introns %>% nrow()
  
  
  df_introns %>% 
    dplyr::filter(gene_name == "SNCA")
  
  
  ## LOAD NOVEL JUNCTIONS DATABASE DATA
  file_name <- paste0(folder_name, "/", cluster, "_db_novel.rds")
  df_novel <- readRDS(file = file_name)
  df_novel %>% nrow()
  
  df_novel <- df_novel %>% 
    dplyr::filter(ref_junID %in% df_introns$ref_junID)
  
  df_novel %>% nrow()
  df_novel %>% head()
  
  
  
  df_novel %>% 
    dplyr::filter(gene_name == "SNCA") %>% distinct(novel_junID)
  
  
  
  
  ## ADD START END POSITION OF THE REF INTRON
  
  
  df_novel_ref <- merge(x = df_novel,
                        y = df_introns[, c("ref_junID", "start", "end")],
                        by = "ref_junID",
                        all.x = T)
  
  df_novel_ref <- df_novel_ref %>%
    dplyr::rename(start = start.x,
                  end = end.x,
                  ref_start = start.y,
                  ref_end = end.y)
  
  df_novel_ref %>% head
  
  
  
  
  ## NOVEL EXON DISCOVERY CRITERIA ------------------------------------------------------------------------------------
  
  x <- df_novel_ref %>% #filter(gene_name == "SNCA") %>%
    filter(novel_n_individuals >= 95, novel_mean_counts >= 5, distance < 0) %>%
    
    mutate(novel_acceptor_start = ifelse(novel_type == "novel_acceptor", start %>% as.integer(), 0)) %>%
    mutate(novel_acceptor_end = ifelse(novel_type == "novel_acceptor", end %>% as.integer(), 0)) %>%
    mutate(novel_donor_start = ifelse(novel_type == "novel_donor", start %>% as.integer(), 0)) %>%
    mutate(novel_donor_end = ifelse(novel_type == "novel_donor", end %>% as.integer(), 0)) %>%
    
    mutate(ref_start_exon = ifelse(strand == "-", novel_donor_end + 1, novel_acceptor_end + 1)) %>%
    
    mutate(ref_end_exon = ifelse(strand == "-", novel_acceptor_start - 1, novel_donor_start - 1)) %>%
    
    mutate(ref_start_exon = ifelse(is.na(ref_start_exon), 0, ref_start_exon)) %>%
    mutate(ref_end_exon = ifelse(is.na(ref_end_exon), 0, ref_end_exon)) 
  
  x %>% 
    filter(gene_name == "SNCA", ref_junID == "67197834") %>% 
    dplyr::select(ref_junID, novel_junID, novel_n_individuals, novel_type,
           ref_start,
           ref_end,
           novel_start = start,
           novel_end = end,
           ref_start_exon, ref_end_exon)
  
  
  
  y <- x %>% 
    group_by(ref_junID) %>%
    tidyr::expand(ref_start_exon, ref_end_exon) %>%
    ungroup() 
  
  
  y %>% 
    filter(ref_junID == "67197834")
  
  
  hits <- y %>% #filter(ref_junID == "67197834") %>%
    filter(ref_start_exon > 1, ref_end_exon > 1) %>%
    filter(ref_start_exon < ref_end_exon) %>%
    mutate(exon_modulo = (ref_start_exon - ref_end_exon) %% 3) %>%
    filter(exon_modulo == 0)
  
  hits %>% 
    filter(ref_junID == "18383134") %>% 
    dplyr::select(ref_junID, ref_start_exon, ref_end_exon)
  
  hits_complete <- x %>%
    filter((ref_junID %in% (hits$ref_junID %>% unique()) & ref_start_exon %in% (hits$ref_start_exon %>% unique())) |
             (ref_junID %in% (hits$ref_junID %>% unique()) & ref_end_exon %in% (hits$ref_end_exon %>% unique())) )
  
  
  
  hits_complete %>% 
    #group_by(ref_junID) %>%
    #mutate(mean_n_individuals = mean(novel_n_individuals)) %>%
    #filter(mean_n_individuals > 90) %>%
    #filter(gene_name == "SNCA") %>% 
    dplyr::select(gene_name, seqnames, ref_start, ref_end ) %>%
    distinct(gene_name, .keep_all = T) %>%
    arrange(gene_name %>% as.character())# %>% pull(gene_name) %>% unlist() %>% unname()
  
  
}

comparing_reference_transcriptome_versions <- function(cluster = "Brain-FrontalCortex_BA9") {
  
  
  
  ###############################
  
  ## Load recount2 original data
  
  if (!exists("all_recount2_split_reads")) {
    all_recount2_split_reads <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/annotated_split_reads.rda")
  } else {
    print("'all_recount2_split_reads' file already loaded!")
  }
  all_recount2_split_reads$junID %>% unique() %>% length()
  all_recount2_split_reads$chr %>% unique()
  
  split_read_counts <- fread(paste0("/data/recount/GTEx_SRP012682/gtex_full_split_read_count_table/", cluster, "/", cluster, ".csv"))
  split_read_counts$junID %>% unique() %>% length()
  
  split_read_counts %>%
    filter(junID %in% all_recount2_split_reads$junID) %>%
    nrow()
  
  
  
  ###############################
  
  folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/"
  
  all_split_reads_details_97 <- readRDS(file = paste0(folder_root, "base_data/", cluster, "/", cluster, "_annotated_SR_details_length_v97.rds"))
  all_split_reads_details_104 <- readRDS(file = paste0(folder_root, "base_data/", cluster, "/", cluster, "_annotated_SR_details_length_v104.rds"))
  db_introns_104 <- readRDS(file = paste0(folder_root, "pipeline3/missplicing-ratio/", cluster, "/", cluster, "_db_introns_tpm.rds"))
  db_novel_104 <- readRDS(file = paste0(folder_root, "pipeline3/missplicing-ratio/", cluster, "/", cluster, "_db_novel.rds"))
  
  ###############################
  
  ## VERSION JULY 2019 - STATISTICS
  all_split_reads_details_97$junID %>% unique() %>% length()
  all_split_reads_details_97$seqnames %>% unique()
  all_split_reads_details_97$type %>% unique()
  
  all_split_reads_details_97 %>%
    filter(type %in% c("annotated")) %>% nrow()
  all_split_reads_details_97 %>%
    filter(type %in% c("novel_donor", "novel_acceptor")) %>% nrow()
  all_split_reads_details_97 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>% nrow()
  all_split_reads_details_97 %>%
    filter(type %in% c("none")) %>% nrow()
  
  
  
  ## VERSION MARCH 2021 - STATISTICS
  all_split_reads_details_104$junID %>% unique() %>% length()
  #all_split_reads_details_97_new_raw$junID %>% unique() %>% length()
  all_split_reads_details_104$seqnames %>% unique()
  all_split_reads_details_104$type %>% unique()
  
  all_split_reads_details_104 %>%
    filter(type %in% c("annotated")) %>% nrow()
  all_split_reads_details_104 %>%
    filter(type %in% c("novel_donor", "novel_acceptor")) %>% nrow()
  all_split_reads_details_104 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>% nrow()
  all_split_reads_details_104 %>%
    filter(type %in% c("none")) %>% nrow()
  
  
  
  ## CROSSING DATA BETWEEN VERSIONS USING THE JUNCTION ID ----------------------------------------------------------------------------------------
  
  annotated_old_junID <- all_split_reads_details_97 %>%
    filter(type == "annotated") %>%
    pull(junID)
  novel_old_junID <- all_split_reads_details_97 %>%
    filter(type %in% c("novel_donor","novel_acceptor")) %>%
    pull(junID)
  none_old_junID <- all_split_reads_details_97 %>%
    filter(type == "none") %>%
    pull(junID)
  other_old_junID <- all_split_reads_details_97 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>%
    pull(junID)
  
  annotated_new_junID <- all_split_reads_details_104 %>%
    filter(type == "annotated") %>% 
    pull(junID)
  novel_new_junID <- all_split_reads_details_104 %>%
    filter(type %in% c("novel_donor","novel_acceptor")) %>% 
    pull(junID)
  none_new_junID <- all_split_reads_details_104 %>%
    filter(type == "none") %>% 
    pull(junID)
  other_new_junID <- all_split_reads_details_104 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>%
    pull(junID)
  
  
  ## Intersecting the 'old' version 4 categories with the 'annotated' new junctions
  intersect(annotated_old_junID, annotated_new_junID) %>% length()
  intersect(novel_old_junID, annotated_new_junID) %>% length()
  intersect(none_old_junID, annotated_new_junID) %>% length()
  intersect(other_old_junID, annotated_new_junID) %>% length()
  
  
  ## Intersecting the 'old' version 4 categories with the 'novel' new junctions
  intersect(annotated_old_junID, novel_new_junID) %>% length()
  intersect(novel_old_junID, novel_new_junID) %>% length()
  intersect(none_old_junID, novel_new_junID) %>% length()
  intersect(other_old_junID, novel_new_junID) %>% length()
  
  
  ## Intersecting the 'old' version 4 categories with the 'none' new junctions
  intersect(annotated_old_junID, none_new_junID) %>% length()
  intersect(novel_old_junID, none_new_junID) %>% length()
  intersect(none_old_junID, none_new_junID) %>% length()
  intersect(other_old_junID, none_new_junID) %>% length()
  
  
  ## Intersecting the 'old' version 4 categories with the 'other' new junctions
  intersect(annotated_old_junID, other_new_junID) %>% length()
  intersect(novel_old_junID, other_new_junID) %>% length()
  intersect(none_old_junID, other_new_junID) %>% length()
  intersect(other_old_junID, other_new_junID) %>% length()
  
  
  
  ## CROSSING DATA BETWEEN VERSIONS USING GENOMIC COORDINATES -----------------------------------------------------------------------------------
  
  
  annotated_old_GR <- all_split_reads_details_97 %>%
    filter(type == "annotated") %>%
    dplyr::select(seqnames, start, end, width, strand) %>% 
    GRanges()
  novel_old_GR <- all_split_reads_details_97 %>%
    dplyr::filter(type %in% c("novel_donor","novel_acceptor")) %>%
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  none_old_GR <- all_split_reads_details_97 %>%
    filter(type == "none") %>%
    dplyr::select(seqnames, start, end, width, strand) %>% 
    GRanges()
  other_old_GR <- all_split_reads_details_97 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>%
    dplyr::select(seqnames, start, end, width, strand) %>% 
    GRanges()
  
  annotated_new_GR <- all_split_reads_details_104 %>%
    filter(type == "annotated") %>%
    dplyr::select(seqnames, start, end, width, strand) %>% 
    GRanges()
  novel_new_GR <- all_split_reads_details_104 %>%
    filter(type %in% c("novel_donor","novel_acceptor"))  %>%
    dplyr::select(seqnames, start, end, width, strand) %>% 
    GRanges()
  none_new_GR <- all_split_reads_details_104 %>%
    filter(type == "none")  %>%
    dplyr::select(seqnames, start, end, width, strand) %>% 
    GRanges()
  other_new_GR <- all_split_reads_details_104 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>%
    dplyr::select(seqnames, start, end, width, strand) %>% 
    GRanges()
  
  
  
  
  ## Intersecting the 'old' version 4 categories with the 'annotated' new junctions
  findOverlaps(query = annotated_old_GR, subject = annotated_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = novel_old_GR, subject = annotated_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = none_old_GR, subject = annotated_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = other_old_GR, subject = annotated_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  
  
  ## Intersecting the 'old' version 4 categories with the 'novel' new junctions
  findOverlaps(query = annotated_old_GR, subject = novel_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = novel_old_GR, subject = novel_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = none_old_GR, subject = novel_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = other_old_GR, subject = novel_new_GR, type=c("equal"), ignore.strand=FALSE)%>% queryHits() %>% length()
  
  
  ## Intersecting the 'old' version 4 categories with the 'none' new junctions
  findOverlaps(query = annotated_old_GR, subject = none_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = novel_old_GR, subject = none_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = none_old_GR, subject = none_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = other_old_GR, subject = none_new_GR, type=c("equal"), ignore.strand=FALSE)%>% queryHits() %>% length()
  
  
  ## Intersecting the 'old' version 4 categories with the 'other' new junctions
  findOverlaps(query = annotated_old_GR, subject = other_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = novel_old_GR, subject = other_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = none_old_GR, subject = other_new_GR, type=c("equal"), ignore.strand=FALSE) %>% queryHits() %>% length()
  findOverlaps(query = other_old_GR, subject = other_new_GR, type=c("equal"), ignore.strand=FALSE)%>% queryHits() %>% length()
  
  
  
  ## FROM THE NOVEL JUNCTIONS IN 2019, COLOUR THEM DEPENDING ON THEIR CATEGORY IN 2021 ---------------------------------------------------------
  
  novel_97 <- all_split_reads_details_97 %>%
    filter(type %in% c("novel_donor","novel_acceptor")) 
  novel_97 %>% distinct(junID) %>% nrow()
  
  
  
  annotation_104 <- db_introns_104 %>%
    filter(ref_junID %in% novel_97$junID) %>% 
    distinct(ref_junID, .keep_all = T) %>%
    mutate(ref_p_individuals = (ref_n_individuals * 100) / 108)
  annotation_104 %>% nrow()
  
  
  novel_104 <- db_novel_104 %>%
    filter(novel_junID %in% novel_97$junID) %>% 
    distinct(novel_junID, .keep_all = T) %>%
    mutate(novel_p_individuals = (novel_n_individuals * 100) / 108)
  novel_104 %>% nrow()
  
  
  
  ggplot() +
    geom_point(data = novel_104,
               mapping = aes(x = novel_p_individuals,
                             y = novel_mean_counts, 
                             color = "black")) + 
    geom_point(data = annotation_104,
               mapping = aes(x = ref_p_individuals,
                             y = ref_mean_counts, 
                             color = "red")) + 
    #geom_abline(mapping = aes(slope = 0, intercept = 25, linetype = factor(3), fill = "red")) +
    ggtitle("Ensembl v104 junctions previously categorised as novel in v97\nFCTX") +
    ylab("mean number of reads") +
    xlab("percentage of individuals (%)") +
    
    scale_y_continuous(breaks = c(0,25,50,75,100,300,600)) +
    scale_x_continuous(breaks = c(0,25,50,75,100)) +
    scale_color_manual(values = c("red","black"),
                       breaks = c("red","black"),
                       labels = c(paste0(annotation_104 %>% distinct(ref_junID) %>% nrow(), " fully annotated"),
                                  paste0(novel_104 %>% distinct(novel_junID) %>% nrow(), " novel"))) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Junction type: ", ncol = 2, nrow = 1))
  
  # SAVE PLOT
  file_name <- paste0(folder_root, "pipeline3/missplicing-ratio/", cluster, "/images/", cluster, "_ensemblv104junc_novelv97.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ## GETTING CONTAMINATION RATES -------------------------------------------------------------------------------------------------------
  
  # MR - however that contamination rate is going to vary at different % individuals AND 
  # at different read counts, but at the moment you can't see that.
  
  # In order to see this, I think you want to deal with % individuals and counts separately, 
  # unless David has a better idea. Thinking about % individuals, bin into 0-10%, 11-20% etc. then in each bin calculate 
  # the proportion of the fully annotated junctions relative to novel according to Ensembl v104.
  
  # You might also want to plot this more as a cumulative distribution plot ie. as the % individual numbers increase,
  # what % of your fully annotated are accounted for. You can repeat but this time binning counts.
  # Does this make sense? The point is not so much that all the annotated junctions are at one end of the plot
  # but that the rate of annotated really rises here.
  
  novel_97 <- all_split_reads_details_97 %>%
    filter(type %in% c("novel_donor","novel_acceptor")) 
  novel_97 %>% distinct(junID) %>% nrow()
  novel_97 %>% head()
  
  
  introns_104 <- db_introns_104 %>%
    filter(ref_junID %in% novel_97$junID) %>% 
    distinct(ref_junID, .keep_all = T) %>%
    mutate(ref_p_individuals = (ref_n_individuals * 100) / 108) %>%
    mutate(ref_p_mean_counts = (ref_mean_counts * 100) / max(ref_mean_counts))
  introns_104 %>% nrow()
  introns_104 %>% head()
  
  
  novel_104 <- db_novel_104 %>%
    filter(novel_junID %in% novel_97$junID) %>% 
    distinct(novel_junID, .keep_all = T) %>%
    mutate(novel_p_individuals = (novel_n_individuals * 100) / 108)
  novel_104 %>% nrow()
  novel_104 %>% head()
  
  ggplot(introns_104) +
    geom_bar(aes(ref_p_individuals)) +
    scale_x_binned() +
    ggtitle(paste0(introns_104 %>% nrow(), " ensembl-v104 introns previously categorised as novel in v97. FCTX tissue.\n",
                   "Introns binned depending on the % of individuals expressed in.")) +
    xlab("percentage of individuals (%)") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") 
  
  file_name <- paste0(folder_root, "pipeline3/missplicing-ratio/", cluster, "/images/", cluster, "_ensemblv104introns_perc_individuals.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ggplot(introns_104) +
    geom_bar(aes(ref_p_mean_counts)) +
    scale_x_binned() +
    ggtitle(paste0(introns_104 %>% nrow(), " ensembl-v104 introns previously categorised as novel in v97. FCTX tissue.\n",
                   "Introns binned depending on the % of counts they expressed.")) +
    xlab("percentage of counts (%)") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") 
  file_name <- paste0(folder_root, "pipeline3/missplicing-ratio/", cluster, "/images/", cluster, "_ensemblv104introns_perc_counts.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  ## GETTING THE 220 BRAND-NEW ANNOTATED JUNCTIONS -------------------------------------------------------------------------------------
  
  brand_new_intronIDs <- setdiff(annotated_new_junID,
                                 c(intersect(annotated_old_junID, annotated_new_junID),
                                   intersect(novel_old_junID, annotated_new_junID),
                                   intersect(none_old_junID, annotated_new_junID),
                                   intersect(other_old_junID, annotated_new_junID)))
  
  brand_new_introns <- all_split_reads_details_104 %>%
    filter(junID %in% brand_new_intronIDs)
  
  
  brand_new_introns$type %>% unique
  brand_new_introns$gene_name_junc %>% unlist() %>% unname() %>% sort() %>% unique()
  
  
  ## GETTING THE NOVEL JUNCTIONS IN 2019 THAT ARE ANNOTATED IN 2021 -----------------------------------------------------------------------------
  
  
  old_novel_now_annotated <- intersect(novel_old_junID,
                                       c(annotated_new_junID, other_new_junID))
  
  old_novel_now_annotated_recount2 <- all_recount2_split_reads %>%
    filter(junID %in% old_novel_now_annotated)
  
  gtex_samples <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_used.rda")
  samples <- gtex_samples$`Brain-FrontalCortex_BA9`
  
  ## Get the count data
  split_read_counts_novel <- split_read_counts %>%
    filter(junID %in% old_novel_now_annotated_recount2$junID) %>%
    dplyr::select(junID, samples %>% as.character()) %>%
    #mutate_all(~replace(., is.na(.), 0)) %>%
    as_tibble()
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  
  split_read_counts_novel[,"novel_n_individuals"] <- matrixStats::rowCounts(split_read_counts_novel[, -c(1)] > 0, na.rm = T)
  split_read_counts_novel[,"novel_mean_counts"] <- rowMeans(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(), 1)], na.rm = T)
  split_read_counts_novel[,"novel_sum_counts"] <- rowSums(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(),
                                                                                      (split_read_counts_novel %>% ncol()) - 1,
                                                                                      1)], na.rm = T)
  
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  if (any(split_read_counts_novel[, "novel_n_individuals"] < 1)) {
    print("Error: some novel junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_novel <- split_read_counts_novel[, c(1, 
                                                         (split_read_counts_novel %>% ncol() - 2),
                                                         (split_read_counts_novel %>% ncol() - 1), 
                                                         split_read_counts_novel %>% ncol())]
  split_read_counts_novel %>% head()
  
  
  
  df_old_novel_now_annotated <- merge(x = old_novel_now_annotated_recount2,
                                      y = split_read_counts_novel,
                                      by.x = "junID", 
                                      by.y = "junID",
                                      all.x = T)
  
  
  
  
  ggplot(data = df_old_novel_now_annotated) +
    geom_point(mapping = aes(x = novel_n_individuals,
                             y = novel_mean_counts)) + 
    geom_abline(mapping=aes(slope=0, intercept=5, linetype=factor(3), color ="red")) +
    scale_y_continuous(breaks = c(0,5,25,50,75,100)) +
    theme_light() +
    ylab("mean number of reads") +
    xlab("number of individuals") +
    ggtitle(paste0("Novel junctions from ensembl 2019 incorporated into 2021 - FCTX\nNovel donor/acceptor in 2019 --> annotated/novel combo/exon skip in 2021.\nJunctions represented: ", 
                   df_old_novel_now_annotated %>% distinct(junID) %>% nrow())) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "none")
  
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/Brain-FrontalCortex_BA9/images/novel2019_annotated2021.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  
  
  ## OLD NOVEL DONOR JUNCTIONS IN 2019 THAT ARE ANNOTATED IN 2021 -----------------------------------------------------------------------------
  
  donor_old_junID <- all_split_reads_details_97 %>%
    filter(type == "novel_donor") %>%
    pull(junID)
  
  old_novel_now_annotated <- intersect(donor_old_junID,
                                       c(annotated_new_junID, other_new_junID))
  
  old_novel_now_annotated_recount2 <- all_recount2_split_reads %>%
    filter(junID %in% old_novel_now_annotated)
  
  gtex_samples <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_used.rda")
  samples <- gtex_samples$`Brain-FrontalCortex_BA9`
  
  ## Get the count data in the current cluster for the current NOVEL junction
  split_read_counts_novel <- split_read_counts %>%
    filter(junID %in% old_novel_now_annotated_recount2$junID) %>%
    dplyr::select(junID, samples %>% as.character()) %>%
    #mutate_all(~replace(., is.na(.), 0)) %>%
    as_tibble()
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  
  split_read_counts_novel[,"novel_n_individuals"] <- matrixStats::rowCounts(split_read_counts_novel[, -c(1)] > 0, na.rm = T)
  split_read_counts_novel[,"novel_mean_counts"] <- rowMeans(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(), 1)], na.rm = T)
  split_read_counts_novel[,"novel_sum_counts"] <- rowSums(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(),
                                                                                      (split_read_counts_novel %>% ncol()) - 1,
                                                                                      1)], na.rm = T)
  
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  if (any(split_read_counts_novel[, "novel_n_individuals"] < 1)) {
    print("Error: some novel junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_novel <- split_read_counts_novel[, c(1, 
                                                         (split_read_counts_novel %>% ncol() - 2),
                                                         (split_read_counts_novel %>% ncol() - 1), 
                                                         split_read_counts_novel %>% ncol())]
  split_read_counts_novel %>% head()
  
  
  
  df_old_novel_now_annotated <- merge(x = old_novel_now_annotated_recount2,
                                      y = split_read_counts_novel,
                                      by.x = "junID", 
                                      by.y = "junID",
                                      all.x = T)
  
  
  
  
  ggplot(data = df_old_novel_now_annotated) +
    geom_point(mapping = aes(x = novel_n_individuals,
                             y = novel_mean_counts)) + 
    geom_abline(mapping=aes(slope=0, intercept=5, linetype=factor(3), color ="red")) +
    scale_y_continuous(breaks = c(0,5,25,50,75,100,130)) +
    #scale_color_manual(values = c("red" = "red", 
    #                              "black" = "black"),
    #                   labels = c(paste0("individuals < 96 & mean_counts < 5 (",df_novel_new %>% filter(color == "black") %>% nrow,")"),  
    #                              paste0("individuals >= 96 & mean_counts >= 5 (", df_novel_new %>% filter(color == "red") %>% nrow,")"))) +
    theme_light() +
    
    ylab("mean number of reads") +
    xlab("number of individuals") +
    ggtitle(paste0("Novel DONOR junctions from ensembl 2019 incorporated into 2021 - FCTX\nNovel donor in 2019 --> annotated/novel combo/exon skip in 2021.\nJunctions represented: ", 
                   df_old_novel_now_annotated %>% distinct(junID) %>% nrow())) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "none")
  
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/Brain-FrontalCortex_BA9/images/noveldonor2019_annotated2021.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ## OLD NOVEL ACCEPTOR JUNCTIONS IN 2019 THAT ARE ANNOTATED IN 2021 -----------------------------------------------------------------------------
  
  
  
  old_annotated_now_novel <- intersect(annotated_old_junID,
                                       c(novel_new_junID, other_new_junID))
  
  
  old_annotated_now_novel_recount2 <- all_recount2_split_reads %>%
    filter(junID %in% old_annotated_now_novel)
  
  
  gtex_samples <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_used.rda")
  samples <- gtex_samples$`Brain-FrontalCortex_BA9`
  
  ## Get the count data in the current cluster for the current NOVEL junction
  split_read_counts_novel <- split_read_counts %>%
    filter(junID %in% old_annotated_now_novel_recount2$junID) %>%
    dplyr::select(junID, samples %>% as.character()) %>%
    #mutate_all(~replace(., is.na(.), 0)) %>%
    as_tibble()
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  
  split_read_counts_novel[,"novel_n_individuals"] <- matrixStats::rowCounts(split_read_counts_novel[, -c(1)] > 0, na.rm = T)
  split_read_counts_novel[,"novel_mean_counts"] <- rowMeans(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(), 1)], na.rm = T)
  split_read_counts_novel[,"novel_sum_counts"] <- rowSums(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(),
                                                                                      (split_read_counts_novel %>% ncol()) - 1,
                                                                                      1)], na.rm = T)
  
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  if (any(split_read_counts_novel[, "novel_n_individuals"] < 1)) {
    print("Error: some novel junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_novel <- split_read_counts_novel[, c(1, 
                                                         (split_read_counts_novel %>% ncol() - 2),
                                                         (split_read_counts_novel %>% ncol() - 1), 
                                                         split_read_counts_novel %>% ncol())]
  split_read_counts_novel %>% head()
  
  
  
  df_old_annotated_now_novel <- merge(x = old_annotated_now_novel_recount2,
                                      y = split_read_counts_novel,
                                      by.x = "junID", 
                                      by.y = "junID",
                                      all.x = T)
  
  
  
  
  ggplot(data = df_old_annotated_now_novel) +
    geom_point(mapping = aes(x = novel_n_individuals,
                             y = novel_mean_counts)) + 
    geom_abline(mapping=aes(slope=0, intercept=5, linetype=factor(3), color ="red")) +
    scale_y_continuous(breaks = c(0,5,10,15,20,50,75,100,130)) +
    theme_light() +
    
    ylab("mean number of reads") +
    xlab("number of individuals") +
    ggtitle(paste0("Annotated junctions ensembl 2019 --> novel/combo/exon skip in 2021.\nFrontal Cortex - Junctions represented: ", 
                   df_old_annotated_now_novel %>% distinct(junID) %>% nrow())) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "none")
  
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/Brain-FrontalCortex_BA9/images/annotated2019_novel2021.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  ## GETTING FROM MY DATABASE OF INTRONS THE NOVEL JUNCTIONS IN 2019 THAT ARE ANNOTATED IN 2021 ---------------------------------------------------------------
  
  
  old_novel_now_annotated <- intersect(novel_old_junID,
                                       c(annotated_new_junID, other_new_junID))
  
  
  split_reads_details_now_annotated <- all_split_reads_details_104 %>%
    filter(junID %in% old_novel_now_annotated)
  
  
  split_reads_details_now_annotated$type %>% unique
  split_reads_details_now_annotated$seqnames %>% unique
  
  
  ## Add counts and number of individuals
  
  
  gtex_samples <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_used.rda")
  samples <- gtex_samples$`Brain-FrontalCortex_BA9`
  
  ## Get the count data in the current cluster for the current NOVEL junction
  split_read_counts_novel <- split_read_counts %>%
    filter(junID %in% split_reads_details_now_annotated$junID) %>%
    dplyr::select(junID, samples %>% as.character()) %>%
    #mutate_all(~replace(., is.na(.), 0)) %>%
    as_tibble()
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  
  split_read_counts_novel[,"novel_n_individuals"] <- matrixStats::rowCounts(split_read_counts_novel[, -c(1)] > 0, na.rm = T)
  split_read_counts_novel[,"novel_mean_counts"] <- rowMeans(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(), 1)], na.rm = T)
  split_read_counts_novel[,"novel_sum_counts"] <- rowSums(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(),
                                                                                      (split_read_counts_novel %>% ncol()) - 1,
                                                                                      1)], na.rm = T)
  
  split_read_counts_novel %>% as.data.frame() %>% head()
  
  
  if (any(split_read_counts_novel[, "novel_n_individuals"] < 1)) {
    print("Error: some novel junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_novel <- split_read_counts_novel[, c(1, 
                                                         (split_read_counts_novel %>% ncol() - 2),
                                                         (split_read_counts_novel %>% ncol() - 1), 
                                                         split_read_counts_novel %>% ncol())]
  split_read_counts_novel %>% head()
  
  
  
  df <- merge(x = split_reads_details_now_annotated,
              y = split_read_counts_novel,
              by.x = "junID", 
              by.y = "junID",
              all.x = T)
  
  df %>%
    filter(gene_id_start %>% as.character() != gene_id_end %>% as.character())
  
  df %>%
    filter(novel_n_individuals >= 54, novel_mean_counts >= 5) %>%
    filter(gene_id_start %>% as.character() == gene_id_end %>% as.character())
  
  
  
  ## LOAD SOURCE DATA
  file_name <- paste0(folder_root, "/pipeline3/missplicing-ratio/", cluster, "/", cluster, "_db_introns.rds")
  df_introns <- readRDS(file = file_name)
  
  df_introns %>%
    filter(ref_junID %in% old_novel_now_annotated)
  
  
  ggplot(data = df_introns %>%
           filter(ref_junID %in% old_novel_now_annotated)) +
    geom_point(mapping = aes(x = ref_n_individuals,
                             y = ref_mean_counts)) + 
    geom_abline(mapping=aes(slope=0, intercept=0, linetype=factor(3), color ="red")) +
    #scale_color_manual(values = c("red" = "red", 
    #                              "black" = "black"),
    #                   labels = c(paste0("individuals < 96 & mean_counts < 5 (",df_novel_new %>% filter(color == "black") %>% nrow,")"),  
    #                              paste0("individuals >= 96 & mean_counts >= 5 (", df_novel_new %>% filter(color == "red") %>% nrow,")"))) +
    theme_light() +
    xlim(c(0,108)) +
    ylim(c(0,130)) +
    ylab("mean number of reads") +
    xlab("number of individuals") +
    ggtitle(paste0("Annotated junctions ensembl 2019 --> novel/combo/exon skip in 2021.\nFrontal Cortex - Junctions represented: ", 
                   df_old_annotated_now_novel %>% distinct(junID) %>% nrow())) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "none")
}










##########################
## PLOTS #################
##########################


## DISTANCES

plot_distances <- function(cluster = "Brain - Frontal Cortex (BA9)",
                           folder_name,
                           limit_bp = 30,
                           save_plot = F) {
  
  # cluster <- "Brain-FrontalCortex_BA9"
  # folder_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/"
  # folder_name <- paste0(folder_name, cluster, "/v104/")
  # print(stringr::str_c(Sys.time(), " - Plotting distances for '", cluster, "' samples..."))
  
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  if (file.exists(paste0(folder_root, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_root, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_root, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
  }
  
  
 
  
  #############################
  ## PLOT THE DISTANCES GRAPH 
  #############################
  
  
  df_novel_tidy <- df_novel %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  
  #title <- paste0("Distances - ", cluster)
  

  #################################
  ## GENERATE PLOT ################
  #################################

  plot <- ggplot(data = df_novel_tidy) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggforce::facet_col(vars(novel_type)) +
    #ggtitle(paste0(title)) +
    #ylim(y_axes) +
    xlab("Distance to the reference intron (in bp)") +
    ylab("Number of unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("novel donor", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
                               override.aes = list(size = 3),
                               ncol = 3 )) +
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
  
  if (save_plot) {
    plot
    file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/distances.png"
    ggplot2::ggsave(filename = file_name,
                    width = 183, height = 183, units = "mm", dpi = 300)
  } else {
    return(plot)
  }
 
  
  
  
}

plot_distances_PC <- function(cluster = "Brain - Frontal Cortex (BA9)",
                              folder_name,
                              limit_bp = 30,
                              stats = F,
                              save_results = F) {
  
  # cluster <- "Brain-FrontalCortex_BA9"
  # folder_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/"
  # folder_name <- paste0(folder_name, cluster, "/v104/")
  # print(stringr::str_c(Sys.time(), " - Plotting distances for '", cluster, "' samples..."))
  
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  
  if (file.exists(paste0(folder_root, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_root, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_root, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
  }
  
  

  
  #############################
  ## PLOT THE DISTANCES GRAPH 
  #############################
  
  
  df_novel_tidy <- df_novel %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  df_novel_tidy$type_PC = factor(df_novel_tidy$type_PC, 
                                 levels = c("protein coding (PC)", "non PC"))
  
  # title <- paste0("Distances - ", cluster)
  

  #################################
  ## GENERATE PLOT ################
  #################################

  
  
  plot_PC <- ggplot(data = df_novel_tidy %>% filter(type_PC == "protein coding (PC)")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggforce::facet_col(vars(novel_type)) +
    ggtitle(paste0("Protein coding (PC)")) +
    #ylim(y_axes) +
    xlab("Distance (bp) to the reference intron)") +
    ylab("Number of unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("novel donor", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
                               override.aes = list(size = 3),
                               ncol = 3 )) +
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
  
  plot_NPC <- ggplot(data = df_novel_tidy %>% filter(type_PC == "non PC")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggforce::facet_col(vars(novel_type)) +
    ggtitle(paste0("Non-protein coding (NPC)\n")) +
    #ylim(y_axes) +
    xlab("Distance (bp) to the reference intron") +
    ylab("Number of unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = seq(from = (limit_bp * -1), to = limit_bp, by = 6)) +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("novel donor", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
                               override.aes = list(size = 3),
                               ncol = 3 )) +
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
  
  
  plot <- ggpubr::ggarrange(plot_PC, 
                    plot_NPC,
                    common.legend = T,
                    labels = c("A","B"),
                    align = "hv",
                    ncol = 2, 
                    nrow = 1)
  
  #plot <- annotate_figure(plot, top = text_grob(title, face = "bold", size = 14))
  
  if (save_results) {
    plot
    file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/distancespc.png"
    ggplot2::ggsave(file_name, width = 183, height = 143, units = "mm", dpi = 300)
  } else {
    return(plot)
  }
  
  
  ## RELEASE SOME MEMORY
  rm(df_tidy)
  rm(title)
  rm(file_name)
  rm(y_axes)
}

## MIS-SPLICING RATIO
# cluster <- gtex_tissues[11]
# folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", cluster, "/v104/")

plot_missplicing_ratio <- function(cluster,
                                   folder_name,
                                   #PC = NULL,
                                   #intron_type = NULL,
                                   save_results = F)  {
  
  
  
  df <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>% as_tibble()
  df %>% head()
  df %>% nrow()
  
  
  ##########################
  ## PLOT MIS-SPLICING RATIO 
  ##########################
  
  y_axes_values <- density(x = df %>% pull(ref_missplicing_ratio_tissue_ND))$y
  y_axes <- c(0, y_axes_values %>% max() + 
                y_axes_values %>% max() / 4)
  
  plot <- ggplot(data = df) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"),
                 alpha = 0.8) +
    geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
                 alpha = 0.8) +
    
    ggtitle(paste0("Mis-splicing ratio (MSR) - ", cluster)) +
    xlab("Mis-splicing ratio") +
    #ylim(y_axes) +
    ggforce::facet_zoom(xlim = c(0,0.15)) +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text.x =  element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  
    return(plot)
  
  
  
  
}

plot_missplicing_ratio_PC <- function(cluster,
                                      folder_name,
                                      save_results = F)  {
  
  
  
  df <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>% as_tibble()
  df %>% head()
  df %>% nrow()
  
  df_tidy <- df %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) 
  
  ##########################
  ## PLOT MIS-SPLICING RATIO 
  ##########################
  
  y_axes_values <- density(x = df_tidy %>% filter(protein_coding == 100) %>% pull(ref_missplicing_ratio_tissue_ND))$y
  y_axes <- c(0, y_axes_values %>% max() + 
                y_axes_values %>% max() / 4)
  
  plot_PC <- ggplot(data = df_tidy %>% filter(protein_coding == 100)) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"),
                 alpha = 0.8) +
    geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
                 alpha = 0.8) +
    
    ggtitle(paste0("Protein coding (PC)\n")) +
    xlab("Mis-splicing ratio") +
    #ylim(y_axes) +
    ggforce::facet_zoom(xlim = c(0,0.15)) +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text.x =  element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  plot_NPC <- ggplot(data = df_tidy %>% filter(protein_coding == 0)) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"),
                 alpha = 0.8) +
    geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
                 alpha = 0.8) +
    ggtitle(paste0("Non-protein coding (NPC)\n")) +
    xlab("Mis-splicing ratio") +
    ggforce::facet_zoom(xlim = c(0,0.15)) +
    #ylim(y_axes) +
    ylab("") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  plot <- ggarrange(plot_PC, 
                    plot_NPC,
                    common.legend = T,
                    labels = c("A","B"),
                    align = "hv",
                    ncol = 2, 
                    nrow = 1)
  title <- paste0("Mis-splicing ratio (MSR) - ", cluster)
  plot <- annotate_figure(plot, top = text_grob(title, face = "bold", size = 14))
  
  if (save_results) {
    plot
    folder_img <- paste0(folder_name, "/images/")
    dir.create(file.path(folder_img), showWarnings = F)
    file_name <- paste0(folder_img, cluster, "_missplicingratio.png")
    ## Save the results
    ggplot2::ggsave(filename = file_name, 
                    width = 183, height = 143, units = "mm", dpi = 300)
    print(paste0(Sys.time(), " - density plot saved!"))
  } else {
    return(plot)
  }
  
  
  
}





## MODULO

plot_distances_modulo <- function(cluster,
                                     folder_name,
                                     #PC = NULL,
                                     save_results = F) {
  
  cluster <- "Brain - Frontal Cortex (BA9)"
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/results/pipeline3/missplicing-ratio/",
                        cluster, "/v105/")
  
  if (file.exists(paste0(folder_root, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_root, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_root, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
  }
  ## Load novel junction database for current tissue
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  df_novel %>% head()
  
  
  df_novel_tidy <- df_novel %>%
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

  plot <- ggplot(data = df_novel_tidy, 
                 aes(x = factor(module), y = perc*100, fill = factor(type_p))) + 
    geom_bar(stat = "identity", position = "dodge") +
    ggtitle(paste0(title)) +
    # geom_text(aes( label = scales::percent(..prop.., accuracy = 0.1),
    #                y= ..prop.. ), stat= "count", vjust = 2, colour = "red") +
    #scale_x_continuous(breaks = c(0, 1, 2)) +
    #scale_y_continuous(labels = scales::percent) +
    scale_fill_viridis_d() +
    ylab("% of novel junctions") +
    xlab("modulo 3") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(colour = "black",size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          strip.text.y = element_blank(),
          plot.caption = element_text(colour = "black",size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 2))
  
  return(plot)
  
  
  
  
}


plot_distances_modulo_PC <- function(cluster,
                                  folder_name,
                                  #PC = NULL,
                                  save_results = F) {
  
  
  ## Load novel junction database for current tissue
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  df_novel %>% head()
  
  
  df_novel_tidy <- df_novel %>%
    filter(protein_coding %in% c(0, 100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " ")) %>%
    mutate(type_p = ifelse(distance < 0, paste0(novel_type," intron"), paste0(novel_type," exon"))) %>% 
    mutate(module = abs(distance) %% 3)
  
  df_novel_tidy$novel_type = factor(df_novel_tidy$novel_type, 
                                    levels = c("novel donor", "novel acceptor"))
  df_novel_tidy$type_PC = factor(df_novel_tidy$type_PC, 
                                 levels = c("protein coding (PC)", "non PC"))
  df_novel_tidy$type_p = factor(df_novel_tidy$type_p, 
                                levels = c("novel acceptor intron", 
                                           "novel acceptor exon",
                                           "novel donor intron", 
                                           "novel donor exon"))
  title <- paste0("Modulo 3 - ", cluster)
  
  
  plot <- ggplot(data = df_novel_tidy, 
                 aes(x = module, group = type_p, fill = type_p)) + 
    geom_bar(aes(y = ..prop..), stat="count") +
    ggtitle(paste0(title)) +
    geom_text(aes( label = scales::percent(..prop.., accuracy = 0.1),
                   y= ..prop.. ), stat= "count", vjust = 2, colour = "red") +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_viridis_d() +
    facet_grid(type_p~type_PC) + 
    ylab("% of novel junctions") +
    xlab("modulo 3") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(colour = "black",size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          strip.text.y = element_blank(),
          plot.caption = element_text(colour = "black",size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 4, 
                               nrow = 1))
  
  if (save_results) {
    plot
    folder_img <- paste0(folder_name, "/images/")
    dir.create(file.path(folder_img), showWarnings = T)
    file_name <- paste0(folder_img, cluster, "_distances_24bp_module.png")
    
    ## SAVE THE RESULTS
    ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  } else {
    return(plot)
  }
  
  
  
}


plot_delta_maxentscanscore <- function(cluster,
                                       folder_name,
                                       save_results = F) {
  
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))  %>% 
    as_tibble() %>%
    filter(u2_intron == T)
  
  df_introns <- df_introns %>%
    as.data.frame() %>%
    dplyr::select(ref_junID, 
                  ref_ss5score, ref_ss3score, 
                  ref_type)
  df_introns %>% head()
  df_introns %>% nrow()

  
  
  ## Load the novel junction database ---------------------------------------------
  
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  df_novel <- df_novel %>%
    dplyr::select(ref_junID, 
                  novel_junID, 
                  novel_ss5score, novel_ss3score, 
                  novel_type, gene_id)
  df_novel %>% head()
  df_novel %>% nrow()
  
  
  ## Merge both databases ---------------------------------------------------------
  
  df_merged <- merge(x = df_novel %>% data.table::as.data.table(),
                     y = df_introns %>% data.table::as.data.table(),
                     by = "ref_junID",
                     all.x = T)
  
  df_merged %>% nrow()
  df_merged %>% head()
  
  df_merged <- df_merged %>%
    mutate(diff_ss5score = ref_ss5score - novel_ss5score,
           diff_ss3score = ref_ss3score - novel_ss3score)
  
  
  ## Gather
  df_all_tidy <- df_merged %>% 
    dplyr::select(diff_ss5score, diff_ss3score, novel_type, gene_id) %>%
    tidyr::gather(key = "type", value = "mean", -novel_type, -gene_id) %>%
    dplyr::filter(mean != 0)
  
  df_all_tidy <- df_all_tidy %>%
    mutate(type = type %>% as.factor()) %>%
    mutate(type = relevel(type, ref = "diff_ss5score"))
  
  title <- paste0("Delta MaxEntScanScore - ", cluster)
  
  ## Plot
  plot <- ggplot(data = df_all_tidy) +
    geom_density(mapping = aes(x = mean, fill = type), 
                 alpha = 0.8) +
    ggtitle(paste0(title)) +
    geom_vline(xintercept = 0) +
    xlab("Delta MaxEntScan score") +
    theme_light() +
    scale_color_viridis_d() +
    scale_fill_manual(values =  c("#35B779FF","#440154FF"),
                      breaks = c("diff_ss5score", "diff_ss3score"),
                      labels = c("Delta MES 5'ss (novel donor)", "Delta MES 3'ss (novel acceptor)")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text.x = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  return(plot)

}

plot_delta_maxentscanscore_PC <- function(cluster,
                                          folder_name,
                                          save_results = F) {
  
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))  %>% 
    as_tibble() %>%
    filter(u2_intron == T)
  
  df_introns <- df_introns %>%
    as.data.frame() %>%
    dplyr::select(ref_junID, 
                  ref_ss5score, ref_ss3score, 
                  protein_coding, ref_type)
  df_introns %>% head()
  df_introns %>% nrow()
  
  df_introns <- df_introns %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) 
  
  df_introns$type_PC = factor(df_introns$type_PC, 
                              levels = c("protein coding (PC)", "non PC"))
  
  
  
  ## Load the novel junction database ---------------------------------------------
  
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  df_novel <- df_novel %>%
    dplyr::select(ref_junID, 
                  novel_junID, 
                  novel_ss5score, novel_ss3score, 
                  novel_type, gene_id)
  df_novel %>% head()
  df_novel %>% nrow()
  
  
  ## Merge both databases ---------------------------------------------------------
  
  df_merged <- merge(x = df_novel,
                     y = df_introns,
                     by = "ref_junID",
                     all.x = T)
  
  df_merged %>% nrow()
  df_merged %>% head()
  
  df_merged <- df_merged %>%
    mutate(diff_ss5score = ref_ss5score - novel_ss5score,
           diff_ss3score = ref_ss3score - novel_ss3score)
  
  
  ## Gather
  df_all_tidy <- df_merged %>% 
    dplyr::select(diff_ss5score, diff_ss3score, novel_type, gene_id,protein_coding,type_PC) %>%
    tidyr::gather(key = "type", value = "mean", -novel_type, -gene_id, -protein_coding, -type_PC ) %>%
    dplyr::filter(mean != 0)
  
  df_all_tidy <- df_all_tidy %>%
    mutate(type = type %>% as.factor()) %>%
    mutate(type = relevel(type, ref = "diff_ss5score"))
  
  title <- paste0("Delta MaxEntScanScore - ", cluster)
  
  ## Plot
  plot <- ggplot(data = df_all_tidy) +
    geom_density(mapping = aes(x = mean, fill = type), 
                 alpha = 0.8) +
    ggtitle(paste0(title)) +
    geom_vline(xintercept = 0) +
    xlab("Delta MaxEntScan score") +
    theme_light() +
    scale_color_viridis_d() +
    facet_grid(~type_PC) +
    scale_fill_manual(values =  c("#35B779FF","#440154FF"),
                      breaks = c("diff_ss5score", "diff_ss3score"),
                      labels = c("Delta MES 5'ss (novel donor)", "Delta MES 3'ss (novel acceptor)")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text.x = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
  
  if (save_results) {
    plot
    file_name <- "/home/sruiz/PROJECTS/splicing-project/results/paper/supplementary_figure51.png"
    ggplot2::ggsave(filename = file_name,
                    width = 183, height = 183, units = "mm", dpi = 300)
  } else {
    return(plot)
  }
  
}


plot_cons_scores_tissue <- function(cluster,
                                    folder_name,
                                    save_results = F) {
  
  
  db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))   %>% 
    as_tibble() %>%
    filter(u2_intron == T)
  
  
  #################################
  ## CONSERVATION
  #################################
  
  db_introns_cons <- db_introns %>%
    dplyr::select(cons_5ss_mean = phastCons20way_5ss_mean, 
                  cons_3ss_mean = phastCons20way_3ss_mean, 
           #protein_coding,
           #cons_5ss_median = phastCons20way_5ss_median,
           #cons_3ss_median = phastCons20way_3ss_median,
           ref_type) %>%
    gather(key = cons_type, value = cons_value, -ref_type) %>%
    mutate(cons_type = cons_type %>% as.factor())
  
  
  db_introns_cons$cons_type = factor(db_introns_cons$cons_type, 
                                     levels = c("cons_5ss_mean", "cons_3ss_mean"))
  
  title <- paste0("Conservation score - ", cluster)
  plot <- ggplot(db_introns_cons) +
    geom_density(aes(x = cons_value, 
                     fill = ref_type), alpha = 0.6) +
    ggtitle(paste0(title)) +
    scale_fill_viridis_d() +
    theme_light() +
    xlab("Conservation score") +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Intron type: ",
                               ncol = 4, 
                               nrow = 1))
  
  plot %>% return()

}


plot_cons_scores_tissue_PC <- function(cluster,
                                       folder_name,
                                       save_results = F) {
  
  
  db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))   %>% 
    as_tibble() %>%
    filter(u2_intron == T)
  
  
  #################################
  ## CONSERVATION
  #################################
  
  db_introns_cons <- db_introns %>%
    dplyr::select(cons_5ss_mean = phastCons20way_5ss_mean, 
                  cons_3ss_mean = phastCons20way_3ss_mean, 
                  protein_coding,
           #cons_5ss_median = phastCons20way_5ss_median,
           #cons_3ss_median = phastCons20way_3ss_median,
           ref_type) %>%
    gather(key = cons_type, value = cons_value, -ref_type,-protein_coding) %>%
    mutate(cons_type = cons_type %>% as.factor()) %>%
    #filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) 
  
  
  db_introns_cons$cons_type = factor(db_introns_cons$cons_type, 
                                     levels = c("cons_5ss_mean", "cons_3ss_mean"))
  db_introns_cons$type_PC = factor(db_introns_cons$type_PC, 
                                   levels = c("protein coding (PC)", "non PC"))
  
  title <- paste0("Conservation score - ", cluster)
  plot <- ggplot(db_introns_cons) +
    geom_density(aes(x = cons_value, 
                     fill = ref_type), alpha = 0.6) +
    facet_grid(cons_type~type_PC) + 
    ggtitle(paste0(title)) +
    scale_fill_viridis_d() +
    theme_light() +
    xlab("Conservation score") +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Intron type: ",
                               ncol = 4, 
                               nrow = 1))
  
  
  if (save_results) {
    plot
    file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/conservation.png"
    ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  } else {
    plot %>% return()
  }
  
  
  
  
}

plot_CDTS_scores_tissue <- function(cluster,
                                    folder_name,
                                    save_results = F) {
  
  #folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", cluster, "/v104/")
  db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))   %>% 
    as_tibble() %>%
    filter(u2_intron == T)
  
  #################################
  ## CONSTRAINT
  #################################
  
  db_introns_CDTS <- db_introns %>%
    dplyr::select(CDTS_5ss_mean, 
           #CDTS_5ss_median, 
           CDTS_3ss_mean, 
           #CDTS_3ss_median,
           protein_coding,
           ref_type) %>%
    gather(key = CDTS_type, value = CDTS_value, -ref_type, -protein_coding) %>%
    mutate(CDTS_type = CDTS_type %>% as.factor()) %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(type_PC = ifelse(protein_coding == 100, "protein coding (PC)", "non PC")) 
  
  
  db_introns_CDTS$CDTS_type = factor(db_introns_CDTS$CDTS_type, 
                                     levels = c("CDTS_5ss_mean", "CDTS_3ss_mean"))
  db_introns_CDTS$type_PC = factor(db_introns_CDTS$type_PC, 
                                   levels = c("protein coding (PC)", "non PC"))
  
  title <- paste0("CDTS score - ", cluster)
  
  
  plot <- ggplot(db_introns_CDTS) +
    geom_density(aes(x = CDTS_value, 
                     fill = ref_type), alpha = 0.6) +
    facet_grid(CDTS_type~type_PC)+ 
    theme_light() +
    scale_fill_viridis_d() +
    ggtitle(paste0(title)) +
    xlab("CDTS score") +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Intron type: ",
                               ncol = 4, 
                               nrow = 1))
  
  
  if (save_results) {
    plot
    file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/CDTS.png"
    ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  } else {
    plot %>% return()
  }
}

plot_lm_single_tissue <- function(cluster,
                                  folder_name,
                                  GTEx = F,
                                  save_results = F) {
  
  #########################
  ## LOAD AND TIDY DATA
  #########################
  
  idb <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  
  # idb <- idb %>%
  #   as.data.frame() %>%
  #   dplyr::distinct(ref_junID, .keep_all = T) %>%
  #   dplyr::rename(intron_length = width,
  #                 intron_5ss_score = ref_ss5score,
  #                 intron_3ss_score = ref_ss3score,
  #                 gene_length = gene_width,
  #                 gene_tpm = tpm_median_rct3,
  #                 gene_num_transcripts = n_transcripts,
  #                 CDTS_5ss = CDTS_5ss_mean,
  #                 CDTS_3ss = CDTS_3ss_mean,
  #                 mean_phastCons20way_5ss = phastCons20way_5ss_mean,
  #                 mean_phastCons20way_3ss = phastCons20way_3ss_mean)
  idb <- idb %>%
    as.data.frame() %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    dplyr::rename(intron_length = width,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score,
                  gene_length = gene_width,
                  gene_tpm = tpm_median_rct3,
                  gene_num_transcripts = n_transcripts,
                  CDTS_5ss = CDTS_5ss_mean,
                  CDTS_3ss = CDTS_3ss_mean,
                  protein_coding = protein_coding,
                  mean_phastCons20way_5ss = phastCons20way_5ss_mean,
                  mean_phastCons20way_3ss = phastCons20way_3ss_mean)
  
  #########################
  ## LINEAR MODELS
  #########################
  

  fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ 
                    intron_length +
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
  
  fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ 
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
  
  plot <- jtools::plot_summs(fit_donor, 
                             fit_acceptor,
                             #scale = TRUE, 
                             robust = T,
                             #inner_ci_level = .75,
                             #n.sd = 2,
                             pvals = TRUE,
                             legend.title = "Model:",
                             #plot.distributions = TRUE,
                             ci_level = 0.95,
                             #coefs = coef_names,
                             colors = c("#35B779FF","#440154FF"),
                             model.names = model_names) + 
    theme_minimal() + 
    theme(axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = "black", size = "13"),
          axis.text.y = element_text(colour = "black", size = "13"),
          axis.title = element_text(colour = "black", size = "13"),
          legend.text = element_text(colour = "black", size = "13"),
          legend.title = element_text(colour = "black", size = "13"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 2)) +  
    ylab("Predictors") +
    geom_hline(yintercept = seq(from = 0.5,
                                to = length((fit_donor$coefficients %>% names)[-1]) + .5,
                                by = 1)) +
    guides(colour = guide_legend(ncol = 2, 
                                 nrow = 1))
  
  
  if (save_results) {
    plot
    filename <- paste0(folder_name, "/images/lm_MSR.png")
    ggplot2::ggsave(filename = filename, 
                    width = 183, height = 183, units = "mm", dpi = 300)
  } else {
    plot %>% return()
  }
  
}


plot_modulo_tissues <- function(type) {
  
  if (type == "exon") {
    df_modulo_prop <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/modulo_exon_all_tissues.rds")
  } else {
    df_modulo_prop <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/modulo_intron_all_tissues.rds")
  }
  
  
  df_modulo_prop$novel_type = factor(df_modulo_prop$novel_type, 
                                     levels = c( "donor","acceptor"))
  
  df_modulo_prop <- df_modulo_prop %>% 
    ungroup() %>%
    arrange(type , desc(modulo)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df_modulo_prop$tissue), pattern = "Brain"), "red", "black")
  
  ggplot(data = df_modulo_prop) +
    geom_bar(mapping = aes(x = tissue, y = modulo, fill = factor(type)), 
             stat = "identity", position = position_dodge()) +
    xlab("") +
    ylab("Percentage of junctions (%)") +
    theme_light() +
    facet_col(~novel_type) +
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
          axis.text.x = element_text(color = colours,
                                     angle = 70, 
                                     vjust = 1,
                                     hjust = 1),
          legend.text = element_text(size = "14"),
          legend.title = element_text(size = "14"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) %>% return()
}


plot_estimate_variance <- function(project_id = "BRAIN",
                                   gtf_version = 105,
                                   brain_tissues = F,
                                   save_results = F) {
  
  if (brain_tissues) {
    df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                         "/results/pipeline3/variance_estimate.rds"))
    graph_title <- "Distribution of the estimate values across 11 brain tissues"
  } else {
    df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/variance_estimate.rds"))
    graph_title <- "Distribution of the estimate values across 54 GTEx tissues"
  }
  
  
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
    gather(tissue, feature, -type, -tissue) %>%
    mutate(type = "MSR_Donor")
  
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
  
  plot <- ggplot(MSR_tidy, aes(tissue, feature)) + 
    geom_boxplot() +
    facet_grid(vars(type)) +
    ggtitle(graph_title) +
    ylab("Distribution of the estimate") +
    xlab("Predictor") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    theme(axis.text.x = element_text(angle = 70, 
                                     vjust = 1,
                                     hjust = 1)) +
    geom_hline(yintercept = 0,linetype='dotted')
  
  
  if (save_results) {
    plot
    
    ## SAVE THE RESULTS
    if (brain_tissues) {
      file_name <- "/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_brain_tissues.png"
    } else {
      file_name <- "/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_all_tissues.png"
    }
    
    ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
    
  } else {
    plot %>% return()
  }
}


##########################
## SUPPLEMENTARY        ##
##########################



get_TPM_analysis <- function(cluster,
                             folder_name) {
  
  ## LOAD TPM DATA FOR THE CURRENT TISSUE
  
  library(GSRI)
  
  gtex_tpm <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_gene_median_tpm_tidy.rds")
  
  gtex_tpm <- gtex_tpm %>% 
    dplyr::select(gene_id, cluster) %>%
    distinct(gene_id, .keep_all = T) 
  
  gtex_tpm %>% head()
  gtex_tpm %>% nrow()
  
  
  ## LOAD MIS-SPLICING RATIO DATA FOR THE CURRENT TISSUE
  
  df <- readRDS(file = paste0(folder_name, "/", cluster, "_missplicingratio_tidy_v2.rds"))
  df %>% head()
  df %>% 
    filter(type == "acceptor") %>% 
    dplyr::select(c(gene_id_start, gene_id_end)) %>%
    unname() %>% unlist() %>% unique() %>%
    length()
  
  genes <- df %>%
    dplyr::select(c(gene_id_start, gene_id_end)) %>%
    unname() %>% unlist() %>% unique()
  genes %>% length()
  
  ## GET TPM VALUES FOR GENES MIS-SPLICED IN THE CURRENT TISSUE
  
  fctx_tpm <- gtex_tpm %>%
    filter(gene_id %in% genes)
  
  fctx_tpm %>% nrow()  
  genes %>% length()
  
  
  ## PLOT TPM VALUES
  plot(density(gtex_tpm$tpm), xlab = "Median TPM values", main = "Median TPM values.\nAll FCTX genes vs. filtered FCTX genes.")
  lines(density(fctx_tpm$tpm), col = "red")
  
  legend(y = 2, x = 40000, 
         legend=c("All FCTX genes", "filtered FCTX genes"),
         col=c("black", "red"), lty=1, cex=0.8)
  
  
  plot(density(fctx_tpm$tpm), xlab = "Median TPM values", main = "Median TPM values.\nFiltered FCTX genes only.", col = "red")
  
  legend(y = 0.03, x = 1500, 
         legend=c("Filtered FCTX genes"),
         col=c("red"), lty=1, cex=0.8)
  
  
  
  ## General stats
  
  gtex_tpm$tpm %>% summary(digits = 6)
  gtex_tpm$tpm %>% get_mode()
  
  ## FCTX stats
  fctx_tpm$tpm %>% summary(digits = 8)
  fctx_tpm$tpm %>% get_mode()
  
}

##############################################################################

get_novel_annotation_numbers <- function(tissue) {
  
  all_people_tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_people_used_tissue.rda")
  
  folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue, "/v104/")
  people_tissue <- all_people_tissues[[tissue]] %>% length()
  #df_intron <- readRDS(file = paste0(folder_name, "/", tissue, "_db_introns.rds")) %>% as_tibble()
  df_novel <- readRDS(file = paste0(folder_name, "/", tissue, "_db_novel.rds")) %>% 
    as_tibble() %>%
    mutate(novel_p_individuals = (novel_n_individuals * 100) / people_tissue)
  
  
  print(paste0("Number of novel junctions shared by <= 95% of samples: ", df_novel %>%
                 filter(novel_p_individuals <= 95) %>%
                 nrow()))
  
  
  print(paste0("Number of novel junctions shared by > 95% of samples: ", df_novel %>%
                 filter(novel_p_individuals > 95) %>%
                 nrow()))
  
}

get_novel_annotation_plot <- function(tissue) {
  
  all_people_tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_people_used_tissue.rda")
  
  folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue, "/v104/")
  people_tissue <- all_people_tissues[[tissue]] %>% length()
  
  #df_intron <- readRDS(file = paste0(folder_name, "/", tissue, "_db_introns.rds")) %>% as_tibble()
  df_novel <- readRDS(file = paste0(folder_name, "/", tissue, "_db_novel.rds")) %>% 
    as_tibble() %>%
    mutate(novel_p_individuals = (novel_n_individuals * 100) / people_tissue)
  
  ggplot(data = df_novel) +
    geom_point(aes(x = novel_p_individuals, 
                   y = novel_mean_counts)) +
    ylab("mean number of reads\n(across individuals)") + 
    xlab("% of individuals") +
    geom_vline(xintercept = 95, 
               linetype = "dashed", 
               color = "red") + 
    scale_x_continuous(breaks = c(0, 25, 50, 75, 95, 100))    %>%
    return()
}
