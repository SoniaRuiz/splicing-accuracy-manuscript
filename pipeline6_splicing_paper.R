library(tidyverse)
library(GenomicRanges)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)

####################################################
## PAPER ###########################################
####################################################


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
    arrange(type , prop) %>%
    mutate(tissue = fct_inorder(tissue))
  colours <- ifelse(str_detect(string = as.factor(df2$tissue), pattern = "Brain"), "red", "black")
  
  
  ggplot(df2, aes(x = tissue, 
                  y = prop, 
                  group = type, 
                  fill = type)) +
    geom_col(position = position_dodge()) +
    #geom_text(aes(label = round(x = prop, digits = 2)), 
    #          position = position_stack(vjust = 0.5, reverse = TRUE), colour = "white", size = 1.5) +
    coord_flip() +
    theme_light() +
    ylab("unique novel junctions") +
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
                      breaks = c( "acceptor", "donor" ),
                      labels = c( "novel acceptor", "novel donor")) +
    guides(fill = guide_legend(title = NULL,
                               ncol = 1, 
                               nrow = 2))
  
  
  ## Save the figure 3
  file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/percent_unique_donor_acceptor.svg"
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
    arrange(type , desc(mean)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df_mean_counts_tidy$tissue), pattern = "Brain"), "red", "black")
  
  
  ggplot(data = df_mean_counts_tidy) +
    geom_col(mapping = aes(x = tissue, y = mean, fill = type), 
             #position = position_dodge(width = 0.6), 
             alpha = 0.9) +
    coord_flip() +
    xlab("") +
    ylab("percent of avg. read count across samples (%)") +
    theme_light() +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = c(0,10,25,50,75,100)) +
    #coord_cartesian(ylim=c(1.1,50)) +
    scale_fill_manual(values =  c("#0D0887FF","#EF7F4FFF","#A41F9AFF"),#c( "#440154FF", "#35B779FF", "#39558CFF"),
                      breaks = c("annotated", "acceptor", "donor"),
                      labels = c("annotated", "novel acceptor", "novel donor")) +
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



## SECTION 2 ---------------------------------------------

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



## SECTION 3 ---------------------------------------------

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
    mutate(type_PC = ifelse(protein_coding == 100, "PC", "non PC")) 
  
  df_tidy$type_PC = factor(df_tidy$type_PC, 
                                    levels = c("PC", "non PC"))
  
  
  df_tidy <- df_tidy %>%
    dplyr::select(type_PC, 
                  MSR_D = ref_missplicing_ratio_tissue_ND,
                  MSR_A = ref_missplicing_ratio_tissue_NA) %>%
    gather(key = "MSR_type", value = "MSR", -type_PC)
  
  ##########################
  ## PLOT MIS-SPLICING RATIO 
  ##########################
  
  #y_axes_values <- density(x = df_tidy %>% filter(protein_coding == 100) %>% pull(ref_missplicing_ratio_tissue_ND))$y
  #y_axes <- c(0, y_axes_values %>% max() + 
  #              y_axes_values %>% max() / 4)
  
  plotMSR <- ggplot(data = df_tidy) + 
    geom_density(aes(x = MSR, fill = MSR_type), alpha = 0.8) +
    facet_grid(vars(type_PC)) +
    xlab("Mis-splicing ratio") +
    #ylim(y_axes) +
    #ggforce::facet_zoom(xlim = c(0,0.15)) +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("MSR_D","MSR_A"),
                      labels = c("MSR_D","MSR_A")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          strip.text =  element_text(colour = "black", size = "9"),
          plot.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  # plot_NPC <- ggplot(data = df_tidy %>% filter(protein_coding == 0)) + 
  #   geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"),
  #                alpha = 0.8) +
  #   geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
  #                alpha = 0.8) +
  #   ggtitle(paste0("Non-protein coding (NPC)\n")) +
  #   xlab("Mis-splicing ratio") +
  #   ggforce::facet_zoom(xlim = c(0,0.15)) +
  #   #ylim(y_axes) +
  #   ylab("") +
  #   #ylab("Intron count") +
  #   #scale_y_continuous(limits = c(0, 85000)) +
  #   
  #   theme_light() +
  #   scale_fill_manual(values = c("#35B779FF","#440154FF"),
  #                     breaks = c("#35B779FF","#440154FF"),
  #                     labels = c("MSR_D","MSR_A")) +
  #   theme(axis.line = element_line(colour = "black"), 
  #         axis.text = element_text(colour = "black", size = "9"),
  #         axis.title = element_text(colour = "black", size = "9"),
  #         plot.title = element_text(colour = "black", size = "9"),
  #         legend.text = element_text(size = "9"),
  #         legend.title = element_text(size = "9"),
  #         legend.position = "top") +
  #   guides(fill = guide_legend(title = NULL,
  #                              ncol = 2, 
  #                              nrow = 1))
  # 
  # plotMSR <- ggpubr::ggarrange(plot_PC, 
  #                           plot_NPC,
  #                           common.legend = T,
  #                           #labels = c("a.1","a.2"),
  #                           align = "hv",
  #                           ncol = 2, 
  #                           nrow = 1)

  
  
  
  
  
  
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
  
  
  # if (save_results) {
  #   plotLM
  #   filename <- paste0(folder_name, "/images/lm_MSR.png")
  #   ggplot2::ggsave(filename = filename, 
  #                   width = 183, height = 183, units = "mm", dpi = 300)
  # } else {
  #   plot %>% return()
  # }
  
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
    coord_flip() +
    facet_grid(vars(type)) +
    #ggtitle(graph_title) +
    ylab("Distribution of the estimate") +
    xlab(" ") +
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
    theme(axis.text.y = element_text(#angle = 70, 
                                     vjust = 0.5,
                                     hjust = 1)) +
    geom_hline(yintercept = 0,linetype='dotted')
  
  
  ggpubr::ggarrange(ggpubr::ggarrange(plotMSR,           
                                      plotLM,
                                      ncol = 2,
                                      labels = c("a", "b")),
                    plotTissuesLM, 
                    common.legend = F,
                    labels = c("", "c"),
                    #ncol = 2, 
                    nrow = 2,
                    #widths = c(1,1,2,2),
                    align = "v")
  
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel_4.png")
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



## SECTION 4 ---------------------------------------------
## AGE STRATIFICATION AND RBP EXPRESSION -----------------

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

encori <- function(project_id = "BRAIN") {
  
  #############################################################
  ## Load the genes whose MSR values increase with age
  #############################################################
  
  genes_increase_MSRD <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRD.rds")
  genes_increase_MSRA <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRA.rds")
  
  genes_increase_MSR <- intersect(genes_increase_MSRD, genes_increase_MSRA) %>% unique() %>% sort()
  
  write.csv(x = genes_increase_MSR,
            row.names = F,
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/genes_increase_MSR.csv")
  
  
  ####################################################################
  ## Load RBPs whose TPM levels show a negative linear 
  ## relationship with age
  ####################################################################
  
  df_lm_age_tidy <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                                           "/results/pipeline3/rbp/tpm_lm_all.csv")) %>%
    drop_na()
  
  ## PLOT DENSITIES (to decide the 'Estimate' values that means the RBP is affected by age)
  df_lm_age_tidy
  plot(density(RBPs_lm_age %>%
                 filter(type != "other") %>%
                 pull(Estimate)))
  lines(density(RBPs_lm_age %>%
                  filter(type == "other") %>%
                  pull(Estimate)), col = "red")
  
  ## RBPs affected by age
  RBPs_lm_age <- df_lm_age_tidy %>%
    filter(Estimate < -0.05, 
           type != "other",
           pval < 0.05) %>%
    as_tibble()
  
  write.csv(x = RBPs_lm_age$name %>% unique,
            row.names = F,
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_decrease.csv")
  
  
  ## RBPs not affected by age
  RBPs_lm_other <- df_lm_age_tidy %>%
    filter(RBP_ID %in% setdiff(df_lm_age_tidy$RBP_ID, RBPs_lm_age$RBP_ID))
  
  write.csv(x = RBPs_lm_other$name %>% unique,
            row.names = F,
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_other.csv")
  
  
  ####################################################################
  ## Load RBPs whose level of expression (i.e. TPM) decrease with age
  ####################################################################
  
  df_analysis_corrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                                                  "/results/pipeline3/rbp/tpm_age_spread.csv"), check.names = F) %>% as_tibble()
  RBPs <- df_analysis_corrected %>%
    filter(`20-39` > `40-59`,
           `40-59` > `60-79`) 
  
  ## Add RBPs the gene symbol
  ensembl106 <- biomaRt::useEnsembl(biomart = 'genes', 
                                    dataset = 'hsapiens_gene_ensembl',
                                    version = 106)
  RBP_symbols <- biomaRt::getBM(
    attributes = c('hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = RBPs$RBP,
    mart = ensembl106
  ) %>% pull(hgnc_symbol)
  
  write.csv(x = RBP_symbols,
            row.names = F,
            file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_decrease.csv")
  
  
  
  ####################################################################
  ## INTERSECT BOTH DATASETS OF RBPs
  ####################################################################
  RBPs_splicing$hgnc_symbol %>% unique %>% length()
  RBPs_other$hgnc_symbol %>% unique %>% length()
  RBP_symbols %>% unique %>% length()
  
  intersect(RBPs_splicing$hgnc_symbol, RBP_symbols)
  intersect(RBPs_other$hgnc_symbol, RBP_symbols)
  
  # curl 'https://starbase.sysu.edu.cn/api/RBPTarget/?assembly=hg19&geneType=mRNA&RBP=GEMIN5&clipExpNum=5&pancancerNum=0&target=APOE&cellType=all' > ENCORI_hg19_RBPTarget_GEMIN5_APOE.txt
  
}

analise_encori <- function(project_id = "BRAIN") {
  
  ## LOAD MSR
  
  df_RBP_result <- map_df(c("RBPs_affected_age", "RBPs_notaffected_age"), function(type) {
    
    # type <- "RBPs_notaffected_age"
    # type <- "RBPs_affected_age"
    
    ## LOAD RPB list
    if (type == "RBPs_affected_age") {
      RBPs <- read.csv(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_decrease.csv",
                       header = F)
    } else {
      RBPs <- read.csv(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_other.csv", 
                       header = F)
    }
    
    ## Check iCLIP results
    map_df(RBPs$V1, function (RBP) {
      
      # RBP <- RBPs$V1[1]
      file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/results/", 
                          type, "/ENCORI_hg19_", RBP, "_allgenes_copy.txt")
      
      #print(file_name)
      if (file.exists(file_name)) {
        
        RBP_result <- read.delim(file = file_name, header = T,sep = "\t") %>%
          as_tibble()
        
        if (RBP_result %>% nrow() > 1) {
          
          print(RBP)
        
          RBP_result_GR <- RBP_result %>% 
            dplyr::rename(start = narrowStart, end = narrowEnd) %>% 
            GRanges()
          
          
          # liftover
          
          # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
          chain <- rtracklayer::import.chain(con = "data/hg19ToHg38.over.chain")
  
          RBP_result_GRh38 <- rtracklayer::liftOver(x = RBP_result_GR, 
                                                    chain = chain) %>% 
            unlist() %>% 
            diffloop::rmchr()
          
          ## Overlaps - MSR_D
          
          genes_MSRD_increasing <- genes_MSRD_increasing %>%
            GRanges()
          
          overlaps_MSRD <- GenomicRanges::findOverlaps(query = genes_MSRD_increasing,
                                                       subject = RBP_result_GRh38,
                                                       type = "any",
                                                       ignore.strand = F)
          
          
          # genes_MSRD_increasing_RBP <- genes_MSRD_increasing[queryHits(overlaps),] %>%
          #   as_tibble() %>%
          #   mutate(RBP_name = RBP_result_GRh38[subjectHits(overlaps),]$RBP,
          #          RBP_gene.name = RBP_result_GRh38[subjectHits(overlaps),]$geneName,
          #          RBP_cellline.tissue = RBP_result_GRh38[subjectHits(overlaps),]$cellline.tissue) %>%
          #   unnest(gene_name)
          
          
          ## Overlaps - MSR_A
          
          genes_MSRA_increasing <- genes_MSRA_increasing %>%
            GRanges()
          
          overlaps_MSRA <- GenomicRanges::findOverlaps(query = genes_MSRA_increasing,
                                                       subject = RBP_result_GRh38 ,
                                                       type = "any",
                                                       ignore.strand = F)
          
          # genes_MSRA_increasing_RBP <- genes_MSRA_increasing[queryHits(overlaps),] %>%
          #   as_tibble() %>%
          #   mutate(RBP_name = RBP_result_GRh38[subjectHits(overlaps),]$RBP,
          #          RBP_gene.name = RBP_result_GRh38[subjectHits(overlaps),]$geneName,
          #          RBP_cellline.tissue = RBP_result_GRh38[subjectHits(overlaps),]$cellline.tissue) %>%
          #   unnest(gene_name)
          
          
          ## Return result
          
          # rbind(genes_MSRD_increasing_RBP,
          #       genes_MSRA_increasing_RBP) %>% return()
          #print("File exists!")
          
          data.frame(RBP_name = RBP,
                     RBP_type = type,
                     ovlps_MSRD = genes_MSRD_increasing[queryHits(overlaps_MSRD),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRA = genes_MSRA_increasing[queryHits(overlaps_MSRA),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRD = genes_MSRD_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRA = genes_MSRA_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRD_perc = (((genes_MSRD_increasing[queryHits(overlaps_MSRD),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                       genes_MSRD_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow()),
                     ovlps_MSRA_perc = (((genes_MSRA_increasing[queryHits(overlaps_MSRA),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                          genes_MSRA_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow())) %>% return()
        }
        
       
      
      }
      
    })
    
    # file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/results/", 
    #                     type, "/RBP_MSR.rds")
    # saveRDS(object = df_RBP_result, file = file_name)
    
  })
  
  write_csv(x = df_RBP_result,
            file = "paper_figures/iCLIP/results/RBPs_intronsMSR.csv", col_names = T)
}
 





  