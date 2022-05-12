library(tidyverse)
library(data.table)
library(GenomicRanges)
library(ggforce)
# library(data.table)
# library(DBI)
# library(stringr)
# library(rtracklayer)
# library(ggforce)
# library(GenomicRanges)
# library(diffloop)
# library(ggpubr)
# library(matrixStats)
# library(xlsx)

##  source("/home/sruiz/PROJECTS/splicing-project/pipeline3_age_stratification.R")


################################
## FUNCTIONS       
################################

get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}


age_stratification_init_data <- function (project_id) {
  
  # project_id <- "BRAIN"
  
  if (!file.exists(paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/raw_data/sample_metadata.rds"))) {
    
    rse <- recount3::create_rse_manual(
      project = project_id,
      project_home = "data_sources/gtex",
      organism = "human",
      annotation = "gencode_v29",
      type = "gene")
    
    
    sample_metadata <- rse %>% 
      SummarizedExperiment::colData() %>%
      as_tibble() %>%
      filter(gtex.smafrze != "EXCLUDE")
    
    
    saveRDS(object = sample_metadata,
            file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/raw_data/sample_metadata.rds"))
    
    rm(rse)
    gc()
    
  } else {
    
    sample_metadata <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                                             "/raw_data/sample_metadata.rds"))
  }
  
  age_samples_clusters <- age_stratification_clustering(sample_metadata)
  
  age_samples_clusters_tidy <- age_samples_clusters %>%
    mutate(age_group = "")
  
  for (age_cluster in (age_samples_clusters$age %>% unique())) {
    
    # age_cluster = (age_samples_clusters$age %>% unique())[1]
    
    switch(age_cluster, 
           '20-29' = {
             indx <- which(age_samples_clusters_tidy$age == age_cluster)
             age_samples_clusters_tidy[indx, "age_group"] <- "20-39"
           },
           '30-39' = {
             indx <- which(age_samples_clusters_tidy$age == age_cluster)
             age_samples_clusters_tidy[indx, "age_group"] <- "20-39"
           },
           '40-49' = {
             indx <- which(age_samples_clusters_tidy$age == age_cluster)
             age_samples_clusters_tidy[indx, "age_group"] <- "40-59"
           },
           '50-59' = {
             indx <- which(age_samples_clusters_tidy$age == age_cluster)
             age_samples_clusters_tidy[indx, "age_group"] <- "40-59"
           },
           '60-69' = {
             indx <- which(age_samples_clusters_tidy$age == age_cluster)
             age_samples_clusters_tidy[indx, "age_group"] <- "60-79"
           },
           '70-79' = {
             indx <- which(age_samples_clusters_tidy$age == age_cluster)
             age_samples_clusters_tidy[indx, "age_group"] <- "60-79"
           },
           {
             print('not found')
           }
    )
  }
  
 return(age_samples_clusters_tidy)
  
}


age_stratification_clustering <- function (sample_metadata) {
  

  
  ## 1. Get key metadata
  
  clusters <- sample_metadata$gtex.smtsd %>% unique()
  age_groups <- sample_metadata$gtex.age %>% unique()
  samples <- sample_metadata$external_id %>% unique()
  
  ## 2. Loop through the individuals
  
  gtex_age_clustering <- map_df(samples, function(individual) {
    
    # individual <- (samples)[1]
    
    person <- sample_metadata %>%
      filter(external_id == individual)
    
    # print(person$smtsd %>% unique)
    
    if (person$gtex.age %>% unique() %in% age_groups) {
      
      if (any(person$gtex.smtsd %in% clusters)) {
        
        sex <- ifelse(person$gtex.sex == 1, "male", "female")
        
        data.frame(individual = person$external_id %>% unique(),
                   age = person$gtex.age %>% unique(),
                   rin = person$gtex.smrin,
                   mapped_read_count = person$recount_qc.star.all_mapped_reads,
                   sex = sex,
                   sample_recount_id = person$gtex.sampid,
                   region = person$gtex.smtsd) %>% return()
      }
      
    }
    
  }) 
  
  gtex_age_clustering %>% 
    return()

}


age_stratification_annotate <- function (age_groups,
                                         project_id = "BRAIN",
                                         age_samples_clusters) {
  
  ## Loop
  for (age_cluster in age_groups) {
    
    # age_cluster <- age_groups[1]
    
    ## SET RESULTS FOLDER
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/results/pipeline3/distances/age/", age_cluster, "/")
    dir.create(file.path(folder_name), showWarnings = T, recursive = T)
    
    #######################################################################
    ## GENERATE THE COUNTS AND ANNOTATED SPLIT READ COUNTS BY AGE CLUSTER
    #######################################################################
    
    # age_cluster <- age_groups[1]
    # age_cluster <- age_groups[2]
    
    all_split_reads_details_105_age <- NULL
    split_read_counts_age <- NULL
    samples_age <- NULL
    i <- 1
    
    for (cluster in age_samples_clusters$region %>% unique()) {
      
      # cluster <- (age_samples_clusters$region %>% unique())[1]
      # cluster <- (age_samples_clusters$region %>% unique())[2]
      
      print(paste0(Sys.time(), " - ", cluster, ", ", age_cluster))
      
      ## LOAD SAMPLES DATA
      samples <- age_samples_clusters %>% 
        filter(age_group == age_cluster, region == cluster) %>% 
        pull(individual)
      samples_age <- c(samples_age, samples)
      
      print(paste0(samples %>% length(), " samples... "))
      
      ## Load the annotated split reads
      all_split_reads_details_105 <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                                           project_id, "/results/base_data/",
                                                           cluster, "/", cluster, "_annotated_SR_details_length_105.rds"))
      
      # Load split read counts
      split_read_counts <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                                 project_id, "/results/base_data/",
                                                 cluster, "/", cluster, "_split_read_counts.rds")) %>%
        data.table::as.data.table() %>%
        select(c("junID", all_of(samples)))
      
      if (i == 1) {
        
        split_read_counts_age <- split_read_counts
        all_split_reads_details_105_age <- all_split_reads_details_105
        
      } else {
        
        all_split_reads_details_105_age <- data.table::rbindlist(list(all_split_reads_details_105,
                                                                      all_split_reads_details_105_age)) %>%
          distinct(junID, .keep_all = T)
        
        split_read_counts_age <- merge(x = split_read_counts %>% data.table::as.data.table(),
                                       y = split_read_counts_age %>% data.table::as.data.table(),
                                       by = "junID",
                                       all.x = T,
                                       all.y = T)
      }
      
      i <- i + 1
      print(paste0(split_read_counts_age %>% names() %>% length(), " samples processed in total."))
      
      rm(all_split_reads_details_105)
      rm(split_read_counts)
      rm(samples)
      gc()
      
      
      
    }
    
    ## SAVE RESULTS
    
    split_read_counts_age <- split_read_counts_age[rowSums(is.na(split_read_counts_age[,2:ncol(split_read_counts_age)])) != 
                                                     (ncol(split_read_counts_age)-1),]
    all_split_reads_details_105_age <- all_split_reads_details_105_age %>%
      filter(junID %in% (split_read_counts_age$junID %>% unique()))
    
    
    saveRDS(object = samples_age, file = paste0(folder_name, "/samples_used.rds"))
    saveRDS(object = split_read_counts_age %>% data.table::as.data.table(), file = paste0(folder_name, "/split_read_counts.rds"))
    saveRDS(object = all_split_reads_details_105_age, file = paste0(folder_name, "/all_split_reads_details_105.rds"))
    
    
    rm(all_split_reads_details_105_age)
    rm(split_read_counts_age)
    rm(samples_age)
    gc()
    
    
  }
}


age_stratification_IDB <- function (age_groups,
                                    project_id,
                                    age_samples_clusters) {
  
  source("/home/sruiz/PROJECTS/splicing-project/pipeline3_methods.R")
  
  
  
  ## Loop
  for (age_cluster in age_groups) {
    
    # age_cluster <- age_groups[1]
    # age_cluster <- age_groups[3]
    
    print(paste0(Sys.time(), " - Starting analysis for '", age_cluster, "' samples..."))
    
    ## SET RESULTS FOLDER
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                          project_id, "/results/pipeline3/distances/age/", age_cluster, "/")
    
    dir.create(file.path(folder_name), showWarnings = T, recursive = T)
    
    
    #######################################################################
    ## GENERATE THE IDB
    #######################################################################
    
    
    ## Load data ----------------------------------------------------------------
    
    samples_age <- readRDS(file = paste0(folder_name, "/samples_used.rds"))
    
    split_read_counts_age <- readRDS(file = paste0(folder_name, "/split_read_counts.rds")) %>%
      as_tibble()
    
    all_split_reads_details_105_age <- readRDS(file = paste0(folder_name, "/all_split_reads_details_105.rds")) %>%
      as_tibble()
    
    
    
    ## Call functions ----------------------------------------------------------------
    
    get_distances(cluster = age_cluster,
                  samples = samples_age,
                  split_read_counts = split_read_counts_age,
                  all_split_reads_details = all_split_reads_details_105_age,
                  folder_name = folder_name)
    

    extract_distances(cluster = age_cluster,
                      samples = samples_age,
                      split_read_counts = split_read_counts_age,
                      folder_name = folder_name)

    get_never_misspliced(cluster = age_cluster,
                         split_read_counts = split_read_counts_age,
                         all_split_reads_details = all_split_reads_details_105_age,
                         samples = samples_age,
                         folder_name = folder_name)
    
    extract_never_misspliced(cluster = age_cluster,
                             split_read_counts = split_read_counts_age,
                             samples = samples_age,
                             folder_name = folder_name)
    
    add_never_misspliced_to_df(cluster = age_cluster, 
                               all_split_reads_details = all_split_reads_details_105_age,
                               samples = samples_age,
                               folder_name = folder_name)
    
    
    ## Get mis-splicing ratio
    
    folder_idb_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/" ,
                              project_id, "/results/pipeline3/missplicing-ratio/age/", age_cluster)
    dir.create(file.path(folder_idb_name), showWarnings = F, recursive = T)
    
    get_missplicing_ratio(cluster = age_cluster,
                          split_read_counts = split_read_counts_age,
                          samples = samples_age,
                          folder_name = folder_name,
                          folder_save_name = folder_idb_name)
    
    add_missplicing_class_to_df(cluster = age_cluster,
                                folder_name = folder_idb_name)
    
    ############################################
    ## ADD FEATURES TO THE IDB
    ############################################
    
    remove_MT_genes(cluster = age_cluster,
                    folder_name = folder_idb_name)
    
    add_intron_type(cluster = age_cluster,
                    folder_name = folder_idb_name)
    
    clinvar_analysis(cluster = age_cluster,
                     folder_name = folder_idb_name)
    
    add_MANE_info(cluster = age_cluster,
                  folder_name = folder_idb_name)
    
    
    ############################################
    ## QC
    ############################################
    
    get_missplicing_QC(cluster = age_cluster,
                       split_read_counts = split_read_counts_age,
                       all_split_reads_details = all_split_reads_details_105_age,
                       samples = samples_age,
                       folder_name = folder_idb_name)
    
    
    
    # }
  }
}


age_stratification_plot_distances <- function(df = NULL,
                                              age_groups = NULL,
                                              common = T,
                                              age_levels = c("60-79", "40-59", "20-39"),
                                              distance_limit = 30,
                                              QC = F) {
  
  folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/age/"
  
  
  
  
  if (is.null(df)) {
    ## LOAD IDB across age supergroups ----------------------------------------------
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      
      #print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
      readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_novel.rds")) %>%
        mutate(sample_type = age_group) %>%
        return()
      
    }) 
    
    
    ## GET common junctions across age supergroups ----------------------------------
    
    if (common) {
      
      common_junctions <- df_age_groups %>%
        group_by(sample_type) %>%
        distinct(ref_junID, .keep_all = T) %>%
        ungroup() %>%
        dplyr::count(ref_junID) %>%
        filter(n == age_groups %>% length()) %>%
        pull(ref_junID)
      
      
      
      df_age_groups_tidy <- df_age_groups %>%
        filter(ref_junID %in% common_junctions)
      
      
      title <- paste0(title, "\n", common_junctions %>% unique() %>% length(), " common reference introns used.")
      
      # QC
      if (QC) {
        
        ## Each age group should have the same reference junction IDs
        
        for (age in age_groups) {
          
          print(age)
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(ref_junID) %>% 
            nrow() %>% 
            print()
          
          print("Novel junctions:")
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(novel_junID) %>% 
            nrow() %>% 
            print()
          
        }
        
        
        ## Second QC - novel junction IDs
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[2]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[2]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[3]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[3]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        
        
        
      }
    } else {
      df_age_groups_tidy <- df_age_groups 
    }
  } else {
    df_age_groups_tidy <- df
  }
  
  
  n_ref_jun <- df_age_groups_tidy %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type) %>%
    ungroup() %>%
    distinct(n) %>%
    pull
  
  
  title <- paste0("Distances to the reference junction\n", n_ref_jun, " common ref introns.")
  
  df_age_groups_tidy <- df_age_groups_tidy %>%
    mutate(novel_type = factor(novel_type, levels = c("novel_donor", "novel_acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = age_levels))
  
  
  
  ggplot(data = df_age_groups_tidy) + 
    geom_histogram(aes(x = distance, fill = sample_type),
                   alpha = 0.9,
                   bins = distance_limit * 2,
                   binwidth = 1,
                   position = "identity"
    ) +
    facet_grid(vars(novel_type)) +
    ggtitle(title) +
    xlab("Distance to the reference intron (in bp)") +
    ylab("Number of unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((distance_limit * -1), distance_limit),
                       breaks = c((distance_limit * -1), (round(distance_limit / 2) * -1), 0, round(distance_limit / 2), distance_limit)) +
    
    scale_fill_manual(values =  c("#440154FF","#FDE725FF","#21908CFF"),
                      breaks = age_levels) +
    
    # scale_fill_manual(breaks = c("20-39", "40-59", "60-79"),
    #                   labels = c("20-39", "40-59", "60-79")) +
    guides(fill = guide_legend(title = NULL, 
                               ncol = 4, nrow = 1 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          strip.text = element_text(colour = "black", size = "16"), 
          legend.text = element_text(colour = "black", size = "14"),
          plot.caption = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "16"),
          legend.title = element_text(colour = "black", size = "14"),
          legend.position = "top") %>% 
    return()
  
  
}

age_stratification_plot_distances_proportion <- function(df = NULL,
                                                         age_groups = NULL,
                                                         common = T,
                                                         age_levels = c("60-79", "40-59", "20-39"),
                                                         distance_limit = 30,
                                                         QC = F) {
  
  if (is.null(df)) {
    
    folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/age/"
    
    title <- paste0("Proportion of splicing noise at each distance point to the reference junction")
    
    ## LOAD IDBs
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      
      #print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
      readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_novel.rds")) %>%
        mutate(sample_type = age_group) %>%
        return()
      
    }) 
    
    
    ## GET common junctions
    
    if (common) {
      
      
      common_junctions <- df_age_groups %>%
        group_by(sample_type) %>%
        distinct(ref_junID, .keep_all = T) %>%
        ungroup() %>%
        dplyr::count(ref_junID) %>%
        filter(n == age_groups %>% length()) %>%
        pull(ref_junID)
      
      
      
      df_age_groups_tidy <- df_age_groups %>%
        filter(ref_junID %in% common_junctions)
      
      
      title <- paste0(title, "\n", common_junctions %>% unique() %>% length(), " common reference introns used.")
      
      # QC
      if (QC) {
        
        ## Each age group should have the same reference junction IDs
        
        for (age in age_groups) {
          
          print(age)
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(ref_junID) %>% 
            nrow() %>% 
            print()
          
          print("Novel junctions:")
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(novel_junID) %>% 
            nrow() %>% 
            print()
          
        }
        
        
        ## Second QC - novel junction IDs
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[2]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[2]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[3]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[3]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        
        
        
      }
    } else {
      
      df_age_groups_tidy <- df_age_groups 
    }
    
  } else {
    df_age_groups_tidy <- df 
  }
  
  
  n_ref_jun <- df_age_groups_tidy %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type) %>%
    ungroup() %>%
    distinct(n) %>%
    pull
  
  
  title <- paste0("Proportion of noise per distance locus.\n", n_ref_jun, " common ref introns.")
  
  
  df_age_groups_tidy <- df_age_groups_tidy %>%
    mutate(novel_type = factor(novel_type, levels = c("novel_donor", "novel_acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = age_levels))
  
  
  
  df_prop_noise <- df_age_groups_tidy %>%
    filter(abs(distance) <= distance_limit) %>%
    group_by(sample_type, novel_type) %>%
    mutate(N = n()) %>%
    ungroup() %>%
    group_by(sample_type, novel_type, distance) %>%
    mutate(n_distance = n()) %>%
    ungroup() %>%
    mutate(p_noise = n_distance/N)%>%
    as.data.frame()
  
  
  # df_prop_noise %>% 
  #   filter(distance == 15) %>% head
  # 
  # 
  # df_prop_noise %>%
  #   filter(sample_type == "20-39",
  #          novel_type == "novel_acceptor",
  #          distance == 20) %>% 
  #   nrow()
  # df_prop_noise %>%
  #   filter(sample_type == "20-39") %>% 
  #   nrow()
  
  
  
  
  ggplot(data = df_prop_noise) + 
    geom_col(aes(x = distance, y = p_noise, fill = sample_type),
             position = "dodge") +
    facet_grid(vars(novel_type)) +
    ggtitle(title) +
    xlab("Distance to the reference intron (in bp)") +
    ylab("Proportion of splicing noise") +
    #scale_fill_manual(values =  c("#FDE725FF", "#21908CFF", "#440154FF"),
    #                  breaks = c("20-39", "40-59", "60-79")) +
    guides(fill = guide_legend(title = NULL, 
                               ncol = 4, nrow = 1 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          strip.text = element_text(colour = "black", size = "16"), 
          legend.text = element_text(colour = "black", size = "14"),
          plot.caption = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "16"),
          legend.title = element_text(colour = "black", size = "14"),
          legend.position = "top") %>% 
    return()
  
  
}

age_stratification_mode_distances <- function(age_groups,
                                              common = T) {
  
  # age_groups = c("20-29", "40-49", "50-59", "60-69")
  folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/age/"
  
  
  ## LOAD IDBs
  
  df_age_groups <- map_df(age_groups, function(age_group) {
    
    print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
    readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_novel.rds")) %>%
      mutate(sample_type = age_group) %>%
      return()
    
  }) 
  
  
  df_age_modes <- map_df(age_groups, function(age_group) {
    
    # age_group <- age_groups[1]
    
    map_df(c("novel_acceptor", "novel_donor"), function(type) {
      
      # type <- "novel_acceptor"
      
      mode_in <- df_age_groups %>%
        filter(sample_type == age_group,
               novel_type == type,
               distance < 0) %>%
        pull(distance) %>%
        get_mode() 
      
      mode_ex <- df_age_groups %>%
        filter(sample_type == age_group,
               novel_type == type,
               distance > 0) %>%
        pull(distance) %>%
        get_mode() 
      
      
      return(data.frame(age = age_group,
                        novel_type = type,
                        mode_intron = mode_in,
                        mode_exon = mode_ex))
      
      
    })
    
  })
  
  return(df_age_modes)
  
}



age_stratification_Wilcoxon <- function(age_groups,
                                        common = T) {
  
  folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/age/"
  
  title <- paste0("MSR distribution - comparison across age supergroups")
  
  ## LOAD IDBs
  
  df_age_groups <- map_df(age_groups, function(age_group) {
    
    #print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
    readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_introns.rds")) %>%
      mutate(sample_type = age_group) %>%
      return()
    
  }) 
  
  
  
  
  
  ## GET common junctions
  
  if (common) {
    
    
    common_junctions <- df_age_groups %>%
      group_by(sample_type) %>%
      distinct(ref_junID, .keep_all = T) %>%
      ungroup() %>%
      dplyr::count(ref_junID) %>%
      filter(n == age_groups %>% length()) %>%
      pull(ref_junID)
    
    
    
    df_age_groups_tidy <- df_age_groups %>%
      filter(ref_junID %in% common_junctions)
    
    
    title <- paste0(title, "\n", common_junctions %>% unique() %>% length(), " common reference introns used.")
    
  } else {
    df_age_groups_tidy <- df_age_groups 
  }
  
  
  
  
  # df_age_groups_tidy <- df_age_groups_tidy %>%
  #   #mutate(novel_type = factor(novel_type, levels = c("novel_donor", "novel_acceptor"))) %>%
  #   mutate(sample_type = factor(sample_type, levels = c("60-79", "40-59", "20-39")))
  # 
  # limit_bp <- 60
  # 
  # ggplot(data = df_age_groups_tidy) + 
  #   geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = sample_type),
  #                position = "identity"
  #   ) +
  #   facet_zoom(xlim = c(0,0.05))
  # 
  # ggplot(data = df_age_groups_tidy) + 
  #   geom_boxplot(aes(x = sample_type, y = ref_missplicing_ratio_tissue_ND, fill = sample_type),
  #                position = "identity"
  #   ) 
  # 
  # ggplot(data = df_age_groups_tidy) + 
  #   geom_boxplot(aes(x = sample_type, y = ref_missplicing_ratio_tissue_NA, fill = sample_type),
  #                position = "identity"
  #   ) 
  
  ####################################
  
  
  ## MSR at donor and acceptor values increase in introns with age
  
  ## TESTING MSR_D
  
  df_correlation <- df_age_groups_tidy %>%
    select(ref_junID, sample_type, ref_missplicing_ratio_tissue_ND) %>%
    spread(sample_type, ref_missplicing_ratio_tissue_ND)
  
  
  wilcox.test(x = df_correlation$`20-39`,
              y = df_correlation$`40-59`,
              alternative = "less")
  
  wilcox.test(x = df_correlation$`20-39`,
              y = df_correlation$`60-79`,
              alternative = "less")
  
  wilcox.test(x = df_correlation$`40-59`,
              y = df_correlation$`60-79`,
              alternative = "less")
  
  
  ## TESTING MSR_A
  
  df_correlation <- df_age_groups_tidy %>%
    select(ref_junID, sample_type, ref_missplicing_ratio_tissue_NA) %>%
    spread(sample_type, ref_missplicing_ratio_tissue_NA)
  
  
  wilcox.test(x = df_correlation$`20-39`,
              y = df_correlation$`40-59`,
              alternative = "less")
  
  wilcox.test(x = df_correlation$`20-39`,
              y = df_correlation$`60-79`,
              alternative = "less")
  
  wilcox.test(x = df_correlation$`40-59`,
              y = df_correlation$`60-79`,
              alternative = "less")
  
  
  
  bg_genes <- df_age_groups_tidy %>% unnest(gene_name) %>% distinct(gene_name) %>% pull()
  
  #####################################
  ## MSR_D
  #####################################
  
  df_correlation_MSRD <- df_age_groups_tidy %>%
    select(ref_junID, sample_type, ref_missplicing_ratio_tissue_ND,gene_name) %>%
    spread(sample_type, ref_missplicing_ratio_tissue_ND)
  
  
  introns_always_less_MSRD <- df_correlation_MSRD %>%
    rowwise() %>%
    filter(`20-39` < `40-59`,
           `20-39` < `60-79`,
           `40-59` < `60-79`)
  genes_always_less_MSRD <- introns_always_less_MSRD %>% unnest(gene_name) %>% distinct(gene_name) %>% drop_na()
  
  
  
  introns_always_more_MSRD <- df_correlation_MSRD %>%
    rowwise() %>%
    filter(`20-39` > `40-59`,
           `20-39` > `60-79`,
           `40-59` > `60-79`)
  genes_always_more_MSRD <- introns_always_more_MSRD %>% unnest(gene_name) %>% distinct(gene_name) %>% drop_na()
  
  
  
  genes_less_MSRD <- setdiff(genes_always_less_MSRD$gene_name, genes_always_more_MSRD$gene_name)
  genes_more_MSRD <- setdiff(genes_always_more_MSRD$gene_name, genes_always_less_MSRD$gene_name)
  
  
  ## GO analysis
  
  gene_enrichment <- gprofiler2::gost(query = genes_less_MSRD,
                                      custom_bg = bg_genes,
                                      organism = "hsapiens",
                                      ordered_query = F,
                                      correction_method = "bonferroni",
                                      significant = T)
  gene_enrichment$result
  gene_enrichment$result %>%
    filter(source == "KEGG")
  
  
  gene_enrichment <- gprofiler2::gost(query = genes_more_MSRD,
                                      custom_bg = bg_genes,
                                      organism = "hsapiens",
                                      ordered_query = F,
                                      correction_method = "bonferroni",
                                      significant = T)
  gene_enrichment$result %>%
    filter(source == "KEGG")
  
  
  
  #####################################
  ## MSR_A
  #####################################
  
  df_correlation_MSRA <- df_age_groups_tidy %>%
    select(ref_junID, sample_type, ref_missplicing_ratio_tissue_NA,gene_name) %>%
    spread(sample_type, ref_missplicing_ratio_tissue_NA)
  
  introns_always_less_MSRA <- df_correlation_MSRA %>%
    rowwise() %>%
    filter(`20-39` < `40-59`,
           `20-39` < `60-79`,
           `40-59` < `60-79`)
  genes_always_less_MSRA <- introns_always_less_MSRA %>% unnest(gene_name) %>% distinct(gene_name) %>% drop_na()
  
  
  introns_always_more_MSRA <- df_correlation_MSRA %>%
    rowwise() %>%
    filter(`20-39` > `60-79`,
           `20-39` > `40-59`,
           `40-59` > `60-79`)
  genes_always_more_MSRA <- introns_always_more_MSRA %>% unnest(gene_name) %>% distinct(gene_name) %>% drop_na()
  
  
  genes_less_MSRA <- setdiff(genes_always_less_MSRA$gene_name, genes_always_more_MSRA$gene_name)
  genes_more_MSRA <- setdiff(genes_always_more_MSRA$gene_name, genes_always_less_MSRA$gene_name)
  
  
  
  ## GO analysis
  
  gene_enrichment <- gprofiler2::gost(query = genes_less_MSRA,
                                      custom_bg = bg_genes,
                                      organism = "hsapiens",
                                      ordered_query = F,
                                      correction_method = "bonferroni",
                                      significant = T)
  gene_enrichment$result 
  gene_enrichment$result %>% filter(source == "KEGG")
  
  
  gene_enrichment <- gprofiler2::gost(query = genes_more_MSRA,
                                      custom_bg = bg_genes,
                                      organism = "hsapiens",
                                      ordered_query = F,
                                      correction_method = "bonferroni",
                                      significant = T)
  gene_enrichment$result %>% filter(source == "KEGG")
  
  
  ## COMMON
  genes_common_less <- intersect(genes_less_MSRD, genes_less_MSRA)
  gene_enrichment <- gprofiler2::gost(query = genes_common_less,
                                      custom_bg = bg_genes,
                                      organism = "hsapiens",
                                      ordered_query = F,
                                      correction_method = "bonferroni",
                                      significant = T,
                                      sources = "KEGG")
  gene_enrichment$result
  any(genes_common_less == "PINK1")
  any(genes_common_less == "APOE")
  any(genes_common_less == "MAPT")
  any(genes_common_less == "SNCA")
  
  genes_common_more <- intersect(genes_more_MSRD, genes_more_MSRA)
  gene_enrichment <- gprofiler2::gost(query = genes_common_more,
                                      custom_bg = bg_genes,
                                      organism = "hsapiens",
                                      ordered_query = F,
                                      correction_method = "bonferroni",
                                      significant = T,
                                      sources = "KEGG")
  gene_enrichment$result
  
  
  
}


age_stratification_plot_MSRD <- function(df = NULL,
                                         age_groups = NULL,
                                         age_levels = c("60-79", "40-59", "20-39"),
                                         common = T,
                                         QC = F) {
  
  
  if (is.null(df)) {
    folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/age/"
    
    title <- paste0("Distances to the reference junction")
    
    
    ## LOAD IDB across age supergroups ----------------------------------------------
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      
      #print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
      readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_novel.rds")) %>%
        mutate(sample_type = age_group) %>%
        return()
      
    }) 
    
    
    ## GET common junctions across age supergroups ----------------------------------
    
    if (common) {
      
      common_junctions <- df_age_groups %>%
        group_by(sample_type) %>%
        distinct(ref_junID, .keep_all = T) %>%
        ungroup() %>%
        dplyr::count(ref_junID) %>%
        filter(n == age_groups %>% length()) %>%
        pull(ref_junID)
      
      
      
      df_age_groups_tidy <- df_age_groups %>%
        filter(ref_junID %in% common_junctions)
      
      
      title <- paste0(title, "\n", common_junctions %>% unique() %>% length(), " common reference introns used.")
      
      # QC
      if (QC) {
        
        ## Each age group should have the same reference junction IDs
        
        for (age in age_groups) {
          
          print(age)
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(ref_junID) %>% 
            nrow() %>% 
            print()
          
          print("Novel junctions:")
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(novel_junID) %>% 
            nrow() %>% 
            print()
          
        }
        
        
        ## Second QC - novel junction IDs
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[2]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[2]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[3]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[3]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        
        
        
      }
    } else {
      df_age_groups_tidy <- df_age_groups 
    }
  } else {
    df_age_groups_tidy <- df
  }
  
  n_ref_jun <- df_age_groups_tidy %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type) %>%
    ungroup() %>%
    distinct(n) %>%
    pull
  
  
  title <- paste0("MSR_Donor from common reference introns\n", n_ref_jun, " common ref introns.")
  
  
  
  df_age_groups_tidy <- df_age_groups_tidy %>%
    # mutate(ref_type = factor(ref_type, levels = c("novel_donor", "novel_acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = age_levels))
  
  ggplot(data = df_age_groups_tidy) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = sample_type), alpha = 0.8) +
    ggtitle(title) +
    xlab("MSR_Donor") +
    facet_zoom(xlim = c(0,0.2)) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text.x =  element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "14"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 3,  nrow = 1)) %>% return()
  
  
  
}

age_stratification_plot_MSRA <- function(df = NULL,
                                         age_groups = NULL,
                                         age_levels = c("60-79", "40-59", "20-39"),
                                         common = T,
                                         QC = F) {
  
  
  if (is.null(df)) {
    folder_root <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/age/"
    
    title <- paste0("Distances to the reference junction")
    
    
    ## LOAD IDB across age supergroups ----------------------------------------------
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      
      #print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
      readRDS(file = paste0(folder_root, "/", age_group, "/", age_group, "_db_novel.rds")) %>%
        mutate(sample_type = age_group) %>%
        return()
      
    }) 
    
    
    ## GET common junctions across age supergroups ----------------------------------
    
    if (common) {
      
      common_junctions <- df_age_groups %>%
        group_by(sample_type) %>%
        distinct(ref_junID, .keep_all = T) %>%
        ungroup() %>%
        dplyr::count(ref_junID) %>%
        filter(n == age_groups %>% length()) %>%
        pull(ref_junID)
      
      
      
      df_age_groups_tidy <- df_age_groups %>%
        filter(ref_junID %in% common_junctions)
      
      
      title <- paste0(title, "\n", common_junctions %>% unique() %>% length(), " common reference introns used.")
      
      # QC
      if (QC) {
        
        ## Each age group should have the same reference junction IDs
        
        for (age in age_groups) {
          
          print(age)
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(ref_junID) %>% 
            nrow() %>% 
            print()
          
          print("Novel junctions:")
          
          df_age_groups_tidy %>%
            filter(sample_type == age) %>%
            distinct(novel_junID) %>% 
            nrow() %>% 
            print()
          
        }
        
        
        ## Second QC - novel junction IDs
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[2]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[2]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                                subject = df_age_groups_tidy %>% filter(sample_type == age_groups[3]) %>% GenomicRanges::GRanges(),
                                                type = "equal")
        
        if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                       (df_age_groups_tidy %>% filter(sample_type == age_groups[3]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
          print("Massive error: some overlapping junctions don't have the same junction ID!")
        }
        
        
        
        
      }
    } else {
      df_age_groups_tidy <- df_age_groups 
    }
  } else {
    df_age_groups_tidy <- df
  }
  
  
  n_ref_jun <- df_age_groups_tidy %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type) %>%
    ungroup() %>%
    distinct(n) %>%
    pull
  
  
  title <- paste0("MSR_Acceptor from common reference introns\n", n_ref_jun, " common ref introns.")
  
  
  df_age_groups_tidy <- df_age_groups_tidy %>%
    # mutate(ref_type = factor(ref_type, levels = c("novel_donor", "novel_acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = age_levels))
  
  plot <- ggplot(data = df_age_groups_tidy) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = sample_type), alpha = 0.8) +
    ggtitle(title) +
    xlab("MSR_Acceptor") +
    ggforce::facet_zoom(xlim = c(0,0.2)) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text.x =  element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "14"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1))
  
  return(plot)
  
  
}

age_stratification_GO_enrichment <- function(genes,
                                             bg,
                                             ordered = F) {
  
  gene_enrichment <- gprofiler2::gost(query = genes,
                                      custom_bg = bg,
                                      organism = "hsapiens",
                                      ordered_query = ordered,
                                      correction_method = "bonferroni",
                                      significant = T)
  
  
  #gene_enrichment$result %>% as_tibble()
  gene_enrichment$result %>% 
    as_tibble() %>% 
    arrange(p_value) %>%
    mutate(p_value = p_value %>% formatC(format = "e", digits = 2)) %>%
    select(p_value, term_size, query_size, intersection_size, source, term_name) %>% 
    return()
}


age_stratification_RBP_corrected_TPM_analysis <- function(project_id = "BRAIN") {
  
  
  folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                        "/results/pipeline3/missplicing-ratio/age/")
  
  
  
  ## LOAD TPM CORRECTED VALUES
  tpm_corrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                          project_id, "/results/pipeline3/rbp/tpm_residuals.csv"), header = T)
  
  ## LOAD THE REFERENCE TRANSCRIPTOME
  ref <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  
  ## INIT AGE SUPERGROUPS
  age_samples_clusters_tidy <- age_stratification_init_data(project_id)


  
  ## FIlter TPM by age group and merge
  tpm_20_39 <- tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "20-39") %>% 
                     pull(individual))) %>%
    gather(key = "RBP", value = "TPM", -X ) %>%
    mutate(age_text = paste("20-39_", X)) %>%
    mutate(age = "20-39") %>% 
    as_tibble()
  
  
  tpm_40_59 <-tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "40-59") %>% 
                     pull(individual)))  %>%
    gather(key = "RBP", value = "TPM", -X ) %>%
    mutate(age_text = paste("40-59_", X)) %>%
    mutate(age = "40-59") %>% as_tibble()
  
  tpm_60_79 <- tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "60-79") %>% 
                     pull(individual))) %>%
    gather(key = "RBP", value = "TPM", -X ) %>%
    mutate(age_text = paste("60-79_", X)) %>%
    mutate(age = "60-79") %>% as_tibble()
  
  gattered_matrix <- do.call("rbind", list(tpm_20_39, tpm_40_59, tpm_60_79)) %>%
    dplyr::rename(sample = X) %>% as_tibble() %>%
    group_by(age, RBP) %>%
    mutate(mean_tpm = (TPM %>% mean())) %>%
    distinct(mean_tpm, .keep_all = T)  %>%
    ungroup() %>%
    arrange(age , mean_tpm) %>%
    mutate(RBP = fct_inorder(RBP))
  
  
  
  ## SAVE RESULTS
  spread_matrix <- gattered_matrix %>%
    dplyr::select(RBP, age, mean_tpm) %>%
    spread(key = age, value = mean_tpm)
  spread_matrix
  write_csv(x = spread_matrix,
            file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/results/pipeline3/rbp/tpm_age_spread.csv"))
  
  
  ## PLOT HEATMAP
  ggplot(data = gattered_matrix, 
         aes(x = age, y = RBP, fill = mean_tpm)) +
    geom_tile() + 
    scale_fill_gradient(low = "green", high = "red") +
    ggtitle("Mean RBP level of expression across samples\nfrom each age cluster.\nTPM values have been covariate corrected.")+
    ylab("RBP") +
    xlab("age cluster") +
    theme(axis.text.y = element_blank()) 
  
}

age_stratification_RBP_uncorrected_TPM_lm <- function(project_id = "BRAIN") {

  
  ref <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  
  ## Load and tidy the uncorrected TPMs
  ## TPMs here should be uncorrected as the covariates are going to be included in the linear models
  tpm_batch_corrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                                project_id, "/results/pipeline3/rbp/tpm.csv")) ##tpm_residuals.csv
  tpm_batch_corrected <- tpm_batch_corrected %>%
    dplyr::rename(gene = "X") %>%
    as_tibble()
  # tpm_batch_corrected_t <- t(tpm_batch_corrected) %>% as.data.frame() 
  # colnames(tpm_batch_corrected_t) <- tpm_batch_corrected_t[1,]
  # tpm_batch_corrected_t <- tpm_batch_corrected_t[-1,] %>% as_tibble(rownames = "gene")
  

  
  
  ## Load and tidy the sample metadata
  sample_metadata <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                              project_id, "/results/pipeline3/rbp/covariates.csv"), header = T,
                              fileEncoding = "UTF-8") %>% as_tibble()
  
  # sample_metadata_t <- t(sample_metadata) %>% as_tibble(rownames = "sample")
  # colnames(sample_metadata_t) <- sample_metadata_t[1,]
  # sample_metadata_t <- sample_metadata_t[-1,] %>%
  #   dplyr::rename(sample = "covariates") %>%
  #   mutate(gtex.age = gtex.age %>% as.double(),
  #          gtex.sex = gtex.sex %>% as.double(),
  #          gtex.smrin = gtex.smrin %>% as.double())
  
  
  
  ## Check the samples are the same
  identical(x = names(tpm_batch_corrected)[-1],
            y = names(sample_metadata)[-1])
  
  
  
  RBPs <- (tpm_batch_corrected$gene)
  
  ## Obtain values and run lm
  df_lm_output <- map_df(RBPs, function(RBP) {
    
    # RBP <- RBPs[1]
    
    print(RBP)
    
    tpm <- tpm_batch_corrected %>%
      filter(gene == RBP) %>%
      mutate(gene = "tpm")
    
    df <- rbind(tpm %>% 
                  dplyr::rename(covariates = "gene"),
                sample_metadata) %>%
      gather(sample, value, -covariates)  %>%
      spread(key = covariates, value) %>%
      dplyr::mutate(gtex.age = gtex.age %>% as.double(),
                    gtex.dthhrdy = gtex.dthhrdy %>% as.double(),
                    gtex.sex = gtex.sex %>% as.double(),
                    gtex.smcenter = gtex.smcenter %>% as.double(),
                    gtex.smgebtch = gtex.smgebtch %>% as.double(),
                    gtex.smgebtchd = gtex.smgebtchd %>% as.double(),
                    gtex.smnabtch = gtex.smnabtch %>% as.double(),
                    gtex.smnabtchd = gtex.smnabtchd %>% as.double(),
                    gtex.smnabtcht = gtex.smnabtcht %>% as.double(),
                    gtex.smrin = gtex.smrin %>% as.double(),
                    tpm = tpm %>% as.double()) %>%
      as_tibble()
    
    # print(df)
    lm_output <- lm(tpm ~ 
                      gtex.smrin + 
                      gtex.smcenter +
                      gtex.smgebtch +
                      gtex.smgebtchd +
                      gtex.smnabtch +
                      gtex.smnabtchd +
                      gtex.smnabtcht +
                      gtex.dthhrdy +
                      gtex.sex +
                      gtex.smtsd +
                      gtex.age,
                    data = df)
    
    lm_output <- lm_output %>% summary()
    df_lm_output <- lm_output$coefficients %>% as_tibble(rownames = "covariate")
    
    # cor(df)
    df_lm_output <- df_lm_output %>%
      mutate(RBP_ID = RBP) %>%
      dplyr::select(covariate, Estimate, pval = `Pr(>|t|)`, RBP_ID)
    
    return(df_lm_output)
    
  }) 
  
  
  # any(df_lm_output$pval > 0.05)
  
 
  ## Filter by age covariate and order
  df_lm_output_age <- df_lm_output %>%
    filter( pval <= 0.05,
            str_detect(string = covariate, pattern = "gtex.age")) %>%
    arrange(Estimate)
  
  ## Add gene SYMBOL info
  df_lm_age_tidy <- left_join(x = df_lm_output_age, 
                              y = RBPs_tidy, 
                              by = c("RBP_ID" = "ensembl_gene_id")) %>%
    group_by(RBP_ID) %>%
    distinct(pval, .keep_all = T)
  
  
  write_csv(x = df_lm_age_tidy,
            file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/results/pipeline3/rbp/tpm_lm.csv"))
  
  ###############################################
  ## Get RBPs that decrease expression with age
  ###############################################
  
  
  RBP_decrease <- df_lm_age_tidy %>%
    filter(Estimate < 0) %>%
    arrange(desc(abs(Estimate)))
  
  RBP_increase <- df_lm_age_tidy %>%
    filter(Estimate > 0) %>%
    arrange(desc(Estimate))
  
  intersect(RBP_decrease$hgnc_symbol,
            RBP_increase$hgnc_symbol)
  
  RBP_decrease$hgnc_symbol %>% unique() %>% sort()
  RBP_increase$hgnc_symbol %>% unique() %>% sort()
  
  # ###############################################
  # ## GO ENRICHMENT
  # ###############################################
  # 
  # ## Get RBPs that increase expression with age
  # 
  # gene_enrichment <- gprofiler2::gost(query = RBP_decrease$hgnc_symbol %>% unique(),
  #                                     custom_bg = RBP_bg,
  #                                     organism = "hsapiens",
  #                                     ordered_query = T,
  #                                     correction_method = "bonferroni",
  #                                     significant = T)
  # 
  # if (!is.null(gene_enrichment)) {
  # GO_enrichment <- gene_enrichment$result %>% 
  #   filter(str_detect(source, pattern = "GO")) %>%
  #   mutate(go_type = str_sub(source, start = 4, end = 5))
  # 
  # GO_enrich_reduc <- rutils::go_reduce(pathway_df = data.frame(go_type = GO_enrichment$go_type,
  #                                                              go_id = GO_enrichment$term_id,
  #                                                              go_term = GO_enrichment$term_name),
  #                                      orgdb = "org.Hs.eg.db",
  #                                      #threshold = 0.7,
  #                                      scores = NULL,
  #                                      measure = "Wang")
  # 
  # df_result <- merge(x = gene_enrichment$result %>% 
  #                      filter(source %in% c("GO:BP", "GO:MF", "GO:CC", "KEGG")),
  #                    y = GO_enrich_reduc,
  #                    by.x = "term_id",
  #                    by.y = "go_id",
  #                    all.x = T)
  # 
  # rbind(df_result %>%
  #         as_tibble() %>%
  #         dplyr::select(parent_term, p_value, query_size, intersection_size, source) %>%
  #         group_by(parent_term, source) %>%
  #         mutate(p_value = p_value %>% max) %>%
  #         ungroup() %>%
  #         distinct(parent_term, .keep_all = T) %>%
  #         drop_na(), 
  #       df_result %>%
  #         as_tibble() %>%
  #         dplyr::select(parent_term = term_name, p_value, query_size, intersection_size, source) %>%
  #         filter(source == "KEGG") %>%
  #         group_by(parent_term, source) %>%
  #         mutate(p_value = p_value %>% max) %>%
  #         ungroup() %>%
  #         distinct(parent_term, .keep_all = T))%>%
  #   arrange(p_value) %>%
  #   mutate(p_value = p_value %>% formatC(format = "e", digits = 2)) %>%
  #   as.data.frame()
  # }
  # 
  # 
  # ##############################################
  # ## Get RBPs that increase expression with age
  # ##############################################
  # 
  # 
  # 
  # gene_enrichment <- gprofiler2::gost(query = RBP_increase$hgnc_symbol %>% unique,
  #                                     custom_bg = RBP_bg,
  #                                     organism = "hsapiens",
  #                                     ordered_query = T,
  #                                     correction_method = "bonferroni",
  #                                     significant = T)
  # if (!is.null(gene_enrichment)) {
  # gene_enrichment$result%>%
  #   dplyr::select(parent_term = term_name, p_value, query_size, intersection_size, source) %>% 
  #   filter(str_detect(source, pattern = c("GO","WP")))
  # }
  
  
  
  
  ##############################################################################
  ##############################################################################
  
  
  
 
}

# 
# tpm_age_spread <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                                          "/results/pipeline3/rbp/tpm_age_spread.csv")) %>%
#   as_tibble() %>%
#   dplyr::rename(RBP_id = RBP)
# 
# tpm_lm <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
#                                  "/results/pipeline3/rbp/tpm_lm.csv")) %>%
#   filter(pval < 0.05) %>% 
#   dplyr::select(-pval) %>%
#   spread(key = covariate, value = "Estimate") %>%
#   replace(is.na(.),0) %>%
#   
#   as_tibble()
# 
# merge(x = tpm_age_spread,
#       y = tpm_lm,
#       by = "RBP_id")
#   

################################
## CALLS
################################


# age_samples_clusters_tidy <- age_samples_clusters %>% dplyr::rename(age_group = age)

# age_stratification_annotate(age_groups = clusters,
#                             age_samples_clusters = age_samples_clusters_tidy,
#                             project_id = project_id)
# 
# age_stratification_IDB(age_groups = clusters,
#                        project_id = project_id,
#                        age_samples_clusters = age_samples_clusters_tidy)


 


