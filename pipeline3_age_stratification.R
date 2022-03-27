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
  
  
  title <- paste0("Proportion of noise per distance locus\n", n_ref_jun, " common ref introns.")
  
  
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


age_stratification_RBP_analysis <- function(age_samples_clusters_tidy,
                                            project_id) {
  
  
  ###############################
  # # 10. Load output from Aine's function
  ###############################
  tpm_corrected <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                          project_id, "/results/pipeline3/rbp/tpm_residuals.csv"),
                            header = T)
  tpm_corrected[1,]
  age_samples_clusters_tidy[1,]
  
  tpm_20_39 <- tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "20-39") %>% pull(individual)))
  
  
  tpm_40_59 <-tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "40-59") %>% pull(individual))) 
  
  tpm_60_79 <- tpm_corrected %>%
    filter(X %in% (age_samples_clusters_tidy %>%
                     filter(age_group == "60-79") %>% pull(individual)))
  
  intersect(tpm_20_39$X, tpm_60_79$X)
  
  #####################
  ## U2AF1
  #####################
  
  wilcox.test(x = tpm_20_39$ENSG00000160201, 
              y = tpm_40_59$ENSG00000160201, 
              correct = T,
              alternative = "greater")
  
  wilcox.test(x = tpm_20_39$ENSG00000160201, 
              y = tpm_60_79$ENSG00000160201, 
              alternative = "greater")
  
  wilcox.test(x = tpm_40_59$ENSG00000160201, 
              y = tpm_60_79$ENSG00000160201,
              alternative = "greater")
  
}




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


 


