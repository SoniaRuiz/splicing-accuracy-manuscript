
## Connect to the database -----------------------------------------------------
#source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline3_idb_SQL_generation.R")

gtf_version <- 105
main_project <- "age"
database_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v",
                        gtf_version, "/", main_project, "/", main_project, ".sqlite")

##################################
## CALLS - PLOTS
##################################

age_stratification_plot_distances <- function(df = NULL,
                                              age_groups = c("60-79", "40-59", "20-39"),
                                              distance_limit = 30,
                                              QC = F) {
  
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  
  df_age_distances <- map_df((df_metadata$SRA_project %>% unique())[1], function(project_id) {
    
    # project_id <- (df_metadata$SRA_project %>% unique())[1]
    
    print(paste0(Sys.time(), " --> ", project_id))
    
    age_groups <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      # age_group <- age_groups[1]
      print(paste0(Sys.time(), " --> ", age_group))
      
      query <- paste0("SELECT novel_junID FROM '", age_group, "_", project_id, "_misspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                      paste(introns$novel_junID, collapse = ","),")")
      df_novel <- introns %>%
        inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                   by = "novel_junID") %>%
        mutate(sample_type = age_group,
               project_id = project_id)
      
      return(df_novel)
      
    })
    
  })
  
  ## GET common junctions across age supergroups ----------------------------------
  
  common_junctions <- df_age_distances %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  
  df_age_distances_common <- df_age_distances %>%
    group_by(ref_junID, sample_type) %>%
    distinct(novel_coordinates, .keep_all = T)  %>%
    ungroup() %>%
    filter(ref_junID %in% common_junctions) %>%
    as_tibble()
  
  
  
  if (QC) {
    
    ## Each age group should have the same reference junction IDs
    for (age in age_groups) {
      
      print(age)
      
      df_age_distances %>%
        filter(sample_type == age) %>%
        distinct(ref_junID) %>% 
        nrow() %>% 
        print()
      
      print("Novel junctions:")
      
      df_age_distances %>%
        filter(sample_type == age) %>%
        distinct(novel_junID) %>% 
        nrow() %>% 
        print()
      
    }
    
    
    ## Second QC - novel junction IDs
    
    overlaps <- GenomicRanges::findOverlaps(query = df_age_distances %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                            subject = df_age_distances %>% filter(sample_type == age_groups[2]) %>% GenomicRanges::GRanges(),
                                            type = "equal")
    
    if (!identical((df_age_distances %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                   (df_age_distances %>% filter(sample_type == age_groups[2]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
      print("Massive error: some overlapping junctions don't have the same junction ID!")
    }
    
    overlaps <- GenomicRanges::findOverlaps(query = df_age_distances %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
                                            subject = df_age_distances %>% filter(sample_type == age_groups[3]) %>% GenomicRanges::GRanges(),
                                            type = "equal")
    
    if (!identical((df_age_distances %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
                   (df_age_distances %>% filter(sample_type == age_groups[3]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
      print("Massive error: some overlapping junctions don't have the same junction ID!")
    }
    
    
    
    
  }
  
  
  n_ref_jun <- df_age_distances_common %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type) %>%
    ungroup() %>%
    distinct(n) %>%
    pull
  
  df_age_distances_common <- df_age_distances_common %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  
  df_age_distances_common <- df_age_distances_common %>%
    mutate(novel_type = factor(novel_type, levels = c("novel donor", "novel acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = c("60-79", "40-59", "20-39")))
  
  
  
  plot_distances <- ggplot(data = df_age_distances_common) + 
    geom_histogram(aes(x = distance, fill = sample_type),
                   alpha = 0.9,
                   bins = distance_limit * 2,
                   binwidth = 1,
                   position = "identity"
    ) +
    facet_grid(vars(novel_type)) +
    #ggtitle(title) +
    xlab("Distance to the reference intron (in bp)") +
    ylab("Number of unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((distance_limit * -1), distance_limit),
                       breaks = c((distance_limit * -1), (round(distance_limit / 2) * -1), 0, round(distance_limit / 2), distance_limit)) +
    
    scale_fill_manual(values =  c("#21908CFF","#FDE725FF","#440154FF"),
                      labels = c("20-39", "40-59", "60-79"),
                      breaks = c("20-39", "40-59", "60-79")) +
    
    # scale_fill_manual(breaks = c("20-39", "40-59", "60-79"),
    #                   labels = c("20-39", "40-59", "60-79")) +
    guides(fill = guide_legend(title = NULL, 
                               ncol = 4, nrow = 1 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text = element_text(colour = "black", size = "12"), 
          legend.text = element_text(colour = "black", size = "12"),
          plot.caption = element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.title = element_text(colour = "black", size = "12"),
          legend.position = "top") %>% 
    return()
  
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = distance_limit, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
    geom_rect(aes(xmin = (distance_limit)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
    theme_void()
  
  
  plot_distances / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", 
                      main_project, "/figures/panel4_age_distances_blood.svg")
  ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
}

age_stratification_plot_distances_proportion <- function(age_levels = c("60-79", "40-59", "20-39"),
                                                         distance_limit = 40,
                                                         QC = F) {
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  
  con <- dbConnect(RSQLite::SQLite(), "./database/age-stratification.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  
  df_age_distances <- map_df((df_metadata$SRA_project %>% unique())[2], function(project_id) {
    
    # project_id <- (df_metadata$SRA_project %>% unique())[1]
    
    print(paste0(Sys.time(), " --> ", project_id))
    
    age_groups <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      pull(cluster)
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      # age_group <- age_groups[1]
      
      query <- paste0("SELECT novel_junID FROM '", age_group, "_", project_id, "_misspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                      paste(introns$novel_junID, collapse = ","),")")
      df_novel <- merge(x = introns,
                        y = dbGetQuery(con, query) %>% as_tibble(),
                        by = "novel_junID",
                        all.x = T) %>%
        mutate(sample_type = age_group,
               project_id = project_id)
      
      return(df_novel)
      
    })
    
  })
  
  ## GET common junctions across age supergroups ----------------------------------
  
  
  
  common_junctions <- df_age_distances %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  
  common_junctions %>% length()
  
  
  df_age_distances_tidy <- df_age_distances %>%
    group_by(ref_junID, sample_type) %>%
    distinct(novel_coordinates, .keep_all = T)  %>%
    ungroup() %>%
    filter(ref_junID %in% common_junctions) %>%
    as_tibble()
  
  
  
  n_ref_jun <- df_age_distances_tidy %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type) %>%
    ungroup() %>%
    distinct(n) %>%
    pull
  
  df_age_distances_tidy <- df_age_distances_tidy %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  
  
  df_age_distances_tidy <- df_age_distances_tidy %>%
    mutate(novel_type = factor(novel_type, levels = c("novel donor", "novel acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = age_levels))
  
  
  
  df_prop_noise <- df_age_distances_tidy %>%
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
  
  
  
  
  plot_distances <- ggplot(data = df_prop_noise) + 
    geom_col(aes(x = distance, y = p_noise, fill = sample_type),
             position = "identity") +
    facet_grid(vars(novel_type)) +
    xlab("Distance to the reference intron (in bp)") +
    ylab("Proportion of splicing noise") +
    theme_light() +
    scale_fill_manual(values =  c("#21908CFF","#FDE725FF","#440154FF"),
                      labels = c("20-39", "40-59", "60-79"),
                      breaks = c("20-39", "40-59", "60-79")) +
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
          legend.position = "top")
  
  distance_rectangle <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = distance_limit, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
    geom_rect(aes(xmin = (distance_limit)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
    theme_void()
  
  
  plot_distances / distance_rectangle + patchwork::plot_layout(heights = c(8, 1))  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/panel4_age_prop_noise.svg")
  ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
}

age_stratification_mode_distances <- function(age_groups,
                                              common = T) {
  
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  
  con <- dbConnect(RSQLite::SQLite(), "./database/age-stratification.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  
  df_age_distances <- map_df((df_metadata$SRA_project %>% unique())[2], function(project_id) {
    
    # project_id <- (df_metadata$SRA_project %>% unique())[1]
    
    print(paste0(Sys.time(), " --> ", project_id))
    
    age_groups <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      pull(cluster)
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      # age_group <- age_groups[1]
      
      query <- paste0("SELECT novel_junID FROM '", age_group, "_", project_id, "_misspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                      paste(introns$novel_junID, collapse = ","),")")
      df_novel <- merge(x = introns,
                        y = dbGetQuery(con, query) %>% as_tibble(),
                        by = "novel_junID",
                        all.x = T) %>%
        mutate(sample_type = age_group,
               project_id = project_id)
      
      return(df_novel)
      
    })
    
  })
  
  ## GET common junctions across age supergroups ----------------------------------
  
  
  
  common_junctions <- df_age_distances %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  
  
  
  df_age_distances_tidy <- df_age_distances %>%
    group_by(ref_junID, sample_type) %>%
    distinct(novel_coordinates, .keep_all = T)  %>%
    ungroup() %>%
    filter(ref_junID %in% common_junctions) %>%
    as_tibble()
  
  
  
  df_age_modes <- map_df(age_groups, function(age_group) {
    
    # age_group <- age_groups[1]
    
    map_df(c("novel_acceptor", "novel_donor"), function(type) {
      
      # type <- "novel_acceptor"
      
      mode_in <- df_age_distances %>%
        filter(sample_type == age_group,
               novel_type == type,
               distance < 0) %>%
        pull(distance) %>%
        get_mode() 
      
      mode_ex <- df_age_distances %>%
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
  
  con <- dbConnect(RSQLite::SQLite(), "./database/age-stratification.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  
  df_age_groups <- map_df((df_metadata$SRA_project %>% unique())[2], function(project_id) {
    
    # project_id <- (df_metadata$SRA_project %>% unique())[1]
    
    print(paste0(Sys.time(), " --> ", project_id))
    
    age_groups <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      pull(cluster)
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      # age_group <- age_groups[1]
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_type, MSR_D, MSR_A, gene_id FROM '", age_group, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT DISTINCT ref_junID, ref_type, MSR_D, MSR_A, gene_id FROM '", age_group, "_", project_id, "_misspliced'")
      introns <- plyr::rbind.fill(introns, dbGetQuery(con, query)) %>% as_tibble()
      
      query <- paste0("SELECT id, gene_name FROM 'gene' WHERE id IN (", paste(introns$gene_id, collapse = ","),")")
      genes <- dbGetQuery(con, query) %>% as_tibble()
      
      
      introns <- merge(x = introns,
                       y = genes,
                       by.x = "gene_id",
                       by.y = "id", all.x = T) %>%
        mutate(sample_type = age_group) %>%
        as_tibble()
      
      return(introns)
      
    })
    
  })
  
  common_junctions <- df_age_groups %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  
  
  
  df_age_groups_tidy <- df_age_groups %>%
    filter(ref_junID %in% common_junctions)
  
  
  
  ####################################
  ## TESTING MSR_D
  ####################################
  
  df_correlation <- df_age_groups_tidy %>%
    select(ref_junID, sample_type, MSR_D) %>%
    spread(sample_type, MSR_D)
  
  
  wilcox.test(x = df_correlation$`20-39`,
              y = df_correlation$`40-59`,
              alternative = "less")
  
  wilcox.test(x = df_correlation$`20-39`,
              y = df_correlation$`60-79`,
              alternative = "less")
  
  wilcox.test(x = df_correlation$`40-59`,
              y = df_correlation$`60-79`,
              alternative = "less")
  
  
  ####################################
  ## TESTING MSR_A
  ####################################
  
  df_correlation <- df_age_groups_tidy %>%
    select(ref_junID, sample_type, MSR_A) %>%
    spread(sample_type, MSR_A)
  
  
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
    select(ref_junID, sample_type, MSR_D ,gene_name) %>%
    spread(sample_type, MSR_D)
  
  
  introns_always_less_MSRD <- df_correlation_MSRD %>%
    rowwise() %>%
    filter(`20-39` < `40-59`,
           `20-39` < `60-79`,
           `40-59` < `60-79`)
  genes_always_less_MSRD <- introns_always_less_MSRD %>% 
    unnest(gene_name) %>% 
    distinct(gene_name) %>% 
    drop_na()
  
  
  
  introns_always_more_MSRD <- df_correlation_MSRD %>%
    rowwise() %>%
    filter(`20-39` > `40-59`,
           `20-39` > `60-79`,
           `40-59` > `60-79`)
  genes_always_more_MSRD <- introns_always_more_MSRD %>% 
    unnest(gene_name) %>% 
    distinct(gene_name) %>% 
    drop_na()
  
  
  
  genes_less_MSRD <- setdiff(genes_always_less_MSRD$gene_name, 
                             genes_always_more_MSRD$gene_name)
  genes_more_MSRD <- setdiff(genes_always_more_MSRD$gene_name, 
                             genes_always_less_MSRD$gene_name)
  
  
  ## GO analysis
  
  gene_enrichment_less_MSRD <- gprofiler2::gost(query = genes_less_MSRD,
                                                custom_bg = bg_genes,
                                                organism = "hsapiens",
                                                ordered_query = F,
                                                correction_method = "bonferroni",
                                                significant = T)
  gene_enrichment_less_MSRD$result
  gene_enrichment_less_MSRD$result %>%
    filter(source == "KEGG")
  
  
  gene_enrichment_more_MSRD <- gprofiler2::gost(query = genes_more_MSRD,
                                                custom_bg = bg_genes,
                                                organism = "hsapiens",
                                                ordered_query = F,
                                                correction_method = "bonferroni",
                                                significant = T)
  gene_enrichment_more_MSRD$result %>%
    filter(source == "KEGG")
  
  
  
  #####################################
  ## MSR_A
  #####################################
  
  df_correlation_MSRA <- df_age_groups_tidy %>%
    select(ref_junID, sample_type, MSR_A, gene_name) %>%
    spread(sample_type, MSR_A)
  
  introns_always_less_MSRA <- df_correlation_MSRA %>%
    rowwise() %>%
    filter(`20-39` < `40-59`,
           `20-39` < `60-79`,
           `40-59` < `60-79`)
  genes_always_less_MSRA <- introns_always_less_MSRA %>% 
    unnest(gene_name) %>% 
    distinct(gene_name) %>% 
    drop_na()
  
  
  introns_always_more_MSRA <- df_correlation_MSRA %>%
    rowwise() %>%
    filter(`20-39` > `60-79`,
           `20-39` > `40-59`,
           `40-59` > `60-79`)
  genes_always_more_MSRA <- introns_always_more_MSRA %>% 
    unnest(gene_name) %>% 
    distinct(gene_name) %>% 
    drop_na()
  
  
  genes_less_MSRA <- setdiff(genes_always_less_MSRA$gene_name, 
                             genes_always_more_MSRA$gene_name)
  genes_more_MSRA <- setdiff(genes_always_more_MSRA$gene_name, 
                             genes_always_less_MSRA$gene_name)
  
  
  
  ## GO analysis
  
  gene_enrichment_less_MSRA <- gprofiler2::gost(query = genes_less_MSRA,
                                                custom_bg = bg_genes,
                                                organism = "hsapiens",
                                                ordered_query = F,
                                                correction_method = "bonferroni",
                                                significant = T)
  gene_enrichment_less_MSRA$result 
  gene_enrichment_less_MSRA$result %>% filter(source == "KEGG")
  
  
  gene_enrichment_more_MSRA <- gprofiler2::gost(query = genes_more_MSRA,
                                                custom_bg = bg_genes,
                                                organism = "hsapiens",
                                                ordered_query = F,
                                                correction_method = "bonferroni",
                                                significant = T)
  gene_enrichment_more_MSRA$result %>% filter(source == "KEGG")
  
  
  
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
  any(genes_common_more == "PINK1")
  any(genes_common_more == "APOE")
  any(genes_common_more == "MAPT")
  any(genes_common_more == "SNCA")
  
  
}


age_stratification_plot_MSR <- function(df = NULL,
                                        age_groups = c("60-79", "40-59", "20-39"),
                                        age_levels = c("60-79", "40-59", "20-39"),
                                        common = T,
                                        QC = F) {
  
  
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  

  
  
  
  df_age_groups <- map_df((df_metadata$SRA_project %>% unique())[3], function(project_id) {
    
    # project_id <- (df_metadata$SRA_project %>% unique())[1]
    
    print(paste0(Sys.time(), " --> ", project_id))
    
    age_groups <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      # age_group <- age_groups[1]
      print(paste0(Sys.time(), " --> ", age_group))
      query <- paste0("SELECT ref_junID, ref_type, MSR_D, MSR_A FROM '", age_group, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT ref_junID, ref_type, MSR_D, MSR_A FROM '", age_group, "_", project_id, "_misspliced'")
      introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
      
      return(introns %>%
               mutate(sample_type = age_group))
      
    })
    
  })
  
  
  ###############################################
  ## QC
  ###############################################
  
  any(df_age_groups %>% filter(sample_type == "never") %>% pull(MSR_D) > 0)
  any(df_age_groups %>% filter(sample_type == "never") %>% pull(MSR_A) > 0)
  
  intersect(df_age_groups %>% filter(sample_type == "never") %>% pull(ref_junID),
            df_age_groups %>% filter(sample_type != "never") %>% pull(ref_junID))
  
  ###############################################
  ## GET ONLY COMMON JUNCTIONS ACROSS AGE GROUPS
  ###############################################
  
  common_introns <- df_age_groups %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  
  common_introns %>% length()
  
  ## Filter the INTRONS table by the common mis-spliced introns
  df_age_groups_tidy <- merge(x = df_age_groups %>% data.table::as.data.table(),
                              y = data.table::data.table(ref_junID = common_introns),
                              by = "ref_junID",
                              all.y = T) %>%
    distinct(ref_junID, sample_type, .keep_all = T)
  
  
  
  #######################################
  ## EXTRA - CHECK INTRONS 
  ## MSR INCREASING WITH AGE
  #######################################
  
  df_MSRD <- df_age_groups_tidy %>%
    distinct(ref_junID, sample_type, .keep_all = T) %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_A) %>%
    mutate(MSR_A = MSR_A %>% round(digits = 4)) %>%
    spread(sample_type, MSR_A)
  
  df_MSRD %>%
    filter(`20-39` < `40-59`,
           `40-59` < `60-79`)
  df_MSRD %>%
    filter(`20-39` > `40-59`,
           `40-59` > `60-79`)
  
  
  #######################################
  ## CONTINUE PLOT
  #######################################
  
  n_ref_jun <- df_age_groups_tidy %>%
    group_by(sample_type) %>%
    distinct(ref_junID) %>%
    dplyr::count(sample_type) %>%
    ungroup() %>%
    distinct(n) %>%
    pull
  
  
  df_age_groups_tidy <- df_age_groups_tidy %>%
    distinct(ref_junID, sample_type, .keep_all = T) %>%
    # mutate(ref_type = factor(ref_type, levels = c("novel_donor", "novel_acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = age_levels)) %>%
    as_tibble()
  
  
  ggplot(data = df_age_groups_tidy %>% 
           dplyr::select(MSR_A, sample_type) %>%
           gather(key = "MSR_type", value = "MSR", -sample_type) %>%
           mutate(MSR_type = factor(MSR_type, levels = c("MSR_A")))) + 
    geom_density(aes(x = MSR, fill = sample_type), alpha = 0.8) +
    #ggtitle(title) +
    xlab("MSR") +
    facet_grid(vars(MSR_type)) +
    ggforce::facet_zoom(xlim = c(0,0.2)) +
    theme_light() +
    scale_fill_manual(values =  c("#21908CFF","#FDE725FF","#440154FF"),
                      labels = c("20-39", "40-59", "60-79"),
                      breaks = c("20-39", "40-59", "60-79")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          strip.text =  element_text(colour = "black", size = "9"),
          plot.title = element_text(colour = "black", size = "9"),
          legend.text = element_text(size = "9"),
          legend.title = element_text(size = "9"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 3,  nrow = 1)) %>% return()
  
  
  ########################################
  ## STATISTICAL TEST - MSR_DONOR
  ########################################
  
  wilcox.test(x = df_age_groups_tidy %>% filter(sample_type == "20-39") %>% pull(MSR_D),
              y = df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_D),
              alternative = "less",
              paired = T)
  wilcox.test(x = df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_D),
              y = df_age_groups_tidy %>% filter(sample_type == "60-79") %>% pull(MSR_D),
              alternative = "less",
              paired = T)
  
  t.test(x = df_age_groups_tidy %>% filter(sample_type == "20-39") %>% pull(MSR_D),
         y = df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_D),
         alternative = "less",
         paired = T)
  
  ########################################
  ## STATISTICAL TEST - MSR_ACCEPTOR
  ########################################
  
  wilcox.test(x = df_MSRD$`20-39`,
              y = df_MSRD$`40-59`,
              alternative = "less",
              paired = T)
  wilcox.test(x = df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_A),
              y = df_age_groups_tidy %>% filter(sample_type == "60-79") %>% pull(MSR_A),
              alternative = "less",
              correct = T)
  
  
  
  df_age_groups_tidy %>% filter(sample_type == "20-39") %>% pull(MSR_A) %>% summary()
  df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_A) %>% summary()
  df_age_groups_tidy %>% filter(sample_type == "60-79") %>% pull(MSR_A) %>% summary()
  
}


####################################
# MSR CHANGING WITH AGE
####################################

age_stratification_MSR_changing_with_age <- function(project_id = "BRAIN") {
  
  
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################

  
  
  df_age_groups <- map_df((df_metadata$SRA_project %>% unique())[3], function(project_id) {
    
    # project_id <- (df_metadata$SRA_project %>% unique())[1]
    
    print(paste0(Sys.time(), " --> ", project_id))
    
    age_groups <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      # age_group <- age_groups[1]
      
      print(paste0(age_group))
      
      query <- paste0("SELECT ref_junID, ref_type, MSR_D, MSR_A, gene_id FROM '", age_group, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      
      query <- paste0("SELECT ref_junID, ref_type, MSR_D, MSR_A, gene_id FROM '", age_group, "_", project_id, "_misspliced'")
      introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
      
      query <- paste0("SELECT id, gene_name FROM 'gene' WHERE id IN (", paste(introns$gene_id %>% unique(), collapse =","), ")")
      introns <- merge(x = introns,
                       y = dbGetQuery(con, query) %>% as_tibble(),
                       by.x = "gene_id",
                       by.y = "id")
      
      return(introns %>%
               mutate(sample_type = age_group))
      
    })
    
  })
  
  
  ###############################################
  ## QC
  ###############################################
  
  any(df_age_groups %>% filter(sample_type == "never") %>% pull(MSR_D) > 0)
  any(df_age_groups %>% filter(sample_type == "never") %>% pull(MSR_A) > 0)
  
  intersect(df_age_groups %>% filter(sample_type == "never") %>% pull(ref_junID),
            df_age_groups %>% filter(sample_type != "never") %>% pull(ref_junID))
  
  ###############################################
  ## GET ONLY COMMON JUNCTIONS ACROSS AGE GROUPS
  ###############################################
  
  common_introns <- df_age_groups %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  
  common_introns %>% length()
  
  ## Filter the INTRONS table by the common mis-spliced introns
  df_age_groups_tidy <- merge(x = df_age_groups %>% data.table::as.data.table(),
                              y = data.table::data.table(ref_junID = common_introns),
                              by = "ref_junID",
                              all.y = T) %>%
    distinct(ref_junID, sample_type, .keep_all = T)
  
  
  
  ## MSR_D -------------------------------------------------
  df_MSRD <- df_age_groups_tidy %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_D,
                  gene_name) %>%
    mutate(MSR_D = MSR_D %>% round(digits = 4)) %>%
    spread(sample_type, MSR_D)
  
  # genes_MSRD_discard <- df_MSRD %>%
  #   rowwise() %>%
  #   filter(`20-39` > `40-59` |
  #            `20-39` > `60-79`) %>%
  #   distinct(gene_name) %>% 
  #   pull() %>% unlist()
  
  genes_MSRD_increasing <- df_MSRD %>%
    filter((`20-39` < `40-59` &
              `40-59` < `60-79`) | (`20-39` < `40-59` &
                                      `20-39` < `60-79`))
  genes_MSRD_decreasing <- df_MSRD %>%
    filter((`20-39` > `40-59` & 
              `40-59` > `60-79`) | (`20-39` > `40-59` & 
                                      `20-39` > `60-79`))
  
  
  saveRDS(object = genes_MSRD_increasing %>% distinct(ref_junID, .keep_all = T),
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRD.rds")
  saveRDS(object = genes_MSRD_decreasing %>% distinct(ref_junID, .keep_all = T),
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_decrease_MSRD.rds")
  
  ## MSR_A -------------------------------------------------
  df_MSRA <- df_age_groups_tidy %>%
    #filter(ref_junID %in% df_age_groups_novel_tidy$ref_junID) %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_A,
                  gene_name) %>%
    mutate(MSR_A = MSR_A %>% round(digits = 4)) %>%
    spread(sample_type, MSR_A)
  
  
  genes_MSRA_increasing <- df_MSRA %>%
    filter((`20-39` < `40-59` & 
              `40-59` < `60-79`) | (`20-39` < `40-59` & 
                                      `20-39` < `60-79`))
  genes_MSRA_decreasing <- df_MSRA %>%
    filter((`20-39` > `40-59` &
              `40-59` > `60-79`) | (`20-39` > `40-59` &
                                      `20-39` > `60-79`))
  
  
  saveRDS(object = genes_MSRA_increasing %>% distinct(ref_junID, .keep_all = T),
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRA.rds")
  saveRDS(object = genes_MSRA_decreasing %>% distinct(ref_junID, .keep_all = T),
          file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_decrease_MSRA.rds")
  
  
  
  
  #######################################
  ## STATISTICAL TEST - MSR_DONOR
  ########################################
  
  wilcox.test(x = df_MSRD$`20-39`,
              y = df_MSRD$`40-59`,
              alternative = "less",
              paired = T)
  wilcox.test(x = df_MSRD$`40-59`,
              y = df_MSRD$`60-79`,
              alternative = "less",
              paired = T)
  
  t.test(x = df_MSRD$`20-39`,
         y = df_MSRD$`40-59`,
         alternative = "less",
         paired = T)
  t.test(x = df_MSRD$`20-39`,
         y = df_MSRD$`60-79`,
         alternative = "less",
         paired = T)
  
  wilcox.test(x = df_MSRA$`20-39`,
              y = df_MSRA$`40-59`,
              alternative = "less",
              paired = T)
  wilcox.test(x = df_MSRA$`40-59`,
              y = df_MSRA$`60-79`,
              alternative = "less",
              paired = T)
}

####################################
# SPLICEOSOMAL RBPs AND THEIR TPMs
####################################

age_stratification_RBP_TPM <- function (project_id = "BRAIN") {
  
  ################################################################################
  # 0. Prepare the necessary data ------------------------------------------------
  ################################################################################
  
  
  project_id <- "BRAIN"
  
  ## Load reference gtf and get only the genes
  ensembl105 <- biomaRt::useEnsembl(biomart = 'genes', 
                                    dataset = 'hsapiens_gene_ensembl',
                                    version = 105)
  
  
  ## Load the RBPs and add the ensemblID
  all_RBPs <- xlsx::read.xlsx(file = '/home/sruiz/PROJECTS/splicing-project-recount3/markdowns/41586_2020_2077_MOESM3_ESM.xlsx', 
                              sep = '\t', header = TRUE,
                              sheetIndex = 1) %>%  as_tibble() 
  
  RBPs_annotated <- biomaRt::getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
    filters = 'hgnc_symbol',
    values = all_RBPs$name %>% unique,
    mart = ensembl105
  )
  RBPs_annotated <- RBPs_annotated %>%
    distinct(hgnc_symbol, .keep_all = T)
  
  
  ## Load the recount3 data
  rse <- recount3::create_rse_manual(
    project = project_id,
    project_home = "data_sources/gtex",
    organism = "human",
    annotation = "gencode_v29",
    type = "gene")
  
  ## Transform counts
  SummarizedExperiment::assays(rse)$counts <- recount3::transform_counts(rse)
  ## See that now we have two assayNames()
  SummarizedExperiment::assayNames(rse)
  
  
  ################################################################################
  # 1. Get the RBP TPM expression using GTEX data --------------------------------
  ################################################################################
  
  dds_tpm <- get_GTEx_gene_expression(rse = rse, ensembl = ensembl105)
  dds_tpm %>% head()
  
  dds_tpm <- dds_tpm %>%
    filter(tpm > 0)
  
  # Tidy the 'dds_tpm' object
  dds_tpm_pv <- dds_tpm %>%
    dplyr::filter(gene %in% RBPs_annotated$ensembl_gene_id) %>%
    dplyr::select(gene, sample, tpm) %>% 
    group_by(gene) %>%
    distinct(sample, .keep_all = T) %>%
    pivot_wider(names_from = sample, values_from = tpm, values_fill = 0)
  
  write.csv(x = dds_tpm_pv %>%
              tibble::column_to_rownames("gene"),
            file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/results/pipeline3/rbp/tpm_ensembl105.csv"))
}

age_stratification_RBP_TPM_LM <- function(projects_id = c("BRAIN")) {
  
  ################################
  ## LOAD RBPs OF INTEREST
  ################################
  
  all_RBPs <- xlsx::read.xlsx(file = '/home/sruiz/PROJECTS/splicing-project-recount3/markdowns/41586_2020_2077_MOESM3_ESM.xlsx', 
                              sep = '\t', header = TRUE,
                              sheetIndex = 1) %>%  
    as_tibble() %>%
    mutate(type = ifelse(Splicing.regulation == 1, "splicing regulation", "other")) %>%
    dplyr::select(name, ensembl_gene_id = id, type) %>%
    distinct(name, .keep_all = T) 
  
  ## Load reference gtf and get only the genes
  # ensembl105 <- biomaRt::useEnsembl(biomart = 'genes', 
  #                                   dataset = 'hsapiens_gene_ensembl',
  #                                   version = 105, mirror = "asia")
  ensembl105 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf") %>%
    as_tibble()
  # RBPs_annotated <- biomaRt::getBM(
  #   attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  #   filters = 'hgnc_symbol',
  #   values = all_RBPs$name %>% unique,
  #   mart = ensembl105
  # )
  RBPs_annotated <- ensembl105 %>%
    as_tibble() %>%
    filter(gene_id %in% all_RBPs$ensembl_gene_id) %>%
    distinct(gene_name, .keep_all = T)
  
  
  RBPs_annotated_tidy <- left_join(x = RBPs_annotated %>% select(gene_id, gene_name) %>% as_tibble(), 
                                   y = all_RBPs %>% as_tibble(), 
                                   by = c("gene_id" = "ensembl_gene_id")) 
  
  
  
  
  for (project_id in projects_id) {
    
    ## project_id <- "BRAIN"
    
    ## Load and tidy the uncorrected TPMs
    ## TPMs here should be uncorrected as the covariates are going to be included in the linear models
    rbp_recount_tpm <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                              project_id, "/results/pipeline3/rbp/tpm_ensembl105.csv"))
    
    
    tpm_uncorrected <- rbp_recount_tpm %>%
      dplyr::rename(gene = "X") %>%
      as_tibble() %>%
      distinct(gene, .keep_all = T)
    
    
    ## Load and tidy the sample metadata
    sample_metadata <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                              project_id, "/results/pipeline3/rbp/covariates.csv"), header = T,
                                fileEncoding = "UTF-8") %>% as_tibble()
    
    print(sample_metadata %>% nrow)
    ## Check the samples are the same
    if ( !( identical(x = names(tpm_uncorrected)[-1],
                      y = names(sample_metadata)[-1]) ) ) {
      print("ERROR!")
      break;
    }
    
    
    
    RBPs <- (tpm_uncorrected$gene)
    
    # tpm_uncorrected[rowSums(tpm_uncorrected[, c(2:ncol(tpm_uncorrected))]),]
    
    ## Obtain values and run lm
    df_lm_output <- map_df(RBPs, function(RBP) {
      
      # RBP <- RBPs[1]
      # RBP <- "ENSG00000277564"
      print(RBP)
      
      tpm <- tpm_uncorrected %>%
        filter(gene == RBP) %>%
        mutate(gene = "tpm") %>%
        as_tibble()
      
      if (rowSums(tpm[, c(2:ncol(tpm))]) != 0) {
        
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
        
        if (project_id == "MUSCLE") {
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
                            #gtex.smtsd +
                            gtex.age,
                          data = df)
        } else {
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
        }
        
        
        lm_output <- lm_output %>% summary()
        df_lm_output <- lm_output$coefficients %>% as_tibble(rownames = "covariate")
        
        # cor(df)
        df_lm_output <- df_lm_output %>%
          mutate(RBP_ID = RBP) %>%
          dplyr::select(covariate, Estimate, pval = `Pr(>|t|)`, RBP_ID)
        
        return(df_lm_output)
        
      } else {
        
        return(NULL)
        
      }
      
    }) 
    
    
    # any(df_lm_output$pval > 0.05)
    
    
    ## Filter by age covariate and order
    df_lm_output_age <- df_lm_output %>%
      group_by(covariate) %>%
      mutate(pval_corrected = p.adjust(pval)) %>%
      ungroup() %>%
      filter(str_detect(string = covariate, pattern = "gtex.age")) %>%
      arrange(Estimate)
    
    
    
    
    ## Add gene SYMBOL info
    df_lm_age_tidy <- left_join(x = df_lm_output_age, 
                                y = RBPs_annotated_tidy %>% drop_na(), 
                                by = c("RBP_ID" = "gene_id")) %>%
      group_by(RBP_ID) #%>%
    #distinct(pval, .keep_all = T)
    
    
    write_csv(x = df_lm_age_tidy,
              file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                            "/results/pipeline3/rbp/tpm_lm_all_ensembl105.csv"))
    
  }
}

age_stratification_RBPs_affected_age <- function (project_id = "BRAIN") {
  
  ####################################
  # LOAD RBPs ORIGINAL PAPER AND TIDY
  ####################################
  
  all_RBPs_tidy <- xlsx::read.xlsx(file = '/home/sruiz/PROJECTS/splicing-project-recount3/markdowns/41586_2020_2077_MOESM3_ESM.xlsx', 
                                   sep = '\t', header = TRUE,
                                   sheetIndex = 1) %>%  
    as_tibble() %>%
    mutate(type = ifelse(Splicing.regulation == 1, "splicing regulation", "other")) %>%
    dplyr::select(name, ensembl_gene_id = id, type) %>%
    distinct(name, .keep_all = T) %>%
    dplyr::rename(RBP_name = name,
                  RBP_type = type)
  
  
  write.table(x = all_RBPs_tidy$RBP_name,
              file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs.csv"),
              row.names = F, col.names = F, quote = F )
  
  ####################################
  # LOAD RBPs LM RESULTS
  ####################################
  
  RBPs_age_lm <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                                        "/results/pipeline3/rbp/tpm_lm_all_ensembl105.csv")) %>%
    drop_na() %>% 
    as_tibble()
  
  
  
  ##########################################
  # CLASSIFY AFFECTED / NOT AFFECTED BY AGE
  ##########################################
  
  RBP_age <- RBPs_age_lm %>%
    filter(pval_corrected  <= 0.05,
           Estimate < 0) %>%
    drop_na() 
  
  write.table(x = RBP_age$name,
              file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_decrease.csv"),
              row.names = F, col.names = F, quote = F )
  
  RBP_notage <- RBPs_age_lm %>%
    filter(!(name %in% RBP_age$name    )) %>%
    drop_na()
  
  write.table(x = RBP_notage$name,
              file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_other.csv"),
              row.names = F, col.names = F,quote = F)
  
  
  plot(density(RBP_notage$Estimate))
  lines(density(RBP_age$Estimate), col = "red")
  
  RBPs_age_tidy <- RBPs_age_lm %>%
    drop_na() %>%
    mutate(age_type = ifelse(pval_corrected <= 0.05,
                             "affected_age",
                             "not_affected_age"))
  RBPs_age_tidy %>%
    dplyr::count(age_type)
  
  write.csv(x = RBPs_age_tidy,
            file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id,
                          "/results/pipeline3/rbp/tpm_lm_all_ensembl105_tidy.csv"))
}


####################################
# ENCORI
####################################

overlap_ENCORI_introns_MSR <- function(project_id = "BRAIN") {
  
  
  con <- dbConnect(RSQLite::SQLite(), "./database/age-stratification.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  chain <- rtracklayer::import.chain(con = "data/hg19ToHg38.over.chain")
  
  
  ###########################
  ## TIDY INTRONS
  ###########################
  
  query <- paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query)
  
  indx <- str_locate_all(string = df_intron$ref_coordinates, pattern = ":")
  chr_positions <- NULL
  start_positions <- NULL
  end_positions <- NULL
  strand_positions <- NULL
  
  for (i in c(1:(indx %>% length()))) {
    
    chr_j <- df_intron$ref_coordinates[i] %>%
      str_sub(start = 1,
              end = str_locate_all(string = df_intron$ref_coordinates[i], pattern = ":")[[1]][1,2]-1)
    start_j <- df_intron$ref_coordinates[i] %>%
      str_sub(start = str_locate_all(string = df_intron$ref_coordinates[i], pattern = ":")[[1]][1,2]+1,
              end = str_locate_all(string = df_intron$ref_coordinates[i], pattern = "-")[[1]][1,2]-1) %>% as.integer()
    end_j <- df_intron$ref_coordinates[i] %>%
      str_sub(start = str_locate_all(string = df_intron$ref_coordinates[i], pattern = "-")[[1]][1,2]+1,
              end = str_locate_all(string = df_intron$ref_coordinates[i], pattern = ":")[[1]][2,2]-1) %>% as.integer()
    strand_j <- df_intron$ref_coordinates[i] %>%
      str_sub(start = str_locate_all(string = df_intron$ref_coordinates[i], pattern = ":")[[1]][2,2]+1,
              end = df_intron$ref_coordinates[i] %>% stringr::str_count())
    
    chr_positions <- c(chr_positions, chr_j)
    start_positions <- c(start_positions, start_j)
    end_positions <- c(end_positions, end_j)
    strand_positions <- c(strand_positions, strand_j)
    
    print(paste0(i, " out of ", indx %>% length()))
  }
  
  df_intron_tidy <- df_intron %>%
    mutate(start = start_positions,
           end = end_positions,
           strand = strand_positions,
           seqnames = chr_positions)
  
  ###########################
  ## MSR_D
  ###########################
  
  genes_MSRD_increasing <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRD.rds")
  genes_MSRD_increasing <- merge(genes_MSRD_increasing,
                                 y = df_intron_tidy %>% select(ref_junID, start, end, seqnames, strand),
                                 by.x = "ref_junID",
                                 by.y = "ref_junID")
  genes_MSRD_increasing <- genes_MSRD_increasing %>%
    mutate(start = start - 25,
           end = end + 25) %>% 
    GRanges()
  
  
  genes_MSRD_decreasing <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_decrease_MSRD.rds")
  genes_MSRD_decreasing <- merge(genes_MSRD_decreasing,
                                 y = df_intron_tidy %>% select(ref_junID, start, end, seqnames, strand),
                                 by.x = "ref_junID",
                                 by.y = "ref_junID")
  genes_MSRD_decreasing <- genes_MSRD_decreasing %>%
    mutate(start = start - 25,
           end = end + 25) %>% 
    GRanges()
  
  
  
  ###########################
  ## MSR_A
  ###########################
  
  genes_MSRA_increasing <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_increase_MSRA.rds")
  genes_MSRA_increasing <- merge(genes_MSRA_increasing,
                                 y = df_intron_tidy %>% select(ref_junID, start, end, seqnames, strand),
                                 by.x = "ref_junID",
                                 by.y = "ref_junID")
  
  genes_MSRA_increasing <- genes_MSRA_increasing %>%
    mutate(start = start - 25,
           end = start + 25) %>% 
    GRanges()
  
  genes_MSRA_decreasing <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/genes_decrease_MSRA.rds")
  genes_MSRA_decreasing <- merge(genes_MSRA_decreasing,
                                 y = df_intron_tidy %>% select(ref_junID, start, end, seqnames, strand),
                                 by.x = "ref_junID",
                                 by.y = "ref_junID")
  
  genes_MSRA_decreasing <- genes_MSRA_decreasing %>%
    mutate(start = start - 25,
           end = start + 25) %>% 
    GRanges()
  
  ######################################################################
  ## GET OVERLAPS - ENCORI AND INTRONS WITH INCREASING MSR_D AND MSR_A
  ######################################################################
  
  df_RBP_ENCORI_MSR_result <- map_df(c("RBPs_affected_age", "RBPs_notaffected_age"), function(type) {
    
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
    
    ##############################
    ## Check ENCORI iCLIP results
    ##############################
    
    map_df(RBPs$V1 %>% sort(), function (RBP) {
      
      # RBP <- "ADAR"
      file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/results/ENCORI/",
                          "/ENCORI_hg19_", RBP, "_allgenes.txt")
      
      #print(file_name)
      if (file.exists(file_name)) {
        
        
        ENCORI_RBP_result <- read.delim(file = file_name, header = T, skip = 3, sep = "\t") %>%
          as_tibble()
        
        if (ENCORI_RBP_result %>% nrow() > 1) {
          
          print(paste0(RBP, " - ", type))
          
          ENCORI_RBP_result <- ENCORI_RBP_result %>% 
            dplyr::rename(start = broadStart, end = broadEnd) %>% 
            GRanges()
          
          #######################
          # Liftover
          #######################
          
          # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
          
          ENCORI_RBP_result_GRh38 <- rtracklayer::liftOver(x = ENCORI_RBP_result, 
                                                           chain = chain) %>% 
            unlist() 
          
          #######################
          ## Overlaps - MSR_D
          #######################
          
          overlaps_MSRD <- GenomicRanges::findOverlaps(query = genes_MSRD_increasing,
                                                       subject = ENCORI_RBP_result_GRh38,
                                                       type = "any",
                                                       ignore.strand = F)
          
          overlaps_MSRDd <- GenomicRanges::findOverlaps(query = genes_MSRD_decreasing,
                                                        subject = ENCORI_RBP_result_GRh38,
                                                        type = "any",
                                                        ignore.strand = F)
          
          
          #######################
          ## Overlaps - MSR_A
          #######################
          
          overlaps_MSRA <- GenomicRanges::findOverlaps(query = genes_MSRA_increasing,
                                                       subject = ENCORI_RBP_result_GRh38,
                                                       type = "any",
                                                       ignore.strand = F)
          
          overlaps_MSRAd <- GenomicRanges::findOverlaps(query = genes_MSRA_decreasing,
                                                        subject = ENCORI_RBP_result_GRh38,
                                                        type = "any",
                                                        ignore.strand = F)
          
          #######################
          ## Return result
          #######################
          
          
          data.frame(RBP_name = RBP,
                     RBP_type = type,
                     ovlps_MSRD = genes_MSRD_increasing[queryHits(overlaps_MSRD),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRA = genes_MSRA_increasing[queryHits(overlaps_MSRA),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRD = genes_MSRD_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRA = genes_MSRA_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRD_perc = (((genes_MSRD_increasing[queryHits(overlaps_MSRD),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                          genes_MSRD_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow()),
                     ovlps_MSRA_perc = (((genes_MSRA_increasing[queryHits(overlaps_MSRA),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                          genes_MSRA_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow()),
                     
                     ovlps_MSRDd = genes_MSRD_decreasing[queryHits(overlaps_MSRDd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRAd = genes_MSRA_decreasing[queryHits(overlaps_MSRAd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRDd = genes_MSRD_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRAd = genes_MSRA_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRD_percd = (((genes_MSRD_decreasing[queryHits(overlaps_MSRDd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                           genes_MSRD_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow()),
                     ovlps_MSRA_percd = (((genes_MSRA_decreasing[queryHits(overlaps_MSRAd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                           genes_MSRA_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow())
          ) %>% return()
        }
      }
      
    })
    
  })
  
  write_csv(x = df_RBP_ENCORI_MSR_result,
            file = "paper_figures/iCLIP/results/RBPs_ENCORI_MSR.csv", col_names = T)
  
}

ENCORI_test_results <- function () {
  
  
  df_RBP_result <- read.csv(file = "paper_figures/iCLIP/results/RBPs_ENCORI_MSR.csv") %>% as_tibble()
  df_RBP_result
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRA_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRA_perc),
              paired = F,
              alternative = "greater")
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRD_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRD_perc),
              paired = F,
              alternative = "greater")
  
  
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRD_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRD_percd),
              paired = T,
              alternative = "greater")
  
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRA_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRA_percd),
              paired = T,
              alternative = "less")
  
  
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRD_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRD_percd),
              paired = T,
              alternative = "less")
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRA_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRA_percd),
              paired = T,
              alternative = "less")
  
  
  
  
}

################################
## CALLS
################################




for (project_id in c("BLOOD", "MUSCLE", "BRAIN")) {
  
  print(paste0(Sys.time(), " --> ", project_id))
  # project_id <- "BLOOD"
  project_init <- age_stratification_init_data(project_id)
  
  age_stratification_IDB(age_groups = project_init$age_group %>% unique(),
                         project_id = project_id,
                         age_samples_clusters = project_init)
  
}
