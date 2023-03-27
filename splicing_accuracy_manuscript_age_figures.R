library(tidyverse)
library(DBI)

# source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline6_splicing_paper_age.R")


setwd(normalizePath("."))


##################################
## CONNECT TO THE DATABASE
##################################

gtf_version <- 105
main_project <- "age_subsampled"
database_path <- paste0(getwd(), "/database/v", gtf_version, "/", main_project, "/", main_project, ".sqlite")
con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)

age_projects <- readRDS(file = paste0(getwd(), "/results/",main_project,"_final_projects_used.rds"))

custom_ggtheme <-  theme(text = element_text(size = 7, family="Arial", colour = "black"),
                         
                         axis.ticks = element_line(colour = "black", linewidth = 2),
                         axis.text = element_text(size = 7, family="Arial", colour = "black"),
                         axis.line = element_line(colour = "black"),
                         axis.title = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.y = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.x = element_text(size = 7, family="Arial", colour = "black", 
                                                    hjust = 0.5, vjust = 0.5),
                         strip.text = element_text(size = 7, family="Arial", colour = "black"),
                         
                         legend.text = element_text(size = "7", family="Arial", colour = "black"),
                         legend.title = element_blank(),
                         legend.position = "top",
                         legend.box = "vertical")


##################################
## STATS
##################################

age_stratification_get_stats <- function() {
  
  
  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'master'")
  
  ## Total number of samples considered
  db_metadata <- dbGetQuery(con, query) %>% as_tibble()
  db_metadata %>% nrow() 
  
  ## Number of samples per age group
  db_metadata %>%
    dplyr::count(cluster)
  
  ## Number of tissues considered
  db_metadata %>%
    dplyr::count(SRA_project)
  
  
  db_metadata %>%
    group_by(SRA_project) %>%
    distinct(cluster, .keep_all = T) %>%
    dplyr::select(SRA_project,cluster) %>%
    #print(n = 60) %>%
    write.csv(file = paste0(getwd(), "/results/_paper/results/", main_project, "_body_sites_tissues.csv"), row.names = FALSE)
  
  db_metadata %>%
    group_by(SRA_project) %>%
    mutate(nsamples = n()) %>%
    filter(nsamples >= 70) %>%
    distinct(SRA_project, .keep_all = T)
  
  db_metadata %>%
    group_by(SRA_project) %>%
    mutate(nsamples = n()) %>%
    filter(nsamples >= 70) %>%
    distinct(SRA_project, .keep_all = T) %>%
    pull(nsamples) %>% sum
  
  db_metadata %>%
    group_by(SRA_project) %>%
    mutate(nsamples = n()) %>%
    filter(nsamples >= 70) %>%
    distinct(cluster, .keep_all = T) 
  
  
  
  query <- paste0("SELECT * from 'intron'")
  db_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  ## This database contained a total of 308,717 annotated introns
  db_introns %>% distinct(ref_junID) %>% nrow()
  
  ## From which 228,534 presented evidence of at least one type of mis-splicing event. 
  db_introns %>% dplyr::count(misspliced)
  

  
  query <- paste0("SELECT * from 'novel'")
  novel <- dbGetQuery(con, query) 
  ## It also included 719,069 novel donor and 999,041 novel acceptor junctions
  novel %>% distinct(novel_junID) %>% nrow()
  novel %>% dplyr::count(novel_type)
  
  ## covering 199,1991 transcripts, 30,580 genes and 6,111 samples from 40 tissues and 18 body sites. 
  query <- paste0("SELECT * from 'gene'")
  gene <- dbGetQuery(con, query) 
  gene %>% nrow()
  
  query <- paste0("SELECT * from 'transcript'")
  transcript <- dbGetQuery(con, query) 
  transcript %>% distinct(transcript_id) %>% nrow()

}


age_stratification_plot_metadata <- function() {
  
  custom_theme <- theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "9"),
          axis.text.x = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "9"),
          strip.text = element_text(colour = "black", size = "9"), 
          legend.text = element_text(colour = "black", size = "9"),
          plot.caption = element_text(colour = "black", size = "9"),
          plot.title = element_text(colour = "black", size = "9"),
          legend.title = element_text(colour = "black", size = "9"),
          legend.position = "top",
          ## This is to plot the two legends in two rows
          legend.box="vertical")  
  
  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'master'")
  db_metadata <- dbGetQuery(con, query) %>% 
    as_tibble() %>% 
    mutate(gender = gender %>% as.integer())
  
  # db_metadata <- project_init %>% dplyr::rename(cluster = age_group, SRA_project = project, gender = sex)

  ## Num samples
  plot_num_samples <- ggplot(db_metadata %>% 
           dplyr::count(SRA_project,cluster) %>%
           arrange(cluster , n) %>%
           mutate(SRA_project = fct_inorder(SRA_project)) ) +
    geom_bar(aes(y = SRA_project, x = n, fill = cluster),
             stat = "identity", position = position_dodge()) + 
    theme_light() +
    scale_fill_hue() +
    labs(y = ("Body site"), x = "Num. samples" ) +
    guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_theme
  
  plot_num_samples
  
  
  ## Gender
  plot_gender <- ggplot(db_metadata %>% 
           mutate(gender = gender %>% as.character()) %>%
           dplyr::count(SRA_project, gender) %>%
           arrange(gender , n) %>%
           mutate(SRA_project = fct_inorder(SRA_project)) ) +
    geom_bar(aes(y = SRA_project, x = n, group = gender, fill = gender),
             stat = "identity", position = "dodge") + 
    theme_light() +
    labs(y = ("Body site"), x = "Num. samples" ) +
    scale_fill_manual(labels = c("Male", "Female"), values = c("1", "2"), palette=scales::hue_pal()) +
    guides(fill = guide_legend(title = "Gender: ", ncol = 2, nrow = 1)) +
    custom_theme
  
  plot_gender
  
  
  ## RIN
  plot_rin <- ggplot(db_metadata ) +
    geom_density(aes(x = rin, fill = cluster), alpha = 0.5) + 
    theme_light() +
    labs(x = "RIN" ) +
    scale_fill_hue() +
    guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_theme
  
  plot_rin
  
  
  ## Mapped read depth
  plot_read_depth <- ggplot(db_metadata ) +
    geom_density(aes(x = mapped_read_count, fill = cluster), alpha = 0.5) + 
    theme_light() +
    scale_fill_hue() +
    labs(x = "Mapped Read Count" ) +
    guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_theme
  
  plot_read_depth
  
  ggpubr::ggarrange(plot_num_samples,plot_gender,
                    plot_rin,plot_read_depth,
                    labels = c("a","b","c","d"))
  
  
  figures_path <- paste0(getwd(), "/results/_paper/figures/")
  file_name <- paste0(figures_path, "/age_metadata")
  ggplot2::ggsave(paste0(file_name, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(paste0(file_name, ".png"), width = 210, height = 190, units = "mm", dpi = 300)
  
  
}




## 9. Increasing age is associated with increasing levels of splicing errors


get_common_introns_across_age_groups <- function() {
  
  ## CONNECT TO THE DATABASE
  con <- dbConnect(RSQLite::SQLite(), database_path) 
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  
  
  ## Load metadata
  df_metadata <- dbGetQuery(con, query) %>%
    group_by(SRA_project) %>%
    mutate(nsamples = n()) %>%
    filter(nsamples >= 70)
  
  
  folder_results <- paste0(getwd(), "/results/_paper/results/")
  dir.create(file.path(folder_results), recursive = TRUE, showWarnings = T)
  
  
  df_age_groups_all_introns <- map_df(age_projects, function(project_id) {
      
    # project_id <- age_projects[1]
    
    print(paste0(Sys.time(), " --> ", project_id))
    
    ## Get the age clusters
    age_groups <- df_metadata %>%
      filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull(cluster)
      
    
    df_age_groups <- map_df(age_groups, function(age_group) {
      # age_group <- age_groups[1]
      
      print(paste0(age_group))
      
      query <- paste0("SELECT name 
                      FROM sqlite_master 
                      WHERE type='table' AND name='", 
                      age_group, "_", project_id, "_nevermisspliced';")
        
      if ( (dbGetQuery(con, query) %>% nrow()) > 0 ) {
        
        ###########################################
        ## GET DATA FROM THE DATABASE
        
        query <- paste0("SELECT ref_junID, ref_type, MSR_D, MSR_A, transcript_id, 
                        novel_junID, ref_sum_counts, novel_sum_counts 
                        FROM '", age_group, "_", project_id, "_misspliced'")
        introns <- dbGetQuery(con, query) %>% as_tibble()
        
        
        ## Merge with 'master' novel table
        query <- paste0("SELECT novel_junID, novel_type 
                        FROM 'novel' 
                        WHERE novel_junID IN (", 
                        paste(introns$novel_junID %>% unique(), collapse =","), ")")
        introns <- introns %>%
          inner_join(y = dbGetQuery(con, query) %>% as_tibble(),  by = "novel_junID")
        
        
        ## Add never-misspliced info
        query <- paste0("SELECT ref_junID, ref_type, MSR_D, MSR_A, transcript_id 
                        FROM '", 
                        age_group, "_", project_id, "_nevermisspliced'")
        introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
        
        
        ## Add transcript info
        query <- paste0("SELECT id, gene_id FROM 'transcript' WHERE id IN (", 
                        paste(introns$transcript_id %>% unique(), collapse =","), ")")
        introns %>% nrow()
        introns <- introns %>%
          inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                     by = c("transcript_id" = "id"))
        
        
        ## Add gene info
        query <- paste0("SELECT id, gene_name FROM 'gene' WHERE id IN (", 
                        paste(introns$gene_id %>% unique(), collapse =","), ")")
        introns <- introns %>%
          inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                     by = c("gene_id" = "id"))
        
        return(introns %>% 
                 distinct(ref_junID, .keep_all = T) %>%
                 dplyr::select(-transcript_id, -gene_id) %>%
                 mutate(sample_type = age_group,
                        project_id = project_id))
      } else {
        return(NULL)
      }
    })
    
  })
    
  ## Get the age clusters
  age_groups <- df_metadata$cluster %>% unique()
  age_projects <- df_age_groups_all_introns$project_id %>% unique()
  
  
  ## Get common introns across age groups per tissue
  df_common_introns <- df_age_groups_all_introns %>%
    as_tibble() %>%
    group_by(project_id, sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    group_by(ref_junID) %>%
    dplyr::count() %>%
    filter(n == (age_groups %>% length()) * (age_projects %>% length())) %>%
    ungroup()
  df_common_introns$ref_junID %>% unique() %>% length()
    
  
  ## Filter by common introns
  df_age_groups_tidy <- df_age_groups_all_introns %>% 
    inner_join(y = df_common_introns %>% distinct(ref_junID),
               by = "ref_junID") %>%
    group_by(project_id, sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup()
  df_age_groups_tidy %>% distinct(ref_junID) %>% nrow()
  
  
  ## Save data
  file_name <- paste0(folder_results, "/", main_project, "_common_introns_all_age_groups.rds")
  saveRDS(object = df_age_groups_tidy, file = file_name)
  
}


get_effsize_MSR_with_age <- function() {
  

  folder_results <- paste0(getwd(), "/results/_paper/results/")
  dir.create(file.path(folder_results), recursive = TRUE, showWarnings = T)
  
  ## Load common introns across age groups and tissues
  file_name <- paste0(folder_results, "/", main_project, "_common_introns_all_age_groups.rds")
  df_age_groups_tidy <- readRDS(file = file_name)
  
  
  ## MSR_D data preparation
  df_MSRD <- df_age_groups_tidy %>%
    dplyr::select(ref_junID, sample_type, MSR_D, gene_name, project_id)  %>% 
    #distinct(project_id, sample_type, ref_junID, .keep_all = T) %>%
    mutate(MSR_D = MSR_D %>% round(digits = 4)) %>%
    spread(sample_type, MSR_D) %>%
    as_tibble()
    
  
  ## MSR_A data preparation
  df_MSRA <- df_age_groups_tidy %>%
    dplyr::select(ref_junID, sample_type, MSR_A, gene_name, project_id)  %>% 
    #distinct(project_id, sample_type, ref_junID, .keep_all = T) %>%
    mutate(MSR_A = MSR_A %>% round(digits = 4)) %>%
    spread(sample_type, MSR_A) %>%
    as_tibble()
    
    
  df_wilcoxon <- map_df( c("MSR Donor", "MSR Acceptor"), function(MSR_type) { #, "MSR_A"
      
    # MSR_type <- "MSR Donor"
    print( paste0( "Getting data from: '", MSR_type, "'..." ) )

    map_df(age_projects, function(proj_id) {
      
      # proj_id <- age_projects[2]
      # proj_id <- "KIDNEY"
      print( paste0( proj_id ) )
      
      ###############################
      ## MSR_D INCREASING WITH AGE
      ###############################
      
      df_introns <- NULL
      
      if (MSR_type == "MSR Donor") {
        df_introns <- df_MSRD 
      } else {
        df_introns <- df_MSRA
      }
      
      
      df_introns <- df_introns %>%
        filter(project_id == proj_id) 
      
      if ( nrow(df_introns) > 0 ) {
        
        ## Wilcoxon text
        wilcox_pval <- data.frame(pval = c((wilcox.test(x = df_introns$`20-39`,
                                                        y = df_introns$`60-79`,
                                                        alternative = "less",
                                                        paired = T, 
                                                        conf.int = T))$p.value))
        
        # Zstat <- qnorm(wilcox_pval$pval/2)
        # abs(Zstat)/sqrt(df_introns %>% nrow())
        
        ## Wilcoxon effect size
        r_effect1 <- rstatix::wilcox_effsize(data = df_introns %>%
                                               dplyr::select(ref_junID, `20-39`,`60-79`) %>%
                                               gather(key = age, value = MSR, -ref_junID) %>%
                                               mutate(age = age %>% as.factor()),
                                             formula = MSR ~ age,
                                             paired = T,
                                             p.adjust.method = "fdr") %>%
          mutate(MSR_type = paste0(MSR_type), 
                 tissue = proj_id)
        
        
        return(rbind(r_effect1) %>% cbind(wilcox_pval))
        
      } else {
        return(NULL)
      }
      
    })
    
  })
  
  file_name <- paste0(folder_results, "/", main_project, "_effsize_MSR_with_age.rds")
  saveRDS(object = df_wilcoxon, file = file_name)
  
}


plot_effsize_MSR_with_age <- function() {
  
  
  folder_results <- paste0(getwd(), "/results/_paper/results/")
  file_name <- paste0(folder_results, "/", main_project, "_df_wilcoxon_effsize_paired_all_projects.rds")
  df_wilcoxon <- readRDS(file = file_name)
  
  #######################################################
  ## PLOTS
  #######################################################
  
  df_wilcoxon_tidy <- df_wilcoxon %>%
    group_by(MSR_type) %>%
    mutate(q = p.adjust(pval, method = "fdr"))%>%
    ungroup()
  
  df_wilcoxon_tidy$MSR_type = factor(df_wilcoxon_tidy$MSR_type, 
                                     levels = c("MSR Donor",
                                                "MSR Acceptor"))
  
  df_wilcoxon_tidy_final <- df_wilcoxon_tidy %>%
    filter(tissue %in% (df_wilcoxon_tidy %>%
                          group_by(tissue) %>%
                          filter(pval <= 0.05) %>%
                          dplyr::count(pval) %>%
                          #filter(q < 0.05) %>%
                          ungroup() %>% 
                          pull(tissue)) ) %>%
    mutate(tissue = str_replace(tissue, 
                                pattern = "_",
                                replacement = " ")) %>%
    group_by(MSR_type) %>%
    mutate(tissue = fct_reorder(tissue, plyr::desc(effsize))) %>%
    ungroup()  
  
  
  
  write.csv(x = df_wilcoxon_tidy,
            file = paste0(getwd(), "/results/_paper/results/age_wilcoxon_MSR_all_tissues.csv"))
  
  
  
  
  ## PLOT 1
  ggplot(data = df_wilcoxon_tidy_final ,
         aes(x = effsize, y = tissue, color = MSR_type, size = q)) +
    geom_point(alpha=.7) +
    facet_wrap(vars(MSR_type)) +
    theme_light() +
    ylab("") +
    xlab("Probability of superior MSR in\n60-79yrs compared to 20-39yrs") +
    scale_color_manual(values = c("#35B779FF","#64037d"),
                       breaks = c("MSR Donor", "MSR Acceptor"),
                       labels = c("MSR Donor", "MSR Acceptor")) +
    custom_ggtheme + 
    scale_size(range = c(4, 1))+
    theme(legend.position="top", 
          legend.box="horizontal", 
          
          plot.margin = margin(0,0,0,0),
          legend.box.margin=margin(b = -11)) +
    guides(color = guide_legend(title = ""))+
    scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
  
  ## Save the figure 3
  folder_figures <- paste0(getwd(), "/results/_paper/figures/")
  dir.create(file.path(folder_figures), recursive = TRUE, showWarnings = T)
  
  ggplot2::ggsave(filename = paste0(folder_figures, "/panel7a.svg"), 
                  width = 180, height = 55, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(folder_figures, "/panel7a.png"), 
                  width = 180, height = 55, units = "mm", dpi = 300)
  
  
  ###################################
  ## STATS - DONOR
  ###################################
  
  df_wilcoxon_tidy_final %>%
    filter(MSR_type == "MSR Donor", 
           q < 0.05) %>%
    arrange(effsize)
  
  df_wilcoxon_tidy_final %>%
    filter(MSR_type == "MSR Donor", 
           q < 0.05) %>%
    arrange(q)
  
  
  ###################################
  ## STATS - ACCEPTOR
  ###################################
  
  df_wilcoxon_tidy_final %>%
    filter(MSR_type == "MSR Acceptor", 
           q < 0.05) %>%
    arrange(effsize)
  
  df_wilcoxon_tidy_final %>%
    filter(MSR_type == "MSR Acceptor", 
           q < 0.05) %>%
    arrange(q)
  
  
}


stats_effsize_MSR_with_age <- function() {
  
  ## CONNECT TO THE DATABASE
  con <- dbConnect(RSQLite::SQLite(), database_path) 
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  
  ## Load metadata
  df_metadata <- dbGetQuery(con, query) %>%
    group_by(SRA_project) %>%
    mutate(nsamples = n()) %>%
    filter(nsamples >= 70)
  
  ## Get the age clusters
  age_groups <- df_metadata$cluster %>% unique()
  
  
  
  #######################################
  ## STATISTICAL TEST 
  ## WILCOXON MSR
  ########################################
  
  folder_results <- paste0(getwd(), "/results/_paper/results/")
  if ( (age_projects %>% length()) > 1) {
    file_name <- paste0(folder_results, "/", main_project, "_df_wilcoxon_effsize_paired_all_projects.rds")
  } else {
    file_name <- paste0(folder_results, "/", main_project, "_df_wilcoxon_effsize_paired_", project_id,".rds")
  }
  
  if ( file.exists(file_name) ) {
    
    df_wilcoxon <- readRDS(file = file_name)
    
  } 
  
  
  
  #######################################################
  ## PAPER STATS
  #######################################################
  
  
  ## Focusing on the donor splice site, we found that MSRD values in the “60-79” age group were significantly 
  ## higher than those in the “20-39” age group in 9 of the 18 body sites analysed 
  df_wilcoxon_tidy <- df_wilcoxon %>%
    group_by(MSR_type) %>%
    mutate(pval_sig = ifelse(pval < 0.05, "yes", "no")) %>%
    #count(pval_sig) %>%
    ungroup() 
  
  
  ## (Wilcoxon signed rank test with continuity correction, 2e-16 < p-value < 1.09e-43)
  df_wilcoxon_tidy %>%
    group_by(MSR_type) %>%
    mutate(pval_sig = ifelse(pval < 0.05, "yes", "no")) %>%
    filter(pval_sig == "yes") %>%
    ungroup() %>%
    arrange(MSR_type, pval) 
  
  
  ## and that the effect of age was highest in XXXXX
  df_wilcoxon_tidy %>%
    group_by(MSR_type) %>%
    mutate(pval_sig = ifelse(pval < 0.05, "yes", "no")) %>%
    filter(pval_sig == "yes") %>%
    ungroup() %>%
    arrange(MSR_type, effsize)
  
  
  
  
  
  #######################################################
  ## OTHER STATS
  #######################################################
  
  df_wilcoxon_tidy <- df_wilcoxon_tidy %>%
    mutate(group = paste0(group1, "<", group2)) %>%
    as_tibble() %>%
    as.data.frame() 
  
  df_wilcoxon_tidy %>%
    filter(MSR_type == "MSR Donor",
           pval <= 0.05) %>%
    arrange(MSR_type, desc(pval))
  
  df_wilcoxon_tidy %>%
    filter(MSR_type == "MSR Donor",
           pval <= 0.05) %>%
    pull(effsize) %>%
    summary()
  
  
  
  df_wilcoxon_tidy %>%
    filter(MSR_type == "MSR Acceptor",
           pval <= 0.05) %>%
    arrange(MSR_type, desc(pval))
  df_wilcoxon_tidy %>%
    filter(MSR_type == "MSR Acceptor",
           pval <= 0.05) %>%
    pull(effsize) %>%
    summary()
  
  #######################################################
  ## PLOTS
  #######################################################
  
  df_wilcoxon_tidy <- df_wilcoxon_tidy %>%
    group_by(MSR_type) %>%
    mutate(q = p.adjust(pval, method = "fdr"))%>%
    ungroup()
  
  df_wilcoxon_tidy$MSR_type = factor(df_wilcoxon_tidy$MSR_type, 
                                     levels = c("MSR Donor",
                                                "MSR Acceptor"))
  
  df_wilcoxon_tidy_final <- df_wilcoxon_tidy %>%
    filter(tissue %in% (df_wilcoxon_tidy %>%
                          mutate(tissue = str_replace(tissue, 
                                                      pattern = "_",
                                                      replacement = " "))  %>%
                          group_by(tissue) %>%
                          filter(pval_sig == "yes") %>%
                          dplyr::count(pval_sig) %>%
                          #filter(q < 0.05) %>%
                          ungroup() %>% pull(tissue))) %>%
    group_by(MSR_type) %>%
    mutate(#effsize = format(effsize, digits = 2, scientific = T),
      tissue = fct_reorder(tissue, effsize)) %>%
    ungroup() 
  
  write.csv(x = df_wilcoxon_tidy,
            file = paste0(getwd(), "/results/_paper/results/age_wilcoxon_MSR_all_tissues.csv"))
  
  ## PLOT 1
  ggplot(data = df_wilcoxon_tidy_final ,
         aes(x = effsize, y = tissue, color = MSR_type, size = q)) +
    geom_point(alpha=.7) +
    facet_wrap(vars(MSR_type)) +
    theme_light() +
    ylab("") +
    xlab("Median MSR difference between\n20-39yrs vs. 60-79yrs") +
    scale_color_manual(values = c("#35B779FF","#64037d"),
                       breaks = c("MSR Donor", "MSR Acceptor"),
                       labels = c("MSR Donor", "MSR Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "10"),
          axis.text.x = element_text(colour = "black", size = "9"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          strip.text = element_text(colour = "black", size = "12")
    ) + 
    scale_size(range = c(7, 1))+
    theme(legend.position="top", 
          legend.box="vertical", 
          legend.margin=margin()) +
    guides(color = guide_legend(title = ""))+
    scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
  
  ## Save the figure 3
  folder_figures <- paste0(getwd(), "/results/_paper/figures/")
  dir.create(file.path(folder_figures), recursive = TRUE, showWarnings = T)
  
  ggplot2::ggsave(filename = paste0(folder_figures, "/effsize_common_introns_all_tissues.svg"), 
                  width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(folder_figures, "/effsize_common_introns_all_tissues.png"), 
                  width = 183, height = 100, units = "mm", dpi = 300)
  
  
  
  
}


age_stratification_brain_GO <- function() {
  
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  
  #############################################
  ## COMMON INTRONS ACROSS AGE GROUPS
  ## ONLY IN BRAIN
  #############################################
  
  project_id <- "BRAIN"
  print(paste0(Sys.time(), " --> ", project_id))
  
  age_groups <- df_metadata %>%
    filter(SRA_project == project_id) %>%
    distinct(cluster) %>%
    pull()
  
  df_age_groups <- map_df(age_groups, function(age_group) {
    # age_group <- age_groups[1]
    
    print(paste0(Sys.time(), " --> ", age_group))
    
    query <- paste0("SELECT DISTINCT ref_junID, ref_type, MSR_D, MSR_A, transcript_id 
                    FROM '", age_group, "_", project_id, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    query <- paste0("SELECT DISTINCT ref_junID, ref_type, MSR_D, MSR_A, transcript_id 
                    FROM '", age_group, "_", project_id, "_misspliced'")
    introns <- plyr::rbind.fill(introns, dbGetQuery(con, query)) %>% as_tibble()
    
    introns <- introns %>%
      distinct(ref_junID, .keep_all = T)
    
    query <- paste0("SELECT transcript.transcript_id, transcript.id, gene.gene_name, gene.gene_id
                    FROM 'transcript' INNER JOIN 'gene' ON gene.id = transcript.gene_id
                    WHERE transcript.id IN (", paste(introns$transcript_id, collapse = ","),")")
    genes <- dbGetQuery(con, query) %>% as_tibble()
    
    introns <- introns %>%
      inner_join(y = genes,
                 by = c("transcript_id" = "id")) %>%
      mutate(sample_type = age_group) %>%
      as_tibble()
    
    return(introns)
  })
 
  
  common_junctions <- df_age_groups %>%
    group_by(sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  common_junctions %>% unique() %>% length()
  
  df_age_groups_tidy <- df_age_groups %>%
    filter(ref_junID %in% common_junctions) %>%
    drop_na()
  
  
  #####################################
  ## TIDY THE DATAFRAME OF COMMON
  ## INTRONS TO GET INCREASING
  ## MSR_D AND MSR_A VALUES WITH AGE
  #####################################
  
  
  df_MSRD <- df_age_groups_tidy %>%
    dplyr::select(ref_junID, sample_type, MSR_D, gene_name, gene_id) %>%
    spread(sample_type, MSR_D)
  
  df_MSRA <- df_age_groups_tidy %>%
    dplyr::select(ref_junID, sample_type, MSR_A, gene_name, gene_id) %>%
    spread(sample_type, MSR_A)
  
  ## Introns increasing MSR with age
  df_MSRD_increasing <- df_MSRD %>%
    filter((`20-39` < `40-59` & `40-59` < `60-79`) | 
             (`20-39` == `40-59` & `40-59` < `60-79`) |
             (`20-39` < `40-59` & `40-59` == `60-79`))
  df_MSRA_increasing <- df_MSRA %>%
    filter((`20-39` < `40-59` & `40-59` < `60-79`) | 
             (`20-39` == `40-59` & `40-59` < `60-79`) |
             (`20-39` < `40-59` & `40-59` == `60-79`))
  
  c(df_MSRD_increasing %>% distinct(ref_junID) %>% pull(),
    df_MSRA_increasing %>% distinct(ref_junID) %>% pull()) %>% unique() %>% length()
  
  # df_MSRD_increasing <- df_MSRD %>%
  #   filter((`20-39` < `40-59` & `40-59` < `60-79`))
  # 
  # df_MSRA_increasing <- df_MSRA %>%
  #   filter((`20-39` < `40-59` & `40-59` < `60-79`))
  
  
  ## Genes with introns increasing noise with age
  genes_increasing <- c(df_MSRD_increasing %>% distinct(gene_name) %>% pull(),
                        df_MSRA_increasing %>% distinct(gene_name) %>% pull()) %>% unique()
  genes_increasing %>% length()

  ## GENES BACKGROUND
  bg_genes <- df_age_groups_tidy %>% unnest(gene_name) %>% distinct(gene_name) %>% pull()
  bg_genes %>% length()
  
  
  #############################################
  ## MSR OF GENES INCREASING NOISE AT DONOR
  ## AND ACCEPTOR
  #############################################
  

  ego_MSR <- clusterProfiler::enrichGO(
    gene          = genes_increasing,
    universe      = bg_genes,
    keyType       = "SYMBOL",
    OrgDb         = "org.Hs.eg.db", ##Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
    ont           = "ALL",
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
    )
  
  
  library('org.Hs.eg.db')
  mapIds(org.Hs.eg.db, genes_increasing, 'ENTREZID', 'SYMBOL')
 
  ekegg_MSR <- clusterProfiler::enrichKEGG(
    gene          = mapIds(x = org.Hs.eg.db, keys = genes_increasing, column = 'ENTREZID', keytype = 'SYMBOL'),
    organism      = "hsa",
    keyType       = "kegg",
    universe      = mapIds(x = org.Hs.eg.db, keys = bg_genes, column = 'ENTREZID', keytype = 'SYMBOL'), 
    pAdjustMethod = "fdr")
  
  xlsx::write.xlsx2(x = genes_increasing, 
                    file = paste0(getwd(), "/results/_paper/results/genes.xlsx"), 
                    sheetName = "genes_increasing", row.names = F, append = T)
  
  xlsx::write.xlsx2(x = bg_genes, 
                    file = paste0(getwd(), "/results/_paper/results/genes.xlsx"), 
                    sheetName = "bg_genes", row.names = F, append = T)
  
  #############################################
  ## PLOTS
  #############################################
  
  ## DOTPLOT
  
  #edox2 <- enrichplot::pairwise_termsim(ego_MSR, showCategory = 30)
  #clusterProfiler::heatplot(ego_MSR,showCategory = 30)
  
  
  plotGO <- clusterProfiler::dotplot(ego_MSR %>%
                                       filter(Description != "RNA splicing, via transesterification reactions with bulged adenosine as nucleophile",
                                              Description != "plasma membrane bounded cell projection cytoplasm",
                                              ONTOLOGY == "CC"), 
                                     x = "GeneRatio", 
                                     showCategory = 30, 
                                     split="ONTOLOGY") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
    xlab("Gene Ratio") +
    ggforce::facet_row(ONTOLOGY~., scales = "free_x", space = "free") +
    custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top",
          legend.box="horizontal",
          plot.margin = margin(0,0,0,0),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(b = -9)) + 
    scale_size(range = c(1, 5))+
    coord_flip() +
    guides(size = guide_legend(title = "Gene Count: "),
           colour = guide_legend(title = "q: "))               
  
  plotGO
  
  plotKEGG <-  clusterProfiler::dotplot(ekegg_MSR %>%
                                          mutate(ONTOLOGY = "KEGG") %>%
                                          filter(Description %in% c("Lysosome", 
                                                                    "Thyroid hormone signaling pathway", 
                                                                    "Huntington disease"     ,                          
                                                                    "Synaptic vesicle cycle",                         
                                                                    "Spliceosome",                 
                                                                    "Pathways of neurodegeneration - multiple diseases", 
                                                                    "Parkinson disease",
                                                                    "Glutamatergic synapse", 
                                                                    "Ubiquitin mediated proteolysis",
                                                                    "Amyotrophic lateral sclerosis",                    
                                                                    "Neurotrophin signaling pathway",
                                                                    "Dopaminergic synapse",                            
                                                                    "Nucleocytoplasmic transport",
                                                                    "RNA degradation")), 
                                      showCategory = 20, 
                                      split="ONTOLOGY") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    xlab("Gene Ratio") +
    ggforce::facet_row(ONTOLOGY~., scales = "free", space = "free") +
    coord_flip() +
    custom_ggtheme +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
      legend.position = "top",
      legend.box="horizontal",
      #legend.margin=margin(0,0,0,0),
      plot.margin = margin(0,0,0,0),
      legend.box.margin=margin(b = -9)) + 
    scale_size(range = c(1, 5))+
    guides(colour = guide_legend(title = "q: "))
  
  plotKEGG
  
  plot_legend <- ggpubr::get_legend(plotKEGG)
  
  ggpubr::ggarrange(plotGO + ggpubr::rremove("legend"),
                    plotKEGG + ggpubr::rremove("ylab")+ ggpubr::rremove("legend"), 
                    common.legend = T,
                    ncol = 2,
                    nrow = 1,
                    widths = c(1.5,1),
                    legend.grob = plot_legend)
  
  ggplot2::ggsave(filename = paste0(getwd(), 
                                    "/results/_paper/figures/panel7b.svg"), 
                  width = 180, height = 80, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(getwd(), 
                                    "/results/_paper/figures/panel7b.png"), 
                  width = 180, height = 80, units = "mm", dpi = 300)
  
  # ############################################
  # ## BAR PLOT -------------------------------
  # ############################################
  # 
  # plot1 <- barplot(ego_MSR, 
  #                  x = "Count", 
  #                  order=T, 
  #                  showCategory = 40) +
  #   scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  #   xlab("Gene Count") +
  #   ggforce::facet_row(facets = vars(ONTOLOGY),
  #                      scales = "free_x", space = "free",
  #                      strip.position = "top") +
  #   coord_flip() +
  #   theme_light() + 
  #   theme(text = element_text(colour = "black",size = 12),
  #         axis.line = element_line(colour = "black"),
  #         
  #         strip.text = element_text(colour = "black", size = "12"),
  #         legend.position = "top",
  #         axis.ticks = element_line(colour = "black", linewidth = 2),
  #         axis.text.x = element_text(colour = "black",
  #                                    angle = 90,
  #                                    vjust = 0.3,
  #                                    hjust = 1),
  #         axis.text.y = element_text(colour = "black",
  #                                    vjust = 0.3,
  #                                    hjust = 1)) +
  #   guides(fill = guide_legend(title = "q",
  #                              label.position = "bottom") ) #+xlim(c(0,0.072))
  # plot1
  # 
  # 
  # plot2 <- barplot(ekegg_MSR %>%
  #                    mutate(ONTOLOGY = "KEGG"), 
  #                  x = "Count", 
  #                  order=T, 
  #                  color = "qvalue") +
  #   scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  #   xlab("Gene Ratio") +
  #   ggforce::facet_row(facets = vars(ONTOLOGY),
  #                      scales = "free_x", space = "free",
  #                      strip.position = "top") +
  #   theme_light() + 
  #   theme(text = element_text(colour = "black",size = 12),
  #         axis.line = element_line(colour = "black"),
  #         
  #         strip.text = element_text(colour = "black", size = "12"),
  #         legend.position = "top",
  #         #panel.grid.major.x = element_blank(),
  #         #panel.grid.major.y = element_blank(),
  #         axis.ticks = element_line(colour = "black", size = 2),
  #         axis.text.x = element_text(colour = "black",
  #                                    angle = 90,
  #                                    vjust = 0.3,
  #                                    hjust = 1),
  #         axis.text.y = element_text(colour = "black",
  #                                    vjust = 0.3,
  #                                    hjust = 1)) +
  #   guides(fill = guide_legend(title = "q",
  #                              label.position = "bottom") ) +
  #   coord_flip() #+xlim(c(0,0.072))
  # 
  # 
  # ggpubr::ggarrange(plot0,
  #                   plot2+ ggpubr::rremove("ylab"), 
  #                   common.legend = T,
  #                   widths = c(3,1))
  # 
  # ggplot2::ggsave(filename = paste0(getwd(), 
  #                                   "/results/_paper/figures/go_barplot_ALL_age_brain_acceptor.png"), 
  #                 width = 260, height = 160, units = "mm", dpi = 300)
  
  
  #############################################
  ## SAVE TERMS AS SUPPLEMENTARY DATA 
  #############################################
  
  
  write.csv(x = ego_MSR %>% as_tibble(),
            file = paste0(getwd(), "/results/_paper/results/go_msr_brain.csv"),
            row.names = F)
  write.csv(x = ekegg_MSR %>% as_tibble(),
            file = paste0(getwd(), "/results/_paper/results/kegg_msr_brain.csv"),
            row.names = F)
  
  
  
  #############################################
  ## EVALUATE TERMS
  #############################################
  
  
  
  query <- paste0("SELECT DISTINCT ref_junID, ref_coordinates FROM 'intron' WHERE ref_junID = '88181'")
  TARDBP_introns <- dbGetQuery(con, query) %>% as_tibble() 
  
  df_MSRD_increasing %>%
    filter(gene_name == "TARDBP") %>%
    inner_join(y = TARDBP_introns,
               by = "ref_junID") %>%
    dplyr::select(-ref_junID) %>%
    dplyr::rename(intron_coordinates = ref_coordinates) %>%
    dplyr::relocate(intron_coordinates) 
  
  
  ## All enrichment
  ego_MSR %>% 
    as_tibble() %>%
    dplyr::select(Description, p.adjust) %>%
    arrange(p.adjust) %>%
    print(n=140)
  
  
  ## Find terms with TARDBP (the gene that encodes for the TDP43 gene)
  ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = geneID,
                      pattern = "TARDBP")) %>%
    dplyr::select( ONTOLOGY, ID,Description,GeneRatio,p.adjust,geneID)
  
  
  ## Get genes with "splic" term
  ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = Description,
                      pattern = "splic")) %>%
    arrange(p.adjust)
  
  
  ## Get genes with "tau" term
  ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = Description,
                      pattern = "tau")) %>%
    arrange(p.adjust)
  
  ## Get genes with "neur" term
  ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = Description,
                      pattern = "neur")) %>%
    arrange(p.adjust)
  
  ego_MSR_neur_genes <- ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = Description,
                      pattern = "neur")) %>%
    mutate(geneID = str_replace_all(string = geneID, pattern = "/", replacement = " "))
  
  ego_MSR_neur_genes_tidy <- paste(ego_MSR_neur_genes$geneID,collapse = ",") %>%
    str_replace_all(pattern = "/", replacement = " ") %>%
    str_split_1(pattern = " ") %>%
    unique() %>%
    sort()
  
  any(ego_MSR_neur_genes_tidy == "MAPT")
  any(ego_MSR_neur_genes_tidy == "SNCA")
  any(ego_MSR_neur_genes_tidy == "PSEN1")
  any(ego_MSR_neur_genes_tidy == "PSEN2")
  
  
  ## Get genes with "dendritic" term
  ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = Description,
                      pattern = "dend")) %>%
    arrange(p.adjust)
  
  ## Get genes with "proteosome" term
  ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = Description,
                      pattern = "prote")) %>%
    arrange(p.adjust)
  
  
  ## Get genes with "splic" term
  ego_MSR_splic_genes <- ego_MSR %>% 
    as_tibble() %>%
    filter(str_detect(string = Description,
                      pattern = "splic")) %>%
    mutate(geneID = str_replace_all(string = geneID, pattern = "/", replacement = " "))
  
  ego_MSR_splic_genes_tidy <- paste(ego_MSR_splic_genes$geneID,collapse = ",") %>%
    str_replace_all(pattern = "/", replacement = " ") %>%
    str_split_1(pattern = " ") %>%
    unique() %>%
    sort()
  
  ego_MSR_splic_genes_tidy <- data.frame(gene_name = ego_MSR_splic_genes_tidy)
  
  RBPs <- xlsx::read.xlsx(file = paste0(getwd(), "/data/RBPs_subgroups.xlsx"),sheetIndex = 1,header = T)
  
  ego_MSR_splic_genes_tidy_RBPs <-  ego_MSR_splic_genes_tidy %>%
    left_join(y = RBPs,
              by = c("gene_name" = "name"))
  
  
  (ego_MSR_splic_genes_tidy_RBPs %>%
      drop_na() %>%
      distinct(gene_name) %>%
      nrow() * 100) / (RBPs %>% 
                         drop_na() %>%
                         distinct(name) %>%
                         nrow())
  
  
  ## splicing regulator
  (ego_MSR_splic_genes_tidy_RBPs %>% dplyr::count(Splicing.regulation) %>% filter(Splicing.regulation == 1) %>% pull(n) * 100) /
    (RBPs %>% dplyr::count(Splicing.regulation) %>% filter(Splicing.regulation == 1) %>% pull(n))
  
  ego_MSR_splic_genes_tidy_RBPs %>% filter(Splicing.regulation == 1) %>% pull(gene_name)
  
  ## spliceosome
  (ego_MSR_splic_genes_tidy_RBPs %>% dplyr::count(Spliceosome) %>% filter(Spliceosome == 1) %>% pull(n) * 100) /
    (RBPs %>% dplyr::count(Spliceosome) %>% filter(Spliceosome == 1) %>% pull(n))
  
  ego_MSR_splic_genes_tidy_RBPs %>% filter(Spliceosome == 1) %>% pull(gene_name)
  RBPs %>% filter(Spliceosome == 1) %>% nrow()
  
  ## exon-junction complex
  (ego_MSR_splic_genes_tidy_RBPs %>% dplyr::count(Exon.Junction.Complex) %>% filter(Exon.Junction.Complex == 1) %>% pull(n) * 100) /
    (RBPs %>% dplyr::count(Exon.Junction.Complex) %>% filter(Exon.Junction.Complex == 1) %>% pull(n))
  
  ego_MSR_splic_genes_tidy_RBPs %>% filter(Exon.Junction.Complex == 1) %>% pull(gene_name)
  RBPs %>% filter(Exon.Junction.Complex == 1) %>% nrow()
  
  
  
}






####################################
# SPLICEOSOMAL RBPs AND THEIR TPMs
####################################

age_stratification_RBP_uncorrected_TPM_lm <- function() {
  
  
  # Only measure RBPs that are splicing regulator, spliceosome, and Exon.Junction.Complex
  all_RBPs <- xlsx::read.xlsx(file = '/home/sruiz/PROJECTS/splicing-project-recount3/data/RBPs_subgroups.xlsx', 
                              header = TRUE,
                              sheetIndex = 1) %>% 
    as_tibble() %>%
    drop_na()
  

  
  
  # 1. The age stratification database has been corrected by RIN number
  
  # 2. Per RBP, get their TPM values across the age groups
  
  project_id <- "BRAIN"
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) %>%
    filter(SRA_project == project_id)
  
    
  ## Get the MSR values of the RBPs across the age groups
  brain_rbp_MSR_age <- map_df(df_metadata$cluster %>% unique(), function(age_group) {
    
    # age_group <- (df_metadata$cluster %>% unique())[1]
      
    print(paste0(age_group))
    
    ###########################################
    ## GET DATA FROM THE DATABASE
        
    query <- paste0("SELECT DISTINCT ref_junID, ref_type, gene_tpm, transcript_id
                    FROM '", age_group, "_", project_id, "_misspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
        
        
    
    ## Add never-misspliced info
    query <- paste0("SELECT ref_junID, ref_type, gene_tpm, transcript_id 
                        FROM '", 
                    age_group, "_", project_id, "_nevermisspliced'")
    introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
    
    
    
    ## Add transcript info
    query <- paste0("SELECT id, gene_id FROM 'transcript' WHERE id IN (", 
                    paste(introns$transcript_id %>% unique(), collapse =","), ")")
    introns %>% nrow()
    introns <- introns %>%
      inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                 by = c("transcript_id" = "id"))
    
    
    ## Add gene info
    query <- paste0("SELECT id, gene_name FROM 'gene' WHERE id IN (", 
                    paste(introns$gene_id %>% unique(), collapse =","), ")")
    introns <- introns %>%
      inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                 by = c("gene_id" = "id"))
    
    
    introns %>%
      filter(gene_name %in% all_RBPs$name) %>%
      mutate(age_cluster = age_group) %>%
      dplyr::select(-c("transcript_id", "gene_id")) %>%
      return()
  })
  
  brain_rbp_MSR_age <- brain_rbp_MSR_age %>% as_tibble()
  
  # 3. Get theannotated introns present in the 3 age groups
  common_introns <- brain_rbp_MSR_age %>%
    dplyr::count(ref_junID) %>%
    filter(n == 3) %>%
    pull(ref_junID)
  
  brain_rbp_MSR_age_common <- brain_rbp_MSR_age %>%
    filter(ref_junID %in% common_introns) 
  
  # 4. Compare the TPM values
  ggplot(data = brain_rbp_MSR_age_common %>%
           mutate(gene_tpm = log10(gene_tpm))) +
    geom_tile(aes(x = age_cluster, y = gene_name, fill = gene_tpm))
  
  
}

age_stratification_RBP_uncorrected_TPM_lm <- function(project_id = c("BRAIN")) {
  
  
  ################################
  ## LOAD RBPs OF INTEREST
  ################################
  
  if ( !exists("all_RBPs") ) {
    
    # Only measure RBPs that are splicing regulator, spliceosome, and Exon.Junction.Complex
    all_RBPs <- xlsx::read.xlsx(file = '/home/sruiz/PROJECTS/splicing-project-recount3/data/41586_2020_2077_MOESM3_ESM.xlsx', 
                                sep = '\t', header = TRUE,
                                sheetIndex = 1) %>% 
      as_tibble() 
  }
  
  
  all_RBPs_tidy <- all_RBPs %>%
    filter(Splicing.regulation == 1 |
             Spliceosome == 1 | 
             Exon.Junction.Complex == 1) %>%
    dplyr::select(name, id,
                  Splicing.regulation,
                  Spliceosome,
                  Exon.Junction.Complex)
  
  write.table(x = all_RBPs_tidy$name,
              file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs.csv"),
              row.names = F, col.names = F, quote = F )
  saveRDS(object = all_RBPs_tidy,
          file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/data/all_RBPs_tidy.rds"))
  

  ################################
  ## CONNECT TO THE DATABASE
  ################################
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) %>%
    filter(SRA_project == project_id)
  
  
  ## Get the TPM of the RBPs selected
  
  for (project_id in project_id) {
    ## project_id <- "BRAIN"
    
    ## Load and tidy the uncorrected TPMs
    tpm_uncorrected <- map_df((df_metadata$region %>% unique()), function(cluster) {
      
      print(paste0(Sys.time(), " - ", cluster, " ... "))
      if (file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/EXPRESSION_ANALYSIS/RBP/",
                             cluster, "/tpm_ensembl105.csv"))) {
        read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/EXPRESSION_ANALYSIS/RBP/",
                               cluster, "/tpm_ensembl105.csv"), header = T, fileEncoding = "UTF-8") 
      } else {
        return(NULL)
      }   
    })
    
    tpm_uncorrected_tidy <- tpm_uncorrected %>%
      dplyr::rename(gene = "X") %>%
      as_tibble() %>%
      distinct(gene, .keep_all = T)
    
    
    ## Load and tidy the sample metadata
    sample_metadata <- map_df((df_metadata$region %>% unique()), function(cluster) {
      
      print(paste0(Sys.time(), " - ", cluster, " ... "))
      if (file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/EXPRESSION_ANALYSIS/RBP/",
                             cluster, "/covariates.csv"))) {
        read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/EXPRESSION_ANALYSIS/RBP/",
                               cluster, "/covariates.csv"), header = T, fileEncoding = "UTF-8") 
      } else {
        return(NULL)
      }   
    })
    
    
    
    print(sample_metadata %>% nrow)
    ## Check the samples are the same
    if ( !( identical(x = names(tpm_uncorrected_tidy)[-1],
                      y = names(sample_metadata)[-1]) ) ) {
      print("ERROR!")
      break;
    }
    
    # ## Select only the samples from the current database
    # tpm_uncorrected_tidy <- tpm_uncorrected_tidy %>%
    #   dplyr::select(gene, all_of(str_replace_all(string = df_metadata$individual,
    #                                              pattern = "-",
    #                                              replacement = "."))) 
    # sample_metadata <- sample_metadata %>%
    #   dplyr::select(covariates, all_of(str_replace_all(string = df_metadata$individual,
    #                                                    pattern = "-",
    #                                                    replacement = "."))) 
    # 
    # ## Check the samples are the same after subsampliing by the samples of the current database
    # if ( !( identical(x = names(tpm_uncorrected_tidy)[-1],
    #                   y = names(sample_metadata)[-1]) ) ) {
    #   print("ERROR!")
    #   break;
    # }
    
    # RBPs <- (tpm_uncorrected_tidy$gene)
    
    
    
    
    # tpm_uncorrected[rowSums(tpm_uncorrected[, c(2:ncol(tpm_uncorrected))]),]
    
    ## Obtain values and run lm
    df_lm_output <- map_df(all_RBPs_tidy$id, function(RBP) {
      
      
      # RBP <- all_RBPs_tidy$id[1]
      # RBP <- "ENSG00000169045"
      print(RBP)
      
      ## Filter by the current RBP
      tpm <- tpm_uncorrected_tidy %>%
        filter(gene == RBP) %>%
        mutate(gene = "tpm") %>%
        select_if(~ !any(is.na(.))) %>%
        as_tibble()
      
      RBP_sample_metadata <- sample_metadata %>%
        dplyr::select(covariates, all_of(names(tpm %>% dplyr::select(-gene))))%>%
        drop_na() 
      
      if ( nrow(tpm) > 0 && rowSums(tpm[, c(2:ncol(tpm))], na.rm = T) != 0) {
        
        
        df <- rbind(tpm %>%  dplyr::rename(covariates = "gene"),
                    RBP_sample_metadata) %>%
          gather(sample, value, -covariates)  %>%
          drop_na() %>%
          spread(key = covariates, value = value) %>%
          dplyr::mutate(gtex.age = ifelse(gtex.age %>% as.double() == 1 | gtex.age %>% as.double() == 2, 30, gtex.age %>% as.double()),
                        gtex.age = ifelse(gtex.age %>% as.double() == 3 | gtex.age %>% as.double() == 4, 50, gtex.age %>% as.double()),
                        gtex.age = ifelse(gtex.age %>% as.double() == 5 | gtex.age %>% as.double() == 6, 70, gtex.age %>% as.double()),
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
                          #gtex.smnabtcht +
                          gtex.dthhrdy +
                          gtex.sex +
                          #gtex.smtsd +
                          gtex.age,
                        data = df)
        
        
        
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
    

    ## Filter by age covariate and order
    df_lm_output_age <- df_lm_output %>%
      filter(str_detect(string = covariate, pattern = "gtex.age")) %>%
      arrange(Estimate)
    
    if ( exists("RBPs_annotated") ) {
      
      RBPs_annotated_tidy <- left_join(x = RBPs_annotated %>% as_tibble(), 
                                       y = all_RBPs_tidy %>% as_tibble(), 
                                       by = c("ensembl_gene_id" = "ensembl_gene_id"))
    } else {
      RBPs_annotated_tidy <- all_RBPs_tidy
    }
    
    ## Add gene SYMBOL info
    df_lm_age_tidy <- left_join(x = df_lm_output_age, 
                                y = RBPs_annotated_tidy %>% drop_na(), 
                                by = c("RBP_ID" = "id")) %>%
      group_by(RBP_ID) #%>%
    #distinct(pval, .keep_all = T)
    
    
    write_csv(x = df_lm_age_tidy,
              file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/EXPRESSION_ANALYSIS/RBP/",
                            project_id,"_tpm_lm_all_ensembl105.csv"))
    
  }
  
  ###############################################
  ## Get RBPs that decrease expression with age
  ###############################################
  
  
  # RBP_decrease <- df_lm_age_tidy %>%
  #   filter(Estimate < 0) %>%
  #   arrange(desc(abs(Estimate)))
  # 
  # RBP_increase <- df_lm_age_tidy %>%
  #   filter(Estimate > 0) %>%
  #   arrange(desc(Estimate))
  # 
  # intersect(RBP_decrease$hgnc_symbol,
  #           RBP_increase$hgnc_symbol)
  # 
  # RBP_decrease$hgnc_symbol %>% unique() %>% sort()
  
  
}

age_stratification_RBPs_affected_age <- function (project_id = "BRAIN") {
  
  ####################################
  # LOAD RBPs ORIGINAL PAPER AND TIDY
  ####################################
  
  all_RBPs_tidy <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/data/RBPs_subgroups.xlsx"))  %>%
    dplyr::select(RBP_name = name, ensembl_gene_id = id) %>%
    distinct(RBP_name, .keep_all = T) 
  
  ####################################
  # LOAD RBPs LM RESULTS
  ####################################
  
  RBPs_age_lm <- read.csv(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/EXPRESSION_ANALYSIS/RBP/", project_id,
                                        "_tpm_lm_all_ensembl105.csv")) %>%
    drop_na() %>% 
    as_tibble() 
  
  
  ##########################################
  # CLASSIFY AFFECTED / NOT AFFECTED BY AGE
  ##########################################
  
  RBP_age <- RBPs_age_lm %>%
    filter(pval <= 0.05) %>%
    drop_na() 
  
  write.table(x = RBP_age$name  %>% sort(),
              file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_decrease.csv"),
              row.names = F, col.names = F, quote = F )
  
  RBP_notage <- RBPs_age_lm %>%
    filter( !(name  %in% RBP_age$name) ) %>%
    drop_na() 
  
  write.table(x = RBP_notage$name ,
              file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_other.csv"),
              row.names = F, col.names = F,quote = F)
  
  
  plot(density(RBP_notage$Estimate))
  lines(density(RBP_age$Estimate), col = "red")
  
  RBPs_age_tidy <- RBPs_age_lm %>%
    drop_na() %>%
    mutate(age_type = ifelse(pval <= 0.05,
                             "affected_age",
                             "not_affected_age"))
  RBPs_age_tidy %>%
    dplyr::count(age_type)
  
  write.csv(x = RBPs_age_tidy,
            file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/EXPRESSION_ANALYSIS/RBP/", project_id,
                          "_tpm_lm_all_ensembl105_tidy.csv"))
}


####################################
# ENCORI
####################################

overlap_ENCORI_introns_MSR <- function(project_id = "BRAIN") {
  
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  chain <- rtracklayer::import.chain(con = "data/hg19ToHg38.over.chain")
  
  
  ###########################
  ## TIDY INTRONS
  ###########################
  
  query <- paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query)

 
  df_positions <- map_df (df_intron$ref_coordinates, function (coordinate) {
    
    # coordinate <- "chr7:69596434-69597233:-"
    
    print(coordinate)
    
    chr_j <- str_sub(string = coordinate,
                     start = 1,
                     end = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]-1)
    start_j <- str_sub(string = coordinate,
                       start = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]+1,
                       end = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]-1) %>% as.integer()
    end_j <- str_sub(string = coordinate,
                     start = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]+1,
                     end = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]-1) %>% as.integer()
    strand_j <- str_sub(string = coordinate,
                        start = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]+1,
                        end = coordinate %>% stringr::str_count())
    print(coordinate)
    
    return(data.frame(ref_coordinates = coordinate,
                      seqnames = chr_j %>% as.character(),
                      start = start_j %>% as.integer(),
                      end = end_j %>% as.integer(),
                      strand = strand_j %>% as.character()))
    
  })
  
  df_intron_tidy <- df_intron %>%
    left_join(y = df_positions,
              by = "ref_coordinates")
  
  ###########################
  ## MSR_D
  ###########################

  
  df_MSRD <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v",
                                   gtf_version,"/",main_project,"/results/df_MSRD_common_introns_all_projects.rds")) %>%
    dplyr::rename("body_site" = project_id) %>%
    filter(body_site == project_id) %>%
    left_join( y = df_intron_tidy %>% dplyr::select(ref_junID, start, end, seqnames, strand),
               by.x = "ref_junID",
               by.y = "ref_junID" )
  
  df_MSRD_introns_increasing <- df_MSRD %>%
    filter((`20-39` < `40-59` |
             `20-39` < `60-79`) | 
             `40-59` < `60-79`) %>%
    mutate(start = start - 25) %>% 
    GRanges()
  
  
  df_MSRD_introns_decreasing <- df_MSRD %>%
    filter( !(ref_junID %in% df_MSRD_introns_increasing$ref_junID) ) 
  df_MSRD_introns_decreasing <- df_MSRD_introns_decreasing[sample(nrow(df_MSRD_introns_decreasing), 
                                                                  (df_MSRD_introns_increasing$ref_junID %>% length()) ), ]
  df_MSRD_introns_decreasing <- df_MSRD_introns_decreasing %>%
    mutate(start = start - 25) %>% 
    GRanges()
  
  # genes_MSRD_increasing <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, 
  #                                                "/", main_project, "/results/genes_increase_MSRD_", project_id, ".rds"))
  # genes_MSRD_increasing <- merge(genes_MSRD_increasing,
  #                                y = df_intron_tidy %>% dplyr::select(ref_junID, start, end, seqnames, strand),
  #                                by.x = "ref_junID",
  #                                by.y = "ref_junID")
  # genes_MSRD_increasing <- genes_MSRD_increasing %>%
  #   mutate(start = start - 25, end = end + 25) %>% 
  #   GRanges()
  # 
  # 
  # genes_MSRD_decreasing <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, 
  #                                                "/", main_project, "/results/genes_decrease_MSRD_", project_id, ".rds"))
  # genes_MSRD_decreasing <- merge(genes_MSRD_decreasing,
  #                                y = df_intron_tidy %>% dplyr::select(ref_junID, start, end, seqnames, strand),
  #                                by.x = "ref_junID",
  #                                by.y = "ref_junID")
  # genes_MSRD_decreasing <- genes_MSRD_decreasing %>%
  #   mutate(start = start - 25, end = end + 25) %>% 
  #   GRanges()
  
  
  
  ###########################
  ## MSR_A
  ###########################
  
  df_MSRA <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v",
                                   gtf_version,"/",main_project,"/results/df_MSRA_common_introns_all_projects.rds")) %>%
    dplyr::rename("body_site" = project_id) %>%
    filter(body_site == project_id) %>%
    left_join( y = df_intron_tidy %>% dplyr::select(ref_junID, start, end, seqnames, strand),
               by.x = "ref_junID",
               by.y = "ref_junID" )
  
  df_MSRA_introns_increasing <- df_MSRA %>%
    filter( (`20-39` < `40-59` &
             `20-39` < `60-79`) | 
             `40-59` < `60-79`) %>%
    mutate(end = end + 25) %>% 
    GRanges()
  
  df_MSRA_introns_decreasing <- df_MSRA %>%
    filter( !(ref_junID %in% df_MSRA_introns_increasing$ref_junID) ) 
  df_MSRA_introns_decreasing <- df_MSRA_introns_decreasing[sample(nrow(df_MSRA_introns_decreasing), 
                                                                  (df_MSRA_introns_increasing$ref_junID %>% length()) ), ]
  df_MSRA_introns_decreasing <- df_MSRA_introns_decreasing %>%
    mutate(end = end + 25) %>% 
    GRanges()
  # df_MSRA_introns_decreasing <- df_MSRA %>%
  #   filter((`20-39` > `40-59` #&
  #          #  `40-59` > `60-79`
  #          )) %>%
  #   mutate(end = end + 25) %>% 
  #   GRanges()
  
  # genes_MSRA_increasing <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, 
  #                                                "/", main_project, "/results/genes_increase_MSRA_", project_id, ".rds"))
  # genes_MSRA_increasing <- merge(genes_MSRA_increasing,
  #                                y = df_intron_tidy %>% dplyr::select(ref_junID, start, end, seqnames, strand),
  #                                by.x = "ref_junID",
  #                                by.y = "ref_junID")
  # 
  # genes_MSRA_increasing <- genes_MSRA_increasing %>%
  #   mutate(start = start - 25, end = end + 25) %>% 
  #   GRanges()
  # 
  # genes_MSRA_decreasing <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, 
  #                                                "/", main_project, "/results/genes_decrease_MSRA_", project_id, ".rds"))
  # genes_MSRA_decreasing <- merge(genes_MSRA_decreasing,
  #                                y = df_intron_tidy %>% dplyr::select(ref_junID, start, end, seqnames, strand),
  #                                by.x = "ref_junID",
  #                                by.y = "ref_junID")
  # 
  # genes_MSRA_decreasing <- genes_MSRA_decreasing %>%
  #   mutate(start = start - 25, end = end + 25) %>% 
  #   GRanges()
  
  ######################################################################
  ## GET OVERLAPS - ENCORI AND INTRONS WITH INCREASING MSR_D AND MSR_A
  ######################################################################
  
  
 
  
  
  df_RBP_ENCORI_MSR_result <- map_df( c("RBPs_affected_age", "RBPs_notaffected_age"), function(type) {
    
    # type <- "RBPs_notaffected_age"
    # type <- "RBPs_affected_age"
    
    if ( !exists("ensembl105") ) {
      ensembl105 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")  
    } 
    
    ## LOAD RPB list
    if ( type == "RBPs_affected_age" ) {
      
      RBPs <- read.csv(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_decrease.csv",
                       header = F)
      
      
      # RBPs <- spread_matrix %>%
      #   filter(`20-39` > `40-59` ,
      #          `40-59` > `60-79`) %>%
      #   left_join(y = ensembl105 %>%
      #               as_tibble() %>%
      #               dplyr::select(gene_id, gene_name) %>%
      #               distinct(gene_id, .keep_all = T),
      #             by = c("RBP" = "gene_id")) %>%
      #   dplyr::select(V1 = gene_name) 

      
    } else {
      
      RBPs <- read.csv(file = "/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/RBPs_other.csv", 
                       header = F)
      
      # RBPs <-  spread_matrix %>%
      #   filter(!(RBP %in%  (spread_matrix %>%
      #                         filter(`20-39` > `40-59` ,
      #                                `40-59` > `60-79`) %>%
      #                         pull(RBP)))) %>%
      #   left_join(y = ensembl105 %>%
      #               as_tibble() %>%
      #               dplyr::select(gene_id, gene_name) %>%
      #               distinct(gene_id, .keep_all = T),
      #             by = c("RBP" = "gene_id")) %>%
      #   dplyr::select(V1 = gene_name) 
      
    }
    
    ##############################
    ## Check ENCORI iCLIP results
    ##############################
    
    map_df( RBPs$V1 %>% sort(), function (RBP) {
      
      # RBP <- "ADAR"
      file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/paper_figures/iCLIP/results/ENCORI/",
                          "/ENCORI_hg19_", RBP, "_allgenes.txt")
      
      #print(file_name)
      
      if ( file.exists(file_name) ) {
        
        ENCORI_RBP_result <- read.delim(file = file_name, header = T, skip = 3, sep = "\t") 
        
        if ( ENCORI_RBP_result %>% nrow() > 1) {
          
          print(paste0(RBP, " - ", type))
          
          ENCORI_RBP_result <- ENCORI_RBP_result %>% 
            as.data.frame() %>%
            mutate(chromosome = chromosome %>% as.factor(),
                   strand = strand %>% as.factor()) %>%
            distinct(geneID, .keep_all = T) %>%
            dplyr::select(RBP, seqnames = chromosome, start = broadStart, end = broadEnd, strand) %>% 
            GenomicRanges::GRanges()
          
          
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
          
          overlaps_MSRD <- GenomicRanges::findOverlaps(query = df_MSRD_introns_increasing,
                                                       subject = ENCORI_RBP_result_GRh38,
                                                       type = "any",
                                                       ignore.strand = F)
          
          overlaps_MSRDd <- GenomicRanges::findOverlaps(query = df_MSRD_introns_decreasing,
                                                        subject = ENCORI_RBP_result_GRh38,
                                                        type = "any",
                                                        ignore.strand = F)
          
          
          #######################
          ## Overlaps - MSR_A
          #######################
          
          overlaps_MSRA <- GenomicRanges::findOverlaps(query = df_MSRA_introns_increasing,
                                                       subject = ENCORI_RBP_result_GRh38,
                                                       type = "any",
                                                       ignore.strand = F)
          
          overlaps_MSRAd <- GenomicRanges::findOverlaps(query = df_MSRA_introns_decreasing,
                                                        subject = ENCORI_RBP_result_GRh38,
                                                        type = "any",
                                                        ignore.strand = F)
          
          #######################
          ## Return result
          #######################
          
          
          data.frame(RBP_name = RBP,
                     RBP_type = type,
                     ovlps_MSRD = df_MSRD_introns_increasing[queryHits(overlaps_MSRD),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRA = df_MSRA_introns_increasing[queryHits(overlaps_MSRA),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRD = df_MSRD_introns_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRA = df_MSRA_introns_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRD_perc = (((df_MSRD_introns_increasing[queryHits(overlaps_MSRD),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                          df_MSRD_introns_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow()),
                     ovlps_MSRA_perc = (((df_MSRA_introns_increasing[queryHits(overlaps_MSRA),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                          df_MSRA_introns_increasing %>% as_tibble %>% distinct(ref_junID) %>% nrow()),
                     
                     ovlps_MSRDd = df_MSRD_introns_decreasing[queryHits(overlaps_MSRDd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRAd = df_MSRA_introns_decreasing[queryHits(overlaps_MSRAd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRDd = df_MSRD_introns_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     total_MSRAd = df_MSRA_introns_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow,
                     ovlps_MSRD_percd = (((df_MSRD_introns_decreasing[queryHits(overlaps_MSRDd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                           df_MSRD_introns_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow()),
                     ovlps_MSRA_percd = (((df_MSRA_introns_decreasing[queryHits(overlaps_MSRAd),] %>% as_tibble %>% distinct(ref_junID) %>% nrow) * 100) / 
                                           df_MSRA_introns_decreasing %>% as_tibble %>% distinct(ref_junID) %>% nrow())
          ) %>% return()
        }
      }
      
    })
    
  })
  
  write_csv(x = df_RBP_ENCORI_MSR_result,
            file = "paper_figures/iCLIP/results/RBPs_ENCORI_MSR.csv", col_names = T)
  
  df_RBP_ENCORI_MSR_result %>% as_tibble()
  
}

ENCORI_test_results <- function () {
  
  df_RBP_result <- read.csv(file = "paper_figures/iCLIP/results/RBPs_ENCORI_MSR.csv") %>% as_tibble()
  df_RBP_result
  
  plot(density(df_RBP_result %>% 
                 filter(RBP_type == "RBPs_affected_age") %>%
                 pull(ovlps_MSRD_perc)))
  lines( density ( df_RBP_result %>% 
                     filter(RBP_type == "RBPs_notaffected_age") %>%
                     pull(ovlps_MSRD_perc)), col = "red")
  
  df_RBP_result %>% 
    filter(RBP_type == "RBPs_affected_age") %>%
    pull(ovlps_MSRA_perc) %>% summary  
  df_RBP_result %>% 
    filter(RBP_type == "RBPs_notaffected_age") %>%
    pull(ovlps_MSRA_perc) %>% summary
  
  
  ###########################################
  ## RBPs_affected_age vs RBPs_notaffected_age
  ###########################################
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRA_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRA_perc),
              alternative = "greater",
              correct = T)
  
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_affected_age") %>%
                pull(ovlps_MSRA_percd),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRA_percd),
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
                pull(ovlps_MSRD_percd),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRD_percd),
              paired = F,
              alternative = "greater")
 
  
  
  ###########################################
  ## RBPs_affected_age vs RBPs_notaffected_age
  ###########################################
  
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
              alternative = "greater")
  
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRD_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRD_percd),
              paired = T,
              alternative = "greater")
  
  wilcox.test(x = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRA_perc),
              y = df_RBP_result %>% 
                filter(RBP_type == "RBPs_notaffected_age") %>%
                pull(ovlps_MSRA_percd),
              paired = T,
              alternative = "greater")
}

##################################
## BASIC PLOTS
##################################

age_stratification_plot_distances <- function(distance_bp = 30) {
  
  
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  #all_projects <- df_metadata$SRA_project %>% unique()
  all_projects <- "BRAIN"
  
  df_age_distances <- map_df(all_projects, function(project_id) {
    
    # project_id <- (all_projects)[1]
    
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
                   alpha = 0.6,
                   bins = distance_bp * 2,
                   binwidth = 1,
                   position = "identity"
    ) +
    facet_grid(vars(novel_type)) +
    #ggtitle(title) +
    xlab("Distance to the reference intron (in bp)") +
    ylab("Number of unique novel junctions") +
    theme_light() +
    scale_x_continuous(limits = c((distance_bp * -1), distance_bp),
                       breaks = c((distance_bp * -1), (round(distance_bp / 2) * -1), 0, round(distance_bp / 2), distance_bp)) +
    
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
    geom_rect(aes(xmin = 0, xmax = distance_bp, ymin = 1, ymax = 100),
              fill = "grey", color = "black") +
    geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
    geom_rect(aes(xmin = (distance_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
              fill = "grey", alpha = 1, color = "black") +
    geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
    theme_void()
  
  
  plot_distances / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
  
  
  folder_image <- paste0(getwd(), "/recount3/_paper/figures/")
  dir.create(file.path(folder_image), recursive = TRUE, showWarnings = T)
  ggplot2::ggsave(filename = paste0(folder_image, "/", main_project, "_distances_", project_id, ".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  ggplot2::ggsave(filename = paste0(folder_image, "/", main_project, "_distances_", project_id, ".png"), width = 183, height = 183, units = "mm", dpi = 300)
  
}


age_stratification_plot_distances_across_tissues <- function(df = NULL,
                                                             age_projects = c("BRAIN", "BLOOD","MUSCLE"),
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
  
  
  
  df_age_distances <- map_df(age_projects, function(project_id) {
    
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
    group_by(project_id, sample_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    ungroup() %>%
    dplyr::count(ref_junID) %>%
    filter(n == age_groups %>% length()) %>%
    pull(ref_junID)
  
  df_age_distances_common <- df_age_distances %>%
    group_by(project_id, ref_junID, sample_type) %>%
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
  
  
  df_age_distances_common <- df_age_distances_common %>%
    mutate(novel_type = str_replace(string = novel_type,
                                    pattern = "_",
                                    replacement = " "))
  df_age_distances_common <- df_age_distances_common %>%
    mutate(novel_type = factor(novel_type, levels = c("novel donor", "novel acceptor"))) %>%
    mutate(sample_type = factor(sample_type, levels = c("60-79", "40-59", "20-39")))
  
  
  
  
  df_age_distances_common_tidy <- df_age_distances_common %>%
    filter(abs(distance) <= 30) %>%
    group_by(project_id, sample_type, novel_type) %>%
    mutate(median_distance = distance %>% median()) %>%
    ungroup()
  
  df_age_distances_common_tidy_donor <- df_age_distances_common_tidy %>% 
    filter(novel_type == "novel donor") %>%
    dplyr::select(sample_type, project_id, median_distance) %>%
    distinct(sample_type, project_id,median_distance) %>%
    spread(key = sample_type , value = median_distance)
  
  df_age_distances_common_tidy_donor %>%
    column_to_rownames("project_id" ) %>%
    as.matrix()
  
  
  heatmap(x =   df_age_distances_common_tidy_donor %>%
            column_to_rownames("project_id" ) %>%
            as.matrix())
  
  
  
  
  
  
  folder_image <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", main_project, "/figures/")
  dir.create(file.path(folder_image), recursive = TRUE, showWarnings = T)
  ggplot2::ggsave(filename = paste0(folder_image, "/panel4_age_distances_",project_id,".svg"), width = 183, height = 183, units = "mm", dpi = 300)
  
}


age_stratification_plot_distances_proportion <- function(age_levels = c("60-79", "40-59", "20-39"),
                                                         distance_limit = 40,
                                                         QC = F) {
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
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
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
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






age_stratification_plot_MSR <- function(df = NULL,
                                        project_id = "MUSCLE",
                                        age_groups = c("60-79", "40-59", "20-39"),
                                        common = T,
                                        QC = F) {
  
  
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  df_age_groups <- map_df(project_id, function(project_id) {
    
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
  df_age_groups_tidy <- df_age_groups %>%
    inner_join(y = data.table::data.table(ref_junID = common_introns),
               by = "ref_junID") %>%
    distinct(ref_junID, sample_type, .keep_all = T)
  
  
  
  
  
  
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
    mutate(sample_type = factor(sample_type, levels = c( "60-79","20-39", "40-59" ))) %>%
    as_tibble()
  
  
  ggplot(data = df_age_groups_tidy %>% 
           dplyr::select(MSR_D, sample_type) %>%
           gather(key = "MSR_type", value = "MSR", -sample_type) %>%
           mutate(MSR_type = factor(MSR_type, levels = c("MSR_A", "MSR_D")))) + 
    geom_density(aes(x = MSR, fill = sample_type), alpha = 0.8) +
    #ggtitle(title) +
    xlab("MSR") +
    facet_wrap(vars(MSR_type)) +
    ggforce::facet_zoom(xlim = c(0,0.005)) +
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
    guides(fill = guide_legend(title = NULL, ncol = 3,  nrow = 1)) 
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", 
                      main_project, "/figures/panel4_MSR_D_", project_id,".svg")
  ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ggplot(data = df_age_groups_tidy %>% 
           dplyr::select(MSR_A, sample_type) %>%
           gather(key = "MSR_type", value = "MSR", -sample_type) %>%
           mutate(MSR_type = factor(MSR_type, levels = c("MSR_A", "MSR_D")))) + 
    geom_density(aes(x = MSR, fill = sample_type), alpha = 0.8) +
    #ggtitle(title) +
    xlab("MSR_A") +
    facet_wrap(vars(MSR_type)) +
    ggforce::facet_zoom(xlim = c(0,0.005)) +
    theme_light() +
    scale_fill_manual(values =  c("#21908CFF","#FDE725FF","#440154FF"),
                      labels = c("20-39", "40-59", "60-79"),
                      breaks = c("20-39", "40-59", "60-79")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          strip.text =  element_text(colour = "black", size = "12"),
          plot.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL, ncol = 3,  nrow = 1)) 
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", 
                      main_project, "/figures/panel4_MSR_A_",project_id,".svg")
  ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  ########################################
  ## STATISTICAL TEST - MSR_DONOR
  ########################################
  
  
  wilcox.test(x = df_age_groups_tidy %>% filter(sample_type == "20-39") %>% pull(MSR_D),
              y = df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_D),
              alternative = "less",
              correct = T,
              paired = T)
  wilcox.test(x = df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_D),
              y = df_age_groups_tidy %>% filter(sample_type == "60-79") %>% pull(MSR_D),
              alternative = "less",
              paired = T)
  
  t.test(x = df_age_groups_tidy %>% filter(sample_type == "20-39") %>% pull(MSR_D),
         y = df_age_groups_tidy %>% filter(sample_type == "40-59") %>% pull(MSR_D),
         alternative = "less",
         paired = T)
  
  
  ## EFFECT SIZE 
  
  df_MSRD <- df_age_groups_tidy %>%
    distinct(ref_junID, sample_type, .keep_all = T) %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_D) %>%
    mutate(MSR_D = MSR_D %>% round(digits = 4)) %>%
    spread(sample_type, MSR_D)
  
  df_MSRD %>%
    filter(`20-39` < `40-59`,
           `40-59` < `60-79`)
  df_MSRD %>%
    filter(`20-39` > `40-59`,
           `40-59` > `60-79`)
  
  
  rstatix::wilcox_effsize(data = df_MSRD %>%
                            tidyr::gather(key = ref_junID, value = year),
                          formula = ref_junID  ~ year )
  
  ########################################
  ## STATISTICAL TEST - MSR_ACCEPTOR
  ########################################
  
  ## EFFECT SIZE 
  
  df_MSRA <- df_age_groups_tidy %>%
    distinct(ref_junID, sample_type, .keep_all = T) %>%
    dplyr::select(ref_junID,
                  sample_type,
                  MSR_A) %>%
    mutate(MSR_A = MSR_A %>% round(digits = 4)) %>%
    spread(sample_type, MSR_A)
  
  df_MSRA %>%
    filter(`20-39` < `40-59`,
           `40-59` < `60-79`)
  df_MSRA %>%
    filter(`20-39` > `40-59`,
           `40-59` > `60-79`)
  
  df_MSRA_tidy <- df_MSRA %>%
    tidyr::gather(key = ref_junID, value = year) %>%
    dplyr::rename(MSR = year) %>%
    mutate(year = as.numeric(factor(as.matrix(ref_junID))) ) %>%
    as_tibble() %>%
    dplyr::select(MSR, year)
  
  rstatix::wilcox_effsize(data = df_MSRA_tidy,
                          formula = MSR  ~ year )
  
  wilcox.test(x = df_MSRA$`20-39`,
              y = df_MSRA$`40-59`,
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


age_stratification_plot_MSR_across_tissues <- function(age_projects = c("BRAIN", "BLOOD","MUSCLE"),
                                                       age_groups = c("60-79", "40-59", "20-39")) {
  
  
  #######################################
  ## CONNECT TO THE DATABASE
  #######################################
  
  
  con <- dbConnect(RSQLite::SQLite(), database_path)
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  df_age_groups_MSR <- map_df(age_projects, function(project_id) {
    
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
               mutate(sample_type = age_group,
                      project_id = project_id))
      
    })
    
    
    ## GET ONLY COMMON JUNCTIONS ACROSS AGE GROUPS
    
    common_introns <- df_age_groups %>%
      group_by(sample_type) %>%
      distinct(ref_junID, .keep_all = T) %>%
      ungroup() %>%
      dplyr::count(ref_junID) %>%
      filter(n == age_groups %>% length() ) %>% 
      dplyr::select(-n)
    
    
    ## Filter the INTRONS table by the common mis-spliced introns
    df_age_groups_tidy <- df_age_groups %>%
      inner_join(y = common_introns ,
                 by = c("ref_junID" = "ref_junID")) %>%
      distinct(sample_type, ref_junID, .keep_all = T) %>%
      mutate(sample_type = factor(sample_type, levels = c( "20-39", "40-59", "60-79" ))) %>%
      group_by(sample_type) %>%
      mutate(mean_MSRD = MSR_D  %>% mean()) %>%
      mutate(mean_MSRA = MSR_A  %>% mean()) %>%
      mutate(median_MSRD = MSR_D  %>% median()) %>%
      mutate(median_MSRA = MSR_A  %>% median()) %>%
      ungroup()
    
    
    df_age_groups_tidy %>% 
      return()
  })
  
  
  ###############################################
  ## QC
  ###############################################
  
  any(df_age_groups_MSR %>% filter(sample_type == "never") %>% pull(MSR_D) > 0)
  any(df_age_groups_MSR %>% filter(sample_type == "never") %>% pull(MSR_A) > 0)
  
  intersect(df_age_groups_MSR %>% filter(sample_type == "never") %>% pull(ref_junID),
            df_age_groups_MSR %>% filter(sample_type != "never") %>% pull(ref_junID))
  
  
  
  
  
  
  #######################################
  ## PLOT - DONOR
  #######################################
  
  
  
  
  df_age_groups_tidy_donor <- df_age_groups_MSR %>% 
    dplyr::select(sample_type, project_id, mean_MSRD) %>%
    distinct(project_id, sample_type, mean_MSRD) %>%
    mutate(sample_type = factor(sample_type, levels = c( "20-39", "40-59", "60-79" ))) 
  
  
  ggplot(df_age_groups_tidy_donor %>%
           mutate(mean_MSRD = mean_MSRD %>% log10()), 
         aes(project_id, sample_type, fill= mean_MSRD)) + 
    geom_tile()+
    scale_fill_gradient(low="white", high="blue")
  
  
  
  
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", 
                      main_project, "/figures/heatmap_MSR_D_all_tissues.svg")
  ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  #######################################
  ## PLOT - ACCEPTOR
  #######################################
  
  
  df_age_groups_tidy_acceptor <- df_age_groups_MSR %>% 
    dplyr::select(sample_type, project_id, mean_MSRA) %>%
    distinct(project_id, sample_type, mean_MSRA)  %>%
    mutate(sample_type = factor(sample_type, levels = c( "20-39", "40-59", "60-79" ))) 
  
  ggplot(df_age_groups_tidy_acceptor, 
         aes(project_id, sample_type, fill= mean_MSRA)) + 
    geom_tile()+
    scale_fill_gradient(low="white", high="blue")
  
  file_name <- paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", 
                      main_project, "/figures/heatmap_MSR_A_all_tissues.svg")
  ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
}


################################
## CALLS
################################

age_stratification_MSR_changing_with_age(database_path,
                                         gtf_version)



