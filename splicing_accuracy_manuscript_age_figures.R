library(tidyverse)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(Biostrings)
library(tidyverse)

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/splicing_accuracy_manuscript_age_figures.R")

##################################
## CONNECT TO THE DATABASE
##################################

base_folder <- here::here()
# base_folder <- "/mnt/PROJECTS/splicing-accuracy-manuscript"

gtf_version <- 105

supporting_reads <- 1
project_name <-  paste0("splicing_", supporting_reads, "read")
database_name <- paste0(project_name, "_age")


database_folder <- paste0(base_folder, "/database/", database_name, "/", gtf_version, "/")

figures_folder <- paste0(base_folder, "/results/", project_name, "/", gtf_version, "/_paper_review/figures/")
results_folder <- paste0(base_folder, "/results/", project_name, "/", gtf_version, "/_paper_review/results/")


database_path <- paste0(database_folder,  "/", database_name, ".sqlite")

con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)

query = paste0("SELECT * FROM 'metadata'")
df_metadata <- dbGetQuery(con, query) %>% distinct(external_id, .keep_all = T) %>% as_tibble()
all_projects <- df_metadata$SRA_project %>% unique

age_projects <- readRDS(file = paste0(base_folder, "/results/",  project_name, 
                                      "/", gtf_version, "/all_final_projects_used.rds"))


custom_ggtheme <-  theme(text = element_text(size = 7, family="Arial", colour = "black"),
                         axis.ticks = element_line(colour = "black", linewidth = 2),
                         axis.text = element_text(size = 7, family="Arial", colour = "black"),
                         axis.line = element_line(colour = "black"),
                         axis.title = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.y = element_text(size = 7, family="Arial", colour = "black"),
                         axis.text.x = element_text(size = 7, family="Arial", colour = "black", hjust = 0.5, vjust = 0.5),
                         strip.text = element_text(size = 7, family="Arial", colour = "black"),
                         legend.text = element_text(size = 7, family="Arial", colour = "black"),
                         #legend.title = element_blank(),
                         legend.position = "top",
                         legend.box = "vertical")


##################################
## STATS
##################################

age_stratification_get_stats <- function() {
  
  
  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'metadata'")
  
  ## Total number of samples considered
  db_metadata <- dbGetQuery(con, query) %>% as_tibble()
  db_metadata %>% nrow() 
  
  ## Number of samples
  db_metadata %>%
    distinct(external_id, .keep_all = T)
  
  
  ## Number of samples per age group
  db_metadata %>%
    distinct(external_id, .keep_all = T)%>%
    dplyr::count(cluster)
  
  
  ## Number of tissues considered
  db_metadata %>%
    distinct(external_id, .keep_all = T)%>%
    dplyr::count(SRA_project,cluster)
  
  
  
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
  

  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'metadata'")
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
    custom_ggtheme
  
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
    custom_ggtheme
  
  plot_gender
  
  
  ## RIN
  plot_rin <- ggplot(db_metadata ) +
    geom_density(aes(x = rin, fill = cluster), alpha = 0.5) + 
    theme_light() +
    labs(x = "RIN" ) +
    scale_fill_hue() +
    guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_ggtheme
  
  plot_rin
  
  
  ## Mapped read depth
  plot_read_depth <- ggplot(db_metadata ) +
    geom_density(aes(x = all_mapped_reads, fill = cluster), alpha = 0.5) + 
    theme_light() +
    scale_fill_hue() +
    labs(x = "Mapped Read Count" ) +
    guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
    custom_ggtheme
  
  plot_read_depth
  
  ggpubr::ggarrange(plot_num_samples,plot_gender,
                    plot_rin,plot_read_depth,
                    labels = c("a","b","c","d"))
  
  

  file_name <- paste0(figures_folder, "/age_metadata")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 210, height = 190, units = "mm", dpi = 300)
  
  
}




## 9. Increasing age is associated with increasing levels of splicing errors


get_common_introns_across_age_groups <- function() {
  
  ## CONNECT TO THE DATABASE
  con <- dbConnect(RSQLite::SQLite(), database_path) 
  dbListTables(con)
  query <- paste0("SELECT * FROM 'metadata'")
  
  
  ## Load metadata
  df_metadata <- dbGetQuery(con, query) %>%
    group_by(SRA_project) %>%
    mutate(nsamples = n()) %>%
    filter(nsamples >= 70)
  
  
  #folder_results <- paste0(getwd(), "/results/_paper/results/")
  #dir.create(file.path(folder_results), recursive = TRUE, showWarnings = T)
  
  
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
                        WHERE novel_junID IN (", paste(introns$novel_junID %>% unique(), collapse =","), ")")
        introns <- introns %>%
          inner_join(y = dbGetQuery(con, query) %>% as_tibble(),  by = "novel_junID")
        
        
        ## Add never-misspliced info
        query <- paste0("SELECT ref_junID, ref_type, MSR_D, MSR_A, transcript_id 
                        FROM '", 
                        age_group, "_", project_id, "_nevermisspliced'")
        introns <- plyr::rbind.fill(introns, dbGetQuery(con, query) %>% as_tibble())
        
        
        ## Add transcript info
        query <- paste0("SELECT id, gene_id 
                        FROM 'transcript' WHERE id IN (", 
                        paste(introns$transcript_id %>% unique(), collapse =","), ")")
        introns %>% nrow()
        introns <- introns %>%
          inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
                     by = c("transcript_id" = "id"))
        
        
        ## Add gene info
        query <- paste0("SELECT id, gene_name 
                        FROM 'gene' WHERE id IN (", 
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
  file_name <- paste0(results_folder, "/common_introns_all_age_groups.rds")
  saveRDS(object = df_age_groups_tidy, file = file_name)
  
}


get_effsize_MSR_with_age <- function() {
  

  ## Load common introns across age groups and tissues
  file_name <- paste0(results_folder, "/common_introns_all_age_groups.rds")
  df_common_introns_age <- readRDS(file = file_name)
  
  df_common_introns_age %>% as_tibble()
  

  
  ## MSR_D data preparation
  df_MSRD <- df_common_introns_age %>%
    dplyr::select(ref_junID, sample_type, MSR_D, gene_name, project_id)  %>% 
    mutate(MSR_D = MSR_D %>% round(digits = 4)) %>%
    spread(sample_type, MSR_D) %>%
    as_tibble()
  
  ## MSR_A data preparation
  df_MSRA <- df_common_introns_age %>%
    dplyr::select(ref_junID, sample_type, MSR_A, gene_name, project_id)  %>% 
    mutate(MSR_A = MSR_A %>% round(digits = 4)) %>%
    spread(sample_type, MSR_A) %>%
    as_tibble()

  
  if ( !file.exists(paste0(results_folder, "/effsize_MSR_with_age.rds")) ) {
    
    df_wilcoxon <- map_df( c("MSR Donor", "MSR Acceptor"), function(MSR_type) { #, "MSR_A"
        
      # MSR_type <- "MSR Donor"
      message( "Getting data from: '", MSR_type, "'..." ) 
  
      doParallel::registerDoParallel(10)
      all_samples_age_used <- foreach(i = seq(length(age_projects)), .combine = "rbind") %dopar%{
  
        proj_id <- age_projects[i]
        
        message( proj_id ) 
        
        # proj_id <- age_projects[2]
        # proj_id <- "BRAIN"
        print(proj_id)
        
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
        
      }
      
    })
  
    file_name <- paste0(results_folder, "/effsize_MSR_with_age.rds")
    saveRDS(object = df_wilcoxon, file = file_name)
    
  } else {
    message("Loading data...")
    df_wilcoxon <- readRDS(file = paste0(results_folder, "/effsize_MSR_with_age.rds"))
  }
  
}




stats_effsize_MSR_with_age <- function() {
  
  ## CONNECT TO THE DATABASE
  con <- dbConnect(RSQLite::SQLite(), database_path) 
  dbListTables(con)
  query <- paste0("SELECT * FROM 'metadata'")
  
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

  if ( (age_projects %>% length()) > 1) {
    file_name <- paste0(results_folder, "/effsize_MSR_with_age.rds")
  } else {
    file_name <- paste0(results_folder, "/", main_project, "_df_wilcoxon_effsize_paired_", project_id,".rds")
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
  
  df_wilcoxon_tidy %>%
    dplyr::count(MSR_type)
  
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


age_brain_enrichment <- function() {
  
  
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
    # age_group <- age_groups[3]
    
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
                    FROM 'transcript' 
                    INNER JOIN 'gene' ON gene.id = transcript.gene_id
                    WHERE transcript.id IN (", paste(introns$transcript_id, collapse = ","),")")
    genes <- dbGetQuery(con, query) %>% as_tibble()
    
    introns <- introns %>%
      inner_join(y = genes,
                 by = c("transcript_id" = "id")) %>%
      mutate(sample_type = age_group) %>%
      as_tibble()
    
    return(introns)
  })
 
  ## Get common junctions across age groups in BRAIN
  
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
  
  df_age_groups_tidy$ref_junID %>% unique %>% length()
  
  
  #####################################
  ## TIDY THE DATAFRAME OF COMMON
  ## INTRONS TO GET INCREASING
  ## MSR_D AND MSR_A VALUES WITH AGE
  #####################################
  
  
  ## We identified 37,743 annotated introns of interest based on increasing MSRD or MSRA values with age. 
  
  df_MSRD <- df_age_groups_tidy %>%
    dplyr::select(ref_junID, sample_type, MSR_D, gene_name, gene_id) %>%
    spread(sample_type, MSR_D)
  
  df_MSRA <- df_age_groups_tidy %>%
    dplyr::select(ref_junID, sample_type, MSR_A, gene_name, gene_id) %>%
    spread(sample_type, MSR_A)
  
  ## Introns increasing MSR with age
  df_MSRD_increasing <- df_MSRD %>%
    filter((`20-39` < `40-59` & `40-59` < `60-79`)) #| 
             #(`20-39` == `40-59` & `40-59` < `60-79`) |
             #(`20-39` < `40-59` & `40-59` == `60-79`))
  df_MSRA_increasing <- df_MSRA %>%
    filter((`20-39` < `40-59` & `40-59` < `60-79`))# | 
             #(`20-39` == `40-59` & `40-59` < `60-79`) |
             #(`20-39` < `40-59` & `40-59` == `60-79`))
  
  c(df_MSRD_increasing %>% distinct(ref_junID) %>% pull(),
    df_MSRA_increasing %>% distinct(ref_junID) %>% pull()) %>% 
    unique() %>% 
    length()
  
  ##After assigning these introns to their unique genes (n=12,408), 
  
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
  
  ego_MSR %>%
    as_tibble() %>%
    filter(Description == "dendritic spine")
  ego_MSR %>%
    as_tibble() %>%
    filter(str_detect(Description, pattern = "neuron to neuron synapse")) 
  ego_MSR %>%
    as_tibble() %>%
    filter(str_detect(Description, pattern = "tau protein binding")) 
  
  ## Save result
  write.csv(x = ego_MSR %>% as.data.frame() %>% dplyr::select(-"geneID"),
            file = paste0(results_folder, "SuppTable9.csv"),
            row.names = F)
  
}


plot_effsize_MSR_with_age <- function(effect.size.file.path = paste0(results_folder, "/effsize_MSR_with_age.rds"),
                                      figure.name = paste0(figures_folder, "/age_MSR_effectsize.png"),
                                      plot.stats = F) {
  

  # file_name <- paste0(results_folder, "/effsize_MSR_with_age.rds")
  
  df_wilcoxon_age <- readRDS(file = effect.size.file.path) %>%
    group_by(MSR_type) %>%
    mutate(q = p.adjust(pval, method = "fdr"))%>%
    ungroup()
  
  #######################################################
  ## PLOTS
  #######################################################
  
  df_wilcoxon_age$MSR_type = factor(df_wilcoxon_age$MSR_type, 
                                    levels = c(df_wilcoxon_age$MSR_type %>% unique))
  
  df_wilcoxon_tidy_final <- df_wilcoxon_age %>%
    filter(tissue %in% (df_wilcoxon_age %>%
                           group_by(tissue) %>%
                           #filter(q <= 0.05) %>%
                           ungroup() %>% 
                           pull(tissue)) ) %>%
    mutate(tissue = str_replace(tissue, pattern = "_",replacement = " ")) %>%
    group_by(MSR_type) %>%
    mutate(tissue = fct_reorder(tissue, plyr::desc(effsize))) %>%
    ungroup()  
  
  
  
  # write.csv(x = df_wilcoxon_tidy_final %>%
  #             mutate(statistical_test = "Wilcoxon Rank text: rstatix::wilcox_test(data, formula, paired = TRUE, correct = TRUE, alternative = 'less')",
  #                    H0 = "The MSR_D observations from the '20-39' & '60-79' distributions are symmetric about their median value.",
  #                    H1 = "The MSR_D observations from the '20-39' distribution are smaller at their median value than the MSR_D observations from the 60-79' distribution"),
  #           file = paste0(results_folder, "/age_wilcoxon_MSR_all_tissues.csv"),
  #           row.names = F)
  
  
  
  ####################################
  ## PLOT 
  ####################################
  
  ggplot(data = df_wilcoxon_tidy_final,
         aes(x = effsize, y = tissue, color = MSR_type, size = q)) +
    geom_point(alpha=.7) +
    facet_grid(~MSR_type, scales = "free") +
    theme_light() +
    ylab("") +
    xlab("Probability of superior MSR in 60-79yrs compared to 20-39yrs") +
    scale_color_manual(values = c("#35B779FF","#64037d"),
                       breaks = c("MSR Donor", "MSR Acceptor"),
                       labels = c("MSR Donor", "MSR Acceptor")) +
    custom_ggtheme + 
    theme( plot.margin = margin(0,10,0,0),
           legend.box.margin=margin(b = -11),
           legend.position="top", 
           legend.box="horizontal") +
    scale_size(name = "q:",
               #trans="log10",
               range=c(5, 1), 
               breaks=c(2.2e-16, 0.5)) +
    
    #guides(size = guide_legend(title = "q:"))+
    guides(color = guide_legend(title = ""))+
    scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
  
  ## Save the figure 
  ggplot2::ggsave(filename = figure.name, width = 180, height = 280, units = "mm", dpi = 300)
  
  
  
  if (plot.stats) {
    
    ## STATS - DONOR
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Donor", 
             q < 0.05) %>%
      arrange(effsize)
    
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Donor", 
             q < 0.05) %>%
      pull(effsize) %>% summary
    
    
    
    ## STATS - ACCEPTOR
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Acceptor", 
             q < 0.05) %>%
      arrange(effsize)
    
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Acceptor", 
             q < 0.05) %>%
      pull(effsize) %>% summary
  }
  
  
  
  
  
}


####################################
# SPLICEOSOMAL RBPs AND THEIR TPMs
####################################


age_stratification_age_effsize_MSR_normalised_by_TPM <- function(project.list = c("BRAIN"),
                                                                 do_paired_test = T,
                                                                 get_median_TPM = T,
                                                                 replace = F) {
  
  
  
  local_results_folder <- paste0(results_folder, "/MSR_normalisation_by_TPM")
  dir.create(path = local_results_folder)
  
  file_name <- paste0(results_folder, "/common_introns_all_age_groups.rds")
  df_common_introns_age <- readRDS(file = file_name)
  
  
  
  
  
  ################################
  ## LOAD RBPs OF INTEREST
  ################################
  
  # Only get genes with splicing regulator, spliceosome, Exon.Junction.Complex and NMD functions

  all_RBPs <- (xlsx::read.xlsx(file = paste0(base_folder,'/dependencies/RBPs_subgroups.xlsx'),
                              header = TRUE, sheetIndex = 1) %>% as_tibble() %>% mutate(NMD = 0))[-116,]
  
  all_NMDs <- data.frame(name = c("SMG1", "SMG5", "SMG6", "SMG7", 
                                  "UPF1", "UPF2", "UPF3A"),
                         id = c("ENSG00000157106", "ENSG00000198952", "ENSG00000070366", "ENSG00000116698", 
                                "ENSG00000005007", "ENSG00000151461", "ENSG00000169062"),
                         Splicing.regulation = c(0,0,0,0,0,0,0),
                         Spliceosome = c(0,0,0,0,0,0,0),
                         Exon.Junction.Complex = c(0,0,0,0,0,0,0),
                         NMD = c(1,1,1,1,1,1,1))

  all_RBPs_tidy <- rbind(all_NMDs,all_RBPs)

  # write_csv(x = all_RBPs_tidy %>% arrange(name) %>%
  #             dplyr::rename(gene_symbol = name,
  #                           gene_ensembl = id),
  #           col_names = T,
  #           file = paste0(results_folder,"/list_of_RBPs_NMD_vanNostrand.csv" ))
  
  ################################
  ## PER TISSUE TO THE DATABASE
  ################################
  
  #doParallel::registerDoParallel(cores = 4)
  #foreach(j = seq(length(local_age_projects)), .combine = "rbind") %dopar% {
  for (project_id in project.list) {
    
    ## project_id = project.list[1]
  
    ## project_id <- "ADIPOSE_TISSUE"
    ## project_id <- "ADRENAL_GLAND"
    ## project_id <- "BRAIN"
    
    message(project_id, "...")
    
    
    tissue_local_results_folder <- paste0(local_results_folder, "/", project_id, "/", project_id, "_paired", do_paired_test)
    dir.create(path = tissue_local_results_folder, recursive = T)
    
    
    if ( !file.exists(paste0(tissue_local_results_folder,  
                             "/", project_id, "_effect_size_paired", do_paired_test, "_median", get_median_TPM,".rds")) || 
         replace ) {
      
      
      ################################
      ## CALCULATE MSR VALUES 
      ## NORMALISED BY NMD/TPM VALUES
      ################################
      
      MSR_normalised <- age_stratification_get_normalised_MSR_values_by_gene_TPM(project.id = project_id,
                                                                                 get.median = get_median_TPM,
                                                                                 gene.list = all_RBPs_tidy)
      
      
      ################################
      ## CALCULATE AGE EFFECT SIZES 
      ## USING MSR NORMALISED DATA
      ################################
      
      message("Calculating age effect sizes using MSR normalised data....")
      
      
      
      gene_list <- c(MSR_normalised$gene_normalised) %>% unique() %>% sort()
      
      
     
      doParallel::registerDoParallel(cores = 10)
      effect_size_age_normalised_TPM <- foreach(j = seq(length(gene_list)), .combine = "rbind") %dopar% {

        
        # gene <- gene_list[1]
        
        
        gene <- gene_list[j]
 
       
        message(gene, "....")
        

        map_df(c("MSR Donor","MSR Acceptor"), function(MSR_type) {
          
          tryCatch(
            {
              # MSR_type = "MSR Donor"
              
              message(gene, " --> ", MSR_type, " --> Wilcoxon paired=", do_paired_test, " effect-size test ....")
              
              if (MSR_type == "MSR Donor") {
                
                MSR_normalised_data <- MSR_normalised %>%
                  dplyr::select(ref_junID, MSR_normalised = MSR_D_normalised, age, gene_normalised) 
                
              } else {
                MSR_normalised_data <- MSR_normalised %>%
                  dplyr::select(ref_junID, MSR_normalised = MSR_A_normalised, age, gene_normalised) 
              }
              
              MSR_normalised_data <- MSR_normalised_data  %>%
                filter(ref_junID %in% (df_common_introns_age$ref_junID %>% unique)) %>%
                filter(gene_normalised == gene) %>%
                spread(key = age, MSR_normalised)
              
              
              ## Wilcoxon one-tailed test
              wilcox_pval <- data.frame(pval = c((wilcox.test(x = MSR_normalised_data$`20-39`,
                                                              y = MSR_normalised_data$`60-79`,
                                                              alternative = "less",
                                                              paired = do_paired_test))$p.value))
              
              
              ## Wilcoxon effect size
              w_effsize <- rstatix::wilcox_effsize(data = MSR_normalised_data %>%
                                                     dplyr::select(ref_junID, `20-39`,`60-79`) %>%
                                                     gather(key = age, value = MSR, -ref_junID) %>%
                                                     mutate(age = age %>% as.factor()),
                                                   formula = MSR ~ age,
                                                   paired = do_paired_test) %>%
                mutate(MSR_type = MSR_type, 
                       project = project_id,
                       gene_normalised = gene)
              
              
              ## Return results
              return(w_effsize %>% cbind(wilcox_pval))
              
            },
            error = function(cond) {
              message("Here's the original error message:")
              message(conditionMessage(cond))
              return(NULL)
            }
          )
        })
      } 
      #})

      
      saveRDS(object = effect_size_age_normalised_TPM %>% 
                inner_join(y = all_RBPs_tidy, 
                           by = c("gene_normalised" = "name")),
              file = paste0(tissue_local_results_folder, 
                            "/", project_id, "_effect_size_paired", do_paired_test, "_median", get_median_TPM,".rds"))
      
      
      
    } else {
      
      effect_size_age_normalised_TPM <- readRDS(file = paste0(tissue_local_results_folder, 
                                                              "/", project_id, "_effect_size_paired", 
                                                              do_paired_test, "_median", get_median_TPM,".rds")) %>%
        as_tibble()
    }
    
    
    ################################
    ## PLOT AGE EFFECT SIZES 
    ## USING MSR NORMALISED DATA
    ################################
    
    plot_age_stratification_effsize_normalised_MSR(effect.size.file.path = paste0(tissue_local_results_folder, 
                                                                                  "/", project_id, "_effect_size_paired", 
                                                                                  do_paired_test, "_median", get_median_TPM,".rds"),
                                                   figure.path = paste0(figures_folder, "/age_normalised_MSR_by_TPM/", project_id),
                                                   paired.test = do_paired_test,
                                                   get.median = get_median_TPM,
                                                   project.id = project_id)
  }

  
 
  
}


age_stratification_get_normalised_MSR_values_by_gene_TPM <- function(project.id,
                                                                     get.median = T,
                                                                     gene.list) {

  
  
  
  local_results_folder <- paste0(results_folder, "/MSR_normalisation_by_TPM/")
  
  # if ( !file.exists(paste0(local_results_folder, "/",project.id,"/",project.id,"_age_MSR_D_normalised_by_TPM.rds")) || 
  #      !file.exists(paste0(local_results_folder, "/",project.id,"/",project.id,"_age_MSR_D_normalised_by_TPM.rds")) ) {
    
    median_TPM_values_age_groups <- map_df(c("20-39","60-79"), function(cluster_id) {
      
      # cluster_id <- c("20-39","60-79")[1]
      
      message(cluster_id, "...")
      
      # cluster_id <- (df_metadata$cluster %>% unique())[1]
      
      ## 1. Calculate median TPM values corresponding to each RBP/NMD gene across the samples of each age group
      print(paste0(Sys.time(), " - ", cluster_id, " ... "))
      
      if (file.exists(paste0(base_folder,"/results/splicing_1read/", gtf_version, "/", project.id, "/tpm/",
                             project.id, "_", cluster_id, "_tpm.rds"))) {
        
        
        message("Calculating median TPM values across sample groups...")
        
        local_TPM <- readRDS(file = paste0(base_folder,"/results/splicing_1read/", gtf_version, "/", project.id, "/tpm/",
                                           project.id, "_", cluster_id, "_tpm.rds")) %>%
          filter(gene_id %in% (gene.list$id) )
        
        
        if (get.median) {
          local_TPM_w_median <- local_TPM %>% 
            rowwise() %>% 
            mutate(TPM = median( c_across(where(is.numeric))) )
        } else {
          local_TPM_w_median <- local_TPM %>% 
            rowwise() %>% 
            mutate(TPM = mean( c_across(where(is.numeric))) )
        }
        local_TPM_w_median <- local_TPM_w_median %>%
          dplyr::select(gene_id, TPM)  %>%
          inner_join(y = gene.list, by = c("gene_id"="id")) %>%
          dplyr::rename(gene_name = name) %>%
          dplyr::relocate(gene_name, .after = "gene_id")
        
        
        # 2. Query database to get MSR values across age groups
        con <- dbConnect(RSQLite::SQLite(), database_path)
        query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_misspliced'")
        tissue_introns <- dbGetQuery(con, query) %>% as_tibble() %>% distinct(ref_junID, .keep_all = T)
        
        query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_nevermisspliced' " )
        tissue_introns <- rbind(tissue_introns, dbGetQuery(con, query)) %>% as_tibble() %>% 
          distinct(ref_junID, .keep_all = T) %>%
          arrange(ref_junID)
        DBI::dbDisconnect(con)
      
        
        message("Calculating normalised MSR values by median TPMs...")
        
        # 3. Normalise MSR values per RBP/NMD gene in the age group and tissue
        normalised_TPM <- map_df(local_TPM_w_median$gene_name, function(gene) {
          
          # message(gene, "...")
          # gene <- local_TPM_w_median$gene_name[1]
          
          TPM <- local_TPM_w_median %>% filter(gene_name == gene) %>% pull(TPM)
          
          tissue_introns %>%
            mutate(MSR_D_normalised = MSR_D/TPM,
                   MSR_A_normalised = MSR_A/TPM) %>%
            mutate(age = cluster_id,
                   gene_normalised = gene) %>%
            return()
          
        })
        
        
        normalised_TPM[which(!is.finite(normalised_TPM$MSR_D_normalised)),"MSR_D_normalised"] <- 0
        normalised_TPM[which(!is.finite(normalised_TPM$MSR_A_normalised)),"MSR_A_normalised"] <- 0
        
        normalised_TPM %>% 
          return()
        
        
      } else {
        return(NULL)
      }   
    })
    
    
    # ## 4. Get only the annotated introns overlapping the two age groups
    # common_ref_junID <- median_TPM_values_age_groups %>%
    #   dplyr::count(ref_junID) %>%
    #   filter(n == ((median_TPM_values_age_groups$gene_normalised %>% unique %>% length) * 2)) %>%
    #   pull(ref_junID)
    # 
    # 
    # ## 5. Filter the results by common introns between the two age groups
    # median_TPM_values_age_groups_common_introns <- median_TPM_values_age_groups %>%
    #   filter(ref_junID %in% common_ref_junID) 
    
    
    # ## 6. Calculate MSR_D normalised values & save results
    # MSR_D_normalised <- median_TPM_values_age_groups_common_introns %>%
    #   dplyr::select(-c(MSR_D, MSR_A, MSR_A_normalised)) %>%
    #   spread(key = age, MSR_D_normalised)
    # saveRDS(object = MSR_D_normalised,
    #         file = paste0(local_results_folder, "/",project.id,"/", project.id,"_age_MSR_D_normalised_by_TPM.rds"))
    # 
    # 
    # ## 7. Calculate MSR_A normalised values & save results
    # MSR_A_normalised <- median_TPM_values_age_groups_common_introns %>%
    #   dplyr::select(-c(MSR_D,MSR_A,MSR_D_normalised)) %>%
    #   spread(key = age, MSR_A_normalised)
    # saveRDS(object = MSR_A_normalised,
    #         file = paste0(local_results_folder, "/",project.id,"/", project.id,"_age_MSR_A_normalised_by_TPM.rds"))
    
    
  # } else {
  #   message(project.id, " - loading normalised MSR values....")
  #   MSR_D_normalised <- readRDS(file = paste0(local_results_folder, 
  #                                             "/", project.id, "/", project.id,"_age_MSR_D_normalised_by_TPM.rds"))
  #   MSR_A_normalised <- readRDS(file = paste0(local_results_folder, 
  #                                             "/", project.id, "/", project.id,"_age_MSR_A_normalised_by_TPM.rds"))
  # }
  
  return(median_TPM_values_age_groups)
}




age_stratification_get_representative_effsize <- function(project.list = c("BLOOD","BLOOD_VESSEL","BRAIN","COLON","MUSCLE","SKIN"),
                                                          paired.test = T,
                                                          get.median = T) {
  
  
  for (project_id in project.list) {
    
    # project_id = project.list[1]
    
    message(project_id, "...")
    
    ref_effsize <- readRDS(file = paste0(results_folder, "/effsize_MSR_with_age.rds")) %>%
      filter(tissue == project_id) %>%
      arrange(desc(MSR_type)) %>%
      dplyr::select(ref_effsize = effsize, MSR_type)
    
    
    effect.size.folder.path <- paste0(results_folder, "/MSR_normalisation_by_TPM/", project_id, "/", project_id, "_paired", paired.test, "/")
    
    df_wilcoxon_age <- readRDS(file = paste0(effect.size.folder.path,
                                             "/", project_id, "_effect_size_paired", paired.test, "_median", get.median,".rds")) %>%
      group_by(MSR_type) %>%
      mutate(q = p.adjust(pval, method = "fdr"))%>%
      ungroup() %>%
      left_join(y = ref_effsize,
                by = "MSR_type") #%>%
      #filter(q<0.05)
    
    ### CALCULATE REFERENCE EFFSIZE PER RBP TYPE
    
    effsize_summary_categories <- map_df(c("Splicing.regulation","Spliceosome","Exon.Junction.Complex","NMD"), function(functional_category) {
      
      # functional_category <- "Splicing.regulation"
      
      df_wilcoxon_age %>%
        filter(if_any(.cols = all_of(functional_category),
                      .fns = function(x){ x == 1})) %>%
       
        group_by(MSR_type) %>%
        mutate(adding_effsize = sum(effsize),
               IQR_effsize = IQR(effsize),
               mean_effsize = mean(effsize),
               effsize_sd = sd(effsize),
               signed_effsize = effsize - ref_effsize) %>%
        ungroup() %>%
        group_by(MSR_type) %>%
        mutate(signed_adding_effsize = sum(signed_effsize))%>%
        ungroup() %>%
        
        group_by(MSR_type) %>%
        mutate(n_genes = n())%>%
        ungroup() %>%
        
        dplyr::select(MSR_type,
                      effsize, ref_effsize, signed_effsize, signed_adding_effsize,
                      gene_normalised,
                      adding_effsize, IQR_effsize, mean_effsize, effsize_sd,
                      all_of(functional_category) ,n_genes) %>%
        #distinct(MSR_type, .keep_all = T) %>%
        gather(key = type, value = functional_category,-c(effsize, ref_effsize, signed_effsize, signed_adding_effsize,
                                                          gene_normalised,
                                                          MSR_type,adding_effsize, IQR_effsize, mean_effsize, effsize_sd, n_genes)) %>%
        dplyr::select(-functional_category) %>%
        return()
      
    })
    
    

    
    data_to_plot <- effsize_summary_categories %>%
      group_by(type) %>%
      mutate(type_tidy = paste0(type, " (", max(n_genes), " genes)")) %>%
      ungroup() %>%
      mutate(MSR_type = MSR_type %>% as.factor(),
             type_tidy = type_tidy %>% as.factor()) %>% 
      replace(is.na(.), 0)
    
    
    
    ##################################
    ## PLOTS
    ##################################
    
    figure_path = paste0(figures_folder, "/age_normalised_MSR_by_TPM/", project_id)
    dir.create(path = figure_path, recursive = T)
    
    ## Distribution of effect sizes
    ggplot2::ggplot(data = data_to_plot) +
      geom_boxplot(mapping = aes(x = effsize, y  = type_tidy, color = MSR_type), show.legend = FALSE) +
      facet_grid(~ fct_rev(MSR_type)) +
      labs(y = '',
           x = "Wilcoxon effect size",
           #title = paste0(project_id,  
            #              '\nSummary of distribution of Wilcoxon effect sizes per functional category'),
           caption = 'Effect sizes measure differences in splicing noise between samples aged 20-39 vs 60-79yrs-old.') +
      theme_light()
    
    
    ggplot2::ggsave(filename = paste0(figure_path, "/", project_id, "_IQR_effsize_paired",paired.test,"_median",get.median,".png"), 
                    width = 180, height = 90, units = "mm", dpi = 300)
    
    
    
    
    ## Adding effect sizes
    
    ggplot2::ggplot(data = data_to_plot %>% 
                      group_by(type) %>%
                      distinct(MSR_type, .keep_all = T) %>%
                      ungroup()) +
      geom_bar(mapping = aes(y = type_tidy, x  = signed_adding_effsize, fill = MSR_type), show.legend = FALSE,
               stat = "identity" ) +
      facet_grid(~ fct_rev(MSR_type)) +
      labs(y = '',
           x = "Cumulative Wilcoxon Effect Size",
           #title = paste0(project_id,  
           #              '\nSummary of distribution of Wilcoxon effect sizes per functional category'),
           caption = 'Effect sizes measure differences in splicing noise between samples aged 20-39 vs 60-79yrs-old.') +
      theme_light()
    
    
    saveRDS(object = data_to_plot %>% 
              filter(type == "NMD") %>%
              #group_by(type) %>%
              #distinct(MSR_type, .keep_all = T) %>%
              #ungroup() %>%
              dplyr::select(MSR_type, gene = gene_normalised, effsize, ref_effsize, signed_effsize, signed_adding_effsize, type_tidy) %>%
              arrange(desc(MSR_type)),
            file = paste0(effect.size.folder.path, 
                          "/", project_id, "_CUMULATIVE_effsize_paired",paired.test,"_median",get.median,".rds"))
    
    ggplot2::ggsave(filename = paste0(figure_path, "/", project_id, "_CUMULATIVE_effsize_paired",paired.test,"_median",get.median,".png"), 
                    width = 180, height = 90, units = "mm", dpi = 300)
    
    # ggplot(data = data_to_plot,
    #        aes(x = ref_effsize,
    #            y  = type_tidy,
    #            color = MSR_type)) +
    #   geom_point(shape = 15, size  = 4) +
    #   
    #   geom_errorbar(aes(xmin = ref_effsize - ref_effsize_sd,
    #                     xmax = ref_effsize + ref_effsize_sd,
    #                     width = 0.15)) +
    #   facet_grid(~ fct_rev(MSR_type)) +
    #   theme_light() +
    #   ylab("") +
    #   xlab("Mean effect-size across the RBPs from each group") +
    #   scale_color_manual(values = c("#34b275", "#64037d"),
    #                      breaks = c("MSR Donor", "MSR Acceptor"),
    #                      labels = c("MSR Donor", "MSR Acceptor")) +
    #   custom_ggtheme + 
    #   theme( plot.margin = margin(0,10,0,0),
    #          legend.box.margin=margin(b = -11),
    #          legend.position="top", 
    #          legend.box="horizontal")  +
    #   
    #   #guides(size = guide_legend(title = "q:"))+
    #   guides(color = guide_legend(title = ""))+
    #   scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
    
    
    
    # ggplot(data = data_to_plot,
    #        aes(x = type,
    #            y = ref_effsize,
    #            colour = MSR_type)) +
    #   geom_point(alpha=.7, size = 3) +
    #   
    #   facet_grid(~ MSR_type) +
    #   theme_light() +
    #   ylab("") +
    #   xlab("Mean effect-size across the RBPs from each group") +
    #   scale_color_manual(values = c("#34b275", "#64037d"),
    #                      breaks = c("MSR Donor", "MSR Acceptor"),
    #                      labels = c("MSR Donor", "MSR Acceptor")) +
    #   custom_ggtheme + 
    #   theme( plot.margin = margin(0,10,0,0),
    #          legend.box.margin=margin(b = -11),
    #          legend.position="top", 
    #          legend.box="horizontal")  +
    #   
    #   #guides(size = guide_legend(title = "q:"))+
    #   guides(color = guide_legend(title = ""))+
    #   scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
    
  }
  
}

# project.id = "ADIPOSE_TISSUE"
# project.id = "ADRENAL_GLAND"
# project.id = "BRAIN"
# effect.size.file.path = paste0(results_folder, "/MSR_normalisation_by_TPM/",project.id, "/", project.id,"_pairedTRUE/", project.id,"_effect_size_normalised_MSR_with_age_pairedTRUE.rds")
# figure.path = paste0(figures_folder, "/age_normalised_MSR_by_TPM/",project.id)

plot_age_stratification_effsize_normalised_MSR <- function(effect.size.file.path,
                                                           figure.path,
                                                           plot.stats = F,
                                                           paired.test,
                                                           get.median,
                                                           project.id) {
  
  
  
  # MSR_type_to_plot <- "MSR Donor"
  
  #for (MSR_type_to_plot in c("MSR Donor", "MSR Acceptor")) {
    
  # file_name <- paste0(results_folder, "/effsize_MSR_with_age.rds")
  
  # if (MSR_type_to_plot == "MSR Donor") {
  #   plot_color = "#34b275"
  # } else {
  #   plot_color = "#64037d"
  # }
  
  
  df_wilcoxon_age <- readRDS(file = effect.size.file.path) %>%
    group_by(MSR_type) %>%
    mutate(q = p.adjust(pval, method = "fdr"))%>%
    ungroup()
  
  
  #######################################################
  ## PLOTS
  #######################################################
  
  df_wilcoxon_age$MSR_type = factor(df_wilcoxon_age$MSR_type, 
                                    levels = c(df_wilcoxon_age$MSR_type %>% unique))
  
  df_wilcoxon_tidy_final <- df_wilcoxon_age %>%
    mutate(project = str_replace(project, pattern = "_",replacement = " ")) %>%
    group_by(MSR_type) %>%
    mutate(gene_normalised = fct_reorder(gene_normalised, plyr::desc(effsize))) %>%
    ungroup()  #%>%
    #filter(MSR_type == MSR_type_to_plot)
  
  
  
  # write.csv(x = df_wilcoxon_tidy_final %>%
  #             mutate(statistical_test = "Wilcoxon Rank text: rstatix::wilcox_test(data, formula, paired = TRUE, correct = TRUE, alternative = 'less')",
  #                    H0 = "The MSR_D observations from the '20-39' & '60-79' distributions are symmetric about their median value.",
  #                    H1 = "The MSR_D observations from the '20-39' distribution are smaller at their median value than the MSR_D observations from the 60-79' distribution"),
  #           file = paste0(results_folder, "/age_wilcoxon_MSR_all_tissues.csv"),
  #           row.names = F)
  
  
  
  ####################################
  ## PLOT 
  ####################################
  
  reference_effsize <- readRDS(file = paste0(results_folder, "/effsize_MSR_with_age.rds")) %>%
    filter(tissue == project.id) %>%
    arrange(desc(MSR_type)) %>%
    pull(effsize) %>% mean()
  
  ggplot(data = df_wilcoxon_tidy_final %>%
           dplyr::select(effsize, gene_normalised, MSR_type, q, 
                         Splicing.regulation, Spliceosome, Exon.Junction.Complex, NMD ) %>%
           mutate(type = ifelse(Splicing.regulation == 1, "Splicing.regulation",
                                ifelse(Spliceosome == 1, "Spliceosome",
                                       ifelse(Exon.Junction.Complex == 1, "Exon.Junction.Complex","NMD")))) %>%
           mutate(color = ifelse(q > 0.05, "grey", ifelse(MSR_type == "MSR Donor", "#34b275", "#64037d"))),
         aes(x = effsize, y = gene_normalised, #size = q, 
             colour = color)) +
    geom_point(alpha=.7, size = 3) +
    geom_vline(mapping = aes(xintercept = reference_effsize), linetype="dotted") +
    facet_grid(~type) +
    theme_light() +
    ylab("") +
    xlab("Probability of superior MSR in 60-79yrs compared to 20-39yrs") +
    scale_color_manual(values = c("#666666", "#34b275", "#64037d"),
                       breaks = c("grey", "#34b275", "#64037d"),
                       labels = c("Non-significant FDR", "MSR Donor", "MSR Acceptor")) +
    custom_ggtheme + 
    theme( plot.margin = margin(0,10,0,0),
           legend.box.margin=margin(b = -11),
           legend.position="top", 
           legend.box="horizontal") +
    scale_size(name = "FDR pval:",
               #trans="log10",
               range=c(5), 
               breaks=c(2.2e-16, 0.5)) +
    
    #guides(size = guide_legend(title = "q:"))+
    guides(color = guide_legend(title = ""))+
    scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
  
  ## Save the figure 
  dir.create(path = figure.path, recursive = T)
  ggplot2::ggsave(filename = paste0(figure.path, "/", project.id, "_age_effect_size_paired",paired.test,"_median",get.median,".png"), 
                  width = 180, height = 280, units = "mm", dpi = 300)
    
  #}
  
  
  
  
  
}


################################
## OTHER FUNCTIONS
################################


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



################################
## CALLS
################################


age_stratification_age_effsize_MSR_normalised_by_TPM(project.list = c("BLOOD_VESSEL"),
                                                     do_paired_test = T,
                                                     get_median_TPM = T,
                                                     replace = F)



