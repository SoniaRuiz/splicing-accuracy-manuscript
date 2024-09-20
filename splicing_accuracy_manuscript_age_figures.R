library(tidyverse)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(Biostrings)
library(tidyverse)

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/splicing_accuracy_manuscript_age_figures.R")

##############################################
## CONNECT TO THE AGE STRATIFICATION DATABASE
##############################################

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

#' Title
#' Extracts simple stats from the 'Age stratification' database for the manuscript
#' Manuscript Methods section 'Age stratification and sample clustering'
#' @return
#' @export
#'
#' @examples
age_stratification_get_stats <- function() {
  
  
  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'metadata'")
  
  ## Total number of samples considered
  db_metadata <- dbGetQuery(con, query) %>% as_tibble()
  db_metadata %>%
    distinct(sample_id, .keep_all = T)
  
  
  ## Number of individuals
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
  
  
  ## ANNOTATED INTRONS 
  
  query <- paste0("SELECT * from 'intron'")
  db_introns <- dbGetQuery(con, query) %>% as_tibble()
  
  ## This database contained a total of 321,663 annotated introns
  db_introns %>% distinct(ref_junID) %>% nrow()
  
  ## From which 254,416 presented evidence of at least one type of mis-splicing event. 
  db_introns %>% dplyr::count(misspliced)
  
  
  
  
  ## NOVEL JUNCTIONS
  
  query <- paste0("SELECT novel_junID, novel_type from 'novel'")
  novel <- dbGetQuery(con, query) 
  
  ## Collectively, it included 2,848,776 novel junctions (novel_acceptor = 1,664,788; novel_donor = 1,183,988)
  novel %>% distinct(novel_junID) %>% nrow()
  novel %>% dplyr::count(novel_type)
  
  ## covering 200,837 transcripts, 31,544 genes
  query <- paste0("SELECT gene_id from 'gene'")
  gene <- dbGetQuery(con, query) 
  gene %>% nrow()
  
  query <- paste0("SELECT transcript_id from 'transcript'")
  transcript <- dbGetQuery(con, query) 
  transcript %>% distinct(transcript_id) %>% nrow()

}


#' Title
#' Produces a figure containing summary metadata of the samples included in the "Age Stratification" intron database.
#' Corresponding to 'Supplementary Figure 23'
#' @return
#' @export
#'
#' @examples
age_stratification_plot_metadata <- function() {
  

  tables <- dbListTables(con)
  query <- paste0("SELECT * from 'metadata'")
  db_metadata <- dbGetQuery(con, query) %>% 
    distinct(external_id, .keep_all = T) %>%
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


#' Title
#' Function to obtain the annotated introns that are overlapping age sample clusters and tissues.
#' This method was performed to reduce biases in comparisons across age groups and tissues.
#' @return
#' @export
#'
#' @examples
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


#' Title
#' Calculates the wilcoxon effect size in MSRs from overlapping introns in samples aged '60-79' vs '20-39'yrs.
#' Corresponding to the Results section 'Increasing age is associated with increasing levels of splicing inaccuracies'
#' @return
#' @export
#'
#' @examples
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
    
  }
  
}


#' Title
#' Plots the output of the function 'get_effsize_MSR_with_age()'.
#' In other words, it plots the age effect size between the distributions of MSR values between the overlapping introns in '60-79' vs '20-39'yrs per tissue.
#' Corresponding to the Main Figure (Fig. 7. a)
#' @param effect.size.file.path 
#' @param figure.name 
#' @param plot.stats 
#'
#' @return
#' @export
#'
#' @examples
plot_MSR_age_effsize <- function(effect.size.file.path = paste0(results_folder, "/effsize_MSR_with_age.rds"),
                                 figure.name = paste0(figures_folder, "/age_MSR_effectsize.png"),
                                 plot.stats = F) {
  

  
  #######################################################
  ## LOAD AND TIDY THE DATA
  #######################################################
  
  
  df_wilcoxon_age <- readRDS(file = effect.size.file.path) %>%
    group_by(MSR_type) %>%
    mutate(q = p.adjust(pval, method = "fdr"))%>%
    ungroup()
  
  
  
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


#' Title
#' Obtains the GO and KEGG enrichment of the most mis-spliced genes with age in BRAIN tissue.
#' Manuscript Results Section: Increasing age is associated with increasing levels of splicing inaccuracies
#' Creates Main Figure (Fig. 7. b) 
#' @return
#' @export
#'
#' @examples
age_brain_GO_KEGG_enrichment <- function() {
  
  
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
    
    ## Load never mis-spliced introns
    query <- paste0("SELECT DISTINCT ref_junID, ref_type, MSR_D, MSR_A, transcript_id 
                    FROM '", age_group, "_", project_id, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    
    ## Load mis-spliced introns
    query <- paste0("SELECT DISTINCT ref_junID, ref_type, MSR_D, MSR_A, transcript_id 
                    FROM '", age_group, "_", project_id, "_misspliced'")
    introns <- plyr::rbind.fill(introns, dbGetQuery(con, query)) %>% as_tibble()
    introns <- introns %>% distinct(ref_junID, .keep_all = T)
    
    
    ## Add transcript and gene info
    query <- paste0("SELECT transcript.transcript_id, transcript.id, gene.gene_name, gene.gene_id
                    FROM 'transcript' 
                    INNER JOIN 'gene' ON gene.id = transcript.gene_id
                    WHERE transcript.id IN (", paste(introns$transcript_id, collapse = ","),")")
    genes <- dbGetQuery(con, query) %>% as_tibble()
    
    ## Return
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
  
  
  
  df_age_groups_tidy <- df_age_groups %>% filter(ref_junID %in% common_junctions) %>% drop_na()
  df_age_groups_tidy$ref_junID %>% unique %>% length()
  
  
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
    filter((`20-39` < `40-59` & `40-59` < `60-79`)) 
  df_MSRA_increasing <- df_MSRA %>%
    filter((`20-39` < `40-59` & `40-59` < `60-79`))# | 
  
  
  ## We identified 37,743 annotated introns of interest based on increasing MSRD or MSRA values with age. 
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
  
  ego_MSR %>% as_tibble() %>% filter(Description == "dendritic spine")
  ego_MSR %>% as_tibble() %>% filter(str_detect(Description, pattern = "neuron to neuron synapse")) 
  ego_MSR %>% as_tibble() %>% filter(str_detect(Description, pattern = "tau protein binding")) 
  
  ## Save result
  write.csv(x = ego_MSR %>% as.data.frame() %>% dplyr::select(-"geneID"),
            file = paste0(results_folder, "SuppTable9.csv"),
            row.names = F)
  
}


#################################################################
# CALCULATE AGE EFFECT SIZE IN NORMALISED MSR VALUES
# MSR VALUES WERE NORMALISED BY EXPRESSION LEVELS OF RBP/NMD FACTORS
#################################################################


#' Title
#' Main function designated to calculate the age effect size in MSR values between overlapping introns in '60-79' vs '20-39'yrs.
#' MSR values are normalised by the inverse fold-change TPM of each RBP/NMD gene across the samples of each age group and tissue.
#' This function aims to assess the degree of effect that the expression of each RBP/NMD gene produces over the levels of splicing noise.
#' @param project.list 
#' @param do_paired_test 
#' @param get_median_TPM 
#' @param replace 
#'
#' @return
#' @export
#'
#' @examples
age_stratification_age_effsize_MSR_normalised_by_TPM <- function(project.list = c("BRAIN"),
                                                                 do_paired_test = T,
                                                                 get_median_TPM = T,
                                                                 replace = F) {
  
  
  
  local_results_folder <- paste0(results_folder, "/MSR_normalisation_by_TPM")
  dir.create(path = local_results_folder, recursive = T, showWarnings = F)
  
  file_name <- paste0(results_folder, "/common_introns_all_age_groups.rds")
  df_common_introns_age <- readRDS(file = file_name)
  
  
  ################################
  ## LOAD RBP/NMD GENES OF INTEREST
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


  
  
  
  ################################
  ## PER GTEX TISSUE
  ################################

  
  for (project_id in project.list) {
    
    ## project_id = project.list[1]
    ## project_id <- "BRAIN"
    
    message(project_id, "...")
    
    
    tissue_local_results_folder <- paste0(local_results_folder, "/", project_id, "/", project_id, "_paired", do_paired_test)
    dir.create(path = tissue_local_results_folder, recursive = T, showWarnings = F)
    
    
    if ( !file.exists(paste0(tissue_local_results_folder,  
                             "/", project_id, "_effect_size_paired", do_paired_test, "_median", get_median_TPM,".rds")) || 
         replace ) {
      
      
      #####################################
      ## CALCULATE INVERSE FOLD-CHANGE TPM
      #####################################
      
      if (file.exists(paste0(tissue_local_results_folder, "/", project_id, "_inverse_fold-change_TPM.rds"))) {
        
        fold_change_TPM <- readRDS(file = paste0(tissue_local_results_folder, "/", project_id, "_inverse_fold-change_TPM.rds"))
        
      } else {
        fold_change_TPM <- age_stratification_calculate_fold_change_TPM(project.id = project_id,
                                                                        RBP.list = all_RBPs_tidy,
                                                                        get.median = get_median_TPM)
        
        fold_change_TPM %>% arrange(t_test) %>% as_tibble()
        saveRDS(object = fold_change_TPM,
                file = paste0(tissue_local_results_folder, "/", project_id, "_inverse_fold-change_TPM.rds"))
      }

      
      #####################################
      ## NORMALISE MSR VALUES USING
      ## INVERSE FOLD-CHANGE TPM
      #####################################
      
      
      MSR_normalised <- age_stratification_get_normalised_MSR_values_by_gene_TPM(project.id = project_id,
                                                                                 fold.change.TPM = fold_change_TPM,
                                                                                 gene.list = all_RBPs_tidy)
      
      
      
      ################################
      ## CALCULATE AGE EFFECT SIZES 
      ## USING MSR NORMALISED DATA
      ################################
      
      
      message("Calculating age effect sizes using MSR normalised data....")
      
      
      gene_list <- (fold_change_TPM %>% pull(name))
      

      doParallel::registerDoParallel(cores = 3)
      effect_size_age_normalised_TPM <- foreach(j = seq(length(gene_list)), .combine = "rbind") %dopar% {

        
        # gene <- gene_list[1]
        # gene <- "KHDRBS2"
        # gene <- "UPF1"
        # gene <- "U2AF1"
        
        
        gene <- (gene_list)[j]
 
       
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
                                                   alternative = "less",
                                                   ref.group = "20-39",
                                                   conf.level = .05,
                                                   paired = do_paired_test) %>%
                mutate(MSR_type = MSR_type,
                       project = project_id,
                       gene_normalised = gene)
            
              ## Return results
              return(w_effsize %>%
                       cbind(wilcox_pval))
              
            },
            error = function(cond) {
              message("Here's the original error message:")
              message(conditionMessage(cond))
              return(NULL)
            }
          )
        })
      } 

      
      
      saveRDS(object = effect_size_age_normalised_TPM %>% 
                inner_join(y = all_RBPs_tidy, 
                           by = c("gene_normalised" = "name")),
              file = paste0(tissue_local_results_folder, 
                            "/", project_id, "_effect_size_paired", 
                            do_paired_test, "_median", get_median_TPM,"_inverse_fold-change_TPM.rds"))
      
      
      
    } 
    
    

  }

  
 
  
}

#' Title
#' Auxiliary function to 'age_stratification_age_effsize_MSR_normalised_by_TPM()'
#' Calculates the inverse fold-change TPM of each RBP/NMD gene between the sample group '60-79' as compared to '20-39'yrs.
#' @param project.id 
#' @param RBP.list 
#' @param get.median 
#'
#' @return
#' @export
#'
#' @examples
age_stratification_calculate_fold_change_TPM <- function(project.id,
                                                         RBP.list,
                                                         get.median = T) {
  
  

  
  fold_change_TPM_values_age_groups <- map_df(RBP.list$id, function(RBP) {
    
    # RBP = RBP.list$id[1]
    # RBP = "ENSG00000271885"
    
    message(RBP,"...")
    
    young_group_TPM <- readRDS(file = paste0(base_folder,"/results/splicing_1read/", gtf_version, "/", project.id, "/tpm/",
                                       project.id, "_20-39_tpm.rds")) %>%
      filter(gene_id %in% RBP ) %>% 
      rowwise() %>% 
      mutate(TPM = mean( c_across(where(is.numeric))) ) %>%
      ungroup()
    
    
    eldest_group_TPM <- readRDS(file = paste0(base_folder,"/results/splicing_1read/", gtf_version, "/", project.id, "/tpm/",
                                       project.id, "_60-79_tpm.rds")) %>%
      filter(gene_id %in% RBP ) %>% 
      rowwise() %>% 
      mutate(TPM = mean( c_across(where(is.numeric))) ) %>%
      ungroup() 
    
    
    
    if (young_group_TPM %>% nrow == 1 &&
        eldest_group_TPM %>% nrow == 1 ) {
      
      
      
      fold_change = eldest_group_TPM$TPM/young_group_TPM$TPM
      
      data.frame(gene = RBP,
                 mean_youngest = young_group_TPM$TPM,
                 mean_eldest = eldest_group_TPM$TPM,
                 fold_change_eldest = fold_change,
                 inverse_fold_change = 1/fold_change,
                 t_test = t.test(young_group_TPM[,-1] %>% gather() %>% pull(value),
                                 eldest_group_TPM[,-1] %>% gather() %>% pull(value))$p.value,
                 type = ifelse(fold_change > 1, "upregulation with age", 
                               ifelse(fold_change < 1, "downregulation with age", 
                                      "no change in expression w age"))) %>%
        return()
      
    } else {
      return(NULL)
    }
    
    
    
    
  })
  
  fold_change_TPM_values_age_groups %>%
    inner_join(y = RBP.list ,
               by = c("gene" = "id")) %>%
    arrange(desc(fold_change_eldest)) %>%
    return()
  
}




#' Title
#' Auxiliary function to 'age_stratification_age_effsize_MSR_normalised_by_TPM()'.
#' Normalises MSR values of the overlapping introns between age groups using 
#' the inverse fold-change TPM of each RBP/NMD gene between the sample group '60-79' as compared to '20-39'yrs.
#' @param project.id 
#' @param fold.change.TPM 
#' @param gene.list 
#'
#' @return
#' @export
#'
#' @examples
age_stratification_get_normalised_MSR_values_by_gene_TPM <- function(project.id,
                                                                     fold.change.TPM,
                                                                     gene.list) {

  
  
  
  local_results_folder <- paste0(results_folder, "/MSR_normalisation_by_TPM/")
  
  median_TPM_values_age_groups <- map_df(c("20-39","60-79"), function(cluster_id) {
    
    # cluster_id <- c("20-39","60-79")[1]
    
    message(cluster_id, "...")
    
    # cluster_id <- (df_metadata$cluster %>% unique())[1]
    
    ## 1. Calculate median TPM values corresponding to each RBP/NMD gene across the samples of each age group
    message(cluster_id, " - getting introns from database .... ")
    
   
    
    # 2. Query database to get MSR values across age groups
    con <- dbConnect(RSQLite::SQLite(), database_path)
    query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_misspliced'")
    tissue_introns <- dbGetQuery(con, query) %>% as_tibble() %>% distinct(ref_junID, .keep_all = T)
    
    query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_nevermisspliced' " )
    tissue_introns <- rbind(tissue_introns, dbGetQuery(con, query)) %>% as_tibble() %>% 
      distinct(ref_junID, .keep_all = T) %>%
      arrange(ref_junID)
    DBI::dbDisconnect(con)
  
    
    message("Calculating normalised MSR values by inverse fold-change TPM...")
    
    
    # 3. Normalise MSR values per RBP/NMD gene in the age group and tissue
    normalised_TPM <- map_df(fold.change.TPM$name, function(gene_name) {
      
      # message(gene_name, "...")
      # gene_name <- fold.change.TPM$name[1]
      
      if (cluster_id == "60-79") {
        TPM <- fold.change.TPM %>% filter(name == gene_name) %>% pull(inverse_fold_change)
      } else {
        #TPM <- fold.change.TPM %>% filter(name == gene_name) %>% pull(mean_youngest)
        TPM <- 1
      }
      
      
      
      tissue_introns %>%
        mutate(MSR_D_normalised = MSR_D/TPM,
               MSR_A_normalised = MSR_A/TPM) %>%
        mutate(TPM_group = fold.change.TPM %>% filter(name == gene_name) %>% pull(fold_change_eldest),
               inverse_TPM_group = TPM,
               age = cluster_id,
               gene_normalised = gene_name) %>%
        return()
      
    })
 
    
    normalised_TPM %>% 
      return()
      
      

  })
  
  
  return(median_TPM_values_age_groups)

}


#' Title
#' Auxiliary function to 'age_stratification_age_effsize_MSR_normalised_by_TPM()'.
#' Plots the age effect size using normalised MSR values by the inverse fold-change TPM of each RBP/NMD gene.
#' Corresponding to the supplementary figures 19-21
#' @param effect.size.file.path 
#' @param figure.path 
#' @param plot.stats 
#' @param paired.test 
#' @param get.median 
#' @param project.id 
#'
#' @return
#' @export
#'
#' @examples
plot_age_stratification_effsize_normalised_MSR <- function(effect.size.file.path,
                                                           figure.path,
                                                           plot.stats = F,
                                                           paired.test,
                                                           get.median,
                                                           project.id) {
  
  
  # SET VARIABLES
  # project.id = "BLOOD_VESSEL"
  # project.id = "BLOOD"
  # project.id = "BRAIN"
  # paired.test = T
  # get.median = T
  # effect.size.file.path = paste0(results_folder, "/MSR_normalisation_by_TPM/", project.id, "/", project.id, "_paired",paired.test,"/",
  #                                project.id, "_effect_size_paired",
  #                                paired.test, "_median", get.median,"_inverse_fold-change_TPM.rds")
  # figure.path = paste0(figures_folder, "/age_normalised_MSR_by_TPM/", project.id)
  
  
  df_wilcoxon_age <- readRDS(file = effect.size.file.path) %>%
    group_by(MSR_type) %>%
    mutate(q = 0)%>%#p.adjust(pval, method = "fdr"))%>%
    ungroup()
  
  
  #######################################################
  ## PLOTS
  #######################################################
  
  reference_effsize <- readRDS(file = paste0(results_folder, "/effsize_MSR_with_age.rds")) %>%
    filter(tissue == project.id) %>%
    arrange(desc(MSR_type)) %>%
    pull(effsize) %>%
    mean()
  
  df_wilcoxon_age$MSR_type = factor(df_wilcoxon_age$MSR_type, 
                                    levels = c(df_wilcoxon_age$MSR_type %>% unique))
  
  df_wilcoxon_tidy_final <- df_wilcoxon_age %>%
    mutate(project = str_replace(project, pattern = "_",replacement = " ")) %>%
    group_by(MSR_type) %>%
    mutate(gene_normalised = fct_reorder(gene_normalised, plyr::desc(effsize))) %>%
    ungroup()  %>%
    mutate(delta_effsize = reference_effsize - effsize) %>%
    dplyr::select(delta_effsize, 
                  effsize,
                  gene_normalised,
                  MSR_type, q, 
                  Splicing.regulation, Spliceosome, Exon.Junction.Complex, NMD ) %>%
    mutate(type = ifelse(Splicing.regulation == 1, "Splicing Regulation",
                         ifelse(Spliceosome == 1, "Spliceosome",
                                ifelse(Exon.Junction.Complex == 1, 
                                       "Exon Junction Complex", "NMD")))) %>%
    mutate(color = ifelse(q > 0.05, "grey", ifelse(MSR_type == "MSR Donor", "#34b275", "#64037d"))) %>%
    mutate(type = factor(type, levels=c('NMD','Splicing Regulation','Spliceosome',
                                        'Exon Junction Complex')))
  
  
  
  # write.csv(x = df_wilcoxon_tidy_final %>%
  #             mutate(statistical_test = "Wilcoxon Rank text: rstatix::wilcox_test(data, formula, paired = TRUE, correct = TRUE, alternative = 'less')",
  #                    H0 = "The MSR_D observations from the '20-39' & '60-79' distributions are symmetric about their median value.",
  #                    H1 = "The MSR_D observations from the '20-39' distribution are smaller at their median value than the MSR_D observations from the 60-79' distribution"),
  #           file = paste0(results_folder, "/age_wilcoxon_MSR_all_tissues.csv"),
  #           row.names = F)
  
  
  
  ####################################
  ## PLOT 
  ####################################
  
  
  # reference_effsize = 0.0388
  
  age_effsize_plot  <- ggplot(data = df_wilcoxon_tidy_final,
         aes(x = effsize, y = gene_normalised, colour = color)) +
    geom_point(alpha=.7) +
    geom_vline(mapping = aes(xintercept = reference_effsize), linetype="dotted") +
    facet_grid(~type) +
    theme_light() +
    ylab("") +
    xlab("Probability of superior MSR after controlling for individual RBP/NMD activity in '60-79'yrs compared to '20-39'yrs") +
    scale_color_manual(values = c("#666666", "#34b275", "#64037d"),
                       breaks = c("grey", "#34b275", "#64037d"),
                       labels = c("Non-significant FDR", "MSR Donor", "MSR Acceptor")) +
    #custom_ggtheme + 
    theme( plot.margin = margin(0,0,0,0),
           legend.box.margin=margin(l = -11),
           legend.position="right", 
           legend.box="horizontal", 
           axis.text.y = element_text(size = 6, 
                                      family="Arial", colour = "black"), 
           axis.text.x = element_text(size = 6, 
                                      family="Arial", colour = "black"),
           axis.title.x = element_text(size = 7, 
                                       family="Arial", colour = "black"),
           strip.text.x = element_text(size = 7, 
                                       family="Arial", colour = "black"),
           legend.text = element_text(size = 7, 
                                      family="Arial", colour = "black")) +
    scale_size(name = "FDR pval:",
               #trans="log10",
               
               range=c(5), 
               breaks=c(2.2e-16, 0.5)) +
    
    #guides(size = guide_legend(title = "q:"))+
    guides(color = guide_legend(title = ""))+
    scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
  
  
  
  # Set the colors to be used
  fill_colors <- c("#cccccc","#999999","#999999","#999999")
  
  # Find strips glob
  gt<-ggplot_gtable(ggplot_build(age_effsize_plot))
  strips <- which(startsWith(gt$layout$name,'strip'))
  
  # Change the fill color of each strip
  for (s in seq_along(strips)) {
    gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- fill_colors[s]
  }
  
  dir.create(path = figure.path, recursive = T, showWarnings = F)
  png(filename=paste0(figure.path, "/", project.id, "_age_effect_size_paired",paired.test,"_median",get.median,"_fold-change_TPM.png"), 
      width = 180, height = 240, 
      units = "mm", res = 300)
  plot(gt)
  dev.off()
  
}


################################
## CALLS
################################

# # c("BLOOD_VESSEL", "BLOOD", "BRAIN", "COLON", "MUSCLE", "SKIN")
# age_stratification_age_effsize_MSR_normalised_by_TPM(project.list = c("MUSCLE", "SKIN"),
#                                                      do_paired_test = T,
#                                                      get_median_TPM = T,
#                                                      replace = T)



