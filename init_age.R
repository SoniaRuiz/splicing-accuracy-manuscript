library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(Biostrings)
library(tidyverse)
library(protr)

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/init_age.R")

setwd(normalizePath("."))

dependencies_folder <- paste0(here::here(), "/dependencies/")
data_source <- "data_sources/gtex"
base_folder <- here::here()

#####################################
## LOAD SOURCE SCRIPTS
#####################################


setwd(file.path(base_folder,"scripts"))
files.sources = list.files()
sapply(files.sources, source)

setwd(file.path(base_folder))


#####################################
## SET MAIN VARIABLES
#####################################


gtf_version <- 105

supportive_reads <- 1
replace <- T
project_name <- paste0("splicing_",supportive_reads,"read")
database_name <- paste0(project_name, "_age")

database_folder <- paste0(getwd(), "/database/", database_name, "/", gtf_version, "/")
results_folder <- file.path(here::here("results"), paste0(project_name, "/", gtf_version, "/"))

age_projects <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))




## AGE SAMPLE CLUSTERING ---------------------------------------------------------------------

age_stratification_init_data <- function (projects.id,
                                          gtf.version,
                                          project.name,
                                          results.folder) {
  
  # project_id <- "BRAIN"
  map_df(projects.id, function(project_id) {
    
    message(Sys.time(), " --> ", project_id, " - getting sample's metadata.")
    
    ## Samples are loaded from the splicing project
    sample_metadata <- readRDS(file = file.path(results.folder, project_id, "/base_data/", 
                                                paste0(project_id, "_samples_metadata.rds")) )
    ## Cluster the samples by age
    age_samples_clusters_tidy <- age_stratification_clustering(sample_metadata, project_id) %>%
      mutate(age_group = "")
    
    for (age_cluster in (age_samples_clusters_tidy$age %>% unique())) {
      
      # age_cluster = (age_samples_clusters_tidy$age %>% unique())[1]
      #message(age_cluster, "...")
      
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
               indx <- which(age_samples_clusters_tidy$age == age_cluster)
               age_samples_clusters_tidy[indx, "age_group"] <- NA
               print('not found')
             }
      )
    }
    
    age_samples_clusters_tidy <- age_samples_clusters_tidy %>%
      mutate(rin_transformed = rin / 10) %>%
      mutate(age_group_trans = (age_group %>% as.factor %>% as.integer())/3)
    
    ## If the tissue has the three age groups with at least 25 samples
    if ( (age_samples_clusters_tidy %>%
          dplyr::count(age_group) %>%
          drop_na() %>%
          filter(n >= 25) %>%
          pull(n)) %>% length() == 3 ) {
    
      
      # "40-59" ----------------------------------------------

      data_combined <- age_samples_clusters_tidy %>%
        filter(age_group %in% c("20-39","40-59")) 

      m.out0 <- MatchIt::matchit(age_group_trans ~ rin_transformed,
                                 data = data_combined,
                                 method = "optimal",
                                 ratio = 1, replace = F)
      
      st_match <- MatchIt::match.data(m.out0) %>%
        dplyr::select(-c(distance,weights, subclass))
      
      # m.out <- MatchIt::matchit(formula = age_group ~ rin_transformed, 
      #                           data = age_samples_clusters_tidy, 
      #                           distance = age_samples_clusters_tidy$rin_transformed,
      #                           method = "nearest", 
      #                           caliper = c(rin_transformed = 0.005), 
      #                           std.caliper = FALSE, 
      #                           ratio = 1)
      # 
      # st_match <- MatchIt::match.data(m.out) %>%
      #   dplyr::select(-c(distance,weights, subclass))
      
      
      
      
      # "60-79" ----------------------------------------------
      
      df_youngest <- st_match %>%
        filter(age_group == "20-39") 
      
      df_eldest <- age_samples_clusters_tidy %>%
        filter(age_group == "60-79")
      
      data_combined <- rbind(df_youngest, df_eldest)
      

      m.out <- MatchIt::matchit(age_group_trans ~ rin_transformed,
                                data = data_combined,
                                method = "optimal",
                                ratio = 1, 
                                replace = F)

      nd_match <- MatchIt::match.data(m.out) %>%
        dplyr::select(-c( distance, weights, subclass))
      
      

      # Join both ----------------------------------------------
     

      age_samples_clusters_subsampled <- rbind(st_match, nd_match) %>%
        distinct(individual, .keep_all = T) 
      
      message(project_id, " --> ", age_samples_clusters_subsampled %>%
                dplyr::count(age_group))
      

      

      return(age_samples_clusters_subsampled)
      
    } else {
      return(NULL)
    }
  })
  
}

age_stratification_clustering <- function (sample_metadata, 
                                           project.id) {
  
  print(paste0(Sys.time(), " --> ", project.id, " - performing the sample clustering."))
  
  ## 1. Get key metadata
  
  clusters <- sample_metadata$cluster %>% unique()
  age_groups <- sample_metadata$age %>% unique()
  samples <- sample_metadata$external_id %>% unique()
  
  ## 2. Loop through the individuals
  
  gtex_age_clustering <- map_df(samples, function(individual) {
    
    # individual <- (samples)[1]
    
    person <- sample_metadata %>%
      filter(external_id == individual)
    
    # print(person$smtsd %>% unique)
    
    if ( person$age %>% unique() %in% age_groups ) {
      
      if ( any(person$cluster %in% clusters) ) {
        
        sex <- ifelse(person$gender == 1, "male", "female")
        
        data.frame(individual = person$external_id %>% unique(),
                   age = person$age %>% unique(),
                   rin = person$rin,
                   mapped_read_count = person$all_mapped_reads,
                   sex = sex,
                   #sample_recount_id = person$gtex.sampid,
                   cluster = person$cluster,
                   project = project.id) %>% return()
      }
      
    }
    
  }) 
  
  gtex_age_clustering %>% 
    return()

}



## JUNCTION ANNOTATION -----------------------------------------------------------------

age_stratification_annotate <- function (age.groups,
                                         project.name,
                                         gtf.version,
                                         project.id,
                                         age.samples.clusters,
                                         results.folder,
                                         supportive.reads) {
  
  
  all_samples_used <- NULL
  local_results_folder <- file.path(results.folder, "/base_data/")
  
  ## Loop per age cluster
  doParallel::registerDoParallel(3)
  all_samples_age_used <- foreach(i = seq(length(age.groups))) %dopar%{
    
    age_cluster <- age.groups[i]
    message(age_cluster)
    
    # age_cluster <- age.groups[1]
    # age_cluster <- age.groups[2]
    
    
    split_read_counts_age <- NULL
    all_split_reads_age <- NULL
    samples_age <- NULL
    j <- 1
    
    ## GENERATE THE COUNTS AND ANNOTATED SPLIT READ COUNTS BY AGE CLUSTER ---------------------------------------------------
    
    for ( cluster_id in age.samples.clusters$cluster %>% unique() ) {
      
      # cluster_id <- (age.samples.clusters$cluster %>% unique())[1]
      
      print(paste0(Sys.time(), " - ", cluster_id, ", ", age_cluster))
      
      ## Get the samples from the current tissue according to the current age cluster
      samples <- age.samples.clusters %>% 
        filter(age_group == age_cluster, cluster == cluster_id) %>% 
        pull(individual)
      
      
      if ( samples %>% length() > 0 ) {
        
        samples_age <- c(samples_age, samples)
        print(paste0(samples %>% length(), " samples... "))
        
        ## Load the annotated split reads (the splicing project contains all samples QC'ed and checked with recount3 origin)
        all_split_reads_local <- readRDS(file = paste0(local_results_folder, "/", project.id,  
                                                       "_", cluster_id, "_all_split_reads.rds"))

        # Load split read counts
        split_read_counts_local <- readRDS(file = paste0(local_results_folder, "/", project.id,  
                                                         "_", cluster_id, "_split_read_counts.rds"))

        
        if ( is.null(names(split_read_counts_local)) ) {
          split_read_counts_local <- split_read_counts_local %>%
            as_tibble(rownames = "junID")
        }
        
        if ( split_read_counts_local %>%
             mutate(total=rowSums(select_if(., is.numeric))) %>% 
             filter(total < supportive.reads) %>% 
             nrow() != 0 ) {
          paste0(project.id, " - There are split reads with less than the number of supportive reads indicated by parameter!") %>% print
        }

        split_read_counts_local <- split_read_counts_local %>%
          dplyr::select(c("junID", any_of(samples)))
        
        ## Only split reads with at least 2 supportive reads after filtering using the local cluster of samples
        split_read_counts_local <- split_read_counts_local %>%
          mutate(total = rowSums(select_if(., is.numeric))) %>% 
          filter(total >= supportive.reads) %>%
          dplyr::select(-total)
        
        
        all_split_reads_local <- all_split_reads_local %>%
          filter(junID %in% split_read_counts_local$junID)

        
        if ( j == 1 ) {
          split_read_counts_age <- split_read_counts_local
          all_split_reads_age <- all_split_reads_local
          j <- 2

        } else {

          ## The age stratification analysis is done per body site category,
          ## thus, only annotated introns overlapping the clusters from each
          ## body site are considered

          all_split_reads_age <- data.table::rbindlist(list(all_split_reads_local,
                                                            all_split_reads_age)) %>%
            distinct(junID, .keep_all = T)

          split_read_counts_age <- split_read_counts_local %>%
            inner_join(y = split_read_counts_age,
                       by = "junID")
        }

        
        #print(paste0(split_read_counts_age %>% names() %>% length(), " samples processed in total."))
        
        rm(all_split_reads_local)
        rm(split_read_counts_local)
        rm(samples)
        gc()
      
      }
      
    }
    

    
    print(paste0(Sys.time(), " - SAVING RESULTS ... "))
    
    split_read_counts_age <- split_read_counts_age[rowSums(is.na(split_read_counts_age[,2:ncol(split_read_counts_age)])) !=
                                                     (ncol(split_read_counts_age)-1),]

    all_split_reads_age <- all_split_reads_age %>%
      distinct(junID, .keep_all = T) %>%
      filter( junID %in% (split_read_counts_age$junID %>% unique()) )
    
    
    
    ############################################
    ## SAVE RESULTS FOR THE CURRENT AGE CLUSTER
    ############################################
    
    saveRDS(object = samples_age,
            file = paste0(local_results_folder, "/", project.id, "_", age_cluster, "_samples_used.rds"))

    ## Split read counts for age
    saveRDS(object = split_read_counts_age %>% data.table::as.data.table(),
            file = paste0(local_results_folder, "/", project.id, "_", age_cluster, "_split_read_counts.rds"))

    ## All split reads for age
    saveRDS(object = all_split_reads_age %>% data.table::as.data.table(),
            file = paste0(local_results_folder, "/", project.id, "_", age_cluster, "_all_split_reads.rds"))
    
    print(paste0(Sys.time(), " - ", age_cluster, " - ", project.id, " - RESULTS SAVED!"))
    
    
    return(samples_age)
    
  }
  
  all_samples_age_used %>% head()
  saveRDS(object = all_samples_age_used %>% unlist,
          file = paste0(local_results_folder, "/", project.id, "_age_all_samples_used.rds"))
  
  
  metadata_project <- readRDS(file = paste0(local_results_folder, "/", project.id, "_samples_metadata.rds") )
  
  ## Save age-related metadata
  saveRDS(object = metadata_project %>% 
            filter(external_id %in% (all_samples_age_used %>% unlist)) %>%
            inner_join(y = age.samples.clusters %>% dplyr::select(individual, age_group),
                       by = c("external_id" = "individual")) %>%
            dplyr::select(-cluster) %>%
            dplyr::rename(cluster = age_group),
          file = paste0(local_results_folder, "/", project.id, "_age_samples_metadata.rds"))
  
  
  ## Save age-related clusters
  saveRDS(object = age.groups %>% unique(),
          file = paste0(local_results_folder, "/", project.id, "_age_clusters_used.rds"))
  

  gc()
  
}



## JUNCTION PAIRING --------------------------------------------------------------------

age_stratification_junction_pairing <- function (age.groups,
                                                 project.id,
                                                 project.name,
                                                 gtf.version,
                                                 age.samples.clusters,
                                                 results.folder,
                                                 replace) {

  

  ## Loop per age cluster
  doParallel::registerDoParallel(3)
  
  foreach(i = seq(length(age.groups)), .combine = "rbind") %dopar%{
    
    age_cluster <- age.groups[i]
    
    # age_cluster <- age.groups[1]
    # age_cluster <- age.groups[2]
    
    print(paste0(Sys.time(), " - Starting analysis for '", project.id, "' and '", age_cluster, "' samples..."))
    
    #######################################################################
    ## GENERATE THE IDB
    #######################################################################
    
    
    ## Load data ----------------------------------------------------------------
    
    samples_age <- readRDS(file = paste0(results.folder, "/base_data/", 
                                         project.id, "_", age_cluster, "_samples_used.rds"))
    
    split_read_counts_age <- readRDS(file = paste0(results.folder, "/base_data/", 
                                                   project.id, "_", age_cluster, "_split_read_counts.rds")) 
    
    if ( is.null(names(split_read_counts_age)) ) {
      split_read_counts_age <- split_read_counts_age %>%
        as_tibble("junID")
    }
    
    all_split_reads_age <- readRDS(file = paste0(results.folder, "/base_data/",
                                                 project.id, "_", age_cluster, "_all_split_reads.rds")) %>% 
      as_tibble()
    
    if (!identical(split_read_counts_age$junID,
                               all_split_reads_age$junID)) {
      paste0(project.id, " - some jxn have not been considered!") %>% print()
    }
    
    ## Call functions ----------------------------------------------------------------
    
    folder_results <- paste0(results.folder, "/junction_pairing/", age_cluster, "/")
    dir.create(file.path(folder_results), showWarnings = T, recursive = T)
    
    get_distances(cluster = age_cluster,
                  samples = samples_age,
                  split_read_counts = split_read_counts_age,
                  all_split_reads_details = all_split_reads_age,
                  folder_name = folder_results,
                  replace = replace)
    
    
    extract_distances(cluster = age_cluster,
                      samples = samples_age,
                      folder_name = folder_results,
                      replace = replace)
    
    
    get_never_misspliced(cluster = age_cluster,
                         samples = samples_age,
                         split_read_counts = split_read_counts_age,
                         all_split_reads_details = all_split_reads_age,
                         folder_name = folder_results,
                         replace = replace)
    
  }
  
  rm(split_read_counts_age)
  rm(all_split_reads_age)
  rm(samples_age)
  gc()
}


##################################
## CALLS - SQL GENERATION
##################################




## Only if the database does not exist, we run this code to generate it from scratch

#if ( !file.exists(paste0(database_folder, "/", database_name , ".sqlite")) ) {

  if (!exists("project_init")) {
    project_init <- age_stratification_init_data(projects.id = age_projects,
                                                 gtf.version = gtf_version,
                                                 project.name = project_name,
                                                 results.folder = results_folder) %>%
      as_tibble()
  }

  age_projects <- project_init$project %>% unique()

  # for (project_id in age_projects) {
  # 
  #   # project_id <- age_projects[2]
  #   print(paste0(Sys.time(), " --> ", project_id))
  #   project_init_local <- project_init %>% filter(project == project_id)
  #   project_results_folder <- file.path(results_folder, project_id)
  #   
  #   
  #   
  #   message(project_id)
  #   
  #   age_stratification_annotate(age.groups = project_init_local$age_group %>% unique(),
  #                               project.id = project_id,
  #                               gtf.version = gtf_version,
  #                               project.name = project_name,
  #                               age.samples.clusters = project_init_local,
  #                               results.folder = project_results_folder,
  #                               supportive.reads = supportive_reads)
  #   
  #   
  #   age_stratification_junction_pairing(age.groups = project_init_local$age_group %>% unique(),
  #                                       project.id = project_id,
  #                                       gtf.version = gtf_version,
  #                                       project.name = project_name,
  #                                       age.samples.clusters = project_init_local,
  #                                       results.folder = project_results_folder,
  #                                       replace)
  # }
  

  # get_all_annotated_split_reads(recount3.project.IDs = age_projects,
  #                               all.clusters = project_init$age_group %>% unique(),
  #                               database.folder = database_folder,
  #                               results.folder = results_folder)
  # 
  # 
  # get_all_raw_jxn_pairings(recount3.project.IDs = age_projects,
  #                          all.clusters = project_init$age_group %>% unique(),
  #                          database.folder = database_folder,
  #                          results.folder = results_folder)
  
   
   
  age_projects <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))
   
   
  # generate_recount3_tpm(recount3.project.IDs = age_projects,
  #                       data.source = data_source,
  #                       results.folder = results_folder)
  # 
  # 
  # 
  # tidy_data_pior_sql(recount3.project.IDs = age_projects,
  #                    project.name = project_name,
  #                    all.clusters = project_init$age_group %>% unique(),
  #                    database.folder = database_folder,
  #                    results.folder = results_folder)
  


  database_path <- paste0(database_folder,  "/", database_name, ".sqlite")

  sql_database_generation(database.path = database_path,
                          recount3.project.IDs = age_projects,
                          project.name = project_name,
                          gtf.version = gtf_version,
                          remove.all = F,
                          database.folder = database_folder,
                          results.folder = results_folder,
                          supportive.reads = supportive_reads)


#}


  ####################################
  ## Visualise sample pairing metadata
  ####################################
  
  # ggplot(project_init %>% 
  #          dplyr::count(project,age_group)%>%
  #          arrange(age_group , n) %>%
  #          mutate(project = fct_inorder(project)))+
  #   geom_bar(aes(y = project, x = n, fill = age_group),
  #            stat = "identity", position = position_dodge()) +
  #   custom_theme
  # 
  # ggplot(project_init %>% 
  #          mutate(sex = sex %>% as.character()) %>%
  #          dplyr::count(project, sex) %>%
  #          arrange(sex , n) %>%
  #          mutate(project = fct_inorder(project)) ) +
  #   geom_bar(aes(y = project, x = n, group = sex, fill = sex),
  #            stat = "identity", position = "dodge") + 
  #   theme_light()+
  #   custom_theme
  # 
  # 
  # ggplot(project_init ) +
  #   geom_density(aes(x = rin, fill = age_group), alpha = 0.5) + 
  #   theme_light() +
  #   labs(x = "RIN" ) +
  #   scale_fill_hue() +
  #   guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
  #   custom_theme
  # 
  # 
  # ggplot(project_init ) +
  #   geom_density(aes(x = mapped_read_count, fill = age_group), alpha = 0.5) + 
  #   theme_light() +
  #   scale_fill_hue() +
  #   labs(x = "Mapped Read Count" ) +
  #   guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
  #   custom_theme
  