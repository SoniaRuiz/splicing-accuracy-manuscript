library(tidyverse)
library(data.table)
library(GenomicRanges)
library(ggforce)
library(DBI)

## source("/home/sruiz/PROJECTS/splicing-project-recount3/pipeline4-2_age_stratification.R")


source("/home/sruiz/PROJECTS/splicing-project-recount3/init.R")


gtf_version <- 105
main_project <- "age_subsampled"
database_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v",
                        gtf_version, "/", main_project, "/", main_project, ".sqlite")

all_projects <- readRDS(file = "~/PROJECTS/splicing-project-results/splicing-recount3-projects/all_projects_used.rds")
age_projects <- all_projects


##################################
## CALLS - PREPARE AGE DATA
##################################

age_stratification_init_data <- function (projects_id,
                                          gtf_version,
                                          main_project) {
  
  # project_id <- "BRAIN"
  map_df(projects_id, function(project_id) {
    
    print(paste0(Sys.time(), " --> ", project_id, " - getting sample's metadata."))
    
    if ( !file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", project_id,
                             "/v", gtf_version, "/",
                             main_project, "_project/raw_data/samples_metadata.rds")) ) {
      
      rse <- recount3::create_rse_manual(
        project = project_id,
        project_home = "data_sources/gtex",
        organism = "human",
        annotation = "gencode_v29",
        type = "gene")
      
      
      sample_metadata <- rse %>% 
        SummarizedExperiment::colData() %>%
        as_tibble() %>%
        filter(gtex.smrin >= 6.0,
               gtex.smafrze != "EXCLUDE")
      
      folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", project_id, 
                            "/v", gtf_version, "/",
                            main_project, "_project/raw_data/")
      
      dir.create(file.path(folder_root), recursive = TRUE, showWarnings = T)
      saveRDS(object = sample_metadata,
              file = paste0(folder_root, "/samples_metadata.rds"))
      
      rm(rse)
      gc()
      
    } else {

      sample_metadata <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", project_id,
                                               "/v", gtf_version, "/",
                                               main_project, "_project/raw_data/samples_metadata.rds"))
    }
    
    age_samples_clusters <- age_stratification_clustering(sample_metadata, project_id)
    
    age_samples_clusters_tidy <- age_samples_clusters %>%
      mutate(age_group = "")
    
    for (age_cluster in (age_samples_clusters$age %>% unique())) {
      
      # age_cluster = (age_samples_clusters$age %>% unique())[1]
      print(age_cluster)
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
    
    
    if ( (age_samples_clusters_tidy %>%
          dplyr::count(age_group) %>%
          drop_na() %>%
          filter(n >= 25) %>%
          pull(n)) %>% length() == 3 ) {
    
      df_subsampled <- map_df(c("40-59", "60-79"), function(group) {
        
        df_youngest <- age_samples_clusters_tidy %>%
          filter(age_group == "20-39")
        
        df_samples <- age_samples_clusters_tidy %>%
          filter(age_group == group)
        
        data_combined <- rbind(df_youngest, df_samples)
        
        m.out <- MatchIt::matchit(age_group ~ rin, 
                                  data = data_combined, 
                                  distance = data_combined$rin,
                                  method = "nearest", 
                                  caliper = c(rin = 0), 
                                  std.caliper = FALSE)
        MatchIt::match.data(m.out) %>%
          return()
        
      })
      
      age_samples_clusters_tidy <- df_subsampled %>%
        distinct(individual, .keep_all = T) 
      
      # kruskal.test(mapped_read_count ~ age_group, data = age_samples_clusters_tidy) 
      # kruskal.test(sex ~ age_group, data = age_samples_clusters_tidy) 
      # kruskal.test(rin ~ age_group, data = age_samples_clusters_tidy) 
      
      return(age_samples_clusters_tidy)
    } else {
      return(NULL)
    }
  })
  
  
}

age_stratification_clustering <- function (sample_metadata, project_id) {
  
  print(paste0(Sys.time(), " --> ", project_id, " - performing the sample clustering."))
  
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
                   region = person$gtex.smtsd,
                   project = project_id) %>% return()
      }
      
    }
    
  }) 
  
  gtex_age_clustering %>% 
    return()

}


##################################
## CALLS - JUNCTION ANNOTATION
##################################

age_stratification_annotate <- function (age_groups,
                                         main_project,
                                         project_id,
                                         age_samples_clusters) {
  
  ## Loop
  for (age_cluster in age_groups) {
    
    # age_cluster <- age_groups[1]
    # age_cluster <- age_groups[2]
    
    ## SET RESULTS FOLDER
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", project_id, 
                          "/v", gtf_version, "/",
                          main_project, "_project/")
    folder_name <- paste0(folder_root, "/results/", age_cluster, "/")
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
        filter(age_group == age_cluster, 
               region == cluster) %>% 
        pull(individual)
      
      if (samples %>% length() > 0) {
        
        samples_age <- c(samples_age, samples)
        print(paste0(samples %>% length(), " samples... "))
        
        ## Load the annotated split reads
        all_split_reads_details_105 <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                                                             project_id, "/v", gtf_version, "/splicing_project/raw_data/",
                                                             project_id,  "_", cluster, "_all_split_reads_sample_tidy.rds"))
        
        # Load split read counts
        split_read_counts <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
                                                   project_id, "/v", gtf_version, "/splicing_project/raw_data/",
                                                   project_id,  "_", cluster, "_split_read_counts_sample_tidy.rds"))
        
        if ( is.null(names(split_read_counts)) ) {
          split_read_counts <- split_read_counts %>%
            as_tibble(rownames = "junID")
        }
        
        split_read_counts <- split_read_counts %>%
          dplyr::select(c("junID", any_of(samples)))
        
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
      
    }
    
    ## SAVE RESULTS
    print(paste0(Sys.time(), " - SAVING RESULTS ... "))
    
    split_read_counts_age <- split_read_counts_age[rowSums(is.na(split_read_counts_age[,2:ncol(split_read_counts_age)])) != 
                                                     (ncol(split_read_counts_age)-1),]
    

    all_split_reads_details_105_age <- all_split_reads_details_105_age %>%
      data.table::as.data.table() %>%
      distinct(junID, .keep_all = T) %>%
      filter(junID %in% (split_read_counts_age$junID %>% unique()))
    
    
    saveRDS(object = samples_age, file = paste0(folder_root, "/raw_data/",
                                                project_id, "_", age_cluster,
                                                "_samples_used.rds"))
    saveRDS(object = split_read_counts_age %>% 
              data.table::as.data.table(), 
            file = paste0(folder_root, "/raw_data/",
                          project_id, "_", age_cluster, 
                          "_split_read_counts_sample_tidy.rds"))
    
    saveRDS(object = all_split_reads_details_105_age %>% data.table::as.data.table(), 
            file = paste0(folder_root, "/raw_data/",
                          project_id, "_", age_cluster, "_all_split_reads_sample_tidy.rds"))
    
    print(paste0(Sys.time(), " - RESULTS SAVED!"))
    
    rm(all_split_reads_details_105_age)
    rm(split_read_counts_age)
    rm(samples_age)
    gc()
    
  }
  
  rm(all_split_reads_details_105_age)
  rm(all_split_reads_details_105)
  rm(split_read_counts)
  rm(split_read_counts_age)
  rm(samples_age)
  gc()
  
}


##################################
## CALLS - JUNCTION PAIRING
##################################

age_stratification_junction_pairing <- function (age_groups,
                                                 project_id,
                                                 age_samples_clusters) {
  
  
  
  ## Loop
  for (age_cluster in age_groups) {
    
    # age_cluster <- age_groups[1]
    # age_cluster <- age_groups[3]
    
    print(paste0(Sys.time(), " - Starting analysis for '", age_cluster, "' samples..."))
    
    ## SET RESULTS FOLDER
    folder_root <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/", 
                          project_id, "/v", gtf_version, "/", main_project, "_project/")
    folder_name <- paste0(folder_root, "/results/", age_cluster)
    dir.create(file.path(folder_name), showWarnings = T, recursive = T)
    
    
    #######################################################################
    ## GENERATE THE IDB
    #######################################################################
    
    
    ## Load data ----------------------------------------------------------------
    
    samples_age <- readRDS(file = paste0(folder_root, "/raw_data/",
                                         project_id, "_", age_cluster, "_samples_used.rds"))
    
    split_read_counts_age <- readRDS(file = paste0(folder_root, "/raw_data/",
                                                   project_id, "_", age_cluster, 
                                                   "_split_read_counts_sample_tidy.rds")) 
    
    if ( is.null(names(split_read_counts_age)) ) {
      split_read_counts_age <- split_read_counts_age %>%
        as_tibble("junID")
    }
    
    all_split_reads_details_105_age <- readRDS(file = paste0(folder_root, "/raw_data/",
                                                             project_id, "_", age_cluster, 
                                                             "_all_split_reads_sample_tidy.rds")) %>% 
      as_tibble()
    
    
    ## Call functions ----------------------------------------------------------------
    
    get_distances(cluster = age_cluster,
                  samples = samples_age,
                  split_read_counts = split_read_counts_age,
                  all_split_reads_details = all_split_reads_details_105_age,
                  folder_name)
    
    
    extract_distances(cluster = age_cluster,
                      samples = samples_age,
                      split_read_counts = split_read_counts_age,
                      folder_name = folder_name)
    
    
    get_never_misspliced(cluster = age_cluster,
                         samples = samples_age,
                         split_read_counts = split_read_counts_age,
                         all_split_reads_details = all_split_reads_details_105_age,
                         folder_name = folder_name,
                         save_results = T)
    
    # get_distances(cluster = age_cluster,
    #               samples = samples_age,
    #               split_read_counts = split_read_counts_age,
    #               all_split_reads_details = all_split_reads_details_105_age,
    #               folder_name = folder_name)
    # 
    # 
    # extract_distances(cluster = age_cluster,
    #                   samples = samples_age,
    #                   split_read_counts = split_read_counts_age,
    #                   folder_name = folder_name)
    # 
    # get_never_misspliced(cluster = age_cluster,
    #                      split_read_counts = split_read_counts_age,
    #                      all_split_reads_details = all_split_reads_details_105_age,
    #                      samples = samples_age,
    #                      folder_name = folder_name)

    # extract_never_misspliced(cluster = age_cluster,
    #                          split_read_counts = split_read_counts_age,
    #                          samples = samples_age,
    #                          folder_name = folder_name)
    # 
    # add_never_misspliced_to_df(cluster = age_cluster,
    #                            all_split_reads_details = all_split_reads_details_105_age,
    #                            samples = samples_age,
    #                            folder_name = folder_name)
    # 
    # 
    # ## Get mis-splicing ratio
    # 
    # folder_idb_name <- paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/" ,
    #                           project_id, "/results/pipeline3/missplicing-ratio/age/", age_cluster)
    # dir.create(file.path(folder_idb_name), showWarnings = F, recursive = T)
    # 
    # get_missplicing_ratio(cluster = age_cluster,
    #                       split_read_counts = split_read_counts_age,
    #                       samples = samples_age,
    #                       folder_name = folder_name,
    #                       folder_save_name = folder_idb_name)
    # 
    # add_missplicing_class_to_df(cluster = age_cluster,
    #                             folder_name = folder_idb_name)
    
    ############################################
    ## ADD FEATURES TO THE IDB
    ############################################
    
    # remove_MT_genes(cluster = age_cluster,
    #                 folder_name = folder_idb_name)
    # 
    # add_intron_type(cluster = age_cluster,
    #                 folder_name = folder_idb_name)
    # 
    # clinvar_analysis(cluster = age_cluster,
    #                  folder_name = folder_idb_name)
    
    # add_MANE_info(cluster = age_cluster,
    #               folder_name = folder_idb_name)
    
    
    ############################################
    ## QC
    ############################################

    
    # get_missplicing_QC(cluster = age_cluster,
    #                    split_read_counts = split_read_counts_age,
    #                    all_split_reads_details = all_split_reads_details_105_age,
    #                    samples = samples_age,
    #                    folder_name = folder_idb_name)
    
    
    
    # }
  }
  
  rm(samples_age)
  rm(split_read_counts_age)
  rm(all_split_reads_details_105_age)
  gc()
}


##################################
## CALLS - SQL GENERATION
##################################




## Only if the database does not exist, we run this code to generate it from scratch

# if ( !file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/",
#                          main_project, "/", main_project , ".sqlite")) ) {
#   
#   project_init <- age_stratification_init_data(projects_id = age_projects,
#                                                gtf_version = gtf_version,
#                                                main_project = main_project) %>%
#     as_tibble()
#   
#   age_projects <- project_init$project %>% unique()
#   
#   for (project_id in age_projects) {
# 
#     
#     # project_id <- age_projects[1]
#     # project_id <- "BRAIN"
#     # project_id <- "MUSCLE"
#     # project_id <- "LIVER"
#     # project_id <- "KIDNEY"
#     
#     print(paste0(Sys.time(), " --> ", project_id))
# 
#     project_init_local <- project_init %>%
#       filter(project == project_id)
# 
#     if ( !file.exists(paste0("/home/sruiz/PROJECTS/splicing-project-results/splicing-recount3-projects/",
#                              project_id, "/v", gtf_version, "/", main_project, "_project/results/",
#                              (project_init_local$age_group %>% unique())[1], "/",
#                              (project_init_local$age_group %>% unique())[1], "_raw_distances_tidy.rds")) ) {
# 
# 
# 
#       age_stratification_annotate(age_groups = project_init_local$age_group %>% unique(),
#                                   project_id = project_id,
#                                   main_project = main_project,
#                                   age_samples_clusters = project_init_local)
# 
# 
#       age_stratification_junction_pairing(age_groups = project_init_local$age_group %>% unique(),
#                                           project_id = project_id,
#                                           age_samples_clusters = project_init_local)
#       
#       
#     } else {
#       print(paste0("File '", (project_init_local$age_group %>% unique())[1], "_raw_distances_tidy.rds' exists!"))
#     }
# 
# 
# 
#   }
#   
#   
#   
#   if ( ! exists("project_init") ) {
#     project_init <- age_stratification_init_data(projects_id = age_projects,
#                                                  gtf_version = gtf_version,
#                                                  main_project = main_project) %>%
#       as_tibble()
#   }
# 
#   get_all_annotated_split_reads(all_projects = age_projects,
#                                 gtf_version = gtf_version,
#                                 all_clusters = project_init$age_group %>% unique(),
#                                 main_project = main_project)
# 
#   get_all_raw_distances_pairings(all_projects = age_projects,
#                                  gtf_version = gtf_version,
#                                  all_clusters = project_init$age_group %>% unique(),
#                                  main_project = main_project)
# 
#   tidy_data_pior_sql(projects_used = age_projects,
#                      gtf_version = gtf_version,
#                      all_clusters = project_init$age_group %>% unique(),
#                      main_project = main_project)
#   
#   
# 
# }

project_init <- age_stratification_init_data(projects_id = age_projects,
                                             gtf_version = gtf_version,
                                             main_project = main_project) %>%
  as_tibble()

age_projects <- project_init$project %>% unique()

sql_database_generation(database_path = database_path,
                        projects_used = age_projects,
                        main_project = main_project,
                        gtf_version = gtf_version,
                        remove_all = T)





