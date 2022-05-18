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

summarise_lm_tissues <- function(project_id = "BRAIN",
                                 gtf_version = 105,
                                 all_tissues = F) {
  
  if (all_tissues) {
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
                      gene_tpm = tpm,
                      gene_tpm_rct = tpm_median_rct3,
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
                    gene_tpm_rct +
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
                              gene_tpm_rct = model$coefficients["gene_tpm_rct"] %>% unname(),
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
                    gene_tpm_rct +
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
                                 gene_tpm_rct = model$coefficients["gene_tpm_rct"] %>% unname(),
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
  
  all_data <- map_df(all_projects, function(project_id) {
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    
    map_df(all_clusters, function(cluster) { 
      
      folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                            project_id, "/results/pipeline3/missplicing-ratio/",
                            cluster, "/v", gtf_version, "/")
      
      readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
        filter(ref_junID %in% common_introns$ref_junID) %>%
        unnest(gene_name) %>%
        select(ref_junID, ref_type, gene_name) %>%
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
  
  saveRDS(object = df_all, file = file_name)
  
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
  
  df_estimates <- map_df(all_projects, function(project_id) {
    
    all_clusters <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
    map_df(clusters_ID, function(cluster) {
    
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
      
      
      ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
      model$coefficients[-ind_sign] <- 0
      
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
      
      ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
      model$coefficients[-ind_sign] <- 0
      
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
    
  saveRDS(object = df_estimates,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/variance_estimate.rds"))
  
  library(xlsx)
  
  ## Save the .xlsx object
  xlsx::write.xlsx(MSR_Donor, 
                   file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                 "/results/pipeline3/variance_estimate.xlsx"), 
                   sheetName = "MSR_Donor")
  
  xlsx::write.xlsx(MSR_Acceptor, 
                   file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                 "/results/pipeline3/variance_estimate.xlsx"),   
                   sheetName = "MSR_Acceptor",
                   append = T)


  
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
      select(ref_junID, tpm, ref_missplicing_ratio_tissue_ND, ref_missplicing_ratio_tissue_NA) %>%
      distinct(ref_junID, .keep_all = T)
    
    df_introns %>% nrow %>% print
    db_common <- rbind(db_common,
                       df_introns %>% 
                         mutate(tissue = cluster) %>%
                         select(ref_junID, tpm, MSR_D = ref_missplicing_ratio_tissue_ND, MSR_A=  ref_missplicing_ratio_tissue_NA, tissue))
    
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
  
  df_never %>% head() %>% select(gene_id)
  
  
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
  
  df_donor %>% head() %>% select(gene_id)
  
  
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
  
  df_acceptor %>% head() %>% select(gene_id)
  
  
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
  gene_enrichment_both$result %>% select(p_value, term_name)
  
  ## 'Donor' type genes
  gene_enrichment_donor <- gprofiler2::gost(query = donor_genes_unique,
                                            custom_bg = all_genes,
                                            organism = "hsapiens",
                                            ordered_query = F,
                                            correction_method = "bonferroni",
                                            sources = c("GO"))
  gene_enrichment_donor$result %>% select(p_value, term_name)
  
  
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
  gene_enrichment_never$result %>% select(p_value, term_name)
  
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
    select(mean_phastCons20way_5ss, 
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
    select(CDTS_5ss, CDTS_3ss, final_type) %>%
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
    select(seqnames,
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
    select(chrom = seqnames, chromStart = start, chromEnd = end, strand, name = novel_junID)
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
    select(ref_junID, novel_junID, novel_n_individuals, novel_type,
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
    select(ref_junID, ref_start_exon, ref_end_exon)
  
  hits_complete <- x %>%
    filter((ref_junID %in% (hits$ref_junID %>% unique()) & ref_start_exon %in% (hits$ref_start_exon %>% unique())) |
             (ref_junID %in% (hits$ref_junID %>% unique()) & ref_end_exon %in% (hits$ref_end_exon %>% unique())) )
  
  
  
  hits_complete %>% 
    #group_by(ref_junID) %>%
    #mutate(mean_n_individuals = mean(novel_n_individuals)) %>%
    #filter(mean_n_individuals > 90) %>%
    #filter(gene_name == "SNCA") %>% 
    select(gene_name, seqnames, ref_start, ref_end ) %>%
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
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  novel_old_GR <- all_split_reads_details_97 %>%
    filter(type %in% c("novel_donor","novel_acceptor")) %>%
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  none_old_GR <- all_split_reads_details_97 %>%
    filter(type == "none") %>%
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  other_old_GR <- all_split_reads_details_97 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>%
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  
  annotated_new_GR <- all_split_reads_details_104 %>%
    filter(type == "annotated") %>%
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  novel_new_GR <- all_split_reads_details_104 %>%
    filter(type %in% c("novel_donor","novel_acceptor"))  %>%
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  none_new_GR <- all_split_reads_details_104 %>%
    filter(type == "none")  %>%
    select(seqnames, start, end, width, strand) %>% 
    GRanges()
  other_new_GR <- all_split_reads_details_104 %>%
    filter(type %in% c("novel_exon_skip", "novel_combo")) %>%
    select(seqnames, start, end, width, strand) %>% 
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


get_contamination_rates <- function(cluster = "Brain-FrontalCortex_BA9",
                                    folder_root = "/home/sruiz/PROJECTS/splicing-project/results/") {
  
  # cluster <- gtex_tissues[17]
  # all_split_reads_details_81 <- readRDS(file = paste0(folder_root, "base_data/", 
  #                                                     cluster, "/", cluster, "_annotated_SR_details_length_v81.rds"))
  # all_split_reads_details_90 <- readRDS(file = paste0(folder_root, "base_data/", 
  #                                                      cluster, "/", cluster, "_annotated_SR_details_length_v90.rds"))
  
  versions <- c("76", "81", "90", "97")
  dates <- c("18-Jul-2014", "07-Jul-2015", "28-Jul-2017", "26-May-2019")
  
  
  df_contamination <- map_df(versions, function(version) {
    
    
    print(version)
    
    db_introns_old <- readRDS(file = paste0(folder_root, "pipeline3/missplicing-ratio/", 
                                            cluster, "/v", version, "/", cluster, "_db_introns.rds"))
    db_novel_old <- readRDS(file = paste0(folder_root, "pipeline3/missplicing-ratio/", 
                                          cluster, "/v", version, "/", cluster, "_db_novel.rds"))
    
    
    db_introns_new <- readRDS(file = paste0(folder_root, "pipeline3/missplicing-ratio/", 
                                            cluster, "/v104/", cluster, "_db_introns.rds"))
    db_novel_new <- readRDS(file = paste0(folder_root, "pipeline3/missplicing-ratio/", 
                                          cluster, "/v104/", cluster, "_db_novel.rds"))
    
    in_annotation <- db_introns_new %>%
      filter(ref_junID %in% db_novel_old$novel_junID) %>% 
      distinct(ref_junID, .keep_all = T)
    in_annotation %>% nrow()
    
    
    
    out_annotation <- db_novel_new %>%
      filter(novel_junID %in% db_introns_old$ref_junID) %>% 
      distinct(novel_junID, .keep_all = T) 
    out_annotation %>% nrow()
    
    label <- NULL
    
    if (version == "76") {
      label <- paste0("Ensembl_v", version, " (", dates[1], ")")
    } else if (version == "81") {
      label <- paste0("Ensembl_v", version, " (", dates[2], ")")
    } else if (version == "90") {
      label <- paste0("Ensembl_v", version, " (", dates[3], ")")
    } else {
      label <- paste0("Ensembl_v", version, " (", dates[4], ")")
    } 
    
    
    return(data.frame(tissue = cluster,
                      contamination_rates = label,
                      in_annotation = (in_annotation %>% nrow() * 100) / db_novel_old %>% nrow(),
                      out_annotation = (out_annotation %>% nrow() * 100) / db_introns_old %>% nrow()))
    
  })
  
  
  
  ## GETTING CONTAMINATION RATES - % OF INDIVIDUALS
  ggplot(data = df_contamination) +
    #geom_bar(data = novel_104, mapping = aes(x = novel_p_individuals, y = ..count..), stat = "count", color = "black") +
    geom_bar(aes(x = contamination_rates, y = in_annotation), stat = "identity") + 
    xlab(" ") +
    ylab("contamination rates (%)") +
    ggtitle(paste0("Contamination rates across different Ensembl versions compared\nwith Ensembl v104 - ", cluster)) +
    #scale_x_binned() +
    #scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100))+
    #scale_color_manual(values = c("red","black"),
    #                   breaks = c("red","black"),
    #                   labels = c(paste0(new_annotation_104 %>% distinct(ref_junID) %>% nrow(), " fully annotated"),
    #                              paste0(novel_104 %>% distinct(novel_junID) %>% nrow(), " novel"))) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top",
          axis.text.x = element_text(angle = 50, 
                                     vjust = 1,
                                     hjust = 1)) +
    guides(color = guide_legend(title = "Junction type: ", ncol = 2, nrow = 1))
  
  
  
  file_name <- paste0(folder_root, "pipeline3/missplicing-ratio/", cluster, "/images/", cluster, "_contaminationrates.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  ## GETTING CONTAMINATION RATES - NORMALISED COUNTS
  
  
  
  
  
  
}







##########################
## PLOTS #################
##########################

## DISTANCES

plot_distances <- function(cluster,
                           folder_name,
                           limit_bp = 30) {
  
  # cluster <- "Brain-FrontalCortex_BA9"
  # folder_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/"
  # folder_name <- paste0(folder_name, cluster, "/v104/")
  # print(stringr::str_c(Sys.time(), " - Plotting distances for '", cluster, "' samples..."))
  
  
  if (file.exists(paste0(folder_name, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_name, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
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
  
  title <- paste0("Distances - ", cluster)
  

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
    ggtitle(paste0(title)) +
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
  
  
  return(plot)
 
  
  
  
}

plot_distances_PC <- function(cluster,
                              folder_name,
                              limit_bp = 30,
                              stats = F,
                              save_results = F) {
  
  # cluster <- "Brain-FrontalCortex_BA9"
  # folder_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/"
  # folder_name <- paste0(folder_name, cluster, "/v104/")
  # print(stringr::str_c(Sys.time(), " - Plotting distances for '", cluster, "' samples..."))
  
  
  if (file.exists(paste0(folder_name, "/", cluster, "_db_novel.rds"))) {
    df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  } else {
    print(paste0(paste0(folder_name, "/", cluster, "_db_novel.rds"), " file doesn't exist."))
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
  
  title <- paste0("Distances - ", cluster)
  

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
  
  plot_NPC <- ggplot(data = df_novel_tidy %>% filter(type_PC == "non PC")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    ggforce::facet_col(vars(novel_type)) +
    ggtitle(paste0("Non-protein coding (NPC)\n")) +
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
  
  
  plot <- ggarrange(plot_PC, 
                    plot_NPC,
                    common.legend = T,
                    labels = c("A","B"),
                    align = "hv",
                    ncol = 2, 
                    nrow = 1)
  
  plot <- annotate_figure(plot, top = text_grob(title, face = "bold", size = 14))
  
  if (save_results) {
    plot
    ## SAVE THE RESULTS
    file_name <- paste0(cluster, "_distances", limit_bp, "bp.png")
    folder_img <- paste0(folder_name, "/distances/", tissue, "/", "/images/")
    dir.create(file.path(folder_img), showWarnings = T)
    
    file_name <- paste0(folder_img, file_name)
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
    select(cons_5ss_mean = phastCons20way_5ss_mean, 
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
    select(cons_5ss_mean = phastCons20way_5ss_mean, 
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
    select(CDTS_5ss_mean, 
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
                                   brain_tissues = T,
                                   save_results = F) {
  
  if (brain_tissues) {
    df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                         "/results/pipeline3/variance_estimate.rds"))
    graph_title <- "Distribution of the estimate values across 11 brain tissues"
  } else {
    df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_all_tissues.rds"))
    graph_title <- "Distribution of the estimate values across 40 GTEx tissues"
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
  
  
  MSR_tidy <- MSR_tidy %>%
    filter(!(tissue %in% c("CDTS_5ss",
                           "CDTS_3ss",
                           "mean_phastCons20way_5ss",
                           "mean_phastCons20way_3ss")))
  
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
    select(gene_id, cluster) %>%
    distinct(gene_id, .keep_all = T) 
  
  gtex_tpm %>% head()
  gtex_tpm %>% nrow()
  
  
  ## LOAD MIS-SPLICING RATIO DATA FOR THE CURRENT TISSUE
  
  df <- readRDS(file = paste0(folder_name, "/", cluster, "_missplicingratio_tidy_v2.rds"))
  df %>% head()
  df %>% 
    filter(type == "acceptor") %>% 
    select(c(gene_id_start, gene_id_end)) %>%
    unname() %>% unlist() %>% unique() %>%
    length()
  
  genes <- df %>%
    select(c(gene_id_start, gene_id_end)) %>%
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
