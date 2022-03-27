####################################################
## RUNNING ANALYSES USING THE INTRON DATABASE ######
####################################################

## DISTANCES ----------------------------------------

project_id <- "BRAIN"

summarise_distances_tissues <- function() {
  
  if (!exists(gtex_tissues)) {
    gtex_tissues <-  readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                                           project_id, "/raw_data/all_clusters_used.rds"))
  }
  
  all_data <- map_df(gtex_tissues, function(tissue) {
    
    print(tissue)
    
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
                          project_id, "/results/pipeline3/missplicing-ratio/",
                          tissue, "/v105/")
    df <- readRDS(file = paste0(folder_name, "/", tissue, "_db_novel.rds")) %>%
      as_tibble()
    
    df_tidy_donor <- df %>%
      filter(novel_type == "novel_donor") %>%
      pull(distance) %>%
      summary() %>% as.array %>%as.data.frame() %>%
      mutate(tissue = tissue) %>%
      mutate("NovelType" = "Novel Donor")
    
    df_tidy_donor <- rbind(df_tidy_donor,
                           data.frame(Var1 = "Mode (+)",
                                      Freq = get_mode(df %>%
                                                        filter(novel_type == "novel_donor", distance > 0) %>%
                                                        pull(distance)),
                                      tissue = tissue,
                                      "NovelType" = "Novel Donor"))
    
    df_tidy_donor <- rbind(df_tidy_donor,
                           data.frame(Var1 = "Mode (-)",
                                      Freq = get_mode(df %>%
                                                        filter(novel_type == "novel_donor", distance < 0) %>%
                                                        pull(distance)),
                                      tissue = tissue,
                                      "NovelType" = "Novel Donor"))
    ##############
    
    df_tidy_acceptor <- df %>%
      filter(novel_type == "novel_acceptor") %>%
      pull(distance) %>%
      summary() %>% 
      as.array %>%
      as.data.frame() %>%
      mutate(tissue = tissue) %>%
      mutate("NovelType" = "Novel Acceptor")
    
    df_tidy_acceptor <- rbind(df_tidy_acceptor,
                              data.frame(Var1 = "Mode (+)",
                                         Freq = get_mode(df %>%
                                                           filter(novel_type == "novel_acceptor", distance > 0) %>%
                                                           pull(distance)),
                                         tissue = tissue,
                                         "NovelType" = "Novel Acceptor"))
    
    df_tidy_acceptor <- rbind(df_tidy_acceptor,
                              data.frame(Var1 = "Mode (-)",
                                         Freq = get_mode(df %>%
                                                           filter(novel_type == "novel_acceptor", distance < 0) %>%
                                                           pull(distance)),
                                         tissue = tissue,
                                         "NovelType" = "Novel Acceptor"))
    
    
    
    return(rbind(df_tidy_donor, df_tidy_acceptor))
  })
  
  saveRDS(object = all_data,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/",
          project_id, "/results/pipeline3/summary_distances.rds"))
}

## MODULO -------------------------------------------

summarise_modulo_tissues <- function() {
  
  folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/"
  ## Load GTEx tissues
  tissues <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")
  tissues_tidy <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used_tidy.rda")
  samples <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_used.rda")
  
  ## Import MANE annotation
  mane_transcripts <- rtracklayer::import(con = "/home/sruiz/PROJECTS/splicing-project/data/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf") %>% 
    as.data.frame()  %>%
    filter(tag == "MANE_Select" | tag == "MANE_Plus_Clinical") %>%
    mutate(transcript_id = str_remove(string = transcript_id, pattern = "\\..*") ) %>%
    distinct(transcript_id)
  
  
  df_modulo_prop <- #map_df(c(-1,1), function(type) {
    
    #print(type)
    
    map_df(gtex_tissues, function(cluster) {
      
      # cluster <- gtex_tissues[1]
      
      i <- which(cluster == gtex_tissues)
      print(tissues_tidy[i])
      
      folder_name <- paste0(folder_root, "/", cluster, "/v104/")
      
      ## Load IDB intron data -------------------------------------------------
      file_path <- paste0(folder_name, "/", cluster, "_db_introns.rds")
      df_introns <- readRDS(file = file_path) 
      df_introns[1,]
      
      df_introns_mane <- df_introns %>%
        unnest(tx_id_junction) %>%
        filter(tx_id_junction %in% mane_transcripts$transcript_id)
      
      ## Load IDB novel data -------------------------------------------------
      file_path <- paste0(folder_name, "/", cluster, "_db_novel.rds")
      db_novel <- readRDS(file = file_path)
      df_novel_mane <- db_novel %>%
        filter(ref_junID %in% df_introns_mane$ref_junID)
      
      
      ## Get Modulo ----------------------------------------------------------
      
      df_novel_mane <- df_novel_mane %>%
        dplyr::rowwise() %>%
        filter(abs(distance) <= 100) %>%
        mutate(modulo = abs(distance) %% 3) %>%
        distinct(novel_junID, .keep_all = T) 
      
      num_novel_donor <- df_novel_mane %>% filter(novel_type == "novel_donor") %>% nrow()
      num_novel_acceptor <- df_novel_mane %>% filter(novel_type == "novel_acceptor") %>% nrow()
      
      ## Number of modulo 0 unique introns
      modulo0_donor <- (df_novel_mane %>%
                          filter(modulo == 0, novel_type == "novel_donor") %>%
                          nrow()) * 100 / num_novel_donor
      
      modulo0_aceptor <- (df_novel_mane %>%
                            filter(modulo == 0, novel_type == "novel_acceptor") %>%
                            nrow()) * 100 / num_novel_acceptor
      
      ## Number of modulo 1 or 2 unique introns
      modulo1_donor <- (df_novel_mane %>%
                          filter(modulo == 1, novel_type == "novel_donor") %>%
                          nrow()) * 100 / num_novel_donor
      
      modulo1_aceptor <- (df_novel_mane %>%
                            filter(modulo == 1, novel_type == "novel_acceptor") %>%
                            nrow()) * 100 / num_novel_acceptor
      
      modulo2_donor <- (df_novel_mane %>%
                          filter(modulo == 2, novel_type == "novel_donor") %>%
                          nrow()) * 100 / num_novel_donor
      
      modulo2_aceptor <- (df_novel_mane %>%
                            filter(modulo == 2, novel_type == "novel_acceptor") %>%
                            nrow()) * 100 / num_novel_acceptor
      
      df <- data.frame(tissue = tissues_tidy[i],
                       modulo = modulo0_donor,
                       novel_type = "donor",
                       type = "modulo0")
      
      df <- rbind(df,
                  data.frame(tissue = tissues_tidy[i],
                             modulo = modulo1_donor,
                             novel_type = "donor",
                             type = "modulo1"))
      
      df <- rbind(df,
                  data.frame(tissue = tissues_tidy[i],
                             modulo = modulo2_donor,
                             novel_type = "donor",
                             type = "modulo2"))
      
      df <- rbind(df,
                  data.frame(tissue = tissues_tidy[i],
                             modulo = modulo0_aceptor,
                             novel_type = "acceptor",
                             type = "modulo0"))
      
      df <- rbind(df,
                  data.frame(tissue = tissues_tidy[i],
                             modulo = modulo1_aceptor,
                             novel_type = "acceptor",
                             type = "modulo1"))
      
      df <- rbind(df,
                  data.frame(tissue = tissues_tidy[i],
                             modulo = modulo2_aceptor,
                             novel_type = "acceptor",
                             type = "modulo2"))
      return(df)
      
    })
  
  #if (type == -1) {
  #  saveRDS(df_modulo_prop, file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/modulo_intron_all_tissues.rds")
  #} else {
  #  saveRDS(df_modulo_prop, file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/modulo_exon_all_tissues.rds")
  #}
  
  saveRDS(df_modulo_prop, file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/modulo_all_tissues.rds")
  
  #})
}

## MSR ----------------------------------------

summarise_MSR_tissues <- function() {
  
  if (!exists(gtex_tissues)) {
    gtex_tissues <-  readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")
  }
  
  all_data <- map_df(gtex_tissues, function(tissue) {
    
    print(tissue)
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/",
                          tissue, "/v104/")
    
    
    df <- readRDS(file = paste0(folder_name, "/", tissue, "_db_introns.rds")) %>% as_tibble()
    df %>% head()
    df %>% nrow()
    
    
    df_tidy_donor <- df %>%
      pull(ref_missplicing_ratio_tissue_ND) %>%
      summary() %>% 
      as.array %>%
      as.data.frame() %>%
      mutate(tissue = tissue) %>%
      mutate("MSR_type" = "MSR_Donor")
    
    
    ##############
    
    df_tidy_acceptor <- df %>%
      pull(ref_missplicing_ratio_tissue_NA) %>%
      summary() %>% 
      as.array() %>%
      as.data.frame() %>%
      mutate(tissue = tissue) %>%
      mutate("MSR_type" = "MSR_Acceptor")
    
    
    return(rbind(df_tidy_donor, df_tidy_acceptor))
  })
  
  saveRDS(object = all_data,
          file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/summary_MSR_tissues.rds")
}

## LINEAR MODELS ------------------------------------

summarise_lm_tissues <- function() {
  
  if (!exists(gtex_tissues)) {
    gtex_tissues <-  readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")
  }
  
  
  all_data <- map_df(gtex_tissues, function(tissue) {
    
    # tissue <- gtex_tissues[11]
    print(tissue)
    
    ## Load the IDB
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", 
                          tissue, "/v104/", tissue, "_db_introns.rds")
    
    idb <- readRDS(file = folder_name) 
    
    idb <- idb %>%
      as.data.frame() %>%
      dplyr::distinct(ref_junID, .keep_all = T) %>%
      dplyr::rename(intron_length = width,
                    intron_5ss_score = ref_ss5score,
                    intron_3ss_score = ref_ss3score,
                    gene_length = gene_width,
                    gene_tpm = tpm,
                    gene_num_transcripts = transcript_number,
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
    
    
    
    
    MSR_Donor <- data.frame(tissue = tissue,
                            "MSR_Type" = "Donor",
                            intron_length = model$coefficients["intron_length"] %>% unname(),
                            intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                            intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                            gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
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
    
    MSR_Acceptor <- data.frame(tissue = tissue,
                               "MSR_Type" = "Acceptor",
                               intron_length = model$coefficients["intron_length"] %>% unname(),
                               intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                               intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                               gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
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
  
  saveRDS(object = all_data,
          file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/summary_tissues_lm.rds")
}

get_estimate_variance <- function(clusters = gtex_tissues[7:17],
                                  save_results = T) {
  
  
  MSR_Donor <- data.frame(tissue = as.character(),
                          intron_length = as.double(),
                          intron_5ss_score = as.double(),
                          intron_3ss_score = as.double(),
                          gene_tpm = as.double(),
                          gene_length = as.double(),
                          # clinvarTRUE = as.double(),
                          protein_coding = as.double(),
                          gene_num_transcripts = as.double(),
                          CDTS_5ss = as.double(),
                          CDTS_3ss = as.double(),
                          mean_phastCons20way_5ss = as.double(),
                          mean_phastCons20way_3ss = as.double(),
                          intron_5ss_scoreintron_3ss_score = as.double())
  
  
  MSR_Acceptor <- data.frame(tissue = as.character(),
                             intron_length = as.double(),
                             intron_5ss_score = as.double(),
                             intron_3ss_score = as.double(),
                             gene_tpm = as.double(),
                             gene_length = as.double(),
                             # clinvarTRUE = as.double(),
                             protein_coding = as.double(),
                             gene_num_transcripts = as.double(),
                             CDTS_5ss = as.double(),
                             CDTS_3ss = as.double(),
                             mean_phastCons20way_5ss = as.double(),
                             mean_phastCons20way_3ss = as.double(),
                             intron_5ss_scoreintron_3ss_score = as.double())
  
  # MSR_Donor_list <- list()
  # MSR_Acceptor_list <- list()
  
  
  #or (protein_percentage in c(0, 100)) {
  
  #print(paste0("Protein-percentage: ", protein_percentage))
  
  # protein_percentage = 0
  for (tissue in clusters) {
    # tissue <- clusters[1]
    print(tissue)
    
    ## Load the IDB
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", 
                          tissue, "/v104/", tissue, "_db_introns.rds")
    
    idb <- readRDS(file = folder_name) #%>%
    #filter(protein_coding == protein_percentage)
    
    common_introns <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/all_introns_brain_tidy.rds")
    
    idb <- idb %>%
      as.data.frame() %>%
      dplyr::distinct(ref_junID, .keep_all = T) %>%
      filter(u2_intron == T) %>%
      dplyr::rename(intron_length = width,
                    intron_5ss_score = ref_ss5score,
                    intron_3ss_score = ref_ss3score,
                    gene_length = gene_width,
                    gene_tpm = tpm,
                    gene_num_transcripts = transcript_number,
                    CDTS_5ss = CDTS_5ss_mean,
                    CDTS_3ss = CDTS_3ss_mean,
                    mean_phastCons20way_5ss = phastCons20way_5ss_mean,
                    mean_phastCons20way_3ss = phastCons20way_3ss_mean)  %>%
      filter(gene_tpm > 0)  %>% 
      filter(intron_length < gene_length) %>%
      filter(ref_junID %in% common_introns$ref_junID)
    
    # idb %>% nrow() %>% print()
    
    
    ## Donor
    model <- lm(ref_missplicing_ratio_tissue_ND ~ 
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
    #MSR_Donor_list[[tissue]] <- model
    
    
    ind_sign <- which(((model %>% summary())$coefficients  %>% as.data.frame())[,4] < 0.05)
    model$coefficients[-ind_sign] <- 0
    
    MSR_Donor <- rbind(MSR_Donor,
                       data.frame(tissue = tissue,
                                  intron_length = model$coefficients["intron_length"] %>% unname(),
                                  intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                                  intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                                  gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                                  gene_length = model$coefficients["gene_length"] %>% unname(),
                                  #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),
                                  protein_coding = model$coefficients["protein_coding"] %>% unname(),
                                  gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                                  intron_5ss_scoreintron_3ss_score = model$coefficients["intron_5ss_score:intron_3ss_score"] %>% unname(),
                                  CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                                  CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                                  mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                                  mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname()))
    
    
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
    
    MSR_Acceptor <- rbind(MSR_Acceptor,
                          data.frame(tissue = tissue,
                                     intron_length = model$coefficients["intron_length"] %>% unname(),
                                     intron_5ss_score = model$coefficients["intron_5ss_score"] %>% unname(),
                                     intron_3ss_score = model$coefficients["intron_3ss_score"] %>% unname(),
                                     gene_tpm = model$coefficients["gene_tpm"] %>% unname(),
                                     gene_length = model$coefficients["gene_length"] %>% unname(),
                                     #clinvarTRUE = model$coefficients["clinvarTRUE"] %>% unname(),
                                     protein_coding = model$coefficients["protein_coding"] %>% unname(),
                                     gene_num_transcripts = model$coefficients["gene_num_transcripts"] %>% unname(),
                                     intron_5ss_scoreintron_3ss_score = model$coefficients["intron_5ss_score:intron_3ss_score"] %>% unname(),
                                     CDTS_5ss = model$coefficients["CDTS_5ss"] %>% unname(),
                                     CDTS_3ss  = model$coefficients["CDTS_3ss"] %>% unname(),
                                     mean_phastCons20way_5ss  = model$coefficients["mean_phastCons20way_5ss"] %>% unname(),
                                     mean_phastCons20way_3ss = model$coefficients["mean_phastCons20way_3ss"] %>% unname()))
  }
  
  if (save_results) {
    saveRDS(object = rbind(MSR_Donor %>% mutate(type = "MSR_Donor"), 
                           MSR_Acceptor %>% mutate(type = "MSR_Acceptor")),
            file = paste0("/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_brain_tissues.rds"))
    
    
    #############################
    ## ADD DATA TO AN EXCEL FILE
    #############################
    
    library(xlsx)
    ## Save the .xlsx object
    write.xlsx(MSR_Donor, 
               file = paste0("/home/sruiz/PROJECTS/splicing-project/results/paper/tissues_variance_features.xlsx"), 
               sheetName="MSR_Donor")
    write.xlsx(MSR_Acceptor, 
               file = paste0("/home/sruiz/PROJECTS/splicing-project/results/paper/tissues_variance_features.xlsx"),   
               sheetName="MSR_Acceptor",
               append = T)
  }
  #}
  
  
}


correlation_core_spliceosome_unique_junc <- function(clusters = gtex_tissues) {
  
  ## 0. Load the mean number of reads
  
  # df_mean_counts <- readRDS("/home/sruiz/PROJECTS/splicing-project/results/paper/figure1_data.rds")
  # df_mean_counts %>% head()
  
  ## 1. Load the unique number of junctions across samples of each tissue and normalised by the number of samples
  
  df_unique_junc <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/paper/figure3_data.rds")
  df_unique_junc %>% head()
  
  df_unique_junc <- df_unique_junc %>%
    rowwise() %>%
    mutate(annotated_junc = annotated_junc / samples,
           donor_junc = donor_junc / samples,
           acceptor_junc = acceptor_junc / samples)
  # df_unique_junc <- df_unique_junc %>%
  #   mutate(annotated_prop = annotated_junc/(annotated_junc + donor_junc + acceptor_junc),
  #          donor_prop = donor_junc/(annotated_junc + donor_junc + acceptor_junc),
  #          acceptor_prop = acceptor_junc/(annotated_junc + donor_junc + acceptor_junc))
  
  df_unique_junc %>% head()
  
  ## 2. Load the TPM values of the core spliceosome genes
  
  # SF3B1 (ENSG00000115524): recognition of branch point
  # U2AF1 (ENSG00000160201): recognition of AG
  # SRSF2 (ENSG00000161547): recognition of ESE
  # ZRSR2 (ENSG00000169249): splicing of u12-introns
  
  tpm <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_gene_median_tpm_tidy.rds") 
  
  tpm_SF3B1 <- tpm %>%
    filter(gene_id == "ENSG00000177628") %>%
    gather(key = tissue, value = gene_id) %>% 
    rename(SF3B1 = gene_id)
  tpm_SF3B1 %>% head()
  
  tpm_U2AF1 <- tpm %>%
    filter(gene_id == "ENSG00000160201") %>%
    gather(key = tissue, value = gene_id) %>% 
    rename(U2AF1 = gene_id)
  tpm_U2AF1 %>% head()
  
  tpm_SRSF2 <- tpm %>%
    filter(gene_id == "ENSG00000161547") %>%
    gather(key = tissue, value = gene_id) %>% 
    rename(SRSF2 = gene_id)
  tpm_SRSF2 %>% head()
  
  tpm_core <- merge(tpm_SF3B1, tpm_U2AF1, by = "tissue")
  tpm_core <- merge(tpm_core, tpm_SRSF2, by = "tissue")
  
  
  # ## 3. Load the TPM values of the core NMD genes
  # NMD_genes <- c("ENSG00000005007", "ENSG00000151461", "ENSG00000169062", "ENSG00000125351")
  # 
  # # UPF1 (ENSG00000005007)
  # # UPF2 (ENSG00000151461)
  # # UPF3A (ENSG00000169062)
  # # UPF3B (ENSG00000125351)
  # 
  # tpm_UPF1 <- tpm %>%
  #   filter(gene_id == "ENSG00000005007") %>%
  #   gather(key = tissue, value = gene_id) %>% 
  #   rename(UPF1 = gene_id)
  # tpm_UPF1 %>% head()
  # 
  # tpm_UPF2 <- tpm %>%
  #   filter(gene_id == "ENSG00000151461") %>%
  #   gather(key = tissue, value = gene_id) %>% 
  #   rename(UPF2 = gene_id)
  # tpm_UPF2 %>% head()
  # 
  # tpm_UPF3A <- tpm %>%
  #   filter(gene_id == "ENSG00000169062") %>%
  #   gather(key = tissue, value = gene_id) %>% 
  #   rename(UPF3A = gene_id)
  # tpm_UPF3A %>% head()
  # 
  # tpm_UPF3B <- tpm %>%
  #   filter(gene_id == "ENSG00000125351") %>%
  #   gather(key = tissue, value = gene_id) %>% 
  #   rename(UPF3B = gene_id)
  # tpm_UPF3B %>% head()
  # 
  # 
  # tpm_core <- merge(tpm_core, tpm_UPF1, by = "tissue")
  # tpm_core <- merge(tpm_core, tpm_UPF2, by = "tissue")
  # tpm_core <- merge(tpm_core, tpm_UPF3A, by = "tissue")
  # tpm_core <- merge(tpm_core, tpm_UPF3B, by = "tissue")
  # 
  
  tpm_core %>% head()
  
  ## 3. Merge TPM and unique junctions datasets
  
  df_merged <- merge(df_unique_junc,
                     tpm_core, 
                     by = "tissue")
  
  df_merged %>% head()
  
  # df_merged <- merge(df_merged,
  #                    df_mean_counts, 
  #                    by = "tissue")
  # 
  # df_merged %>% head()
  
  ## 4. Run correlation tests
  
  # ## Annotated - Core Spliceosome
  cor.test(x = df_merged$annotated_junc,
           y = df_merged$SF3B1)
  cor.test(x = df_merged$annotated_junc,
           y = df_merged$U2AF1)
  cor.test(x = df_merged$annotated_junc,
           y = df_merged$SRSF2)
  
  # 
  # ## Annotated - Core NMD
  # cor.test(x = df_merged$annotated_junc,
  #          y = df_merged$UPF1)
  # cor.test(x = df_merged$annotated_junc,
  #          y = df_merged$UPF2)
  # cor.test(x = df_merged$annotated_junc,
  #          y = df_merged$UPF3A)
  # cor.test(x = df_merged$annotated_junc,
  #          y = df_merged$UPF3B)
  
  ## Acceptor - Core Spliceosome
  cor.test(x = (df_merged$acceptor_junc + df_merged$donor_junc),
           y = df_merged$SF3B1)
  cor.test(x = (df_merged$acceptor_junc/(df_merged$acceptor_junc+df_merged$donor_junc)),
           y = df_merged$U2AF1)
  cor.test(x = (df_merged$acceptor_junc/(df_merged$acceptor_junc+df_merged$donor_junc)),
           y = df_merged$SRSF2)
  
  
  # ## Acceptor - Core NMD
  # cor.test(x = df_merged$acceptor_junc,
  #          y = df_merged$UPF1)
  # cor.test(x = df_merged$acceptor_junc,
  #          y = df_merged$UPF2)
  # cor.test(x = df_merged$acceptor_junc,
  #          y = df_merged$UPF3A)
  # cor.test(x = df_merged$acceptor_junc,
  #          y = df_merged$UPF3B)
  # 
  ## Donor - Core Spliceosome
  cor.test(x = df_merged$donor_junc,
           y = df_merged$SF3B1)
  cor.test(x = df_merged$donor_junc,
           y = df_merged$U2AF1)
  cor.test(x = df_merged$donor_junc,
           y = df_merged$SRSF2)
  # 
  # 
  # ## Donor - Core NMD
  # cor.test(x = df_merged$donor_junc,
  #          y = df_merged$UPF1)
  # cor.test(x = df_merged$donor_junc,
  #          y = df_merged$UPF2)
  # cor.test(x = df_merged$donor_junc,
  #          y = df_merged$UPF3A)
  # cor.test(x = df_merged$donor_junc,
  #          y = df_merged$UPF3B)
  
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

get_common_junctions <- function(clusters = gtex_tissues,
                                 folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/",
                                 folder_save_name,
                                 only_brain = F) {
  
  if (only_brain) {
    clusters <- gtex_tissues[c(7:17)]
  } else {
    clusters <- gtex_tissues
  }
  
  print(paste0(Sys.time(), " - getting unique and common junctions across clusters..."))
  
  
  ## Getting all introns that are common across all GTEx tissues -------------------------------------------
  
  all_introns <- list()
  
  for(cluster in clusters) { 
    
    folder_name <- paste0(folder_root, "/", cluster)
    
    df <- readRDS(file = paste0(folder_name, "/v104/", cluster, "_db_introns.rds"))
    
    all_introns[[cluster]] <- df$ref_junID %>% 
      unique()
    
    print(paste0(Sys.time(), " - intron IDs collected from '", cluster, "'"))
    
  }
  
  common_introns <- data.frame(ref_junID = Reduce(intersect,  all_introns))
  common_introns %>% head()
  common_introns %>% nrow()
  
  
  ## QC - all junctions commonly found across all tissues should be all also found in 
  ## the original dataset of junctions of each tissue
  for (cluster in names(all_introns)) {
    
    if (intersect(all_introns[[cluster]], common_introns$ref_junID) %>% 
        unique %>% 
        length() != common_introns$ref_junID %>% length()) {
      
      print(paste0("Cluster ", cluster, " presents a mismatch of common junctions!"))
      break;
      
    } else {
      
      print(paste0(cluster, " -> OK"))
      
    }
    
  }
  
  ## Obtaining all info from the common introns --------------------------------------------------------------
  
  i <- 1
  all_introns_details <- NULL
  print(paste0(Sys.time(), " - getting mis-splicing ratio data from all introns ..."))
  
  for(cluster in clusters) { 
    # cluster <- clusters[5]
    
    folder_name <- paste0(folder_root, "/", cluster)
    df <- readRDS(file = paste0(folder_name, "/v104/", cluster, "_db_introns.rds")) %>%
      as_tibble() %>%
      dplyr::filter(ref_junID %in% (common_introns$ref_junID %>% unique()))
    
    if (i == 1) {
      if ((df %>% nrow() == common_introns %>% nrow()) &&
          (identical(sort(df$ref_junID), sort(common_introns$ref_junID)))) {
        
        all_introns_details <- df
      }
      
    } else {
      all_introns_details <- rbind(all_introns_details, df)
      
    }
    
    print(paste0(Sys.time(), " - mis-splicing info collected from '", cluster, "'"))
    i <- i + 1
    
  }
  
  ## Save results -------------------------------------------------
  
  all_introns_details %>% head()
  all_introns_details %>% nrow()
  
  if (only_brain) {
    file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_brain_raw.rds"
  } else{
    file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_all_tissues_raw.rds"
  }
  
  saveRDS(object = all_introns_details, file = file_name)
  
  
  ###############################################
  ## ASSIGN JUNCTION CATEGORY TO COMMON JUNCTIONS
  ###############################################
  
  if (!exists("all_introns_details")) {
    if (only_brain) {
      file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_brain_raw.rds"
    } else{
      file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_all_tissues_raw.rds"
    }
    all_introns_details <- readRDS(file = file_name)
  }
  
  
  ## Tidy up the data and define the intron type -------------------------------------------------
  
  identical(all_introns_details$ref_junID %>% unique(),
            common_introns$ref_junID)
  
  all_introns_tidy <- all_introns_details %>% 
    group_by(ref_junID) %>%
    mutate(ref_missplicingratio_ND_tissues = ref_missplicing_ratio_tissue_ND %>% mean(na.rm = T)) %>%
    mutate(ref_missplicingratio_NA_tissues = ref_missplicing_ratio_tissue_NA %>% mean(na.rm = T)) %>%
    mutate(final_type = case_when(
      (ref_type %>% unique() %>% length() == 1) ~ ref_type,
      (!(any(ref_type == "never")) & ref_type %>% unique() %>% length() > 1) ~ "both",
      (any(ref_type == "never") & ref_type %>% unique() %>% length > 2) ~ "both",
      (ref_missplicingratio_ND_tissues == 0 & ref_missplicingratio_NA_tissues == 0) ~ "never",
      (any(ref_type == "acceptor") & any(ref_type == "never") & ref_type %>% unique() %>% length() == 2) ~ "acceptor",
      (any(ref_type == "donor") && any(ref_type == "never") & ref_type %>% unique() %>% length == 2) ~ "donor",
      (any(ref_type == "both") && any(ref_type == "never") & ref_type %>% unique() %>% length == 2) ~ "both",
      TRUE ~ "NA"
    )) %>%
    ungroup() %>%
    distinct(ref_junID, .keep_all = T) 
  
  
  ## Save results -------------------------------------------------
  
  all_introns_tidy %>% head()
  all_introns_tidy %>% nrow()
  
  if (only_brain) {
    file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_brain_tidy.rds"
  } else{
    file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_all_tissues_tidy.rds"
  }
  
  all_introns_tidy <- all_introns_tidy %>%
    filter(u2_intron == T)
  all_introns_tidy %>% nrow()
  #df_introns<-all_introns_tidy
  #all_introns_tidy<- df_merged
  saveRDS(object = all_introns_tidy, file = file_name)
  
}


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
                              #PC = NULL,
                              gene = NULL,
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
  
  
  if (!is.null(gene)) {
    df_aux <- merge(x = df_novel, 
                    y = all_split_reads_details_97, 
                    by.x = "ref_junID", 
                    by.y = "junID")
    df_aux <- df_aux %>%
      filter(gene_name_start == gene, gene_name_end == gene) 
    
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
  
  # y_axes <- c(0, (df_novel_tidy %>% 
  #                   filter(novel_type == "novel_acceptor", distance == (df_novel_tidy$distance %>% get_mode()))  %>% 
  #                   nrow()))
  #   
  
  
  
  #################################
  ## GENERATE PLOT ################
  #################################
  
  # print("Novel donor:")
  # # print(df_novel_tidy %>%
  # #                        filter(novel_type == "novel_donor") %>%
  # #                        pull(distance) %>% summary() %>% as.table())
  # print(paste0("Mode (+): ", df_novel_tidy %>%
  #                filter(novel_type == "novel_donor") %>%
  #                filter(distance > 0) %>%
  #                pull(distance) %>% get_mode(), " - "))
  # print(paste0("Mode (-): ", df_novel_tidy %>%
  #                filter(novel_type == "novel_donor") %>%
  #                filter(distance < 0) %>%
  #                pull(distance) %>% get_mode(), " - "))
  # 
  # print("Novel acceptor:")
  # print(paste0("Mode (+): ", df_novel_tidy %>%
  #                filter(novel_type == "novel_acceptor") %>%
  #                filter(distance > 0) %>%
  #                pull(distance) %>% get_mode(), " - "))
  # print(paste0("Mode (-): ", df_novel_tidy %>%
  #                filter(novel_type == "novel_acceptor") %>%
  #                filter(distance < 0) %>%
  #                pull(distance) %>% get_mode(), " - "))
  
  
  plot_PC <- ggplot(data = df_novel_tidy %>% filter(type_PC == "protein coding (PC)")) + 
    geom_histogram(aes(x = distance, fill = novel_type),
                   bins = limit_bp * 2,
                   binwidth = 1,
                   position = "stack"
    ) +
    facet_col(vars(novel_type)) +
    ggtitle(paste0(title, "\nprotein coding (PC)")) +
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
    facet_col(vars(novel_type)) +
    ggtitle(paste0(title, "\nnon PC")) +
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
    facet_zoom(xlim = c(0,0.15)) +
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
                                   #PC = NULL,
                                   #intron_type = NULL,
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
    facet_zoom(xlim = c(0,0.15)) +
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
    facet_zoom(xlim = c(0,0.15)) +
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
    mutate(ref_mean_counts = ref_sum_counts/ref_n_individuals) %>%
    dplyr::select(ref_junID, 
                  ref_ss5score, ref_ss3score, 
                  ref_mean_counts, ref_type)
  df_introns %>% head()
  df_introns %>% nrow()

  
  
  ## Load the novel junction database ---------------------------------------------
  
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  df_novel <- df_novel %>%
    mutate(novel_mean_counts = novel_sum_counts/novel_n_individuals) %>%
    dplyr::select(ref_junID, 
                  novel_junID, 
                  novel_ss5score, novel_ss3score, 
                  novel_mean_counts, novel_type, gene_id)
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
    dplyr::select(diff_ss5score, diff_ss3score, novel_mean_counts, ref_mean_counts, novel_type, gene_id) %>%
    tidyr::gather(key = "type", value = "mean", -novel_mean_counts, -novel_type, -ref_mean_counts, -gene_id) %>%
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
                  ref_mean_counts, protein_coding, ref_type, protein_coding)
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
                  novel_mean_counts, novel_type, gene_id)
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
    dplyr::select(diff_ss5score, diff_ss3score, novel_mean_counts, ref_mean_counts, novel_type, gene_id,protein_coding,type_PC) %>%
    tidyr::gather(key = "type", value = "mean", -novel_mean_counts, -novel_type, -ref_mean_counts, -gene_id, -protein_coding, -type_PC ) %>%
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
    filter(protein_coding %in% c(0,100)) %>%
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
  #                 gene_tpm = tpm,
  #                 gene_num_transcripts = transcript_number,
  #                 CDTS_5ss = CDTS_5ss_mean,
  #                 CDTS_3ss = CDTS_3ss_mean,
  #                 mean_phastCons20way_5ss = phastCons20way_5ss_mean,
  #                 mean_phastCons20way_3ss = phastCons20way_3ss_mean)
  idb <- idb %>%
    as.data.frame() %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    dplyr::rename(intron_length = width,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score)
  
  #########################
  ## LINEAR MODELS
  #########################
  
  # file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/common_introns_brain_raw.rds"
  # common_introns_brain <- readRDS(file = file_name)
  # idb <- idb %>%
  #   filter(ref_junID %in% common_introns_brain$ref_junID)
  fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ intron_length +
                    intron_length + 
                    intron_5ss_score + 
                    intron_3ss_score, 
                  data = idb) #%>%)
  # fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ intron_length +
  #                   intron_length + 
  #                   intron_5ss_score + 
  #                   intron_3ss_score + 
  #                   gene_tpm + 
  #                   gene_length + 
  #                   CDTS_5ss + 
  #                   CDTS_3ss + 
  #                   mean_phastCons20way_5ss +
  #                   mean_phastCons20way_3ss +
  #                   #clinvar +
  #                   protein_coding +
  #                   gene_num_transcripts, 
  #                 data = idb) #%>%)
  fit_donor %>% summary()
  
  fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ 
                       intron_length + 
                       intron_5ss_score + 
                       intron_3ss_score, 
                     data = idb)
  fit_acceptor %>% summary()
  
  model_names <- c("MSR_Donor", "MSR_Acceptor")
  
  
  coef_names <- c("Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  "Inter. 5'ss & 3'ss MES score" = "intron_5ss_score:intron_3ss_score", 
                  "Intron Type u2" = "u2_intronTRUE",
                  "Intron ClinVar mutation" = "clinvarTRUE",
                  "Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts", 
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

plot_estimate_variance <- function(brain_tissues = F,
                                   save_results = F) {
  
  if (brain_tissues) {
    df_estimate <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/paper/variance_estimate_brain_tissues.rds"))
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
                  "Inter. 5'ss & 3'ss MES score" = "intron_5ss_scoreintron_3ss_score", 
                  #"Intron Type u2" = "u2_intronTRUE",
                  #"Intron ClinVar mutation" = "clinvarTRUE",
                  "Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  #"CDTS 5'ss" = "CDTS_5ss",
                  #"CDTS 3'ss" = "CDTS_3ss",
                  #"Conservation 5'ss" = "mean_phastCons20way_5ss",
                  #"Conservation 3'ss" = "mean_phastCons20way_3ss", 
                  "Protein coding" = "protein_coding") %>%
    gather(tissue, feature, -type, -tissue) %>%
    mutate(type = "MSR_Donor")
  
  MSR_Acceptor_tidy <- MSR_Acceptor %>%
    dplyr::rename("Intron Length" = "intron_length",
                  "Intron 5'ss MES score" = "intron_5ss_score",
                  "Intron 3'ss MES score" = "intron_3ss_score", 
                  "Inter. 5'ss & 3'ss MES score" = "intron_5ss_scoreintron_3ss_score", 
                  # "Intron Type u2" = "u2_intronTRUE",
                  # "Intron ClinVar mutation" = "clinvarTRUE",
                  "Gene Length" = "gene_length",
                  "Gene TPM" = "gene_tpm",
                  "Gene num. transcripts" = "gene_num_transcripts",
                  #"CDTS 5'ss" = "CDTS_5ss",
                  #"CDTS 3'ss" = "CDTS_3ss",
                  #"Conservation 5'ss" = "mean_phastCons20way_5ss",
                  #"Conservation 3'ss" = "mean_phastCons20way_3ss", 
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

plot_missplicing_ratio_vs_reads <- function(cluster,
                                            folder_name) {
  
  df <- readRDS(file = paste0(folder_name, "/", cluster, "_missplicingratio_tidy_v2.rds"))
  df %>% head()
  df %>% nrow()
  
  print(paste0(Sys.time(), " - preparing the data ..."))
  
  folder_img <- paste0(folder_name, "/images/")
  dir.create(file.path(folder_img), showWarnings = T)
  
  
  ggplot(df, aes(x = missplicingratio_ND_tissue, y = ref_counts)) +
    geom_point(size = 1) + 
    facet_zoom(ylim = c(1,1000)) +
    scale_color_viridis_d() +
    ggtitle(paste0("Mis-splicing ratio at donor positions vs. reads of the reference intron.\nFrontalCortex tissue.")) +
    labs(x = "mis-splicing ratio at donor positions",
         y = "intron reads") +
    theme_light(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  ggplot(df, aes(x = missplicingratio_NA_tissue, y = ref_counts)) +
    geom_point(size = 1) + 
    facet_zoom(ylim = c(1,1000)) +
    scale_color_viridis_d() +
    ggtitle(paste0("Mis-splicing ratio at acceptor positions vs. reads of the reference intron.\nFrontalCortex tissue.")) +
    labs(x = "mis-splicing ratio at acceptor positions",
         y = "intron reads") +
    theme_light(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

constituive_exons <- function() {
  
  # constitutive_exons <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/constitutive_exons.rds")
  constitutive_exons <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/constitutive_exons_ranges.rds") %>%
    addchr()
  
  
  
  
  ## Load never-misspliced junctions
  ## Convert them into GRanges object
  never_misspliced_junc <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/never_misspliced_junctions.rds") %>%
    distinct(ref_junID, .keep_all = T) %>%
    GRanges() %>%
    addchr()
  
  
  ## Load sometimes mis-spliced junctions
  ## Convert them into GRanges object
  always_misspliced_junc <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/always_misspliced_junctions.rds") %>%
    distinct(ref_junID, .keep_all = T) %>%
    GRanges() %>%
    addchr()
  
  
  
  ## Constitutive exons and forward strand -------------------------------------------------------------------------------------------------------------
  
  # constitutive_exons_forward <- constitutive_exons[constitutive_exons@strand == "+"]
  # never_misspliced_forward <- never_misspliced_junc[never_misspliced_junc@strand == "+"]
  
  ## Find upstream constitutive exons in the forward strand
  overlaps_upstream <- GenomicRanges::findOverlaps(GenomicRanges::GRanges(seqnames = constitutive_exons %>% seqnames(),
                                                                          ranges = IRanges(start = constitutive_exons %>% end(), 
                                                                                           end = constitutive_exons %>% end())),
                                                   GenomicRanges::GRanges(seqnames = never_misspliced_junc %>% seqnames(),
                                                                          ranges = IRanges(start = never_misspliced_junc %>% start() - 1, 
                                                                                           end = never_misspliced_junc %>% start() - 1))
  )
  
  constitutive_upstream <- constitutive_exons[queryHits(overlaps_upstream),]
  never_upstream <- never_misspliced_junc[subjectHits(overlaps_upstream),]
  
  df_never_upstream <- data.frame(constitutive_exon = constitutive_upstream$exon_id,
                                  never_misspliced = never_upstream$ref_junID)
  
  ## Find downstream constitutive exons in the forward strand
  overlaps_downstream <- GenomicRanges::findOverlaps(GenomicRanges::GRanges(seqnames = constitutive_exons %>% seqnames(),
                                                                            ranges = IRanges(start = constitutive_exons %>% start(), 
                                                                                             end = constitutive_exons %>% start())),
                                                     GenomicRanges::GRanges(seqnames = never_misspliced_junc %>% seqnames(),
                                                                            ranges = IRanges(start = never_misspliced_junc %>% end() + 1, 
                                                                                             end = never_misspliced_junc %>% end() + 1))
  )
  
  constitutive_downstream <- constitutive_exons[queryHits(overlaps_downstream),]
  never_downstream <- never_misspliced_junc[subjectHits(overlaps_downstream),]
  
  df_never_downstream <- data.frame(constitutive_exon = constitutive_downstream$exon_id,
                                    never_misspliced = never_downstream$ref_junID)
  
  
  ## Join results
  df <- rbind(df_never_upstream,
              df_never_downstream)
  
  df %>%
    filter(never_misspliced == "11607971")
  
  df <- df %>%
    group_by(never_misspliced) %>% 
    dplyr::count() %>%
    ungroup()
  
  df$n %>% unique
  ## Constitutive exons and reverse strand -------------------------------------------------------------------------------------------------------------
  
  # constitutive_exons_reverse <- constitutive_exons[constitutive_exons@strand == "-"]
  # never_misspliced_reverse <- never_misspliced_junc[never_misspliced_junc@strand == "-"]
  # 
  # ## Find downstream constitutive exons in the reverse strand
  # overlaps3 <- GenomicRanges::findOverlaps(GenomicRanges::GRanges(seqnames = constitutive_exons_reverse %>% seqnames(),
  #                                                                ranges = IRanges(start = constitutive_exons_reverse %>% end(), 
  #                                                                                 end = constitutive_exons_reverse %>% end())),
  #                                         GenomicRanges::GRanges(seqnames = never_misspliced_reverse %>% seqnames(),
  #                                                                ranges = IRanges(start = never_misspliced_reverse %>% start() - 1, 
  #                                                                                 end = never_misspliced_reverse %>% start() - 1)),
  #                                         type = "equal"
  # )
  # ## Find upstream constitutive exons in the reverse strand
  # overlaps4 <- GenomicRanges::findOverlaps(GenomicRanges::GRanges(seqnames = constitutive_exons_reverse %>% seqnames(),
  #                                                                 ranges = IRanges(start = constitutive_exons_reverse %>% start(), 
  #                                                                                  end = constitutive_exons_reverse %>% start())),
  #                                          GenomicRanges::GRanges(seqnames = never_misspliced_reverse %>% seqnames(),
  #                                                                 ranges = IRanges(start = never_misspliced_reverse %>% end() + 1, 
  #                                                                                  end = never_misspliced_reverse %>% end() + 1)),
  #                                          type = "equal"
  # )
  
  
  
  upstream_const_exon_reverse <- constitutive_exons_reverse[queryHits(overlaps3),]
  never_reverse <- never_misspliced_reverse[subjectHits(overlaps3),]
  downstream_const_exon_reverse <- constitutive_exons_reverse[queryHits(overlaps4),]
  never_reverse2 <- never_misspliced_reverse[subjectHits(overlaps4),]
  
  df_reverse_upstream <- data.frame(constitutive_exon = upstream_const_exon_reverse$exon_id,
                                    never_misspliced = never_reverse$ref_junID)
  
  df_reverse_downstream <- data.frame(constitutive_exon = downstream_const_exon_reverse$exon_id,
                                      never_misspliced = never_reverse2$ref_junID)
  
  df2 <- rbind(df_reverse_upstream,
               df_reverse_downstream)
  
  
  
  df3 <- rbind(df, df2)
  
  
  
  
  df3 %>% 
    group_by(never_misspliced)  %>% 
    summarise(n = n())
  
  
  nearest_constitutive <- follow(x = GenomicRanges::GRanges(seqnames = never_misspliced_junc %>% seqnames(),
                                                            ranges = never_misspliced_junc %>% ranges,
                                                            strand = never_misspliced_junc %>% strand()),
                                 subject = GenomicRanges::GRanges(seqnames = constitutive_exons %>% seqnames,
                                                                  ranges = constitutive_exons %>% ranges,
                                                                  strand = constitutive_exons %>% strand),
                                 type = "start")
  never_misspliced_junc
  constitutive_exons[nearest_constitutive,]
  
  
  
  
  always_misspliced_junc$gene_id %>% unique
  intersect(constitutive_exons$gene_id, always_misspliced_junc$gene_id)
  
  
  
  # never_misspliced_junc <- never_misspliced_junc[never_misspliced_junc@strand == "+"]
  # always_misspliced_junc <- always_misspliced_junc[always_misspliced_junc@strand == "+"]
  
}


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

correlate_TPM_NMD_genes <- function(clusters,
                                    gtex_tissues,
                                    folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/") {
  
  ## Load TPM file and format it -------------------------------------------------------------------------------------
  
  
  library(GSRI)
  
  gtex_tpm <- GSRI::readGct(file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct") %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_id") 
  
  
  gtex_tpm <- gtex_tpm  %>% 
    separate(gene_id, into = c("gene_id", "post"), sep = "\\.") %>%
    #select(gene_id, tpm = Brain...Frontal.Cortex..BA9.) %>%
    distinct(gene_id, .keep_all = T) 
  
  gtex_tpm %>% head()
  gtex_tpm %>% nrow()
  
  
  gtex_tpm <- gtex_tpm %>% 
    select(-post, 
           -Bladder,
           -Cervix...Ectocervix,
           -Cervix...Endocervix,
           -Brain...Cerebellum,
           -Brain...Cortex,
           -Breast...Mammary.Tissue,
           -Fallopian.Tube,
           -Kidney...Cortex,
           -Kidney...Medulla,
           -Ovary,
           -Prostate,
           -Testis,
           -Uterus,
           -Vagina)
  gtex_tpm %>% head()
  gtex_tpm %>% nrow()
  
  names(gtex_tpm) <- c("gene_id", gtex_tissues)
  gtex_tpm %>% head
  
  
  saveRDS(object = gtex_tpm,
          file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_gene_median_tpm_tidy.rds")
  
  ## Filter the TPM file to obtain only the expression levels for the NMD genes -----------------------------
  
  NMD_genes <- c("ENSG00000005007", "ENSG00000151461", "ENSG00000169062", "ENSG00000125351")
  df_NMD_tpm <- gtex_tpm %>%
    filter(gene_id %in% NMD_genes) 
  
  
  
  df_all <- map_df(clusters, function(cluster) { 
    
    
    df_missplicing_all <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/df_missplicing_all.rds")
    
    # folder_name <- paste0(folder_root, "/", cluster)
    # df <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
    
    missplicing_local <- df_missplicing_all %>%
      filter(tissue == cluster) %>%
      select( prop_acceptor ) %>%
      sum()
    
    
    NMD_tpm_local <- df_NMD_tpm %>%
      select(all_of(cluster)) %>%
      sum()
    
    data.frame(tissue = cluster,
               missplicing = missplicing_local,
               TPM = NMD_tpm_local)
  })
  
  cor.test(x = df_all$missplicing,
           y = df_all$TPM,
           method = "pearson")
  
  
  
}

plot_TPM_NMD_gene <- function(genes = c("ENSG00000005007",
                                        "ENSG00000151461",
                                        "ENSG00000169062",
                                        "ENSG00000125351"),
                              gtex_tissues,
                              folder_name) {
  
  ## GET TPM DATA
  
  library(GSRI)
  
  gtex_tpm <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/data/GTEx_gene_median_tpm_tidy.rds")
  
  
  gtex_tpm %>% head()
  gtex_tpm %>% nrow()
  
  
  gtex_tpm <- gtex_tpm %>%
    gather(key = tissue, value = "TPM", -gene_id )
  
  
  
  
  
  
  gtex_tpm %>% head
  gtex_tpm %>% nrow()
  
  ## UPF1 = ENSG00000005007
  ## UPF2 = ENSG00000151461
  ## UPF3A = ENSG00000169062
  ## UPF3B = ENSG00000125351
  
  
  ## GET TPM VALUES FOR NMD GENES
  NMD_genes <- c("ENSG00000005007", "ENSG00000151461", "ENSG00000169062", "ENSG00000125351")
  df_gtex_tpm <- gtex_tpm %>%
    filter(gene_id %in% NMD_genes) 
  
  # Add gene name
  ind <- which(df_gtex_tpm$gene_id == "ENSG00000005007")
  df_gtex_tpm[ind,"gene_name"] <- "UPF1"
  ind <- which(df_gtex_tpm$gene_id == "ENSG00000151461")
  df_gtex_tpm[ind,"gene_name"] <- "UPF2"
  ind <- which(df_gtex_tpm$gene_id == "ENSG00000169062")
  df_gtex_tpm[ind,"gene_name"] <- "UPF3A"
  ind <- which(df_gtex_tpm$gene_id == "ENSG00000125351")
  df_gtex_tpm[ind,"gene_name"] <- "UPF3B"
  df_gtex_tpm
  
  
  
  
  df_gtex_tpm$gene_name = factor(df_gtex_tpm$gene_name, 
                                 levels = c(   "UPF2","UPF3A",  "UPF3B", "UPF1"))
  
  
  df2 = df_gtex_tpm %>% 
    ungroup() %>%
    arrange(gene_name , desc(TPM)) %>%
    mutate(tissue = fct_inorder(tissue))
  
  colours <- ifelse(str_detect(string = as.factor(df2$tissue), pattern = "Brain"), "red", "black")
  
  
  
  ggplot(df2, aes(x = tissue, 
                  y = TPM, 
                  group = gene_name, 
                  fill = gene_name)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    #coord_flip() +
    
    
    
    # ggplot(data = df_gtex_tpm_ordered) +
    # geom_bar(aes(x = reorder(tissue,TPM,FUN="mean",order = T),y = TPM), stat = 'identity') +
    #geom_point()+ 
    #geom_col(position = "stack") +
    #geom_bar(stat = "fill")  +
    theme_light() +
    #ylab("proportion of introns") +
    #xlab(xlabel) +
    scale_fill_viridis_d(option = "plasma") +
    #ggtitle("Proportion of intron type per GTEx tissue.") +
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
    guides(fill = guide_legend(title = "NMD gene: ",
                               ncol = 4, 
                               nrow = 1)) 
  
  
  
  file_name <- "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/TPM_NMD_genes.png"
  ggplot2::ggsave(filename = file_name,
                  width = 183, height = 183, units = "mm", dpi = 300)
  
}

gene_name_to_gene_ID <- function(gene_list) {
  
  library("biomaRt")
  
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  gene_ids <- getBM(attributes = 'ensembl_gene_id', 
                    filters = 'hgnc_symbol', 
                    values = gene_list, 
                    mart = ensembl)
  
  
  return(gene_ids)
  # df_data %>%
  #   filter(str_detect(string = gene_id, pattern = "ENSG00000188603"))
  # 
  # 
  # df_mina <- df_data %>%
  #   group_by(gene_id) %>%
  #   mutate(min_missplicing_d = missplicingratio_ND_tissue %>% min()) %>%
  #   mutate(max_missplicing_d = missplicingratio_ND_tissue %>% max()) %>%
  #   mutate(min_missplicing_a = missplicingratio_NA_tissue %>% min()) %>%
  #   mutate(max_missplicing_a = missplicingratio_NA_tissue %>% max()) %>%
  #   ungroup() %>%
  #   distinct(gene_id, .keep_all = T) %>%
  #   dplyr::select(c(gene_id,min_missplicing_d,max_missplicing_d,min_missplicing_a,max_missplicing_a))
  # df_mina %>% nrow()
  # df_mina %>% head()
  # 
  # library(xlsx)
  # ## Save the .xlsx object
  # write.xlsx(df_mina, file="/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/genes_missplicing_summary_FCTX.xlsx", 
  #            sheetName="frontal_cortex", 
  #            row.names=T)
  # 
  # df_mina %>% 
  #   filter(gene_id == "ENSG00000188603")
}


gene_ID_to_gene_name <- function(gene_list) {
  
  library("biomaRt")
  
  
  
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  gene_ids <- getBM(attributes = 'hgnc_symbol', 
                    filters = 'ensembl_gene_id', 
                    values = gene_list, 
                    mart = ensembl)
  
  
  return(gene_ids)
  
}


get_MT_gene_id <- function() {
  
  ## Load HUMAN REFERENCE data
  
  if (!exists("homo_sapiens_v104")) {
    
    homo_sapiens_v104 <- ensemblGenome()
    basedir(homo_sapiens_v104) <- dirname("/data/references/ensembl/gtf_gff3/v104/")
    read.gtf(homo_sapiens_v104, "/v104/Homo_sapiens.GRCh38.104.gtf")
    print("'homo_sapiens' v104 file loaded!")
    
    ## Get all transcripts
    homo_sapiens_v104_gtf <- getGtf(homo_sapiens_v104)
    transcripts_v97 <- homo_sapiens_v104_gtf %>%
      filter(feature == "transcript") 
    
  } else {
    print("'homo_sapiens' v104 file already loaded!")
  }
  
  MT_genes <- homo_sapiens_v104_gtf %>% 
    filter(str_starts(gene_name, pattern = "MT-")) 
  
  if (MT_genes %>% 
      distinct(gene_name, .keep_all = T) %>% 
      nrow() == 37) {
    saveRDS(object = MT_genes, file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/MT_genes.rds")
  } else {
    print("Error: it has been found more than 37 MT genes!")
  }
}


reduce_GO_terms <- function(query, 
                            custom_bg = NULL,
                            ordered = T,
                            significant = T) {
  
  ## GO analysis
  gene_enrichment <- gprofiler2::gost(query = query,
                                      custom_bg = custom_bg,
                                      organism = "hsapiens",
                                      ordered_query = ordered,
                                      correction_method = "bonferroni",
                                      sources = c("GO:BP", "GO:MF", "GO:CC"),
                                      significant = significant)
  
  
  if (!is.null(gene_enrichment)) {
    
    pathway_df <- data.frame(go_type = (substr(gene_enrichment$result$source, start = 4, stop = 5)),
                             go_id = gene_enrichment$result$term_id,
                             go_term = gene_enrichment$result$term_name)
    
    pathway_df <- pathway_df %>% 
      mutate(go_type = go_type %>% as.character()) %>%
      mutate(go_id = go_id %>% as.character()) %>%
      mutate(go_term = go_term %>% as.character()) %>%
      as_tibble()
    
    
    gene_enrichment_reduced <- rutils::go_reduce(pathway_df = pathway_df,
                                                 threshold = 0.95)
    
    
    go_desc_reduced <- merge(x = gene_enrichment_reduced,
                             y = gene_enrichment$result,
                             by.x = "go_id",
                             by.y = "term_id")
    
    
    
    return(go_desc_reduced %>%
             distinct(parent_term, .keep_all = T) %>% 
             dplyr::select(go_id, go_type, p_value, parent_term, term_size, intersection_size) %>%
             arrange(p_value))
  } else {
    print("No gene enrichment results to show!")
  }
}


reduce_df_to_xlsx <- function(folder_name = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/") {
  
  df_distances_all <- readRDS(file = paste0(folder_name, "/df_distances_all.rds"))
  df_missplicing_all <- readRDS(file = paste0(folder_name, "/df_missplicing_all.rds"))
  df_modulo_all <- readRDS(file = paste0(folder_name, "/df_modulo0_all.rds"))
  df_clinvar_all <- readRDS(file = paste0(folder_name, "/df_clinvar_all.rds"))
  
  
  
  library(xlsx)
  ## Save the .xlsx object
  write.xlsx(df_distances_all, file = paste0(folder_name, "/all_tissues_summary.xlsx"), sheetName = "distances", row.names = T)
  write.xlsx(df_modulo_all, file = paste0(folder_name, "/all_tissues_summary.xlsx"), sheetName = "modulo", row.names = T, append = T)
  write.xlsx(df_missplicing_all, file = paste0(folder_name, "/all_tissues_summary.xlsx"), sheetName = "mis-splicing", row.names = T, append = T)
  write.xlsx(df_clinvar_all, file = paste0(folder_name, "/all_tissues_summary.xlsx"), sheetName = "clinvar", row.names = T, append = T)
  
  
}

u2intron_scores <- function() {
  
  
  load(file = "/home/egust/Projects/Alu_exonisation/results/CNC_CDTS_CONS_gr.rda")
  
  df_mcols <- CNC_CDTS_CONS_gr %>% 
    mcols() %>% 
    attr(which = "listData")
  
  
  
  ## Load intron type files
  u12_introns <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/intron-type/minor_introns_tidy.rds")
  u2_introns <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/intron-type/major_introns_tidy.rds")
  
  
  u12_introns <- u12_introns %>%
    GRanges() %>%
    addchr()
  u2_introns <- u2_introns %>%
    GRanges() %>%
    addchr()
  
  
  u12_introns %>% width() %>% get_mode()
  u12_introns %>% width() %>% summary()
  u2_introns %>% width() %>% get_mode()
  u2_introns %>% width() %>% summary()
  
  # Check the length of each intron type
  # ggplot() +
  #   geom_freqpoly(aes(x = u12_introns %>% width(), color = "black")) +
  #   geom_freqpoly(aes(x = u2_introns %>% width(), color = "red")) +
  #   xlim(c(1,200)) 
  
  
  u12_introns_start <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                   ranges = CNC_CDTS_CONS_gr@ranges,
                                                                   strand = CNC_CDTS_CONS_gr@strand,
                                                                   mcols = df_mcols),
                                        ranges = GenomicRanges::GRanges(seqnames = u12_introns %>% seqnames(),
                                                                        ranges = IRanges(start = u12_introns %>% start() - 31, 
                                                                                         end = u12_introns %>% start() - 1),
                                                                        strand = u12_introns %>% strand()))
  u12_introns_start2 <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                    ranges = CNC_CDTS_CONS_gr@ranges,
                                                                    strand = CNC_CDTS_CONS_gr@strand,
                                                                    mcols = df_mcols),
                                         ranges = GenomicRanges::GRanges(seqnames = u12_introns %>% seqnames(),
                                                                         ranges = IRanges(start = u12_introns %>% start() + 10, 
                                                                                          end = u12_introns %>% start() + 40),
                                                                         strand = u12_introns %>% strand()))
  
  u12_introns_end <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                 ranges = CNC_CDTS_CONS_gr@ranges,
                                                                 strand = CNC_CDTS_CONS_gr@strand,
                                                                 mcols = df_mcols),
                                      ranges = GenomicRanges::GRanges(seqnames = u12_introns %>% seqnames(),
                                                                      ranges = IRanges(start = u12_introns %>% end() - 35, 
                                                                                       end = u12_introns %>% end() -5),
                                                                      strand = u12_introns %>% strand()))
  u12_introns_end2 <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                  ranges = CNC_CDTS_CONS_gr@ranges,
                                                                  strand = CNC_CDTS_CONS_gr@strand,
                                                                  mcols = df_mcols),
                                       ranges = GenomicRanges::GRanges(seqnames = u12_introns %>% seqnames(),
                                                                       ranges = IRanges(start = u12_introns %>% end() + 1, 
                                                                                        end = u12_introns %>% end() + 31),
                                                                       strand = u12_introns %>% strand()))
  # u12_introns_start <- c(u12_introns_start, u12_introns_start2)
  # u12_introns_end <- c(u12_introns_end, u12_introns_end2)
  
  ##
  
  u2_introns_start <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                  ranges = CNC_CDTS_CONS_gr@ranges,
                                                                  strand = CNC_CDTS_CONS_gr@strand,
                                                                  mcols = df_mcols),
                                       ranges = GenomicRanges::GRanges(seqnames = u2_introns %>% seqnames(),
                                                                       ranges = IRanges(start = u2_introns %>% start() - 31, 
                                                                                        end = u2_introns %>% start() - 1),
                                                                       strand = u2_introns %>% strand()))
  u2_introns_start2 <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                   ranges = CNC_CDTS_CONS_gr@ranges,
                                                                   strand = CNC_CDTS_CONS_gr@strand,
                                                                   mcols = df_mcols),
                                        ranges = GenomicRanges::GRanges(seqnames = u2_introns %>% seqnames(),
                                                                        ranges = IRanges(start = u2_introns %>% start() + 10, 
                                                                                         end = u2_introns %>% start() + 40),
                                                                        strand = u2_introns %>% strand()))
  
  u2_introns_end <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                ranges = CNC_CDTS_CONS_gr@ranges,
                                                                strand = CNC_CDTS_CONS_gr@strand,
                                                                mcols = df_mcols),
                                     ranges = GenomicRanges::GRanges(seqnames = u2_introns %>% seqnames(),
                                                                     ranges = IRanges(start = u2_introns %>% end() - 35, 
                                                                                      end = u2_introns %>% end() - 5),
                                                                     strand = u2_introns %>% strand()))
  
  u2_introns_end2 <- subsetByOverlaps(x = GenomicRanges::GRanges(seqnames = CNC_CDTS_CONS_gr@seqnames,
                                                                 ranges = CNC_CDTS_CONS_gr@ranges,
                                                                 strand = CNC_CDTS_CONS_gr@strand,
                                                                 mcols = df_mcols),
                                      ranges = GenomicRanges::GRanges(seqnames = u2_introns %>% seqnames(),
                                                                      ranges = IRanges(start = u2_introns %>% end() + 1, 
                                                                                       end = u2_introns %>% end() + 31),
                                                                      strand = u2_introns %>% strand()))
  
  
  # u2_introns_start <- c(u2_introns_start, u2_introns_start2)
  # u2_introns_end <- c(u2_introns_end, u2_introns_end2)
  
  
  
  ## Plot - CDTS (constraint scores)  -------------------------------------------------------------------------------------------------------------
  
  
  ggplot() +
    geom_density(aes(x = u2_introns_start$mcols.CDTS, color = "black")) +
    geom_density(aes(x = u2_introns_start2$mcols.CDTS, color = "blue")) +
    geom_density(aes(x = u12_introns_start$mcols.CDTS, color = "red")) +
    geom_density(aes(x = u12_introns_start2$mcols.CDTS, color = "darkgreen")) +
    ggtitle("CDTS scores - 5' intronic region") +
    xlab("CDTS") +
    theme_light() +
    scale_color_manual(values = c("black", "blue", "red", "darkgreen"),
                       breaks = c("black", "blue", "red", "darkgreen"),
                       labels = c("u2-introns 5' exon", "u2-introns 5' intron", "u12-introns 5' exon", "u12-introns 5' intron")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Intron type: ",
                                ncol = 2, 
                                nrow = 2))
  
  ggplot2::ggsave(filename = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/minorintron3-constraint-5ss.png", 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  ggplot() +
    geom_density(aes(x = u2_introns_start$mcols.CDTS, color = "black")) +
    geom_density(aes(x = u2_introns_start2$mcols.CDTS, color = "blue")) +
    geom_density(aes(x = u12_introns_start$mcols.CDTS, color = "red")) +
    geom_density(aes(x = u12_introns_start2$mcols.CDTS, color = "darkgreen")) +
    ggtitle("CDTS scores - 3' intronic region") +
    xlab("CDTS") +
    theme_light() +
    scale_color_manual(values = c("black", "blue", "red", "darkgreen"),
                       breaks = c("black", "blue", "red", "darkgreen"),
                       labels = c("u2-introns 5' exon", "u2-introns 5' intron", "u12-introns 5' exon", "u12-introns 5' intron")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Intron type: ",
                                ncol = 2, 
                                nrow = 2))
  ggplot2::ggsave(filename = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/minorintron3-constraint-3ss.png", 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  ## Plot - Conservation --------------------------------------------------------------------------------------------
  
  
  # conservation_u12intron_s <- u12_introns_start$mcols.mean_phastCons20way
  # conservation_u12intron_e <- u12_introns_end$mcols.mean_phastCons20way
  # 
  # 
  # conservation_u2intron_s <- u2_introns_start$mcols.mean_phastCons20way
  # if(which(conservation_u2intron_s %>% is.na()) %>% length > 0) {
  #   conservation_u2intron_s <- u2_introns_start[-which(conservation_u2intron_s %>% is.na())]$mcols.mean_phastCons20way
  # }
  # 
  # conservation_u2intron_e <- u2_introns_end$mcols.mean_phastCons20way
  # if(which(conservation_u2intron_e %>% is.na()) %>% length > 0) {
  #   conservation_u2intron_e <- u2_introns_end[-which(conservation_u2intron_e %>% is.na())]$mcols.mean_phastCons20way
  # }
  
  
  
  
  ggplot() +
    geom_density(aes(x = u2_introns_start$mcols.mean_phastCons20way, color = "black")) +
    geom_density(aes(x = u2_introns_start2$mcols.mean_phastCons20way, color = "blue")) +
    geom_density(aes(x = u12_introns_start$mcols.mean_phastCons20way, color = "red")) +
    geom_density(aes(x = u12_introns_start2$mcols.mean_phastCons20way, color = "green")) +
    ggtitle("mean phastCons20 - 5' intronic region") +
    xlab("mean phastCons20") +
    theme_light() +
    scale_color_manual(values = c("black", "blue", "red", "green"),
                       breaks = c("black", "blue", "red", "green"),
                       labels = c("u2-introns 5' exon", "u2-introns 5' intron", "u12-introns 5' exon", "u12-introns 5' intron")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Intron type: ",
                                ncol = 2, 
                                nrow = 2))
  
  ggplot2::ggsave(filename = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/minorintron-conservation-5ss.png", 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  ggplot() +
    geom_density(aes(x = u2_introns_end$mcols.mean_phastCons20way, color = "black")) +
    geom_density(aes(x = u2_introns_end2$mcols.mean_phastCons20way, color = "blue")) +
    geom_density(aes(x = u12_introns_end$mcols.mean_phastCons20way, color = "red")) +
    geom_density(aes(x = u12_introns_end2$mcols.mean_phastCons20way, color = "green")) +
    ggtitle("mean phastCons20 - 3' intronic region") +
    xlab("mean phastCons20") +
    theme_light() +
    scale_color_manual(values = c("black", "blue", "red", "green"),
                       breaks = c("black", "blue", "red", "green"),
                       labels = c("u2-introns 3' intron", "u2-introns 3' exon", "u12-introns 3' intron", "u12-introns 3' exon")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Intron type: ",
                                ncol = 2, 
                                nrow = 2))
  
  ggplot2::ggsave(filename = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/minorintron-conservation-3ss.png", 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  ## Plot - CNCRs -------------------------------------------------------------------------------------------------------------
  
  ggplot() +
    geom_density(aes(x = u2_introns_start$mcols.CNC, color = "black")) +
    geom_density(aes(x = u2_introns_start2$mcols.CNC, color = "blue")) +
    geom_density(aes(x = u12_introns_start$mcols.CNC, color = "red")) +
    geom_density(aes(x = u12_introns_start2$mcols.CNC, color = "darkgreen")) +
    ggtitle("CNCR scores - 5' intronic region") +
    xlab("CNCRs") +
    theme_light() +
    scale_color_manual(values = c("black", "blue", "red", "darkgreen"),
                       breaks = c("black", "blue", "red", "darkgreen"),
                       labels = c("u2-introns 3' intron", "u2-introns 3' exon", "u12-introns 3' intron", "u12-introns 3' exon")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Intron type: ",
                                ncol = 2, 
                                nrow = 2))
  
  
  ggplot() +
    geom_density(aes(x = u2_introns_start$mcols.CNC, color = "black")) +
    geom_density(aes(x = u2_introns_start2$mcols.CNC, color = "blue")) +
    geom_density(aes(x = u12_introns_start$mcols.CNC, color = "red")) +
    geom_density(aes(x = u12_introns_start2$mcols.CNC, color = "darkgreen")) +
    ggtitle("CNCR scores - 3' intronic region") +
    xlab("CNCRs") +
    theme_light() +
    scale_color_manual(values = c("black", "blue", "red", "darkgreen"),
                       breaks = c("black", "blue", "red", "darkgreen"),
                       labels = c("u2-introns 3' intron", "u2-introns 3' exon", "u12-introns 3' intron", "u12-introns 3' exon")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Intron type: ",
                                ncol = 2, 
                                nrow = 2))
  
  
}

add_coordinates_to_novel_junc_distances <- function(clusters) {
  
  
  for (cluster in clusters) {
    
    print(cluster)
    
    all_split_reads_details_97 <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/base_data/", 
                                                        cluster, "/", cluster, "_annotated_SR_details_length.rds")) %>%
      as.data.frame()
    
    all_split_reads_details_97_novel <- all_split_reads_details_97 %>%
      filter(type %in% c("novel_donor", "novel_acceptor")) %>%
      select(seqnames,start,end,width, strand, junID)
    
    
    
    distances <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/distances/", cluster, "/", cluster, "_distances_tidy_all.rds"))
    
    
    
    distances_gr <- merge(x = distances %>% select(-seqnames,-start,-end,-strand),
                          y = all_split_reads_details_97_novel,
                          by.x = "novel_junID",
                          by.y = "junID",
                          all.x = T)
    
    
    saveRDS(object = distances_gr,
            file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/distances/", cluster, "/", cluster, "_distances_tidy_all_gr.rds"))
  }
  
  
  
  
}

check_TPM_changes <- function(cluster,
                              cluster2,
                              folder_root) {
  
  ## LOAD SOURCE DATA FROM CLUSTER 1
  file_name <- paste0(folder_root, "/", cluster, "/", cluster, "_db_introns.rds")
  df_introns_cluster1 <- readRDS(file = file_name)
  
  df_introns_cluster1 %>% head()
  df_introns_cluster1 %>% nrow()
  
  
  ## LOAD SOURCE DATA FROM CLUSTER 2
  file_name <- paste0(folder_root, "/", cluster2, "/", cluster2, "_db_introns.rds")
  df_introns_cluster2 <- readRDS(file = file_name)
  
  df_introns_cluster2 %>% head()
  df_introns_cluster2 %>% nrow()
  
  
  ## CHECK DIFFERENCES IN TPM
  df_introns_cluster1 %>%
    filter(gene_name == "SNCA") %>%
    dplyr::select(ref_missplicing_ratio_tissue_ND, ref_missplicing_ratio_tissue_NA, tpm)
  
  
  df_introns_cluster2 %>%
    filter(gene_name == "SNCA") %>%
    dplyr::select(ref_missplicing_ratio_tissue_ND, ref_missplicing_ratio_tissue_NA, tpm = TPM)
}


get_wilcoxon_test <- function(cluster,
                              folder_name) {
  
  print(paste0(cluster))
  
  
  df_stats <- readRDS(file = paste0(folder_name, "/", cluster, "_missplicingratio_tidy_v2.rds"))
  df_stats %>% head()
  df_stats %>% nrow()
  
  
  print(paste0("MSR Novel donor: ", 
               df %>%
                 filter(missplicingratio_ND_tissue > 0 ) %>%
                 pull(ref_junID) %>%
                 unique() %>%
                 length()))
  
  print(paste0("MSR Novel acceptor: ", 
               df %>%
                 filter(missplicingratio_NA_tissue > 0 ) %>%
                 pull(ref_junID) %>%
                 unique() %>%
                 length()))
  
  
  print(paste0("Reference junctions type U12-intron: ", 
               df %>%
                 filter(u12_intron == T) %>%
                 pull(ref_junID) %>%
                 unique() %>%
                 length()))
  
  print(paste0("Reference junctions type U2-intron: ", 
               df %>%
                 filter(u2_intron == T) %>%
                 pull(ref_junID) %>%
                 unique() %>%
                 length()))
  
  
  ######################################################################
  ######################################################################
  ######################################################################
  
  df_discard <- df_stats %>% 
    distinct(ref_junID, .keep_all = T)  %>%
    filter(u2_intron == F, u12_intron == F) 
  df_discard %>% nrow()
  df_discard %>%
    filter(type == "none") %>%
    nrow()
  df_discard %>%
    filter(type == "both") %>%
    nrow()
  df_discard %>%
    filter(type == "donor") %>%
    nrow()
  df_discard %>%
    filter(type == "acceptor") %>%
    nrow()
  
  
  df_data <- df_stats %>% 
    distinct(ref_junID, .keep_all = T)  %>%
    filter(u2_intron == T | u12_intron == T)
  df_data %>% nrow()
  
  
  
  
  
  ## Replace NA's by zeros
  #df_data[is.na(df_data$missplicingratio_ND_tissue),]$missplicingratio_ND_tissue <- 0#1/(df_data[is.na(df_data$missplicingratio_ND_tissue),]$ref_counts + 1)
  #df_data[is.na(df_data$missplicingratio_NA_tissue),]$missplicingratio_NA_tissue <- 0#1/(df_data[is.na(df_data$missplicingratio_NA_tissue),]$ref_counts + 1)
  
  
  df_data %>% head()
  
  
  
  
  ## Summary of the data
  
  summary(df_data$missplicingratio_ND_tissue)
  summary(df_data$missplicingratio_NA_tissue)
  
  which(!is.na(df_data$missplicingratio_ND_tissue)) %>% length()
  which(!is.na(df_data$missplicingratio_NA_tissue)) %>% length()
  
  which(df_data$missplicingratio_ND_tissue == 0) %>% length()
  which(df_data$missplicingratio_ND_tissue > 0) %>% length()
  
  which(df_data$missplicingratio_NA_tissue == 0) %>% length()
  which(df_data$missplicingratio_NA_tissue > 0) %>% length()
  
  
  ## Statistical tests (missing values are ignored)
  
  wilcox.test(x = df_data$missplicingratio_ND_tissue,
              y = df_data$missplicingratio_NA_tissue,
              paired = T,
              correct = T,
              alternative = "less")
  
  
  
  
  
  ######################################################
  ## WILCOXON TESTS
  ######################################################
  
  ## Reference junctions from protein-coding transcripts present more stringent (smaller) mis-splicing values
  
  df_protein <- df_data %>% 
    filter(protein_coding == 100)
  df_protein %>% nrow()
  
  df_notprotein <- df_data %>% 
    filter(protein_coding == 0)
  df_notprotein %>% nrow()
  
  ## TEST AT DONOR SPLICE SITES
  
  wilcox.test(x = df_protein$missplicingratio_ND_tissue,
              y =  df_notprotein$missplicingratio_ND_tissue,
              correct = T,
              paired = F,
              alternative = "less")
  
  
  ggplot() + 
    geom_density(aes(x = df_protein$missplicingratio_ND_tissue, fill = "#440154FF"), 
                 alpha = 0.8) +
    geom_density(aes(x = df_notprotein$missplicingratio_ND_tissue, fill = "#35B779FF"), 
                 alpha = 0.8) +
    
    ggtitle(paste0("Mis-splicing ratio mean at donor splice sites.\n108 samples - ", tissue)) +
    xlab("Mis-splicing ratio mean at donor sites") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    xlim(c(0, 1)) +
    ylim(c(0, 17)) + 
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("non-protein-coding","protein-coding")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  ## Save the results
  folder_img <- paste0(folder_name, "/images/")
  file_name <- paste0(folder_img, cluster, "_missplicingratio_donor_PC.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  ## TEST AT ACCEPTOR SPLICE SITES
  
  wilcox.test(x = df_protein$missplicingratio_NA_tissue,
              y = df_notprotein$missplicingratio_NA_tissue,
              correct = T,
              paired = F,
              alternative = "less")
  
  
  ggplot() + 
    geom_density(aes(x = df_protein$missplicingratio_NA_tissue, fill = "#440154FF"), 
                 alpha = 0.8) +
    geom_density(aes(x = df_notprotein$missplicingratio_NA_tissue, fill = "#35B779FF"), 
                 alpha = 0.8) +
    
    ggtitle(paste0("Mis-splicing ratio mean at acceptor splice sites.\n108 samples - ", tissue)) +
    xlab("Mis-splicing ratio mean at acceptor sites") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    xlim(c(0, 1)) +
    ylim(c(0, 17)) + 
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("non-protein-coding","protein-coding")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  ## Save the results
  folder_img <- paste0(folder_name, "/images/")
  file_name <- paste0(folder_img, cluster, "_missplicingratio_acceptor_PC.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  ######################################################
  ######################################################
  
  
  ## Reference junctions from u12-introns present more stringent (smaller) mis-splicing values
  ## than u2-introns
  df_data %>% nrow()
  
  df_u12intron <- df_data %>% 
    filter(u12_intron == T)
  df_u12intron %>% nrow()
  
  
  df_u2intron <- df_data %>% 
    filter(u2_intron == T)
  df_u2intron %>% nrow()
  
  
  ## NOVEL DONOR TEST
  
  wilcox.test(x = df_u12intron$missplicingratio_ND_tissue,
              y = df_u2intron$missplicingratio_ND_tissue,
              correct = T,
              paired = F,
              alternative = "less")
  
  
  ggplot() + 
    geom_density(aes(x = df_u2intron$missplicingratio_ND_tissue, fill = "#35B779FF"), 
                 alpha = 0.8) +
    geom_density(aes(x = df_u12intron$missplicingratio_ND_tissue, fill = "#440154FF"), 
                 alpha = 0.8) +
    
    
    ggtitle(paste0("Mis-splicing ratio mean at donor splice sites.\n108 samples - ", tissue)) +
    xlab("Mis-splicing ratio mean at donor sites") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    xlim(c(0, 1)) +
    ylim(c(0, 20)) + 
    theme_light() +
    scale_fill_manual(values = c("#440154FF","#35B779FF"),
                      breaks = c("#440154FF","#35B779FF"),
                      labels = c("u12-intron","u2-intron")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  ## Save the results
  folder_img <- paste0(folder_name, "/images/")
  file_name <- paste0(folder_img, cluster, "_missplicingratio_donor_introntype.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  
  ## NOVEL ACCETOR TEST
  
  
  wilcox.test(x = df_u12intron$missplicingratio_NA_tissue,
              y = df_u2intron$missplicingratio_NA_tissue,
              #exact = T,
              correct = T,
              paired = F,
              alternative = "less")
  
  
  ggplot() + 
    geom_density(aes(x = df_u2intron$missplicingratio_NA_tissue, fill = "#35B779FF"), 
                 alpha = 0.8) +
    geom_density(aes(x = df_u12intron$missplicingratio_NA_tissue, fill = "#440154FF"), 
                 alpha = 0.8) +
    
    
    ggtitle(paste0("Mis-splicing ratio mean at acceptor splice sites.\n108 samples - ", tissue)) +
    xlab("Mis-splicing ratio mean at acceptor sites") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    xlim(c(0, 1)) +
    ylim(c(0, 20)) + 
    theme_light() +
    scale_fill_manual(values = c("#440154FF","#35B779FF"),
                      breaks = c("#440154FF","#35B779FF"),
                      labels = c("u12-intron","u2-intron")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  ## Save the results
  folder_img <- paste0(folder_name, "/images/")
  file_name <- paste0(folder_img, cluster, "_missplicingratio_acceptor_introntype.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  ######################################################
  ######################################################
  
  
  ## Reference junctions from u12-introns & protein-coding transcripts present more stringent (smaller) mis-splicing values
  
  
  df_u2intron_p <- df_data %>% 
    filter(u12_intron == F, protein_coding == 100)
  df_u2intron_p %>% nrow()
  
  
  df_u12intron_p <- df_data %>% 
    filter(u12_intron == T, protein_coding == 100)
  df_u12intron_p %>% nrow()
  
  
  
  
  ## NOVEL DONOR TEST
  
  wilcox.test(x = df_u12intron_p$missplicingratio_ND_tissue,
              y = df_u2intron_p$missplicingratio_ND_tissue,
              paired = F,
              exact = F,
              alternative = "less")
  
  ggplot() + 
    geom_density(aes(x = df_u2intron_p$missplicingratio_ND_tissue, fill = "#35B779FF"), 
                 alpha = 0.8) +
    geom_density(aes(x = df_u12intron_p$missplicingratio_ND_tissue, fill = "#440154FF"), 
                 alpha = 0.8) +
    
    
    ggtitle(paste0("Mis-splicing ratio mean at donor splice sites.\nOnly protein-coding transcripts.\n108 samples - ", tissue)) +
    xlab("Mis-splicing ratio mean at donor sites") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    xlim(c(0, 1)) +
    ylim(c(0, 20)) + 
    theme_light() +
    scale_fill_manual(values = c("#440154FF","#35B779FF"),
                      breaks = c("#440154FF","#35B779FF"),
                      labels = c("u12-intron","u2-intron")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  ## Save the results
  folder_img <- paste0(folder_name, "/images/")
  file_name <- paste0(folder_img, cluster, "_missplicingratio_donor_introntype_PC.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  
  ## NOVEL ACCEPTOR TEST
  
  
  wilcox.test(x = df_u12intron_p$missplicingratio_NA_tissue,
              y =  df_u2intron_p$missplicingratio_NA_tissue,
              paired = F,
              exact = F,
              alternative = "less")
  
  
  ggplot() + 
    geom_density(aes(x = df_u2intron_p$missplicingratio_NA_tissue, fill = "#35B779FF"), 
                 alpha = 0.8) +
    geom_density(aes(x = df_u12intron_p$missplicingratio_NA_tissue, fill = "#440154FF"), 
                 alpha = 0.8) +
    
    
    ggtitle(paste0("Mis-splicing ratio mean at acceptor splice sites.\nOnly protein-coding transcripts.\n108 samples - ", tissue)) +
    xlab("Mis-splicing ratio mean at acceptor sites") +
    #ylab("Intron count") +
    #scale_y_continuous(limits = c(0, 85000)) +
    xlim(c(0, 1)) +
    ylim(c(0, 20)) + 
    theme_light() +
    scale_fill_manual(values = c("#440154FF","#35B779FF"),
                      breaks = c("#440154FF","#35B779FF"),
                      labels = c("u12-intron","u2-intron")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1))
  
  
  ## Save the results
  folder_img <- paste0(folder_name, "/images/")
  file_name <- paste0(folder_img, cluster, "_missplicingratio_acceptor_introntype_PC.png")
  ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
  
  ######################################################
  ######################################################
  
  # 
  # ## Reference junctions from u12-introns & non-protein-coding transcripts present more stringent (smaller) mis-splicing values
  # ## than u2-introns
  # 
  # 
  # df_u12intron_p <- df_data %>% 
  #   filter(u12_intron == T, protein_coding == 0)
  # df_u12intron_p %>% nrow()
  # 
  # 
  # df_u2intron_p <- df_data %>% 
  #   filter(u2_intron == T, protein_coding == 0)
  # df_u2intron_p %>% nrow()
  # 
  # #df_nd_p_u2 <- sample_n(df_nd_p_u2, df_nd_p_u12 %>% nrow())
  # 
  # wilcox.test(x = df_u12intron_p$missplicingratio_ND_tissue,
  #             y = df_u2intron_p$missplicingratio_ND_tissue,
  #             paired = F,
  #             exact = F,
  #             alternative = "less")
  # 
  # ##
  # 
  # 
  # wilcox.test(x = df_u12intron_p$missplicingratio_NA_tissue,
  #             y = df_u2intron_p$missplicingratio_NA_tissue,
  #             paired = F,
  #             exact = F,
  #             alternative = "less")
  
}

get_missplicing_ratio_Mina <- function(cluster,
                                       samples,
                                       split_read_counts,
                                       folder_name,
                                       folder_save_name,
                                       db = NULL) {
  
  
  df_tidy <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_tidy_all.rds"))
  
  print(paste0(Sys.time(), " - Getting mis-splicing ratio for mis-spliced junctions."))
  
  
  df_tidy %>% head()
  df_tidy %>% nrow()
  
  print(paste0(Sys.time(), " - ", df_tidy %>%
                 distinct(ref_junID) %>%
                 nrow(), " number of reference junctions"))
  
  df_tidy %>%
    filter(is.na(ref_counts)) %>%
    nrow()
  
  ## Threshold for the minimum number of reads
  threshold <- 10
  
  
  ## 50% of the total number of samples
  print(paste0("Threshold samples: ", round(samples %>% length() * 0.5)))
  print(paste0("Threshold counts: ", (round(samples %>% length() * 0.5) * threshold)  / samples %>% length()))
  
  
  ## Obtain only reference junction expressed in 50% of samples presenting at least 5 reads
  df <- df_tidy %>%
    filter(ref_n_individuals >= round(samples %>% length() * 0.5), 
           ref_mean_counts >= (round(samples %>% length() * 0.5) * threshold) / samples %>% length(),
           !is.na(ref_counts))
  
  
  
  ## Get novel junctions ---------------------------------------------------------------------------------------------------------------------------
  
  df_novel <- df %>%
    group_by(novel_junID, sample) %>% 
    dplyr::mutate(novel_missplicing_ratio_sample = missplicing_ratio %>% mean()) %>% 
    ungroup() %>%
    group_by(novel_junID) %>% 
    dplyr::mutate(novel_missplicing_ratio_tissue = sum(novel_missplicing_ratio_sample)/ref_n_individuals) %>% 
    ungroup()
  
  df_novel$novel_missplicing_ratio_tissue %>% summary
  
  df_novel %>% head
  df_novel %>% nrow
  
  df_novel %>%
    filter(novel_junID == "67197814") %>% 
    select(novel_junID, sample, novel_counts, ref_counts, ref_n_individuals, missplicing_ratio, novel_missplicing_ratio_sample, novel_missplicing_ratio_tissue) %>% 
    as.data.frame()
  
  
  db_novel <- df_novel %>% 
    as.data.frame() %>%
    select(-ref_counts,
           -ref_width,
           -ref_ss5score,
           -ref_ss3score,
           -missplicing_ratio, 
           -ref_n_individuals,
           -ref_mean_counts,
           -gene_name_end,
           -start,
           -end,
           -u12_intron,
           -u2_intron) %>%
    distinct(novel_junID, .keep_all = T) %>% 
    dplyr::rename(novel_type = type,
                  gene_name = gene_name_start)
  
  db_novel %>%
    filter(novel_junID == "67197814")
  
  
  saveRDS(object = db_novel ,
          file = paste0(folder_save_name, "/", cluster, "_db_novel.rds"))
  
  print(paste0(Sys.time(), " - Novel junctions DB saved!"))
  
  ## Get novel junctions details ----------------------------------------------------------------------------
  
  
  split_read_counts_novel <- split_read_counts %>%
    filter(junID %in% (df_novel$novel_junID %>% unique())) %>%
    dplyr::select(junID, samples %>% as.character()) %>%
    #mutate_all(~replace(., is.na(.), 0)) %>%
    as_tibble()
  
  split_read_counts_novel_details <- split_read_counts_novel %>% #head() %>%
    gather(sample, count, samples %>% as.character()) %>%
    #dplyr::filter(count > 0) %>%
    dplyr::rename(novel_junID = junID,
                  novel_counts = count)
  
  split_read_counts_novel_details %>%
    as.data.frame() %>%
    filter(novel_junID == "67197814")
  
  # db_novel_details <- df_novel %>% 
  #   as.data.frame() %>%
  #   select(novel_junID,
  #          sample,
  #          novel_type = type,
  #          novel_counts) %>%
  #   mutate(age = "",
  #          sex = "",
  #          mapped_read_count = 0,
  #          avg_read_length = 0)
  # 
  # 
  # df_novel %>%
  #   filter(novel_junID == "67197814") %>%
  #   select(novel_junID,
  #          sample,
  #          novel_n_individuals,
  #          novel_type = type,
  #          novel_counts)
  # split_read_counts_novel_details %>%
  #   filter(junID == "67197814") %>%
  #   as.data.frame()
  
  for (sample in split_read_counts_novel_details$sample %>% unique()) { # sample <- db_novel_details$sample[1]
    
    if (!is.null(db)) {
      
      ind <- which(db$run %>% as.character() == sample)
      people <- db[ind, ]
      
      if (str_detect(people$characteristics[1], "female")) {
        sex <- "female"
      }else if (str_detect(people$characteristics[1], "male")) {
        sex <- "male"
      }
      
      data <- people$characteristics[1]
      position <- str_locate(data, "age at death")
      age <- substr(data, position[2] + 3, position[2] + 4) %>% 
        as.integer()
      
      mapped_read_count <- people$mapped_read_count
      avg_read_length <- people$avg_read_length
      
    } else {
      
      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/data/splicing_tolerance/splicing_tolerance.sqlite")
      query <- DBI::dbSendQuery(con, paste0("SELECT * FROM GTEX_info WHERE sample_recount_id == '", sample,"';")) #  WHERE smafrze == 'USE ME'
      people <- DBI::dbFetch(query)
      DBI::dbClearResult(query)
      DBI::dbDisconnect(con)
      
      if (str_detect(people$bigwig_file[1], "female")) {
        sex <- "female"
      } else if (str_detect(people$bigwig_file[1], "male")) {
        sex <- "male"
      }
      
      age <- people$age
      mapped_read_count <- people$mapped_read_count
      avg_read_length <- people$avg_read_length
    }
    
    
    
    ind <- which(split_read_counts_novel_details$sample == sample)
    split_read_counts_novel_details[ind,"age"] <- age
    split_read_counts_novel_details[ind,"sex"] <- sex
    split_read_counts_novel_details[ind,"mapped_read_count"] <- mapped_read_count
    split_read_counts_novel_details[ind,"avg_read_length"] <- avg_read_length
    
  }
  
  
  
  saveRDS(object = split_read_counts_novel_details,
          file = paste0(folder_save_name, "/", cluster, "_db_novel_details.rds"))
  
  print(paste0(Sys.time(), " - Details of novel junctions DB saved!"))
  
  
  
  ## Generate reference introns data base ---------------------------------------------------------------------------
  
  
  df_ref <- df %>%
    group_by(ref_junID, sample, type) %>% 
    dplyr::mutate(ref_missplicing_ratio_sample = sum(novel_counts)/(sum(novel_counts) + ref_counts)) %>% 
    distinct(ref_missplicing_ratio_sample, .keep_all = T) %>%
    ungroup() 
  
  
  
  
  df_ref <- df_ref %>%
    spread(key = type, value = ref_missplicing_ratio_sample) %>% 
    group_by(ref_junID) %>% 
    dplyr::mutate(ref_missplicing_ratio_tissue_ND = (sum(novel_donor, na.rm = T) / ref_n_individuals)) %>%
    dplyr::mutate(ref_missplicing_ratio_tissue_NA = (sum(novel_acceptor, na.rm = T) / ref_n_individuals)) %>% 
    ungroup() 
  
  
  
  db_ref <- df_ref %>% 
    as.data.frame() %>%
    distinct(ref_junID, .keep_all = T) %>%
    select(-novel_junID, 
           -sample, 
           -distance,
           -novel_counts,
           -novel_width,
           -novel_ss5score,
           -novel_ss3score,
           -ref_counts,
           -missplicing_ratio,
           -novel_n_individuals,
           -novel_mean_counts,
           -gene_name_end,
           -novel_start,
           -novel_end,
           -novel_donor,
           -novel_acceptor) %>%
    dplyr::rename(gene_name = gene_name_start)
  
  
  
  
  #print(paste0(Sys.time(), " - Getting mis-splicing ratio for never mis-spliced junctions."))
  
  ## APPLY THRESHOLD TO NEVER MIS-SPLICED
  
  df_never <- df_tidy %>% 
    filter(is.na(ref_counts), 
           ref_n_individuals >= round(samples %>% length() * 0.5),
           ref_mean_counts >= (round(samples %>% length() * 0.5) * threshold)  / samples %>% length())
  
  db_never <- df_never %>% 
    select(ref_junID,
           strand,
           ref_width,
           ref_ss5score,
           ref_ss3score,
           ref_n_individuals,
           ref_mean_counts,
           seqnames,
           start,
           end,
           gene_id,
           gene_name = gene_name_start,
           u12_intron,
           u2_intron)
  
  
  
  
  db_ref_joined <- plyr::rbind.fill(db_ref, db_never)
  db_ref_joined %>% nrow()
  db_ref_joined %>% head()
  
  
  
  #print(paste0(Sys.time(), " - Tidying results..."))
  
  if (db_ref_joined %>%
      filter(is.na(ref_missplicing_ratio_tissue_ND),
             is.na(ref_missplicing_ratio_tissue_NA)) %>%
      distinct(ref_junID) %>%
      nrow() != df_never %>%
      #filter(is.na(novel_junID)) %>%
      distinct(ref_junID) %>%
      nrow()) {
    print("Error: disimilar number of never mis-spliced junctions")
  }
  
  
  # db_ref_joined %>% 
  #    distinct(ref_junID, .keep_all = T)
  
  
  db_ref_joined[is.na(db_ref_joined$ref_missplicing_ratio_tissue_ND),"ref_missplicing_ratio_tissue_ND"] <- 0
  db_ref_joined[is.na(db_ref_joined$ref_missplicing_ratio_tissue_NA),"ref_missplicing_ratio_tissue_NA"] <- 0
  db_ref_joined %>% head()
  db_ref_joined %>% nrow()
  
  
  db_ref_joined$ref_missplicing_ratio_tissue_ND %>% summary()
  db_ref_joined$ref_missplicing_ratio_tissue_NA %>% summary()
  
  
  print(paste0(Sys.time(), " - ", db_ref_joined %>% nrow(), " final number of ref junctions."))
  #print(paste0(Sys.time(), " - ", db_ref_joined %>% distinct(ref_junID) %>% nrow(), " final number of unique ref junctions."))
  
  
  
  saveRDS(object = db_ref_joined,
          file = paste0(folder_save_name, "/", cluster, "_db_introns.rds"))
  
  #print(paste0(Sys.time(), " - Results saved!"))
  print(paste0(Sys.time(), " - Introns DB saved!"))
  
  
  
  
  ## Get reference introns details ----------------------------------------------------------------------------
  
  split_read_counts_ref <- split_read_counts %>%
    filter(junID %in% (db_ref_joined$ref_junID %>% unique())) %>%
    dplyr::select(junID, samples %>% as.character()) %>%
    #mutate_all(~replace(., is.na(.), 0)) %>%
    as_tibble()
  
  split_read_counts_ref_details <- split_read_counts_ref %>% #head() %>%
    gather(sample, count, samples %>% as.character()) %>%
    #dplyr::filter(count > 0) %>%
    dplyr::rename(ref_junID = junID)
  
  
  
  for (sample in split_read_counts_ref_details$sample %>% unique()) { # sample <- split_read_counts_ref_details$sample[1]
    
    
    
    if (!is.null(db)) {
      
      ind <- which(db$run %>% as.character() == sample)
      people <- db[ind, ]
      
      if (str_detect(people$characteristics[1], "female")) {
        sex <- "female"
      }else if (str_detect(people$characteristics[1], "male")) {
        sex <- "male"
      }
      
      data <- people$characteristics[1]
      position <- str_locate(data, "age at death")
      age <- substr(data, position[2] + 3, position[2] + 4) %>% 
        as.integer()
      
      mapped_read_count <- people$mapped_read_count
      avg_read_length <- people$avg_read_length
      
    } else {
      
      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/data/splicing_tolerance/splicing_tolerance.sqlite")
      query <- DBI::dbSendQuery(con, paste0("SELECT * FROM GTEX_info WHERE sample_recount_id == '", sample,"';")) #  WHERE smafrze == 'USE ME'
      people <- DBI::dbFetch(query)
      DBI::dbClearResult(query)
      DBI::dbDisconnect(con)
      
      if (str_detect(people$bigwig_file[1], "female")) {
        sex <- "female"
      }else if (str_detect(people$bigwig_file[1], "male")) {
        sex <- "male"
      }
      
      age <- people$age
      mapped_read_count <- people$mapped_read_count
      avg_read_length <- people$avg_read_length
    }
    
    # con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/data/splicing_tolerance/splicing_tolerance.sqlite")
    # query <- DBI::dbSendQuery(con, paste0("SELECT * FROM GTEX_info WHERE sample_recount_id == '", sample,"';")) #  WHERE smafrze == 'USE ME'
    # people <- DBI::dbFetch(query)
    # DBI::dbClearResult(query)
    # DBI::dbDisconnect(con)
    # 
    # if (str_detect(people$bigwig_file[1], "female")) {
    #   sex <- "female"
    # }else if (str_detect(people$bigwig_file[1], "male")) {
    #   sex <- "male"
    # }
    
    ind <- which(split_read_counts_ref_details$sample == sample)
    split_read_counts_ref_details[ind,"age"] <- age
    split_read_counts_ref_details[ind,"sex"] <- sex
    split_read_counts_ref_details[ind,"mapped_read_count"] <- mapped_read_count
    split_read_counts_ref_details[ind,"avg_read_length"] <- avg_read_length
    
    #split_read_counts_ref_details
    print(sample)
  }
  
  
  saveRDS(object = split_read_counts_ref_details,
          file = paste0(folder_save_name, "/", cluster, "_db_introns_details.rds"))
  
  print(paste0(Sys.time(), " - Details of introns DB saved!"))
  
  # ggplot(data = db_ref_joined %>% 
  #          distinct(ref_junID, .keep_all = T)) + 
  #   geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"), 
  #                alpha = 0.8) +
  #   geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
  #                alpha = 0.8) +
  #   
  #   ggtitle(paste0("Mean mis-splicing ratio. ", cluster, "\nEach intron presents a mean of ", 
  #                  (round(samples %>% length() * 0.5) * threshold)  / samples %>% length()," reads in half of the samples (54 samples)")) +
  #   xlab("Mis-splicing ratio mean") +
  #   #facet_zoom(xlim(c(0,0.1))) +
  #   #ylab("Intron count") +
  #   #scale_y_continuous(limits = c(0, 85000)) +
  #   xlim(c(0, 1)) +
  #   #ylim(c(0, 17)) + 
  #   theme_light() +
  #   scale_fill_manual(values = c("#35B779FF","#440154FF"),
  #                     breaks = c("#35B779FF","#440154FF"),
  #                     labels = c("novel donor","novel acceptor")) +
  #   theme(axis.line = element_line(colour = "black"), 
  #         axis.text = element_text(colour = "black", size = "12"),
  #         axis.title = element_text(colour = "black", size = "12"),
  #         legend.text = element_text(size = "12"),
  #         legend.title = element_text(size = "12"),
  #         legend.position = "top") +
  #   guides(fill = guide_legend(title = NULL,
  #                              ncol = 2, 
  #                              nrow = 1))
  # 
  # ## Save the results
  # ggplot2::ggsave(filename = paste0(folder_save_name, "/images/", cluster, "_missplicingratio_",
  #                                   (round(samples %>% length() * 0.5) * threshold)  / samples %>% length(), "reads.png"), 
  #                 width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  
  
}

#############################
## ACROSS MULTIPLE TISSUES ##
#############################

## DISTANCES
summarise_distances_clusters <- function(clusters,
                                         folder_root,
                                         folder_results = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/") {
  
  
  df_distances_all <- map_df(clusters, function(cluster) {
    
    print(paste0(Sys.time(), " - loading ", cluster, " data..."))
    
    folder_name <- paste0(folder_root, "/", cluster, "/")
    df_tidy <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_tidy.rds"))
    
    df_donor_p <- df_tidy %>%
      filter(type == "novel_donor", distance > 0)
    
    df_donor_n <- df_tidy %>%
      filter(type == "novel_donor", distance < 0)
    
    df_acceptor_p <- df_tidy %>%
      filter(type == "novel_acceptor", distance > 0 )
    
    df_acceptor_n <- df_tidy %>%
      filter(type == "novel_acceptor", distance < 0 )
    
    
    return(data.frame(tissue = cluster,
                      
                      donor_mode_p = df_donor_p %>%
                        pull(distance) %>%
                        get_mode(),
                      
                      donor_mode_n = df_donor_n %>%
                        pull(distance) %>%
                        get_mode(),
                      
                      donor_median_p = df_donor_p %>%
                        pull(distance) %>%
                        median(),
                      
                      donor_median_n = df_donor_n %>%
                        pull(distance) %>%
                        median(),
                      
                      donor_mean_p = df_donor_p %>%
                        pull(distance) %>%
                        mean(),
                      
                      donor_mean_n = df_donor_n %>%
                        pull(distance) %>%
                        mean(),
                      
                      donor_min = df_donor_n %>%
                        select(distance) %>%
                        min(),
                      
                      donor_max = df_donor_p %>%
                        select(distance) %>%
                        max(),
                      
                      acceptor_mode_p = df_acceptor_p %>%
                        pull(distance) %>%
                        get_mode(),
                      
                      acceptor_mode_n = df_acceptor_n %>%
                        pull(distance) %>%
                        get_mode(),
                      
                      acceptor_median_p = df_acceptor_p %>%
                        pull(distance) %>%
                        median(),
                      
                      acceptor_median_n = df_acceptor_n %>%
                        pull(distance) %>%
                        median(),
                      
                      acceptor_mean_p = df_acceptor_p %>%
                        pull(distance) %>%
                        mean(),
                      
                      acceptor_mean_n = df_acceptor_n %>%
                        pull(distance) %>%
                        mean(),
                      
                      acceptor_min = df_acceptor_n %>%
                        pull(distance) %>%
                        min(),
                      
                      acceptor_max = df_acceptor_p %>%
                        pull(distance) %>%
                        max()))
    
  })
  
  
  saveRDS(object = df_distances_all,
          file = paste0(folder_results, "/df_distances_all.rds"))
  
}

## CLINVAR
summarise_clinvar_clusters <- function(clusters,
                                       folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/",
                                       folder_results = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/") {
  
  
  ## Load clinvar data
  clinvar_tidy <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/clinvar/clinvar-tidy.rda")
  clinvar_tidy %>% head()
  
  
  df_clinvar_all <- map_df(clusters, function(cluster) { # cluster = clusters[1]
    
    
    print(paste0(Sys.time(), " - loading ", cluster, " clinvar data..."))
    
    folder_name <- paste0(folder_root, "/", cluster, "/")
    df <- readRDS(file = paste0(folder_name, "/", cluster, "_missplicingratio_tidy_v2.rds")) %>%
      as_tibble() %>%
      filter(u2_intron == T) 
    
    
    ## Find overlaps between clinvar mutations and mis-splicing ratios
    overlaps <- GenomicRanges::findOverlaps(GenomicRanges::GRanges(seqnames = clinvar_tidy$chr,
                                                                   ranges = IRanges(start = clinvar_tidy$start_clinvar, 
                                                                                    end = clinvar_tidy$start_clinvar),
                                                                   strand = clinvar_tidy$strand),
                                            GenomicRanges::GRanges(seqnames = df$seqnames,
                                                                   ranges = IRanges(start = df$start, 
                                                                                    end = df$end),
                                                                   strand = df$strand),
                                            type = "within"
    )
    
    query <- S4Vectors::queryHits(overlaps)
    hits <- as(overlaps, "List")
    
    print(paste0(Sys.time(), " - obtaining clinvar mutations overlaps..."))
    
    df_clinvar <- map_df(unique(query), function(i) { # i <- 16
      
      clinvar_mut <- clinvar_tidy[i,]
      reference_junction <- df[hits[[i]],]
      
      return(data.frame(ref_junID = reference_junction$ref_junID,
                        missplicingratio_ND = reference_junction$missplicingratio_ND_tissue,
                        missplicingratio_NA = reference_junction$missplicingratio_NA_tissue,
                        clinvar_mut = clinvar_mut$hgvs,
                        closer_start_end = clinvar_mut$closer_start_end,
                        strand = clinvar_mut$strand,
                        type = reference_junction$type
                        
      ))
      
    })
    
    df_clinvar$clinvar <- 
      case_when(df_clinvar$closer_start_end == "start" & df_clinvar$strand == "+" ~ "acceptor", 
                df_clinvar$closer_start_end == "start" & df_clinvar$strand == "-" ~ "donor", 
                df_clinvar$closer_start_end == "end" & df_clinvar$strand == "+" ~ "donor", 
                df_clinvar$closer_start_end == "end" & df_clinvar$strand == "-" ~ "acceptor")
    
    
    
    ind <- which(df$ref_junID %in% df_clinvar$ref_junID)
    df_not_clinvar <- df[-ind,]
    
    
    ## WILCOXON TESTS
    
    donor_wilcoxon <- wilcox.test(x = df_clinvar %>% filter(clinvar == "donor") %>% pull(missplicingratio_ND),
                                  y = df_not_clinvar$missplicingratio_ND_tissue,
                                  paired = F,
                                  correct = T,
                                  alternative = "less")
    
    acceptor_wilcoxon <- wilcox.test(x = df_clinvar %>% filter(clinvar == "acceptor") %>% pull(missplicingratio_NA),
                                     y = df_not_clinvar$missplicingratio_NA_tissue,
                                     paired = F,
                                     correct = T,
                                     alternative = "less")
    
    
    return(data.frame(tissue = cluster,
                      donor_clinvar = which(df_clinvar$clinvar == "donor") %>% length(),
                      acceptor_clinvar = which(df_clinvar$clinvar == "acceptor") %>% length(),
                      donor_pvalue = donor_wilcoxon$p.value,
                      acceptor_pvalue = acceptor_wilcoxon$p.value))
    
    
  })
  
  
  
  df_clinvar_all
  df_clinvar_all %>% filter(tissue == "Muscle-Skeletal")
  
  saveRDS(object = df_clinvar_all,
          file = paste0(folder_results, "/df_clinvar_all.rds"))
  
}

## MODULO
summarise_modulo3_clusters <- function(clusters,
                                       folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/",
                                       folder_results = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/") {
  
  
  
  df_modulo_all <- map_df(clusters, function(cluster) { # cluster = clusters[1]
    
    
    print(paste0(Sys.time(), " - loading ", cluster, " modulo3 data..."))
    
    
    folder_name <- paste0(folder_root, "/", cluster, "/")
    
    df <- readRDS(file = paste0(folder_name, "/", cluster, "_missplicingratio_tidy_v2.rds")) %>%
      as_tibble()
    
    
    
    donorPC <- df %>%
      filter(missplicingratio_ND_tissue > 0, protein_coding == 100) %>%
      mutate(modulo = distance %% 3) 
    
    acceptorPC <- df %>%
      filter(missplicingratio_NA_tissue > 0, protein_coding == 100) %>%
      mutate(modulo = distance %% 3) 
    
    
    
    donorNPC <- df %>%
      filter(missplicingratio_ND_tissue > 0, protein_coding == 0) %>%
      mutate(modulo = distance %% 3)
    
    acceptorNPC <- df %>%
      filter(missplicingratio_NA_tissue > 0, protein_coding == 0) %>%
      mutate(modulo = distance %% 3)
    
    
    
    
    return(data.frame(tissue = cluster,
                      donorPC = ((donorPC %>% filter(modulo == 0) %>% nrow()) * 100) / donorPC %>% nrow(),
                      acceptorPC = ((acceptorPC %>% filter(modulo == 0) %>% nrow()) * 100) / acceptorPC %>% nrow(),
                      donorNPC = ((donorNPC %>% filter(modulo == 0) %>% nrow()) * 100) / donorNPC %>% nrow(),
                      acceptorNPC = ((acceptorNPC %>% filter(modulo == 0) %>% nrow()) * 100) / acceptorNPC %>% nrow()))
    
    
  })
  
  
  
  df_modulo_all
  df_modulo_all %>% filter(tissue == "Muscle-Skeletal")
  
  saveRDS(object = df_modulo_all,
          file = paste0(folder_results, "/df_modulo0_all.rds"))
  
}

## MIS-SPLICING RATIO
plot_missplicing_ratio_clusters <- function(clusters = c("Brain-FrontalCortex_BA9", "Liver", "WholeBlood"),
                                            folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/",
                                            folder_results = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/all_tissues/") {
  
  df <- data.frame(missplicing_donor = as.double(),
                   missplicing_acceptor = as.double(),
                   tissue = as.character())
  
  for (cluster in clusters) { # cluster <- clusters[1]
    
    folder_name <- paste0(folder_root, "/", cluster, "/")
    
    df <- rbind(df,
                readRDS(file = paste0(folder_name, "/", cluster, "_missplicingratio_tidy_v2.rds")) %>%
                  dplyr::select(missplicing_donor = missplicingratio_ND_tissue,
                                missplicing_acceptor = missplicingratio_NA_tissue) %>%
                  mutate(tissue = cluster))
    
    df %>% head()
    df %>% nrow()
    
  }
  
  ggplot(data = df) + 
    geom_density(aes(x = missplicing_donor, color = tissue)) +
    
    ggtitle("Comparison of mis-splicing ratio at donor splice sites") +
    xlab("Mis-splicing ratio mean") +
    xlim(c(0, 1)) +
    theme_light() +
    # scale_fill_manual(values = c("#35B779FF","#440154FF"),
    #                   breaks = c("#35B779FF","#440154FF"),
    #                   labels = c("novel donor","novel acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") 
  
  
  file_name <- paste0(folder_results, "/missplicing_donor_comparison.png")
  ggplot2::ggsave(filename = file_name, 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ggplot(data = df) + 
    geom_density(aes(x = missplicing_acceptor, color = tissue)) +
    
    ggtitle("Comparison of mis-splicing ratio at acceptor splice sites") +
    xlab("Mis-splicing ratio mean") +
    xlim(c(0, 1)) +
    theme_light() +
    # scale_fill_manual(values = c("#35B779FF","#440154FF"),
    #                   breaks = c("#35B779FF","#440154FF"),
    #                   labels = c("novel donor","novel acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") 
  file_name <- paste0(folder_results, "/missplicing_acceptor_comparison.png")
  ggplot2::ggsave(filename = file_name, 
                  width = 183, height = 183, units = "mm", dpi = 300)
  
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
# plot_distances <- function(clusters = clusters_ID,
#                            project_id,
#                            folder_root = paste0("PROJECTS/splicing-project/splicing-recount3-projects/", 
#                                                 project_id, "/results/pipeline3/missplicing-ratio/"),
#                            common = F,
#                            QC = F) {
#   
#   
#   title <- paste0("Distances to the reference junction")
#   
#   ## LOAD IDBs
#   
#   df_distances_clusters <- map_df(clusters, function(cluster) {
#     
#     print(paste0(Sys.time(), " - loading IDB for '", cluster, "' samples ..."))
#     readRDS(file = paste0(folder_root, "/", cluster, "/v105/", cluster, "_db_novel.rds")) %>%
#       mutate(sample_type = cluster) %>%
#       return()
#     
#   }) 
#   
#   
#   
#   
#   
#   ## GET common junctions
#   
#   if (common) {
#     
#     
#     common_junctions <- df_distances_clusters %>%
#       group_by(sample_type) %>%
#       distinct(ref_junID, .keep_all = T) %>%
#       ungroup() %>%
#       dplyr::count(ref_junID) %>%
#       filter(n == clusters %>% length()) %>%
#       pull(ref_junID)
#     
#     
#     
#     df_distances_tidy <- df_distances_clusters %>%
#       filter(ref_junID %in% common_junctions)
#     
#     
#     title <- paste0(title, "\n", common_junctions %>% unique() %>% length(), " common reference introns used.")
#     
#     # QC
#     # if (QC) {
#     #   
#     #   ## Each age group should have the same reference junction IDs
#     #   
#     #   for (cluster in df_clusters) {
#     #     
#     #     print(cluster)
#     #     
#     #     df_age_groups_tidy %>%
#     #       filter(sample_type == cluster) %>%
#     #       distinct(ref_junID) %>% 
#     #       nrow() %>% 
#     #       print()
#     #     
#     #     print("Novel junctions:")
#     #     
#     #     df_age_groups_tidy %>%
#     #       filter(sample_type == cluster) %>%
#     #       distinct(novel_junID) %>% 
#     #       nrow() %>% 
#     #       print()
#     #     
#     #   }
#     #   
#     #   
#     #   ## Second QC - novel junction IDs
#     #   
#     #   overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
#     #                                           subject = df_age_groups_tidy %>% filter(sample_type == age_groups[2]) %>% GenomicRanges::GRanges(),
#     #                                           type = "equal")
#     #   
#     #   if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
#     #                  (df_age_groups_tidy %>% filter(sample_type == age_groups[2]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
#     #     print("Massive error: some overlapping junctions don't have the same junction ID!")
#     #   }
#     #   
#     #   overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
#     #                                           subject = df_age_groups_tidy %>% filter(sample_type == age_groups[3]) %>% GenomicRanges::GRanges(),
#     #                                           type = "equal")
#     #   
#     #   if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
#     #                  (df_age_groups_tidy %>% filter(sample_type == age_groups[3]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
#     #     print("Massive error: some overlapping junctions don't have the same junction ID!")
#     #   }
#     #   
#     #   
#     #   
#     #   
#     # }
#   } else {
#     df_distances_tidy <- df_distances_clusters
#   }
#   
#   limit_bp <- 60
#   
#   
#   df_distances_tidy <- df_distances_tidy %>%
#     mutate(novel_type = factor(novel_type, levels = c("novel_donor", "novel_acceptor")))
#   
#   df_prop_noise <- df_distances_tidy %>%
#     filter(abs(distance) <= limit_bp) %>%
#     group_by(sample_type) %>%
#     mutate(N = n()) %>%
#     ungroup() %>%
#     group_by(sample_type, novel_type, distance) %>%
#     mutate(n_distance = n()) %>%
#     ungroup() %>%
#     mutate(p_noise = n_distance/N)%>%
#     as.data.frame()
#   
#   
#   df_prop_noise <- df_prop_noise %>%
#     mutate(sample_type = factor(sample_type, levels = clusters_ID))
#   
#   ggplot(data = df_prop_noise) + 
#     geom_col(aes(x = distance, y = p_noise, fill = sample_type),
#                    position = "dodge") +
#     facet_grid(vars(novel_type)) +
#     theme(legend.position = "top") 
#     
#   
#   
#   
#   ggplot(data = df_distances_tidy) + 
#     geom_histogram(aes(x = distance, fill = sample_type),
#                    bins = limit_bp * 2,
#                    binwidth = 1,
#                    position = "identity"
#     ) +
#     facet_grid(vars(novel_type)) +
#     ggtitle(title) +
#     xlab("Distance to the reference intron (in bp)") +
#     ylab("Number of unique novel junctions") +
#     theme_light() +
#     scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
#                        breaks = c((limit_bp * -1), (round(limit_bp / 2) * -1), 0, round(limit_bp / 2), limit_bp)) +
#     
#     # scale_fill_manual(values =  c("#FDE725FF", "#21908CFF", "#440154FF"),
#     #                   breaks = c("20-39", "40-59", "60-79")) +
#     
#     # scale_fill_manual(breaks = c("20-39", "40-59", "60-79"),
#     #                   labels = c("20-39", "40-59", "60-79")) +
#     guides(fill = guide_legend(title = NULL, 
#                                ncol = clusters %>% length(), nrow = 1 )) +
#     theme(axis.line = element_line(colour = "black"), 
#           axis.text = element_text(colour = "black", size = "14"),
#           axis.title = element_text(colour = "black", size = "14"),
#           strip.text = element_text(colour = "black", size = "16"), 
#           legend.text = element_text(colour = "black", size = "14"),
#           plot.caption = element_text(colour = "black", size = "14"),
#           plot.title = element_text(colour = "black", size = "16"),
#           legend.title = element_text(colour = "black", size = "14"),
#           legend.position = "top") %>% 
#     return()
#   
#   
# }
# 
# 
# plot_MSR <- function(clusters = clusters_ID,
#                      project_id,
#                      folder_root = paste0("PROJECTS/splicing-project/splicing-recount2-projects/", 
#                                           project_id, "/results/pipeline3/missplicing-ratio/"),
#                      common = F,
#                      QC = F) {
#   
#   
#   title <- paste0("MSR")
#   
#   ## LOAD IDBs
#   
#   df_clusters <- map_df(clusters, function(cluster) {
#     
#     #print(paste0(Sys.time(), " - loading IDB for '", age_group, "' samples ..."))
#     readRDS(file = paste0(folder_root, "/", cluster, "/v105/", cluster, "_db_introns.rds")) %>%
#       mutate(sample_type = cluster) %>%
#       return()
#     
#   }) 
#   
#   
#   
#   
#   
#   ## GET common junctions
#   
#   if (common) {
#     
#     
#     common_junctions <- df_clusters %>%
#       group_by(sample_type) %>%
#       distinct(ref_junID, .keep_all = T) %>%
#       ungroup() %>%
#       dplyr::count(ref_junID) %>%
#       filter(n == clusters %>% length()) %>%
#       pull(ref_junID)
#     
#     
#     
#     df_clusters_tidy <- df_clusters %>%
#       filter(ref_junID %in% common_junctions)
#     
#     
#     title <- paste0(title, "\n", common_junctions %>% unique() %>% length(), " common reference introns used.")
#     
#     # QC
#     # if (QC) {
#     #   
#     #   ## Each age group should have the same reference junction IDs
#     #   
#     #   for (cluster in df_clusters) {
#     #     
#     #     print(cluster)
#     #     
#     #     df_age_groups_tidy %>%
#     #       filter(sample_type == cluster) %>%
#     #       distinct(ref_junID) %>% 
#     #       nrow() %>% 
#     #       print()
#     #     
#     #     print("Novel junctions:")
#     #     
#     #     df_age_groups_tidy %>%
#     #       filter(sample_type == cluster) %>%
#     #       distinct(novel_junID) %>% 
#     #       nrow() %>% 
#     #       print()
#     #     
#     #   }
#     #   
#     #   
#     #   ## Second QC - novel junction IDs
#     #   
#     #   overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
#     #                                           subject = df_age_groups_tidy %>% filter(sample_type == age_groups[2]) %>% GenomicRanges::GRanges(),
#     #                                           type = "equal")
#     #   
#     #   if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
#     #                  (df_age_groups_tidy %>% filter(sample_type == age_groups[2]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
#     #     print("Massive error: some overlapping junctions don't have the same junction ID!")
#     #   }
#     #   
#     #   overlaps <- GenomicRanges::findOverlaps(query = df_age_groups_tidy %>% filter(sample_type == age_groups[1]) %>% GenomicRanges::GRanges(),
#     #                                           subject = df_age_groups_tidy %>% filter(sample_type == age_groups[3]) %>% GenomicRanges::GRanges(),
#     #                                           type = "equal")
#     #   
#     #   if (!identical((df_age_groups_tidy %>% filter(sample_type == age_groups[1]))[S4Vectors::queryHits(overlaps),]$novel_junID,
#     #                  (df_age_groups_tidy %>% filter(sample_type == age_groups[3]))[S4Vectors::subjectHits(overlaps),]$novel_junID)){
#     #     print("Massive error: some overlapping junctions don't have the same junction ID!")
#     #   }
#     #   
#     #   
#     #   
#     #   
#     # }
#   } else {
#     df_clusters_tidy <- df_clusters
#   }
#   
#   
#   
#   
#   df_clusters_tidy <- df_clusters_tidy %>%
#     mutate(novel_type = factor(novel_type, levels = c("novel_donor", "novel_acceptor")))
#   
#   
#   limit_bp <- 60
#   
#   ggplot(data = df_clusters_tidy) + 
#     geom_boxplot(aes(y = ref_missplicing_ratio_tissue_ND, fill = sample_type)) 
#   ggplot(data = df_clusters_tidy) + 
#     geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = ref_type)) 
#    
#   
# }