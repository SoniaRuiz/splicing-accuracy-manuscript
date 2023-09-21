

#' Title
#' Creates the child tables for each of the recount3 ID projects indicated
#' @param database.path Local path to the .sqlite database
#' @param recount3.project.IDs Vector with the ID of the recount3 projects to work with
#' @param database.folder Local path to the folder that contains the database
#' @param results.folder Local path to the folder that contains the result files
#'
#' @return
#' @export
#'
#' @examples
sql_create_child_tables <- function(database.path,
                                    recount3.project.IDs = NULL,
                                    database.folder,
                                    results.folder,
                                    supportive.reads) {
  
  all_split_reads_details_105 <- readRDS(file = paste0(database.folder, "/all_split_reads_qc_level2.rds"))
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  tables <- DBI::dbListTables(conn = con)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  print( paste0(Sys.time(), " --> SQL connection stablished!") )
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'metadata'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## GET FROM INTRON TABLE
  query = paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query) %>% as_tibble()
  df_intron %>% nrow()
  df_intron %>% 
    dplyr::count(misspliced)
  
  ## GET FROM NOVEL JUNCTION TABLE
  query = paste0("SELECT * FROM 'novel'")
  df_novel <- dbGetQuery(con, query) %>% as_tibble()
  df_novel %>% nrow()
  df_novel %>% 
    dplyr::count(novel_type) %>%
    print()
  
  ## GET FROM GENE TABLE
  query = paste0("SELECT * FROM 'gene'")
  df_gene <- dbGetQuery(con, query) %>% as_tibble()
  df_gene %>% nrow()
  
  ## GET FROM TRANSCRIPT TABLE
  query = paste0("SELECT * FROM 'transcript'")
  df_transcript <- dbGetQuery(con, query) %>% as_tibble()
  df_transcript %>% nrow()
  
  DBI::dbDisconnect(conn = con) 
  
  if ( is.null(recount3.project.IDs) ){
    recount3.project.IDs <- (df_metadata$SRA_project %>% unique())
  }
  
  ## Loop through the projects parallely
  # doParallel::registerDoParallel(10)
  # foreach(i = seq(length(recount3.project.IDs))) %dopar%{
  #   
  #   project_id <- recount3.project.IDs[i]
  
  for (project_id in recount3.project.IDs) { 
    
    # project_id <- recount3.project.IDs[1]
    
    message(Sys.time(), " --> Working with '", project_id, "' ...")
    results_folder_local <- paste0(results.folder, "/", project_id, "/")
    
    clusters <- df_metadata %>%
      dplyr::filter(SRA_project == project_id) %>%
      distinct(cluster) %>%
      pull()
    
    for (cluster_id in clusters) { 
      
      # cluster_id <- clusters[1]
      message(Sys.time(), " --> ", cluster_id)
      
      ###############################
      ## PREPARE DATA
      ###############################
    
      if ( file.exists(paste0(results_folder_local, "/junction_pairing/", 
                              cluster_id, "/", cluster_id, "_raw_distances_tidy.rds")) ) {
        
        
        ## LOAD BASE DATA ONLY FOR THE CURRENT CLUSTER ID -------------------------------------------------
        
        message(Sys.time(), " --> ", cluster_id, " load base data ... ")
        
        ## Load split read counts
        split_read_counts <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                                   project_id, "_", cluster_id, "_split_read_counts.rds")) 
        print(paste0(Sys.time(), " --> ", cluster_id, " split read counts loaded!"))
        stopifnot(
          "Still there are split reads with less than N number of supportive reads" =
            split_read_counts %>% 
            mutate(sumCounts = rowSums(select(., !contains("junID")))) %>%
            filter(sumCounts < supportive.reads) %>% 
            nrow() == 0
        )
        
        
        if ( is.null(names(split_read_counts)) ) {
          split_read_counts <- split_read_counts %>%
            as_tibble(rownames = "junID")
        }
        
        ## Load samples
        samples <- readRDS(file = paste0(results_folder_local, "/base_data/", 
                                         project_id, "_", cluster_id, "_samples_used.rds"))
        
        
        ## LOAD INTRONS AND NOVEL JUNCTIONS ------------------------------------
        df_cluster_distances <- readRDS(file = paste0(results_folder_local, "/junction_pairing/", cluster_id, "/", 
                                                      cluster_id, "_raw_distances_tidy.rds")) %>% as_tibble()
        
        
        
        ## INTRONS ---------------------------------------------------
        df_introns_gr <- df_cluster_distances %>%
          distinct(ref_junID, .keep_all = T) %>%
          dplyr::select(ref_junID, seqnames = ref_seq, start = ref_start,
                        end = ref_end, strand = ref_strand)
        
        
        
        ## Add reads detected for the introns in the current tissue
        ## It might be the case that, after subdividing the junctions across the samples of a given cluster, that
        ## a junction only presents a read within that cluster. However, that junction for sure has at least two supportive reads across all samples from that cluster.
        split_read_counts_intron <- generate_coverage(split_read_counts = split_read_counts,
                                                      samples = samples,
                                                      junIDs = df_introns_gr$ref_junID) %>%
          dplyr::rename(ref_n_individuals = n_individuals,
                        ref_sum_counts = sum_counts)
        
        split_read_counts_intron %>% head()
        
        # if ( split_read_counts_intron$ref_sum_counts %>% unique() %>% min() == 1 ) {
        #   message(cluster_id, " - ERROR! Minimum number of supporting reads is 1!")
        #   break;
        # }
        
        df_introns_merged <- df_introns_gr %>%
          left_join(y = split_read_counts_intron,
                    by = c("ref_junID" = "junID"))
        
        
        ## NOVEL JUNCTIONS -----------------------------------------------------
        df_novel_gr <- df_cluster_distances %>%
          distinct(novel_junID, .keep_all = T) %>%
          dplyr::select(novel_junID, ref_junID, seqnames = novel_seq, 
                        start = novel_start, end = novel_end, strand = novel_strand)
        
        
        split_read_counts_novel <- generate_coverage(split_read_counts = split_read_counts,
                                                     samples = samples,
                                                     junID = df_novel_gr$novel_junID) %>%
          dplyr::rename(novel_n_individuals = n_individuals,
                        novel_sum_counts = sum_counts)
        
        
        df_novel_merged <- df_novel_gr %>%
          left_join(y = split_read_counts_novel,
                    by = c("novel_junID" = "junID"))
        
        df_novel_merged %>% as_tibble()
        
        
        ## QC data - NO '*' within the ID ---------------------------------------------
        
        df_novel_merged <- df_novel_merged %>% 
          rowwise() %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"), 
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID)) 
        
        df_novel_merged <- df_novel_merged %>% 
          rowwise() %>%
          mutate(novel_junID = ifelse(str_detect(string = novel_junID, pattern = "\\*"),
                                      str_replace(string = novel_junID, pattern = "\\*", strand ),
                                      novel_junID)) 
        
        df_introns_merged <- df_introns_merged %>% 
          rowwise() %>%
          mutate(strand = strand %>% as.character()) %>%
          mutate(ref_junID = ifelse(str_detect(string = ref_junID, pattern = "\\*"),
                                    str_replace(string = ref_junID, pattern = "\\*", strand ),
                                    ref_junID))
        
        
        ## Merge tidy data -----------------------------------------------------
        
        print(paste0(Sys.time(), " --> merge introns and novel junctions ... "))
        
        df_intron_novel_merged <- df_novel_merged %>%
          dplyr::select(-c(seqnames, start, end, strand)) %>%
          inner_join(y = df_introns_merged %>% 
                       dplyr::select(-c(seqnames, start, end, strand)),
                     by = "ref_junID")
        
        
        
        
        ####################################
        ## QC
        ####################################
        
        if (any(str_detect(df_intron_novel_merged$ref_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(str_detect(df_intron_novel_merged$novel_junID,pattern = "//*"))) {
          print("ERROR: some IDs contain '*'")
        }
        if (any(df_intron_novel_merged$ref_junID %>% is.na())) {
          print("ERROR: There are missing reference introns!!")
        }
        
        
        
        ## JOIN data with MASTER NOVEL table
        ## Only novel junctions previously paired and QC'ed are considered
        
        print(paste0(Sys.time(), " --> merge data with MASTER NOVEL table ... "))
        
        df_all_misspliced <- df_intron_novel_merged %>% 
          inner_join(y = df_novel %>% 
                       dplyr::select(novel_junID, novel_coordinates, novel_type),
                     by = c("novel_junID" = "novel_coordinates")) %>%
          dplyr::rename(novel_coordinates = novel_junID) %>%
          dplyr::rename(novel_junID = novel_junID.y)
        
        
        if (setdiff(df_all_misspliced$ref_junID, df_intron$ref_coordinates) %>% length() > 0) {
          print("ERROR! Some introns detected in this tissue are not stored in the master intron table.")
          break;
        }
        
        ## JOIN data with MASTER INTRON table
        df_all_misspliced <- df_all_misspliced %>% 
          left_join(y = df_intron %>%
                      dplyr::filter(misspliced == T) %>%
                      dplyr::select(ref_junID, ref_coordinates, transcript_id),
                    by = c("ref_junID" = "ref_coordinates")) %>%
          dplyr::select(-ref_junID) %>%
          dplyr::rename(ref_junID = ref_junID.y)
        
        
        if (which(str_detect(df_all_misspliced$novel_coordinates,pattern = "//*")) %>% length() > 0) {
          print("ERROR: some IDs contain '*'")
        }
        
        
        ## PREPARE DATA PRIOR POPULATING THE TABLE
        df_all_misspliced <- df_all_misspliced %>%
          relocate(ref_junID, novel_junID)
        
        
        
        #######################################
        ## CHECK INTEGRITY WITH PARENT TABLES
        #######################################
        
        print(paste0(Sys.time(), " --> checking integrity with parent tables ... "))
        
        master_novel <- df_novel %>%
          dplyr::filter(novel_junID %in% 
                          (df_all_misspliced %>%
                             pull(novel_junID))) %>% 
          dplyr::select(novel_coordinates) 
        
        if ( !(identical(df_all_misspliced$novel_coordinates %>% sort(), 
                         master_novel$novel_coordinates %>% sort())) ) {
          print("ERROR! Tables not identical")
          if ( all(intersect(setdiff(df_all_misspliced$novel_coordinates, 
                                     master_novel$novel_coordinates),
                             ambiguous_novel_junc$novel_junID) == setdiff(df_all_misspliced$novel_coordinates, 
                                                                          master_novel$novel_coordinates)) == T ) {
            df_all_misspliced <- df_all_misspliced %>%
              as.data.table() %>%
              dplyr::filter(!(novel_coordinates %in% ambiguous_novel_junc$novel_junID))
          }
        }
        
        if ( identical(df_all_misspliced$novel_coordinates %>% sort(), 
                       master_novel$novel_coordinates %>% sort()) ) {
          
          df_all_misspliced <- df_all_misspliced %>%
            dplyr::select(-novel_coordinates)
          
          
          
          
          ## CHECK INTEGRITY WITH PARENT TABLE
          
          df <- df_novel %>%
            dplyr::select(novel_junID, ref_junID) %>%
            arrange(novel_junID) %>% 
            inner_join(df_all_misspliced %>%
                         dplyr::select(novel_junID, ref_junID) %>%
                         arrange(novel_junID),
                       by = "novel_junID")
          
          diff <- df %>% 
            dplyr::filter(ref_junID.x != ref_junID.y)
          
          
          if (diff %>% nrow() > 0) {
            df_all_misspliced <- df_all_misspliced %>%
              dplyr::filter(!(ref_junID %in% diff$ref_junID.y)) 
            print(paste0("ERROR!: ", diff, " --> mismatch junctions."))
            break;
          } 
          
          
          #####################################
          ## CALCULATE MSR MEASURES
          #####################################
          
          print(paste0(Sys.time(), " --> calculating MSR measures ... "))
          
          # df_all_misspliced %>% dplyr::filter(ref_junID == "17") %>% as.data.frame()
          
          db_introns <- df_all_misspliced %>%
            group_by(ref_junID, novel_type) %>%
            mutate(MSR = sum(novel_sum_counts)/(sum(novel_sum_counts) + ref_sum_counts)) %>%
            ungroup()
          
          # db_introns %>% dplyr::filter(ref_junID == "17") %>% as.data.frame()
          
          db_introns <- db_introns %>% 
            spread(key = novel_type, value = MSR, fill = 0) %>%
            group_by(ref_junID) %>% 
            ## If the annotated intron is mis-spliced at both splice sites, it will have some rows MSR=0 and some rows MSR>0 at the 5' 
            ## and at the 3'ss. We select the highest value
            dplyr::mutate(MSR_Donor = max(novel_donor, na.rm = T)) %>%
            dplyr::mutate(MSR_Acceptor = max(novel_acceptor, na.rm = T)) %>%
            ungroup()
          
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
          if ( file.exists( paste0(results.folder, "/", project_id, "/tpm/", project_id, "_", cluster_id, "_tpm.rds")) ) {
            
            tpm <- readRDS(file = paste0(results.folder,  "/", project_id, "/tpm/", project_id, "_", cluster_id, "_tpm.rds")) %>% 
              dplyr::select(gene_id, all_of(samples))
            
            tpm <- tpm  %>%
              dplyr::mutate(tpm_median = matrixStats::rowMedians(x = as.matrix(.[2:(ncol(tpm))]))) %>%
              dplyr::select(gene_id, tpm_median) 
            
            ## In case there are any duplicates, take the genes with the maximum tpms
            tpm <- tpm %>% 
              distinct(gene_id, .keep_all = T) %>%
              group_by(gene_id) %>% 
              summarize_all(max) %>%
              ungroup()
            
            
            df_gene %>%
              filter(gene_id == "ENSG00000000419")
            
            
            tpm_tidy <- tpm %>%
              inner_join(y = df_gene %>% as_tibble(),
                         by = c("gene_id" = "gene_id")) %>%
              inner_join(y = df_transcript %>% as_tibble(),
                         by = c("id" = "gene_id"),
                         multiple = "all") %>%
              dplyr::select(transcript_id = id.y, 
                            tpm_median)
            
            db_introns <- db_introns %>%
              left_join(y = tpm_tidy,
                        by = c("transcript_id" = "transcript_id")) %>% 
              dplyr::rename(gene_tpm = tpm_median)
            
          }
          
          
          #####################################
          ## ADD THE TYPE OF INTRON
          ## i.e. mis-spliced at the donor, at the 
          ## acceptor or at both splice sites
          #####################################
          
          db_introns[, "ref_type"] <- ""
          db_introns %>% head()
          
          ## TYPE 'BOTH'
          indx <- which(db_introns$MSR_Donor > 0 & db_introns$MSR_Acceptor > 0)
          db_introns[indx, "ref_type"] <- "both"
          print(paste0("Introns type 'both': ", indx %>% length()))
          
          
          ## TYPE 'DONOR'
          indx <- which(db_introns$MSR_Donor > 0 & db_introns$MSR_Acceptor == 0)
          db_introns[indx, "ref_type"] <- "donor"
          print(paste0("Junctions type 'donor': ", indx %>% length()))
          
          
          ## TYPE 'ACCEPTOR'
          indx <- which(db_introns$MSR_Donor == 0 & db_introns$MSR_Acceptor > 0)
          db_introns[indx, "ref_type"] <- "acceptor"
          print(paste0("Junctions type 'acceptor': ", indx %>% length()))
          
        
          db_introns %>% nrow()
          db_introns %>% distinct(ref_junID) %>% nrow()
          
          
          #########################################################
          ## CREATE AND POPULATE CHILD 'MIS-SPLICED' INTRON TABLE
          #########################################################
          
          print(paste0(Sys.time(), " --> creating 'intron' table ... "))
          
          # dbRemoveTable(conn = con, paste0(cluster_id, "_", project_id))
          query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster_id, "_", project_id, "_misspliced"), "'", 
                          "(ref_junID INTEGER NOT NULL,
                          novel_junID INTEGER NOT NULL,
                          ref_n_individuals INTEGER NOT NULL,
                          ref_sum_counts INTEGER NOT NULL,
                          ref_type TEXT NOT NULL, 
                          novel_n_individuals INTEGER NOT NULL, 
                          novel_sum_counts INTEGER NOT NULL, 
                          MSR_D DOUBLE NOT NULL, 
                          MSR_A DOUBLE NOT NULL, 
                          gene_tpm DOUBLE, 
                          transcript_id INTEGER NOT NULL, 
                          FOREIGN KEY (ref_junID, novel_junID) REFERENCES novel (ref_junID, novel_junID),
                          FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
          
          ## Connect the database
          con <- dbConnect(RSQLite::SQLite(), database.path)
          DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
          ## Create the child table
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          ## POPULATE THE TABLE
          db_introns_final <- db_introns %>%
            dplyr::select(-novel_acceptor,-novel_donor)%>%
            dplyr::rename(MSR_D = MSR_Donor, MSR_A = MSR_Acceptor)
          
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster_id, "_", project_id, "_misspliced"), 
                             value = db_introns_final)
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", paste0(cluster_id, "_", project_id, "_misspliced"), "' ON '",
                          paste0(cluster_id, "_", project_id, "_misspliced"), "'(ref_junID,novel_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          ## Disconnect the database
          DBI::dbDisconnect(conn = con)
          
          
          print(paste0(Sys.time(), ". '", paste0(cluster_id, "_", project_id, "_misspliced"), 
                       "' table populated!"))
          
          ####################################
          ## START WITH NEVER MISSPLICED INTRONS
          ####################################
          
          print(paste0(Sys.time(), " --> getting never mis-spliced introns ... "))
          
          ## TYPE 'NONE'
          introns_never <- readRDS(file = paste0(results_folder_local, "/junction_pairing/", cluster_id, "/not-misspliced/", 
                                                 cluster_id, "_all_notmisspliced.rds"))
          ## The introns not misspliced in this tissue, should have not been detected as spliced.
          ## Thus, this should be zero
          if (intersect(introns_never, df_intron %>%
                        dplyr::filter(ref_junID %in% db_introns_final$ref_junID) %>%
                        pull(ref_coordinates) %>% 
                        unique()) %>% length() > 0 ){
            print("ERROR! The introns not misspliced in this tissue, should have not been detected as spliced.")
            break;
          }
          df_introns_never <- data.frame(ref_junID = introns_never)
          
          
          
          if (any(str_detect(string = df_introns_never$ref_junID, pattern = "\\*"))) {
            print("ERROR! * in the IDs")
            break;
          }
          
          # split_read_counts <- split_read_counts %>%
          #   left_join(y = all_split_reads_details_105 %>% 
          #               dplyr::select(junID, seqnames, start, end, strand),
          #             by = "junID") %>%
          #   dplyr::select(-junID) %>%
          #   mutate(junID = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>%
          #   dplyr::relocate(junID)
          
          split_read_counts_intron_never <- generate_coverage(split_read_counts = split_read_counts,
                                                              samples = samples,
                                                              junID = df_introns_never$ref_junID) %>%
            dplyr::rename(ref_n_individuals = n_individuals,
                          ref_sum_counts = sum_counts) 
          
          if (any(str_detect(split_read_counts_intron_never$junID,pattern = "\\*"))) {
            print("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          
          
          df_never_merged <- df_introns_never %>%
            inner_join(y = split_read_counts_intron_never,
                       by = c("ref_junID" = "junID")) %>%
            mutate(MSR_D = 0, MSR_A = 0) %>%
            mutate(ref_type = "never")
          df_never_merged <- df_never_merged %>% as_tibble()
          
          if ( any(df_never_merged$ref_n_individuals %>% is.na()) ) {
            print("ERROR! some never mis-spliced junctions without the number of individuals")
            break;
          }
          
          ## QC
          if (any(str_detect(df_never_merged$ref_junID, pattern = "\\*"))) {
            print("ERROR! * in the IDs")
            break;
          }
          
          ## TYPE 'NONE'
          print(paste0("Junctions type 'never': ", df_never_merged %>% distinct(ref_junID) %>% nrow()))
          
          ##################################################
          ## ADD REFERENCE KEY TO THE MASTER INTRON TABLE 
          ##################################################
          
          print(paste0(Sys.time(), " --> adding the intron reference key to the never mis-spliced introns ... "))
          
          df_never_merged <- df_never_merged %>%
            inner_join(df_intron %>% 
                         dplyr::select(ref_junID, ref_coordinates, transcript_id),
                       by = c("ref_junID" = "ref_coordinates")) %>%
            dplyr::filter(!is.na(ref_junID)) %>%
            dplyr::select(-ref_junID) %>% 
            dplyr::rename(ref_junID = ref_junID.y) %>%
            relocate(ref_junID)
          
          
          if (df_never_merged %>% dplyr::filter(is.na(ref_junID)) %>% nrow > 0) {
            print("ERROR! IDs are NA")
            break;
          }
          if (any(df_never_merged$ref_type != "never")) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
            break;
          }
          if ((intersect(df_never_merged$ref_junID, db_introns$ref_junID) %>% length()) > 0) {
            print("Error! Some never mis-spliced junctions have been stored as mis-spliced.")
            break;
          }
          if (any(duplicated(df_never_merged$ref_junID))) {
            print("Error! Some never mis-spliced junctions are duplicated.")
            break;
          }
          
          ## Remove introns and novel junctions with not TPM data
          # df_never <- df_never %>%
          #   dplyr::filter(!(gene_tpm %>% is.na()))    
          
          #####################################
          ## GET THE GENE TPM
          #####################################
          
          
          if (exists("tpm_tidy")) {
            df_never_merged <- df_never_merged %>%
              left_join(y = tpm_tidy,
                        by = "transcript_id") %>% 
              dplyr::rename(gene_tpm = tpm_median) 
          }
          
          ## QC --------------------------------------------------------
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "donor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
            break;
          }
          
          
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "aceptor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
            break;
          }
          
          if (intersect(df_never_merged %>%
                        dplyr::filter(ref_type == "never") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "both") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
            break;
          }
          
          if (intersect(db_introns %>%
                        dplyr::filter(ref_type == "acceptor") %>%
                        distinct(ref_junID) %>%
                        pull(), db_introns %>%
                        dplyr::filter(ref_type == "donor") %>%
                        distinct(ref_junID) %>%
                        pull()) %>% length() > 0) {
            print("Error!")
            break;
          }
          
          
          #############################################################
          ## CREATE AND POPULATE CHILD 'NEVER MIS-SPLICED' INTRON TABLE
          #############################################################
          
          
          query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster_id, "_", project_id, "_nevermisspliced"), "'", 
                          "(ref_junID INTEGER NOT NULL,
                          ref_n_individuals INTEGER NOT NULL, 
                          ref_sum_counts INTEGER NOT NULL,
                          MSR_D DOUBLE NOT NULL, 
                          MSR_A DOUBLE NOT NULL, 
                          ref_type TEXT NOT NULL, 
                          gene_tpm DOUBLE,
                          transcript_id INTEGER NOT NULL,
                          FOREIGN KEY (ref_junID) REFERENCES intron (ref_junID),
                          FOREIGN KEY (transcript_id) REFERENCES 'transcript'(id))")
          
          
          ## Connect the database
          con <- dbConnect(RSQLite::SQLite(), database.path)
          DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
          
          ## Create the child table 
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          ## POPULATE TABLE
          DBI::dbAppendTable(conn = con,
                             name = paste0(cluster_id, "_", project_id, "_nevermisspliced"), 
                             value = df_never_merged)
          
          ## CREATE INDEX
          query <- paste0("CREATE UNIQUE INDEX 'index_", 
                          paste0(cluster_id, "_", project_id, "_nevermisspliced"), "' ON '",
                          paste0(cluster_id, "_", project_id, "_nevermisspliced"),"'(ref_junID)");
          res <- DBI::dbSendQuery(conn = con, statement = query)
          DBI::dbClearResult(res)
          
          ## Disconnect the database
          DBI::dbDisconnect(conn = con)
          
          print(paste0(Sys.time(), ". '", 
                       paste0(cluster_id, "_", project_id, "_nevermisspliced"), 
                       "' table populated!"))
          
        } else {
          print("Error: novel junctions are distinct!")
          break;
        }
       
      }
    }
    gc()
  }
}
