
####################################################
## GENERATION OF THE INTRON DATABASE ###############
####################################################


## DISTANCES ---------------------------------------

get_distances <- function(cluster, 
                          samples,
                          # recount.info,
                          split_read_counts,
                          all_split_reads_details,
                          folder_name) {
  
  ## DECLARE SOME VARIABLES
  
  num_sample <- 1
  
  ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
  for (sample in samples) { 
    
    # sample <- samples[1]
    # sample <- samples[2]
    # sample <- "54285"

    if ( !file.exists(paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")) ) {
      
      if ( length(which((colnames(split_read_counts) == sample) == TRUE)) == 1 ) { 
        
        # Every sampleID is unique, so the result should be equal to 1
        
        print(cluster)
        print(paste0(Sys.time(), " - sample '", num_sample, "'"))
        
        # indx <- which(recount.info$sample == sample)
        # sample_mapped_read_count <- recount.info[indx,]$mapped_read_count
        ############################################################
        ## OBTAINING ALL JUNCTIONS AND COUNTS FOR THE CURRENT SAMPLE
        ############################################################
        
        
        split_read_counts_sample <- split_read_counts %>%
          as_tibble() %>%
          dplyr::select(junID, all_of(sample %>% as.character())) %>%
          drop_na() 
        
        split_read_counts_sample <- split_read_counts_sample[(split_read_counts_sample[,2] >0),]
          
          
        split_read_counts %>% nrow()
        split_read_counts_sample %>% nrow()
        split_read_counts_sample %>% head()
        
        split_reads_details_97_sample <- merge(all_split_reads_details, 
                                               split_read_counts_sample, 
                                               by = "junID", 
                                               sort = TRUE)  %>%
          dplyr::rename(counts = all_of(sample %>% as.character())) %>%
          as_tibble()
        
        split_reads_details_97_sample %>% head()
        
        ## Do some QC
        index <- runif(n = 1, 1, split_reads_details_97_sample %>% nrow()) %>% as.integer()
        data <- split_reads_details_97_sample[index, c(1, split_reads_details_97_sample %>% ncol())]
        data_col <- split_read_counts_sample %>%
          dplyr::filter(junID == data$junID)
        if (data_col[, colnames(data_col) == sample] != data$counts) {
          print("Error: QC failed!")
          break;
        }
        
        ###################################################
        ##########    ANNOTATED JUNC   ####################
        ###################################################
        all_annotated <- split_reads_details_97_sample %>%
          dplyr::filter(type == "annotated") %>% 
          GenomicRanges::GRanges()
        
        all_annotated_forward <- all_annotated[all_annotated@strand == "+"]
        all_annotated_reverse <- all_annotated[all_annotated@strand == "-"]
        
        
        
        ###################################################
        ##########    NOVEL DONOR      ####################
        ###################################################
        all_donor <- split_reads_details_97_sample %>%
          dplyr::filter(type == "novel_donor") %>% 
          GRanges()
        
        all_donor_forward <- all_donor[all_donor@strand == "+"]
        all_donor_reverse <- all_donor[all_donor@strand == "-"]
        
        
        
        ###################################################
        ##########    NOVEL ACCEPTOR      #################
        ###################################################
        all_acceptor <- split_reads_details_97_sample%>%
          dplyr::filter(type == "novel_acceptor") %>% 
          GRanges()
        
        all_acceptor_forward <- all_acceptor[all_acceptor@strand == "+"]
        all_acceptor_reverse <- all_acceptor[all_acceptor@strand == "-"]
        
        
        
        ####################################################
        ########  NOVEL-DONOR & FORWARD STRAND  ############
        ####################################################
        
        overlaps <- GenomicRanges::findOverlaps(query = all_donor_forward,
                                                subject = all_annotated_forward,
                                                ignore.strand = FALSE,
                                                type = "end")
        
        novel_junctions <- all_donor_forward[queryHits(overlaps),]
        ref_junctions <- all_annotated_forward[subjectHits(overlaps),]
        
        
        distance <- start(ref_junctions) - start(novel_junctions)
        
        # distance <- ifelse(distance < 0, distance + 1, distance - 1)
        if (distance %>% length() > 0) {
        
          df_nd_f <- data.frame(sample = sample,
                                type = "novel_donor", 
                                distance = distance,
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                #novel_counts_norm = novel_junctions$counts / sample_mapped_read_count,
                                novel_seq = novel_junctions %>% seqnames() %>% as.character(),
                                novel_start = novel_junctions %>% start() %>% as.character(),
                                novel_end = novel_junctions %>% end() %>% as.character(),
                                novel_strand = novel_junctions %>% strand() %>% as.character(),
                                novel_width = novel_junctions %>% width() %>% as.character(),
                                #novel_ss5score = novel_junctions$ss5score,
                                #novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_end = ref_junctions %>% end(),
                                ref_width = ref_junctions %>% width(),
                                #ref_ss5score = ref_junctions$ss5score,
                                #ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          ## In case the novel junction has been attached to many ref junctions, the ref junction chosen will be 
          ## the one with the highest number of counts. In case many ref junctions present the same number of counts,
          ## the closest one will be chosen.
          
          df_nd_f <- df_nd_f %>% 
            group_by(novel_junID) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            ungroup() 
        
        } else {
          df_nd_f <- NULL
          
        }
        
        # novel <- GenomicRanges::GRanges(seqnames = df_nd_f[1,]$novel_seq,
        #                                 ranges = IRanges(start = df_nd_f[1,]$novel_start %>% as.integer(),
        #                                                  end = df_nd_f[1,]$novel_start %>% as.integer()),
        #                                 strand = df_nd_f[1,]$novel_strand)
        # ref <-  GenomicRanges::GRanges(seqnames = df_nd_f[1,]$ref_seq,
        #                                ranges = IRanges(start = df_nd_f[1,]$ref_start %>% as.integer(),
        #                                                 end = df_nd_f[1,]$ref_start %>% as.integer()),
        #                                strand = df_nd_f[1,]$ref_strand)
        # GenomicRanges::distance(x = novel, y = ref)
        
        # df_nd_f %>%
        #   dplyr::select(novel_strand, novel_start, novel_end, ref_strand, ref_start, ref_end, distance)
        
        ####################################################
        ########  NOVEL-DONOR & REVERSE STRAND  ############
        ####################################################
        
        overlaps <- GenomicRanges::findOverlaps(query = all_donor_reverse,
                                                subject = all_annotated_reverse,
                                                ignore.strand = FALSE,
                                                type = "start")
        
        
        
        novel_junctions <- all_donor_reverse[queryHits(overlaps),]
        ref_junctions <- all_annotated_reverse[subjectHits(overlaps),]
        
        
        distance <- end(novel_junctions) - end(ref_junctions)
        # distance <- ifelse(distance < 0, distance + 1, distance - 1)
        
        if (distance %>% length() > 0) {
          df_nd_r <- data.frame(sample = sample,
                                type = "novel_donor", 
                                distance = distance,
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                #novel_counts_norm = novel_junctions$counts / sample_mapped_read_count,
                                novel_seq = novel_junctions %>% seqnames(),
                                novel_start = novel_junctions %>% start(),
                                novel_end = novel_junctions %>% end(),
                                novel_strand = novel_junctions %>% strand(),
                                novel_width = novel_junctions %>% width(),
                                #novel_ss5score = novel_junctions$ss5score,
                                #novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                #ref_ss5score = ref_junctions$ss5score,
                                #ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          
          
          df_nd_r <- df_nd_r %>% 
            group_by(novel_junID) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            ungroup() 
        } else {
          df_nd_r <- NULL
        }
        
        
        
        ####################################################
        ########  NOVEL-ACCEPTOR & FORWARD STRAND  #########
        ####################################################
        
        
        overlaps <- GenomicRanges::findOverlaps(query = all_acceptor_forward,
                                                subject = all_annotated_forward,
                                                ignore.strand = FALSE,
                                                type = "start")
        
        
        novel_junctions <- all_acceptor_forward[queryHits(overlaps),]
        ref_junctions <- all_annotated_forward[subjectHits(overlaps),]
        
        
        distance = end(novel_junctions) - end(ref_junctions)
        # distance = ifelse(distance < 0, distance + 1, distance - 1)
        
        if (distance %>% length() > 0) {
          df_na_f <- data.frame(sample = sample,
                                type = "novel_acceptor", 
                                distance = distance,
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                #novel_counts_norm = novel_junctions$counts / sample_mapped_read_count,
                                novel_seq = novel_junctions %>% seqnames(),
                                novel_start = novel_junctions %>% start(),
                                novel_end = novel_junctions %>% end(),
                                novel_strand = novel_junctions %>% strand(),
                                novel_width = novel_junctions %>% width(),
                                #novel_ss5score = novel_junctions$ss5score,
                                #novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                #ref_ss5score = ref_junctions$ss5score,
                                #ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          
          df_na_f <- df_na_f %>% 
            group_by(novel_junID) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            ungroup() 
          
        } else {
          df_na_f <- NULL
        }
        
        
        
        
        ####################################################
        ########  NOVEL-ACCEPTOR & REVERSE STRAND  #########
        ####################################################
        
        overlaps <- GenomicRanges::findOverlaps(query = all_acceptor_reverse,
                                                subject = all_annotated_reverse,
                                                ignore.strand = FALSE,
                                                type = "end")
        
        
        novel_junctions <- all_acceptor_reverse[queryHits(overlaps),]
        ref_junctions <- all_annotated_reverse[subjectHits(overlaps),]
        
        
        distance = start(ref_junctions) - start(novel_junctions)
        # distance = ifelse(distance < 0, distance + 1, distance - 1)
        
        if (distance %>% length() > 0) {
          df_na_r <- data.frame(sample = sample,
                                type = "novel_acceptor", 
                                distance = distance,
                                novel_junID = novel_junctions$junID,
                                novel_counts = novel_junctions$counts,
                                #novel_counts_norm = novel_junctions$counts / sample_mapped_read_count,
                                novel_seq = novel_junctions %>% seqnames(),
                                novel_start = novel_junctions %>% start(),
                                novel_end = novel_junctions %>% end(),
                                novel_strand = novel_junctions %>% strand(),
                                novel_width = novel_junctions %>% width(),
                                #novel_ss5score = novel_junctions$ss5score,
                                #novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                #ref_ss5score = ref_junctions$ss5score,
                                #ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          
          df_na_r <- df_na_r %>% 
            group_by(novel_junID) %>%
            dplyr::filter(ref_counts == max(ref_counts)) %>%
            dplyr::filter(abs(distance) == min(abs(distance))) %>%
            ungroup() 
        } else {
          df_na_r <- NULL
        }
        
        #################################
        ########  SAVE RESULTS  #########
        #################################
        
        distances <- rbind(df_nd_f, 
                           df_nd_r,
                           df_na_f,
                           df_na_r)
        
        
  
        distances <- distances %>% 
          dplyr::select(-novel_width, -ref_width, -sample)
        
        
        saveRDS(object = distances,
                file = paste0(folder_name, "/", cluster, "_", sample, "_distances.rds"))
        
        
        
        
        num_sample <- num_sample + 1
        
        ## FREE-UP SOME MEMORY
        rm(split_reads_details_97_sample)
        rm(split_read_counts_sample)
        
        
        rm(all_annotated)
        rm(all_annotated_forward)
        rm(all_annotated_reverse)
        rm(all_donor)
        rm(all_donor_forward)
        rm(all_donor_reverse)
        rm(all_acceptor)
        rm(all_acceptor_forward)
        rm(all_acceptor_reverse)
        rm(overlaps)
        rm(distances)
        rm(df_nd_f)
        rm(df_nd_r)
        rm(df_na_f)
        rm(df_na_r)
        rm(novel_junctions)
        rm(ref_junctions)
        gc()
      
      } else {
        print(paste0("Error: sample '", sample, "' not found in the split-reads file!"))
        break;
      }
      
    } else {
      print(paste0("Sample '", sample, "' exists!"))
    }
  }
}


extract_distances <- function(cluster,
                              samples,
                              #recount.info,
                              split_read_counts,
                              folder_name) {

  if ( !file.exists(paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) ) {
    
    ## Obtain the distances across all samples
    df_all <- map_df(samples, function(sample) { 
      
      # sample <- samples[1]
      # sample <- "GTEX-ZVT2-0426-SM-5E44S.1"
      
      file_name <- paste0(folder_name, "/", cluster, "_", sample, "_distances.rds")
      
      if (file.exists(file_name)) {
        
        print(paste0("sample: ", sample))
        
        df <- readRDS(file = file_name)
        
        return(df)
      } else {
        print(paste0("sample: ", sample, " doesn't exist!"))
        break;
      }
      
    })
    saveRDS(object = df_all %>%
              distinct(novel_junID, ref_junID, .keep_all = T) %>%
              mutate(tissue = cluster),
            file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))
    print("Samples loaded!")
    
    
    
    ## RELEASE SOME MEMORY
    rm(df_all)
    rm(df_tidy)
    rm(df_merged)
    rm(df_notambig)
    rm(split_read_counts_novel)
    rm(split_read_counts_ref)
    rm(split_read_counts)
    gc()
  } else {
    print(paste0("File '", cluster, "_raw_distances_tidy.rds' already created!"))
  }
}





## NEVER MIS-SPLICED --------------------------------


# gtf_version <- 105
# project_id <- "BRAIN"
# cluster <- "Brain - Frontal Cortex (BA9)"
# folder_root <- "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN/"
# samples <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster,  "/", project_id, "_", cluster,  "_samples.rds"))
# folder_name <- paste0(folder_root, "/results/pipeline3/distances/", cluster, "/v105/")
# split_read_counts <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/", project_id, "_", cluster, "_split_read_counts_", gtf_version, ".rds"))
# all_split_reads_details <- readRDS(file = paste0(folder_root, "/results/base_data/", cluster, "/",cluster, "_annotated_SR_details_length_", gtf_version, ".rds"))
get_never_misspliced <- function(cluster, 
                                 samples,
                                 split_read_counts,
                                 all_split_reads_details = NULL,
                                 all_not_misspliced = NULL,
                                 folder_name,
                                 save_results = T) {
  
  
  if (is.null(all_not_misspliced)) {
    ## Get Never Mis-spliced!
    
    print(paste0(Sys.time(), " - filtering junctions that are potentially not mis-spliced..."))
    
    if (file.exists(paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds"))) {
      
      # ## This should be zero
      # if (setdiff(all_split_reads_details$junID,split_read_counts$junID) %>% length() > 0) {
      #   print("ERROR")
      #   break;
      # }
      # The filter above should not be zero as the 'split_read_counts' file contains
      # split reads shorter than 25 bp etc
      
      
      ## Remove all '*' from split_read_counts
      all_split_reads_details$junID %>% length()
      
      split_read_counts <- split_read_counts %>%
        inner_join(y = all_split_reads_details%>% dplyr::select(seqnames, start, end, width, strand, junID),
                   by=c("junID" = "junID")) 
      
      ############################
      ## QC 
      ############################

      ## Split_Read_Counts file
      ind <- which(str_detect(string = split_read_counts$junID, pattern = "\\*"))
      if (ind %>% length() > 0) {
        split_read_counts[ind, "junID"] <- str_replace(string = split_read_counts[ind, "junID"]$junID, 
                                                       pattern = "\\*", 
                                                       replacement = split_read_counts[ind, "strand"]$strand %>% as.character() )
      }
      split_read_counts <- split_read_counts %>%
        dplyr::select(-seqnames, -start, -end, -width, -strand)
      any(str_detect(string = split_read_counts$junID, pattern = "\\*"))
      
      ## Annotated_split_reads file
      ind <- which(str_detect(string = all_split_reads_details$junID, pattern = "\\*"))
      if (ind %>% length() > 0) {
        all_split_reads_details[ind, "junID"] <- str_replace(string = all_split_reads_details[ind, "junID"]$junID, 
                                                               pattern = "\\*", 
                                                               replacement = all_split_reads_details[ind, "strand"]$strand %>% as.character() )
      }
      
      ## QC - 2: this should be zero
      ind <- which(str_detect(string = split_read_counts$junID, pattern = "\\*"))
      if (ind %>% length() > 0) {
        print("ERROR")
        break;
      }
      
      
      ## Get mis-spliced and junctions that were not mis-spliced in the first round
      df_all_misspliced <- readRDS(file = paste0(folder_name, "/", cluster, "_raw_distances_tidy.rds")) %>%
        data.table::as.data.table() %>%
        dplyr::distinct(ref_junID, .keep_all = T)
      
      all_not_misspliced <- all_split_reads_details %>%
        data.table::as.data.table() %>%
        dplyr::filter(!(junID %in% df_all_misspliced$ref_junID))
      
      if (intersect(all_not_misspliced$junID, df_all_misspliced$ref_junID) %>% length() > 0){
        print("Error: some of the never misspliced introns appear within the 'distances' file!")
      }
      if (any(df_all_misspliced$ref_junID %>% duplicated()) || any(all_not_misspliced$junID %>% duplicated())) {
        print("Error: none of the mis-spliced introns should within the never mis-spliced introns!")
      }
    } else {
      print(paste0("File: '", folder_name, "/", cluster, "_distances_raw.rds' doesn't exist."))
    }
    
  }
  
  
  num_sample <- 1
  junc_ignore <- NULL
  junc_not_misspliced <- NULL
  
  
  
  if ( !file.exists(paste0(folder_name, "/not-misspliced/", cluster, "_all_notmisspliced.rds")) || 
       !file.exists(paste0(folder_name, "/not-misspliced/", cluster, "_all_misspliced_not_paired.rds")) ) {
 
    ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
    for (sample in samples) { 
      
      # sample <- samples[1]
      # sample <- samples[16]
      
      if ( length(which((colnames(split_read_counts) == sample) == TRUE)) == 1 ) { # Every sampleID is unique, so the result should be equal to 1
        
        #distances <- list()
        split_read_counts_sample <- split_read_counts %>%
          as_tibble()%>%
          dplyr::select(junID, all_of(sample %>% as.character())) 
        
        split_read_counts_sample[split_read_counts_sample==0] <- NA
        
        split_read_counts_sample <- split_read_counts_sample %>%
          drop_na()
        
        split_reads_details_97_sample <- merge(all_not_misspliced, 
                                               split_read_counts_sample, 
                                               by = "junID", 
                                               sort = TRUE) %>%
          dplyr::rename(counts = all_of(sample %>% as.character()))
        
        split_reads_details_97_sample <- split_reads_details_97_sample %>% as_tibble()
        all_not_misspliced <- all_not_misspliced %>% as_tibble()
        
        
        ## Do some QC
        index <- runif(n = 1, 1, split_reads_details_97_sample %>% nrow()) %>% as.integer()
        data <- split_reads_details_97_sample[index, c(1, split_reads_details_97_sample %>% ncol())]
        data_col <- split_read_counts_sample %>%
          dplyr::filter(junID == data$junID)
        if (data_col[[2]] != data$counts) {
          print("Error: QC failed!")
          break;
        }
        
        ## Print an informative message
        print(cluster)
        
        ###################################################
        ##########    ANNOTATED JUNC   ####################
        ###################################################
        all_annotated <- split_reads_details_97_sample %>%
          dplyr::filter(type == "annotated") %>% 
          GenomicRanges::GRanges()
        
        print(paste0("sample '", num_sample, "': ", 
                     length(all_annotated$junID), " initial ref junctions."))
        
        all_annotated_forward <- all_annotated[all_annotated@strand == "+"]
        all_annotated_reverse <- all_annotated[all_annotated@strand == "-"]
        
        
        ###################################################
        ##########    NOVEL DONOR      ####################
        ###################################################
        all_donor <- split_reads_details_97_sample %>%
          dplyr::filter(type == "novel_donor") %>% 
          GenomicRanges::GRanges()
        
        all_donor_forward <- all_donor[all_donor@strand == "+"]
        all_donor_reverse <- all_donor[all_donor@strand == "-"]
        
        ###################################################
        ##########    NOVEL ACCEPTOR      #################
        ###################################################
        all_acceptor <- split_reads_details_97_sample%>%
          dplyr::filter(type == "novel_acceptor") %>% 
          GRanges()
        
        all_acceptor_forward <- all_acceptor[all_acceptor@strand == "+"]
        all_acceptor_reverse <- all_acceptor[all_acceptor@strand == "-"]
        
        if (all_donor_forward %>% length() > 0 &&
            all_donor_reverse %>% length() > 0 &&
            all_annotated_forward %>% length() > 0 &&
            all_annotated_reverse %>% length() > 0) {
          
          ####################################################
          ########  NOVEL-DONOR & FORWARD STRAND  ############
          ####################################################
          
          
          overl_df <- GenomicRanges::findOverlaps(query = all_donor_forward,
                                                  subject = all_annotated_forward,
                                                  ignore.strand = FALSE,
                                                  type = "end")
          ignore_df <- all_annotated_forward[subjectHits(overl_df),]$junID
          ref_annotated_d_forward <- all_annotated_forward[-subjectHits(overl_df),]$junID
          
          
          
          ####################################################
          ########  NOVEL-DONOR & REVERSE STRAND  ############
          ####################################################
          
          overl_dr <- GenomicRanges::findOverlaps(query = all_donor_reverse,
                                                  subject = all_annotated_reverse,
                                                  type = "start",
                                                  ignore.strand = FALSE)
          ignore_dr <- all_annotated_reverse[subjectHits(overl_dr),]$junID
          ref_annotated_d_reverse <- all_annotated_reverse[-subjectHits(overl_dr),]$junID
        
        } else {
          ignore_df <- NULL
          ref_annotated_d_forward <- NULL
          ignore_dr <- NULL
          ref_annotated_d_reverse <- NULL
        }
        
        ####################################################
        ########  NOVEL-ACCEPTOR & FORWARD STRAND  #########
        ####################################################
        
        if (all_acceptor_forward %>% length() > 0 &&
            all_acceptor_reverse %>% length() > 0 &&
            all_annotated_forward %>% length() > 0 &&
            all_annotated_reverse %>% length() > 0) {
        
          overl_af <- GenomicRanges::findOverlaps(query = all_acceptor_forward,
                                                  subject = all_annotated_forward,
                                                  type = "start",
                                                  ignore.strand = FALSE)
          ignore_af <- all_annotated_forward[subjectHits(overl_af),]$junID
          ref_annotated_a_forward <- all_annotated_forward[-subjectHits(overl_af),]$junID
          
          
          ####################################################
          ########  NOVEL-ACCEPTOR & REVERSE STRAND  #########
          ####################################################
          
          overl_ar <- GenomicRanges::findOverlaps(query = all_acceptor_reverse,
                                                  subject = all_annotated_reverse,
                                                  type = "end",
                                                  ignore.strand = FALSE)
          ignore_ar <- all_annotated_reverse[subjectHits(overl_ar),]$junID
          ref_annotated_a_reverse <- all_annotated_reverse[-subjectHits(overl_ar),]$junID
          
        } else {
          ignore_af <- NULL
          ref_annotated_a_forward <- NULL
          ignore_ar <- NULL
          ref_annotated_a_reverse <- NULL
        }
        
        junc_ignore <- c(junc_ignore,
                         ignore_af,
                         ignore_ar,
                         ignore_df,
                         ignore_dr) %>% unique()
        
        
        junc_not_misspliced <- c(junc_not_misspliced,
                                 ref_annotated_d_forward, 
                                 ref_annotated_d_reverse,
                                 ref_annotated_a_forward, 
                                 ref_annotated_a_reverse) %>% unique()
        
        
        
        
        print(paste0(Sys.time(), " --> Sample '", num_sample, "' processed // ", c(ref_annotated_d_forward, 
                                                                                                 ref_annotated_d_reverse,
                                                                                                 ref_annotated_a_forward, 
                                                                                                 ref_annotated_a_reverse) %>% unique() %>% length(), 
                     " junctions not mis-spliced // ", c(ignore_af,
                                                        ignore_ar,
                                                        ignore_df,
                                                        ignore_dr) %>% unique %>% length(), " mis-spliced junctions."))
        
        
        if (!save_results) {
          if (junc_ignore %>% length() != 0) {
            print("Error: some never mis-spliced introns have been paired to a novel junction.")
            break;
          }
        }
        
        
        num_sample <- num_sample + 1
        
        rm(ignore_af)
        rm(ignore_ar)
        rm(ignore_df)
        rm(ignore_dr)
        rm(ref_annotated_d_forward)
        rm(ref_annotated_d_reverse)
        rm(ref_annotated_a_forward) 
        rm(ref_annotated_a_reverse)
        rm(overl_df)
        rm(overl_dr)
        rm(overl_af)
        rm(overl_ar)
        rm(split_read_counts_sample)
        rm(split_reads_details_97_sample)
        rm(all_annotated)
        rm(all_annotated_forward)
        rm(all_annotated_reverse)
        rm(all_donor)
        rm(all_donor_forward)
        rm(all_donor_reverse)
        rm(all_acceptor)
        rm(all_acceptor_forward)
        rm(all_acceptor_reverse)
        rm(data)
        rm(data_col)
        gc()
        
      }
    }
    
    if ((which(junc_not_misspliced %in% junc_ignore) %>% length()) > 0) {
      junc_not_misspliced <- junc_not_misspliced[-which(junc_not_misspliced %in% junc_ignore)]  
    }
    
    if (save_results) {
      
      if (intersect(junc_not_misspliced,junc_ignore) %>% length() > 0) {
        print(paste0(Sys.time(), " - Error! There are overlapping not-misspliced and misspliced junctions!"))
        
      } else { 
        
        print(paste0(Sys.time(), " - ", junc_not_misspliced %>% length(), " not mis-spliced junctions!"))
        folder_name <- paste0(folder_name, "/not-misspliced/")
        dir.create(file.path(folder_name), showWarnings = F)
        
        saveRDS(object = junc_not_misspliced, 
                file = paste0(folder_name, "/", cluster, "_all_notmisspliced.rds"))
        
        saveRDS(object = junc_ignore, 
                file = paste0(folder_name, "/", cluster, "_all_misspliced_not_paired.rds"))
        
        print(paste0(Sys.time(), " - results saved!"))
      }
      
    } else {
      return(junc_not_misspliced)
    }
    
    
    ## RELEASE SOME MEMORY
    rm(num_sample)
    rm(junc_ignore)
    rm(junc_not_misspliced)
    #rm(df_junc_counts)
    rm(split_reads_details_97_sample)
    rm(all_annotated)
    rm(all_annotated_forward)
    rm(all_annotated_reverse)
    rm(all_donor)
    rm(all_donor_forward)
    rm(all_donor_reverse)
    rm(all_acceptor)
    rm(all_acceptor_forward)
    rm(all_acceptor_reverse)
    rm(overl_df)
    rm(ignore_df)
    rm(ref_annotated_d_forward)
    rm(overl_dr)
    rm(ignore_dr)
    rm(ref_annotated_d_reverse)
    rm(overl_af)
    rm(ignore_af)
    rm(ref_annotated_a_forward)
    rm(overl_ar)
    rm(ignore_ar)
    rm(ref_annotated_a_reverse)
    gc()
  }
  
}






get_missplicing_QC <- function(cluster,
                               samples,
                               split_read_counts,
                               all_split_reads_details,
                               folder_name) {
  
  # cluster = "PD"
  # cluster = "control"
  #folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/", project_id, "/results/pipeline3/missplicing-ratio/", cluster, "/")
  
  ################################################################################################################
  ## Load the intron database
  
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>% as_tibble()
  df_introns %>% head()
  df_introns %>% nrow()
  
  
  
  
  
  ################################################################################################################
  ## Double-check again that the never mis-spliced junctions can't actually be paired with any novel donor or 
  ## novel acceptor from the current tissue
  
  
  ## Get the annotation details from the never mis-spliced introns
  all_not_misspliced <- all_split_reads_details %>%
    dplyr::filter(junID %in% (df_introns %>%
                                dplyr::filter(ref_type == "never") %>%
                                distinct(ref_junID) %>%
                                pull()))
  
  
  ## Add the rest of the novel donor and acceptor junctions
  all_not_misspliced <- rbind(all_not_misspliced,
                              all_split_reads_details %>%
                                dplyr::filter(type %in% c("novel_donor", "novel_acceptor")))
  
  
  
  ## Call the function. The result returned should be of length zero
  never_misspliced_pairings <- get_never_misspliced(cluster = cluster,
                                                    samples = samples,
                                                    split_read_counts = split_read_counts,
                                                    all_not_misspliced = all_not_misspliced,
                                                    folder_name = folder_name,
                                                    save_results = F)
  
  
  if (never_misspliced_pairings %>% length() > 0) {
    print("Error: some introns catalogued as 'never mis-spliced' have been found to be mis-spliced.")
    break;
  }
  
  
  
  ################################################################################################################
  ## Load novel junctions
  
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds")) %>% as_tibble()
  df_novel %>% head()
  df_novel %>% nrow()
  
  
  ## Novel junctions from the novel database are only novel donor or novel acceptor
  if (any(!(df_novel$novel_type %>% unique()) %in% c("novel_acceptor", "novel_donor"))) {
    print("Error: some novel junctions from the novel database are not classified as 'novel donor' or 'novel acceptor'.")
    break;
  }
  
  
  ################################################################################################################
  ## All novel junctions should be attached to exact the same number of introns stored within the intron database
  
  if (df_introns %>%
      dplyr::filter(ref_type != "never") %>% 
      distinct(ref_junID) %>%
      nrow() != df_novel %>%
      distinct(ref_junID) %>%
      nrow()) {
    print("Error: some novel junctions aren't attached to the same introns stored within the intron database.")
    break;
  }
  
  
  
  if (intersect(df_novel %>%
                distinct(ref_junID) %>%
                pull(), 
                df_introns %>%
                dplyr::filter(ref_type != "never") %>%
                distinct(ref_junID) %>%
                pull()) %>% length() != df_introns %>% 
      dplyr::filter(ref_type != "never") %>% nrow()){
    print(paste0("Error: some ", cluster, " novel junctions are attached to introns that are not stored on the database!"))
  }
  
  if (intersect(df_novel %>%
                distinct(novel_junID) %>%
                pull() %>% sort(), all_split_reads_details %>%
                dplyr::filter(type %in% c("novel_donor", "novel_acceptor")) %>%
                distinct(junID) %>%
                pull() %>% sort()) %>% length() != df_novel %>% nrow()){
    
    
    diff <-setdiff(df_novel %>%
                     distinct(novel_junID) %>%
                     pull(), all_split_reads_details %>%
                     dplyr::filter(type %in% c("novel_donor", "novel_acceptor")) %>%
                     distinct(junID) %>%
                     pull())
    
    all_split_reads_details %>%
      dplyr::filter(junID == diff[1])
    
    
    all_split_reads_details %>%
      dplyr::filter(junID == "JUNC00006215")
    
    print(paste0("Error: some ", cluster, " novel junctions have not been found within the original annotation!"))
  }
  
  
  
  
  ## FREE SOME MEMORY
  
  rm(never_misspliced_pairings)
  rm(all_not_misspliced)
  rm(df_introns)
  rm(df_novel)
  gc()
  
}




add_cdts_cons_scores <- function(cluster = NULL,
                                 db_introns = NULL,
                                 folder_name = NULL) {
  
  
  
  
  
  
  if (is.null(db_introns) && !is.null(cluster)) {
    ## Load the IDB 
    db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
      distinct(ref_junID, .keep_all = T)
  }
  
  db_introns <- db_introns %>%
    mutate(CDTS_5ss_mean = 0.0,
           CDTS_3ss_mean = 0.0,
           phastCons20way_5ss_mean = 0.0,
           phastCons20way_3ss_mean = 0.0) %>%
    GRanges()
  
  
  print(paste0(Sys.time(), " - getting 5' scores assigned to introns from IntroVerse ..."))
  
  ## https://www.nature.com/articles/nature09000
  ## Scores from the 5'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% start() - 5,
                                                                                            end = db_introns %>% start() + 35),
                                                                           strand = db_introns %>% strand()),
                                          ignore.strand = FALSE,
                                          type = "any")
  
  overlaps_tidy <- overlaps %>%
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_5ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_5ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  
  
  print(paste0(Sys.time(), " - getting 3' scores assigned to introns from IntroVerse ..."))
  
  ## Scores from the 3'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% end() - 35,
                                                                                            end = db_introns %>% end() + 5),
                                                                           strand = db_introns %>% strand()),
                                          
                                          ignore.strand = FALSE,
                                          type = "any")
  
  overlaps_tidy <- overlaps %>% 
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    as_tibble() %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_3ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_3ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  #####################
  ## SAVE RESULTS
  #####################
  
  if (!is.null(cluster)) {
    file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
    saveRDS(object = db_introns %>% data.table::as.data.table(),
            file = file_name)
    
    print(paste0(Sys.time(), " - CDTS and Conservation scores added! IDB updated!"))
    
    
    rm(overlaps_tidy)
    rm(db_introns)
    rm(file_name)
    rm(overlaps)
    
  } else {
    return(db_introns)
  }
  
  
  
  
  #gc()
}

####################################################
## ADDING FEATURES TO THE INTRON DATABASE ##########
####################################################



## MITOCHONDRIAL GENES --------------------------

remove_MT_genes <- function(cluster,
                            folder_name) {
  
  
  
  ## Load mis-splicing ratios df
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  df_introns %>% head()
  df_introns %>% nrow()
  
  
  ## Load MT genes
  
  print(paste0(Sys.time(), " - checking the existance of MT genes..."))
  
  MT_genes <- readRDS(file = "/data/references/MT_genes/MT_genes.rds")
  MT_geneID <- MT_genes %>% distinct(gene_id) %>% pull(gene_id)
  
  if (any(df_introns$gene_id %in% MT_geneID)) {
    
    ## Remove MT genes
    df_introns <- df_introns %>%
      dplyr::filter(!(gene_id %in% MT_geneID)) 
    
    print(paste0(Sys.time(), " - MT genes removed! ", df_introns %>% nrow(), " final number of junctions."))
    
    
    ## SAVE RESULTS
    saveRDS(object = df_introns,
            file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
    
    print(paste0(Sys.time(), " - results saved!"))
    
  } else {
    print(paste0(Sys.time(), " - the dataset doesn't contain any MT gene!"))
  }
  
  
  ## FREE SOME MEMORY
  rm(df_introns)
  rm(MT_geneID)
  rm(MT_genes)
  gc()
  
}





## GENE TPM AND GENE LENGTH -----------------------

# folder_name <- "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/ADIPOSE_TISSUE//results/pipeline3/missplicing-ratio/Adipose - Subcutaneous/v105/"
# cluster <- "Adipose - Subcutaneous"

add_gene_tpm_length <- function(cluster,
                                tpm_file,
                                folder_name,
                                GTEx = T) {
  
  print(cluster)
  
  ## LOAD SOURCE DATA
  file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
  df_introns <- readRDS(file = file_name) %>%
    as_tibble()
  
  ## PRINT SOME STATS
  df_introns %>% head()
  df_introns %>% nrow()
  df_introns <- df_introns %>%
    distinct(ref_junID, .keep_all = T) # %>%
  df_introns %>% nrow()
  
  
  #########################
  ## GENE TPM 
  
  if (GTEx) {
    
    print(paste0(Sys.time(), " - Adding gene TPM..."))
    
    ## Add gene TPM 
    
    # library(GSRI)
    
    tpm <- readRDS(file = tpm_file)

    tpm %>% head()
    tpm %>% nrow()
    
    tpm <- tpm  %>% 
      dplyr::select(gene_id = gene, 
                    tpm_median = TPM_median,
                    tpm_mean = TPM_mean)
    
    tpm %>% head()
    tpm %>% nrow()
    
    
  } else {
    
    tpm <- readRDS(file = tpm_folder) %>%
      dplyr::select(-tpm)
  }
  
  
  
  
  
  df_introns <- df_introns %>%
    unnest(gene_id)
  
  
  df_merged <- merge(x = df_introns %>% data.table::as.data.table(),
                     y = tpm %>% data.table::as.data.table(),
                     by = "gene_id",
                     all.x = T) %>% 
    dplyr::rename(tpm_mean_rct = tpm_mean) %>% 
    dplyr::rename(tpm_median_rct = tpm_median)
  
  
  df_merged %>% head()
  df_merged %>% nrow()
  

  
  ############################
  ## ADD NUMBER OF TRANSCRIPTS
  ############################
  
  
  print(paste0(Sys.time(), " - Adding number of transcripts..."))
  
  
  if (!exists("homo_sapiens_v105")) {
    
    homo_sapiens_v105 <- rtracklayer::import("/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf") %>%
      as.data.frame()
    
    hsv105_transcripts <- homo_sapiens_v105 %>%
      dplyr::count(gene_id, type) %>%
      dplyr::filter(type == "transcript") %>%
      unnest(gene_id)
    
  }
  
  df_merged <- merge(x = df_merged %>% data.table::as.data.table(),
                     y = hsv105_transcripts[,c("gene_id", "n")] %>% data.table::as.data.table(),
                     by = "gene_id",
                     all.x = T) %>%
    dplyr::rename(n_transcripts = n)
  
  df_merged %>% head()
  df_merged %>% nrow()
  
  
  if (any(is.na(df_merged$n_transcripts))) {
    print("Error: some genes don't have their number of transcript attached!")
  }
  
  
  
  #########################
  ## ADD GENE LENGTH 
  #########################
  
  print(paste0(Sys.time(), " - Adding gene length..."))
  
  ## The matching is going to be donr
  hsv105_genes <- homo_sapiens_v105 %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_width = width)
  
  
  df_merged %>% head()
  df_merged %>% nrow()
  
  df_merged <- merge(x = df_merged %>% data.table::as.data.table(),
                     y = hsv105_genes %>% data.table::as.data.table(),
                     by = "gene_id",
                     all.x = T)
  
  
  df_merged %>% head()
  df_merged %>% nrow()
  
  
  if (any(is.na(df_merged$gene_width))) {
    print("Error: some genes don't have their length attached!")
  }
  
  if (any(df_merged$width > df_merged$gene_width)) {
    print(paste0(Sys.time(), " - removing introns longer than their assigned gene!"))
    
    df_merged %>% nrow()
    
    df_merged <- df_merged %>%
      dplyr::filter(width < gene_width)
    
    df_merged %>% nrow()
  }
  
  df_merged %>% head()
  df_merged %>% nrow()
  
  ## Remove duplicated ref introns (this is due to the 'unnest' function applied before) ----------------
  
  df_merged <- df_merged %>%
    distinct(ref_junID, .keep_all = T)
  
  
  ## Matching results with novel introns ---------------------------------------
  
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))%>%
    dplyr::filter(ref_junID %in% df_merged$ref_junID)
  
  if (!identical(df_novel$ref_junID %>% unique() %>% sort(), 
                 df_merged %>% dplyr::filter(ref_type != "never") %>% dplyr::select(ref_junID) %>% pull()%>% unique() %>% sort())) {
    print("ERROR!")
  }
  
  ## Save results --------------------------------------------------------------
  
  file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
  saveRDS(object = df_merged, file = file_name)
  
  file_name <- paste0(folder_name, "/", cluster, "_db_novel.rds")
  saveRDS(object = df_novel, file = file_name)
  
  print(paste0(Sys.time(), " - file saved!"))
  
  
  ## Free some memory ----------------------------------------------------------
  rm(hsv105_transcripts)
  rm(homo_sapiens_v105)
  rm(hsv105_transcripts)
  rm(df_merged)
  rm(df_introns)
  rm(tpm)
  gc()
  
}


QC_IDB_junctions <- function() {
  
  for (cluster in gtex_tissues) {
    
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", cluster, "/v104/")
    file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
    
    ## IDB - Introns
    df_introns <- readRDS(file = file_name) %>% as.data.frame()
    
    ## IDB - Novel junctions
    df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds")) %>% as_tibble()
    
    ## Compare
    if (!identical(df_novel$ref_junID %>% unique() %>% sort(), 
                   df_introns %>% dplyr::filter(ref_type != "never") %>% dplyr::select(ref_junID) %>% pull()%>% unique() %>% sort())) {
      print("ERROR!")
      break;
    } else {
      print(paste0(cluster, " correct!"))
    }
    
    if (any(db_introns$u2_intron == F)){
      print("ERROR!")
      break;
    }
    
    if (any(db_introns$width > db_introns$gene_width)) {
      print("ERROR!")
      break;
    }
    
  }
  
  # cluster <- "Brain-FrontalCortex_BA9"
  # version <- "104"
  # 
  
  ## IDB - Introns
  
  print(cluster)
  
}


## PROTEIN PERCENTAGE  ----------------------------



#' Title
#'
#' @param df_introns 
#' @param cluster 
#' @param folder_root 
#'
#' @return
#' @export
#'
#' @examples
add_protein_percentage_to_database <- function(db_introns = NULL,
                                               cluster = NULL,
                                               folder_root = NULL) {
  
  
  ## Import HUMAN REFERENCE transcriptome
  # homo_sapiens_v104_gtf <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf") %>% 
  #   as.data.frame()
  
  print(paste0(Sys.time(), " - adding biotypes percentage to junctions..."))
  
  df_protein <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_annotated_SR_details_length_104_biotype.rds") %>%
    as_tibble()
  df_protein %>% head()
  df_protein %>% nrow()
  
  
  
  ## Load mis-splicing data for the current cluster
  
  
  if (!is.null(cluster)) {
    
    print(paste0(Sys.time(), " - Adding biotype percentage to '", cluster, "' junctions..."))
    
    if (cluster == "case" || cluster == "control") {
      folder_name <- paste0(folder_root, "/", cluster, "/pipeline3/missplicing-ratio/")
      db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
    } else {
      folder_name <- folder_root #paste0(folder_root, "/results/pipeline3/missplicing-ratio/", cluster, "/v105")
      db_introns <- readRDS(file = paste0(folder_root, "/", cluster, "_db_introns.rds")) %>%
        dplyr::filter(u12_intron == T | u2_intron == T)
    }
    
    ## Load database of novel junctions
    df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))  %>%
      dplyr::filter(ref_junID %in% db_introns$ref_junID)
    
    df_novel %>% head()
    df_novel %>% distinct(ref_junID) %>% nrow()
    
    ## QC - novel junctions must be attached to the same introns stored on the intron database
    if ((db_introns %>% dplyr::filter(ref_type != "never") %>% distinct(ref_junID) %>% nrow()) != 
        (df_novel %>% distinct(ref_junID) %>% nrow())) {
      print("Error: some novel junctions are attached to introns that are not stored on the intron database!")
    }
    
  }
  
  
  # db_introns %>% head()
  # db_introns %>% nrow()
  # db_introns %>% distinct(ref_junID) %>% nrow()
  
  df_protein_local <- df_protein %>%
    dplyr::filter(junID %in% db_introns$ref_junID) 
  
  # df_protein_local %>% head()
  # df_protein_local %>% nrow()
  # db_introns %>% distinct(ref_junID) %>% nrow()
  
  
  
  if ((db_introns %>% distinct(ref_junID) %>% nrow()) == 
      (df_protein_local %>% distinct(junID) %>% nrow())) {
    
    #db_introns %>% nrow()
    
    ## Add protein-coding to the IDB
    
    db_introns <- merge(x = db_introns %>% data.table::as.data.table(),
                        y = df_protein_local %>% data.table::as.data.table(),
                        by.x = "ref_junID",
                        by.y = "junID",
                        all.x = T)
    #db_introns %>% head()
    
    
    
    
    ## Add protein-coding info to the novel junction database
    if (!is.null(cluster)) {
      df_novel <- merge(x = df_novel %>% data.table::as.data.table(),
                        y = df_protein_local %>% data.table::as.data.table(),
                        by.x = "ref_junID",
                        by.y = "junID",
                        all.x = T)
      df_novel %>% head()
      
      if(identical(db_introns %>% dplyr::filter(ref_type != "never") %>% distinct(ref_junID) %>% pull %>% sort(), 
                   df_novel %>% distinct(ref_junID) %>% pull %>% sort())) {
        
        saveRDS(object = db_introns %>% data.table::as.data.table(),
                file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
        
        saveRDS(object = df_novel,
                file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
        
        print(paste0(Sys.time(), " - ", cluster, " finished!"))
        
      } else {
        print(paste0(Sys.time(), " - ", cluster, " error - the two datasets present different refIDs of rows!"))
      }
      
    } else {
      return(db_introns)
    }
    
    
    
    
    
    
    
  } else {
    return(NULL)
  }
  
  rm(df_protein_local)
  rm(db_introns)
  rm(df_novel)
  
  gc()
  
}

add_protein_percentage_to_database_QC <- function(tissues = gtex_tissues,
                                                  folder_root = "/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/",
                                                  protein_folder = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_annotated_SR_details_length_104_biotype.rds") {
  
  ## Load database of introns
  cluster <- tissues[11]
  folder_name <- paste0(folder_root, "/", cluster, "/v104")
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
    dplyr::filter(u12_intron == T | u2_intron == T)
  
  df_protein <- readRDS(file = protein_folder)
  df_protein %>% head()
  df_protein %>% nrow()
  
  homo_sapiens_v104_gtf <- rtracklayer::import(con = "/data/references/ensembl/gtf_gff3/v104/Homo_sapiens.GRCh38.104.gtf") %>% 
    as.data.frame()
  homo_sapiens_v104_gtf[1,]
  
  
  row_number <- sample(1:nrow(df_introns), 1, replace = F)   
  
  tx <- df_introns[row_number,] %>%
    unnest(tx_id_junction) %>%
    pull(tx_id_junction)
  
  df_protein %>%
    dplyr::filter(junID == df_introns[row_number,]$ref_junID)
  
  homo_sapiens_v104_gtf %>%
    dplyr::filter(type == "transcript",
           transcript_id %in% tx)
  
  
} 



