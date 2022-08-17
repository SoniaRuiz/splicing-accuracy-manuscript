
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

    if (!file.exists(paste0(folder_name, "/", cluster, "_", sample, "_distances.rds"))) {
      
      if (length(which((colnames(split_read_counts) == sample) == TRUE)) == 1) { 
        
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
          dplyr::rename(counts = all_of(sample %>% as.character()))
        
        split_reads_details_97_sample %>% head()
        
        ## Do some QC
        index <- runif(n = 1, 1, split_reads_details_97_sample %>% nrow()) %>% as.integer()
        data <- split_reads_details_97_sample[index, c(1, split_reads_details_97_sample %>% ncol())]
        data_col <- split_read_counts_sample %>%
          filter(junID == data$junID)
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
                                novel_ss5score = novel_junctions$ss5score,
                                novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_end = ref_junctions %>% end(),
                                ref_width = ref_junctions %>% width(),
                                ref_ss5score = ref_junctions$ss5score,
                                ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          ## In case the novel junction has been attached to many ref junctions, the ref junction chosen will be 
          ## the one with the highest number of counts. In case many ref junctions present the same number of counts,
          ## the closest one will be chosen.
          
          df_nd_f <- df_nd_f %>% 
            group_by(novel_junID) %>%
            filter(ref_counts == max(ref_counts)) %>%
            filter(abs(distance) == min(abs(distance))) %>%
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
        #   select(novel_strand, novel_start, novel_end, ref_strand, ref_start, ref_end, distance)
        
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
                                novel_ss5score = novel_junctions$ss5score,
                                novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                ref_ss5score = ref_junctions$ss5score,
                                ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          
          
          df_nd_r <- df_nd_r %>% 
            group_by(novel_junID) %>%
            filter(ref_counts == max(ref_counts)) %>%
            filter(abs(distance) == min(abs(distance))) %>%
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
                                novel_ss5score = novel_junctions$ss5score,
                                novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                ref_ss5score = ref_junctions$ss5score,
                                ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          
          df_na_f <- df_na_f %>% 
            group_by(novel_junID) %>%
            filter(ref_counts == max(ref_counts)) %>%
            filter(abs(distance) == min(abs(distance))) %>%
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
                                novel_ss5score = novel_junctions$ss5score,
                                novel_ss3score = novel_junctions$ss3score,
                                ref_junID = ref_junctions$junID,
                                ref_counts = ref_junctions$counts,
                                #ref_counts_norm = ref_junctions$counts / sample_mapped_read_count,
                                ref_seq = ref_junctions %>% seqnames(),
                                ref_start = ref_junctions %>% start(),
                                ref_end = ref_junctions %>% end(),
                                ref_strand = ref_junctions %>% strand(),
                                ref_width = ref_junctions %>% width(),
                                ref_ss5score = ref_junctions$ss5score,
                                ref_ss3score = ref_junctions$ss3score,
                                stringsAsFactors = FALSE)
          
          
          df_na_r <- df_na_r %>% 
            group_by(novel_junID) %>%
            filter(ref_counts == max(ref_counts)) %>%
            filter(abs(distance) == min(abs(distance))) %>%
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
          select(-novel_width, -ref_width, -sample)
        
        
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
  
  print("Samples loaded!")
  
  # if (df_all %>% nrow() == 0) {
  #   
  #   df_all <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_raw.rds"))
  #   
  # } else {
  #   
  #   saveRDS(object = df_all %>% data.table::as.data.table(),
  #           file = paste0(folder_name, "/", cluster, "_distances_raw.rds"))
  #   print(stringr::str_c(Sys.time(), " - raw 'distances' file saved!")) 
  #   
  # }
  
  # df_all <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_raw.rds"))
  print(paste0(Sys.time(), " --> ", df_all$ref_junID %>% unique() %>% length(), " unique reference introns"))
  print(paste0(Sys.time(), " --> ", df_all$novel_junID %>% unique() %>% length(), " unique novel junctions"))
  
  # SNCA David's novel exon:
  # df %>% filter(novel_junID == "67198600")
  # df %>% filter(novel_junID == "67199298")
  
  # df_all %>%
  #   filter(novel_seq == "4",
  #          novel_start == "89729278",
  #          novel_end == "89771396",
  #          novel_strand == "-")
  # df_all %>%
  #   filter(ref_seq == "4",
  #          ref_start == "89729278",
  #          ref_end == "89771396",
  #          ref_strand == "-")
  #4:89729278-89771396:-
  
  # print(paste0(Sys.time(), " --> ", df_all %>% dplyr::distinct(sample) %>% nrow(), " samples in '", cluster, "' tissue"))
  
  ## STATS OF THE RAW DISTANCES
  any(df_all$distance == 0)
  print(paste0("Novel donor (+) mode: ", get_mode(df_all %>% filter(type == "novel_donor") %>% filter(distance > 0) %>% pull(distance))))
  print(paste0("Novel donor (-) mode: ", get_mode(df_all %>% filter(type == "novel_donor") %>% filter(distance < 0) %>% pull(distance))))
  print(paste0("Novel acceptor (+) mode: ", get_mode(df_all %>% filter(type == "novel_acceptor") %>% filter(distance > 0) %>% pull(distance))))
  print(paste0("Novel acceptor (-) mode: ", get_mode(df_all %>% filter(type == "novel_acceptor") %>% filter(distance < 0) %>% pull(distance))))
  
  ## Introns shouldn't be obtained from sex or MT chromosomes
  if (any((df_all$ref_seq %>% unique()) == "MT") ||
      any((df_all$novel_seq %>% unique()) == "MT")) {
    print(paste0("Error: the dataset contains introns and junctions from sex and MT chromosomes!"))
    break;
  }
  
  ## I check that there're multiple novel junctions duplicated
  if (df_all$novel_junID %>% length() >
      df_all$novel_junID %>% unique() %>% length()) {
    
    print(stringr::str_c(Sys.time(), " - removing ambiguous junctions...: "))
    
  }
  
  ## It should be more novel than reference junctions
  if (!(df_all$novel_junID %>% unique() %>% length() > df_all$ref_junID %>% unique() %>% length())) {
    
    print("Error, there are more introns than novel junctions.")
    break;
    
  }
  
  ## Novel donor junctions cannot also be targeted as novel acceptor
  if (intersect(df_all %>%
                filter(type == "novel_donor") %>%
                rowwise() %>%
                distinct(novel_junID) %>%
                pull(novel_junID),
                
                df_all %>%
                filter(type == "novel_acceptor") %>%
                rowwise() %>%
                distinct(novel_junID) %>%
                pull(novel_junID)) %>% length() > 0) {
    print("Error: there are novel junctions classified as both novel donor and novel acceptor.")
  }
  
  
  
  ## 1. Add two columns to store the mean distances and SD distances per each novel junction across samples
  df_tidy <- df_all %>%
    group_by(novel_junID) %>%
    mutate(mean_distances = distance %>% mean() %>% round()) %>%
    mutate(distances_sd = distance %>% sd()) %>%
    mutate(distances_sd = replace_na(distances_sd, 0)) %>%
    ungroup() %>%
    as.data.frame()
  df_tidy %>% head()
  df_tidy %>% nrow()
  df_tidy$ref_junID %>% unique() %>% length()
  df_tidy$novel_junID %>% unique() %>% length()
  
  
  ## There are multiple novel junctions associated to different introns (we call them AMBIGUOUS JUNCTIONS)
  ## We get their percentage
  print(paste0(Sys.time(), " --> ", ((df_tidy %>%
                                        filter(distances_sd > 0) %>%
                                        distinct(novel_junID) %>%
                                        nrow() * 100) / (df_tidy %>%
                                                           filter(distances_sd == 0) %>%
                                                           distinct(novel_junID) %>%
                                                           nrow())) %>% round(digits = 2), "% of ambiguous junctions"))
  
  
  ## Ambiguous junctions shouldn't overlap with unambiguous junctions
  if (intersect(df_tidy %>%
                filter(distances_sd > 0) %>%
                distinct(novel_junID) %>%
                pull(novel_junID),
                df_tidy %>%
                filter(distances_sd == 0) %>%
                distinct(novel_junID) %>%
                pull(novel_junID)) %>% length() > 0) {
    print("Error: there are novel junctions classified as both ambiguous and unambiguous junctions.")
  }
  
  
  ## 2. Only keep unambiguous novel junctions
  ## We only get the junctions that have the same reference junction across samples (i.e the same distance figures)
  print(stringr::str_c(Sys.time(), " --> removing ambiguous novel junctions...: "))
  
  df_notambig <- df_tidy %>%
    filter(distances_sd == 0) # %>%
  #group_by(novel_junID) %>%
  #summarise(ambiguous = n_distinct(ref_junID)) %>%
  #filter(ambiguous == 1)
  
  # df_notambig <- df_notambig %>%
  #   group_by(novel_junID) %>%
  #   summarise(ambiguous = n_distinct(ref_junID)) %>%
  #   filter(ambiguous == 1)
  
  df_tidy %>% nrow()
  df_tidy <- df_tidy %>%
    filter(novel_junID %in% df_notambig$novel_junID)
  df_tidy %>% nrow()
  df_tidy$ref_junID %>% unique() %>% length()
  df_tidy$novel_junID %>% unique() %>% length()
  
  
  ## Unambiguous junctions must have 0 standard deviation value in their mean distances across samples
  if (any(df_tidy$distances_sd > 0)) {
    print("Error: still there are ambiguous junc.")
    break;
  }
  
  ## After removing ambiguous junctions, the mean distances across samples are equal to the distance at a sample level
  if (!identical(df_tidy$mean_distances %>% sort %>% as.integer(),
                 df_tidy$distance %>% sort %>% as.integer())) {
    print("Error: still there are ambiguous junc.")
    break;
  }
  
  ## Removing unuseful columns
  df_tidy <- df_tidy %>%
    dplyr::select(-c(distances_sd, mean_distances))
  
  
  
  # ## 3. Remove duplicated rows, so the same novel junction only present one distance figure across the samples
  # df_tidy$novel_junID %>% length()
  # df_tidy$novel_junID %>% unique() %>% length()
  # df_tidy %>% head
  #
  # df_tidy <- df_tidy %>%
  #   distinct(novel_junID, .keep_all = TRUE)
  # df_tidy %>% nrow()
  
  
  
  ## GET STATS AFTER REMOVING DUPLICATES
  print(paste0("Novel donor (+) mode: ", get_mode(df_tidy %>%
                                                    distinct(novel_junID, .keep_all = TRUE) %>%
                                                    filter(type == "novel_donor") %>%
                                                    filter(distance > 0) %>%
                                                    pull(distance))))
  print(paste0("Novel donor (-) mode: ", get_mode(df_tidy %>%
                                                    distinct(novel_junID, .keep_all = TRUE) %>%
                                                    filter(type == "novel_donor") %>%
                                                    filter(distance < 0) %>%
                                                    pull(distance))))
  print(paste0("Novel acceptor (+) mode: ", get_mode(df_tidy %>%
                                                       distinct(novel_junID, .keep_all = TRUE) %>%
                                                       filter(type == "novel_acceptor") %>%
                                                       filter(distance > 0) %>%
                                                       pull(distance))))
  print(paste0("Novel acceptor (-) mode: ", get_mode(df_tidy %>%
                                                       distinct(novel_junID, .keep_all = TRUE) %>%
                                                       filter(type == "novel_acceptor") %>%
                                                       filter(distance < 0) %>%
                                                       pull(distance))))
  
  
  
  
  ## It should be more novel junctions than introns (because the same intron can be the reference for multiple novel junctions)
  df_tidy %>% nrow()
  
  if (df_tidy %>% pull(novel_junID) %>% unique() %>% length() <
      df_tidy %>% pull(ref_junID) %>% unique() %>% length()) {
    print("Error: there are less novel junctions than reference junctions")
    break;
  }
  
  
  ########################################################################################
  ## BEFORE SAVING THE FILE, we add 3 extra columns:
  ## 1. Column to store the frequency of the novel junction across samples
  ## 2. Column to store the mean counts of the novel junction across samples
  ## 3. Column to store the total number of counts of the novel junction across samples
  ########################################################################################
  
  print(stringr::str_c(Sys.time(), " - adding novel junction frequency..."))
  
  
  ## Get the count data in the current cluster for the current NOVEL junction
  split_read_counts_novel <- split_read_counts %>%
    filter(junID %in% df_tidy$novel_junID) %>%
    dplyr::select(junID, all_of(samples %>% as.character())) 

  
  split_read_counts_novel[, "novel_n_individuals"] <- matrixStats::rowCounts(split_read_counts_novel[, -c(1)] > 0, na.rm = T)
  split_read_counts_novel <- split_read_counts_novel %>% as.data.frame()
  split_read_counts_novel[, "novel_sum_counts"] <- Matrix::rowSums(split_read_counts_novel[,-c(split_read_counts_novel %>% ncol(), 1)], na.rm = T)
  split_read_counts_novel <- split_read_counts_novel %>% as.data.frame()
  
  
  if (any(split_read_counts_novel[, "novel_n_individuals"] < 1)) {
    print("Error: some novel junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_novel <- split_read_counts_novel[, c(1,(split_read_counts_novel %>% ncol() - 1), split_read_counts_novel %>% ncol())]
  
  
  ########################################################
  ## Convert data.frames to datatables to reduce timing
  ## during merge()
  ########################################################
  
  df_tidy <- df_tidy %>%
    #mutate(novel_junID = novel_junID %>% as.integer()) %>%
    data.table::as.data.table()
  
  split_read_counts_novel <- split_read_counts_novel %>%
    #mutate(junID = junID %>% as.integer()) %>% 
    data.table::as.data.table()
  
  df_merged <- merge(x = df_tidy,
                     y = split_read_counts_novel,
                     by.x = "novel_junID",
                     by.y = "junID",
                     all.x = T)
  
  
  ########################################################################################
  ## BEFORE SAVING THE FILE, we do the same for the reference introns, adding 3 extra columns:
  ## 1. Column to store the frequency of the ref junction across samples
  ## 2. Column to store the mean counts of the ref junction across samples
  ## 3. Column to store the total number of counts of the ref junction across samples
  ########################################################################################
  
  print(stringr::str_c(Sys.time(), " - adding frequency of individuals to the intron db..."))
  
  ## Get the count data in the current cluster for the current REFERENCE junction
  
  split_read_counts_ref <- split_read_counts %>%
    filter(junID %in% df_merged$ref_junID) %>%
    dplyr::select(junID, all_of(samples %>% as.character()))
  
  
  split_read_counts_ref[,"ref_n_individuals"] <- (matrixStats::rowCounts(split_read_counts_ref[, -c(1)] > 0, na.rm = T)) 
  split_read_counts_ref <- split_read_counts_ref %>% as.data.frame()
  
  split_read_counts_ref[,"ref_sum_counts"] <- Matrix::rowSums(split_read_counts_ref[,-c(split_read_counts_ref %>% ncol(),1)], na.rm = T)
  split_read_counts_ref <- split_read_counts_ref %>% as.data.frame()
  
  split_read_counts_ref <- split_read_counts_ref[, c(1,(split_read_counts_ref %>% ncol() - 1),(split_read_counts_ref %>% ncol()))]
  
  
  if (any(split_read_counts_ref[, "ref_n_individuals"] < 1)) {
    print("Error: some ref junctions do not present any read across any of the samples.")
    break;
  }
  
  split_read_counts_ref %>% head()
  
  ########################################################
  ## Convert data.frames to datatables to reduce timing
  ## during merge()
  ########################################################
  
  df_merged <- df_merged %>%
    #dplyr::mutate(ref_junID = ref_junID %>% as.integer()) %>%
    data.table::as.data.table()
  
  split_read_counts_ref <- split_read_counts_ref %>%
    #dplyr::mutate(junID = junID %>% as.integer())%>%
    data.table::as.data.table()
  
  
  ## Add data to the final df
  df_merged <- merge(x = df_merged,
                     y = split_read_counts_ref,
                     by.x = "ref_junID",
                     by.y = "junID",
                     all.x = T)
  df_merged %>% head()
  df_merged %>% nrow()
  
  ###############
  ## QC
  ###############
  
  if (intersect(df_merged$ref_junID, df_merged$novel_junID) %>% length() > 0) {
    print(paste0(Sys.time(), " - Error! Some reference junctions have been classified as novel junctions!"))
  }
  
  df_merged %>% distinct(novel_junID, .keep_all = T) %>% nrow()
  df_merged %>% distinct(ref_junID, .keep_all = T) %>% nrow()
  
  
  ###############
  ## SAVE FINAL 
  ###############
  
  saveRDS(object = df_merged %>% data.table::as.data.table(),
          file = paste0(folder_name, "/", cluster, "_distances_tidy.rds"))
  print(stringr::str_c(Sys.time(), " - file '", cluster, "_distances_tidy.rds' saved!"))
  
  
  ## RELEASE SOME MEMORY
  rm(df_all)
  rm(df_tidy)
  rm(df_merged)
  rm(df_notambig)
  rm(split_read_counts_novel)
  rm(split_read_counts_ref)
  rm(split_read_counts)
  gc()
}





## NEVER MIS-SPLICED --------------------------------


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
    
    if (file.exists(paste0(folder_name, "/", cluster, "_distances_tidy.rds"))) {
      
      df_all_misspliced <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_tidy.rds")) %>%
        dplyr::distinct(ref_junID, .keep_all = T)
      
      all_not_misspliced <- all_split_reads_details %>%
        dplyr::filter(!(junID %in% df_all_misspliced$ref_junID)) %>%
        data.table::as.data.table()
      
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
  
  
  ## Per each sample from the current cluster, we obtain all junctions, counts and ratios
  for (sample in samples) { 
    
    # sample <- samples[1]
    # sample <- samples[16]
    
    # sample <- samples[1]
    if (length(which((colnames(split_read_counts) == sample) == TRUE)) == 1) { # Every sampleID is unique, so the result should be equal to 1
      
      #distances <- list()
      split_read_counts_sample <- split_read_counts %>%
        dplyr::select(junID, all_of(sample %>% as.character())) %>%
        drop_na() 
      split_read_counts_sample <- split_read_counts_sample[split_read_counts_sample[,sample] > 0,] %>%
        data.table::as.data.table()
      
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
        filter(junID == data$junID)
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
      
      
      
      
      print(paste0(Sys.time(), " --> '", num_sample, "' samples processed. Total numbers: ", junc_not_misspliced %>% length(), 
                   " junctions not mis-spliced - ", junc_ignore %>% length(), " mis-spliced junctions."))
      
      
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


extract_never_misspliced <- function(cluster,
                                     samples,
                                     split_read_counts,
                                     folder_name) {
  
  
  df_never <- data.frame(junID = readRDS(file = paste0(folder_name, "/not-misspliced/", cluster, "_all_notmisspliced.rds")))
  print(paste0(df_never %>% nrow(), " not mis-spliced junctions in ", cluster))
  df_never %>% head()
  
  print(paste0(Sys.time(), " - obtaining reads from never mis-spliced junctions."))
  
  
  split_read_counts_never <- split_read_counts %>%
    as.data.frame() %>%
    dplyr::select(junID, all_of(samples %>% as.character())) %>%
    filter(junID %in% df_never$junID)
  
  
  
  split_read_counts_never[, "ref_n_individuals"] <- matrixStats::rowCounts(split_read_counts_never[, -1] > 0, na.rm = T)
  split_read_counts_never[, "ref_sum_counts"] <-  Matrix::rowSums(split_read_counts_never[,-c(split_read_counts_never %>% ncol(), 1)], na.rm = T)
  
  split_read_counts_never %>% head()
  
  print(paste0(Sys.time(), " - mean reads obtained for the never mis-spliced junctions."))
  
  split_read_counts_never <- split_read_counts_never %>%
    dplyr::select(-(samples %>% as.character()))
  
  split_read_counts_never %>% head()
  
  
  if (any(split_read_counts_never[, "ref_n_individuals"] < 1)) {
    print("Error: some ref junctions do not present any read across any of the samples.")
    break;
  }
  
  ## Save the not mis-spliced junctions df
  saveRDS(object = split_read_counts_never %>% data.table::as.data.table(),
          file = paste0(folder_name, "/not-misspliced/", cluster, "_all_notmisspliced_counts.rds"))
  
  print(paste0(Sys.time(), " - 'all_notmisspliced_counts.rds' file saved!"))
  
  
  rm(split_read_counts_never)
  rm(split_read_counts)
  rm(folder_name)
  rm(df_never)
  gc()
  
  
  #################################################
  #################################################
  ## FINDING THE ORIGIN OF GENE 'ENSG00000188603'
  
  # df_tidy <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_tidy.rds"))
  # df_tidy %>% head()
  # df_tidy %>% nrow()
  # 
  # 
  # df_stats <- merge(df_tidy, 
  #                   all_split_reads_details_97[,c("junID", "gene_id_start", "gene_id_end", "seqnames", "start", "end")], 
  #                   by.x = "ref_junID", 
  #                   by.y = "junID" )
  # 
  # df_stats %>% nrow()
  # df_stats %>% head()
  # 
  # df_stats <- df_stats %>%
  #   unnest(gene_id_end) %>%
  #   unnest(gene_id_start)
  # 
  # 
  # df_filtered <- df_stats %>%
  #   filter(gene_id_start == "ENSG00000188603" | gene_id_end == "ENSG00000188603") %>%
  #   as.data.frame()
  # 
  # 
  # library(xlsx)
  # ## Save the .xlsx object
  # write.xlsx(df_filtered, 
  #            file="/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/CLN3_missplicing_summary_FCTX.xlsx", 
  #            sheetName="CLN3", 
  #            row.names=T)
}


add_never_misspliced_to_df <- function(cluster,
                                       samples,
                                       all_split_reads_details,
                                       folder_name) {
  
  
  
  ## Load the not mis-spliced junctions df
  df_never <- readRDS(file = paste0(folder_name, "/not-misspliced/", cluster, "_all_notmisspliced_counts.rds")) 
  
  if (any(df_never %>% names() == "juncID")) {
    df_never <- df_never %>%
      dplyr::rename(junID = juncID)
  }
  df_never %>% head()
  df_never %>% nrow()
  
  
  ## Add features to never mis-spliced junctions
  print(paste0(Sys.time(), " - ", cluster, " adding features to never mis-spliced junctions..."))
  df_never <- merge(df_never %>% data.table::as.data.table(), 
                    all_split_reads_details[, c("junID", "ss5score","ss3score", "seqnames", "start", "end", "strand")], #%>%
                      #mutate(junID = junID %>% as.integer()) %>% data.table::as.data.table(), 
                    by.x = "junID",
                    by.y = "junID", 
                    sort = TRUE,
                    all.x = TRUE)
  df_never %>% head()
  df_never %>% nrow()
  
  ## None of the introns should be located within the X, Y or MT chromosomes
  if (any(df_never$seqnames == "MT")) {
    print("Error: some never mis-spliced introns are located on the sex and MT chromosomes!")
    break;
  }
  if (any(df_never$ss5score %>% is.na())) {
    print("Error: some never mis-spliced introns are located on the sex and MT chromosomes!")
    break;
  }
  
  ## Change name of features to match with features from the mis-splicing df
  df_never <- df_never %>% 
    dplyr::rename(ref_junID = junID,
                  ref_seq = seqnames,
                  ref_start = start,
                  ref_end = end,
                  ref_strand = strand,
                  ref_ss5score = ss5score,
                  ref_ss3score = ss3score) 
  df_never %>% head()
  df_never %>% nrow()
  df_never %>% distinct(ref_junID) %>% nrow()
  
  ######################
  
  
  ## Load the dataframe containing the mis-spliced junctions
  print(paste0(Sys.time(), " - Merging never mis-spliced and mis-spliced junctions..."))
  df_misspliced <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_tidy.rds"))
  df_misspliced %>% head()
  df_misspliced %>% nrow()
  df_misspliced %>% distinct(ref_junID) %>% nrow()
  
  
  ## None of the introns should be located within the X, Y or MT chromosomes
  if (any(df_misspliced$seqnames == "MT")) {
    print("Error: some introns are located on the sex and MT chromosomes!")
    break;
  }
  
  
  
  ## Junctions that are never misspliced shouldn't appear to be mis-spliced df
  if (df_misspliced %>%
      filter(ref_junID %in% (df_never %>% 
                             distinct(ref_junID, .keep_all = T) %>%
                             pull(ref_junID))) %>% nrow() != 0) {
    print("Error: some junctions considered as never mis-spliced have been found to be mis-spliced.")
    break;
  }
  
  ## plyr::rbind.fill fills missing columns with NA
  df_all <- plyr::rbind.fill(df_misspliced, df_never)
  df_all %>% nrow()
  df_all %>% head()
  df_all %>% distinct(ref_junID) %>% nrow()
  
  
  
  if ((df_all %>% 
       distinct(ref_junID) %>%
       nrow()) != (df_misspliced %>%
                   distinct(ref_junID) %>%
                   nrow() + 
                   df_never %>% 
                   distinct(ref_junID) %>%
                   nrow())) {
    print("Error: some junctions have been erroneusly discarded")
  }
  
  
  
  
  ## Add the rest of the features to the joined data.frame
  print(paste0(Sys.time(), " - Adding the rest of the features to the joined dataframe..."))
  df_all <- merge(df_all %>% data.table::as.data.table(), 
                  all_split_reads_details[, c("junID", "gene_id_junction", "gene_name_junction", "tx_id_junction")] %>%
                    #mutate(junID = junID %>% as.integer()) %>% 
                    data.table::as.data.table(), 
                  by.x = "ref_junID",
                  by.y = "junID", 
                  sort = TRUE,
                  all.x = TRUE)
  df_all %>% as_tibble()
  df_all %>% distinct(ref_junID) %>% nrow()
  
  
  ## None of the introns should be located within the MT chromosome
  if (any(df_all$seqnames == "MT")) {
    print("Error: some introns are located on the MT chromosome!")
    break;
  }
  
  
  ## Rename GENE columns
  
  df_all <- df_all %>% 
    dplyr::rename(gene_id = gene_id_junction) %>%
    dplyr::rename(gene_name = gene_name_junction)
  
  df_all %>% 
    filter(gene_name %in% c("RBM5","RBM6")) %>% 
    select(gene_name,gene_id) %>%
    unnest(gene_id) %>%
    unnest(gene_name) %>%
    mutate_if(is.list, simplify_all)
  ######################
  ## QC
  ###################### 
  
  df_all %>% distinct(ref_junID) %>% nrow()
  df_all %>% distinct(novel_junID) %>% nrow()
  
  
  df_all %>%
    filter(ref_junID %in% df_never$ref_junID)
  
  
  all_misspliced_not_paired <- readRDS(file = paste0(folder_name, "/not-misspliced/", cluster, "_all_misspliced_not_paired.rds"))
  if (intersect(c(df_all$ref_junID, df_all$novel_junID), all_misspliced_not_paired) %>% length() > 0) {
    print("Error: the dataframe of mis-spliced junctions contain not-paired junctions!")
  }
  
  
  ######################
  ## SAVE RESULTS
  ######################
  
  saveRDS(object = df_all %>% data.table::as.data.table(),
          file = paste0(folder_name, "/", cluster, "_distances_tidy_all.rds"))
  print(stringr::str_c(Sys.time(), " - file '", cluster, "_distances_tidy_all.rds' saved!"))
  
  
  rm(all_split_reads_details)
  rm(df_misspliced)
  rm(df_never)
  rm(df_all)
  gc()
}




## MIS-SPLICING RATIO   --------------------------------


get_missplicing_ratio <- function(cluster,
                                  samples,
                                  split_read_counts,
                                  folder_name,
                                  folder_save_name,
                                  db = NULL) {
  
  ## Load source file containing the distances to the reference intron
  if (file.exists(paste0(folder_name, "/", cluster, "_distances_tidy_all.rds"))) {
    
    df_tidy <- readRDS(file = paste0(folder_name, "/", cluster, "_distances_tidy_all.rds"))
    
  } else {
    
    print("Error: 'distances' file does not exist!")
    break;
    
  }
  
  
  print(paste0(Sys.time(), " - Calulating MSR_Donor and MSR_Acceptor."))
  
  
  df_tidy %>% head()
  
  # SNCA David's novel exon: 
  # df_tidy %>% filter(novel_junID == "67198600")
  
  print(paste0(Sys.time(), " --> ", cluster, ": ", df_tidy %>%
                 distinct(ref_junID) %>%
                 nrow(), " reference junctions."))
  
  # df_tidy %>%
  #   filter(is.na(ref_counts)) %>%
  #   nrow()
  
  # df_tidy %>%
  #   filter(gene_name_start == "SNCA") %>%
  #   select(ref_junID, ref_n_individuals, ref_mean_counts, ref_sum_counts) %>%
  #   distinct(ref_junID, .keep_all = T)
  
  ## Threshold for the minimum number of reads
  # threshold <- 5
  
  
  ## 50% of the total number of samples must present a mean of at least 5 reads
  # print(paste0("Threshold samples: ", round(samples %>% length() * 0.5)))
  # print(paste0("Threshold counts: ", threshold)) #(round(samples %>% length() * 0.5) * threshold)  / samples %>% length()))
  #print(paste0("Threshold counts: ", threshold))
  
  
  ## Obtain only reference junction expressed in 50% of samples presenting at least 5 reads
  # df <- df_tidy %>%
  #   filter(#ref_n_individuals >= round(samples %>% length() * 0.5), 
  #     #ref_mean_counts >= threshold, #(round(samples %>% length() * 0.5) * threshold)  / samples %>% length(),
  #     !is.na(ref_counts))
  
  ## Get only mis-spliced reference introns and novel junctions
  df <- df_tidy %>%
    dplyr::filter(!is.na(ref_counts), 
                  ref_counts > 0)
  
  # df %>%
  #   filter(gene_name_start == "SNCA")%>%
  #   distinct(ref_junID, .keep_all = T)
  
  ## Get novel junctions ---------------------------------------------------------------------------------------------------------------------------
  
  db_novel <- df %>% 
    as_tibble() %>%
    mutate(novel_start = novel_start %>% as.integer(),
           novel_end = novel_end %>% as.integer()) %>%
    dplyr::select(ref_junID,
                  novel_junID,
                  novel_type = type,
                  distance,
                  seqnames = novel_seq,
                  start = novel_start,
                  end = novel_end,
                  strand = novel_strand,
                  novel_ss5score,
                  novel_ss3score,
                  novel_n_individuals,
                  novel_sum_counts,
                  #novel_missplicing_ratio_sample,
                  #novel_missplicing_ratio_tissue,
                  gene_name,
                  gene_id) %>%
    distinct(novel_junID, .keep_all = T)
  
  # db_novel %>%
  #   filter(novel_junID == "67198600")
  # db_novel %>%
  #   filter(gene_name == "SNCA")%>%
  #   distinct(ref_junID, .keep_all = T)
  
  # db_novel %>% filter(gene_name == "SNCA") %>% distinct(ref_junID) %>% nrow
  
  # db_novel %>% distinct(ref_junID) %>% nrow()
  # db_novel %>% nrow()
  # db_novel %>% distinct(novel_junID) %>% nrow()
  
  saveRDS(object = db_novel, file = paste0(folder_save_name, "/", cluster, "_db_novel.rds"))
  
  print(paste0(Sys.time(), " - Novel junctions DB saved!"))
  
  
  
  ## Get novel junctions details ----------------------------------------------------------------------------
  
  # split_read_counts_novel <- split_read_counts %>%
  #   filter(junID %in% (db_novel$novel_junID %>% unique())) %>%
  #   dplyr::select(junID, samples %>% as.character()) %>%
  #   #mutate_all(~replace(., is.na(.), 0)) %>%
  #   as_tibble()
  # 
  # split_read_counts_novel_details <- split_read_counts_novel %>% #head() %>%
  #   gather(sample, count, samples %>% as.character()) %>%
  #   dplyr::rename(novel_junID = junID,
  #                 novel_counts = count) %>%
  #   dplyr::filter(!is.na(novel_counts))
  # 
  # 
  # split_read_counts_novel_details %>%
  #   as.data.frame() %>%
  #   filter(novel_junID == "1001789")
  # 
  # 
  # 
  # for (sample in split_read_counts_novel_details$sample %>% unique()) { # sample <- db_novel_details$sample[1]
  #   
  #   if (!is.null(db)) {
  #     
  #     ind <- which(db$run %>% as.character() == sample)
  #     people <- db[ind, ]
  #     
  #     if (str_detect(people$characteristics[1], "female")) {
  #       sex <- "female"
  #     }else if (str_detect(people$characteristics[1], "male")) {
  #       sex <- "male"
  #     }
  #     
  #     data <- people$characteristics[1]
  #     position <- str_locate(data, "age at death")
  #     age <- substr(data, position[2] + 3, position[2] + 4) %>% 
  #       as.integer()
  #     
  #     mapped_read_count <- people$mapped_read_count
  #     avg_read_length <- people$avg_read_length
  #     
  #   } else {
  #     
  #     con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/data/splicing_tolerance/splicing_tolerance.sqlite")
  #     query <- DBI::dbSendQuery(con, paste0("SELECT * FROM GTEX_info WHERE sample_recount_id == '", sample,"';")) #  WHERE smafrze == 'USE ME'
  #     people <- DBI::dbFetch(query)
  #     DBI::dbClearResult(query)
  #     DBI::dbDisconnect(con)
  #     
  #     if (str_detect(people$bigwig_file[1], "female")) {
  #       sex <- "female"
  #     } else if (str_detect(people$bigwig_file[1], "male")) {
  #       sex <- "male"
  #     }
  #     
  #     age <- people$age
  #     mapped_read_count <- people$mapped_read_count
  #     avg_read_length <- people$avg_read_length
  #   }
  #   
  #   
  #   
  #   ind <- which(split_read_counts_novel_details$sample == sample)
  #   split_read_counts_novel_details[ind,"age"] <- age
  #   split_read_counts_novel_details[ind,"sex"] <- sex
  #   split_read_counts_novel_details[ind,"mapped_read_count"] <- mapped_read_count
  #   split_read_counts_novel_details[ind,"avg_read_length"] <- avg_read_length
  #   
  #   # print(paste0(sample)) 
  # }
  # 
  # saveRDS(object = split_read_counts_novel_details,
  #         file = paste0(folder_save_name, "/", cluster, "_db_novel_details.rds"))
  # 
  # print(paste0(Sys.time(), " - Details of novel junctions DB saved!"))
  
  
  
  ## Generate reference introns data base ---------------------------------------------------------------------------
  
  
  
  
  ## The novel junction itself is not important anymore. 
  ## We calculate the mis-splicing ratio per novel event and reference intron
  df_ref <- df %>% 
    #filter(ref_junID == "1000209") %>%
    group_by(ref_junID, type) %>%
    mutate(novel_mean_counts = novel_sum_counts / novel_n_individuals) %>%
    mutate(ref_mean_counts = ref_sum_counts / ref_n_individuals) %>%
    dplyr::mutate(ref_missplicing_ratio_tissue = sum(novel_mean_counts) / 
                    (sum(novel_mean_counts) + ref_mean_counts)) %>% 
    #as.data.frame()
    #dplyr::mutate(ref_missplicing_ratio_tissue2 = sum(novel_mean_counts)/(sum(novel_mean_counts) + (ref_mean_counts))) %>%
    distinct(ref_missplicing_ratio_tissue, .keep_all = T) %>%
    ungroup() #%>% select(ref_missplicing_ratio_tissue)
  
  
  df_ref_tidy <- df_ref %>% 
    spread(key = type, value = ref_missplicing_ratio_tissue) %>%
    group_by(ref_junID) %>% 
    dplyr::mutate(ref_missplicing_ratio_tissue_ND = mean(novel_donor, na.rm = T)) %>%
    dplyr::mutate(ref_missplicing_ratio_tissue_NA = mean(novel_acceptor, na.rm = T))  %>%
    distinct(ref_junID, .keep_all = T)  %>%
    ungroup()
  
  df_ref_tidy[1,]
  #df_ref_tidy %>% filter(ref_junID == "10086104") %>% as.data.frame()
  #df_ref %>% filter(ref_junID == "10086104") %>% as.data.frame()
  
  ## For the intron database, we only keep columns related to the intron
  db_ref <- df_ref_tidy %>%
    dplyr::select(-novel_junID, 
                  -distance,
                  -novel_counts,
                  -novel_seq,
                  -novel_start,
                  -novel_end,
                  -novel_strand,
                  -novel_ss5score,
                  -novel_ss3score,
                  -novel_n_individuals,
                  -novel_mean_counts,
                  -novel_sum_counts,
                  -ref_counts,
                  -ref_mean_counts,
                  #-gene_name_end,
                  -novel_donor,
                  -novel_acceptor) %>%
    dplyr::rename(seqnames = ref_seq,
                  start = ref_start,
                  end = ref_end,
                  strand = ref_strand) %>%
    as.data.frame()
  
  db_ref %>% nrow()
  db_ref %>% distinct(ref_junID) %>% nrow()
  
  #db_ref %>% filter(ref_junID == "10086104") %>% as.data.frame()
  
  ## QC
  if (!identical(db_ref %>% distinct(ref_junID) %>% arrange(desc(ref_junID)) %>% pull(),
                 db_novel %>% distinct(ref_junID) %>% arrange(desc(ref_junID)) %>% pull())) {
    print("Error: novel junctions database doesn't contain all reference introns!")
    break;
  }
  
  ## Generate never mis-spliced reference introns ---------------------------------------------------------------
  
  df_never <- df_tidy %>% 
    dplyr::filter(is.na(ref_counts))#, 
  #ref_n_individuals >= round(samples %>% length() * 0.5),
  #ref_mean_counts >= threshold)#(round(samples %>% length() * 0.5) * threshold)  / samples %>% length())
  
  db_never <- df_never %>% 
    dplyr::select(-novel_junID, 
                  -type,
                  #-sample, 
                  -distance,
                  -novel_counts,
                  -novel_seq,
                  -novel_start,
                  -novel_end,
                  -novel_strand,
                  #-novel_width,
                  -novel_ss5score,
                  -novel_ss3score,
                  -novel_n_individuals,
                  #-novel_mean_counts,
                  -novel_sum_counts,
                  -ref_counts) %>%
    dplyr::rename(seqnames = ref_seq,
                  start = ref_start,
                  end = ref_end,
                  strand = ref_strand) %>%
    as.data.frame()
  
  
  
  
  db_ref_joined <- plyr::rbind.fill(db_ref, db_never)
  db_ref_joined %>% nrow()
  
  db_ref_joined %>% distinct(ref_junID) %>% nrow()
  db_ref_joined %>% head()
  
  
  
  print(paste0(Sys.time(), " - Tidying up results..."))
  
  if (db_ref_joined %>%
      dplyr::filter(is.na(ref_missplicing_ratio_tissue_ND),
                    is.na(ref_missplicing_ratio_tissue_NA)) %>%
      distinct(ref_junID) %>%
      nrow() != df_never %>%
      #filter(is.na(novel_junID)) %>%
      distinct(ref_junID) %>%
      nrow()) {
    print("Error: disimilar number of never mis-spliced junctions")
  }
  
  
  db_ref_joined[is.na(db_ref_joined$ref_missplicing_ratio_tissue_ND),"ref_missplicing_ratio_tissue_ND"] <- 0
  db_ref_joined[is.na(db_ref_joined$ref_missplicing_ratio_tissue_NA),"ref_missplicing_ratio_tissue_NA"] <- 0
  db_ref_joined %>% head()
  
  
  db_ref_joined %>% nrow()
  db_ref_joined %>% distinct(ref_junID) %>% nrow()
  
  
  db_ref_joined$ref_missplicing_ratio_tissue_ND %>% summary()
  db_ref_joined$ref_missplicing_ratio_tissue_NA %>% summary()
  
  
  print(paste0(Sys.time(), " - ", db_ref_joined %>% nrow(), " final number of introns."))
  
  
  # db_ref_joined <- readRDS(file = paste0(folder_save_name, "/v104/", cluster, "_db_introns.rds"))
  
  saveRDS(object = db_ref_joined, file = paste0(folder_save_name, "/", cluster, "_db_introns.rds"))
  
  
  print(paste0(Sys.time(), " - Introns DB saved!"))
  
  
  rm(db_ref_joined)
  rm(db_never)
  rm(df_never)
  rm(db_ref)
  rm(db_novel)
  rm(df)
  rm(df_tidy)
  gc()
  
  ## Get reference introns details ----------------------------------------------------------------------------
  
  # split_read_counts_ref <- split_read_counts %>%
  #   filter(junID %in% (db_ref_joined$ref_junID %>% unique())) %>%
  #   dplyr::select(junID, samples %>% as.character()) %>%
  #   #mutate_all(~replace(., is.na(.), 0)) %>%
  #   as_tibble()
  # 
  # 
  # split_read_counts_ref_details <- split_read_counts_ref %>% 
  #   gather(sample, count, samples %>% as.character()) %>%
  #   dplyr::rename(ref_junID = junID) %>%
  #   dplyr::filter(!is.na(count))
  # 
  # 
  # 
  # for (sample in split_read_counts_ref_details$sample %>% unique()) { # sample <- split_read_counts_ref_details$sample[1]
  #   
  #   if (!is.null(db)) {
  #     
  #     ind <- which(db$run %>% as.character() == sample)
  #     people <- db[ind, ]
  #     
  #     if (str_detect(people$characteristics[1], "female")) {
  #       sex <- "female"
  #     }else if (str_detect(people$characteristics[1], "male")) {
  #       sex <- "male"
  #     }
  #     
  #     data <- people$characteristics[1]
  #     position <- str_locate(data, "age at death")
  #     age <- substr(data, position[2] + 3, position[2] + 4) %>% 
  #       as.integer()
  #     
  #     mapped_read_count <- people$mapped_read_count
  #     avg_read_length <- people$avg_read_length
  #     
  #   } else {
  #     
  #     con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "/data/splicing_tolerance/splicing_tolerance.sqlite")
  #     query <- DBI::dbSendQuery(con, paste0("SELECT * FROM GTEX_info WHERE sample_recount_id == '", sample,"';")) #  WHERE smafrze == 'USE ME'
  #     people <- DBI::dbFetch(query)
  #     DBI::dbClearResult(query)
  #     DBI::dbDisconnect(con)
  #     
  #     if (str_detect(people$bigwig_file[1], "female")) {
  #       sex <- "female"
  #     }else if (str_detect(people$bigwig_file[1], "male")) {
  #       sex <- "male"
  #     }
  #     
  #     age <- people$age
  #     mapped_read_count <- people$mapped_read_count
  #     avg_read_length <- people$avg_read_length
  #   }
  #   
  #   
  #   ind <- which(split_read_counts_ref_details$sample == sample)
  #   split_read_counts_ref_details[ind,"age"] <- age
  #   split_read_counts_ref_details[ind,"sex"] <- sex
  #   split_read_counts_ref_details[ind,"mapped_read_count"] <- mapped_read_count
  #   split_read_counts_ref_details[ind,"avg_read_length"] <- avg_read_length
  #   
  #   #split_read_counts_ref_details
  #   #print(sample)
  # }
  # 
  # 
  # 
  # saveRDS(object = split_read_counts_ref_details,
  #         file = paste0(folder_save_name, "/", cluster, "_db_introns_details.rds"))
  # 
  # 
  # print(paste0(Sys.time(), " - Intron details DB saved!"))
  
}




add_missplicing_class_to_df <- function(cluster,
                                        folder_name) {
  
  
  df_missplicing <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  df_missplicing %>% nrow()
  df_missplicing %>% head()
  
  
  
  ## Intron classification --------------------------------------------------------------------------------
  
  print(paste0(Sys.time(), " - ", cluster, ", starting junction classification!"))
  
  # print(paste0("All junctions: ", df_missplicing %>% nrow()))
  # print(paste0("All genes: ", df_missplicing %>% 
  #                dplyr::select(c(gene_id)) %>%
  #                unname() %>% unlist() %>% unique() %>% length()))
  
  
  df_missplicing[, "ref_type"] <- ""
  df_missplicing %>% head()
  
  
  
  ## TYPE 'BOTH'
  indx <- which(df_missplicing$ref_missplicing_ratio_tissue_ND > 0 & df_missplicing$ref_missplicing_ratio_tissue_NA > 0)
  df_missplicing[indx, "ref_type"] <- "both"
  print(paste0("Introns type 'both': ", indx %>% length()))
  
  # print(paste0("Genes type 'both': ", df_missplicing[indx,] %>% 
  #                dplyr::select(c(gene_id)) %>%
  #                unname() %>% unlist() %>% unique() %>% length()))
  
  
  
  ## TYPE 'DONOR'
  indx <- which(df_missplicing$ref_missplicing_ratio_tissue_ND > 0 & df_missplicing$ref_missplicing_ratio_tissue_NA == 0)
  df_missplicing[indx, "ref_type"] <- "donor"
  print(paste0("Junctions type 'donor': ", indx %>% length()))
  
  # print(paste0("Genes type 'donor': ", df_missplicing[indx,] %>% 
  #                dplyr::select(c(gene_id)) %>%
  #                unname() %>% unlist() %>% unique() %>% length()))
  
  
  
  ## TYPE 'ACCEPTOR'
  indx <- which(df_missplicing$ref_missplicing_ratio_tissue_ND == 0 & df_missplicing$ref_missplicing_ratio_tissue_NA > 0)
  df_missplicing[indx, "ref_type"] <- "acceptor"
  print(paste0("Junctions type 'acceptor': ", indx %>% length()))
  
  # print(paste0("Genes type 'acceptor': ", df_missplicing[indx,] %>% 
  #                dplyr::select(c(gene_id)) %>%
  #                unname() %>% unlist() %>% unique() %>% length()))
  
  
  
  ## TYPE 'NONE'
  indx <- which(df_missplicing$ref_missplicing_ratio_tissue_ND == 0 & df_missplicing$ref_missplicing_ratio_tissue_NA == 0)
  df_missplicing[indx, "ref_type"] <- "never"
  print(paste0("Junctions type 'never': ", indx %>% length()))
  
  # print(paste0("Type 'never': ", df_missplicing[indx,] %>% 
  #                dplyr::select(c(gene_id)) %>%
  #                unname() %>% unlist() %>% unique() %>% length()))
  
  
  
  
  df_missplicing$ref_type %>% unique()
  
  
  ## QC --------------------------------------------------------
  if (intersect(df_missplicing %>%
                filter(ref_type == "never") %>%
                distinct(ref_junID) %>%
                pull(), df_missplicing %>%
                filter(ref_type == "donor") %>%
                distinct(ref_junID) %>%
                pull()) %>% length() > 0) {
    print("Error!")
  }
  
  
  if (intersect(df_missplicing %>%
                filter(ref_type == "never") %>%
                distinct(ref_junID) %>%
                pull(), df_missplicing %>%
                filter(ref_type == "aceptor") %>%
                distinct(ref_junID) %>%
                pull()) %>% length() > 0) {
    print("Error!")
  }
  
  if (intersect(df_missplicing %>%
                filter(ref_type == "never") %>%
                distinct(ref_junID) %>%
                pull(), df_missplicing %>%
                filter(ref_type == "both") %>%
                distinct(ref_junID) %>%
                pull()) %>% length() > 0) {
    print("Error!")
  }
  
  if (intersect(df_missplicing %>%
                filter(ref_type == "acceptor") %>%
                distinct(ref_junID) %>%
                pull(), df_missplicing %>%
                filter(ref_type == "donor") %>%
                distinct(ref_junID) %>%
                pull()) %>% length() > 0) {
    print("Error!")
  }
  
  df_missplicing %>% nrow()
  df_missplicing %>% distinct(ref_junID) %>% nrow()
  
  saveRDS(object = df_missplicing,
          file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  
  print(paste0(Sys.time(), " - results saved!")) 
  
  
  # ggplot(data = df_missplicing ) + 
  #   geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"), 
  #                alpha = 0.8) +
  #   geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
  #                alpha = 0.8) +
  #   
  #   #ggtitle(title) +
  #   xlab("Mis-splicing ratio mean") +
  #   #facet_wrap(~ type) +
  #   #ylab("Intron count") +
  #   #scale_y_continuous(limits = c(0, 85000)) +
  #   #xlim(c(0, 0.01)) +
  #   
  #   facet_zoom(xlim = c(0, 0.01)) + 
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
    filter(junID %in% (df_introns %>%
                         filter(ref_type == "never") %>%
                         distinct(ref_junID) %>%
                         pull()))
  
  
  ## Add the rest of the novel donor and acceptor junctions
  all_not_misspliced <- rbind(all_not_misspliced,
                              all_split_reads_details %>%
                                filter(type %in% c("novel_donor", "novel_acceptor")))
  
  
  
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
      filter(ref_type != "never") %>% 
      distinct(ref_junID) %>%
      nrow() != df_novel %>%
      distinct(ref_junID) %>%
      nrow()) {
    print("Error: some novel junctions aren't attached to the same introns stored within the intron database.")
    break;
  }
  
  
  
  if (intersect(df_novel %>%
                distinct(ref_junID) %>%
                pull(), df_introns %>%
                filter(ref_type != "never") %>%
                distinct(ref_junID) %>%
                pull()) %>% length() != df_introns %>% 
      filter(ref_type != "never") %>% nrow()){
    print(paste0("Error: some ", cluster, " novel junctions are attached to introns that are not stored on the database!"))
  }
  
  if (intersect(df_novel %>%
                distinct(novel_junID) %>%
                pull() %>% sort(), all_split_reads_details %>%
                filter(type %in% c("novel_donor", "novel_acceptor")) %>%
                distinct(junID) %>%
                pull() %>% sort()) %>% length() != df_novel %>% nrow()){
    
    
    diff <-setdiff(df_novel %>%
                     distinct(novel_junID) %>%
                     pull(), all_split_reads_details %>%
                     filter(type %in% c("novel_donor", "novel_acceptor")) %>%
                     distinct(junID) %>%
                     pull())
    
    all_split_reads_details %>%
      filter(junID == diff[1])
    
    
    all_split_reads_details %>%
      filter(junID == "JUNC00006215")
    
    print(paste0("Error: some ", cluster, " novel junctions have not been found within the original annotation!"))
  }
  
  
  
  
  ## FREE SOME MEMORY
  
  rm(never_misspliced_pairings)
  rm(all_not_misspliced)
  rm(df_introns)
  rm(df_novel)
  gc()
  
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
      filter(!(gene_id %in% MT_geneID)) 
    
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


## INTRON TYPE   --------------------------------

add_intron_type <- function(cluster,
                            folder_name) {
  
  
  ## Load intron type files
  u12_introns <- readRDS(file = "/data/references/IAOD/minor_introns_tidy.rds")
  u2_introns <- readRDS(file = "/data/references/IAOD/major_introns_tidy.rds")
  
  
  ## Load intron database
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>% as_tibble()
  df_introns %>% head()
  df_introns %>% nrow()
  
  # df_introns <- df_introns %>%
  #   rename(transcript_id_start = transcript_id_start.x,
  #          transcript_id_end = transcript_id_end.x,
  #          protein_coding = protein_coding.x,
  #          gene_biotype = gene_biotype.x) %>%
  #   select(-transcript_id_start.y,-transcript_id_end.y,-protein_coding.y,-gene_biotype.y)
  # 
  # saveRDS(object = df_introns,
  #         file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  
  ## Add two new columns to incorporate info about intron type
  df_introns <- df_introns %>%
    mutate(u12_intron = F) %>%
    mutate(u2_intron = F)
  
  df_introns %>% head()
  df_introns %>% nrow()
  df_introns %>% distinct(ref_junID) %>% nrow()
  #df_introns %>% distinct(novel_junID) %>% nrow()
  
  ################
  ## MINOR INTRON
  ################
  
  print(paste0(Sys.time(), " - Getting junctions spliced out by the minor spliceosome."))
  
  overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u12_introns$seqnames,
                                                                         ranges = IRanges(start = u12_introns$start, 
                                                                                          end = u12_introns$end),
                                                                         strand = u12_introns$strand),
                                          subject = GenomicRanges::GRanges(seqnames = df_introns$seqnames,
                                                                           ranges = IRanges(start = df_introns$start, 
                                                                                            end = df_introns$end),
                                                                           strand = df_introns$strand),
                                          ignore.strand = FALSE,
                                          type = "equal")
  
  queryHits(overlaps) %>% length()
  df_introns[subjectHits(overlaps),]$u12_intron <- T
  
  
  
  ################
  ## MAJOR INTRON
  ################
  
  print(paste0(Sys.time(), " - Getting junctions spliced out by the major spliceosome."))
  
  overlaps <- GenomicRanges::findOverlaps(query = GenomicRanges::GRanges(seqnames = u2_introns$seqnames %>% as.character(),
                                                                         ranges = IRanges(start = u2_introns$start,
                                                                                          end = u2_introns$end),
                                                                         strand = u2_introns$strand),
                                          subject = GenomicRanges::GRanges(seqnames = df_introns$seqnames,
                                                                           ranges = IRanges(start = df_introns$start,
                                                                                            end = df_introns$end),
                                                                           strand = df_introns$strand),
                                          ignore.strand = FALSE,
                                          type = "equal")
  
  # print(overlaps)
  subjectHits(overlaps) %>% length()
  df_introns[subjectHits(overlaps),]$u2_intron <- T
  
  
  # ################
  # ## QC
  # ################
  # 
  # df_introns %>% head
  # df_introns %>% nrow
  # 
  # 
  # 
  # ## Number of minor introns
  # df %>%
  #   filter(u12_intron == T) %>%
  #   nrow()
  # 
  # ## Number of major introns
  # df %>%
  #   filter(u2_intron == T) %>%
  #   nrow()
  # 
  # 
  # ## Number of introns from SNCA
  # df %>%
  #   filter(gene_name_start == "SNCA" | gene_name_end == "SNCA") %>%
  #   distinct(ref_junID) %>%
  #   nrow()
  # 
  # ## Number of uncatalogued junctions
  # 
  # # print(paste0(Sys.time(), " - removing uncatalogued junctions..."))
  # 
  # # if (df %>%
  # #     filter(u2_intron == F, u12_intron == F) %>%
  # #     nrow() > 0) {
  # # 
  # #   ## Remove uncatalogued junctions
  # #   df <- df %>%
  # #     filter(u12_intron == T | u2_intron == T)
  # # }
  # 
  # 
  # ## Number of introns from SNCA
  # df %>%
  #   filter(gene_name_start == "SNCA" | gene_name_end == "SNCA") %>%
  #   distinct(ref_junID) %>%
  #   nrow()
  # 
  # 
  # df %>% nrow()
  
  
  df_introns <- df_introns %>%
    filter(u2_intron == T)
  
  ## Matching results with novel introns
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))%>%
    filter(ref_junID %in% df_introns$ref_junID)
  
  
  # df_introns %>% head()
  # df_introns %>% nrow()
  # df_novel %>% nrow()
  
  if (!identical(df_novel$ref_junID %>% unique(), 
                 df_introns %>% filter(ref_type != "never") %>% select(ref_junID) %>% pull()%>% unique())) {
    print("ERROR!")
  }
  
  ################
  ## SAVE RESULTS
  ################
  
  # df_introns %>% nrow()
  # df_introns %>% distinct(ref_junID) %>% nrow()
  
  saveRDS(object = df_introns,
          file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  
  saveRDS(object = df_novel,
          file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  
  
  print(paste0(Sys.time(), " - results saved!"))
  
  rm(df_introns)
  rm(u12_introns)
  rm(u2_introns)
  rm(overlaps)
  gc()
}


## CLINVAR ANALYSIS  -----------------------------

clinvar_analysis <- function(cluster,
                             folder_name, 
                             plot = F) {
  
  
  print(paste0(Sys.time(), " - loading the ClinVar source file..."))
  
  ## Load clinvar data
  clinvar_tidy <- readRDS(file = "/data/references/clinvar/clinvar_intronic_tidy.rda") %>%
    GRanges() %>%
    diffloop::rmchr()
  clinvar_tidy %>% head()
  
  
  
  ## Load mis-splicing ratios data
  df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
    mutate(clinvar_type = "-") %>%
    mutate(clinvar_start = 0) %>%
    mutate(clinvar_end = 0) %>%
    GRanges() 
  df_introns %>% head()
  df_introns %>% length()
  
  
  
  ## Find overlaps between clinvar mutations and mis-splicing ratios
  overlaps <- GenomicRanges::findOverlaps(query = clinvar_tidy,
                                          subject = df_introns,
                                          ignore.strand = F,
                                          type = "any")
  
  
  elementMetadata(df_introns)[subjectHits(overlaps), "clinvar_type"] <- ifelse(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"] == "none",
                                                                               "-", paste0(elementMetadata(clinvar_tidy)[queryHits(overlaps), "splice_site"], "-",
                                                                                           elementMetadata(clinvar_tidy)[queryHits(overlaps), "distance"], "bp"))
  elementMetadata(df_introns)[subjectHits(overlaps), "clinvar_start"] <- clinvar_tidy[queryHits(overlaps), ] %>% start()
  elementMetadata(df_introns)[subjectHits(overlaps), "clinvar_end"] <- clinvar_tidy[queryHits(overlaps), ] %>% end()
  
  
  df_introns %>%
    as.data.frame() %>%
    filter(clinvar_type != "-")
  
  # df_introns[subjectHits(overlaps),] %>% head()
  # 
  # # print(paste0(Sys.time(), " - locating ClinVar mutations within annotated introns..."))
  # 
  # # for (i in unique(query)) { # i <- unique(query)[1]
  # #   
  # #   # print(i)
  # #   
  # #   clinvar_mut <- clinvar_tidy[i,]
  # #   
  # #   df_introns[hits[[i]],]$clinvar <- 
  # #     case_when(clinvar_mut$closer_start_end == "start" & clinvar_mut$strand == "+" ~ "acceptor", 
  # #               clinvar_mut$closer_start_end == "start" & clinvar_mut$strand == "-" ~ "donor", 
  # #               clinvar_mut$closer_start_end == "end" & clinvar_mut$strand == "+" ~ "donor", 
  # #               clinvar_mut$closer_start_end == "end" & clinvar_mut$strand == "-" ~ "acceptor")
  # # }
  # 
  # df_introns %>% head()
  # df_introns %>% length()
  
  
  
  # which(df_introns$clinvar == "donor") %>% length()
  # which(df_introns$clinvar == "acceptor") %>% length()
  
  
  
  ## Getting introns non containing clinvar mutations
  
  if (plot) {
    
    df_not_clinvar <- df_introns %>%
      as.data.frame() %>%
      filter(clinvar_type == "none")
    
    df_clinvar <- df_introns %>%
      as.data.frame() %>%
      filter(clinvar_type != "none")
    
    if ((intersect(df_not_clinvar$ref_junID, df_clinvar$ref_junID) %>% length()) > 0) {
      print("Error: some introns are classified as both having and non-having clinvar mutations.")
    }
    
    
    
    ## Wilcoxon test
    
    wilcox.test(x = df_clinvar %>% filter(clinvar_type == "donor") %>% pull(ref_missplicing_ratio_tissue_ND),
                y = df_not_clinvar$ref_missplicing_ratio_tissue_ND,
                paired = F,
                correct = T,
                alternative = "less")
    
    
    ## Plot results
    
    ggplot() +
      geom_density(aes(x = df_clinvar %>% filter(clinvar_type == "donor") %>% pull(ref_missplicing_ratio_tissue_ND), fill = "#440154FF"),
                   alpha = 0.8) +
      geom_density(aes(x = df_not_clinvar %>% pull(ref_missplicing_ratio_tissue_ND), fill = "#35B779FF"),
                   alpha = 0.8) +
      theme_light() +
      xlim(c(0,1))+
      #ylab("Rumber of reference junctions containing at least one ClinVar mutation") +
      xlab("Mean mis-splicing ratio at donor splice sites") +
      scale_fill_manual(values = c("#440154FF", "#35B779FF"),
                        breaks = c("#440154FF", "#35B779FF"),
                        labels = c("ClinVar","Non-ClinVar")) +
      ggtitle(paste0("Mean mis-splicing ratio at donor splice sites.")) +
      theme(axis.line = element_line(colour = "black"),
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12"),
            legend.text = element_text(size = "12"),
            plot.caption = element_text(size = "12"),
            #legend.title = element_text(size = "12"),
            legend.position = "top") +
      guides(fill = guide_legend(title = NULL,
                                 ncol = 2,
                                 nrow = 1))
    
    
    
    ## Save the results
    folder_img <- paste0(folder_name, "/images/")
    dir.create(file.path(folder_img), showWarnings = T)
    file_name <- paste0(folder_img, cluster, "_missplicingratio_donor_clinvar.png")
    ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
    
    
    
    
    ## TEST NOVEL ACCEPTOR
    
    wilcox.test(x = df_clinvar %>% filter(clinvar == "acceptor") %>% pull(ref_missplicing_ratio_tissue_NA),
                y = df_not_clinvar$ref_missplicing_ratio_tissue_NA,
                paired = F,
                correct = T,
                alternative = "less")
    
    
    
    ggplot() + 
      geom_density(aes(x = df_clinvar %>% filter(type == "acceptor") %>% pull(ref_missplicing_ratio_tissue_NA), fill = "#440154FF"), 
                   alpha = 0.8) +
      geom_density(aes(x = df_not_clinvar$ref_missplicing_ratio_tissue_NA, fill = "#35B779FF"), 
                   alpha = 0.8) +
      theme_light() +
      xlim(c(0,1))+
      #ylab("Rumber of reference junctions containing at least one ClinVar mutation") +
      xlab("Mean mis-splicing ratio at acceptor splice sites") + 
      scale_fill_manual(values = c("#440154FF", "#35B779FF"),
                        breaks = c("#440154FF", "#35B779FF"),
                        labels = c("ClinVar","Non-ClinVar")) +
      ggtitle(paste0("Mean mis-splicing ratio at acceptor splice sites.")) +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12"),
            legend.text = element_text(size = "12"),
            plot.caption = element_text(size = "12"),
            #legend.title = element_text(size = "12"),
            legend.position = "top") +
      guides(fill = guide_legend(title = NULL,
                                 ncol = 2, 
                                 nrow = 1))
    
    
    ## Save the results
    folder_img <- paste0(folder_name, "/images/")
    file_name <- paste0(folder_img, cluster, "_missplicingratio_acceptor_clinvar.png")
    ggplot2::ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
    
  }
  
  #############################
  ## SAVE RESULTS
  #############################
  
  saveRDS(object = df_introns,
          file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  
  print(paste0(Sys.time(), " - results saved!"))
  
  rm(clinvar_tidy)
  rm(df_introns)
  rm(overlaps)
  gc()
  
}


## MANE TRANSCRIPTS ------------------------------

# folder_name <- "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/BRAIN//results/pipeline3/missplicing-ratio/Brain - Hippocampus/v105/"
# cluster <- "Brain - Hippocampus"

add_MANE_info <- function(cluster,
                          folder_name) {
  
  print(paste0(Sys.time(), " - Adding MANE information to '", cluster, "' reference introns!"))
  
  ## Import MANE annotation
  mane_transcripts <- rtracklayer::import(con = "/data/references/MANE/MANE.GRCh38.v1.0.ensembl_genomic.gtf") %>% 
    as_tibble() %>%
    filter(tag == "MANE_Select" | tag == "MANE_Plus_Clinical") %>%
    mutate(transcript_id = str_remove(string = transcript_id, pattern = "\\..*") )# %>%
    #distinct(transcript_id)
  
  ## Load intron data from the IDB
  IDB_path <- paste0(folder_name, "/", cluster, "_db_introns.rds")
  df_introns <- readRDS(file = IDB_path) 
  
  df_introns_mane <- df_introns %>% 
    as_tibble() %>%
    rowwise() %>%
    mutate(MANE_2 = ifelse(any((tx_id_junction %>% unlist() %>% unname()) %in% mane_transcripts$transcript_id), T, F))
  
  
  ####################
  ## SAVE RESULTS
  ####################
  
  if (identical(df_introns_mane$ref_junID, df_introns$ref_junID)) {
    
    saveRDS(object = df_introns_mane, file = IDB_path)
    
  } else {
    print(paste0(Sys.time(), " - Error! After adding the MANE info some introns have been lost!"))
  }
  
  print(paste0(Sys.time(), " - MANE info saved!"))
  
  rm(mane_transcripts)
  rm(df_introns_mane)
  rm(folder_name)
  rm(df_introns)
  rm(IDB_path)
  gc()
  
}


## MODULO 3  ------------------------------------

get_stats_modulo <- function(cluster,
                             folder_name,
                             PC = T,
                             stats = F) {
  
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
  df_novel %>% head()
  df_novel %>% nrow() 
  
  
  
  
  if (stats) {
    
    ## First hypothesis - Protein-coding novel donor junctions produce less frame-shift
    print(paste0(Sys.time(), " - First hypothesis - Protein-coding novel donor junctions produce less frame-shift."))
    
    donor_data <- df_tidy %>%
      filter(novel_type == "novel_donor", protein_coding == 100) %>%
      mutate(modulo = distance %% 3) 
    
    acceptor_data <- df_tidy %>%
      filter(novel_type == "novel_acceptor", protein_coding == 100) %>%
      mutate(modulo = distance %% 3) 
    
    
    # data <- data.frame(group = "novel_donor",
    #                    value = abs(x_data$modulo))
    # data <- rbind(data,
    #               data.frame(group = "novel_acceptor",
    #                          value = abs(y_data$modulo)))
    
    # group_by(data, group) %>%
    #   summarise(
    #     count = n(),
    #     median = median(value, na.rm = TRUE),
    #     mean = mean(value, na.rm = TRUE),
    #     IQR = IQR(value, na.rm = TRUE)
    #   )
    
    wilcox.test(x = donor_data %>% pull(modulo),
                y = acceptor_data %>% pull(modulo),
                exact = FALSE, 
                alternative = "less") %>% print()
    
    print(paste0((acceptor_data %>% filter(modulo == 0) %>% nrow() * 100) / acceptor_data %>% nrow(), "% of modulo3 = 0 protein-coding acceptor junctions."))
    print(paste0((donor_data %>% filter(modulo == 0) %>% nrow() * 100) / donor_data %>% nrow(), "% of modulo3 = 0 protein-coding donor junctions."))
    
    
    
    ###########################################################################################
    ## Second hypothesis - Never protein-coding novel donor junctions produce less frame-shift
    ###########################################################################################
    
    print(paste0(Sys.time(), " - Second hypothesis - Never protein-coding novel donor junctions produce less frame-shift."))
    
    donor_data <- df_tidy %>%
      filter(novel_type == "novel_donor", protein_coding == 0) %>%
      mutate(modulo = distance %% 3) 
    
    acceptor_data <- df_tidy %>%
      filter(novel_type == "novel_acceptor", protein_coding == 0) %>%
      mutate(modulo = distance %% 3) 
    
    
    
    
    # data <- data.frame(group = "novel_donor",
    #                    value = abs(x_data$modulo))
    # data <- rbind(data,
    #               data.frame(group = "novel_acceptor",
    #                          value = abs(y_data$modulo)))
    
    # group_by(data, group) %>%
    #   summarise(
    #     count = n(),
    #     median = median(value, na.rm = TRUE),
    #     mean = mean(value, na.rm = TRUE),
    #     IQR = IQR(value, na.rm = TRUE)
    #   )
    
    wilcox.test(x = donor_data %>% pull(modulo),
                y = acceptor_data %>% pull(modulo),
                exact = FALSE, 
                alternative = "less") %>% print()
    
    print(paste0((acceptor_data %>% filter(modulo == 0) %>% nrow() * 100) / acceptor_data %>% nrow(), "% of modulo3 = 0 non-protein-coding acceptor junctions."))
    print(paste0((donor_data %>% filter(modulo == 0) %>% nrow() * 100) / donor_data %>% nrow(), "% of modulo3 = 0 non-protein-coding donor junctions."))
  }
}






## CDTS - CONSERVATION ----------------------------

add_cdts_cons_scores <- function(cluster,
                                 folder_name) {
  
  
  
  
  ## Load the IDB ---------------------------------------------------------
  
  print(paste0(Sys.time(), " - loading the intron database..."))
  
  db_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds")) %>%
    distinct(ref_junID, .keep_all = T) %>%
    mutate(CDTS_5ss_mean = 0.0,
           CDTS_3ss_mean = 0.0,
           phastCons20way_5ss_mean = 0.0,
           phastCons20way_3ss_mean = 0.0) %>%
    GRanges()
  
  
  print(paste0(Sys.time(), " - getting scores assigned to introns from the IDB ..."))
  
  ## https://www.nature.com/articles/nature09000
  ## Scores from the 5'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% start() - 10,
                                                                                            end = db_introns %>% start() + 35),
                                                                           strand = db_introns %>% strand()),
                                          ignore.strand = FALSE,
                                          type = "within")
  
  overlaps_tidy <- overlaps %>%
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_5ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_5ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  
  ## Scores from the 3'ss ---------------------------------------------------------
  overlaps <- GenomicRanges::findOverlaps(query = CNC_CDTS_CONS_gr %>% diffloop::rmchr(),
                                          subject = GenomicRanges::GRanges(seqnames = db_introns %>% seqnames(),
                                                                           ranges = IRanges(start = db_introns %>% end() - 35,
                                                                                            end = db_introns %>% end() + 10),
                                                                           strand = db_introns %>% strand()),
                                          
                                          ignore.strand = FALSE,
                                          type = "within")
  
  overlaps_tidy <- overlaps %>%
    as.data.frame() %>%
    mutate(CDTS = CNC_CDTS_CONS_gr[queryHits(overlaps),]$CDTS,
           mean_phastCons20way = CNC_CDTS_CONS_gr[queryHits(overlaps),]$mean_phastCons20way) %>%
    group_by(subjectHits) %>%
    mutate(CDTS_mean = CDTS %>% mean(),
           mean_phastCons20way_mean = mean_phastCons20way %>% mean())
  
  db_introns[subjectHits(overlaps),]$CDTS_3ss_mean <- overlaps_tidy$CDTS_mean
  db_introns[subjectHits(overlaps),]$phastCons20way_3ss_mean <- overlaps_tidy$mean_phastCons20way_mean
  
  #####################
  ## SAVE RESULTS
  #####################
  
  file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds") 
  saveRDS(object = db_introns %>% data.table::as.data.table(), 
          file = file_name)
  
  print(paste0(Sys.time(), " - CDTS and Conservation scores added! IDB updated!"))
  
  
  rm(overlaps_tidy)
  rm(db_introns)
  rm(file_name)
  rm(overlaps)
  
  #gc()
}





## EXON LENGTH ------------------------------------

add_exon_length <- function(cluster,
                            folder_name) {
  
  
  ## Load splice table
  splice_table <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/base_data/homosapiens_v104_splicetable.rds"))
  splice_table %>% head()
  
  
  ## Load intron database
  file_stats <- paste0(folder_name, "/", cluster, "_db_introns.rds")
  df_stats <- readRDS(file = file_stats)
  df_stats <- df_stats %>%
    mutate(left_exon_length = 0,
           right_exon_length = 0)
  df_stats %>% head()
  
  
  
  ## Get the length of both flanking exons
  
  splice_table <- splice_table %>%
    mutate(left_exon_length = splice_table$lend - splice_table$lstart,
           right_exon_length = splice_table$rend - splice_table$rstart)
  
  splice_table %>% head()
  
  
  
  ## Get genomic overlaps between the splice table and the intron database
  overlaps <-  GenomicRanges::findOverlaps(query = GRanges(seqnames = splice_table$seqid,
                                                           IRanges(start = splice_table$lend + 1, 
                                                                   end = splice_table$rstart - 1),
                                                           strand = splice_table$strand),
                                           subject = GRanges(seqnames = df_stats$seqnames,
                                                             IRanges(start = df_stats$start, 
                                                                     end = df_stats$end),
                                                             strand = df_stats$strand),
                                           ignore.strand = F,
                                           type = "equal")
  
  ## QC
  if (subjectHits(overlaps) %>% unique %>% length() !=
      df_stats %>% nrow()) {
    print("Error: some introns from the intron database have not been found in annotation!")
  }
  
  
  
  ## Merge both datasets
  df_stats[subjectHits(overlaps),]$left_exon_length <- splice_table[queryHits(overlaps),]$left_exon_length 
  df_stats[subjectHits(overlaps),]$right_exon_length <- splice_table[queryHits(overlaps),]$right_exon_length 
  
  
  # ## Adding the position of the intron within the transcript -- can't do this, as I don't know which specific transcript was mis-spliced.
  # ## Load the ensembl dataset
  # homo_sapiens_v104_gtf <- rtracklayer::import("/data/references/ensembl/gtf_gff3/v104/Homo_sapiens.GRCh38.104.gtf") %>%
  #   as.data.frame()
  # homo_sapiensv104_exons <- homo_sapiens_v104_gtf %>%
  #   filter(type == "exon")
  # 
  # overlaps <-  GenomicRanges::findOverlaps(query = GRanges(seqnames = homo_sapiensv104_exons$seqnames,
  #                                                          IRanges(start = homo_sapiensv104_exons$end + 1, 
  #                                                                  end = homo_sapiensv104_exons$end + 1),
  #                                                          strand = homo_sapiensv104_exons$strand),
  #                                          subject = GRanges(seqnames = df_stats$seqnames,
  #                                                            IRanges(start = df_stats$start, 
  #                                                                    end = df_stats$start),
  #                                                            strand = df_stats$strand),
  #                                          ignore.strand = F,
  #                                          type = "equal")
  # 
  # ## QC
  # if (subjectHits(overlaps) %>% unique %>% length() !=
  #     df_stats %>% nrow()) {
  #   print("Error: some introns from the intron database have not been found in annotation!")
  # }
  # 
  # # homo_sapiensv104_exons[queryHits(overlaps),][2:3,]
  # # df_stats[subjectHits(overlaps),][2:3,]
  # # homo_sapiensv104_exons[queryHits(overlaps),] %>%
  # #   filter( (transcript_id == "ENST00000379236" ||
  # #             transcript_id == "ENST00000453580")) %>%
  # #   filter(end == "1213093")
  #   
  # 
  # ## Merge both datasets
  # df_stats <- df_stats %>%
  #   mutate(intron_number = 0)
  # df_stats[subjectHits(overlaps),]$intron_number <- (homo_sapiensv104_exons[queryHits(overlaps),]$exon_number %>% as.integer()) + 1 
  # 
  # df_stats %>% head()
  
  
  
  file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
  saveRDS(object = df_stats, file = file_name)
  
  print(paste0(Sys.time(), " - file saved!"))
}





## GENE TPM AND GENE LENGTH -----------------------

add_gene_tpm_length <- function(cluster,
                                tpm_file,
                                folder_name,
                                GTEx = T) {
  
  print(cluster)
  
  ## LOAD SOURCE DATA
  file_name <- paste0(folder_name, "/", cluster, "_db_introns.rds")
  df_introns <- readRDS(file = file_name) 
  
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
    dplyr::rename(tpm_mean_rct3 = tpm_mean)%>% 
    dplyr::rename(tpm_median_rct3 = tpm_median)
  
  
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
    filter(type == "gene") %>%
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
      filter(width < gene_width)
    
    df_merged %>% nrow()
  }
  
  df_merged %>% head()
  df_merged %>% nrow()
  
  ## Remove duplicated ref introns (this is due to the 'unnest' function applied before) ----------------
  
  df_merged <- df_merged %>%
    distinct(ref_junID, .keep_all = T)
  
  
  ## Matching results with novel introns ---------------------------------------
  
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))%>%
    filter(ref_junID %in% df_merged$ref_junID)
  
  if (!identical(df_novel$ref_junID %>% unique() %>% sort(), 
                 df_merged %>% filter(ref_type != "never") %>% dplyr::select(ref_junID) %>% pull()%>% unique() %>% sort())) {
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
                   df_introns %>% filter(ref_type != "never") %>% select(ref_junID) %>% pull()%>% unique() %>% sort())) {
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



add_protein_percentage_to_database <- function(cluster,
                                               folder_root) {
  
  
  ## Import HUMAN REFERENCE transcriptome
  # homo_sapiens_v104_gtf <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf") %>% 
  #   as.data.frame()
  
  print(paste0(Sys.time(), " - adding biotypes percentage to junctions..."))
  
  df_protein <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_annotated_SR_details_length_104_biotype.rds")
  df_protein %>% head()
  df_protein %>% nrow()
  
  
  
  ## Load mis-splicing data for the current cluster
  print(paste0(Sys.time(), " - Adding biotype percentage to '", cluster, "' junctions..."))
  
  if (cluster == "case" || cluster == "control") {
    folder_name <- paste0(folder_root, "/", cluster, "/pipeline3/missplicing-ratio/")
    df_introns <- readRDS(file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
  } else {
    folder_name <- folder_root #paste0(folder_root, "/results/pipeline3/missplicing-ratio/", cluster, "/v105")
    df_introns <- readRDS(file = paste0(folder_root, "/", cluster, "_db_introns.rds")) %>%
      filter(u12_intron == T | u2_intron == T)
  }
  
  
  df_introns %>% head()
  df_introns %>% nrow()
  df_introns %>% distinct(ref_junID) %>% nrow()
  
  df_protein_local <- df_protein %>%
    filter(junID %in% df_introns$ref_junID) 
  
  df_protein_local %>% head()
  df_protein_local %>% nrow()
  df_introns %>% distinct(ref_junID) %>% nrow()
  
  
  ## Load database of novel junctions
  df_novel <- readRDS(file = paste0(folder_name, "/", cluster, "_db_novel.rds"))  %>%
    filter(ref_junID %in% df_introns$ref_junID)
  
  df_novel %>% head()
  df_novel %>% distinct(ref_junID) %>% nrow()
  
  
  ## QC - novel junctions must be attached to the same introns stored on the intron database
  if ((df_introns %>% filter(ref_type != "never") %>% distinct(ref_junID) %>% nrow()) != 
      (df_novel %>% distinct(ref_junID) %>% nrow())) {
    print("Error: some novel junctions are attached to introns that are not stored on the intron database!")
  }
  
  
  if ((df_introns %>% distinct(ref_junID) %>% nrow()) == (df_protein_local %>% distinct(junID) %>% nrow())) {
    
    df_introns %>% nrow()
    
    ## Add protein-coding to the IDB
    
    df_introns <- merge(x = df_introns %>% data.table::as.data.table(),
                        y = df_protein_local %>% data.table::as.data.table(),
                        by.x = "ref_junID",
                        by.y = "junID",
                        all.x = T)
    df_introns %>% head()
    
    
    
    
    ## Add protein-coding info to the novel junction database
    
    df_novel <- merge(x = df_novel %>% data.table::as.data.table(),
                      y = df_protein_local %>% data.table::as.data.table(),
                      by.x = "ref_junID",
                      by.y = "junID",
                      all.x = T)
    df_novel %>% head()
    
    
    
    
    if(identical(df_introns %>% filter(ref_type != "never") %>% distinct(ref_junID) %>% pull %>% sort(), 
                 df_novel %>% distinct(ref_junID) %>% pull %>% sort())) {
      
      saveRDS(object = df_introns %>% data.table::as.data.table(),
              file = paste0(folder_name, "/", cluster, "_db_introns.rds"))
      
      saveRDS(object = df_novel,
              file = paste0(folder_name, "/", cluster, "_db_novel.rds"))
      
      print(paste0(Sys.time(), " - ", cluster, " finished!"))
      
    } else {
      print(paste0(Sys.time(), " - ", cluster, " error - the two datasets present different refIDs of rows!"))
    }
    
    
  } else {
    print(paste0(Sys.time(), " - ", cluster, " error - the two datasets present different number of rows!"))
  }
  
  rm(df_protein_local)
  rm(df_introns)
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
    filter(u12_intron == T | u2_intron == T)
  
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
    filter(junID == df_introns[row_number,]$ref_junID)
  
  homo_sapiens_v104_gtf %>%
    filter(type == "transcript",
           transcript_id %in% tx)
  
  
} 



