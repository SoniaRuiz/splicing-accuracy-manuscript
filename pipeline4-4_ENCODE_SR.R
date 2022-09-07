#################################
## LOAD SOURCE DEPENDENCIES
#################################

source("/home/sruiz/PROJECTS/splicing-project/pipeline1_QC_split_reads.R")
source("/home/sruiz/PROJECTS/splicing-project/pipeline3_methods.R")


# source("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/encode.R")


# project_id <- "U2AF1"
# project_id <- "SF3B1"
# project_id <- "U2AF2"
# project_id <- "SRSF1"



################################
## LOAD REFERENCE TRANSCRIPTOME
################################

if (!exists(x = "edb")) {
  ## Load Reference Transcriptome
  edb <- ensembldb::ensDbFromGtf("/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf", 
                                 outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.105.sqlite"))
  edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.105.sqlite"))
  
} else {
  print("'edb - Homo_sapiens.GRCh38.105.sqlite' DB already loaded!")
}


############################
## FUNCTIONS - ENCODE
############################

# 1. PREPARE DATA ------------------------------------------------------------------

set_ENCODE_samples <- function (project_id) {
  
  
  if (project_id == "PTBP1") {
    
    for (type in c("control", "case")) {
      
      samples <- data.frame(samples = as.character(),
                            RIN = as.double(),
                            cell_line = as.character(),
                            read_depth = as.integer())
      
      if (type == "control") {
        
        ## Two samples from K562 cell lines
        samples <- data.frame(samples = c("ENCFF487ZNA", "ENCFF264WCX"),
                              RIN = c(9.2, 9),
                              cell_line = "K562",
                              read_depth = c(28510158,29426768))
        
        
        
        ## Two samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF613CGT", "ENCFF709LHN"),
                                    RIN = c(9.5),
                                    cell_line = "HepG2",
                                    read_depth = c(308023099,40379990)))
        
      } else if (type == "case") {
        
        ## Two samples from K562 cell lines
        samples <- data.frame(samples = c("ENCFF006ZOW", "ENCFF556DIN"),
                              RIN = c(9.1, 9.2),
                              cell_line = "K562",
                              read_depth = c(27522009,28793900))
        
        
        ## Two samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF707AEZ", "ENCFF125RUG"),
                                    RIN = c(9.5, 9.6),
                                    cell_line = "HepG2",
                                    read_depth = c(21205496,35148355)))
        
      }
      
      saveRDS(object = samples, file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                              project_id, "/", type, "/samples.rds"))
      
    }
    
  } else if (project_id == "U2AF1") {
    
    for (type in c("control", "case")) {
      
      samples <- data.frame(samples = as.character(),
                            RIN = as.double(),
                            cell_line = as.character(),
                            read_depth = as.integer())
      
      if (type == "control") {
        
        ## Two samples from K562 cell lines
        samples <- data.frame(samples = c("ENCFF870ZUC", "ENCFF089NKI"),
                              RIN = c(9.2),
                              cell_line = "K562",
                              read_depth = c(23612508,27216994))
        
        
        
        ## Two samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF032IND", "ENCFF304ITC"),
                                    RIN = c(10),
                                    cell_line = "HepG2",
                                    read_depth = c(19538219,32583010)))
        
      } 
      else if (type == "case") {
        
        ## Two samples from K562 cell lines
        samples <- data.frame(samples = c("ENCFF462PAF", "ENCFF323XGR"),
                              RIN = c(9.3,9),
                              cell_line = "K562",
                              read_depth = c(22358519,22672626))
        
        
        ## Two samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF723JZD", "ENCFF470EQC"),
                                    RIN = c(9.9),
                                    cell_line = "HepG2",
                                    read_depth = c(32507151,53747711)))
        
      }
      
      saveRDS(object = samples, file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                              project_id, "/", type, "/samples.rds"))
      
    }
    
  } else if (project_id == "SF3B1") {
    
    for (type in c("control", "case")) {
      
      samples <- data.frame(samples = as.character(),
                            RIN = as.double(),
                            cell_line = as.character(),
                            read_depth = as.integer())
      
      if (type == "control") {
        
        ## Two samples from K562 cell lines (adult, 53)
        samples <- data.frame(samples = c("ENCFF172ZXY", "ENCFF525ALF"),
                              RIN = c(10,9.4),
                              cell_line = "K562",
                              read_depth = c(27142780,27735686))
        
        ##And two more samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF965UPD", "ENCFF958MMI"),
                                    RIN = c(10),
                                    cell_line = "HepG2",
                                    read_depth = c(33700253,36662855)))
        
        
      } else if (type == "case") {
        
        ## Two samples from K562 cell lines (adult, 53)
        samples <- data.frame(samples = c("ENCFF379JAB", "ENCFF946HGK"),
                              RIN = c(9.7),
                              cell_line = "K562",
                              read_depth = c(33623700,31923353))
        
        ## And two more samples from HepG2 cell lines (child, 15)
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF889YFX", "ENCFF400JIZ"),
                                    RIN = c(9.8,10),
                                    cell_line = "HepG2",
                                    read_depth = c(29422854,33173166)))
        
      }
      
      saveRDS(object = samples, file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                              project_id, "/", type, "/samples.rds"))
      
    }
    
  } else if (project_id == "U2AF2") {
    
    for (type in c("control", "case")) {
      
      samples <- data.frame(samples = as.character(),
                            RIN = as.double(),
                            cell_line = as.character(),
                            read_depth = as.integer())
      
      if (type == "control") {
        
        ## Two samples from K562 cell lines (adult, 53)
        
        samples <- data.frame(samples = c("ENCFF089NKI", "ENCFF870ZUC"),
                              RIN = 9.2,
                              cell_line = "K562",
                              read_depth = c(27216994,23612508))
        
        
        ##And two more samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF905ARG", "ENCFF482MPL"),
                                    RIN = c(9.8,9.7),
                                    cell_line = "HepG2",
                                    read_depth = c(24095732,28490748)))
        
      } else if (type == "case") {
        
        ## Two samples from K562 cell lines (adult, 53)
        samples <- data.frame(samples = c("ENCFF397YUK", "ENCFF748VWE"),
                              RIN = c(9.2,9),
                              cell_line = "K562",
                              read_depth = c(23633966,22337802))
        
        
        ## And two more samples from HepG2 cell lines (child, 15)
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF141OTC", "ENCFF341UWM"),
                                    RIN = c(9.8),
                                    cell_line = "HepG2",
                                    read_depth = c(41105739,52030621)))
      }
      
      saveRDS(object = samples, file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                              project_id, "/", type, "/samples.rds"))
      
    }
    
  } else if (project_id == "SRSF1") {
    
    for (type in c("control", "case")) {
      
      samples <- data.frame(samples = as.character(),
                            RIN = as.double(),
                            cell_line = as.character(),
                            read_depth = as.integer())
      
      if (type == "case") {
        ## Two samples from K562 cell lines (adult, 53)
        samples <- data.frame(samples = c("ENCFF767UXH", "ENCFF590VOY"),
                              RIN = 9.2,
                              cell_line = "K562",
                              read_depth = c(23695454,27397149))
        
        ##And two more samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF694NKE", "ENCFF098AXX"),
                                    RIN = c(9.4,9.5),
                                    cell_line = "HepG2",
                                    read_depth = c(35705730,30439904)))
        
      } else if (type == "control") {
        
        ## Two samples from K562 cell lines (adult, 53)
        samples <- data.frame(samples = c("ENCFF487ZNA", "ENCFF264WCX"),
                              RIN = c(9.2,9),
                              cell_line = "K562",
                              read_depth = c(28510158,29426768))
        
        ##And two more samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF709LHN", "ENCFF613CGT"),
                                    RIN = 9.5,
                                    cell_line = "HepG2",
                                    read_depth = c(35705730,30802309)))
      }
      saveRDS(object = samples, file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                              project_id, "/", type,  "/samples.rds"))
    }
  } else if (project_id == "UPF1") {
    
    for (type in c("control", "case")) {
      
      samples <- data.frame(samples = as.character(),
                            RIN = as.double(),
                            cell_line = as.character(),
                            read_depth = as.integer())
      
      if (type == "case") {
        ## Two samples from K562 cell lines (adult, 53)
        samples <- data.frame(samples = c("ENCFF682DLV", "ENCFF385PCE"),
                              RIN = 10,
                              cell_line = "K562",
                              read_depth = c(35405735,44671298))
        
        ##And two more samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF209ZXB", "ENCFF098NOF"),
                                    RIN = c(9.9,9.8),
                                    cell_line = "HepG2",
                                    read_depth = c(30732488,34992506)))
        
      } else if (type == "control") {
        
        ## Two samples from K562 cell lines (adult, 53)
        samples <- data.frame(samples = c("ENCFF877NSE", "ENCFF064JGI"),
                              RIN = c(10),
                              cell_line = "K562",
                              read_depth = c(35712770,30546711))
        
        ##And two more samples from HepG2 cell lines
        samples <- rbind(samples,
                         data.frame(samples = c("ENCFF678DFI", "ENCFF323SGY"),
                                    RIN = 10,
                                    cell_line = "HepG2",
                                    read_depth = c(37076328,37613468)))
      }
      saveRDS(object = samples, file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                              project_id, "/", type,  "/samples.rds"))
    }
  }
}

# extract_ENCODE_junctions()
extract_all_ENCODE_junctions <- function (project_ids, type_ENCODE) {
  
  # project_ids <- c("U2AF1", "SF3B1", "U2AF2", "SRSF1")
  
  all_junc <- map_df(project_ids, function(project_id) {
    
    print(paste0(Sys.time(), " --> starting '", project_id, "' ...!"))
    
    map_df(c("case", "control"), function(type) {
      
      # type <- "control"
      
      print(paste0(Sys.time(), " --> '", type, "' samples:"))
      
      ## Load samples
      samples <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                       project_id, "/", type,  "/samples.rds"))
      
      
      ###########################
      ## EXTRACT JUNCTION DATA
      ###########################
      
      df_junc <- map_df(samples$samples, function(sample) {
        # sample <- samples[1]
        junc <- read.table(paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", project_id, "/", 
                                  type,  "/", sample, ".bam.sort.s0.junc"), #ENCFF870ZUC.bam.junc",
                           header = FALSE, 
                           sep = "\t", 
                           stringsAsFactors = FALSE, quote = "") %>%
          as.data.frame() %>%
          dplyr::select(chr = V1,
                        start = V2,
                        stop = V3,
                        junID = V4,
                        reads = V5,
                        strand = V6,
                        blockSizes = V11) %>%
          dplyr::mutate(strand = ifelse(strand == "?", "*", strand)) %>%
          separate(col = blockSizes, sep =",", c("blockSizesStart", "blockSizesEnd")) %>%
          mutate(sampleID = sample) %>%
          GRanges() %>%
          diffloop::rmchr() 
        
        ranges(junc) <- IRanges(start = start(junc) + (junc$blockSizesStart %>% as.integer()), 
                                end = end(junc) - (junc$blockSizesEnd %>% as.integer()))
        
        if (any(junc %>% width() <= 25)) {
          print("Error: junctions should be longer than 25bp!")
        }
        
        junc_tidy <- keepSeqlevels(x = junc,
                                   value = c( "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                              "2", "20", "21", "22", "3", "4", "5", "6", "7", "8", "9", "X", "Y"), 
                                   pruning.mode = "tidy")
        
        
        if (any(junc_tidy %>% width() <= 25)) {
          print("Error: junctions should be longer than 25bp!")
        }
        
        ## Junctions are in a .bed file format, which has a coordinate system zero-based for the coordinate start and one-based for the coordinate end.
        ## We'll annotate the junctions using a .gtf file, which has a coordinate system one-based at both ends.
        ## Therefore, we need to increase +1 the starting coordinate of the bed, to adjust it from a zero-based to a one-based coordinate.
        ranges(junc_tidy) <- IRanges(start = start(junc_tidy) + 1, 
                                     end = end(junc_tidy))
        
        print(paste0(Sys.time(), " --> sample '", sample, "' finished!"))
        
        return(junc_tidy %>% as.data.frame())
        
      })
      
      return(df_junc)
      
    })
  })
  
  ## Save raw split reads file
  saveRDS(object = all_junc,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/common_data/",
                        type_ENCODE, "/all_raw_split_reads.rds"))
  
  #############################
  ## ASSIGN UNIQUE IDs
  #############################
  
  ## Save raw split reads file
  all_junc <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/common_data/",
                                    type_ENCODE, "/all_raw_split_reads.rds"))
  
  all_junc %>% nrow()
  all_junc %>% head()
  all_junc$sampleID %>% unique
  
  df_junc_tidy <- all_junc %>% #filter(start == "16789") %>%
    select(-junID, -blockSizesStart, -blockSizesEnd) %>%
    group_by(seqnames, start, end) %>%
    mutate(group_id = cur_group_id())
  
  df_junc_tidy %>% nrow()
  df_junc_tidy %>% head()
  
  df_junc_tidy <- df_junc_tidy %>%
    rename(junID = group_id) %>%
    mutate(junID = sprintf("JUNC%08d", junID))
  
  df_junc_tidy %>% head()
  df_junc_tidy %>% filter(start == "14830", end == "185490") %>% as.data.frame()
  
  ## Save split reads file
  saveRDS(object = df_junc_tidy %>%
            dplyr::select(seqnames, start, end, width, strand, junID) %>%
            distinct(junID, .keep_all = T),
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/common_data/",
                        type_ENCODE, "/all_split_reads.rds"))
  
  print(paste0(Sys.time(), " --> 'all_split_reads' file saved!"))
  
  #############################
  ## 'SPLIT_READ_COUNTS' FILE
  #############################
  
  df_junc_tidy %>% nrow()
  df_junc_tidy_d <- distinct(df_junc_tidy)
  df_junc_tidy_d %>% nrow()
  
  split_read_counts <- df_junc_tidy_d %>%
    as.data.frame() %>%
    dplyr::select(reads, sampleID, junID) %>%
    spread(sampleID, reads)
  
  
  ## QC - 1
  all_junc %>%
    filter(start == "498457", end == "498683") %>% as.data.frame() %>% print()
  
  df_junc_tidy %>%
    filter(start == "498457", end == "498683") %>% as.data.frame() %>% print()
  
  split_read_counts %>%
    filter(junID == df_junc_tidy %>%
             filter(start == "498457", end == "498683") %>%
             distinct(junID) %>% pull(junID)) %>% print()
  
  
  ## QC - 2
  all_junc %>%
    filter(seqnames == 1, start == 17369, end == 17605) %>% as.data.frame() %>% print()
  df_junc_tidy %>%
    filter(seqnames == 1, start == 17369, end == 17605) %>% as.data.frame() %>% print()
  
  split_read_counts %>%
    filter(junID == df_junc_tidy %>%
             filter(seqnames == 1, start == 17369, end == 17605) %>%
             distinct(junID) %>% pull(junID)) %>% print()
  
  
  saveRDS(object = split_read_counts,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/common_data/",
                        type_ENCODE, "/all_split_read_counts.rds"))
  
  print(paste0(Sys.time(), " --> 'all_split_read_counts' file saved!"))
  
}


# annotate_junctions_ENCODE(project_id = "U2AF1", type = "control")
# annotate_junctions_ENCODE(project_id = "U2AF1", type = "case")
annotate_junctions_ENCODE <- function (project_id, type, type_ENCODE) {
  
  ##################################
  ## LOCAL 'SPLIT READ COUNTS' FILE
  ##################################
  
  ## Load samples
  samples <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                   project_id, "/", type,  "/samples.rds"))$samples
  
  all_split_read_counts <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/common_data/",
                                                 type_ENCODE, "/all_split_read_counts.rds")) %>%
    select(junID, all_of(samples %>% as.character()))
  
  all_split_read_counts %>% head()
  
  split_read_counts <- all_split_read_counts %>% 
    filter(if_any(.cols = samples %>% as.character() , ~ !is.na(.)))
  
  split_read_counts %>% head()
  
  split_read_counts %>%
    filter(junID == "JUNC00000049")
  
  if (split_read_counts$junID %>% unique() %>% length() != 
      split_read_counts$junID %>% length()) {
    print("Error!! Some junction IDs are repeated!")
  }
  
  ## Save local 'split_read_counts'
  saveRDS(split_read_counts,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                        project_id, "/", type, "/split_read_counts.rds"))
  
  
  ################################
  ## 'ANNOTATE SPLIT READS' FILE
  ################################
  
  all_split_reads <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/common_data/",
                                           type_ENCODE, "/all_split_reads.rds"))
  
  ind <- which(all_split_reads$junID %in% split_read_counts$junID)
  
  all_split_reads_local <- all_split_reads[ind,]
  
  if (!(identical(all_split_reads_local$junID %>% sort(), 
                  split_read_counts$junID %>% sort()))) {
    print("Error!")
  }
  
  ## ANNOTATE DASPER - Ensembl v104
  split_reads_details <- annotate_dasper(recount2_data = all_split_reads_local %>% GRanges(), edb = edb)
  
  split_reads_details %>% as.data.frame() %>% filter(type == "annotated") %>% nrow()
  split_reads_details %>% as.data.frame() %>% filter(type == "novel_acceptor") %>% nrow()
  split_reads_details %>% as.data.frame() %>% filter(type == "novel_donor") %>% nrow()
  
  
  all_split_reads_details <- generate_max_ent_score(junc_tidy = split_reads_details)
  
  
  saveRDS(object = all_split_reads_details,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", project_id, "/", 
                        type, "/", type, "_annotated_SR_details_length_v104.rds"))
  
  
  print(paste0(Sys.time(), " --> 'annotated_SR_details_length_v104' file saved!"))
}


# 2. CREATE DATABASE ---------------------------------------------------------------


# 3. PLOT STUFF FOR THE PAPER ------------------------------------------------------

# QC_unique_junctions(project_id = "U2AF1")
QC_unique_junctions <- function (project_id = "U2AF1"){
  
  ################################
  ## CHECK THE NOVEL JUNCTION IDs
  ################################
  
  df_all_novel <- map_df(c("case", "control"), function(type) {
    
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                          project_id, "/", type, "/pipeline3/missplicing-ratio/")
    
    df_novel <- readRDS(file = paste0(folder_name, "/", type, "_db_novel.rds")) %>%
      mutate(sample_type = type) 
    
    return(df_novel)
  })
  
  df_case <- df_all_novel %>%
    filter(sample_type == "case") %>%
    GRanges()
  
  df_control <- df_all_novel %>%
    filter(sample_type == "control")%>%
    GRanges()
  
  novel_overlaps <- GenomicRanges::findOverlaps(query = df_case,
                                                subject = df_control,
                                                type = "equal",
                                                ignore.strand = FALSE)
  
  df_common_case <- df_case[queryHits(novel_overlaps),]
  df_common_control <- df_control[subjectHits(novel_overlaps),]
  
  if (!identical(df_common_case$novel_junID %>% unique(),
                 df_common_control$novel_junID %>% unique())) {
    print("Error: Some novel junctions have been assigned with different junctions ID!")
  }
  
  if (intersect(df_case[-queryHits(novel_overlaps),]$novel_junID, 
                df_control[-subjectHits(novel_overlaps),]$novel_junID) %>% length() > 0) {
    print("Error: Some independent introns have been assigned with the same ID!")
  }
  
  
  #############################
  ## CHECK THE REF INTRON IDs
  #############################
  
  df_all_introns <- map_df(c("case", "control"), function(type) {
    
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                          project_id, "/", type, "/pipeline3/missplicing-ratio/")
    
    df_introns <- readRDS(file = paste0(folder_name, "/", type, "_db_introns.rds")) %>%
      mutate(sample_type = type) 
    
    return(df_introns)
  })
  
  df_intron_case <- df_all_introns %>%
    filter(sample_type == "case") %>%
    GRanges()
  
  df_intron_control <- df_all_introns %>%
    filter(sample_type == "control")%>%
    GRanges()
  
  overlaps <- GenomicRanges::findOverlaps(query = df_intron_case,
                                          subject = df_intron_control,
                                          type = "equal",
                                          ignore.strand = FALSE)
  
  df_common_case <- df_intron_case[queryHits(overlaps),]
  df_common_control <- df_intron_control[subjectHits(overlaps),]
  
  if (!identical(df_common_case$ref_junID %>% unique(),
                 df_common_control$ref_junID %>% unique())) {
    print("Error: Some introns have been assigned with different junctions ID!")
  }
  
  if (intersect(df_intron_case[-queryHits(overlaps),]$ref_junID, 
                df_intron_control[-subjectHits(overlaps),]$ref_junID) %>% length() > 0) {
    print("Error: Some independent introns have been assigned with the same ID!")
  }
  
  
  ################################
  ## CHECK THE DISTANCES
  ################################
  
  df_case$distance %>% abs %>% summary()
  df_control$distance %>% abs %>% summary()
  
  df_case %>%
    as.data.frame() %>%
    filter(distance == 0)
  
  df_control %>%
    as.data.frame() %>%
    filter(distance == 0)
  
  df_all_introns %>%
    as.data.frame() %>%
    filter(ref_junID == df_control %>%
             as.data.frame() %>%
             filter(distance == 0) %>%
             pull(ref_junID))
}

plot_distances_ENCODE <- function() {
  
  
  limit_bp <- 30
  con <- dbConnect(RSQLite::SQLite(), "./../splicing-project/splicing-encode-projects/siRNA-knockdown/database/encode.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  for (project_id in df_metadata$SRA_project %>% unique()) {
    ###############################
    ## GET DATA FOR FRONTAL CORTEX
    ###############################
  
    
    df_all <- map_df(c("case", "control"), function(type) {# type <- "case"
      
      query <- paste0("SELECT novel_junID FROM '", type, "_", project_id, "_misspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT * FROM 'novel' WHERE novel_junID IN (",
                      paste(introns$novel_junID, collapse = ","),")")
      df_novel <- merge(x = introns,
                       y = dbGetQuery(con, query) %>% as_tibble(),
                       by = "novel_junID",
                       all.x = T) %>%
        mutate(sample_type = type)
      
      return(df_novel)
    })
    
    df_all <- df_all %>%
      mutate(novel_type = str_replace(string = novel_type,
                                      pattern = "_",
                                      replacement = " "))
    
    
    df_all <- df_all %>%
      mutate(novel_type = factor(novel_type, 
                                 levels = c("novel donor", "novel acceptor")))
    
    
    plot_all <- ggplot(data = df_all) + 
      geom_histogram(aes(x = distance, fill = sample_type),
                     bins = limit_bp * 2,
                     binwidth = 1,
                     position = "stack"
      ) +
      facet_grid(vars(novel_type)) +
      xlab("Distance to the reference intron (in bp)") +
      ylab("Number of unique novel junctions") +
      theme_light() +
      scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                         breaks = c((limit_bp * -1), (round(limit_bp / 2) * -1), 0, round(limit_bp / 2), limit_bp)) +
      scale_fill_manual(values = c("#440154FF","#35B779FF"),
                        breaks = c("case", "control"),
                        labels = c("case", "control")) +
      guides(fill = guide_legend(title = NULL, 
                                 override.aes = list(size = 3),
                                 ncol = 3 )) +
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
      geom_rect(aes(xmin = 0, xmax = limit_bp, ymin = 1, ymax = 100),
                fill = "grey", color = "black") +
      geom_text(aes(x = 15, y = 55),  size = 6, label = "exon") +
      geom_rect(aes(xmin = (limit_bp)*-1, xmax = 0, ymin = 49, ymax = 51),
                fill = "grey", alpha = 1, color = "black") +
      geom_text(aes(x = -15,y = 70),  size = 6, label = "intron") +
      theme_void()
    
    
    plot_all / distance_rectangle +  patchwork::plot_layout(heights = c(8, 1))
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/siRNA-knockdown/figures/",
                        project_id, "_distances.svg")
    ggplot2::ggsave(filename = file_name,
                    width = 183, height = 183, units = "mm", dpi = 300)
  }
}

plot_maxentscan_ENCODE <- function() {
  

  con <- dbConnect(RSQLite::SQLite(), "./../splicing-project/splicing-encode-projects/siRNA-knockdown/database/encode.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  for (project_id in df_metadata$SRA_project %>% unique()) {
    
    
    df_all <- map_df(c("case", "control"), function(type) {# type <- "case"
      
      query <- paste0("SELECT ref_junID, ref_mean_counts, ref_type FROM '", type, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT DISTINCT ref_junID, ref_mean_counts, ref_type FROM '", type, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      query <- paste0("SELECT ref_junID, ref_ss5score, ref_ss3score FROM 'intron' WHERE ref_junID IN (",
                      paste(introns$ref_junID, collapse = ","),")")
      introns <- merge(x = introns,
                       y = dbGetQuery(con, query) %>% as_tibble(),
                       by = "ref_junID",
                       all.x = T) %>% as_tibble() 
      
      
      # introns <- readRDS(file = paste0(folder_root, "/results/pipeline3/missplicing-ratio/",
      #                                  cluster, "/v105/", cluster, "_db_introns.rds"))
      
      ###########################
      ## GET THE NOVEL JUNCTIONS
      ###########################
      
      query <- paste0("SELECT ref_junID, novel_junID, novel_mean_counts FROM '", type, "_", project_id, "_misspliced'")
      novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
      query <- paste0("SELECT novel_junID, novel_type, novel_ss5score, novel_ss3score FROM 'novel' WHERE novel_junID IN (", 
                      paste(novel_junctions$novel_junID, collapse = ","), ")")
      novel_junctions <- merge(x = novel_junctions,
                               y = dbGetQuery(con, query) %>% as_tibble(),
                               by = "novel_junID",
                               all.x = T) %>% as_tibble() 
      
      df_merged <- merge(x = introns %>% dplyr::select(ref_junID, ref_ss5score,
                                                       ref_ss3score, ref_type),
                         y = novel_junctions %>% dplyr::select(ref_junID, novel_junID, 
                                                               novel_ss5score, novel_ss3score, novel_type),
                         by = "ref_junID")  %>%
        mutate(diff_ss5score = ref_ss5score - novel_ss5score,
               diff_ss3score = ref_ss3score - novel_ss3score,
               type = type)
      

      
      return(df_merged)
    })

    
    
    ############################################
    ## PLOT
    ############################################
    
    
    ## Gather
    df_all_tidy <- df_all %>% 
      dplyr::select(diff_ss5score, diff_ss3score, type) %>%
      tidyr::gather(key = "delta_type", value = "mean", -type) %>%
      dplyr::filter(mean != 0)
    
    ## Relevel factors
    df_all_tidy <- df_all_tidy %>%
      mutate(delta_type = ifelse(delta_type == "diff_ss5score", "delta MES 5'ss", "delta MES 3'ss")) %>%
      mutate(delta_type = delta_type %>% as.factor()) %>%
      mutate(delta_type = relevel(delta_type, ref = "delta MES 5'ss")) 
    
    # library(scales)
    # scales::show_col(viridis_pal(option = "magma")(20))
    
    
    ## Plot
    ggplot(data = df_all_tidy) +
      geom_density(mapping = aes(x = mean, fill = type), 
                   alpha = 0.6) +
      geom_vline(xintercept = 0) +
      facet_grid(vars(delta_type)) +
      ggtitle(paste0("Delta MaxEntScan scores - ", project_id, " knockdown.")) +
      xlab("Delta MaxEntScan score") +
      theme_light() +
      scale_fill_manual(values = c("#440154FF","#35B779FF"),
                        breaks = c("case", "control"),
                        labels = c("case", "control")) +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "14"),
            axis.title = element_text(colour = "black", size = "14"),
            strip.text = element_text(colour = "black", size = "16"), 
            legend.text = element_text(colour = "black", size = "14"),
            plot.caption = element_text(colour = "black", size = "14"),
            plot.title = element_text(colour = "black", size = "16"),
            legend.title = element_text(colour = "black", size = "14"),
            legend.position = "top") +
      guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1))
    
    
    file_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/siRNA-knockdown/figures/",
                        project_id, "_deltaMES.svg")
    ggplot2::ggsave(filename = file_name,
                    width = 183, height = 183, units = "mm", dpi = 300)
  
  }
  
}


plot_MSR_ENCODE <- function() {
  
  
  con <- dbConnect(RSQLite::SQLite(), "./../splicing-project/splicing-encode-projects/siRNA-knockdown/database/encode.sqlite")
  dbListTables(con)
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query)
  
  
  for (project_id in df_metadata$SRA_project %>% unique()) {
    
    
    df_all <- map_df(c("case", "control"), function(type) {
      
      query <- paste0("SELECT ref_junID, MSR_D,  MSR_A, ref_type FROM '", type, "_", project_id, "_nevermisspliced'")
      introns <- dbGetQuery(con, query) %>% as_tibble()
      query <- paste0("SELECT ref_junID, MSR_D,  MSR_A, ref_type FROM '", type, "_", project_id, "_misspliced'")
      introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
      
      return(introns%>% 
               mutate(sample_type = type))
    })
    

    
    df_all <- df_all %>%
      mutate(sample_type = factor(sample_type, 
                                 levels = c("case", "control")))
    
    df_all <- df_all %>%
      dplyr::select(sample_type, MSR_D, MSR_A) %>%
      gather(key = "MSR_type", value = "MSR", -sample_type)
    
    plotMSR <- ggplot(data = df_all) + 
      geom_density(aes(x = MSR, fill = sample_type), alpha = 0.8) +
      facet_grid(vars(MSR_type)) +
      xlab("Mis-splicing ratio") +
      theme_light() +
      scale_fill_manual(values = c("#35B779FF","#440154FF"),
                        breaks = c("MSR_D","MSR_A"),
                        labels = c("MSR_Donor","MSR_Acceptor")) +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12"),
            strip.text =  element_text(colour = "black", size = "12"),
            plot.title = element_text(colour = "black", size = "12"),
            legend.text = element_text(size = "12"),
            legend.title = element_text(size = "12"),
            legend.position = "top") +
      guides(fill = guide_legend(title = NULL,
                                 ncol = 2, 
                                 nrow = 1))
    plotMSR
    
    wilcox.test(x = df_all %>% filter(sample_type == "case", MSR_type == "MSR_D") %>% pull(MSR),
                y = df_all %>% filter(sample_type == "case", MSR_type == "MSR_A") %>% pull(MSR),
                alternative = "less")
    wilcox.test(x = df_all %>% filter(sample_type == "control", MSR_type == "MSR_D") %>% pull(MSR),
                y = df_all %>% filter(sample_type == "control", MSR_type == "MSR_A") %>% pull(MSR),
                alternative = "less")
    
    wilcox.test(x = df_all %>% filter(sample_type == "control", MSR_type == "MSR_D") %>% pull(MSR),
                y = df_all %>% filter(sample_type == "case", MSR_type == "MSR_D") %>% pull(MSR),
                alternative = "less")
    wilcox.test(x = df_all %>% filter(sample_type == "control", MSR_type == "MSR_A") %>% pull(MSR),
                y = df_all %>% filter(sample_type == "case", MSR_type == "MSR_A") %>% pull(MSR),
                alternative = "less")
    
    df_all %>% filter(sample_type == "control", MSR_type == "MSR_D") %>% pull(MSR) %>% summary()
    df_all %>% filter(sample_type == "case", MSR_type == "MSR_D") %>% pull(MSR) %>% summary()
    df_all %>% filter(sample_type == "control", MSR_type == "MSR_A") %>% pull(MSR) %>% summary()
    df_all %>% filter(sample_type == "case", MSR_type == "MSR_A") %>% pull(MSR) %>% summary()
  }
}
      
get_lm_ENCODE <- function(type = "case") {
  
  
  ## Load the CASE IDB
  type = "case"
  folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                        project_id, "/", type, "/pipeline3/missplicing-ratio/", type, "_db_introns.rds")
  db_introns <- readRDS(file = folder_name)
  idb_case <- db_introns %>%
    as.data.frame() %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    filter(u2_intron == T | u12_intron == T) %>%
    dplyr::rename(intron_length = width,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score)
  
  
  
  
  ## Load the CONTROL IDB
  type = "control"
  folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                        project_id, "/", type, "/pipeline3/missplicing-ratio/", type, "_db_introns.rds")
  db_introns <- readRDS(file = folder_name)
  idb_control <- db_introns %>%
    as.data.frame() %>%
    dplyr::distinct(ref_junID, .keep_all = T) %>%
    filter(u2_intron == T | u12_intron == T) %>%
    dplyr::rename(intron_length = width,
                  intron_5ss_score = ref_ss5score,
                  intron_3ss_score = ref_ss3score)
  
  idb_case %>% nrow()
  idb_control %>% nrow()
  
  # idb_case %>% filter(ref_type == "donor") %>% nrow()
  # idb_control %>% filter(ref_type == "donor") %>% nrow()
  # 
  # idb_case %>% filter(ref_type == "acceptor") %>% nrow()
  # idb_control %>% filter(ref_type == "acceptor") %>% nrow()
  # 
  # idb_case %>% filter(ref_type == "never") %>% nrow()
  # idb_control %>% filter(ref_type == "never") %>% nrow()
  # 
  # idb_case$ref_missplicing_ratio_tissue_NA %>% summary()
  # idb_control$ref_missplicing_ratio_tissue_NA %>% summary()
  
  
  
  common_junc <- intersect(idb_control$ref_junID %>% unique,
                           idb_case$ref_junID %>% unique)
  
  idb_case <- idb_case %>%
    filter(ref_junID %in% common_junc)
  idb_control <- idb_control %>%
    filter(ref_junID %in% common_junc)
  
  idb_case %>% distinct(ref_junID, .keep_all = T) %>% nrow()
  idb_control %>% distinct(ref_junID, .keep_all = T) %>% nrow()
  
  
  idb_case %>% filter(ref_type == "donor") %>% nrow()
  idb_control %>% filter(ref_type == "donor") %>% nrow()
  
  idb_case %>% filter(ref_type == "acceptor") %>% nrow()
  idb_control %>% filter(ref_type == "acceptor") %>% nrow()
  
  idb_case %>% filter(ref_type == "never") %>% nrow()
  idb_control %>% filter(ref_type == "never") %>% nrow()
  
  idb_case %>% filter(ref_type == "both") %>% nrow()
  idb_control %>% filter(ref_type == "both") %>% nrow()
  
  
  idb_case$ref_missplicing_ratio_tissue_ND %>% summary()
  idb_control$ref_missplicing_ratio_tissue_ND %>% summary()
  idb_case$ref_missplicing_ratio_tissue_NA %>% summary()
  idb_control$ref_missplicing_ratio_tissue_NA %>% summary()
  
  idb <- rbind(idb_case %>% mutate(sample_type = "case"),
               idb_control %>% mutate(sample_type = "control"))
  
  ggplot(data = idb) + 
    geom_density(mapping = aes(ref_missplicing_ratio_tissue_ND, colour = sample_type)) +
    facet_zoom(xlim = c(0,0.1)) +
    #facet_grid(vars(sample_type)) +
    ggtitle("Variance of the estimate across 40 tissues") +
    ylab("Variance of the estimate") +
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
  
  lm(ref_missplicing_ratio_tissue_ND ~ 
       intron_length  +
       intron_5ss_score  +
       intron_3ss_score,
     data = idb_case) %>% summary()
  lm(ref_missplicing_ratio_tissue_ND ~ 
       intron_length  +
       intron_5ss_score  +
       intron_3ss_score,
     data = idb_control) %>% summary()
  
  lm(ref_missplicing_ratio_tissue_NA ~ 
       intron_length  +
       intron_5ss_score  +
       intron_3ss_score,
     data = idb_case) %>% summary()
  lm(ref_missplicing_ratio_tissue_NA ~ 
       intron_length  +
       intron_5ss_score  +
       intron_3ss_score,
     data = idb_control) %>% summary()
  
  # get_lm(idb = idb,
  #        folder_name = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue),
  #        GTEx = T)
}

# ----------------------------------------------------------------------------------



############################
## CALLS
############################

# project_ids <- c("U2AF1", "SF3B1", "U2AF2", "SRSF1")

project_ids <- c("PTBP1")
project_type <- "RBP"

# set_ENCODE_samples(project_id = "PTBP1")
# project_id <- project_ids[1]

for (project_id in project_ids) {
  
  extract_all_ENCODE_junctions(project_ids, type_ENCODE = project_type)
  
  for (sample_type in c("case", "control")) {
    
    print(paste0(Sys.time(), " - processing '", sample_type, "' samples ..."))
    
    # sample_type <- "case"
    
    annotate_junctions_ENCODE(project_id, type = sample_type, type_ENCODE = project_type)
    
    
    ############################################
    ## LOAD DATA FOR THE CURRENT PROJECT
    ############################################
    
    ## Load samples
    samples <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                     project_id, "/", sample_type,  "/samples.rds"))$samples
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                          project_id, "/", sample_type,  "/pipeline3/distances/")
    dir.create(file.path(folder_name), recursive = TRUE, showWarnings = F)
    
    
    ## Load split read data
    all_split_reads_details_104 <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                                         project_id, "/", sample_type, "/", sample_type, 
                                                         "_annotated_SR_details_length_v104.rds"))
    ## Load split read counts
    split_read_counts <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                                               project_id, "/", sample_type, "/split_read_counts.rds"))
    
    
    ############################################
    ## DISTANCES SUITE OF FUNCTIONS
    ############################################
    
    get_distances(cluster = sample_type,
                  samples,
                  split_read_counts,
                  all_split_reads_details_104,
                  folder_name)
    
    extract_distances(cluster = sample_type,
                      samples,
                      split_read_counts = split_read_counts,
                      folder_name = folder_name)
    
    
    get_never_misspliced(cluster = sample_type,
                         samples = samples,
                         split_read_counts = split_read_counts,
                         all_split_reads_details_104 = all_split_reads_details_104,
                         folder_name = folder_name,
                         save_results = T)
    
    
    extract_never_misspliced(cluster = sample_type,
                             samples = samples,
                             split_read_counts = split_read_counts,
                             folder_name = folder_name)
    
    add_never_misspliced_to_df(cluster = sample_type,
                               samples = samples,
                               all_split_reads_details_104 = all_split_reads_details_104,
                               folder_name = folder_name)
    
    ############################################
    ## MIS-SPLICING RATIO SUITE OF FUNCTIONS
    ############################################
    
    folder_save_name <-paste0("/home/sruiz/PROJECTS/splicing-project/splicing-encode-projects/", 
                              project_id, "/", sample_type,  "/pipeline3/missplicing-ratio/")
    dir.create(file.path(folder_save_name), recursive = T, showWarnings = F)
    
    get_missplicing_ratio(cluster = sample_type,
                          samples = samples,
                          split_read_counts = split_read_counts,
                          folder_name = folder_name,
                          folder_save_name = folder_save_name)
    
    add_missplicing_class_to_df(cluster = sample_type,
                                folder_name = folder_save_name)
    
    get_missplicing_QC(cluster = sample_type,
                       samples = samples,
                       split_read_counts = split_read_counts,
                       all_split_reads_details_104 = all_split_reads_details_104,
                       folder_name = folder_save_name)
    
    
    ############################################
    ## ADD FEATURES TO THE IDB
    ############################################
    
    remove_MT_genes(cluster = sample_type,
                    folder_name = folder_save_name)
    add_intron_type(cluster = sample_type,
                    folder_name = folder_save_name)
    
  }
}
