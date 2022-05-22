######################################################################
## PIPELINE 1 - QC functions to clean up a set of raw split reads 
## from recount2 project.
## Sonia García-Ruiz - s.ruiz@ucl.ac.uk - 2020
######################################################################
library(GenomicRanges)
# library(lintr)
# lintr::lint("/home/sruiz/PROJECTS/splicing-project/pipeline1_QC_split_reads.R")

#########################
## 0 - Load libraries
#########################

# library(stringr)
# library(data.table)
# library(RSQLite)
# #library(refGenome)
# library(diffloop)
# library(rtracklayer)
# library(tidyverse)
# #library(Xmisc)
# library(dasper)
# library(ggforce)
# library(ggpubr)

#detach("package:Xmisc", unload = TRUE)

## source("/home/sruiz/PROJECTS/splicing-project/pipeline1_QC_split_reads.R")




###############################
#### LOAD SOURCE FILES    #####
###############################
#GTEx = T

# if (GTEx && !exists("all_recount2_split_reads")) {
#   all_recount2_split_reads <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/annotated_split_reads.rda")
# } else {
#   print("'all_recount2_split_reads' file already loaded!")
# }

## Load all GTEx tissues
#gtex_tissues <-  readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")

###############################
#### 1. FUNCTIONS  ############
###############################





get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}



get_base_data_annotated <- function(cluster,
                                    samples,
                                    all_split_reads,
                                    split_read_counts,
                                    edb = NULL,
                                    gtf_path = NULL,
                                    gtf_version,
                                    blacklist_path,
                                    folder_save_path) {
  
  
  ## LOAD GTF DATA  -------------------------------------------------------------------------------------------------
  

  if (is.null(edb)) {
    if (is.null("gtf_path")) {
      ## Download the current GTF version
      if (is.null(gtf_version)) {
        gtf_version <- 105
        edb <-ensembldb::ensDbFromGtf(gtf = URLencode("http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.", gtf_version, ".gtf.gz"), 
                                      outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
        edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
        
      } else { ## Download the version specified
        edb <-ensembldb::ensDbFromGtf(gtf = URLencode(paste0("http://ftp.ensembl.org/pub/release-", gtf_version,
                                                             "/gtf/homo_sapiens/Homo_sapiens.GRCh38.", gtf_version, ".gtf.gz")), 
                                      outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
        edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
      }
      
    } else {
      edb <- ensembldb::ensDbFromGtf(gtf_path, outfile = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
      edb <- ensembldb::EnsDb(x = file.path(tempdir(), "Homo_sapiens.GRCh38.sqlite"))
    }
  }
  
  
  
  
 
  # tissue <- "Cells-Leukemiacellline_CML" 
  # tissue <- "Brain-FrontalCortex_BA9"
  # tissue <- clusters[1]
  
  ## Print the current tissue we are working with
  print(paste0(Sys.time(), " --> All split reads: ", length(all_split_reads$junID)))
  print(paste0(Sys.time(), " --> ", cluster, ". Split reads: ", length(split_read_counts$junID)))
 
  
  ################################################################
  ########### ANNOTATE SPLIT READS AND SAVE THE RESULT ###########
  ################################################################

  ## Filter by the current cluster's junctions and convert into GRanges format
  all_split_reads_cluster <- all_split_reads %>% 
    as_tibble() %>%
    dplyr::filter(junID %in% split_read_counts$junID) %>% 
    GenomicRanges::GRanges()
  
  
  print(paste0(Sys.time(), " --> ", cluster, ". Split reads in recount2: ", length(all_split_reads_cluster)))
  

  
  ## Annnotate the split reads using dasper
  all_split_reads_details <- annotate_dasper(all_split_reads = all_split_reads_cluster, 
                                             edb = edb,
                                             blacklist_path = blacklist_path)
  
  
  all_split_reads_details <- all_split_reads_details %>% as_tibble()
  
  saveRDS(object = all_split_reads_details %>% dplyr::select(junID, type, width), 
          file = paste0(folder_save_path, "/", cluster, "_annotated_SR_details_",
                        gtf_version, ".rds"))
  
  
  print(paste0(nrow(all_split_reads_details), " junctions NOT overlapping ambiguous genes!"))
  
  ## Apply the reads length filter
  all_split_reads_details <- apply_split_reads_length_filter(all_split_reads_details)
  
  ## Only keep annotated, novel donor and novel acceptor
  all_split_reads_details <- all_split_reads_details %>%
    as_tibble() %>%
    dplyr::filter(type %in% c("annotated", "novel_donor", "novel_acceptor")) %>%
    dplyr::select(seqnames, start, end, width, strand,
                  junID, ss5score, ss3score, type,
                  gene_name_junction, gene_id_junction, tx_id_junction)
  
  saveRDS(object = all_split_reads_details, 
          file = paste0(folder_save_path, "/", cluster, "_annotated_SR_details_length_", 
                        gtf_version, ".rds"))
  
  
  
  
  rm(all_split_reads_details)
  rm(all_split_reads_cluster)
  rm(split_read_counts)
  rm(samples)
  gc()
 
}





#' Title annotate_dasper
#' 
#' @param recount2_data recount2 junctions in a Genomic Ranges format
#' @param gtf path to the gtf file
#'
#' @return
#' @export
#' junctions annotated using dasper - it also removes ENCODE Blacklist regions and ambiguous genes
#' @examples
annotate_dasper <- function(all_split_reads,
                            edb,
                            blacklist_path) {
  
  
  
  
  
  # if (length(junc_ids) != length(recount2_data_tidy)) {
  #   
  #   print(paste0(length(junc_ids) - length(recount2_data_tidy), 
  #                " junctions found in haplotype areas!"))
  #   
  # }
  
  ## Removing all FC junctions that overlap with an ENCODE blacklist region
  junc_tidy <- remove_encode_blacklist_regions(GRdata = all_split_reads,
                                               blacklist_path = blacklist_path)
  
  ## CALLING DAVID'S FUNCTION
  all_split_reads_details_104_w_symbol <- dasper::junction_annot(junc_tidy %>% GRanges(), 
                                                                 ref = edb,
                                                                 ref_cols = c("gene_id", "gene_name", "symbol", "tx_id"), 
                                                                 ref_cols_to_merge = c("gene_id", "gene_name", "tx_id"))
  
  
  print(paste0(all_split_reads_details_104_w_symbol$type %>% unique(), " junction categories."))
 
  ## Return tidy data
  return(all_split_reads_details_104_w_symbol)
}



## REQUIRES THE hg ENCODE Blacklist regions

#' Title 'remove_encode_blacklist_regions'
#'
#' @param GRdata 
#' Genomic Ranges (GR) object
#' @return
#' @export 
#' from a list of junctions, it removes all junctions that overlap with an ENCODE Blacklist region
#' @examples
remove_encode_blacklist_regions <- function(GRdata,
                                            blacklist_path) {
  
  
  if (!exists("encode_blacklist_hg38")) {
    encode_blacklist_hg38 <- rtracklayer::import(con = blacklist_path) %>% diffloop::rmchr()
  } else {
    print("'encode_blacklist_hg38' file already loaded!")
  }
  
  overlaped_junctions <- GenomicRanges::findOverlaps(query = encode_blacklist_hg38, 
                                                     subject = GRdata %>% diffloop::rmchr(),
                                                     ignore.strand = F)
  
  ## JuncID indexes to be removed: they overlap with a black region
  indexes <- S4Vectors::subjectHits(overlaped_junctions)
  
  if (length(indexes) > 0) {
    print(paste0(length(unique(indexes)), " junctions overlap with a ENCODE blacklist region"))
    GRdata <- GRdata[-indexes, ]
    print(paste0(length(GRdata), " junctions after removing overlaps with ENCODE BlackList regions!"))
  }else{
    print("No junctions overlapping with an ENCODE blacklist region")
  }
  
  ## Return tidied data
  return(GRdata)
}






#' Title
#' apply_split_reads_length_filter
#' @param all_split_reads_details object containing recount2 dasper-annotated junctions
#' all_split_reads_details <- readRDS(file = "Brain-FrontalCortex_BA9_annotated_SR_details.rds")
#' @param min minimum split read length (min <- 25)
#' @param max maximum split read length
#' @return 
#' removes all junctions shorter or longer than the specified parameters
#' @export
#'
#' @examples
apply_split_reads_length_filter <- function(all_split_reads_details,
                                            min = 25,
                                            max = NULL) {
  
  
  if (is.null(max)) {
    max <- all_split_reads_details$width %>% max()  
  }
  
  all_split_reads_details <- all_split_reads_details %>%
    dplyr::filter(width >= min, width <= max) 
  
  return(all_split_reads_details)
}





#' Title
#' get_intron_length_ref_transcriptome
#' @param gtf 
#' path to the gtf reference transcriptome object
#' gtf = "/data/references/ensembl/gtf_gff3/v104/Homo_sapiens.GRCh38.104.gtf"
#' @param version 
#' version of the reference transcriptome
#' version = "v104"
#' @return returns the minimum and maximum intron length found within the reference transcriptome
#' @export
#'
#' @examples
get_intron_length_ref_transcriptome <- function(dirname = "/data/references/ensembl/gtf_gff3/v104/",
                                                gtf_file_name = "/v104/Homo_sapiens.GRCh38.104.gtf",
                                                version = "v104") {
  
  
  
  
  homo_sapiens <- ensemblGenome()
  basedir(homo_sapiens) <- dirname(dirname)
  read.gtf(object = homo_sapiens, filename = gtf_file_name)
  print(paste0("'homo_sapiens' ", version," file loaded!"))
  
  
  
  ## Obtaining the splicing table from the reference transcriptome
  # splice_table <- refGenome::getSpliceTable(homo_sapiens) %>% 
  #   refGenome::getGtf()
  # saveRDS(object = splice_table,
  #         file = paste0("/home/sruiz/PROJECTS/splicing-project/results/base_data/homosapiens_",version,"_splicetable.rds"))
  
  
  intron_lengths <- (splice_table$rstart) - (splice_table$lend) - 1
  
  df <- data.frame(length = intron_lengths)
  
  
  
  ## Plot the splicing data within a 'ggplot' histogram plot
  ggplot() + 
    geom_histogram(aes(x = length), dplyr::mutate(df, z = FALSE), bins = 40) + 
    geom_histogram(aes(x = length), dplyr::mutate(df %>% filter(length < 100), z = TRUE), binwidth = 5) +
    facet_zoom(xlim = c(0, 100), ylim = c(0, 12000), zoom.data = z, horizontal = FALSE) + 
    theme(zoom.y = element_blank(), validate = FALSE) +
    ggtitle(paste0("Length of the introns obtained from the raw data of the reference transcriptome",
                   "- version ", version))  +
    xlab("Intron length (in base pairs)") +
    ylab("Intron count (in number of introns)") + 
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "11"),
          axis.title = element_text(colour = "black", size = "11"))
  
  ## Save the plot
  ggplot2::ggsave("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/images/intron_length_raw_ensembl_v104.png",
                  width = 183, height = 183, units = "mm", dpi = 300)
  
  
  
  
  
  
  ## TIDY UP THE RAW REFERENCE TRANSCRIPTOME DATA ------------------------------------------------------------------
  
  ## Load the blacklist regions
  if (!exists("encode_blacklist_hg38")) {
    encode_blacklist_hg38 <- rtracklayer::import("/data/references/ENCODE_blacklist_v2/hg38-blacklist.v2.bed") %>% rmchr()
  } else {
    print("'encode_blacklist_hg38' file already loaded!")
  }
  
  
  ## Remove all genomic data from a low-supported transcript
  homo_sapiens_tidy <- homo_sapiens %>% 
    getGtf() %>%
    dplyr::filter(substr(transcript_support_level, 1, 1) %in% c("1", "2", "3")) %>%
    GenomicRanges::GRanges() 
  
  
  
  ## Remove all genomic data overlapping with an ENCONDE blacklist region 
  overlaped_junctions <- GenomicRanges::findOverlaps(query = encode_blacklist_hg38, 
                                                     subject = homo_sapiens_tidy, 
                                                     ignore.strand = F)
  
  
  # Get the indexes pairings for the overlapping elements.
  indexes <- S4Vectors::subjectHits(overlaped_junctions)
  homo_sapiens_tidy_GR <- homo_sapiens_tidy[-indexes, ]
  print(paste0(length(indexes), " Ensembl-v97 elements overlapped with ENCODE Blacklist regions!"))
  
  
  
  # Convert to refGenome format
  homo_sapiens_tidy <- homo_sapiens_tidy_GR %>%
    as.data.frame() %>%
    dplyr::rename(seqid = seqnames)
  
  # refGenome::setGtf(object = homo_sapiens,
  #                   value = homo_sapiens_tidy)
  # 
  # ## Getting the splice table
  # splice_table <- refGenome::getSpliceTable(homo_sapiens) %>% 
  #   refGenome::getGtf()
  
  saveRDS(object = splice_table,
          file = paste0("/home/sruiz/PROJECTS/splicing-project/results/base_data/homosapiens_",version,"_splicetable_v104.rds"))
  
  intron_lengths <- (splice_table$rstart) - (splice_table$lend) - 1
  
  
  
  ## Removing the frameshift introns (all introns with shorter or equal than 5 nucleotides)
  ind <- which(intron_lengths <= 5)
  if (ind %>% length > 0) {
    intron_lengths <- intron_lengths[-ind]
  }
  print(paste0(length(ind), " frameshift introns removed!"))
  
  
  ## Store the info in a dataframe
  df <- data.frame(length = intron_lengths)
  
  # df$length %>% min
  # df$length %>% max
  
  
  ## Plot the splicing data within a 'ggplot' histogram plot
  ggplot() + 
    geom_histogram(aes(x = length), dplyr::mutate(df, z = FALSE), bins = 40) + 
    geom_histogram(aes(x = length), dplyr::mutate(df %>% filter(length <= 100), z = TRUE), binwidth = 5) +
    geom_vline(aes(xintercept = 28), dplyr::mutate(df, z = TRUE), linetype = "dashed", colour = "red") +
    geom_text(aes(x = 28, label = "28 bp", y = 5000), dplyr::mutate(df, z = TRUE), colour = "black", check_overlap = T) +
    facet_zoom(xlim = c(0, 100), ylim = c(0, 8000), zoom.data = z, horizontal = FALSE) + 
    theme(zoom.y = element_blank(), validate = FALSE) +
    #ggtitle("Length of the introns obtained from the pre-processed data of Ensembl-v97")  +
    xlab("Intron length (in base pairs)") +
    ylab("Intron count (in number of introns)") + 
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"))
  
  ## Save the plot
  ggsave("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/images/intron_length_ensembl_v104_processed.png",
         width = 183, height = 183, units = "mm", dpi = 300)
  
  print("Ensembl-v97 statistics (after pre-processing it):")
  print(paste0("Maximum junction length: ", max(intron_lengths), " bp"))
  print(paste0("Minimum junction length: ", min(intron_lengths), " bp"))
  
  
}




#' Title
#' get_split_reads_length
#' @return returns two objects: 
#' 1. one object to contain the length of the split reads classified by split read category across all GTEx tissues.
#' 2. another object to contain the length of the split reads classifed by split read category and GTEx tissue.
#' @export
#'
#' @examples
get_split_reads_length <- function(gtf_version = 105) {
  
  projects_id <- sort(readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/all_projects_used.rds"))

  
  df_all_lengths <- map_df(projects_id, function(project_id) {
    
    # project_id <- projects_id[3]
    
    clusters <- sort(readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                           project_id, "/raw_data/all_clusters_used.rds")))
    
    map_df(clusters, function(cluster) { 
      
      # cluster <- clusters[1]
      
      print(cluster)
      
      ## LOAD ALL SPLIT READS ANNOTATED FOR THE CURRENT TISSUE
      
      all_split_reads_details <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project_id, 
                                                       "/results/base_data/", cluster, "/", cluster, "_annotated_SR_details_", gtf_version, ".rds"))
      
      all_split_reads_details <- all_split_reads_details %>% as_tibble() 
      
      
      ## ANNOTATED 
      annotated <- all_split_reads_details %>%
        filter(type == "annotated") %>%
        dplyr::select(junID, width)
      
    
      ## NOVEL DONOR 
      donor <- all_split_reads_details %>%
        filter(type == "novel_donor") %>%
        dplyr::select(junID, width)
      
      
      ## NOVEL ACCEPTOR
      acceptor <- all_split_reads_details %>%
        filter(type == "novel_acceptor") %>%
        dplyr::select(junID, width)
      
      
      ## NOVEL COMBO
      combo <- all_split_reads_details %>%
        filter(type == "novel_combo") %>%
        dplyr::select(junID, width)
      

      ## COMPLETELY UNANNOTATED
      none <- all_split_reads_details %>%
        filter(type == "unannotated") %>%
        dplyr::select(junID, width)
      
      
      rm(all_split_reads_details)
      
      ## JOIN ALL DATA
      df <- data.frame(junID = none$junID,
                       type = "completely unannotated",
                       length = none$width)
      df <- rbind(df,
                  data.frame(junID = c(donor$junID, acceptor$junID, combo$junID),
                             type = "partially annotated",
                             length = c(donor$width, acceptor$width, combo$width)))
      df <- rbind(df,
                  data.frame(junID = annotated$junID,
                             type = "annotated",
                             length = annotated$width))
      
      return(df %>% mutate(cluster = cluster))
    })
    
  })
  

  saveRDS(object = df_all_lengths %>% distinct(junID, .keep_all = T),
          file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/junction_all_lengths.rds")
  
  print("File saved!")
  
}




#' Title
#'
#' @param plot_acum_all_tissues 
#' @param plot_one_tissue 
#' @param plot_category_all_tissues 
#' @param tissue 
#'
#' @return
#' @export
#'
#' @examples
plot_junc_length <- function() {

    
    ## LOAD SOURCE DATA
    df_all_lengths <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/junction_all_lengths.rds")
    df_all_lengths %>% head()

    
    ########################
    ## STATS              ##
    ########################
    
    ## SHORTER THAN 25 bp ##
    min_pu <- df_all_lengths %>% 
      dplyr::filter(type == "partially annotated") %>%
      dplyr::pull(length) 
    lv <- length(which(min_pu < 25))
    lt <- length(min_pu)
    print(paste0("Total PUJs ", lt, ". Shorter 25bp: ", lv, " = ", (lv * 100) / lt, "%"))
    print(paste0("PUJs length mode: ", get_mode(min_pu)))
    
    min_a <- df_all_lengths %>% 
      filter(type == "annotated") %>%
      pull(length) 
    lv <- length(which(min_a < 25))
    lt <- length(min_a)
    print(paste0("Total annotated: ", lt, ". Shorter 25bp: ", lv, " = ", (lv * 100) / lt, "%"))
    print(paste0("Annotated length mode: ", get_mode(min_a)))
    
    min_cu <- df_all_lengths %>% 
      filter(type == "completely unannotated") %>%
      pull(length) 
    lv <- length(which(min_cu < 25))
    lt <- length(min_cu)
    print(paste0("Total none: ", lt, ". Shorter 25bp: ", lv, " = ", (lv * 100) / lt, "%"))
    print(paste0("None length mode: ", get_mode(min_cu)))
    
    ########################
    ## SHORTER THAN 28 bp ##
    ########################
    
    min_pu <- df_all_lengths %>% 
      dplyr::filter(type == "partially annotated") %>%
      dplyr::pull(length) 
    lv <- length(which(min_pu < 28))
    lt <- length(min_pu)
    print(paste0("Total PUJs ", lt, ". Shorter 28bp: ", lv, " = ", (lv * 100) / lt, "%"))
    print(paste0("PUJs length mode: ", get_mode(min_pu)))
    
    min_a <- df_all_lengths %>% 
      filter(type == "annotated") %>%
      pull(length) 
    lv <- length(which(min_a < 28))
    lt <- length(min_a)
    print(paste0("Total annotated: ", lt, ". Shorter 28bp: ", lv, " = ", (lv * 100) / lt, "%"))
    print(paste0("Annotated length mode: ", get_mode(min_a)))
    
    min_cu <- df_all_lengths %>% 
      filter(type == "completely unannotated") %>%
      pull(length) 
    lv <- length(which(min_cu < 28))
    lt <- length(min_cu)
    print(paste0("Total none: ", lt, ". Shorter 28bp: ", lv, " = ", (lv * 100) / lt, "%"))
    print(paste0("None length mode: ", get_mode(min_cu)))
    
    ## HISTOGRAM PLOT
    df_all_lengths_tidy <- df_all_lengths %>% 
      filter(length <= 200) %>%
      mutate(type = factor(type, levels = c("completely unannotated", "partially annotated", "annotated"))) 
    saveRDS(object = df_all_lengths_tidy,
            file = "/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/junction_tidy_lengths.rds")
    
    ggplot(data = df_all_lengths_tidy) + 
      geom_histogram(aes(x = length, fill = type), 
                     bins = 180,
                     #dplyr::mutate(df_all_lengths, z = FALSE), 
                     alpha = 0.7, 
                     position = "identity") +
      #geom_histogram(aes(x = length, fill = type), binwidth = 5, 
      #               dplyr::mutate(df_all_lengths, z = TRUE), 
      #               alpha = 0.7, position = "identity") +
      ##ggforce::facet_zoom(xlim = c(0, 200), 
      ##                    ylim = c(0, 710000), 
      ##                    zoom.data = z, horizontal = FALSE) + 
      ## ggtitle(paste0("Length of the annotated, partially unannotated and completely unannotated junctions.\nAll samples used - ", tissue)) +
      xlab("Implied intron length (in bp)") +
      ylab("Number of unique split reads") +
      theme_light() +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12")) +
      guides(fill = guide_legend(title = "Split Reads Category: ",
                                 ncol = 3,
                                 nrow = 2)) +
      scale_fill_manual(values = c("#661100", "#009E73", "#999999"), 
                        name = "Category",
                        breaks = c("annotated", "partially annotated", "completely unannotated"),
                        labels = c("annotated", "partially unannotated", "completely unannotated")) +
      theme(legend.position = "top",
            legend.text = element_text(size = 11)) %>%
      return()
    
    ## Save plot
    file_name <- "/home/sruiz/PROJECTS/splicing-project-recount3/images/junction_length.png"
    ggplot2::ggsave(filename = file_name, width = 183, height = 183, units = "mm", dpi = 300)

}




## Supplementary functions -------------------------------------------------------------------------------------------------

get_num_junc_sample_tissue <- function() {
  
  ## Declare some variables
  df <- data.frame("sample" = character(), 
                   "sex" = character(),
                   "age" = character(), 
                   "annotated" = integer(),
                   "novel_donor" = integer(), 
                   "novel_acceptor" = integer(),
                   "novel_combo" = integer(), 
                   "exon_skip" = integer(),
                   "none" = integer())
  
  ## Load all GTEx tissues
  gtex_tissues <-  readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")
  
  ## Load all GTEx samples
  gtex_samples <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_used.rda")
  
  ## Load all GTEx samples details
  gtex_samples_details <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_samples_details.rda")
  
  
  ## Obtain the info per each tissue
  for (tissue in gtex_tissues) { # tissue <- "Brain-FrontalCortex_BA9"
    
    ## Print the current tissue we are working with
    print(tissue)
    
    ## Obtain all samples from the current tissue
    tissue_samples <- gtex_samples[[tissue]]
    
    ## LOADING THE COUNTS FILE
    split_read_counts <- data.table::fread(paste0("/data/recount/GTEx_SRP012682/gtex_full_split_read_count_table/", 
                                                  tissue, "/", tissue, ".csv"))
    
    ## LOADING THE ANNOTATED JUNCTIONS FILE
    all_split_reads_details_97 <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/base_data/",
                                                        tissue, "/", tissue, "_annotated_SR_details_length.rds"))
    
    ## GET ONLY LOW-SHARED JUNCTIONS
    low_shared_junctions <- get_low_shared_puj(tissue)
    
    ## Per each sample from the current tissue, we obtain all junctions, counts and ratios
    for (sample in tissue_samples) { # sample <- tissue_samples[1]
      
      if (length(which((colnames(split_read_counts) == sample) == TRUE)) == 1) { # Every sampleID is unique, so the result should be equal to 1
        
        ########################################
        ## Getting the sex and age of the donor
        ########################################
        donor <- NULL
        sex <- NULL
        age <- NULL
        for (person in names(gtex_samples_details)) {
          if (any(gtex_samples_details[[person]]$sample_recount_id == sample)) {
            donor <- person
          }
        }
        if (str_detect(gtex_samples_details[[donor]]$bigwig_file[1], "female")) {
          sex <- "female"
        }else if (str_detect(gtex_samples_details[[donor]]$bigwig_file[1], "male")) {
          sex <- "male"
        }
        age <- gtex_samples_details[[donor]]$age[1]
        
        indx <- which(split_read_counts[, colnames(split_read_counts) == sample, with = FALSE] >= 1)
        all_juncs <- split_read_counts[indx, "junID"][[1]]
        
        #######################################################################
        ########### FILTERING THE ANNOTATION OBJECT BY THE JUNCID #############
        #######################################################################
        
        all_split_reads_details_97_2 <- all_split_reads_details_97 %>%
          as.data.frame() %>%
          filter(junID %in% all_juncs)
        
        print(paste0(tissue, " - ", sample))
        print(paste0("Initial number of total junctions: ", length(split_read_counts$junID)))
        print(paste0("Final number of total junctions: ", length(all_split_reads_details_97_2$junID)))
        
        #####################
        ## ANNOTATED
        #####################
        
        junc_ann <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "annotated") %>%
          dplyr::pull("junID")
        
        print(paste0("ANNOTATED: ", length(junc_ann)))
        
        ##################
        ## NOVEL DONOR
        ##################
        
        junc_donor <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_donor") %>%
          dplyr::pull("junID")
        
        low_shared_donor <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_donor") %>%
          dplyr::filter(junID %in% low_shared_junctions) %>%
          dplyr::pull("junID")
        
        print(paste0("NOVEL DONOR: ", length(junc_donor), " - ", length(low_shared_donor)))
        
        ##################
        ## NOVEL ACCEPTOR
        ##################
        
        junc_acceptor <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_acceptor") %>%
          dplyr::pull("junID") 
        
        low_shared_acceptor <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_acceptor") %>%
          dplyr::filter(junID %in% low_shared_junctions) %>%
          dplyr::pull("junID")
        
        print(paste0("NOVEL ACCEPTOR: ", length(junc_acceptor), " - ", length(low_shared_acceptor)))
        
        #####################
        ## NOVEL COMBO
        #####################
        
        junc_combo <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_combo") %>%
          dplyr::pull("junID") 
        
        low_shared_combo <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_combo") %>%
          dplyr::filter(junID %in% low_shared_junctions) %>%
          dplyr::pull("junID")
        
        print(paste0("NOVEL COMBO: ", length(junc_combo), " - ", length(low_shared_combo)))
        
        ####################
        ## NOVEL EXON SKIP
        ####################
        
        junc_skip <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_exon_skip") %>%
          dplyr::pull("junID")
        
        low_shared_skip <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "novel_exon_skip") %>%
          dplyr::filter(junID %in% low_shared_junctions) %>%
          dplyr::pull("junID")
        
        print(paste0("EXON SKIP: ", length(junc_skip), " - ",  length(low_shared_skip)))
        
        #################################
        ## COMPLETELY UNANNOTATED
        #################################
        
        junc_none <- all_split_reads_details_97_2 %>%
          dplyr::filter(junc_cat == "none") %>%
          dplyr::pull("junID")
        
        print(paste0("NONE: ", length(junc_none)))
        
        ###########################################
        ## Save the results in a data.frame format
        ###########################################
        
        df <- rbind(df,
                    data.frame("sample" = sample, "sex" = sex,
                               "age" = age, "annotated" = length(junc_ann),
                               "novel_donor" = length(junc_donor),
                               "low_donor" = length(low_shared_donor),
                               "novel_acceptor" = length(junc_acceptor),
                               "low_acceptor" = length(low_shared_acceptor),
                               "novel_combo" = length(junc_combo),
                               "low_combo" = length(low_shared_combo),
                               "exon_skip" = length(junc_skip),
                               "low_exon_skip" = length(low_shared_skip),
                               "none" = length(junc_none)))
      }
      ## Free up some memory by removing no-longer needed variables
      remove(all_split_reads_details_97_2)
    }
    
    ## Remove/update some variables
    remove(split_read_counts)
    remove(all_split_reads_details_97)
    
    ## Create the folder
    folder_name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/junction_counts/", tissue)
    dir.create(file.path(folder_name), showWarnings = T)
    
    ## Save the R objects
    file_name <- paste0(tissue, "_junction_counts.rda")
    saveRDS(object = df, file = paste0(folder_name, "/", file_name))
    xlsx::write.xlsx(df, file = paste0(folder_name, "/", file_name, ".xlsx"), sheetName = "junctions", row.names = FALSE)
  }
}

plot_num_junc_sample_tissue <- function() {
  
  
  
  ## Load all GTEx tissues
  gtex_tissues <-  sort(readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda"))
  
  for (tissue in gtex_tissues) { 
    # tissue <- "Brain-CerebellarHemisphere"
    # tissue <- "Brain-Cerebellum"
    # tissue <- "Brain-FrontalCortex_BA9"
    
    file.name <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/junction_counts/",
                        tissue, "/", tissue, "_Junction_Counts_2.rda")
    
    if (!file.exists(file.name)) {
      next;
    }
    
    df <- readRDS(file = file.name)
    
    df_counts <- data.frame(num_junctions = df$annotated,
                            type = "annotated")#,
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$novel_donor,
                                  type = "novel_donor"))#,
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$novel_acceptor,
                                  type = "novel_acceptor"))#,
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$novel_combo,
                                  type = "novel_combo"))#,
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$exon_skip,
                                  type = "exon_skip"))#,
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$none,
                                  type = "none"))
    
    annotated <- df_counts %>%
      dplyr::filter(type == "annotated") %>%
      pull("num_junctions") %>%
      sum()
    ndonor <- df_counts %>%
      dplyr::filter(type == "novel_donor") %>%
      pull("num_junctions") %>%
      sum()
    nacceptor <- df_counts %>%
      dplyr::filter(type == "novel_acceptor") %>%
      pull("num_junctions") %>%
      sum()
    ncombo <- df_counts %>%
      dplyr::filter(type == "novel_combo") %>%
      pull("num_junctions") %>%
      sum()
    nexonskip <- df_counts %>%
      dplyr::filter(type == "exon_skip") %>%
      pull("num_junctions") %>%
      sum()
    none <- df_counts %>%
      dplyr::filter(type == "none") %>%
      pull("num_junctions") %>%
      sum()
    
    print(annotated)
    print(ndonor)
    print(nacceptor)
    print(ncombo)
    print(nexonskip)
    print(none)
    total <- sum(c(annotated, ndonor, nacceptor, ncombo, nexonskip, none))
    
    print(annotated * 100 / total)
    print(ndonor * 100 / total)
    print(nacceptor * 100 / total)
    print(ncombo * 100 / total)
    print(nexonskip * 100 / total)
    print(none * 100 / total)
    
    print((annotated * 100 / total) +
            (ndonor * 100 / total) +
            (nacceptor * 100 / total) +
            (ncombo * 100 / total) +
            (nexonskip * 100 / total) +
            (none * 100 / total))
    
    df_counts <- data.frame(num_junctions = df$annotated,
                            type = "annotated")
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$novel_donor,
                                  type = "PUJs"))
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$novel_acceptor,
                                  type = "PUJs"))
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$novel_combo,
                                  type = "PUJs"))
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$exon_skip,
                                  type = "PUJs"))
    df_counts <- rbind(df_counts,
                       data.frame(num_junctions = df$none,
                                  type = "completely unannotated"))
    ggplot() +
      geom_boxplot(data = df_counts, aes(x = type, y = num_junctions, fill = type)) +
      ylab("Cumulative number of junctions") +
      xlab("Junction category") +
      theme_light() + 
      scale_fill_manual(#values=c("#ff5050", "#00cc66", "grey"), 
        name = "Junction category",
        values = c("#ff5050", "#00cc66", "grey"), 
        breaks = c("annotated", "PUJs", "completely unannotated"),
        labels = c("annotated", "PUJs", "completely unannotated")) +
      theme(legend.position = "top",
            legend.text = element_text(size = 11)) +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12"),
            axis.title.y.right = element_text(margin = unit(c(1), "cm"))) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    
    ## Save the plot
    folder_path <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/junction_counts/", tissue)
    dir.create(file.path(folder_path), showWarnings = T)
    ggsave(paste0(folder_path, "/", tissue, "_JunctionsPerJunctionCategory3.png"), 
           width = 183, height = 183, units = "mm", dpi = 300)
    
    #########################
    
    df_counts2 <- data.frame(num_junctions = df$novel_donor,
                             type = "novel donor")
    df_counts2 <- rbind(df_counts2,
                        data.frame(num_junctions = df$novel_acceptor,
                                   type = "novel acceptor"))#,
    df_counts2 <- rbind(df_counts2,
                        data.frame(num_junctions = df$novel_acceptor - df$novel_donor,
                                   type = "novel acceptor - novel donor"))
    ggplot() +
      geom_boxplot(data = df_counts2, aes(x = type, y = num_junctions, fill = type)) + 
      geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "red") +
      ylab("Cumulative number of junctions") + 
      xlab("Junction category") +
      theme_light() + 
      scale_fill_manual(#values=c("#ff5050", "#00cc66", "grey"), 
        name = "Junction category",
        values = c("#00cc66", "#00cc66", "#00cc66"), 
        breaks = c("novel donor", "novel acceptor", "novel acceptor - novel donor"),
        labels = c("novel donor", "novel acceptor", "novel acceptor - novel donor")) +
      theme(legend.position="top",
            legend.text = element_text(size = 11)) +
      theme(axis.line = element_line(colour = "black"), 
            axis.text = element_text(colour = "black", size = "12"),
            axis.title = element_text(colour = "black", size = "12"),
            axis.title.y.right = element_text(margin = unit(c(1), "cm"))) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    
    folder_path <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/junction_counts/", tissue)
    dir.create(file.path(folder_path), showWarnings = T)
    
    ggsave(paste0(folder_path, "/", tissue, "_JunctionsPerJunctionCategory_FCTX.png"), 
           width = 183, height = 183, units = "mm", dpi = 300)
  }
}

apply_kruskal_wallis_test <- function() {
  
  library(PMCMR)
  
  ## Read the file that contains all length of all junctions
  junction_length_per_tissue <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline1/junction_all_lengths_per_tissue.rda")
  
  ## Declare some variables
  annotated_length <- list()
  partially_unannotated_length <- list()
  novel_donor_length <- list()
  novel_acceptor_length <- list()
  novel_combo_length <- list()
  exon_skip_length <- list()
  
  ## Extract junction length information from the file
  for (tissue in names(junction_length_per_tissue)) { 
    if (str_detect(string = tissue, "Brain")) { ## tissue <- "Brain-FrontalCortex_BA9"
      print(tissue)
      
      annotated_length[[tissue]] <- junction_length_per_tissue[[tissue]]$annotated
      
      partially_unannotated_length[[tissue]] <- c(junction_length_per_tissue[[tissue]]$novel_donor,
                                                  junction_length_per_tissue[[tissue]]$novel_acceptor,
                                                  junction_length_per_tissue[[tissue]]$novel_combo,
                                                  junction_length_per_tissue[[tissue]]$exon_skip)
      
      novel_donor_length[[tissue]] <- junction_length_per_tissue[[tissue]]$novel_donor
      novel_acceptor_length[[tissue]] <- junction_length_per_tissue[[tissue]]$novel_acceptor
      novel_combo_length[[tissue]] <- junction_length_per_tissue[[tissue]]$novel_combo
      exon_skip_length[[tissue]] <- junction_length_per_tissue[[tissue]]$exon_skip
    }
  }
  
  ##################
  data_a_result <- list(kruskal.test = kruskal.test(annotated_length),
                        posthoc = PMCMR::posthoc.kruskal.dunn.test(annotated_length, p.adjust = "bonf"))
  saveRDS(object = data_a_result, 
          file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline1/kruskal_wallis_test_annotation.rda")
  print("Annotated saved!")
  rm(data_a_result)
  
  #########################
  data_pu_result <- list(kruskal.test = kruskal.test(partially_unannotated_length),
                         posthoc = PMCMR::posthoc.kruskal.dunn.test(partially_unannotated_length, p.adjust = "bonf"))
  saveRDS(object = data_pu_result, 
          file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline1/kruskal_wallis_test_partially_unannotated.rda")
  print("Partially unannotated saved!")
  rm(data_pu_result)
  
  ####################
  data_nd_result <- list(kruskal.test = kruskal.test(novel_donor_length),
                         posthoc = PMCMR::posthoc.kruskal.dunn.test(novel_donor_length, p.adjust = "bonf"))
  saveRDS(object = data_nd_result, 
          file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline1/kruskal_wallis_test_novel_donor.rda")
  print("Novel Donor saved!")
  rm(data_nd_result)
  
  ########################
  data_na_result <- list(kruskal.test = kruskal.test(novel_acceptor_length),
                         posthoc = PMCMR::posthoc.kruskal.dunn.test(novel_acceptor_length, p.adjust = "bonf"))
  print("Novel Acceptor saved!")
  saveRDS(object = data_na_result, 
          file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline1/kruskal_wallis_test_novel_acceptor.rda")
  rm(data_na_result)
  
  ##########
  data_es_result <- list(kruskal.test = kruskal.test(exon_skip_length),
                         posthoc = PMCMR::posthoc.kruskal.dunn.test(exon_skip_length, p.adjust = "bonf"))
  
  saveRDS(object = data_es_result, 
          file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline1/kruskal_wallis_test_exon_skip.rda")
  rm(data_es_result)
}

read_kruskal_wallis_test <- function() {
  
  library(stringr)
  library(ggforce)
  library(gridExtra)
  
  ############ PLOTTING RESULTS PER JUNCTION TYPE AND TISSUE ##################
  junction_length_per_tissue <- readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/pipeline1/junction_all_lengths_per_tissue.rda")
  
  ############################################
  ## PLOT THE ANNOTATED DIFFERENCES ##########
  ############################################
  
  type <- "annotation"
  kruskal_result <- readRDS(paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/kruskal_wallis_test_",
                                   type, ".rda"))
  
  x_label <- "Junction length (in nucleobases)"
  y_label <- "Number of junctions"
  
  df_c1 <- data.frame(tissue = character(), length = integer())
  for (tissue in c("Brain-CerebellarHemisphere", "Brain-FrontalCortex_BA9")) { 
    print(tissue)
    df_c1 <- rbind(df_c1,
                   data.frame(tissue = tissue,
                              length = junction_length_per_tissue[[tissue]]$annotated))
    
  }
  graph_title <- paste0("Length of the ", str_replace(type, "_", " "),
                        " junctions\nCerebellar Hemisphere and FCTX\nDunn's post-hoc test (p-value = 1.6e-06)")
  plot1 <- ggplot(df_c1, aes(x = length, group = tissue, color = tissue)) + 
    geom_freqpoly(binwidth = 10) +
    facet_zoom(xlim = c(0, 500), horizontal = FALSE, zoom.size = 4) +
    ggtitle(graph_title) +
    xlab(x_label) +
    ylab(y_label) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "11"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Tissue",
                                override.aes = list(size = 3)))
  
  
  df_c2 <- data.frame(tissue = character(), length = integer())
  for (tissue in c("Brain-CerebellarHemisphere", "Brain-Hypothalamus") ){
    print(tissue)
    df_c2 <- rbind(df_c2,
                   data.frame(tissue = tissue,
                              length = junction_length_per_tissue[[tissue]]$annotated))
    
  }
  graph_title <- paste0("Length of the ", str_replace(type, "_", " "), 
                        " junctions\nCerebellar Hemisphere and Hypothalamus\nDunn's post-hoc test (p-value = 1.9e-05)")
  plot2 <- ggplot(df_c2, aes(x = length, group = tissue, color = tissue)) + 
    geom_freqpoly(binwidth = 10) +
    facet_zoom(xlim = c(0, 500), horizontal = FALSE, zoom.size = 4) +
    ggtitle(graph_title) +
    xlab(x_label) +
    ylab(y_label) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "11"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Tissue",
                                override.aes = list(size = 3)))
  
  plot_grid <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)
  
  folder_path <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/images/")
  dir.create(file.path(folder_path), showWarnings = T)
  file_name <- paste0(folder_path, "Kruskal_", stringr::str_remove(type, "_"), "_Significance.png")
  ggplot2::ggsave(file_name, plot = plot_grid, width = 183, height = 183, units = "mm", dpi = 300)
  
  
  ############################################
  ## PLOT THE PARTIALLY UNANNOTATED ##########
  ############################################
  
  type <- "partially_unannotated"
  kruskal_result <- readRDS(paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/kruskal_wallis_test_",
                                   type, ".rda"))
  df <- data.frame(tissue = character(), length = integer())
  
  for (tissue in c("Brain-Amygdala", "Brain-Substantianigra", "Brain-Spinalcord_cervicalc-1")) {
    
    print(tissue)
    df <- rbind(df,
                data.frame(tissue = tissue,
                           length = junction_length_per_tissue[[tissue]]$novel_donor))
  }
  
  graph_title <- paste0("Length of the ", stringr::str_replace(type, "_", " "), 
                        " junctions\nAmygdala and Putamen_basalganglia\nDunn's post-hoc test (p-value = 1)")
  x_label <- "Junction length (in nucleobases)"
  y_label <- "Number of junctions"
  
  ggplot(df, aes(x = length, group = tissue, color = tissue)) + 
    geom_freqpoly(binwidth = 10, position = "identity") +
    facet_zoom(xlim = c(0, 500), horizontal = FALSE) +
    ggtitle(graph_title) +
    xlab(x_label) +
    ylab(y_label) +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "11"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(color = guide_legend(title = "Tissue",
                                override.aes = list(size = 3)))
  
  folder_path <- paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline1/images/")
  dir.create(file.path(folder_path), showWarnings = T)
  file_name <- paste0(folder_path, "Kruskal_", stringr::str_remove(type, "_"), "_Significance.png")
  ggsave(file_name, width = 183, height = 183, units = "mm", dpi = 300)
} 



###############################
#### CALLS  ###################
###############################

# get_base_data_annotated(gtf_path = "/data/references/ensembl/gtf_gff3/v76/Homo_sapiens.GRCh38.76.gtf",
#                         gtf_version = "v76")
# get_base_data_annotated(gtf_path = "/data/references/ensembl/gtf_gff3/v81/Homo_sapiens.GRCh38.81.gtf",
#                         gtf_version = "v81")
# get_base_data_annotated(gtf_path = "/data/references/ensembl/gtf_gff3/v90/Homo_sapiens.GRCh38.90.gtf",
#                         gtf_version = "v90")
# get_base_data_annotated(gtf_path = "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf",
#                         gtf_version = "v97")
# get_base_data_annotated(gtf_path = "/data/references/ensembl/gtf_gff3/v104/Homo_sapiens.GRCh38.104.gtf",
#                         gtf_version = "v104")


# get_base_data_annotated(gtf_path = "/data/references/ensembl/gtf_gff3/v104/Homo_sapiens.GRCh38.104.gtf",
#                         gtex_tissues = gtex_tissues[-c(7:17)],
#                         gtf_version = "v104")