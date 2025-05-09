library(tidyverse)
library(GenomicRanges)
library(DBI)
library(doParallel)

####################################################
## CONNECT TO THE SPLICING DATABASE ################
####################################################
## source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/manuscript_AD_figures.R")

## CONNECT TO THE DATABASE ------------------------------

base_folder <- here::here()
# base_folder <- "/mnt/PROJECTS/splicing-accuracy-manuscript"

gtf_version <- 105

main_project_identifier <- "SRP100948"
supportive_reads <- 1
data_subsample <- F

#project_name <- paste0(main_project_identifier, "_", supportive_reads, "read_subsample", data_subsample) #SRP100948_
project_name <- paste0(main_project_identifier, "_", supportive_reads, "read")

args <-
  list(
    project_id = main_project_identifier,
    case_type = "AD",
    control_type = "control",
    database_folder = file.path(base_folder, "database",  project_name, gtf_version),
    results_folder = file.path(base_folder, "results", project_name, gtf_version, "_paper_review", "results"),
    figures_folder = file.path(base_folder, "results", project_name, gtf_version, "_paper_review", "figures"),
    data_folder = file.path(base_folder, "results", project_name, gtf_version, "_paper_review", "data")
  )


dir.create(file.path(args$results_folder), recursive = TRUE, showWarnings = F)
dir.create(file.path(args$figures_folder), recursive = TRUE, showWarnings = F)
dir.create(file.path(args$data_folder), recursive = TRUE, showWarnings = F)



database_path <- file.path(args$database_folder, "/", paste0(project_name, ".sqlite"))
con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)


## QUERY MASTER TABLES --------------------------------------------------------------------

query = paste0("SELECT * FROM 'metadata'")
df_metadata <- dbGetQuery(con, query) %>% distinct(sample_id, .keep_all = T) %>% as_tibble()
all_projects <- df_metadata$SRA_project %>% unique
 
query <- paste0("SELECT * FROM 'intron'")
master_introns <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'novel'")
master_novel_junctions <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'transcript'")
master_transcripts <- dbGetQuery(con, query) %>% as_tibble()

query <- paste0("SELECT * FROM 'gene'")
master_gene <- dbGetQuery(con, query) %>% as_tibble()

all_projects <- df_metadata$SRA_project %>% unique
all_clusters <- df_metadata$cluster %>% unique()


## UTILS

get_mode <- function(data) {
  uniqv <- unique(data)
  uniqv[which.max(tabulate(match(data, uniqv)))]
}


custom_ggtheme <-  theme(text = element_text(size = 9, family="Arial", colour = "black"),
                         axis.ticks = element_line(colour = "black"),
                         axis.text = element_text(size = 9, family="Arial", colour = "black"),
                         axis.line = element_line(colour = "black"),
                         axis.title = element_text(size = 9, family="Arial", colour = "black"),
                         axis.text.y = element_text(size = 9, family="Arial", colour = "black"),
                         axis.text.x = element_text(size = 9, family="Arial", colour = "black", hjust = 0.5, vjust = 0.5),
                         strip.text = element_text(size = 9, family="Arial", colour = "black"),
                         legend.text = element_text(size = 8, family="Arial", colour = "black"),
                         legend.position = "top",
                         legend.box = "vertical")

########################################
## FUNCTIONS 
########################################


get_database_stats <- function() {
  
  df_metadata %>% nrow() 
  
  if ( any( names(df_metadata) == "rin" ) ) {
    df_metadata %>%
      filter(rin >= 8) %>% nrow()
    df_metadata %>%
      filter(rin >= 7) %>% nrow()
    df_metadata %>%
      filter(rin >= 6) %>% nrow()
  }
  
  ## This database included a set of 245,738 annotated introns (Ensembl-v105=)
  ## 149,649 of them with no evidence of mis-splicing and 96,089 introns with at least one linked novel split read)
  master_introns %>% distinct(ref_junID) %>% nrow()
  master_introns %>% dplyr::count(misspliced)
  
  
  ## and a linked set of 219,658 novel junctions (125,085 novel acceptor and 94,573 novel donor junctions),
  master_novel_junctions %>% nrow() 
  master_novel_junctions %>% dplyr::count(novel_type)
  
  
  ## originating from 23,999 genes and 181,284 transcripts
  master_transcripts %>% distinct(transcript_id) %>% nrow()
  master_gene %>% distinct(gene_id) %>% nrow()
  
}

## SECTION 1 - PAPER FIGURES -----------------------------------------------------------------

## Main figures

#' Title
#' Visualise differences in number of unique novel junctions from common introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
main_figure8_a <- function() {
  
  
  ## Load common introns
  common_introns_subsample <- readRDS(file = file.path(args$results_folder, "/common_subsampled_introns_seed1000.rds")) 
  common_introns_subsample %>% dplyr::count(sample_type)
  
  df_unique_junctions <- map_df(c(args$case_type,args$control_type), function(cluster_id) {
    
    # cluster_id <- all_clusters[1]
    print(cluster_id)
    ####################
    ## GET THE INTRONS
    ####################
    
    query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", args$project_id, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    query <- paste0("SELECT DISTINCT ref_junID FROM '", cluster_id, "_", args$project_id, "_misspliced'")
    introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())
    
    ## Keep only common, subsampled introns
    introns <- introns %>% filter(ref_junID %in% (common_introns_subsample %>% filter(sample_type == cluster_id) %>% pull(ref_junID)))
    
    ############################
    ## GET THE NOVEL JUNCTIONS
    ############################
    
    query <- paste0("SELECT * FROM '", cluster_id, "_", args$project_id, "_misspliced'")
    novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
    
    novel_junctions <- novel_junctions %>%
      inner_join(y = master_novel_junctions %>% dplyr::select(novel_junID,novel_type), by = "novel_junID") %>% 
      filter(ref_junID %in% (common_introns_subsample %>% filter(sample_type == cluster_id) %>% pull(ref_junID)))
    
    
    ###########################
    ## GET THE PROPORTIONS
    ###########################
    
    annotated_junc <- introns %>% distinct(ref_junID) %>% nrow()
    donor_junc <- novel_junctions %>% filter(novel_type == "novel_donor") %>% distinct(novel_junID) %>% nrow()
    acceptor_junc <- novel_junctions %>% filter(novel_type == "novel_acceptor") %>% distinct(novel_junID) %>% nrow()
    
    annotated_prop <- annotated_junc/(annotated_junc + donor_junc + acceptor_junc)
    donor_prop <- donor_junc/(annotated_junc + donor_junc + acceptor_junc)
    acceptor_prop <- acceptor_junc/(annotated_junc + donor_junc + acceptor_junc)
    
    ## Return the data.frame
    return(data.frame(cluster = cluster_id,
                      annotated_junc = annotated_junc ,
                      donor_junc = donor_junc,
                      acceptor_junc = acceptor_junc,
                      annotated_prop = annotated_prop,
                      donor_prop = donor_prop,
                      acceptor_prop = acceptor_prop))
  })
  
  df_unique_junctions_tidy <- df_unique_junctions %>%
    dplyr::select(cluster, donor = donor_prop, acceptor = acceptor_prop, annotated_intron = annotated_prop) %>%
    tidyr::gather(key = "type", value = "prop", -cluster )  %>%
    mutate(prop = prop * 100)
  
  ## Save source data
  write_csv(x = df_unique_junctions_tidy, file = file.path(args$data_folder, "figure8_a.csv"), col_names = T)
  
  ######################
  ## BAR PLOT
  ######################
  
  ggplot(data = df_unique_junctions_tidy %>%
           mutate(type = factor(type, levels = c("donor", "acceptor", "annotated_intron")),
                  prop_label = paste0(round(x = prop, digits = 2), "%")),
         aes(x = cluster, y = prop, fill =  type)) + 
    geom_bar(stat="identity", position = "stack")+
    geom_text(aes(label=prop_label), vjust=2, position = "stack", color="white", size = 3, fontface = "bold") +
    theme_light() +
    ylab("% unique junctions") +
    xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          legend.position = "top") +
    scale_fill_manual(values = c( "#999999","#64037d", "#35B779FF"),
                      breaks = c("annotated_intron","acceptor", "donor" ),
                      label = c("Annotated","Acceptor", "Donor" )) + 
    guides(fill = guide_legend(title = "", ncol = 3, nrow = 1 )) +
    custom_ggtheme  + 
    theme( plot.margin = margin(t = 5, r = 5, l = 5, b = -5),
           legend.box.margin=margin(b = -10, t = -5))#+
  
  
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_a.png"), width = 80, height = 70, units = "mm", dpi = 300)
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_a.svg"), width = 80, height = 70, units = "mm", dpi = 300)
  
  
  
}


#' Title
#' Visualise differences in the cummulative number of novel reads in common introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
main_figure8_b <- function() {
  
  ## Load common introns
  common_introns_subsample <-  readRDS(file = file.path(args$results_folder, "/common_subsampled_introns_seed1000.rds"))
  
  df_mean_counts <- map_df(all_clusters, function(cluster_id) {
    
    # cluster_id <- all_clusters[1]
    
    print(cluster_id)
    
    ####################
    ## GET EXPRESSION LEVELS FROM THE ANNOTATED INTRONS
    ####################
    
    query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", args$project_id, "_nevermisspliced'")
    introns <- dbGetQuery(con, query) %>% as_tibble()
    
    query <- paste0("SELECT DISTINCT ref_junID, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", args$project_id, "_misspliced'")
    introns <- rbind(introns, dbGetQuery(con, query) %>% as_tibble())

    ## Only common, subsampled introns
    introns <- introns %>%
      filter(ref_junID %in% (common_introns_subsample %>% filter(sample_type == cluster_id) %>% pull(ref_junID)))
    
    ###########################
    ## GET THE NOVEL JUNCTIONS
    ###########################
    
    query <- paste0("SELECT ref_junID, novel_junID, novel_sum_counts, novel_n_individuals FROM '", cluster_id, "_", args$project_id, "_misspliced'")
    novel_junctions <- dbGetQuery(con, query) %>% as_tibble() 
    
    
    novel_junctions <- novel_junctions %>%
      inner_join(y = master_novel_junctions %>% dplyr::select(novel_junID, novel_type) %>% as_tibble(), by = "novel_junID") %>% 
      as_tibble() %>%
      ## Only NOVEL READS FROM common, subsampled introns
      filter(ref_junID %in% introns$ref_junID)
     
    
    ###########################
    ## GET THE PROPORTIONS
    ###########################
    
    
    annotated <- introns %>% dplyr::distinct(ref_junID, .keep_all = T) %>% pull(ref_sum_counts) %>% sum()
    
    acceptor <- novel_junctions %>%
      filter(novel_type == "novel_acceptor") %>%
      dplyr::distinct(novel_junID, .keep_all = T) %>%
      pull(novel_sum_counts) %>% 
      sum()
    
    donor <- novel_junctions %>%
      filter(novel_type == "novel_donor") %>%
      dplyr::distinct(novel_junID, .keep_all = T) %>%
      pull(novel_sum_counts) %>% 
      sum()
    
    annotated_p = annotated * 100 / (annotated + acceptor + donor)
    acceptor_p = acceptor * 100 / (annotated + acceptor + donor)
    donor_p = donor * 100 / (annotated + acceptor + donor)
    
    
    return(data.frame(cluster = cluster_id,
                      type = c("annotated","acceptor", "donor"),
                      prop = c(annotated_p, acceptor_p, donor_p),
                      sum_counts = c(annotated, acceptor, donor)))
    
  })
  
  df_mean_counts <- df_mean_counts %>% mutate(type = factor(type, levels = c("donor", "acceptor", "annotated")),
                                              prop_label = paste0(round(x = prop, digits = 2), "%"))

  ## Save source data
  write_csv(x = df_mean_counts, file = file.path(args$data_folder, "figure8_b.csv"), col_names = T)
  
  
  ######################
  ## BAR PLOT
  ######################
  
  p <- ggplot(data = df_mean_counts, aes(x = cluster, y = prop, fill =  type)) + 
    geom_bar(stat = "identity") +
    ggforce::facet_zoom(ylim = c(97, 100)) +
    geom_text(aes(label=prop_label), vjust = 1.3,  position = "stack", color="white", size=3, fontface = "bold") +
    theme_light() +
    ylab("% cumulative split read counts") +
    xlab("") +   
    custom_ggtheme  +
    theme(axis.line = element_line(colour = "black"), legend.position = "top")  +
    scale_fill_manual(values = c( "#999999","#64037d", "#35B779FF"),
                      breaks = c( "annotated","acceptor", "donor"),
                      label = c("Annotated","Acceptor", "Donor")) + 
    guides(fill = guide_legend(title = "", ncol = 3, nrow = 1 )) + 
    theme(plot.margin = margin(t = 5, r = 5, l = 5, b = -5), legend.box.margin=margin(b = -10, t = -5))
  
  ## Only annotated zoomed graph
  pb <- ggplot_build(p)
  pb$data[[2]][1:6, 'alpha'] <- 0
  pg <- ggplot_gtable(pb)
  
  
  ## save plot
  png(filename = file.path(args$figures_folder, "figure8_b.png"), width = 100, height = 70,  units = "mm", res = 300)
  plot(pg)
  dev.off()

  
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_b.svg"), plot =  pg, width = 100, height = 70,  units = "mm", dpi = 300)
  
}

#' Title
#' Visualise differences in modulo3 values in introns from AD vs control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
main_figure8_c <- function() {
  
  common_introns_subsample <- readRDS(file = file.path(args$results_folder, "/common_subsampled_introns_seed1000.rds"))
  
  ## Calculate modulo3 from MANE transcripts in 100bp distance
  df_modulo <- map_df(all_clusters, function(cluster_id) {
      
    # cluster_id <- all_clusters[1]
    print(paste0(Sys.time(), " - ", cluster_id))
    
    #####################################
    ## GET NOVEL JUNCTIONS
    #####################################
    
    query <- paste0("SELECT novel_junID FROM '", cluster_id, "_", args$project_id, "_misspliced'")
    novel_jnx <- dbGetQuery(con, query) 
      
    query <- paste0("SELECT ref_junID, novel_junID, distance FROM 'novel' WHERE novel_junID IN (", paste(novel_jnx$novel_junID, collapse = ","),")")
    novel_jnx <- novel_jnx %>% left_join(y = dbGetQuery(con, query), by = "novel_junID") %>% as_tibble() 
    
      
    #####################################
    ## GET INTRON AND TRANSCRIPT INFO
    #####################################
    
    ## Add the transcript and MANE info
    query <- paste0("SELECT intron.ref_junID, intron.protein_coding, transcript.MANE
                    FROM 'intron' INNER JOIN transcript ON intron.transcript_id = transcript.id
                    WHERE ref_junID IN (", paste(novel_jnx$ref_junID, collapse = ","),")")
    all_jxn <- novel_jnx %>%
      inner_join(y = dbGetQuery(con, query), by = "ref_junID") %>% as_tibble() 
    
    
    df_novel_tidy <- all_jxn %>%
      distinct(novel_junID, .keep_all = T) %>%
      filter(abs(distance) <= 100, MANE == 1, ref_junID %in% (common_introns_subsample$ref_junID)) %>% 
      mutate(modulo = abs(distance) %% 3) 
      
    df_novel_tidy %>% distinct(ref_junID) %>% nrow %>% print
    
    df_novel_tidy <- df_novel_tidy %>% 
      group_by(modulo) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      mutate(cluster = cluster_id)
    
    return(df_novel_tidy)
    
  })
  
  
  df_modulo_tidy <- df_modulo %>% mutate(freq = freq * 100)
  df_modulo_tidy$modulo = factor(df_modulo_tidy$modulo, levels = c( "0", "1", "2"))
  
  
  ## Save Source Data
  write_csv(x = df_modulo_tidy, file = file.path(args$data_folder, "figure8_c.csv"), col_names = T)
  
  
  ################
  ## BAR PLOT
  ################

  ggplot(data = df_modulo_tidy %>%
           mutate(freq_label = paste0(round(x = freq, digits = 2), "%"), modulo = factor(modulo, levels = c("2", "1", "0"))),
         aes(x = cluster, y = freq, fill = modulo)) +
    geom_bar(stat = "identity") +
    xlab("") +
    ylab("% of novel junctions") +
    geom_text(aes(label=freq_label), vjust=2, position = "stack", color="white", size=3, fontface = "bold") +
    ggsci::scale_fill_npg(breaks=c('0', '1', '2')) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = "Modulo: ", ncol = 3, nrow = 1 )) +
    theme(plot.margin = margin(t = 5, r = 5, l = 5, b = -3), legend.box.margin=margin(b = -10, t = -5)) 
  
  
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_c.png"), width = 70, height = 70, units = "mm", dpi = 300)
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_c.svg"), width = 70, height = 70, units = "mm", dpi = 300)
  
}

#' Title
#' GO, KEGG and REACTOME ENRICHMENT analysis of the introns showing increasing MSR values in AD compared to control samples
#' Only using common, subsampled introns
#' @return
#' @export
#'
#' @examples
main_figure8_d_e <- function() {
  
  ## Load common introns
  common_introns_subsample <- readRDS(file = file.path(args$results_folder, "/common_subsampled_introns_seed1000.rds")) 
  
  
  #############################
  # Select and arrange necessary columns
  
  MSR_introns <- map(c("MSR_D", "MSR_A"), function(MSR_type) {
    
    # MSR_type <- "MSR_D"
    
    common_introns_subsample_MSR <- common_introns_subsample %>%
      dplyr::select(MSR = all_of(MSR_type), subclass, sample_type) %>%
      arrange(subclass)
    
    # Filter and merge data for 'control' and 'AD' sample types in one step
    common_introns_subsample_MSR_tidy <- common_introns_subsample_MSR %>%
      filter(sample_type %in% c("control", "AD")) %>%
      spread(key = sample_type, value = MSR) %>%
      mutate(MSR_difference = AD - control)
    
    # Separate introns based on MSR difference
    introns_MSR_increasing <- common_introns_subsample_MSR_tidy %>%
      filter(MSR_difference > 0)
    
    introns_MSR_decreasing <- common_introns_subsample_MSR_tidy %>%
      filter(MSR_difference < 0)
    
    list(MSR_type, introns_increasing = introns_MSR_increasing, introns_decreasing = introns_MSR_decreasing)
  })
  
  #############################
  
  ## Get introns with increasing MSR values (these introns were also paired by similar expression levels)
  
  MSR_D_introns_increasing <- MSR_introns[[1]]$introns_increasing %>%
    left_join(y = common_introns_subsample %>% dplyr::select(ref_junID, subclass), by = "subclass")
  MSR_A_introns_increasing <- MSR_introns[[2]]$introns_increasing %>%
    left_join(y = common_introns_subsample %>% dplyr::select(ref_junID, subclass), by = "subclass")
  MSR_introns_increasing <- rbind(MSR_D_introns_increasing, MSR_A_introns_increasing) %>% distinct(ref_junID)
  
  
  ## GET MASTER INTRON DATA
  df_master_intron_tidy <- master_introns %>% 
    dplyr::select(ref_junID, transcript_id) %>%
    inner_join(y = master_transcripts %>% dplyr::select(id, gene_id), by =c("transcript_id" = "id")) %>%
    inner_join(y = master_gene %>% dplyr::select(id, gene_name, gene_id), by =c("gene_id" = "id")) 
  
  
  ## Add gene name to introns with increasing MSR_A values in AD compared to control
  genes_increasing_MSR <- MSR_introns_increasing %>%
    inner_join(y = df_master_intron_tidy, by = c("ref_junID")) %>%
    distinct(gene_name) %>% pull()
  
  genes_increasing_MSR %>% length()
  
  ## Get gene background data
  bg_genes <- common_introns_subsample %>% inner_join(y = df_master_intron_tidy, by =c("ref_junID")) %>% distinct(gene_name) %>% pull()
  
  ################################
  ## KEGG ENRICHMENT
  ################################
  
  library('org.Hs.eg.db')
  
  category_terms <- 50
  entrez_genes_increasing <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, genes_increasing_MSR, 'ENTREZID', 'SYMBOL')#
  entrez_genes_bg <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, bg_genes, 'ENTREZID', 'SYMBOL')
  
  
  # mapIds(org.Hs.eg.db, introns_increasing_msra, 'ENTREZID', 'SYMBOL')
  
  ekegg_MSR <- clusterProfiler::enrichKEGG(
    gene          = entrez_genes_increasing,
    organism      = "hsa",
    keyType       = "kegg",
    universe      = entrez_genes_bg, 
    pAdjustMethod = "fdr")
  
  
  
  plotKEGG <-  clusterProfiler::dotplot(ekegg_MSR %>% mutate(ONTOLOGY = "KEGG") %>%
                                          filter(!(str_detect(string = Description, pattern = "tyrosine"))) %>%
                                          filter(!(str_detect(string = Description, pattern = "Chemical carcino"))) %>%
                                          filter(!(str_detect(string = Description, pattern = "Protein processing in"))) %>%
                                          filter(!(str_detect(string = Description, pattern = "Amino sugar"))), 
                                        showCategory = category_terms, 
                                        split="ONTOLOGY")
  
  plotKEGG +
    scale_y_discrete(labels = 
                       function(x) stringr::str_wrap(x, width = 40)) +
    xlab("Gene Ratio") +
    ggforce::facet_row(ONTOLOGY~., scales = "free", space = "free") +
    #Coord_flip() +
    custom_ggtheme +
    theme(text = element_text(colour = "black", size = 7),
          legend.position = "top",
          legend.box="horizontal",
          plot.margin = margin(t = 0, b =0, r = 5, 0),
          legend.margin = margin(t = -5, b = -5, r = 5, l = 1),
          legend.box.margin=margin(t = -10,b = -5, r = 0,l = 0)) + 
    scale_size(range = c(1, 5))+
    guides(colour = guide_legend(title = "q: "),
           size = guide_legend(title = "Gene count: "),
           nrow = 2, ncol = 2 ) 
  
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_d.png"), width = 120, height = 70, units = "mm", dpi = 300)
  
  ## Save Source Data
  
  write_csv(x = ekegg_MSR %>% as.data.frame(), file = file.path(args$data_folder, "figure8_d.csv"), col_names = T)
  
  
 
  ################################
  ## GO ENRICHMENT
  ################################
  
  ego_MSR <- clusterProfiler::enrichGO(
    gene          = genes_increasing_MSR,
    universe      = bg_genes,
    keyType       = "SYMBOL",
    OrgDb         = "org.Hs.eg.db", ##Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
    ont           = "ALL",
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  category_terms <- 20
  
  clusterProfiler::dotplot(ego_MSR %>% filter(ONTOLOGY == "CC"), x = "GeneRatio", 
                           showCategory = category_terms, split="ONTOLOGY") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
    xlab("Gene Ratio") +
    ggforce::facet_col(ONTOLOGY~., scales = "free_y", space = "free") +
    custom_ggtheme +
    theme(legend.position = "top",legend.box = "horizontal",
          plot.margin = margin(t = -5,b = 0,l = 5,r = 5),
          legend.margin=margin(t = -5,b = -5, r = 5, 0),
          legend.box.margin=margin(t = -5,r = 10,b = -5,l = 0)) + 
    scale_size(range = c(1, 5)) +
    guides(size = guide_legend(title = "Gene Count: "), colour = guide_legend(title = "q: "))    
  
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_e.png"), width = 180, height = 75, units = "mm", dpi = 300)
  ggplot2::ggsave(file.path(args$figures_folder, "figure8_e.svg"), width = 180, height = 75, units = "mm", dpi = 300)
  
 
  
  ## Save Source Data
  write_csv(x = ego_MSR %>% as.data.frame(), file = file.path(args$data_folder, "figure8_e.csv"), col_names = T)
  
  
  
  
  
  ################################
  ## DISEASE ENRICHMENT
  ################################
  
  x <- DOSE::enrichDO(gene    = entrez_genes_increasing,
                      ont           = "DO",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 10,
                      #universe      = names(entrez_genes_bg),
                      qvalueCutoff  = 0.05,
                      readable      = F)
  head(x)
  
 
}

#' Title
#' Test for differences in MSR_D and MSR_A between the introns from AD vs control samples.
#' Using one-tailed paired Wilcoxon test
#' @return
#' @export
#'
#' @examples
paper_stats_MSR <- function() {
  
  common_introns_subsample <- readRDS(file = file.path(args$results_folder, "/common_subsampled_introns_seed1000.rds"))
  
  ## Pair introns by subsampling criteria
  
  ## donor
  common_introns_subsample_MSRD <- common_introns_subsample %>%
    dplyr::select(ref_junID, MSR_D, subclass, mean_coverage, sample_type) %>%
    arrange(subclass)
  
  
  
  ## acceptor
  common_introns_subsample_MSRA <- common_introns_subsample %>%
    dplyr::select(ref_junID, MSR_A, subclass, mean_coverage, sample_type) %>%
    arrange(subclass)
  
  
  #########################################
  ## TEST DONOR 
  #########################################
  
  
  wilcox.test(x = common_introns_subsample_MSRD %>% filter(sample_type == args$case_type) %>% pull(MSR_D),
              y = common_introns_subsample_MSRD %>% filter(sample_type == args$control_type) %>% pull(MSR_D),
              alternative = "greater", paired = T, correct = T)
  
  rstatix::wilcox_effsize(data = common_introns_subsample_MSRD,
                          formula = MSR_D ~ sample_type,
                          paired = T)
  
  ########################################
  ## TEST ACCEPTOR
  ########################################
  
  
  wilcox.test(x = common_introns_subsample_MSRA %>% filter(sample_type == args$case_type) %>% pull(MSR_A),
              y = common_introns_subsample_MSRA %>% filter(sample_type == args$control_type) %>% pull(MSR_A),
              alternative = "greater", paired = T, correct = T)
  
  rstatix::wilcox_effsize(data = common_introns_subsample_MSRA, formula = MSR_A ~ sample_type, paired = T)
  

}


#' Title
#' Get an overview of the metadata of the project
#' @return
#' @export
#'
#' @examples
supplementary_figure27_a_to_d <- function() {
  
  
  
  ## Num samples
  plot_num_samples <- ggplot(df_metadata %>% dplyr::count(SRA_project, cluster)) +
    geom_bar(aes(x = n, y = SRA_project, fill = cluster),
             stat = "identity", position = position_dodge()) + 
    theme_light() +
    scale_fill_hue() +
    labs(y = "", x = "Num. samples" ) +
    guides(fill = guide_legend(title = NULL, ncol = 2, nrow = 1)) +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    scale_fill_manual(values = c("#bfbfbf","#666666"),
                      breaks = c("control", args$case_type),
                      label = c("Control", args$case_type)) +
    custom_ggtheme
  
  plot_num_samples
  
  if ( any(names(df_metadata) == "gender") ) {
    ## Gender
    plot_gender <- ggplot(df_metadata %>% 
                            mutate(gender = gender %>% as.character()) %>%
                            dplyr::count(cluster, gender)) +
      geom_bar(aes(y = cluster, x = n, group = gender, fill = gender), alpha = 0.8,
               stat = "identity", position = "dodge") + 
      theme_light() +
      labs(y = "", x = "Num. samples" ) +
      scale_fill_manual(labels = c("Male", "Female"), values = c("1", "2"), palette=scales::hue_pal()) +
      guides(fill = guide_legend(title = "Gender: ", ncol = 2, nrow = 1)) +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
      custom_ggtheme
    
    plot_gender
    
  }
  
  if (  any(names(df_metadata) == "age") ) {
    ## AGE
    plot_age <- ggplot(df_metadata ) +
      geom_density(aes(x = age, fill = cluster), alpha = 0.7) + 
      theme_light() +
      labs(x = "AGE" ) +
      scale_fill_hue() +
      guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) +
      scale_fill_manual(values = c("#bfbfbf","#666666"),
                        breaks = c("control", args$case_type),
                        label = c("Control", args$case_type)) +  
      custom_ggtheme
  }
  
  if ( any(names(df_metadata) == "rin") ) {
    ## RIN
    plot_rin <- ggplot(df_metadata ) +
      geom_density(aes(x = rin, fill = cluster), alpha = 0.7) + 
      theme_light() +
      labs(x = "RIN" ) +
      scale_fill_hue() +
      guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) +
      scale_fill_manual(values = c("#bfbfbf","#666666"),
                        breaks = c("control", args$case_type),
                        label = c("Control", args$case_type)) +  
      custom_ggtheme
  }
  
  
  ## Mapped read depth
  plot_read_depth <- ggplot(df_metadata ) +
    geom_density(aes(x = all_mapped_reads, fill = cluster), alpha = 0.6) + 
    theme_light() +
    scale_fill_hue() +
    labs(x = "Mapped Read Count" ) +
    guides(fill = guide_legend(title = NULL, ncol = 3, nrow = 1)) +
    scale_fill_manual(values = c("#bfbfbf","#666666"),
                      breaks = c("control", args$case_type),
                      label = c("Control", args$case_type)) +   
    custom_ggtheme
  
  
  plot_read_depth
  
  if ( exists("plot_age") ) {
    ggpubr::ggarrange(plot_num_samples, plot_gender, plot_age, plot_read_depth,
                      labels = c("a","b","c","d"))
    
  } else if ( exists("plot_rin") ) {
    
    ggpubr::ggarrange(plot_num_samples, plot_gender, plot_rin, plot_read_depth,
                      labels = c("a","b","c","d"))
  } else {
    
    ggpubr::ggarrange(plot_num_samples, plot_read_depth,
                      labels = c("a","b"))
  }
  
  
  ggplot2::ggsave(file.path(args$figures_folder, "supplementary_figure27_ad.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
}



#' Title
#' Get common introns between cases and controls (ONLY FRONTAL CORTEX) and subsample them by read coverage
#' @return
#' @export
#'
#' @examples
supplementary_figure27_e_f <- function() {
  
  
  print("Getting common introns between sample cluster...")
  
  clusters <- df_metadata$cluster %>% unique()
  
  #################################
  ## GET COMMON INTRONS 
  ## ACROSS CONTROL AND AD SAMPLES
  #################################
  
  df_common_introns <- get_common_introns()
  # df_common_introns %>% group_by(sample_type) %>% distinct(ref_junID) %>% dplyr::count(sample_type)
  
  if ( df_common_introns %>% dplyr::count(ref_junID) %>% filter(n == 1) %>% nrow() > 0 ) {
    print("ERROR. Common annotated introns should appear only once per sample cluster.")
  }
  
  ##########################################################
  ## PLOT COVERAGE OF SHARED INTRONS BEFORE SUBSAMPLING
  ##########################################################
  
  plot_BS <- ggplot(data = df_common_introns) +
    geom_density(mapping = aes(x = mean_coverage, fill = sample_type), alpha = 0.8) +
    ggtitle("Before subsampling") +
    theme_light() +
    custom_ggtheme +
    xlab("log10 mean expression level") +
    scale_fill_manual(name = "Sample group: ",
                      values = c("#bfbfbf","#666666"),
                      breaks = c(args$control_type, args$case_type ),
                      label = c(args$control_type, args$case_type )) +
    ylim(c(0, 1))+
    xlim(c(0, 3))
  
  plot_BS
  
  
  #####################################
  #### START SUBSAMPLIING
  #####################################
  
  
  # "To compare differences in splicing accuracy between the annotated introns found across the 48 samples studied, we made use of their MSR values. 
  # To avoid potential biases derived from differences in mean expression levels, we only considered those annotated introns overlapping both groups of samples 
  # that displayed a maximum difference in their log10 expression levels of 0.005 (matchit() function, MatchIt R package, version 4.5.0). "
  
  
  subsample <- if ( file.exists(file.path(args$results_folder, "common_subsampled_introns_seed1000.rds")) ) {
    readRDS(file = file.path(args$results_folder, "/common_subsampled_introns_seed1000.rds"))
  } else {
    subsample_introns(df_common_introns)
  }
  # subsample %>% group_by(sample_type) %>% distinct(ref_junID) %>% dplyr::count(sample_type)

  #######################################
  ## PLOT COVERAGE AFTER SUBSAMPLE
  #######################################
  
  plot_AS <- ggplot(data = subsample) +
    geom_density(mapping = aes(x = mean_coverage, fill = sample_type), alpha = 0.8) +
    ggtitle("After subsampling") +
    theme_light() +
    custom_ggtheme +
    xlab("log10 mean expression level") +
    scale_fill_manual(name = "Sample group: ",
                      values = c("#bfbfbf","#666666"),
                      breaks = c(args$control_type, args$case_type ),
                      label = c(args$control_type, args$case_type )) +
    ylim(c(0, 1))+
    xlim(c(0, 3))
  
  ggpubr::ggarrange(plot_BS, plot_AS, labels = c("a", "b"), common.legend = T) 
  
  dir.create(path = figures_path, recursive = T,showWarnings = T)
  ggplot2::ggsave(file.path(args$figures_path, "supplementary_figure27_e_f.png"), width = 180, height = 100, units = "mm", dpi = 300)
  
}



##################################
## SUPPLEMENTARY FIGURE 22
##################################

## FUNCTIONS TO GENERATE SUPPLEMENTARY FIGURE 22 & Supplementary Table 13

get_data_supplementary_figure_22 <- function(replace = T) {
  
  #common_introns <- readRDS(file = paste0(results_path, "/common_subsampled_introns_seed1000.rds"))
  
  local_results_folder <- paste0(args$results_folder, "/MSR_normalisation_by_TPM/")
  dir.create(path = local_results_folder, showWarnings = F)
  local_figures_folder <-  paste0(args$figures_folder, "/MSR_normalisation_by_TPM/")
  dir.create(path = local_figures_folder, showWarnings = F)
  
  do_paired_test <- T
  
  ################################
  ## LOAD RBPs OF INTEREST
  ################################
  
  # Only get genes with splicing regulator, spliceosome, Exon.Junction.Complex and NMD functions
  
  all_RBPs <- (xlsx::read.xlsx(file = paste0(base_folder,'/dependencies/RBPs_subgroups.xlsx'), header = TRUE, sheetIndex = 1) %>% as_tibble() %>% mutate(NMD = 0))[-116,]
  
  all_NMDs <- data.frame(name = c("SMG1", "SMG5", "SMG6", "SMG7", "UPF1", "UPF2", "UPF3A"),
                         id = c("ENSG00000157106", "ENSG00000198952", "ENSG00000070366", "ENSG00000116698", "ENSG00000005007", "ENSG00000151461", "ENSG00000169062"),
                         Splicing.regulation = c(0,0,0,0,0,0,0),
                         Spliceosome = c(0,0,0,0,0,0,0),
                         Exon.Junction.Complex = c(0,0,0,0,0,0,0),
                         NMD = c(1,1,1,1,1,1,1))
  
  
  all_RBPs_tidy <- rbind(all_NMDs,all_RBPs)
  


  for (project_id in all_projects) {
    
    # project_id = all_projects[1]

    message(project_id, "...")
    
    if (!file.exists(paste0(local_results_folder, "/supplementary_figure_22.rds")) || replace) {
      
      ################################
      ## CALCULATE MSR VALUES 
      ## NORMALISED BY NMD/TPM VALUES
      ################################
      
      fold_change_TPM <- if (file.exists(paste0(local_results_folder, "/", project_id, "_inverse_fold-change_TPM.rds"))) {
         readRDS(file = paste0(local_results_folder, "/", project_id, "_inverse_fold-change_TPM.rds"))
      } else {
        fold_change_TPM <- AD_control_calculate_fold_change_TPM(project.id = project_id, RBP.list = all_RBPs_tidy, get.median = T )
        fold_change_TPM %>% arrange(t_test) %>% as_tibble()
        saveRDS(object = fold_change_TPM, file = paste0(local_results_folder, "/", project_id, "_inverse_fold-change_TPM.rds"))
        fold_change_TPM
      }
      
      MSR_normalised <- AD_control_normalise_MSR_values_by_gene_TPM(project.id = project_id, fold.change.TPM = fold_change_TPM, gene.list = all_RBPs_tidy)
      
      ## Get common introns between AD and control
      
      introns_control <- MSR_normalised %>% filter(groups == "control") %>% distinct(ref_junID)
      introns_AD <- MSR_normalised %>% filter(groups == "AD") %>% distinct(ref_junID)
      common_introns <- intersect(introns_AD$ref_junID, introns_control$ref_junID) %>% unique()
      common_introns %>% length()
      
      
      ################################
      ## CALCULATE AGE EFFECT SIZES 
      ## USING MSR NORMALISED DATA
      ################################
      
      message("Calculating AD effect sizes using MSR normalised data....")
      
      gene_list <- MSR_normalised$gene_normalised %>% unique %>% sort
      
      doParallel::registerDoParallel(cores = 4)
      effect_size_normalised_TPM <- foreach(j = seq(length(gene_list)), .combine = "rbind") %dopar% {
        
        gene = gene_list[j]
        
        # gene <- gene_list[1]
        
        message(gene, "....")
        
        
        map_df( c("MSR Donor", "MSR Acceptor"), function(MSR_type) {
          
          # MSR_type = "MSR Donor"
          
          message(MSR_type, "....")
          
          
          if (MSR_type == "MSR Donor") {
            MSR_normalised_data = MSR_normalised %>% dplyr::select(ref_junID, MSR_normalised = MSR_D_normalised, disease_group = groups, gene_normalised)
          } else {
            MSR_normalised_data = MSR_normalised %>% dplyr::select(ref_junID, MSR_normalised = MSR_A_normalised, disease_group = groups, gene_normalised)
          }
          
          MSR_normalised_data = MSR_normalised_data   %>%
            filter(ref_junID %in% (common_introns)) %>%
            filter(gene_normalised == gene) %>%
            spread(key = disease_group, MSR_normalised) 
          
          
          ## Wilcoxon one-tailed test
          wilcox_pval <- data.frame(pval = c((wilcox.test(x = MSR_normalised_data$control,
                                                          y = MSR_normalised_data$AD,
                                                          alternative = "less",
                                                          paired = T))$p.value))
          
          
          ## Wilcoxon effect size
          message(MSR_type, " - Wilcoxon paired effect size test ....")
          
          w_effsize <- rstatix::wilcox_effsize(data = MSR_normalised_data %>%
                                                 dplyr::select(ref_junID, control, AD) %>%
                                                 gather(key = disease_group, value = MSR, -ref_junID) %>%
                                                 mutate(disease_group = disease_group %>% as.factor()),
                                               formula = MSR ~ disease_group,
                                               alternative = "less",
                                               ref.group = "control",
                                               conf.level = .05,
                                               paired = T) %>%
            mutate(MSR_type = MSR_type, 
                   project = project_id,
                   gene_normalised = gene)
          
          
          
          
          # ## Save results
          # dir.create(path = paste0(local_results_folder, "/",project_id,"_paired",do_paired_test,"/"),
          #            recursive = T)
          # saveRDS(object = r_effect1 %>% cbind(wilcox_pval),
          #         file = paste0(local_results_folder, "/",project_id,"_paired",do_paired_test,"/", 
          #                       project_id,"_", gene, "_",MSR_type,"_effect_size_age_normalised_TPM.rds"))
          
          ## Return results
          return(w_effsize %>% cbind(wilcox_pval))
          
        })
      } 
      
      
      
      saveRDS(object = effect_size_normalised_TPM %>% inner_join(y = all_RBPs_tidy, by = c("gene_normalised" = "name")),
              file = paste0(local_results_folder, "/supplementary_figure_22.rds"))
      
      
      
    } else {
      effect_size_normalised_TPM <- readRDS(file = paste0(local_results_folder, "/supplementary_figure_22.rds"))
    }
  
  }
}



AD_control_calculate_fold_change_TPM <- function(project.id,
                                                 RBP.list,
                                                 get.median = T) {
  
  
  
  
  fold_change_TPM_values_age_groups <- map_df(RBP.list$id, function(RBP) {
    
    # RBP = RBP.list$id[1]
    # RBP = "ENSG00000271885"
    
    message(RBP,"...")
    
    control_group_TPM <- readRDS(file = paste0(base_folder,"/results/", project_name,"/", gtf_version, "/", project.id, "/tpm/", project.id, "_control_tpm.rds")) %>%
      filter(gene_id %in% RBP ) %>% 
      rowwise() %>% 
      mutate(TPM = mean( c_across(where(is.numeric))) ) %>%
      ungroup()
    
    
    AD_group_TPM <- readRDS(file = paste0(base_folder,"/results/",project_name,"/", gtf_version, "/", project.id, "/tpm/",
                                              project.id, "_AD_tpm.rds")) %>%
      filter(gene_id %in% RBP ) %>% 
      rowwise() %>% 
      mutate(TPM = mean( c_across(where(is.numeric))) ) %>%
      ungroup() 
    
    
    
    if (control_group_TPM %>% nrow == 1 && AD_group_TPM %>% nrow == 1) {
      
      fold_change = AD_group_TPM$TPM/control_group_TPM$TPM
      
      data.frame(gene = RBP,
                 mean_TPM_control = control_group_TPM$TPM,
                 mean_TPM_AD = AD_group_TPM$TPM,
                 fold_change = fold_change,
                 inverse_fold_change = 1/fold_change,
                 t_test = t.test(control_group_TPM[,-1] %>% gather() %>% pull(value),
                                 AD_group_TPM[,-1] %>% gather() %>% pull(value))$p.value,
                 type = ifelse(fold_change > 1, "upregulation with AD", 
                               ifelse(fold_change < 1, "downregulation with AD", 
                                      "no change in expression w AD"))) %>%
        return()
    
    } else {
      return(NULL)
    }
  
  })
  
  fold_change_TPM_values_age_groups %>% inner_join(y = RBP.list, by = c("gene" = "id")) %>% arrange(desc(fold_change)) %>% return()
  
}

AD_control_normalise_MSR_values_by_gene_TPM <- function(project.id,
                                                        fold.change.TPM,
                                                        gene.list) {
  
  
  
  
  local_results_folder <- paste0(args$results_folder, "/MSR_normalisation_by_TPM/")
  # 
  # # 2. Query database to get MSR values across age groups
  # introns_MSR_age_groups <- map_df(c("20-39","60-79"), function(cluster_id) {
  # 
  #   # cluster_id = "20-39"
  #   message(cluster_id, "...")
  #   con <- dbConnect(RSQLite::SQLite(), database_path)
  #   query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_misspliced'")
  #   tissue_introns <- dbGetQuery(con, query) %>% as_tibble() %>% distinct(ref_junID, .keep_all = T)
  #   
  #   query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_nevermisspliced' " )
  #   tissue_introns <- rbind(tissue_introns, dbGetQuery(con, query)) %>% as_tibble() %>% 
  #     distinct(ref_junID, .keep_all = T) %>%
  #     arrange(ref_junID)
  #   DBI::dbDisconnect(con)
  #   
  #   tissue_introns %>%
  #     mutate(age = cluster_id) %>%
  #     return()
  # })
  # 
  # 
  # 
  # 
  # introns_MSR_D_age_groups <- introns_MSR_age_groups %>%
  #   dplyr::select(-MSR_A) %>%
  #   mutate(MSR_D = ifelse(MSR_D == 0, 0.00000001, MSR_D)) %>%
  #   tidyr::spread(key = age, value=MSR_D ) %>%
  #   mutate(MSR_D_fold_change = log2(`60-79`/`20-39`)) %>%
  #   drop_na()
  # 
  # 
  # fold_change_TPM <- fold.change.TPM %>%
  #   mutate(fold_change = fold_change_eldest %>% log2)
  # 
  # 
  # wilcox.test(introns_MSR_D_age_groups$MSR_D_fold_change,
  #             fold_change_TPM %>% filter(Splicing.regulation == 1) %>% pull(fold_change))
  # 
  # 
  # wilcox.test(introns_MSR_D_age_groups$MSR_D_fold_change,
  #             fold_change_TPM %>% filter(NMD == 1) %>% pull(fold_change))
  # 
  # 
  # wilcox.test(introns_MSR_D_age_groups$MSR_D_fold_change,
  #             fold_change_TPM %>% filter(Spliceosome == 1) %>% pull(fold_change))
  # 
  # 
  # wilcox.test(introns_MSR_D_age_groups$MSR_D_fold_change,
  #             fold_change_TPM %>% filter(Exon.Junction.Complex == 1) %>% pull(fold_change))
  # 
  # 
  # introns_MSR_D_age_groups
  # wilcox.test(x = introns_MSR_D_age_groups$`20-39`,
  #             y = introns_MSR_D_age_groups$`60-79`,
  #             alternative = "less",
  #             correct = T,
  #             paired = T)
  # 
  # 
  # introns_MSR_A_age_groups <- introns_MSR_age_groups %>%
  #   dplyr::select(-MSR_D) %>%
  #   tidyr::spread(key = age, value = MSR_A)
  # 
  # wilcox.test(x = introns_MSR_A_age_groups$`20-39`,
  #             y = introns_MSR_A_age_groups$`60-79`,
  #             alternative = "less",
  #             correct = T,
  #             paired = T)
  # 
  
  
  
  
  # if ( !file.exists(paste0(local_results_folder, "/",project.id,"/",project.id,"_age_MSR_D_normalised_by_TPM.rds")) || 
  #      !file.exists(paste0(local_results_folder, "/",project.id,"/",project.id,"_age_MSR_D_normalised_by_TPM.rds")) ) {
  
  median_TPM_values_age_groups <- map_df( c("control","AD"), function(cluster_id) {
    
    # cluster_id <- c("control","AD")[1]
    
    message(cluster_id, "...")
    
    # cluster_id <- (df_metadata$cluster %>% unique())[1]
    
    ## 1. Calculate median TPM values corresponding to each RBP/NMD gene across the samples of each age group
    message(cluster_id, " - getting introns from database .... ")
    
    # if (file.exists(paste0(base_folder,"/results/splicing_1read/", gtf_version, "/", project.id, "/tpm/",
    #                        project.id, "_", cluster_id, "_tpm.rds"))) {
    
    
    # message("Calculating median TPM values across sample groups...")
    # 
    # local_TPM <- readRDS(file = paste0(base_folder,"/results/splicing_1read/", gtf_version, "/", project.id, "/tpm/",
    #                                    project.id, "_", cluster_id, "_tpm.rds")) %>%
    #   filter(gene_id %in% (gene.list$id) )
    # 
    # 
    # if (get.median) {
    #   local_TPM_w_median <- local_TPM %>% 
    #     rowwise() %>% 
    #     mutate(TPM = median( c_across(where(is.numeric))) )
    # } else {
    #   local_TPM_w_median <- local_TPM %>% 
    #     rowwise() %>% 
    #     mutate(TPM = mean( c_across(where(is.numeric))) )
    # }
    # local_TPM_w_median <- local_TPM_w_median %>%
    #   dplyr::select(gene_id, TPM)  %>%
    #   inner_join(y = gene.list, by = c("gene_id"="id")) %>%
    #   dplyr::rename(gene_name = name) %>%
    #   dplyr::relocate(gene_name, .after = "gene_id") %>%
    #   mutate(inverse_TPM = 1/TPM) %>%
    #   dplyr::relocate(inverse_TPM, .after = "TPM")
    
    
    # 2. Query database to get MSR values across age groups
    con <- dbConnect(RSQLite::SQLite(), database_path)
    query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_misspliced'")
    all_introns <- dbGetQuery(con, query) %>% as_tibble() %>% distinct(ref_junID, .keep_all = T)
    
    query <- paste0("SELECT ref_junID, MSR_D, MSR_A FROM '", cluster_id, "_", project.id, "_nevermisspliced' " )
    all_introns <- rbind(all_introns, dbGetQuery(con, query)) %>% as_tibble() %>% 
      distinct(ref_junID, .keep_all = T) %>%
      arrange(ref_junID)
    DBI::dbDisconnect(con)
    
    
    message("Calculating normalised MSR values by inverse fold-change TPM...")
    
    
    # 3. Normalise MSR values per RBP/NMD gene in the age group and tissue
    normalised_TPM <- map_df(fold.change.TPM$name, function(gene_name) {
      
      # message(gene_name, "...")
      # gene_name <- fold.change.TPM$name[1]
      
      if (cluster_id == "AD") {
        TPM <- fold.change.TPM %>% filter(name == gene_name) %>% pull(inverse_fold_change)
      } else {
        #TPM <- fold.change.TPM %>% filter(name == gene_name) %>% pull(mean_youngest)
        TPM <- 1
      }
      
      
      
      all_introns %>%
        mutate(MSR_D_normalised = MSR_D/TPM,
               MSR_A_normalised = MSR_A/TPM) %>%
        mutate(TPM_group = fold.change.TPM %>% filter(name == gene_name) %>% pull(fold_change),
               inverse_TPM_group = TPM,
               groups = cluster_id,
               gene_normalised = gene_name) %>%
        return()
      
    })
    
    
    
    
    #normalised_TPM[which(!is.finite(normalised_TPM$MSR_D_normalised)),"MSR_D_normalised"] <- 0
    #normalised_TPM[which(!is.finite(normalised_TPM$MSR_A_normalised)),"MSR_A_normalised"] <- 0
    
    normalised_TPM %>% 
      return()
    
    
    # } else {
    #   return(NULL)
    # }   
  })
  
  
  
  
  # ## 4. Get only the annotated introns overlapping the two age groups
  # common_ref_junID <- median_TPM_values_age_groups %>%
  #   dplyr::count(ref_junID) %>%
  #   filter(n == ((median_TPM_values_age_groups$gene_normalised %>% unique %>% length) * 2)) %>%
  #   pull(ref_junID)
  # 
  # 
  # ## 5. Filter the results by common introns between the two age groups
  # median_TPM_values_age_groups_common_introns <- median_TPM_values_age_groups %>%
  #   filter(ref_junID %in% common_ref_junID) 
  
  
  # ## 6. Calculate MSR_D normalised values & save results
  # MSR_D_normalised <- median_TPM_values_age_groups_common_introns %>%
  #   dplyr::select(-c(MSR_D, MSR_A, MSR_A_normalised)) %>%
  #   spread(key = age, MSR_D_normalised)
  # saveRDS(object = MSR_D_normalised,
  #         file = paste0(local_results_folder, "/",project.id,"/", project.id,"_age_MSR_D_normalised_by_TPM.rds"))
  # 
  # 
  # ## 7. Calculate MSR_A normalised values & save results
  # MSR_A_normalised <- median_TPM_values_age_groups_common_introns %>%
  #   dplyr::select(-c(MSR_D,MSR_A,MSR_D_normalised)) %>%
  #   spread(key = age, MSR_A_normalised)
  # saveRDS(object = MSR_A_normalised,
  #         file = paste0(local_results_folder, "/",project.id,"/", project.id,"_age_MSR_A_normalised_by_TPM.rds"))
  
  
  # } else {
  #   message(project.id, " - loading normalised MSR values....")
  #   MSR_D_normalised <- readRDS(file = paste0(local_results_folder, 
  #                                             "/", project.id, "/", project.id,"_age_MSR_D_normalised_by_TPM.rds"))
  #   MSR_A_normalised <- readRDS(file = paste0(local_results_folder, 
  #                                             "/", project.id, "/", project.id,"_age_MSR_A_normalised_by_TPM.rds"))
  # }
  
  return(median_TPM_values_age_groups)
  
}


supplementary_figure_22 <- function(plot.stats = F) {
  
  paired.test = T
  get.median = T
  project_id = all_projects[1]  

  local_results_folder <- paste0(args$results_folder, "/MSR_normalisation_by_TPM/")
  dir.create(path = local_results_folder, showWarnings = F)

  effect_size_file_path <- file.path(local_results_folder, "supplementary_figure_22.rds")
  
  
  #######################################################
  ## LOAD AND TIDY DATA
  #######################################################
  
  df_wilcoxon_AD <- readRDS(file = effect_size_file_path) %>% group_by(MSR_type) %>% mutate(q = p.adjust(pval, method = "fdr")) %>% ungroup()
  
  df_wilcoxon_AD$MSR_type = factor(df_wilcoxon_AD$MSR_type, levels = c(df_wilcoxon_AD$MSR_type %>% unique))
  
  df_wilcoxon_tidy_final <- df_wilcoxon_AD %>%
    filter(project %in% (df_wilcoxon_AD %>%
                           group_by(project) %>%
                           ungroup() %>% 
                           pull(project)) ) %>%
    mutate(project = str_replace(project, pattern = "_",replacement = " ")) %>%
    group_by(MSR_type) %>%
    mutate(gene_normalised = fct_reorder(gene_normalised, plyr::desc(effsize))) %>%
    dplyr::select(effsize, gene_normalised, MSR_type, q, 
                  Splicing.regulation, Spliceosome, Exon.Junction.Complex, NMD) %>%
    mutate(type = ifelse(Splicing.regulation == 1, "Splicing Regulation",
                         ifelse(Spliceosome == 1, "Spliceosome",
                                ifelse(Exon.Junction.Complex == 1, "Exon Junction Complex", "NMD")))) %>%
    mutate(type = factor(type, levels=c('NMD','Splicing Regulation','Spliceosome', 'Exon Junction Complex'))) %>%
    ungroup()  
  
  

  
  ####################################
  ## PLOT 
  ####################################
  
  age_effsize_plot <- ggplot(data = df_wilcoxon_tidy_final, aes(x = effsize, y = gene_normalised, color = MSR_type)) +
    geom_point(alpha=.7) +
    ggforce::facet_row(~type) +
    geom_vline(mapping = aes(xintercept = 0.053), linetype="dotted") +
    theme_light() +
    ylab("") +
    xlab("Probability of superior MSR after controlling for individual RBP/NMD activity in '60-79'yrs compared to '20-39'yrs") +
    scale_color_manual(values = c("#35B779FF","#64037d"), breaks = c("MSR Donor", "MSR Acceptor"), labels = c("MSR Donor", "MSR Acceptor")) +
    custom_ggtheme + 
    theme( plot.margin = margin(0,0,0,0),
           legend.box.margin = margin(l = -11),
           legend.position = "right", 
           legend.box = "horizontal", 
           axis.text.y = element_text(size = 6, family="Arial", colour = "black"), 
           axis.text.x = element_text(size = 6, family="Arial", colour = "black"),
           axis.title.x = element_text(size = 7, family="Arial", colour = "black"),
           strip.text.x = element_text(size = 7, family="Arial", colour = "black"),
           legend.text = element_text(size = 7, family="Arial", colour = "black") ) +
    scale_size(name = "q:", range=c(5, 1), breaks=c(2.2e-16, 0.5)) +
    guides(color = guide_legend(title = ""))+
    scale_x_continuous(expand = expansion(add = c(0.025, 0.025))) 
  
  ## Save the figure 
  ## ggplot2::ggsave(filename = figure.name, width = 180, height = 240, units = "mm", dpi = 300)
  
  # Set the colors to be used
  fill_colors <- c("#cccccc","#999999","#999999","#999999")
  
  # Find strips glob
  gt<-ggplot_gtable(ggplot_build(age_effsize_plot))
  strips <- which(startsWith(gt$layout$name,'strip'))
  
  # Change the fill color of each strip
  for (s in seq_along(strips)) {
    gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- fill_colors[s]
  }
  
  #dir.create(path = figure.path, recursive = T, showWarnings = F)
  png(filename = file.path(args$figures_folder, "supplementary_figure22.png"), 
      width = 180, height = 240, 
      units = "mm", res = 300)
  plot(gt)
  dev.off()
  
  
  if (plot.stats) {
    
    ## STATS - DONOR
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Donor", q < 0.05) %>%
      arrange(effsize)
    
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Donor", q < 0.05) %>%
      pull(effsize) %>% summary
    
    ## STATS - ACCEPTOR
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Acceptor", q < 0.05) %>%
      arrange(effsize)
    
    df_wilcoxon_tidy_final %>%
      filter(MSR_type == "MSR Acceptor", q < 0.05) %>%
      pull(effsize) %>% summary
  }
}


##################################
## AUXILIARY FUNCTIONS
##################################

## SECTION 0 - GET COMMON INTRONS WITH SIMILAR EXPRESSION LEVELS -----------------------------

get_common_introns <- function() {
  
  all_introns <- map_df(c(args$control_type, args$case_type), function(cluster_id) {
    
    print(paste0(Sys.time(), " - ", args$project_id, " - ", cluster_id))
    
    query <- paste0("SELECT ref_junID, MSR_D, MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                    FROM '", cluster_id, "_", args$project_id, "_nevermisspliced'")
    introns_accurate_splicing <- dbGetQuery(con, query) %>% as_tibble()
    
    
    query <- paste0("SELECT ref_junID, MSR_D,  MSR_A, ref_type, ref_sum_counts, ref_n_individuals 
                    FROM '", cluster_id, "_", args$project_id, "_misspliced'")
    introns_inaccurate_splicing <- dbGetQuery(con, query) %>% as_tibble()
    
    
    rbind(introns_inaccurate_splicing, introns_accurate_splicing) %>% mutate(sample_type = cluster_id) %>% distinct(ref_junID, .keep_all = T)
    
  })
  
  all_introns %>% 
    dplyr::group_by(ref_junID) %>% 
    filter(n() > 1) %>% 
    ungroup  %>%
    group_by(sample_type) %>%
    mutate(mean_coverage = (ref_sum_counts/ref_n_individuals) %>% log10()) %>%
    ungroup() %>% 
    return()
  
  
}

subsample_introns <- function(df_introns) {
  
  set.seed(1000)
  print(paste0(Sys.time(), " - start subsampling ... "))
  
  ## Subsampling introns to control by similarity in mean read coverage
  m.out <- MatchIt::matchit(sample_type ~ mean_coverage,
                            data = df_introns %>% mutate(sample_type = sample_type %>% as.factor()),
                            distance = df_introns$mean_coverage,
                            method = "nearest",
                            caliper = c(mean_coverage = .005),
                            std.caliper = FALSE)
  subsample <- MatchIt::match.data(m.out)
  subsample %>% dplyr::count(sample_type)
  subsample %>% group_by(sample_type) %>% distinct(ref_junID) %>% ungroup() %>% dplyr::count(sample_type)
  
  
  saveRDS(object = subsample, file = file.path(args$results_folder, "common_subsampled_introns_seed1000.rds"))
  
  print(paste0(Sys.time(), " - subsampling finished!"))
}

## SECTION 1 - GENERAL TESTS -----------------------------------------------------------------

index_database <- function() {
  
  
  ## INTRON TABLE -----------------------------------------------------------------
  
  
  ## verify indexes exist on 'intron' master table
  query <- paste0("SELECT * FROM 'sqlite_master' 
                  WHERE tbl_name = 'intron'
                  AND name = 'index_intron_coord'")
  
  if (nrow(dbGetQuery(con, query)) == 0) {
    
    query <- paste0("CREATE UNIQUE INDEX 'index_intron_coord' ON 'intron'(ref_coordinates)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
  }
  
  
  
  ## NOVEL TABLE -----------------------------------------------------------------
  
  
  ## verify indexes exist on 'novel' master table
  query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'novel' AND name = 'index_novel'")
  
  if ( nrow(dbGetQuery(con, query)) == 0 ) {
    
    
    query <- paste0("CREATE UNIQUE INDEX 'index_novel' ON 'novel'(ref_junID,novel_junID)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
    
  }
  
  
  ## verify indexes exist on 'novel' master table
  query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'novel' AND name = 'index_novel_coord'")
  
  if ( nrow(dbGetQuery(con, query)) == 0 ) {
    
    query <- paste0("CREATE UNIQUE INDEX 'index_novel_coord' ON 'novel'(novel_coordinates)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
    
  }
  
  
  ## TRANSCRIPT TABLE -----------------------------------------------------------------
  
  query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'transcript' AND name = 'index_transcript_id'")
  
  if ( nrow(dbGetQuery(con, query)) == 0 ) {
    
    query <- paste0("CREATE UNIQUE INDEX 'index_transcript_id' ON 'transcript'(id)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
    
  }
  
  
  query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'transcript' AND name = 'index_transcript_ensembl_id'")
  
  if ( nrow(dbGetQuery(con, query)) == 0 ) {
    
    query <- paste0("CREATE UNIQUE INDEX 'index_transcript_ensembl_id' ON 'transcript'(transcript_id)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
    
  }
  
  
  ## GENE TABLE -----------------------------------------------------------------
  
  ## verify indexes exist on 'gene' master tableg
  query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'gene' AND name = 'index_gene_id'")
  
  if ( nrow(dbGetQuery(con, query)) == 0 ) {
    
    query <- paste0("CREATE UNIQUE INDEX 'index_gene_id' ON 'gene'(id)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
    
  }
  
  ## verify indexes exist on 'gene' master tableg
  query <- paste0("SELECT * FROM 'sqlite_master' WHERE tbl_name = 'gene' AND name = 'index_gene_ensembl_id'")
  
  if ( nrow(dbGetQuery(con, query)) == 0 ) {
    
    query <- paste0("CREATE UNIQUE INDEX 'index_gene_ensembl_id' ON 'gene'(gene_id)");
    res <- DBI::dbSendQuery(conn = con, statement = query)
    DBI::dbClearResult(res)
    
  }
  
  
}


# "Analysis of MSR measures also demonstrated significantly higher levels of MSRs in AD samples at both donor and acceptor sites 
# (5’ss, effect-size=0.052, one-tailed paired Wilcoxon signed rank test, P<0.001; 3’ss, effect-size=0.054, 
# one-tailed paired Wilcoxon signed rank test, P<0.001)."

#' Title
#' Test for differences in MSR_D and MSR_A between the introns from AD vs control samples considering 
#' the biotype of the transcripts from which the introns originated (protein-coding and noncoding).
#' Using one-tailed paired Wilcoxon test
#' @return
#' @export
#'
#' @examples
plot_MSR_by_biotype <- function() {
  
  common_introns_subsample <- readRDS(file = paste0(results_path, "/common_subsampled_introns_seed1000.rds"))
  common_junID <- common_introns_subsample %>% pull(ref_junID)
  
  
  ################################
  ## TIDY DATA BEFORE PLOTTING
  ################################
  
  df_introns_biotype <- common_introns_subsample %>%
    filter(ref_junID %in% common_junID) %>%
    inner_join(y = master_introns %>% dplyr::select(ref_junID, protein_coding),
               by = "ref_junID") %>%
    filter(protein_coding %in% c(0,100)) %>%
    mutate(biotype = ifelse(protein_coding == 100, "PC", "non PC"))
  
  
  
  df_introns_biotype$biotype = factor(df_introns_biotype$biotype, 
                                      levels = c("non PC","PC"))
  
  df_introns_biotype <- df_introns_biotype %>%
    dplyr::select(ref_junID, biotype, MSR_D, MSR_A) %>%
    gather(key = "MSR_type", value = "MSR", -biotype, -ref_junID)
  
  
  df_introns_biotype$MSR_type = factor(df_introns_biotype$MSR_type, 
                                       levels = c("MSR_A", "MSR_D"))
  
  
  print(paste0(Sys.time(), " - ", df_introns_biotype %>% nrow(), " - introns after tidying!"))
  
  
  df_introns_biotype %>% 
    filter(biotype == "PC") %>%
    pull(MSR) %>%
    summary()
  
  df_introns_biotype %>% 
    filter(biotype == "non PC") %>%
    pull(MSR) %>%
    summary()
  
  
  ###########################
  ## PLOT BY MSR GROUP
  ###########################
  
  # scales::show_col(scales::viridis_pal(option = "B")(40))
  # scales::show_col(scales::viridis_pal(option = "C")(40))
  # scales::show_col(scales::viridis_pal(option = "D")(40))
  # scales::show_col(scales::viridis_pal(option = "E")(40))
  # scales::show_col(scales::viridis_pal(option = "F")(40))
  # scales::show_col(scales::viridis_pal(option = "G")(40))
  # scales::show_col(scales::viridis_pal(option = "H")(40))
  
  
  df_introns_biotype_tidy <- df_introns_biotype %>%
    group_by(biotype, MSR_type) %>%
    distinct(ref_junID, .keep_all = T) %>%
    mutate(percentile_group = case_when(MSR == 0 ~ "0",
                                        MSR > 0 & MSR <= 0.2 ~ "(0,.2]",
                                        MSR > 0.2 & MSR <= 0.4 ~ "(.2,.4]",
                                        MSR > 0.4 & MSR <= 0.6 ~ "(.4,.6]",
                                        MSR > 0.6 & MSR <= 0.8 ~ "(.6,.8]",
                                        MSR > 0.8 & MSR <= 1 ~ "(.8,1]"))
  
  df_introns_biotype_tidy$MSR_type = factor(df_introns_biotype_tidy$MSR_type, 
                                            levels = c("MSR_D","MSR_A"))
  df_introns_biotype_tidy$percentile_group = factor(df_introns_biotype_tidy$percentile_group, 
                                                    levels = c("0","(0,.2]","(.2,.4]","(.4,.6]","(.6,.8]","(.8,1]"))
  
  df_introns_biotype_tidy$biotype = factor(df_introns_biotype_tidy$biotype, 
                                           levels = c("PC","non PC"))
  
  plotMSR_donor <- ggplot(data = df_introns_biotype_tidy ) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = .5, color = "#333333")+
    #ggtitle("MSR Donor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    #ylim(c(0,23000))+
    ggforce::facet_row(~MSR_type) +
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding   ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(title = "",
                               ncol = 2, nrow = 1))
  plotMSR_donor
  
  
  ## Non-protein-coding
  plotMSR_acceptor <- ggplot(data = df_introns_biotype_tidy %>% 
                               filter(MSR_type == "MSR_A")) + 
    geom_bar(aes(x = percentile_group, fill = biotype),
             position = "dodge", linewidth = 0.5, color = "#333333")+
    ggtitle("MSR Acceptor") +
    xlab("Mis-splicing ratio value group") +
    ylab("Number of annotated introns") +
    #ylim(c(0,23000))+
    
    scale_fill_manual(values = c("#333333","#999999"),
                      breaks = c("PC","non PC"),
                      labels = c("Protein-coding ","Non-protein-coding ")) +
    theme_light() +
    custom_ggtheme +
    guides(fill = guide_legend(ncol = 2, nrow = 1))
  
  
  ggpubr::ggarrange(plotMSR_donor,
                    plotMSR_acceptor,
                    nrow = 1,
                    ncol = 2,
                    common.legend = T)
  
  
  
  
  print(paste0(Sys.time(), " - saving plot..."))
  file_name <- paste0(figures_folder, "/panel4bc")
  ggplot2::ggsave(paste0(file_name, ".png"), width = 180, height = 90, units = "mm", dpi = 300)
  
}




##################################
## CALLS
##################################

# AD_control_effsize_MSR_normalised_by_TPM()

