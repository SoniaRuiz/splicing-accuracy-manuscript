gtf_version <- 105
main_project <- "splicing"
folder_database <- paste0(getwd(), "/database/v", gtf_version, "/")


#################################
## GET ALL SPLIT READS FROM RECOUNT3
#################################

print("Loading tidy split reads from recount3....")
all_split_reads_details_w_symbol_reduced_keep_gr <- 
  readRDS(file = paste0(folder_database, "/all_split_reads_gtex_recount3_", gtf_version, "_tidy.rds"))


#################################
## GET RSE DATA FOR THYROID TISSUE
#################################

## Get data from Thyroid tissue
project_id <- "THYROID"

folder_results <- paste0(getwd(), "/results/", project_id, "/v", gtf_version, "/", 
                         main_project, "/base_data/unique_mapping/")
dir.create(file.path(folder_results), recursive = TRUE, showWarnings = T)

print(paste0(Sys.time(), " - getting junction data from recount3 - '", project_id, "' tissue..."))

rse <- recount3::create_rse_manual(
  project = project_id,
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  jxn_format = "UNIQUE",
  type = "jxn"
)


#################################
## GET METADATA
#################################

metadata.info <- rse %>% 
  SummarizedExperiment::colData() %>%
  as_tibble() %>%
  filter(gtex.smafrze != "EXCLUDE") %>%
  filter(gtex.smrin >= 6.0)

saveRDS(object = metadata.info, 
        file = paste0(folder_results, "/", project_id, "_samples_metadata.rds"))


#################################
## GET SPLIT READS AND COUNTS  
## PER TISSUE
#################################

clusters_ID <- rse$gtex.smtsd %>% unique()
cluster_id <- clusters_ID[1]

minimum_samples <- 70
cluster_samples <- metadata.info %>% 
  filter(gtex.smtsd == cluster_id) %>%
  distinct(external_id) %>% 
  pull()

saveRDS(object = cluster_samples, 
        file = paste0(folder_results, "/", project_id, "_", cluster_id, "_samples_used.rds"))


################
## Get counts
################

print(paste0(Sys.time(), " - getting split read counts..."))

counts <- (rse[, rse$external_id %in% cluster_samples] %>% SummarizedExperiment::assays())[[1]]
counts <- counts[(rownames(counts) %in%
                    all_split_reads_details_w_symbol_reduced_keep_gr$junID),]
counts <- counts[rowSums(counts) > 0, ]
counts <- counts %>% as.matrix()
counts <- counts %>% as_tibble(rownames = "junID")

# counts %>% head
print(object.size(counts), units = "GB")
saveRDS(object = counts,
        file = paste0(folder_results, "/", project_id, "_", cluster_id, "_split_read_counts.rds"))
gc()

#######################
## Get split reads ID
#######################

print(paste0(Sys.time(), " - getting split read IDs..."))
all_split_reads <- counts %>%
  dplyr::select(junID) %>%
  data.table::as.data.table() %>%
  inner_join(y = all_split_reads_details_w_symbol_reduced_keep_gr,
             by = "junID")

#######################
## QC and save results
#######################

if (any(all_split_reads$width < 25) |
    any(str_detect(string = all_split_reads$chr, pattern = "random")) |
    any(str_detect(string = str_to_lower(all_split_reads$chr), pattern = "u")) |
    any(!(all_split_reads$type %in% c("annotated", "novel_acceptor", "novel_donor")))) {
  print("ERROR! The split reads do not meet the minimum filtering criteria.")
  break;
}
# if (setdiff(counts$i, all_split_reads$index) %>% length() > 0) {
#   print("ERROR! The split reads are not identical.")
#   break;
# }

## Check how much memory this object uses
print(object.size(all_split_reads %>% data.table::as.data.table()), units = "GB")

## Save the object
saveRDS(object = all_split_reads %>% data.table::as.data.table(),
        file = paste0(folder_results, "/", project_id, "_", cluster_id, "_all_split_reads.rds"))


#################################
## LOAD DATA FROM THE DATABASE
#################################

database_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", main_project,
                        "/", main_project, ".sqlite")

con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)


query <- paste0("SELECT novel_junID, novel_sum_counts, novel_n_individuals FROM '", cluster_id, "_", project_id, "_misspliced'")
df_misspliced_novel <- dbGetQuery(con, query) %>% as_tibble()
query <- paste0("SELECT novel_junID, novel_coordinates FROM 'novel' WHERE novel_junID IN (", paste(df_misspliced_novel$novel_junID, collapse = ","),")")
df_misspliced_novel_database <- df_misspliced_novel %>%
  inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
             by = "novel_junID") %>% 
  as_tibble() 



query <- paste0("SELECT ref_junID, ref_sum_counts, ref_n_individuals FROM '", cluster_id, "_", project_id, "_misspliced'")
df_misspliced_intron <- dbGetQuery(con, query) %>% as_tibble()
query <- paste0("SELECT ref_junID, ref_coordinates FROM 'intron' WHERE ref_junID IN (", paste(df_misspliced_intron$ref_junID, collapse = ","),")")
df_misspliced_intron_database <- df_misspliced_intron %>%
  inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
             by = "ref_junID") %>% 
  as_tibble() 

#################################
## CHECK IF THE NUMBER OF READS 
## IS THE SAME
#################################

unique_junID_coverage <- get_mean_coverage(split_read_counts = counts,
                                           samples = cluster_samples,
                                           junID = df_misspliced_intron_database$ref_coordinates %>% unique())

gc()
unique_junID_coverage %>% head()


df_misspliced_intron_database %>%
  inner_join(y = unique_junID_coverage,
             by = c("ref_coordinates" = "junID"))


df_diff_counts_unique_mapping <- df_misspliced_intron_database %>%
  inner_join(y = unique_junID_coverage,
             by = c("ref_coordinates" = "junID")) %>% 
  mutate(diff_counts_unique_mapping = ref_sum_counts - sum_counts) 

df_diff_counts_unique_mapping$diff_counts_unique_mapping %>% summary()
df_diff_counts_unique_mapping$diff_counts_unique_mapping %>% get_mode()
df_diff_counts_unique_mapping$diff_counts_unique_mapping %>% table()

     