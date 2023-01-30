library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)

#########################
## FOR BRAIN HIPPOCAMPUS
#########################

project_id <- "BRAIN"
cluster_id <- "Brain - Hippocampus"

## 1. Connect to the database

gtf_version <- 105
main_project <- "splicing_unique"
database_path <- paste0("~/PROJECTS/splicing-project-recount3/database/v", gtf_version, "/", main_project,
                        "/", main_project, ".sqlite")

con <- dbConnect(RSQLite::SQLite(), database_path)
dbListTables(con)



######################
## Reference introns
######################

query <- paste0("SELECT ref_junID, ref_n_individuals, ref_sum_counts FROM '", cluster_id, "_", project_id, "_misspliced'")
df_hipp_introns <- dbGetQuery(con, query) %>% as_tibble()
query <- paste0("SELECT ref_junID, ref_coordinates FROM 'intron' WHERE ref_junID IN (", paste(df_hipp_introns$ref_junID, collapse = ","),")")
df_hipp_introns <- df_hipp_introns %>%
  inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
             by = "ref_junID") %>% 
  as_tibble() 


#####################
## Novel junctions
#####################

query <- paste0("SELECT novel_junID, novel_n_individuals, novel_sum_counts FROM '", cluster_id, "_", project_id, "_misspliced'")
df_hipp_novel <- dbGetQuery(con, query) %>% as_tibble()
query <- paste0("SELECT novel_junID, novel_coordinates FROM 'novel' WHERE ref_junID IN (", paste(df_hipp_novel$novel_junID, collapse = ","),")")
df_hipp_novel <- df_hipp_novel %>%
  inner_join(y = dbGetQuery(con, query) %>% as_tibble(),
             by = "novel_junID") %>% 
  as_tibble() 



###############################
## 2. Load the new UNIQUE counts
###############################


hipp_unique_sampes <- readRDS(paste0("~/PROJECTS/splicing-project-recount3/results/", project_id, "/v", gtf_version, "/splicing/",
                                     "/base_data/unique/", project_id, "_", cluster_id, "_samples_used.rds"))
hipp_unique_counts <- readRDS(paste0("~/PROJECTS/splicing-project-recount3/results/", project_id, "/v", gtf_version, "/splicing/",
                                     "/base_data/unique/", project_id, "_", cluster_id, "_split_read_counts.rds"))

hipp_unique_coverage <- get_mean_coverage(split_read_counts = hipp_unique_counts,
                                          samples = hipp_unique_sampes,
                                          junID = c(df_hipp_novel$novel_coordinates,
                                                    df_hipp_introns$ref_coordinates) %>% unique())


###############################
## UPDATE THE DATABASE
## 3. Update the sum_counts 
## and sum_n_individuals columns
###############################


## 1. Update the introns

hipp_introns_to_update <- df_hipp_introns %>%
  distinct(ref_junID, .keep_all = T) %>%
  inner_join(y = hipp_unique_coverage,
             by = c("ref_coordinates" = "junID"))

for (junID in hipp_introns_to_update$ref_junID) {
  
  # junID <- 71818
  record_to_update <- hipp_introns_to_update %>%
    filter(ref_junID == junID)
  
  query <- paste0("UPDATE '", cluster_id, "_", project_id, "_misspliced'
                  SET ref_n_individuals = ", record_to_update$n_individuals,
                  ", ref_sum_counts = ", record_to_update$sum_counts,
                  " WHERE ref_junID = ", record_to_update$ref_junID, ";")
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  print(paste0(junID))
}


## 2. Update the novel Junctions

hipp_novel_to_update <- df_hipp_novel %>%
  distinct(novel_junID, .keep_all = T) %>%
  inner_join(y = hipp_unique_coverage,
             by = c("novel_coordinates" = "junID"))

for (junID in hipp_novel_to_update$ref_junID) {
  
  # junID <- 71818
  record_to_update <- hipp_novel_to_update %>%
    filter(novel_junID == junID)
  
  query <- paste0("UPDATE '", cluster_id, "_", project_id, "_misspliced'
                  SET novel_n_individuals = ", record_to_update$n_individuals,
                  ", novel_sum_counts = ", record_to_update$sum_counts,
                  " WHERE novel_junID = ", record_to_update$novel_junID, ";")
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  print(paste0(junID))
}




query <- paste0("SELECT ref_junID, ref_n_individuals, ref_sum_counts FROM '", cluster_id, "_", project_id, "_misspliced' 
                WHERE ref_junID = 71818" )
dbGetQuery(con, query) %>% as_tibble()


query <- paste0("UPDATE '", cluster_id, "_", project_id, "_misspliced'
SET ref_n_individuals = ", paste(hipp_introns_to_update$ref_n_individuals, collapse = ","),
#"column_2 = ", paste(df_hipp_novel$novel_junID, collapse = ","),
"WHERE ref_junID IN (", paste(hipp_introns_to_update$ref_junID, collapse = ","), ")");


res <- DBI::dbSendQuery(conn = con, statement = query)
DBI::dbClearResult(res)

## 4. Calculate new TPM