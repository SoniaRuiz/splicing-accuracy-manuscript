library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(DBI)
library(dplyr)
library(here)
library(doParallel)
library(Biostrings)
library(protr)

# source("/home/sruiz/PROJECTS/splicing-accuracy-manuscript/init_age.R")



base_folder <- here::here()
dependencies_folder <- paste0(base_folder, "/dependencies/")

#####################################
## LOAD SOURCE SCRIPTS
#####################################


setwd(file.path(base_folder,"scripts"))
files.sources = list.files()
sapply(files.sources, source)
setwd(file.path(base_folder))


#####################################
## SET MAIN VARIABLES
#####################################


gtf_version <- 110

supportive_reads <- 2
replace <- T

data_source <- "data_sources/gtex"


project_name <- paste0("splicing_",supportive_reads,"read")
database_name <- paste0(project_name, "_age")


levelqc1_folder <- paste0(base_folder, "/database/")
database_folder <- paste0(base_folder, "/database/", database_name, "/", gtf_version, "/")
results_folder <- paste0(base_folder, "/results/", project_name, "/", gtf_version, "/")
tpm_folder <- paste0(base_folder, "/results/tpm/")

age_projects <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))


##################################
## CALLS - SQL GENERATION
##################################


# if ( !exists("project_init") ) {
#   project_init <- age_stratification_init_data(projects.id = age_projects,
#                                                gtf.version = gtf_version,
#                                                project.name = project_name,
#                                                results.folder = results_folder) %>%
#     as_tibble()
# }

# age_projects <- project_init$project %>% unique()


# for (project_id in age_projects) {
#   
# 
#   # project_id <- age_projects[2]
#   print(paste0(Sys.time(), " --> ", project_id))
#   project_init_local <- project_init %>% filter(project == project_id)
#   project_results_folder <- file.path(results_folder, project_id)
#   
#   message(project_id)
#   
#   age_stratification_annotate(age.groups = project_init_local$age_group %>% unique(),
#                               project.id = project_id,
#                               gtf.version = gtf_version,
#                               project.name = project_name,
#                               age.samples.clusters = project_init_local,
#                               results.folder = project_results_folder)
# 
#   age_stratification_junction_pairing(age.groups = project_init_local$age_group %>% unique(),
#                                       project.id = project_id,
#                                       gtf.version = gtf_version,
#                                       project.name = project_name,
#                                       age.samples.clusters = project_init_local,
#                                       results.folder = project_results_folder,
#                                       replace)
# }
# 
# 
# get_all_annotated_split_reads(recount3.project.IDs = age_projects,
#                               all.clusters = project_init$age_group %>% unique(),
#                               database.folder = database_folder,
#                               results.folder = results_folder)
# 
# get_all_raw_jxn_pairings(recount3.project.IDs = age_projects,
#                          all.clusters = project_init$age_group %>% unique(),
#                          database.folder = database_folder,
#                          results.folder = results_folder)

age_projects <- readRDS(file = paste0(results_folder, "/all_final_projects_used.rds"))


# generate_recount3_tpm(recount3.project.IDs = age_projects,
#                       data.source = data_source,
#                       tpm.folder = tpm_folder,
#                       results.folder = results_folder)
# 
# 
# 
# tidy_data_pior_sql(recount3.project.IDs = age_projects,
#                    levelqc1.folder = levelqc1_folder,
#                    all.clusters = project_init$age_group %>% unique(),
#                    database.folder = database_folder,
#                    results.folder = results_folder)


database_path <- paste0(database_folder,  "/", database_name, ".sqlite")

sql_database_generation(database.path = database_path,
                        recount3.project.IDs = age_projects,
                        remove.all = T,
                        database.folder = database_folder,
                        results.folder = results_folder,
                        gtf.version = gtf_version)





####################################
## Visualise sample pairing metadata
####################################

# ggplot(project_init %>% 
#          dplyr::count(project,age_group)%>%
#          arrange(age_group , n) %>%
#          mutate(project = fct_inorder(project)))+
#   geom_bar(aes(y = project, x = n, fill = age_group),
#            stat = "identity", position = position_dodge()) +
#   custom_theme
# 
# ggplot(project_init %>% 
#          mutate(sex = sex %>% as.character()) %>%
#          dplyr::count(project, sex) %>%
#          arrange(sex , n) %>%
#          mutate(project = fct_inorder(project)) ) +
#   geom_bar(aes(y = project, x = n, group = sex, fill = sex),
#            stat = "identity", position = "dodge") + 
#   theme_light()+
#   custom_theme
# 
# 
# ggplot(project_init ) +
#   geom_density(aes(x = rin, fill = age_group), alpha = 0.5) + 
#   theme_light() +
#   labs(x = "RIN" ) +
#   scale_fill_hue() +
#   guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
#   custom_theme
# 
# 
# ggplot(project_init ) +
#   geom_density(aes(x = mapped_read_count, fill = age_group), alpha = 0.5) + 
#   theme_light() +
#   scale_fill_hue() +
#   labs(x = "Mapped Read Count" ) +
#   guides(fill = guide_legend(title = "Age group: ", ncol = 3, nrow = 1)) +
#   custom_theme
  