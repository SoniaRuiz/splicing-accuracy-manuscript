#' Title
#' Create metadata table
#' @param database.path Local path to the .sqlite database
#' @param recount3.project.IDs List of recount3 projects 
#' @param results.folder Path to the local folder where the results to read from are stored
#'
#' @return
#' @export
#'
#' @examples
sql_create_master_table_metadata <- function(database.path,
                                             recount3.project.IDs,
                                             results.folder)  {
  
  message(Sys.time(), " - creating metadata table ... ")
  
  df_metadata <- map_df(recount3.project.IDs, function(project_id) { 
    
    message(Sys.time(), " - getting info from ", project_id)
    # project_id <- recount3.project.IDs[1]
    # project_id <- recount3.project.IDs[2]
    # project_id <- recount3.project.IDs[6]
    
    if ( str_detect(database.path, pattern = "age") ) {
      metadata_file <- paste0(results.folder, "/", project_id, 
                              "/base_data/", project_id, "_age_samples_metadata.rds")
    } else {
      metadata_file <- paste0(results.folder, "/", project_id, 
                              "/base_data/", project_id, "_samples_metadata.rds")  
    }
    
    
    if ( file.exists(metadata_file) ) {
      
      metadata <- readRDS(file = metadata_file)
      
      return(metadata %>%
               mutate(SRA_project = project_id))
    } else {
      
      return(NULL)
    }
    
  })
  
  
  df_metadata %>% as_tibble() %>% dplyr::count(cluster)
  
  con <- dbConnect(RSQLite::SQLite(), database.path)
  DBI::dbWriteTable(conn = con,
                    name = "metadata",
                    value = df_metadata %>% as_tibble(),
                    overwrite = T)
  
  DBI::dbDisconnect(conn = con)
  
  print(paste0("Table: 'metadata' created!"))
}
