#' Title
#' Removes either all or only the child tables from the database.
#' @param database.path Local path to the .sqlite database
#' @param all Whether to remove all tables or only the child tables.
#'
#' @return
#' @export
#'
#' @examples
sql_remove_tables <- function(database.path,
                              all) {
  
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = database.path)
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  tables <- dbListTables(con)
  tables %>% print()
  
  for (table in tables) {
    #if (str_detect(table %>% tolower(), pattern = "kidney")) {
    
    if (!all) {
      
      if ( !(table %in% c("gene", "transcript", "intron", "metadata", "novel")) ) {
        dbRemoveTable(conn = con, table)
        print(paste0("Table: '", table, "' has been removed!"))
      }
      
    } else {
      
      dbRemoveTable(conn = con, table)
      print(paste0("Table: '", table, "' has been removed!"))
      
    }
  }
  
  DBI::dbDisconnect(conn = con)
}
#dbRemoveTable(conn = con, "PD_SRP058181_misspliced")