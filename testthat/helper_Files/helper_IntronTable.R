## Load the necessary tables from the database
con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
df_novel <- dplyr::tbl(con, "novel") %>% dplyr::collect()
df_intron <- dplyr::tbl(con, "intron") %>% dplyr::collect()
df_mane <- dplyr::tbl(con, "mane") %>% dplyr::collect()
df_gene <- dplyr::tbl(con, "gene") %>% dplyr::collect()