###################################################################
##
## 16S analysis
##
###################################################################
seu.obj.umap <- readRDS(paste0("~/PRIME-TR/test/DSP_021_seuratObject_umap.rds"))
semisup <- readRDS(paste0("~/PRIME-TR/test/DSP_021_InSituType_semisupervised_clustering.rds"))
seu.obj.umap[["InSituType_semisup"]] <- semisup$clust

# On SQLite databases in R:
# https://www.datacamp.com/tutorial/sqlite-in-r
# Why we probably shouldn't be using SQLite databases for the file sizes that we deal with:
# https://softwareengineering.stackexchange.com/questions/332069/what-is-a-realistic-real-world-maximum-size-for-a-sqlite-database
# https://lowendtalk.com/discussion/166351/sqlite3-vs-mariadb-which-to-use

# https://duckdb.org/docs/api/r.html
# Start a connection with a DB file.
con <- dbConnect(RSQLite::SQLite(), "/rsrch6/scratch/genomic_med/lwfong/DSP_021_tx.db", read_only = FALSE)
# con_new <- dbConnect(RSQLite::SQLite(), "~/PRIME-TR/test/DSP_021_tx.db", read_only = FALSE)
# We are creating the database in /tmp/ because for some reason, we're running out of disk space when we do it in my personal directory ... 
# dbDisconnect(con)
# Read in the CSV.
# https://www.michaelc-m.com/manual_posts/2022-01-27-big-CSV-SQL.html
# Transcript file colnames and column types.
# [1] "fov"         "cell_ID"     "cell"        "x_local_px"  "y_local_px"  "x_global_px" "y_global_px" "z"           "target"      "CellComp"   
# character character character numeric numeric numeric numeric numeric character character
coltypes <- "cccnnnnccc"
f <- file("~/PRIME-TR/test/flatFiles/TR542AcralTMAL/TR542AcralTMAL_tx_file.csv")
# system.time({
read_csv_chunked(file = f, 
                 callback = function(chunk, dummy) {
                   dbWriteTable(conn = con, name = "transcripts", chunk, append = T)}, 
                 chunk_size = 10000, col_types = coltypes)

# })
# chunk_size = 10000
# user   system  elapsed 
# 1053.902   61.296 1109.351 

# Create necessary columns to do calculations later.
# Note that SQLite apparently does not support certain window functions
# Error in `median()`:
#   ! Window function `median()` is not supported by this database.
# so we will have to work around this somehow.
# See https://github.com/tidyverse/dplyr/issues/2960 for suggestions.
# Also:
# Error in `summarize()` at common/scripts/helper_functions.R:33:3:
#   ℹ In argument: `median_x = median(x_um)`
# Caused by error:
#   ! In dbplyr you cannot use a variable created in the same `summarise()`.
# ✖ `x_um` was created earlier in this `summarise()`.
# ℹ You need an extra `mutate()` step to use it.
# !!! See this on how we're going to write a new table after mutating everything:
# https://stackoverflow.com/questions/62467655/r-and-dplyr-how-can-i-use-compute-to-create-a-persistent-table-from-sql-query
# Custom function to write to database, adapted from https://stackoverflow.com/a/62601688/23532435.
write_to_database <- function(input_tbl, db_connection, tbl_name) {
  # Drop table if it exists.
  DBI::dbExecute(con, glue::glue("DROP TABLE IF EXISTS {tbl_name};"))
  
  # SQL query to create new table.
  sql_query <- glue::glue("CREATE TABLE {tbl_name} AS (\n",
                          "SELECT * FROM (\n",
                          dbplyr::sql_render(input_tbl),
                          "\n)) WITH DATA;")
  # Run query.
  DBI::dbExecute(db_connection, as.character(sql_query))
}
# walltime <- system.time({
#   mem_use <- 
#     profmem({
dplyr::tbl(con, "transcripts") |>
  dplyr::mutate(x_um = x_local_px * 0.12,
                y_um = y_local_px * 0.12) |>
  add_group_summaries(c("cell"),
                      median_x = median(x_um),
                      median_y = median(y_um),
                      median_z = median(z),
                      max_x = max(x_um),
                      max_y = max(y_um),
                      max_z = max(z),
                      min_x = min(x_um),
                      min_y = min(y_um),
                      min_z = min(z)
  ) |>
  # 106.9076 Gb (with collect())
  dplyr::mutate(min_80_x = median_x - (median_x - min_x) * 0.8,
                max_80_x = median_x + (max_x - median_x) * 0.8,
                min_80_y = median_y - (median_y - min_y) * 0.8,
                max_80_y = median_y + (max_y - median_y) * 0.8,
                min_80_z = median_z - (median_z - min_z) * 0.8,
                max_80_z = median_z + (max_z - median_z) * 0.8) |>
  dplyr::slice_min(n = 5, order_by = cell_ID) |>
  collect()
# collect() no longer used because it collects the tibble into memory. 
# write as new table to connection.
show_query() |>
  write_to_database(con, "transcripts_16S")
#     })
# })
# user   system  elapsed 
# 1800.697   79.299 1885.201 
# 1885.201 / 60
dbDisconnect(con)
# dbDisconnect(con_out)

# https://stackoverflow.com/questions/1727772/quickly-reading-very-large-tables-as-dataframes/1820610#1820610
system.time(bigdf <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F)))

system.time({
  transcripts_L <- read_csv("~/PRIME-TR/test/flatFiles/TR542AcralTMAL/TR542AcralTMAL_tx_file.csv") %>%
    #   user  system elapsed 
    # 317.882  37.081 154.023 
    # filter(cell %in% Acral$cell_ID) %>%
    group_by(cell) %>%
    #   user  system elapsed 
    # 5.649   0.944   6.649 
    mutate(
      x_um = x_local_px * 0.12,
      y_um = y_local_px * 0.12,
      median_x = median(x_um), 
      median_y = median(y_um),
      median_z = median(z),
      max_x = max(x_um), 
      max_y = max(y_um),
      max_z = max(z),
      min_x = min(x_um), 
      min_y = min(y_um),
      min_z = min(z),
      min_80_x = median_x - (median_x - min_x) * 0.8,
      max_80_x = median_x + (max_x - median_x) * 0.8,
      min_80_y = median_y - (median_y - min_y) * 0.8,
      max_80_y = median_y + (max_y - median_y) * 0.8,
      min_80_z = median_z - (median_z - min_z) * 0.8,
      max_80_z = median_z + (max_z - median_z) * 0.8,
    )
  # user  system elapsed 
  # 89.451  29.248 119.027 
})
transcripts_R <- read_csv("~/PRIME-TR/test/flatFiles/TR542AcralTMAR/TR542AcralTMAR_tx_file.csv") %>%
  #   user  system elapsed 
  # 317.882  37.081 154.023 
  # filter(cell %in% Acral$cell_ID) %>%
  group_by(cell) %>%
  #   user  system elapsed 
  # 5.649   0.944   6.649 
  mutate(
    x_um = x_local_px * 0.12,
    y_um = y_local_px * 0.12,
    median_x = median(x_um), 
    median_y = median(y_um),
    median_z = median(z),
    max_x = max(x_um), 
    max_y = max(y_um),
    max_z = max(z),
    min_x = min(x_um), 
    min_y = min(y_um),
    min_z = min(z),
    min_80_x = median_x - (median_x - min_x) * 0.8,
    max_80_x = median_x + (max_x - median_x) * 0.8,
    min_80_y = median_y - (median_y - min_y) * 0.8,
    max_80_y = median_y + (max_y - median_y) * 0.8,
    min_80_z = median_z - (median_z - min_z) * 0.8,
    max_80_z = median_z + (max_z - median_z) * 0.8,
  )
transcripts_all <- rbind(transcripts_L, transcripts_R)
rm(transcripts_L, transcripts_R)
gc()

#----------------------------------------------------------------------------------------------------
# Annotate confidence levels
#----------------------------------------------------------------------------------------------------
system.time({
  Transcripts_80_16S <- transcripts_L %>%
    filter(target == "Microbial16S") %>%
    mutate(
      x_confidence = ifelse(min_80_x < x_um & x_um < max_80_x, "High", "Low"),
      y_confidence = ifelse(min_80_y < y_um & y_um < max_80_y, "High", "Low"),
      z_confidence = ifelse(min_80_z < z & z < max_80_z, "High", "Low"),
      High_confidence_80_16S = ifelse(x_confidence == "High" & y_confidence == "High" &
                                        z_confidence == "High", "High", "Low")
    ) %>%
    filter(High_confidence_80_16S == "High")
})
# user  system elapsed 
# 6.893   0.547   7.457 
seu.obj$Intracellular_80 <- ifelse(seu.obj$cell_ID %in% Transcripts_80_16S$cell, "Intracellular", "Extracellular")

#----------------------------------------------------------------------------------------------------
# Annotate 16S expression
#----------------------------------------------------------------------------------------------------
system.time({
  Acral_16S <- Acral@assays$RNA$data["Microbial16S", ] %>%
    as.data.frame() %>%
    rename("Microbial16S" = ".") %>%
    mutate(
      Intracellular_80 = Acral$Intracellular_80,
      cell_ID = Acral$cell_ID
    ) %>%
    mutate(`16S_expression` = case_when(
      (Microbial16S <= quantile(Acral_16S$Microbial16S, 0.960)) ~ "Negative",
      (Microbial16S > quantile(Microbial16S, 0.960) & Microbial16S <= quantile(
        Microbial16S,
        0.990
      )) ~ "Low",
      (Microbial16S > quantile(Microbial16S, 0.990) & Microbial16S <= quantile(
        Microbial16S,
        0.999
      )) ~ "Medium",
      (Microbial16S > quantile(Microbial16S, 0.999)) ~ "High"
    )) %>%
    # Where is `Intracellular80`?
    mutate(Intracellular_80_16S = case_when(
      Acral_16S$Intracellular_80 == "Intracellular" &
        Acral_16S$`16S_expression` != "Negative" ~ "Intracellular_16S",
      Acral_16S$Intracellular_80 == "Extracellular" &
        Acral_16S$`16S_expression` != "Negative" ~ "Extracellular_16S",
      Acral_16S$Intracellular_80 == "Extracellular" &
        Acral_16S$`16S_expression` == "Negative" | Acral_16S$Intracellular_80 == "Intracellular" &
        Acral_16S$`16S_expression` == "Negative" ~ "Negative"
    )) %>%
    view()
})