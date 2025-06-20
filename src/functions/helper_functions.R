# # Setup
# if(Sys.info()['sysname']=="Linux" & base::grep("dragon", Sys.info()['nodename'])) .libPaths("/home/lwfong/R/ubuntu/4.3.1")
# source("../LFI-009_bisantine-cosmo/src/ligand-receptor_interactions.R")

# Functions ----
# https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

# https://stackoverflow.com/a/1395280
getSizeAllObjects <- function() { return(sort(sapply(ls(), function(x){object.size(get(x))}), decreasing = TRUE)) }

#' Export RDS file containing name of latest completed module:
#' 
#' @param output_path Directory to which Rds file should be exported
#' @param module The name of the current module that was just completed
#' @return A Boolean
#' @export
updateLatestModule <- function(output_path, module) {
  file_name <- "latest_rule.Rds"
  latest_rule <- module
  saveRDS(latest_rule, paste0(output_path, file_name))
  
  return(TRUE)
}

matrix_to_named_row_list <- function(mat) {
  # Check that input is a matrix
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }
  
  # Check that rownames are present
  if (is.null(rownames(mat))) {
    stop("Matrix must have row names.")
  }
  
  # Convert to list of row vectors
  row_list <- split(mat, row(mat))
  
  # Set names using rownames
  names(row_list) <- rownames(mat)
  
  return(row_list)
}

# QC ----
calcSignalStrength <- function(n_features, n_features_neg, seu_obj) { # reseg = F
  # nCount_RNA <- ifelse(reseg, seu.obj.reseg@meta.data$nCount_RNA, seu.obj@meta.data$nCount_RNA)
  # nCount_negprobes <- ifelse(reseg, seu.obj.reseg@meta.data$nCount_negprobes, seu.obj@meta.data$nCount_negprobes)
  nCount_RNA <- seu_obj@meta.data$nCount_RNA
  nCount_negprobes <- seu_obj@meta.data$nCount_negprobes
  
  # For cell x, calculate
  # 1) the count:gene ratio.
  gene_count_rat <- nCount_RNA / n_features
  
  # 2) the count:(negative probe) ratio.
  # Subset to include only the negative probes. 
  neg_count_rat <- nCount_negprobes / n_features_neg
  
  # Now calculate the (count:gene):(count:negative probe) ratio.
  return(gene_count_rat / neg_count_rat)
}

# Function to perform t-test for each row (target QC).
gene.t.test <- function(row, control_vector) {
  t_result <- t.test(row, control_vector, alternative = "less")
  return(t_result$p.value)
}

#' Calculate QC stats
#' 
#' Calculate QC stats and create summary table.
#' @param data A data.frame containing the cell metadata. 
#' @param total_counts_col The name of the column in the `data` data.frame containing the total counts for a given cell.
#' @param area_col The name of the column in the `data` data.frame containing the area for a given cell.
#' @param fov_col The name of the column in the `data` data.frame containing the FOV identifier. 
#' @param min_counts_per_cell Minimum number of counts a cell must have to pass QC.
#' @param min_features_per_cell Minimum number of unique features a cell must have to pass QC.
#' @param proportion_neg_counts MAXIMUM proportion of negative probes a cell must have to pass QC.
#' @param count_dist Minimum ratio of counts:unique genes (also called "complexity") a cell must have to pass QC.
#' @param max_area Maximum area a cell must have to pass QC.
#' @param min_signal_strength Minimum signal strength (ratio of [counts:genes]:[counts:negative probes]) a cell must have to pass QC.
#' @param min_cells_per_fov Minimum number of cells an FOV must have to pass QC.
#' @export
calculate_qc_stats <- function(data, total_cells_post_qc, total_counts_col, total_features_col, area_col, fov_col,
                               min_counts_per_cell = 200, min_features_per_cell = 200, max_proportion_neg_counts = 0.1, min_count_dist = 1, max_area = 30000, min_signal_strength = 4,
                               min_cells_per_fov = 50) {
  # Cell QC stats
  ## Total cells
  total_cells <- nrow(data)
  ## Total counts
  avg_total_counts <- mean(data[[total_counts_col]], na.rm = TRUE)
  low_counts <- sum(data[[total_counts_col]] < min_counts_per_cell, na.rm = TRUE)
  ## Total features
  avg_total_features <- mean(data[[total_features_col]], na.rm = TRUE)
  low_features <- sum(data[[total_features_col]] < min_features_per_cell, na.rm = TRUE)
  ## Proportion negative counts
  avg_prop_neg_counts <- mean(data[["ProportionNegative"]], na.rm = TRUE)
  high_neg_counts <- sum(data[["ProportionNegative"]] > max_proportion_neg_counts, na.rm = TRUE)
  ## Count distribution (complexity)
  avg_count_dist <- mean(data[["Complexity"]], na.rm = TRUE)
  low_count_dist <- sum(data[["Complexity"]] < min_count_dist, na.rm = TRUE)
  ## Area
  avg_area <- mean(data[[area_col]], na.rm = TRUE)
  high_area <- sum(data[[area_col]] > max_area, na.rm = TRUE)
  ## Signal strength
  avg_signal_strength <- data[["SignalStrength"]] %>% .[is.finite(.)] %>% mean(na.rm = TRUE)
  low_signal_strength <- sum(data[["SignalStrength"]] < min_signal_strength, na.rm = TRUE)
  
  # FOV QC stats
  cells_per_fov <- table(data[[fov_col]])
  avg_cells_per_fov <- mean(cells_per_fov)
  fovs_low_cells <- sum(cells_per_fov < min_cells_per_fov)
  total_fovs <- length(cells_per_fov)
  
  # Create data frame for table
  stats_df <- data.frame(
    Metric = c(
      # Total cells
      "Total number of cells before QC",
      "Total number of cells after QC",

      # Total counts
      "Average total counts",
      paste0("Cells with total counts < ", min_counts_per_cell),
      "Percent cells with low counts",

      # Total features
      "Average total features",
      paste0("Cells with total features < ", min_features_per_cell),
      "Percent cells with low features",

      # Proportion negative counts
      "Average proportion negative counts",
      paste0("Cells with proportion of negative counts > ", max_proportion_neg_counts),
      "Percent cells with high proportion of negative counts",

      # Count distribution (complexity)
      "Average count distribution (complexity)",
      paste0("Cells with count distribution (complexity) < ", min_count_dist),
      "Percent cells with low count distribution",

      # Area
      "Average area",
      paste0("Cells with area > ", max_area),
      "Percent cells with high area",

      # Signal strength
      "Average signal strength",
      paste0("Cells with signal strength < ", min_signal_strength),
      "Percent cells with low signal strength",

      # FOVs
      "Total number of FOVs",
      "Average cells per FOV",
      paste0("FOVs with < ", min_cells_per_fov, " Cells"),
      "Percent FOVs with low cell count"
    ),
    Value = c(
      # Total cells
      as.integer(total_cells),
      as.integer(total_cells_post_qc),
      
      # Total counts
      round(avg_total_counts, 2),
      as.integer(low_counts),
      round(100 * low_counts / total_cells, 2),
      
      # Total features
      round(avg_total_features, 2),
      as.integer(low_features),
      round(100 * low_features / total_cells, 2),
      
      # Proportion negative counts
      round(avg_prop_neg_counts, 2),
      as.integer(high_neg_counts),
      round(100 * high_neg_counts / total_cells, 2),
      
      # Count distribution (complexity)
      round(avg_count_dist, 2),
      as.integer(low_count_dist),
      round(100 * low_count_dist / total_cells, 2),
      
      # Area
      round(avg_area, 2),
      as.integer(high_area),
      round(100 * high_area / total_cells, 2),
      
      # Signal strength
      round(avg_signal_strength, 2),
      as.integer(low_signal_strength),
      round(100 * low_signal_strength / total_cells, 2),
      
      # FOVs
      as.integer(total_fovs),
      round(avg_cells_per_fov, 2),
      as.integer(fovs_low_cells),
      round(100 * fovs_low_cells / total_fovs, 2)
    )
  )
  
  return(stats_df)
}

# # Function to create and save QC stats table as PDF
# create_qc_table <- function(data, total_counts_col, area_col, fov_col, output_file) {
#   # Calculate stats
#   stats_df <- calculate_qc_stats(data, total_counts_col, area_col, fov_col)
#   
#   # Create table using gt
#   table <- gt::gt(stats_df) %>%
#     gt::tab_header(
#       title = "Quality Control Statistics"
#     ) %>%
#     gt::fmt_number(
#       columns = "Value",
#       decimals = 0 
#     ) %>%
#     gt::cols_align(
#       align = "left",
#       columns = "Metric"
#     ) %>%
#     gt::cols_align(
#       align = "right",
#       columns = "Value"
#     ) %>%
#     gt::tab_style(
#       style = list(
#         gt::cell_text(weight = "bold")
#       ),
#       locations = gt::cells_column_labels()
#     )
#   
#   # Save as PDF
#   gt::gtsave(table, output_file)
#   
#   return(table)
# }

# Differential expression ----
#' Calculate coefficient of variation.
#' 
#' @param x A vector.
#' @examples
#' calc_CV(vector)
calc_CV <- function(x) {sd(x) / mean(x)}

# Spatial analysis ----
# For ComplexHeatmap() when identifying spatial niches
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid::grid.rect(x = x, y = y, width = width *0.99, 
                  height = height *0.99,
                  gp = grid::gpar(col = "grey", 
                                  fill = fill, lty = 1, lwd = 0.5))
}

col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

zscore <- function(x) {
  return((x - mean(x))/sd(x))
}

categorize_cells <- function(cell_type_simple) {
  if(base::grepl("^(M(ac|o)|DC)\\|", cell_type_simple)) {
    return("Myeloid")
  } else if(base::grepl("^CD[4|8]", cell_type_simple)) {
    return("T cell")
  } else if(cell_type_simple == "Tumor|Intracellular_16S") {
    return("16S+ tumor")
  } else {
    return("Other")
  }
}