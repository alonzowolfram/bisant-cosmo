#!/usr/bin/env Rscript
# https://training.nextflow.io/advanced/structure/#bin

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Setup -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Source the pipeline/setup.R file.
# Note that we cannot use source() directly in Nextflow; see https://stackoverflow.com/questions/76067626/how-to-call-a-function-present-in-a-different-r-script-from-an-r-executable-usin
# and https://community.seqera.io/t/source-another-r-script-in-the-bin-directory/1059
# So we will have to use the workaround suggested above if the workflow system is Nextflow.
cl_args <- commandArgs(TRUE)
workflow_system <- cl_args[2]
if(workflow_system=="Nextflow") {
  path <- Sys.getenv("PATH") |> strsplit(":")
  bin_path <- tail(path[[1]], n=1)
  source(file.path(bin_path, "pipeline/setup.R"))
} else {
  bin_path <- ""
  source("src/pipeline/setup.R") # I guess the path it sources from is the current working directory, not the path the R script lives in.
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Data import -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/vignette-basic-analysis/
# https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/vignette-basic-analysis/assets/0.-loading-flat-files.html
# Automatically get slide names:
message("Getting slide names from the flat files.")
slide_names <- dir(flat_file_dir)
# Also get the slide names inside run_summaries_dir, because we will need it for the run ID.
message("Getting slide names from the raw files (inside run_summaries_dir).")
run_summaries_slides <- list.files(run_summaries_dir)

# Lists to collect the counts matrices and metadata, one per slide.
count_list <- vector(mode = 'list', length = length(slide_names)) 
metadata_list <- vector(mode = 'list', length = length(slide_names)) 

message("Importing slides.")
for(i in 1:length(slide_names)) {
  slide_name <- slide_names[i] 
  
  msg <- paste0("Loading slide ", slide_name, ", ", i, "/", length(slide_names), ".")
  message(msg)
  
  # Slide-specific files:
  slide_i_files <- dir(paste0(flat_file_dir, "/", slide_name))
  
  # Load in metadata:
  slide_i_metadata <- slide_i_files[base::grepl("metadata\\_file", slide_i_files)]
  temp_data_table <- data.table::fread(paste0(flat_file_dir, "/", slide_name, "/", slide_i_metadata), fill = TRUE)
  
  # Numeric slide ID 
  slide_ID_numeric <- temp_data_table[1,]$slide_ID 
  
  # Load in counts as a data table:
  slide_i_counts <- slide_i_files[base::grepl("exprMat\\_file", slide_i_files)]
  counts_file <- paste0(flat_file_dir, "/", slide_name, "/", slide_i_counts)
  nonzero_elements_per_chunk <- 5*10**7
  ### Safely read in the dense (0-filled ) counts matrices in chunks.
  ### Note: the file is gzip compressed, so we don't know a priori the number of chunks needed.
  
  total_lines <- R.utils::countLines(counts_file)[1]
  chunk_id <- 1
  skiprows <- 0
  cell_count <- 0
  
  required_cols <- data.table::fread(counts_file, select=c("fov", "cell_ID"), fill = TRUE)
  stopifnot("columns 'fov' and 'cell_ID' are required, but not found in the counts file" =
              all(c("cell_ID", "fov") %in% colnames(required_cols)))
  number_of_cells <- nrow(required_cols)

  number_of_cols <-  ncol(data.table::fread(counts_file, nrows = 2))
  number_of_chunks <- ceiling(as.numeric(number_of_cols) * as.numeric(number_of_cells) / (nonzero_elements_per_chunk)) # Because of INTEGER OVERFLOW, must convert integers to numeric to avoid running out of space when doing integer multiplication (R only does 32-bit for integers I guess.)
  chunk_size <- floor(number_of_cells / number_of_chunks)
  sub_counts_matrix <- vector(mode='list', length=number_of_chunks)
  
  pb <- txtProgressBar(min = 0, max = number_of_chunks, initial = 0, char = "=", style = 3)
  
  while (TRUE) {
    read_header <- (chunk_id == 1)
    actual_skip <- skiprows + (chunk_id > 1)
    
    if (actual_skip >= total_lines) {
      message("Reached end of file: trying to skip ", actual_skip, " but file only has ", total_lines, " lines.")
      break
    }
    
    counts_data_table <- tryCatch({
      data.table::fread(counts_file,
                        nrows = chunk_size,
                        skip = actual_skip,
                        header = read_header,
                        fill = TRUE)
    }, error = function(e) {
      message(sprintf("Chunk %d: fread error: %s", chunk_id, e$message))
      return(NULL)
    })
    message(sprintf("Chunk %d: Read %d rows × %d cols", chunk_id, nrow(counts_data_table), ncol(counts_data_table)))
    
    if (is.null(counts_data_table) || nrow(counts_data_table) == 0) break
    
    # Set the header
    if (chunk_id == 1) header <- colnames(counts_data_table)
    
    # Manually pad columns (of `counts_data_table`) to match header length if necessary before setting the column names of the `counts_data_table`
    if (ncol(counts_data_table) != length(header)) {
      missing_cols <- length(header) - ncol(counts_data_table)
      if (missing_cols > 0) {
        # Add NA columns to pad to full width
        for (k in seq_len(missing_cols)) {
          counts_data_table[[paste0("__NA_pad_", k)]] <- NA_real_
        }
      } else {
        message(sprintf("More columns than expected: got %d, expected %d", 
                     ncol(counts_data_table), length(header), ". Skipping to the next slide"))
        next
      }
    }
    
    # Now that we've checked, we can set the column names of the `counts_data_table`
    if(chunk_id != 1) colnames(counts_data_table) <- header
    
    counts_data_table[is.na(counts_data_table)] <- 0
    
    slide_fov_cell_counts <- paste0("c_", i, "_", slide_ID_numeric, "_", counts_data_table$fov, "_", counts_data_table$cell_ID)
    sub_counts_matrix[[chunk_id]] <- as(counts_data_table[ , -c("fov", "cell_ID"), with = FALSE], "sparseMatrix")
    rownames(sub_counts_matrix[[chunk_id]]) <- slide_fov_cell_counts
    
    cell_count <- cell_count + nrow(counts_data_table)
    skiprows <- skiprows + chunk_size
    setTxtProgressBar(pb, chunk_id)
    chunk_id <- chunk_id + 1
  }
  
  close(pb)
  
  count_list[[i]] <- do.call(rbind, sub_counts_matrix) 
  
  # 2024/10/18: add a column for run ID.
  # This will be taken from the raw files directories.
  # Check that the current slide has run summary data. 
  path_regex <- paste0(slide_name, "\\/.+\\/RunSummary\\/Run_.+_ExptConfig\\.txt$") # Regex for the path to the ExptConfig.txt file. # 2025/02/08: $ added after .txt because sometimes there are >1 files.
  run_sum_path <- list.files(run_summaries_dir, full.names = TRUE, recursive = TRUE) %>% regexPipes::grep(path_regex, value = TRUE)
  if(length(run_sum_path) < 1) {
    message(glue::glue("ExptConfig.txt file not found for slide {slide_name}. A random ID will be generated for this slide"))
    set.seed(random_seed + i)
    run_id <- stringi::stri_rand_strings(1, 37, pattern = "[A-Za-z0-9]")
  } else {
    run_sum_path <- run_sum_path %>% .[1]
    # Read in the first line from the file, then extract the run ID.
    con <- file(run_sum_path, "r")
    run_id <- readLines(con, n = 1) %>% regexPipes::gsub("Run Number: ", "")
    close(con)
  }
  # Add the run ID to the metadata. 
  temp_data_table$run_ID <- run_id
  
  # 2024/10/23: replace the slide_name_var (usually Run_Tissue_name) with the name of the slide taken from the folder name, for consistency (since it may be different in the actual metadata flat file.)
  temp_data_table[[slide_name_var]] <- slide_name
  
  # Ensure that cell order in counts matches cell order in metadata.
  slide_fov_cell_metadata <- paste0("c_", i, "_", slide_ID_numeric, "_", temp_data_table$fov, "_", temp_data_table$cell_ID) # As of 2024/10/18, we are adding i to the identifier, which is the slide # in a given data analysis. slide_ID_numeric is the slide # in a given CosMx _run_ (physical). This will allow us to merge data from slides that come from different runs and may therefore have the same slide_ID_numeric. We did the same with slide_fov_cell_counts above. 
  # Original naming convention: c_[slide ID]_[fov ID]_[cell ID]
  # https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/assets/Pancreas-CosMx-ReadMe.html
  # 2025/06/10: The following code is necessary to check if the metadata and counts have the same number of rows.
  idx <- match(slide_fov_cell_metadata, rownames(count_list[[i]]))
  if (anyNA(idx)) {
    message(sprintf("Some metadata cell IDs were not found in the counts matrix (n = %d)", sum(is.na(idx))))
    # Optionally drop those metadata rows:
    valid_rows <- !is.na(idx)
    slide_fov_cell_metadata <- slide_fov_cell_metadata[valid_rows]
    idx <- idx[valid_rows]
    temp_data_table <- temp_data_table[valid_rows, ]
  }
  count_list[[i]] <- count_list[[i]][idx, , drop = FALSE]
  # End new code 2025/06/10.
  # count_list[[i]] <- count_list[[i]][match(slide_fov_cell_metadata, rownames(count_list[[i]])),] 
  identical(slide_fov_cell_metadata, rownames(count_list[[i]]))
  # 2024/10/18: replace the old cell_id with the new one (slide_fov_cell_metadata).
  temp_data_table$updated_cell_id <- slide_fov_cell_metadata
  metadata_list[[i]] <- temp_data_table 
  
  # Track common genes and common metadata columns across slides.
  if(i == 1){
    shared_genes <- colnames(count_list[[i]]) 
    shared_columns <- colnames(temp_data_table)
  } else {
    shared_genes <- intersect(shared_genes, colnames(count_list[[i]]))
    shared_columns <- intersect(shared_columns, colnames(temp_data_table))
  }
  
}

# Reduce to shared metadata columns and shared genes.
message("Reducing data to shared metadata columns and shared genes.")
for(i in 1:length(slide_names)){
  metadata_list[[i]] <- metadata_list[[i]][, ..shared_columns]
  count_list[[i]] <- count_list[[i]][, shared_genes]
}
counts <- do.call(rbind, count_list)
metadata <- data.table::rbindlist(metadata_list)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Cleaning -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

message("Cleaning data and generating Seurat object.")

# # FOV names per slide.
# slide_1_fovs <- metadata %>% dplyr::filter(slide_ID_numeric==1) %>% .$fov %>% unique
# slide_2_fovs <- metadata %>% dplyr::filter(slide_ID_numeric==2) %>% .$fov %>% unique
# intersect(slide_1_fovs, slide_2_fovs)
# # OK, so FOV names aren't unique to each slide. Good to know. 
# So we'll create unique FOV identifiers.
message("Creating unique FOV identifiers.")
metadata$uniqueFOV <- paste0(metadata[[slide_name_var]], ":", metadata[[fov_name_var]])

# Isolate negative control matrices:
message("Isolating negative control matrices.")
neg_counts <- counts[, base::grepl("Negative", colnames(counts))]
false_counts <- counts[, base::grepl("SystemControl", colnames(counts))]
if(ncol(neg_counts) < 1) warning("Warning: no negative counts.")

# NEW, as of 2024/11/03.
# Isolate bacterial counts:
message("Isolating bacterial counts.")
bac_counts <- counts %>% .[, intersect(colnames(.), bacterial_probes)]
if(length(intersect(colnames(counts), bacterial_probes)) < 1) message("Warning: no bacterial counts.")
  
# Reduce counts matrix to only genes:
# But first, make a copy of the colnames of the original counts matrix.
old_counts_colnames <- colnames(counts)
message("Reducing counts matrix to only genes.")
counts <- counts[, !base::grepl("Negative", colnames(counts)) & 
                   !base::grepl("SystemControl", colnames(counts)) & 
                   !(colnames(counts) %in% bacterial_probes)]

# Build a new Seurat object.
message("Building new Seurat object.")
seu.obj <- CreateSeuratObject(counts = t(counts), assay = 'RNA')
# Check if there are any negative probes before creating a new assay object for it. 
if(ncol(neg_counts) > 0) seu.obj[["negprobes"]] <- CreateAssayObject(counts = t(neg_counts)) # https://rdrr.io/github/igordot/scooter/src/R/import.R, add_seurat_assay
# Check if there are any bacterial probes before creating a new assay object for it.
if(length(intersect(old_counts_colnames, bacterial_probes)) > 0) seu.obj[["bacprobes"]] <- CreateAssayObject(counts = t(bac_counts))
seu.obj@meta.data <- metadata
if(identical(Cells(seu.obj), seu.obj@meta.data$updated_cell_id)) {rownames(seu.obj@meta.data) <- seu.obj@meta.data$updated_cell_id} else {warning("Cells do not match updated_cell_id in metadata. Please correct this.")} 

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Additional metadata -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Only join if both `additional_metadata` and `additional_metadata_join_column` are present. 
if(!flagVariable(additional_metadata_join_column) && !flagVariable(additional_metadata)) {
  message("Joining additional metadata.")
  
  if(file.exists(additional_metadata)) {
    # Load the metadata.
    additional_metadata_df <- read.table(additional_metadata, header = TRUE, sep = ",")
    
    # Check that `additional_metadata_join_column` is present in the metadata. 
    if((additional_metadata_join_column %in% colnames(additional_metadata_df)) && (additional_metadata_join_column %in% colnames(seu.obj@meta.data))) {
      seu.obj@meta.data <- seu.obj@meta.data %>% dplyr::left_join(additional_metadata_df, by = additional_metadata_join_column)
      rownames(seu.obj@meta.data) <- seu.obj@meta.data$updated_cell_id
    }
  }
  
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Filtering -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Filter observations by the criteria given in filter_vars, if applicable.
# filter_vars <- "segment,include,Full ROI,Tumor;Tags,exclude,Stroma"

if(!flagVariable(filter_vars)) {
  message("Filtering data.")
  
  filter_vars <- filter_vars %>% strsplit(";") %>% unlist
  for(var in filter_vars) {
    var <- var %>% strsplit(",") %>% unlist
    if(!(var[1] %in% colnames(seu.obj@meta.data))) {
      warning(paste0("The variable ", var[1], " is not included in the metadata for this data set. Skipping this variable."))
      next
    }
    
    var_values <- var[3:length(var)]
    if((var[2] %>% str_to_lower()) == "include") {
      cells_to_keep <- rownames(seu.obj@meta.data)[seu.obj@meta.data[[var[1]]] %in% var_values]
      seu.obj <- subset(seu.obj, cells = cells_to_keep)
      
    } else if((var[2] %>% str_to_lower()) == "exclude") {
      cells_to_keep <- rownames(seu.obj@meta.data)[!(seu.obj@meta.data[[var[1]]] %in% var_values)]
      seu.obj <- subset(seu.obj, cells = cells_to_keep)
      
    } else {
      warning("Skipping this element of filter_vars. Please enter a value of 'include' or 'exclude' as the second element of filter_vars. And check your spelling!")
      next
      
    }
    
  }
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Neovariable generation -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Create the specified neovariables (if provided) and add to the metadata.
if(!flagVariable(neovariables)) {
  message("Generating neovariables.")
  
  neovariables <- neovariables %>% .[. != ""]
  for(neovariable in neovariables) {
    # Split into its component parts.
    neovariable_comps <- neovariable %>% strsplit("\\+") %>% unlist
    
    # Create the neovariable.
    neovariable_name <- paste0(neovariable_comps, collapse = "_")
    seu.obj@meta.data <- seu.obj@meta.data %>% tidyr::unite(!!as.name(neovariable_name), neovariable_comps, remove = FALSE, na.rm = FALSE)
    # Set the rownames.
    # https://github.com/satijalab/seurat/issues/3695#issuecomment-731367783 
    rownames(seu.obj@meta.data) <- seu.obj@meta.data$updated_cell_id
  }
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

message("Exporting Seurat object to RDS file and cleaning up memory.")

# Save Seurat object.
saveRDS(seu.obj, paste0(output_dir_rdata, "seuratObject_raw.rds"))
# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)

# Clean up.
rm(counts, old_counts_colnames, count_list, metadata, metadata_list)
gc()