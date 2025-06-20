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

# Load the Seurat object created in the last step. 
raw_data_path <- cl_args[5]
seu.obj <- readRDS(raw_data_path)

# Check if pre-QC mode is on. 
# If so, override the other general QC settings.
if(pre_qc_mode) {
  filter_cells <- FALSE
  filter_fovs <- FALSE
  filter_neg_probes <- FALSE
  filter_targets <- FALSE
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 16S stats -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if("bacprobes" %in% names(seu.obj.raw@assays)) {
  
  # Step 1: Extract the counts matrix from the 'bacprobes' assay, if there even is one
  bac_counts <- GetAssayData(seu.obj.raw, assay = "bacprobes", slot = "counts")
  
  if(!flagVariable(bacterial_probes_histogram)) {
    if(any(bacterial_probes_histogram %in% rownames(bac_counts))) {
      # At least one of the provided bacterial probes in `bacterial_probes_histogram` is present
      # Are all of them?
      if(all(bacterial_probes_histogram %in% rownames(bac_counts))) {
        # All of them are present
        message("Subsetting bacterial counts")
        bac_counts <- bac_counts %>% .[bacterial_probes_histogram, , drop = FALSE]
        
      } else {
        # Not all of them are present
        probes_not_present <- setdiff(bacterial_probes_histogram, rownames(bac_counts)) # setdiff(x, y): present in x, not in y
        message(paste("Not all of the probes provided by `bacterial_probes_histogram` are present in the bacterial counts. The following were not present and are excluded:", paste(probes_not_present, sep = ", ")))
        message("Subsetting bacterial counts")
        bac_counts <- bac_counts %>% .[bacterial_probes_histogram, , drop = FALSE]
      }
      
    } else {
      # None of the provided bacterial probes in `bacterial_probes_histogram` are present
      message("None of the probes provided by `bacterial_probes_histogram` are present in the bacterial counts. All present probes will be used for 16S QC and stats")
      
    } # End `if()`: the probes in `bacterial_probes_histogram` are actually in the counts matrix
    
  } # End `if()`: `bacterial_probes_histogram` is not NULL/empty
  
  # Step 2: Summarize total bacterial probe counts per cell
  cell_bac_sums <- Matrix::colSums(bac_counts)
  
  # Step 3: Create a data frame with cell IDs and total counts
  bac_df <- data.frame(
    cell_id = names(cell_bac_sums),
    total_bac_count = as.numeric(cell_bac_sums)
  )
  
  # Step 4: Add slide information from the Seurat object's metadata
  bac_df$Run_Tissue_name <- seu.obj.raw@meta.data[["Run_Tissue_name"]]
  
  # Step 5: Add measurements of central tendency and spread
  slide_stats <- bac_df %>%
    group_by(Run_Tissue_name) %>%
    summarise(
      median = median(total_bac_count), # Median will be more robust for 0-inflated counts
      q25 = quantile(total_bac_count, 0.25),
      q75 = quantile(total_bac_count, 0.75),
      .groups = "drop"
    )
  
  # Step 6: Plot histograms faceted by slide
  bacterial_counts_plot <- ggplot(bac_df, aes(x = total_bac_count)) +
    geom_histogram(bins = 50, fill = "#999", color = "black") +
    facet_wrap(~Run_Tissue_name, scales = "free_y") +
    geom_vline(data = slide_stats, aes(xintercept = median), color = "red", linetype = "dashed", size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    geom_rect(
      data = slide_stats,
      aes(xmin = q25, xmax = q75, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      alpha = 0.2, fill = "orange"
    ) +
    labs(
      x = "Total bacterial probe counts per cell",
      y = "Number of cells",
      title = "Distribution of bacterial probe counts by slide"
    )
  
  # Step 6: Generate table of counts
  bacterial_counts_table <- table(bac_df$Run_Tissue_name, bac_df$total_bac_count)
  
} # End `if()`: there are bacterial probes

# Create the list object to save (will be read in by `report_template.Rmd`)
bacterial_counts_list <- list()
if(exists("bacterial_counts_plot")) {
  bacterial_counts_list[["plot"]] <- bacterial_counts_plot
} else {
  bacterial_counts_list[["plot"]] <- NULL
}
if(exists("bacterial_counts_table")) {
  bacterial_counts_list[["table"]] <- bacterial_counts_table
} else {
  bacterial_counts_list[["table"]] <- NULL
}
if(!flagVariable(bacterial_probes_histogram)) {
  bacterial_counts_list[["probes"]] <- bacterial_probes_histogram
} else {
  bacterial_counts_list[["probes"]] <- NULL
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Cell QC ---------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Performing cell QC")

# The two most important variables to filter on are
# 1) signal (mean transcripts/cell) and
# 2) cell area. 

# Extract cell metadata.
metadata <- seu.obj@meta.data %>% as.data.frame
if(identical(Cells(seu.obj), seu.obj@meta.data$updated_cell_id)) {rownames(metadata) <- metadata$updated_cell_id} else {warning("Double-check the cell IDs.")}
dim(metadata)
head(metadata)

## Minimum total counts ----------------------------------------------------------------
message("Creating flags based on total counts")
# Sum up all the counts column-wise.
cts <- as(seu.obj@assays$RNA$counts, "sparseMatrix")
system.time({
  col_sums <- SparseArray::colSums(seu.obj@assays$RNA$counts %>% as.matrix, na.rm = TRUE)
})
# user  system elapsed 
# 349.291  16.535 374.762 
gc() # Run gc() because it uses a lot of memory.
# Check each element of col_sums if it is >= min_counts_per_cell.
qcFlagsMinRNACounts <- ifelse(col_sums >= min_counts_per_cell, FALSE, TRUE)

# Plot histogram of total RNA counts.
metadata %>% 
  ggplot(aes(x = nCount_RNA + 1)) %>% 
  DSPplotHistogram(x_lab = "RNA counts + 1 (log10 scale)", y_lab = "Number of cells", vline_xintercept = min_counts_per_cell) %>% 
  DSPexportPlot(filename = "histogram")

## Mean counts per cell ----------------------------------------------------------------
message("Creating flags based on mean counts per cell")
# For each tissue type (because the mean [counts/cell] will vary by tissue)
# determine which cells are outliers (either too many or too few counts).
tissue_types <- seu.obj@meta.data[[tissue_var]] %>% unique
outliers <- c()
for(tissue_type in tissue_types) {
  # Get all the cells belonging to that tissue type. 
  if(is.na(tissue_type)) {
    cells_i <- seu.obj@meta.data %>% dplyr::filter(is.na(!!as.name(tissue_var))) %>% rownames()
  } else {
    cells_i <- seu.obj@meta.data %>% dplyr::filter(!!as.name(tissue_var) == tissue_type) %>% rownames()
  }
  
  # Get the total counts for each of those cells.
  total_counts_per_cell <- col_sums %>% .[which(names(.) %in% cells_i)]
  
  # Detect outliers. 
  iqr <- IQR(total_counts_per_cell)
  # Outliers = all cells with counts < median - IQR or counts >= median + IQR.
  outliers_i <- total_counts_per_cell %>% between(median(total_counts_per_cell) - iqr, median(total_counts_per_cell) + iqr) %>% !.
  names(outliers_i) <- names(total_counts_per_cell)
  
  outliers <- c(outliers, outliers_i)
}
# match(x, y): For each element i of x, where (at what index) does i occur in y?
# So the element whose order we're trying to match (the standard) should be SECOND (i.e. y).
if(length(outliers) > 0) {
  outliers <- outliers %>% .[match(rownames(seu.obj@meta.data), names(outliers))]
  if(identical(names(outliers), rownames(seu.obj@meta.data))) {qcFlagsMeanRNACounts <- outliers} else {warning("Double-check the cell IDs.")}
} else {
  qcFlagsMeanRNACounts <- rep(FALSE, length(Cells(seu.obj)))
}

## Negative proportions ----------------------------------------------------------------
# Note that if a cell has 0 total counts, its resulting `ProportionNegative` will be 0/0 -> NaN
message("Creating flags based on proportion of negative probes")
cts_neg <- as(seu.obj@assays$negprobes$counts, "sparseMatrix")
col_sums_neg <- SparseArray::colSums(cts_neg, na.rm = TRUE)
gc()
total_per_cell_counts <- col_sums + col_sums_neg
neg_props <- col_sums_neg / total_per_cell_counts
# Add negative proportions to metadata.
seu.obj@meta.data[["ProportionNegative"]] <- neg_props
# Create flags.
qcFlagsCellPropNeg <- ifelse(neg_props < proportion_neg_counts, FALSE, TRUE)

## Unique genes/cell ----------------------------------------------------------------
message("Creating flags based on unique genes/cell")
# This is not a filtering parameter on AtoMx. It is used for visualization purposes only. 
# Extract scaled gene-expression matrix.
# https://r-charts.com/distribution/histogram-boxplot/
qcFlagsMinFeatures <- ifelse(seu.obj@meta.data$nCount_RNA >= min_features_per_cell, FALSE, TRUE)
# genes_per_cell_path <- paste0("data/analysis/", project_name, "_genes-detected-per-cell.png")
# if(Sys.info()['sysname']=="Linux" & base::grep("dragon", Sys.info()['nodename'])) { savePlot(genes_per_cell_path) } else {quartz.save(genes_per_cell_path) }

## Complexity (ratio of counts:unique genes/cell) ----------------------------------------------------------------
message("Creating flags based on complexity")
complexity <- seu.obj[["nCount_RNA"]][[1]] / seu.obj[["nFeature_RNA"]][[1]]
# Add complexity to metadata.
seu.obj@meta.data[["Complexity"]] <- complexity
# Create flags.
qcFlagsCellComplex <- ifelse(complexity >= count_dist, FALSE, TRUE)

## Signal strength (ratio of [counts:genes]:[counts:negative probes]) ----------------------------------------------------------------
message("Creating flags based on signal strength")
# This is not a filtering parameter on AtoMx; it was recommended to us by NanoString.
n_features <- seu.obj@assays$RNA$counts %>% dim %>% .[1]
n_features_neg <- seu.obj@assays$negprobes$counts %>% dim %>% .[1]
signal_strength <- calcSignalStrength(n_features, n_features_neg, seu.obj)
# Add signal strength to metadata.
seu.obj@meta.data[["SignalStrength"]] <- signal_strength
# Create flags.
qcFlagsSignalStrength <- ifelse(signal_strength >= min_signal_strength, FALSE, TRUE)

## Area ----------------------------------------------------------------
message("Creating flags based on area")
# "The Grubbs" test is designed for assessing one outlier only.If more outliers are suspected, alternative tests, such as the Tietjen-Moore test, are recommended."
# https://rpubs.com/DragonflyStats/Grubbs-Test-Outliers
grubbs_area <- outliers::grubbs.test(seu.obj@meta.data$Area, type = 11, opposite = FALSE, two.sided = TRUE)
if(grubbs_area$p.value < area_outlier_pval) {
  # Get the areas that are outliers.
  outlier_areas <- grubbs_area$alternative %>% 
    regexPipes::gsub("[[:alpha:]]", "") %>% 
    str_split("[[:space:]]") %>% 
    unlist %>% 
    .[. != ""] %>% 
    as.numeric
  qcFlagsCellArea <- ifelse(seu.obj@meta.data$Area %in% outlier_areas, TRUE, FALSE)
} else {
  qcFlagsCellArea <- FALSE
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# FOV QC ----------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Performing FOV QC")

## FOVs flagged ----------------------------------------------------------------
# https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/fov-qc/
# QC will be run on each slide/individual tissue type.
# Barcodes can be found here: https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/FOV-QC/code/FOV%20QC 
# Load barcodes.
if(workflow_system=="Nextflow") {
  project_directory <- appendSlashToPath(cl_args[6])
} else {
  project_directory <- ""
}
allbarcodes <- readRDS(paste0(project_directory, path_to_barcodes)) # cl_args[4] = Nextflow project directory. 
names(allbarcodes)
# Get the barcodes for the panel we want.
barcodemap <- allbarcodes[[probe_panel]]

slide_names <- metadata[[slide_name_var]] %>% unique
# Vector to hold IDs of flagged FOVs.
flagged_fov_unique_ids <- c()
# List to hold FOV QC results across a range of values for max_prop_loss.
max_prop_loss_range <- c(0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
fov_qc_res_list <- list()
for(slide_name in slide_names) {
  fov_qc_res_list[[slide_name]] <- list()
  
  # Subset metadata to include only the current slide. 
  metadata_sub <- metadata %>% dplyr::filter(!!as.name(slide_name_var) == slide_name)
  
  # Get the counts.
  # Cells as rows, genes as columns. 
  counts <- seu.obj@assays$RNA$counts %>% 
    t %>% 
    .[rownames(.) %in% rownames(metadata_sub),]
  # Get the FOVs (that each cell is assigned to).
  fov <- metadata_sub$uniqueFOV
  # Get the XY dimensions?
  xy <- metadata_sub[,dimension_name_vars]
  
  # If pre-QC mode is on, run FOV QC across a range of values for max_prop_loss.
  if(pre_qc_mode) {
    for(max_prop_loss_i in max_prop_loss_range) {
      fov_qc_res_list[[slide_name]][[paste0("max_prop_loss_", max_prop_loss_i)]] <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap, max_prop_loss = max_prop_loss_i)
    }
    
    # Unique identifiers of flagged FOVs.
    flagged_fov_unique_ids <- c(flagged_fov_unique_ids, fov_qc_res_list[[slide_name]][[paste0("max_prop_loss_", 0.3)]]$flaggedfovs)
    # paste0("s", metadata$slide_ID, "_f", metadata$fov)
    
  } else {
    # Run FOV QC for the user-provided max_prop_loss value.
    system.time(fov_qc_res <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap, max_prop_loss = max_prop_loss))
    length(fov_qc_res$flaggedfovs) / length(unique(fov))
    # user  system elapsed 
    # 82.884  13.316  72.550 
    mapFlaggedFOVs(fov_qc_res)
    
    # Unique identifiers of flagged FOVs.
    flagged_fov_unique_ids <- c(flagged_fov_unique_ids, fov_qc_res$flaggedfovs)
    # paste0("s", metadata$slide_ID, "_f", metadata$fov)
  }
  
  # Free up memory. 
  rm(metadata_sub, counts, fov, xy)
  gc()
  
}

# Now, using flagged_fov_unique_ids, set qcFlagsFOV.
qcFlagsFOV <- ifelse(seu.obj$uniqueFOV %in% flagged_fov_unique_ids, TRUE, FALSE)

## Cells flagged across FOVs ----------------------------------------------------------------
# I have no idea what this is, honestly. 
# If we're just going to remove entire FOVs, why does it matter if individual cells are flagged in those FOVs?

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Negative probe QC -----------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Outlier test: Grubbs.
# "The Grubbs test is designed for assessing one outlier only. If more outliers are suspected, alternative tests, such as the Tietjen-Moore test, are recommended."
# https://rpubs.com/DragonflyStats/Grubbs-Test-Outliers
neg_probe_outliers <- c()
neg_probe_geo_mean <- seu.obj@assays$negprobes$counts %>% rowMeans()
grubbs_neg_probe <- outliers::grubbs.test(neg_probe_geo_mean, type = 10, opposite = FALSE, two.sided = FALSE)
if(grubbs_neg_probe$p.value < neg_probe_outlier_pval) {
  # Get the negative probes that are outliers.
  outlier_value <- grubbs_neg_probe$alternative %>% 
    regexPipes::gsub("[[:alpha:]]", "") %>% 
    str_split("[[:space:]]") %>% 
    unlist %>% 
    .[. != ""] %>% 
    as.numeric
  neg_probe_outliers <- neg_probe_geo_mean[neg_probe_geo_mean==outlier_value] %>%
    names()
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Target QC -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Performing target QC")

# We can _inspect_ the (gene:count):(negative probe:count) ratio to see if it is >= 4,
# but in general, genes (targets) should not be _removed_.

# # Parameters
# neg_control_probe_quantile_cutoff <- 0.5
# detection_over_bg_p_value <- 0.01
# filter_targets_by_neg_control_quantile <- TRUE
# filter_targets_by_detection_p_value <- TRUE

flagged_targets <- c()
flagged_targets_neg <- c()
flagged_targets_bg <- c()

# Negative control
# Flag with fewer counts than the nth percentile of the means of the negative probes.
nth_percentile_neg_probes <- quantile(neg_probe_geo_mean, probs = neg_control_probe_quantile_cutoff, na.rm = TRUE)
raw_counts_geo_mean <- seu.obj@assays$RNA$counts %>% rowMeans()
if(filter_targets_by_neg_control_quantile) {
  flagged_targets_neg <- raw_counts_geo_mean %>% 
    .[which(. < nth_percentile_neg_probes)] %>%
    names
}

# Detection over background
# Flag probes not determined by background score test to be above background.
if(filter_targets_by_detection_p_value) {
  # Perform t-test for each gene.
  system.time(bg_score_p_values <- apply(seu.obj@assays$RNA$counts, 1, gene.t.test, neg_probe_geo_mean))
  # user  system elapsed 
  # 285.181  77.020 364.367
  gc()
  
  # Adjust p-values for multiple testing.
  adjusted_bg_score_p_values <- p.adjust(bg_score_p_values, method = "fdr")
  
  flagged_targets_bg <- adjusted_bg_score_p_values %>% 
    .[which(. < detection_over_bg_p_value)] %>%
    names
}

if(filter_targets_by_neg_control_quantile & !filter_targets_by_detection_p_value) {
  flagged_targets <- flagged_targets_neg
} else if(filter_targets_by_neg_control_quantile & !filter_targets_by_detection_p_value) {
  flagged_targets <- flagged_targets_bg
} else if(filter_targets_by_neg_control_quantile & filter_targets_by_detection_p_value) {
  flagged_targets <- intersect(flagged_targets_neg, flagged_targets_bg)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Cell and FOV flagging -------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Performing cell and FOV flagging")

## Cell QC ----------------------------------------------------------------
# Add cell QC flags.
for(flag in c("qcFlagsMinRNACounts", "qcFlagsMeanRNACounts", "qcFlagsMinFeatures", "qcFlagsCellPropNeg", "qcFlagsCellComplex", "qcFlagsSignalStrength", "qcFlagsCellArea")) {
  seu.obj[[flag]] <- get(flag) %>% unname
}
# Create a flag to indicate whether cell should be removed or not. 
cell_flag_sum <- qcFlagsMinRNACounts + qcFlagsCellPropNeg + qcFlagsCellComplex + qcFlagsCellArea
qcCellsFlagged <- ifelse(cell_flag_sum > 0, TRUE, FALSE)
seu.obj$qcCellsFlagged <- qcCellsFlagged %>% unname

## FOV QC ----------------------------------------------------------------
# Add FOV QC flags.
seu.obj$qcFlagsFOV <- qcFlagsFOV

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Filtering -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Filtering data")

# If filter_neg_probes is TRUE, filter out flagged negative probes. 
# https://github.com/satijalab/seurat/issues/321
if(filter_neg_probes) {
  neg_probes_filt <- seu.obj@assays$negprobes$counts %>% .[!(rownames(.) %in% neg_probe_outliers),]
} else {
  neg_probes_filt <- seu.obj@assays$negprobes$counts
}

# If filter_targets is TRUE, filter out flagged targets.
# https://github.com/satijalab/seurat/issues/321
if(filter_targets) {
  counts_filt <- seu.obj@assays$RNA$counts %>% .[!(rownames(.) %in% flagged_targets),]
} else {
  counts_filt <- seu.obj@assays$RNA$counts
}

# Build a new seurat object from the filtered counts matrices.
seu.obj.filt <- CreateSeuratObject(counts = counts_filt, assay = 'RNA')
seu.obj.filt[["negprobes"]] <- CreateAssayObject(counts = neg_probes_filt) # https://rdrr.io/github/igordot/scooter/src/R/import.R, add_seurat_assay
if("bacprobes" %in% names(seu.obj.filt)) seu.obj.filt[["bacprobes"]] <- CreateAssayObject(counts = seu.obj@assays$bacprobes$counts)

# Add the metadata to seu.obj.filt.
seu.obj.filt@meta.data <- seu.obj@meta.data
# Set the rownames.
# https://github.com/satijalab/seurat/issues/3695#issuecomment-731367783 
rownames(seu.obj.filt@meta.data) <- seu.obj.filt@meta.data$updated_cell_id

# If filter_cells is TRUE, filter out flagged cells.
if(filter_cells) {
  seu.obj.filt <- subset(seu.obj.filt, subset = qcCellsFlagged == FALSE)
}

# If filter_fovs is TRUE, filter out flagged FOVs. 
if(filter_fovs) {
  seu.obj.filt <- subset(seu.obj.filt, subset = qcFlagsFOV == FALSE)
}

gc()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# QC table (cells and FOVs) -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Creating QC table")

# Run Brenda's function (calculate_qc_stats(), in helper_functions.R)
# The function calculate_qc_stats() takes as input a data frame
# where each row corresponds to 1 cell
# and that has columns for total counts, cell area, and cells/FOV.
stats_df <- calculate_qc_stats(seu.obj@meta.data, total_cells_post_qc = nrow(seu.obj.filt@meta.data), total_counts_col = total_counts_var, total_features_col = total_features_var, area_col = area_var, fov_col = fov_name_var,
                               min_counts_per_cell = min_counts_per_cell, min_features_per_cell = min_features_per_cell, max_proportion_neg_counts = proportion_neg_counts, min_count_dist = count_dist, max_area = max_area, min_signal_strength = min_signal_strength,
                               min_cells_per_fov = min_cells_per_fov)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Exporting results to disk")

# Save Seurat objects.
saveRDS(seu.obj.filt, paste0(output_dir_rdata, "seuratObject_filtered.rds")) # Filtered data. 
# Save metrics.
saveRDS(bacterial_counts_list, paste0(output_dir_rdata, "bac-probe_stats.rds")) # List containing bacterial probe counts plot, table, and the probes selected
saveRDS(stats_df, paste0(output_dir_rdata, "qc_metrics.rds"))
saveRDS(seu.obj@meta.data, paste0(output_dir_rdata, "flagged_metadata_raw.rds")) # Raw metadata but with flags.
saveRDS(fov_qc_res_list, paste0(output_dir_rdata, "fov_qc_res_list.rds")) # A list containing FOV QC res at different values of max_prop_loss.
# Save negative probes.
saveRDS(neg_probes_filt, paste0(output_dir_rdata, "neg-probes_filtered.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)