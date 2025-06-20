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
filt_data_path <- cl_args[5]
seu.obj.filt <- readRDS(filt_data_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Normalization -------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Normalization --------------------------------------
# system.time(seu.obj <- SCTransform(seu.obj, assay = "RNA"))
# system.time({
#   seu.obj <- NormalizeData(
#   seu.obj,
#   assay = "RNA",
#   normalization.method = "LogNormalize",
#   verbose = TRUE
#   )
#   })

# Per https://satijalab.org/seurat/reference/normalizedata,
# "Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p."
# NOTE: NanoString recommended NOT non-linearly transforming CosMx data, though ... 

# Extract the expression matrices.
counts_host <- seu.obj.filt@assays$RNA$counts %>% as.matrix # as.numeric needed if using .colSums instead of colSums below.
if("bacprobes" %in% names(seu.obj.filt)) { 
  counts_bac <- seu.obj.filt@assays$bacprobes$counts %>% as.matrix
  counts_total <- rbind(counts_host, counts_bac)
} else {
  counts_total <- counts_host
}
counts_neg <- seu.obj.filt@assays$negprobes$counts %>% as.matrix
gc()
# Rows are genes, cells are columns.

# Normalize total counts.
totalcounts <- Matrix::rowSums(counts_total)
scaling_factor <- pmax(totalcounts, min_counts_per_cell)
norm <- sweep(counts_total, 1, scaling_factor, "/")
saveRDS(norm, paste0(output_dir_rdata, "norm-unscaled.rds"))
# Normalize host counts.
totalcounts_host <- Matrix::rowSums(counts_host)
scaling_factor_host <- pmax(totalcounts_host, min_counts_per_cell)
norm_host <- sweep(counts_host, 1, scaling_factor_host, "/")
# Save normalized host counts.
saveRDS(norm_host, paste0(output_dir_rdata, "norm-host-unscaled.rds"))

# Apply scaling factors for host normalization.
metadata <- seu.obj.filt@meta.data
scaling_factors <- list(
  host_16S = mean(metadata$nCount_RNA, na.rm = TRUE), # Mean RNA counts for total scaling.
  host = mean(colSums(counts_total), na.rm = TRUE) # Mean total counts for host scaling.
)
# Scale host normalization using scaling factor.
scaled_norm_host <- norm_host * scaling_factors$host
saveRDS(scaled_norm_host, paste0(output_dir_rdata, "norm-host-scaled-with-host-factor.rds"))

# If there are bacterial (16S) counts, normalize them.
if("bacprobes" %in% names(seu.obj.filt)) {
  total_counts_bac <- Matrix::rowSums(counts_bac, na.rm = TRUE)
  
  # Calculate direct background rate from negative probes
  background_rate <- mean(counts_neg[!is.na(counts_neg)], na.rm = TRUE) # Are there bacteria-specific negative probes?
  print(paste("Background rate (mean of negative probes):", background_rate))
  
  # Normalize bacterial counts by background rate
  norm_bac <- as.matrix(counts_bac) / background_rate # Normalize by dividing by the background rate
  
  # Convert normalized bacterial counts back to sparse matrix
  norm_bac <- Matrix::Matrix(norm_bac, sparse = TRUE)
  
  # Preserve row and column names in the normalized matrix
  rownames(norm_bac) <- rownames(counts_bac)
  colnames(norm_bac) <- colnames(counts_bac)
  
  # Save the normalized bacterial counts
  saveRDS(paste0(output_dir_rdata, "norm-bac-unscaled.rds"))
}

# OLD NORMALIZATION CODE:
# # Divide feature counts for each cell by the total counts for that cell.
# # https://stackoverflow.com/a/70965486 - memory efficiency
# system.time(exprs_norm <- sweep(exprs, 2, .colSums(exprs, m = nrow(exprs), n = ncol(exprs)), `/`)) # https://stackoverflow.com/a/9447920/23532435
# # user  system elapsed 
# # 83.720  35.513 137.341 
# # Warning message:
# #   In asMethod(object) :
# #   sparse->dense coercion: allocating vector of size 20.3 GiB
# gc()
# # Multiply by the scale factor.
# scale_factor <- 10000
# exprs_norm <- exprs_norm * scale_factor
# END OLD NORMALIZATION CODE.

# NOTE: NanoString recommends NOT non-linearly transforming CosMx data.
# Natural-log transform.
scaled_norm_host_log <- log1p(scaled_norm_host)
gc()

## Seurat object creation --------------------------------------
# Create Seurat object from normalised counts.
seu.obj.norm <- CreateSeuratObject(counts = counts_host, assay = 'RNA')
gc()
# Add to the RNA assay slot.
# https://satijalab.github.io/seurat-object/reference/Layers.html
if(log_transform) {
  LayerData(seu.obj.norm, layer = "data", assay = "RNA") <- scaled_norm_host_log # https://github.com/satijalab/seurat/wiki/Assay
} else {
  LayerData(seu.obj.norm, layer = "data", assay = "RNA") <- scaled_norm_host
}
# Add in the metadata and other assays.
seu.obj.norm[["negprobes"]] <- CreateAssayObject(counts = counts_neg) # https://rdrr.io/github/igordot/scooter/src/R/import.R, add_seurat_assay
if("bacprobes" %in% names(seu.obj.filt)) {
  seu.obj.norm[["bacprobes"]] <- CreateAssayObject(counts = counts_bac)
  LayerData(seu.obj.norm, layer = "data", assay = "bacprobes") <- norm_bac
}
seu.obj.norm@meta.data <- seu.obj.filt@meta.data

# Scale and center. 
system.time(seu.obj.norm <- FindVariableFeatures(seu.obj.norm,
                                                 assay = "RNA",
                                                 selection.method = "vst", 
                                                 nfeatures = n_variable_features, 
                                                 verbose = TRUE))
# user  system elapsed 
# 19.308   8.467  27.661
system.time(seu.obj.norm <- ScaleData(seu.obj.norm, assay = "RNA")) # https://github.com/satijalab/seurat/issues/2699 # Also, why does this use so much memory???
# user  system elapsed 
# 26.974  19.868  46.668 
gc()

## Cleanup --------------------------------------
rm(counts_bac, counts_host, counts_neg)
gc()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Export Seurat object.
saveRDS(seu.obj.norm, paste0(output_dir_rdata, "seuratObject_normalized.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)