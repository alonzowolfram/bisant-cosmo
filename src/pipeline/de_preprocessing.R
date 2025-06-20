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
niche_data_path <- cl_args[5]
seu.obj.niches <- readRDS(niche_data_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Data preparation -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Coordinate offset ----
# # Visualize pre-offset coordinates
# seu.obj.niches@meta.data %>% 
#   ggplot() +
#   geom_point(aes(x = CenterX_global_px, y = CenterY_global_px, color = Run_Tissue_name))

# Offset coordinates by slide
message("Offsetting coordinates by slide")
seu.obj.niches@meta.data <- seu.obj.niches@meta.data %>% adjust_spatial_coordinates(slide_col = tissue_var, 
                                                                                    x_col = dimension_name_vars[1], 
                                                                                    y_col = dimension_name_vars[2])

# # Visualize post-offset coordinates
# seu.obj.niches@meta.data %>% ggplot() +
#   geom_point(aes(x = x_adj, y = y_adj, color = Run_Tissue_name))

## Inputs ----
message("Preparing inputs")
cts <- seu.obj.niches@assays$RNA$counts # Counts matrix
norm <- seu.obj.niches@assays$RNA$data # Total counts-normalized counts matrix
metainfo <- seu.obj.niches@meta.data # Metadata
# Clean up `metainfo` to match the expected column names due to what I guess is a bug in the smiDE code?
smide_var_names <- c("cell_type", "Run_Tissue_name", "cell_ID", "sdimx", "sdimy")
lookup <- c(cell_type = "InSituType_clusters", 
            Run_Tissue_name = tissue_var,
            cell_ID = "updated_cell_id",
            sdimx = "x_adj", # Using the offset coordinates instead of the original ones
            sdimy = "y_adj"
)
# Remove entries from `lookup` for which the metadata variable already matches what smiDE is expecting.
lookup <- lookup[smide_var_names != lookup] 
# Rename columns in `metainfo` that match those in `names(lookup)`
tmp_col_names <- names(lookup); names(tmp_col_names) <- paste0(tmp_col_names, "_1")
metainfo <- metainfo %>% dplyr::rename(any_of(tmp_col_names)) %>% dplyr::rename(all_of(lookup))

# Get the neighbor radius in appropriate units.
# https://nanostring.com/wp-content/uploads/2023/09/SMI-ReadMe-BETA_humanBrainRelease.html
# "To convert to microns multiply the pixel value by 0.12 um per pixel. For FOV-local coordinates, the origin (0,0) is the top left of the FOV."
# (See same link for global vs local.)
if(coordinate_units=="um") {
  neighbor_distance_converted <- neighbor_distance_de
} else if(coordinate_units=="mm") {
  # Convert um to mm.
  neighbor_distance_converted <- neighbor_distance_de / 1000
} else if(coordinate_units=="px") {
  # Convert um to px.
  neighbor_distance_converted <- neighbor_distance_de / 0.12
} else {
  stop("Error: please set your coordinate units to be one of 'px', 'um', or 'mm'")
}

## Pre-DE ----
message("Performing pre-DE")
# Filter ref_cell_types to include only the ones that are present in the metainfo.
ref_cell_types <- intersect(ref_cell_types, metainfo$cell_type)
pre_de_obj <- pre_de(counts = cts,
                     normalized_dat = norm,
                     metadat = metainfo,
                     ref_celltype = ref_cell_types,
                     cell_type_metadata_colname = "cell_type",
                     split_neighbors_by_colname = "Run_Tissue_name",
                     mm_radius = neighbor_distance_converted,
                     verbose = TRUE)

## Cell contamination ----
message("Calculating cell contamination")
# system.time({
contam_metrics <- 
  contamination_ratio_metric(assay_matrix = norm[,metainfo$cell_ID],
                             metadata = metainfo,
                             cluster_col = "cell_type",
                             radius = neighbor_distance_converted,
                             split_neighbors_by_colname = "Run_Tissue_name",
                             sdimx_col = "sdimx",
                             sdimy_col = "sdimy",
                             verbose = TRUE
  )
# })
# user   system  elapsed 
# 891.813  326.717 1217.979 

# Narrow the list of genes to be analyzed to the ones that do not show contamination from neighboring cells.
analysis_genes_list <- list()
for(cell_type_i in unique(metainfo$cell_type)) {
  analysis_genes_list[[cell_type_i]] <- contam_metrics[cell_type==cell_type_i & ratio <= 1]$target
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Exporting results to disk")

# Save pre-DE objects.
pre_de <- list()
for(obj in c("cts", "norm", "metainfo", "pre_de_obj", "contam_metrics", "analysis_genes_list")) {
  pre_de[[obj]] <- get(obj)
}
saveRDS(pre_de, paste0(output_dir_rdata, "pre_DE.rds"))
# Save Seurat object with adjusted spatial coordinates.
saveRDS(seu.obj.niches, paste0(output_dir_rdata, "seuratObject_spatial.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)