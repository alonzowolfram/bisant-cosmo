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
harmonized_data_path <- cl_args[5]
seu.obj.harmonized <- readRDS(harmonized_data_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Neighbor networks and clustering --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Find neighbors --------------------------------------
options(future.globals.maxSize = 1000 * 1024^2) # https://satijalab.org/seurat/archive/v3.0/future_vignette.html, https://github.com/satijalab/seurat/issues/1845
reduction <- ifelse(correct_batch_effects, "harmony", "pca")
system.time({
  seu.obj.clustered <- FindNeighbors(seu.obj.harmonized, 
                                     reduction = reduction,
                                     dims = 1:30,
                                     prune.SNN = jaccard_cutoff,
                                     nn.method = "annoy",
                                     annoy.metric = dist_metric_nn
  )
})
# user  system elapsed 
# 138.293   2.288 141.172 

## Clustering --------------------------------------
system.time({
  seu.obj.clustered <- FindClusters(seu.obj.clustered,
                                    algorithm = 4, # 4 = Leiden, 1 = Louvain
                                    method = "igraph", # Default is "matrix" which is fine for small datasets. "igraph" avoids casting to dense matrix.
                                    resolution = cluster_resolution,
                                    random.seed = random_seed
  )
})
# user   system  elapsed 
# 2979.840  121.704 3114.610 
gc()

## Visualization --------------------------------------
# Plot the cell locations in 2D space as in the tissue.
# Extract cell metadata.
cell_coords <- seu.obj.clustered@meta.data %>% .[,c(dimension_name_vars,
                                                   slide_name_var,
                                                   tissue_var,
                                                   "seurat_clusters")]
# head(cell_coords) %>% View

# Plot clusters.
plot <- cell_coords %>% ggplot(aes(x = !!as.name(dimension_name_vars[1]), 
                                   y = !!as.name(dimension_name_vars[2]),
                                   color = seurat_clusters
                                   )) 
plot <- DSPplotScatter(plot,
                       point_size = 0.01, point_alpha = 0.05, geom_jitter = TRUE,
                       panel_border_color = "black", panel_border_linewidth = 1,
                       facet_wrap_var = tissue_var,
                       coord_equal = TRUE,
                       legend_position = "none",
                       x_lab = NULL, y_lab = NULL, title = "Cell coordinates in XY space, colored by cluster")
DSPexportPlot(plot, path_out = output_dir_pubs, filename = "cell_coords_cluster", dpi = 600, png_width = NULL, png_height = NULL, png_unit = "mm")

## Marker identification --------------------------------------
seu.obj.clustered.downsamp <- subset(seu.obj.clustered, downsample = downsampling_nn, seed = random_seed)
system.time(clust.markers <- FindAllMarkers(seu.obj.clustered.downsamp, only.pos = TRUE))
rm(seu.obj.clustered.downsamp)
gc()
# user  system elapsed 
# 507.051 133.762 642.003
clust.markers.filt <- clust.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) 

clust.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Save Seurat objects.
saveRDS(seu.obj.clustered, paste0(output_dir_rdata, "seuratObject_clustered.rds")) # Clustered.
# Save cluster markers.
write.csv(clust.markers, paste0(output_dir_tabular, "leiden_cluster_markers.csv")) # Filtered cluster markers.
write.csv(top10, paste0(output_dir_tabular, "leiden_cluster_markers_top10.csv")) # All cluster markers.

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)