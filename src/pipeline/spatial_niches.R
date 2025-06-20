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
celltyped_data_path <- cl_args[5]
seu.obj.celltyping <- readRDS(celltyped_data_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Spatial clustering (niche detection) --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Defining spatial neighbors")
# Define neighbors
nng_subset_var <- ifelse(flagVariable(tissue_core_var), fov_name_var, tissue_core_var)
# system.time({
  neighbors_sparsematrix <- nearestNeighborGraph(
    x = seu.obj.celltyping@meta.data[[dimension_name_vars[1]]], 
    y = seu.obj.celltyping@meta.data[[dimension_name_vars[2]]], 
    N = 50,
    subset = seu.obj.celltyping@meta.data[[nng_subset_var]] # "Edges will only be created for cells that have the same subset value, usually the slide column id but could also be a slide plus FOV id to only create edges within an FOV."
  )
# })
# user  system elapsed 
# 23.792   0.853  24.741 

message("Creating matrix of mean expression of each cell's neighbors")
# Get matrix of mean expression in each cell's neighbors:
# system.time({
  neighborhoodexpression <- neighbor_colMeans(x = t(seu.obj.celltyping[["RNA"]]$data), 
                                              neighbors = neighbors_sparsematrix)
# })
# user  system elapsed 
# 247.704  14.068 262.437 

message("Reducing matrix to top PCs")
# Reduce the matrix to its top PCs:
# system.time({
  pc <- irlba::prcomp_irlba(neighborhoodexpression, n = 25)$x
# })
# user  system elapsed 
# 557.822  47.211 189.292

message("Clustering top PCs")
# Cluster the top PCs:
# system.time({
  temp <- mclust::Mclust(pc,
                         G = 10,
                         modelNames = c("EEI")
  )
# })
# user  system elapsed 
# 51.150   0.712  51.677 

# Add to metadata
spatialclust <- temp$classification
seu.obj.celltyping@meta.data$niche <- spatialclust

# Plot the results:
cols <- c(
  "#616161", "#4285f4", "#db4437", "#f4b400", "#8DD3C7", "#0f9d58", "#ab47bc", "#00acc1",
  "#ff7043", "#9e9d24", "#5c6bc0", "#f06292", "#00796b", "#c2185b", "#7e57c2",
  "#03a9f4", "#8bc34a", "#fdd835", "#fb8c00", "#8d6e63", "#9e9e9e", "#607d8b",
  "#BEBADA", "#FB8072", "#80B1D3", "#FDB462"
)

message("Plotting")
# ! Plot:
for(slide_name in unique(seu.obj.celltyping@meta.data[[slide_name_var]])) {
  ggplot(seu.obj.celltyping@meta.data %>% dplyr::filter(!!as.name(slide_name_var)==slide_name), aes(x = !!as.name(dimension_name_vars[1]), y = !!as.name(dimension_name_vars[2]), color = factor(niche))) +
    geom_point_rast(size = 0.1, alpha = 0.8) +
    scale_color_manual(values = cols) +
    coord_fixed(ratio = 1) + # This maintains the aspect ratio
    theme_bw(base_size = 10) +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    labs(color = "Niche", x = "x", y = "y") +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  ggsave(
    file = here(output_dir_pubs, paste0("UMAP_InSituType-clustering_", slide_name, ".pdf"))
  )
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# UMAP --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ! niche
# Change Identity
Idents(seu.obj.celltyping) <- "niche"

# Plot UMAP
pdf(here(output_dir_pubs, "UMAP_spatial-niches.pdf"),
    width = 8, height = 6
)
DimPlot(seu.obj.celltyping,
        reduction = "umap", label = TRUE, label.size = 4, repel = TRUE,
        cols = cols
)
dev.off()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Metadata for Napari --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
napari_vars_final <- c("updated_cell_id", 
                       fov_name_var, 
                       ifelse(!flagVariable(tissue_core_var), tissue_core_var, NA), 
                       "InSituType_clusters",
                       "niche") %>% .[!is.na(.)]
if(!flagVariable(napari_vars)) {
  tmp <- union(napari_vars_final, napari_vars)
  vars_not_present <- setdiff(tmp, colnames(seu.obj.celltyping@meta.data))
  if(length(vars_not_present) > 0) warning(paste0("The following variables are not present in the metadata: ", vars_not_present))
  napari_vars_final <- intersect(colnames(seu.obj.celltyping@meta.data), union(napari_vars_final, napari_vars))
}

df_metadata <- as.data.frame(seu.obj.celltyping@meta.data)%>%
  dplyr::select(any_of(napari_vars_final))

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Exporting to disk")

# Save Seurat object.
saveRDS(seu.obj.celltyping, paste0(output_dir_rdata, "seuratObject_spatial-niches.rds"))
# Save Napari metadata.
write.csv(df_metadata, here(output_dir_tabular, "metadata_Napari.csv"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)