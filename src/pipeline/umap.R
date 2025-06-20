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
clustered_data_path <- cl_args[5]
seu.obj.clustered <- readRDS(clustered_data_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# UMAP --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
reduction <- ifelse(correct_batch_effects, "harmony", "pca")
system.time(seu.obj.umap <- RunUMAP(seu.obj.clustered, 
                                    assay = "RNA",
                                    slot = "data",
                                    seed.use = random_seed,
                                    reduction = reduction,
                                    dims = 1:n_dims_umap, # OK, so apparently the # of dimensions is what makes the biggest difference in how much separation we get in the graph ... 
                                    n.neighbors = n_neighbors_umap,
                                    min.dist = min_dist_umap,
                                    spread = spread_umap,
                                    metric = dist_metric_umap
)
)
# user  system elapsed 
# 530.604   4.180 535.683 

# Plot UMAP.
plot <- DimPlot(seu.obj.umap,
        reduction = "umap")
DSPexportPlot(plot, path_out = output_dir_pubs, filename = paste0("umap"))

# # UMAP with batch covariates.
# seu.obj.umap@reductions$umap@cell.embeddings %>% as.data.frame() %>%
#   tibble::rownames_to_column(var = "updated_cell_id") %>%
#   dplyr::left_join(seu.obj.umap@meta.data %>% dplyr::select(updated_cell_id, run_ID)) %>%
#   .[sample(nrow(.)),] %>%
#   ggplot(aes(x = umap_1, y = umap_2, color = run_ID)) +
#   geom_point() +
#   theme_bw()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Save Seurat object.
saveRDS(seu.obj.umap, paste0(output_dir_rdata, "seuratObject_UMAP.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)