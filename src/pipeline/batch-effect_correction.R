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
pca_data_path <- cl_args[5]
seu.obj.pca <- readRDS(pca_data_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Batch effect correction --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# By default, the harmony API works on Seurats PCA cell embeddings and corrects them. You can run Harmony within your Seurat workflow with RunHarmony(). Prior RunHarmony() the PCA cell embeddings need to be precomputed through Seurat's API. For downstream analyses, use the harmony embeddings instead of pca.
# https://github.com/immunogenomics/harmony/

# Visualization note:
# With the high number of points, some points might be occluded (on the z-axis) by others.
# In particular, one group might occlude the other if points' order in the data frame is by group.
# https://stackoverflow.com/questions/15706281/controlling-the-order-of-points-in-ggplot2
# So we should probably randomize the order of rows in the data frame before plotting. 
# https://stackoverflow.com/a/11503439

# Visualize pre-BEC PCs.
# Embeddings = coordinates for each cell in low-dimensional space.
# https://satijalab.org/seurat/archive/v3.0/dim_reduction_vignette#:~:text=embeddings%3A%20stores%20the%20coordinates%20for,each%20dimension%20of%20the%20embedding
set.seed(random_seed)
if(correct_batch_effects) {
  for(batch_covariate in batch_covariates) {
    plot <- seu.obj.pca@reductions$pca@cell.embeddings %>% .[,1:2] %>% as.data.frame %>% 
      rownames_to_column(var = "updated_cell_id") %>% 
      dplyr::left_join(seu.obj.pca@meta.data %>% 
                         dplyr::select(updated_cell_id, !!as.name(batch_covariate))) %>% # https://stackoverflow.com/questions/68749491/dplyr-r-selecting-columns-whose-names-are-in-an-external-vector
      .[sample(nrow(.)),] %>% 
      ggplot(aes(x = PC_1, y = PC_2, color = !!as.name(batch_covariate))) %>%
      DSPplotScatter(title = paste0("Pre-harmonization | ", batch_covariate))
    
    DSPexportPlot(plot, path_out = output_dir_pubs, filename = paste0("pre-harmonization_", batch_covariate))
    
    rm(plot)
  }
}

seu.obj.harmonized <- seu.obj.pca
if(correct_batch_effects) {
  rm(seu.obj.harmonized)
  gc()
  
  # Run Harmony. 
  # Harmony can integrate over multiple covariates. To do this, specify a vector of covariates to integrate.
  system.time({
    seu.obj.harmonized <- RunHarmony(seu.obj.pca, batch_covariates)
  })
  # user  system elapsed 
  # 261.311  16.742 276.277 
  
  # Visualize harmonized data.
  for(batch_covariate in batch_covariates) {
    plot <- seu.obj.harmonized@reductions$harmony@cell.embeddings %>% .[,c("harmony_1", "harmony_2")] %>% as.data.frame %>% 
      rownames_to_column(var = "updated_cell_id") %>% 
      dplyr::left_join(seu.obj.harmonized@meta.data %>% dplyr::select(updated_cell_id, !!as.name(batch_covariate))) %>% 
      .[sample(nrow(.)),] %>% 
      ggplot(aes(x = harmony_1, y = harmony_2, color = !!as.name(batch_covariate))) %>%
      DSPplotScatter(title = paste0("Post-harmonization | ", batch_covariate), geom_jitter = TRUE)
    
    DSPexportPlot(plot, path_out = output_dir_pubs, filename = paste0("post-harmonization_", batch_covariate))
  }
  
  rm(plot)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk -------------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Export Seurat object.
saveRDS(seu.obj.harmonized, paste0(output_dir_rdata, "seuratObject_harmonized.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)