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

# Load the pre-DE objects created in the last step. 
pre_de_path <- cl_args[5]
pre_de <- readRDS(pre_de_path)
for(obj in names(pre_de)) {
  assign(obj, pre_de[[obj]], envir = .GlobalEnv)
}

# Create a subdirectory within `output_dir_logs` to hold the logs for DE.
output_dir_logs_de <- paste0(output_dir_logs, "differential_expression/")
if(!dir.exists(output_dir_logs_de)) dir.create(output_dir_logs_de)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Differential expression analysis -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## First-pass screening with independent-clusters DE ----
# Independent clusters passes k-means clusters as random effect to model.
system.time({
  modeldat <- smiDE::xy_kmeans_clusters(metainfo, 
                                        k_prop_n = k_prop_n, 
                                        split_neighbors_by_colname = slide_name_var) # split_neighbors_by_colname = "Run_Tissue_name"
  # k_prop_n = 0.05 (fewer clusters)
  # k_prop_n = 0.25 (more clusters)
})
# user  system elapsed 
# 86.862   0.629  87.283 

# Notes about the formula:
# 1) Offsets:
# "when fitting a negative binomial we need to specify the offset as log(totalcounts)" (from the smiDE vignette.)
# This essentially normalizes for sequencing depth or total transcript abundance by accounting for differences in total RNA capture across cells or regions.
# "An offset in a generalized linear model (GLM) is a term that adjusts for a known covariate without estimating an additional parameter."
# "log(totalcounts) serves as an offset to normalize for varying sequencing depth or transcript abundance across observations (cells, regions, or ROIs in CosMx)."
# "This ensures that the differential expression model accounts for differences in the total number of detected transcripts, preventing bias due to variations in sequencing depth." (ChatGPT)
# Q: Where do we get `totalcounts`?
# A: "totalcounts is the sum of all detected transcripts in a given observation (e.g., a cell, nucleus, or ROI).
#  "In a CosMx dataset, totalcounts can be computed for each unit (e.g., cell) by summing up the gene expression counts." (ChatGPT)
# So I guess that means it's just the `nCount_RNA` column in `metainfo`/`modeldat`.
# 2) 

# MOVING OUTSIDE `for()` loop.
# Add `otherct_expr` to `modeldat_sub`, a variable indicating the contamination coming from other cell types(?).
# I guess we'll just get it from `contam_metrics`?
# modeldat_sub <- modeldat_sub %>% dplyr::left_join(contam_metrics %>% 
#                                                     dplyr::filter(is.finite(ratio)) %>% 
#                                                     dplyr::group_by(cell_type) %>% 
#                                                     dplyr::summarise(otherct_expr = mean(ratio, na.rm = T))
# )
# Wait, it wouldn't make sense to get it from `contam_metrics`, because we're restricting DE analysis to a single cell type, and `ratio` in `contam_metrics` is averaged by cell type--so all the values would be the same.
# So, I guess we'll have to get it on a per-cell basis.
# OH WAIT ... we can get it from the pre-DE object?? 
# According to ChatGPT:
# neighbor_expr_byct is a list that contains the following matrices:
#   
# `allct`:
#   
# This matrix represents the total contamination from all neighboring cells, regardless of their cell type.
# Each row corresponds to a gene, and each column corresponds to a cell.
# The values represent the sum of gene expression from all neighboring cells within the defined radius (mm_radius).
# `otherct`:
#   
# This matrix represents the contamination from all neighboring cells that are NOT the same cell type as the target cell.
# Again, it is a gene × cell matrix, where each value indicates how much expression is coming from "other" cell types.
# Each specific cell type (e.g., Macrophage, Fibroblast, etc.):
#   
# These matrices represent contamination from only the specified cell type.
# For example, the Macrophage matrix contains only the expression contamination that comes from neighboring macrophages.
# The structure remains the same: genes × cells, where each value indicates how much gene expression in that cell is due to neighboring cells of that specific type.

# Add otherct.
otherct <- colSums(pre_de_obj$nblist$neighbor_expr_byct$otherct) %>% 
  as.data.frame %>% 
  tibble::rownames_to_column(var = "cell_ID") %>% 
  dplyr::rename(otherct_exp = 2)
# modeldat <- modeldat %>% dplyr::left_join(otherct)

# Set up parallel backend (80 cores)
plan(multicore, workers = (parallel::detectCores() - 1))

# # When running in RStudio:
# # Warning message:
# #   In supportsMulticoreAndRStudio(...) :
# #   [ONE-TIME WARNING] Forked processing ('multicore') is not supported when running R from RStudio because it is considered unstable. For more details, how to control forked processing or not, and how to silence this warning in future R sessions, see ?parallelly::supportsMulticore
# library(parallelly)
# # Override soft and hard CPU limits
# options(
#   parallelly.availableCores.override = parallel::detectCores(),  # Force R to recognize all cores
#   parallelly.makeNodePSOCK.setup_strategy = "sequential"  # Ensures stable setup
# )
# parallelly::availableCores()
# # Set up parallel workers
# plan(multisession, workers = 4)

# Run parallelized function
options(future.globals.maxSize = 10 * 1024^3)  # Increase to 10 GiB
system.time({
  de_list_nr <- future_map( # `nr` = no random spatial component
    ref_cell_types,
    function(cell_type_i) {
      library(smiDE)  # Ensure package is loaded in workers
      
      # SUBSET cts TO INCLUDE ONLY CELLS OF TYPE cell_type_i
      # This is necessary for CV; otherwise, other cell types would "contaminate" the CV results.
      cells <- modeldat %>% dplyr::filter(cell_type==cell_type_i) %>% .$cell_ID
      cts_sub <- cts %>% .[,cells] # Rows = genes; columns = cells
      
      # If CV cutoff provided, cull list of genes by cell type.
      set.seed(random_seed)
      CV_dat <- apply(cts_sub, MARGIN = 1, calc_CV) %>% na.omit() 
      
      if (is.null(cv_cutoff) || !is.finite(cv_cutoff) || !(0 <= cv_cutoff & cv_cutoff <= 1)) {
        top_cv_genes <- names(CV_dat)
      } else {
        top_cv_genes <- names(CV_dat[CV_dat > quantile(CV_dat, probs = cv_cutoff)])
      }
      # Make sure there are genes that are both in the top # of genes by CV and pass the contamination threshold.
      final_gene_list <- intersect(top_cv_genes, analysis_genes_list[[cell_type_i]])
      if(length(final_gene_list) < 1) {
        warning(glue::glue("No genes passing your supplied CV cutoff of {cv_cutoff} and having a contamination ratio of < 1 were found for cell type {cell_type_i}. Skipping."))
        return(NULL)
      }
      
      do.call(rbind, lapply(final_gene_list, function(gene_name) {
        process_gene_ic(cell_type_i, gene_name, modeldat, cts_sub)
      }))
    },
    .options = furrr_options(seed = random_seed)  # Ensures reproducibility
  )
})
# For two cell types:
# user  system elapsed 
# 31.115   4.058 905.067 

# Combine all results
de_results_nr <- do.call(rbind, de_list_nr)

# The following code is no longer necessary if we consolidate XY-coordinates from different runs to begin with.
# # Perform meta-analysis to combine results from different runs.
# meta_results_nr <- de_results_nr %>%
#   group_by(gene) %>%
#   summarize(
#     meta_estimate = metafor::rma(yi = "estimate", sei = "SE", method = "REML")$b,
#     meta_se = metafor::rma(yi = "estimate", sei = "SE", method = "REML")$se,
#     meta_pval = metafor::rma(yi = "estimate", sei = "SE", method = "REML")$pval
#   )

# Pick the top genes to move on to Gaussian Process (or psuedo-Gaussian Process) modeling.
de_results_nr <- de_results_nr %>% dplyr::mutate(wald_score = abs(estimate / SE))
top_genes <- de_results_nr %>% 
  dplyr::group_by(cell_type) %>% 
  dplyr::arrange(-wald_score) %>% 
  dplyr::slice_head(n = 100)
top_genes_list <- split(top_genes$gene, top_genes$cell_type)

## Gaussian process DE by Gaussian/pseudo-Gaussian modeling ----
# Set up parallel backend (80 cores)
plan(multicore, workers = (parallel::detectCores() - 1))
# plan(multisession, workers = 4)

# Run parallelized function
options(future.globals.maxSize = 10 * 1024^3)  # Increase to 10 GiB
system.time({
  de_list <- future_map(
    names(top_genes_list), function(cell_type_i) {
      library(smiDE)  # Ensure package is loaded in workers
      do.call(rbind, lapply(top_genes_list[[cell_type_i]], function(gene_name) {
        process_gene_gam(cell_type_i, gene_name, modeldat, cts)
      }))
    },
    .options = furrr_options(seed = random_seed)  # Ensures reproducibility
  )
})
# For two cell types, 100 genes each:
# user  system elapsed 
# 258.831 100.510 243.001 

# system.time({
#   de_list <- future_map(
#     names(top_genes_list), function(cell_type_i) {
#       library(smiDE)  # Ensure package is loaded in workers
#       do.call(rbind, lapply(top_genes_list[[cell_type_i]], function(gene_name) {
#         process_gene_test(cell_type_i, gene_name, modeldat, cts)
#       }))
#     },
#     .options = furrr_options(seed = random_seed)  # Ensures reproducibility
#   )
# })

# Combine all results
de_results <- do.call(rbind, de_list)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Save DE results.
# write.csv(df_metadata, here(output_dir_tabular, "metadata_Napari.csv"))
write.csv(de_results, here(output_dir_tabular, "DE_results.csv"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)
