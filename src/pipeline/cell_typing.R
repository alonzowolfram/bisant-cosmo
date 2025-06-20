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
seu.obj.umap <- readRDS(clustered_data_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Cell typing (InSituType) --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Setup ----------------------------------------------------------------------------------------------------

# Required inputs to InSituType.
### 1) Count matrix ---- 
# cells (n rows) x genes (m columns). [I think they're supposed to be the raw counts--the vignette shows integer counts, so ...]
message("Extracting the count matrix.")
counts <- seu.obj.umap@assays$RNA$counts %>% as.matrix %>% t
dim(counts)
# [1] 405339   6182

### 2) Mean negative control values vector ----
# 1 for each cell (size n x 1).
message("Creating the mean negative control values vector.")
neg <- seu.obj.umap@assays$negprobes$counts %>% as.matrix %>% t %>% .[,!(colnames(.) %in% c("NegAdd"))] # We only want the host negative probes
dim(neg)
# [1] 405339     21
negmean <- Matrix::rowMeans(neg)

### 3) Per-cell expected background counts per gene vector ----
# (Optional). 
# If not provided, insitutype will estimate it from the negative controls.

### 4) Reference profile matrix ----
# genes (o rows) x cell types (p columns). 
# MUST BE IN LINEAR SCALE, NOT LOG. 
# If you suspect the tissue has cell types outside the reference matrix, use semi-supervised clustering. 
# Reference profiles will be read in from the .xlsx sheet inside the loop below. EDIT 2024/10/16: switching to Brenda's method. 
# We got our reference profiles from https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/master. 
final_ref_profiles <- list()

# For CosMx reference profiles:
message("Building reference profile matrix.")
final_ref_profile <- read_csv(ref_matrix_file) %>%
  column_to_rownames("...1") %>%
  as.matrix
final_ref_profiles[["IO"]] <- final_ref_profile

# For non-CosMx reference profiles:
ref_profile_names <- c("Lung_HCA", "Skin_HCA")
ref_profile_list <- list()
# https://rpubs.com/sediaz/createenvironment
for(ref_profile_name in ref_profile_names) {
  cell_ref_env <- new.env()
  load(paste0("ext/", ref_profile_name, ".RData"), envir = cell_ref_env)
  attach(cell_ref_env)
  ref_profile_list[[ref_profile_name]] <- profile_matrix %>% as.matrix()
  mode(ref_profile_list[[ref_profile_name]]) <- "numeric"
  rm(list=ls(cell_ref_env), envir = cell_ref_env)
  detach(cell_ref_env)
  cell_ref_env <- NULL
}

# Align genes and normalize the scale.
shared_genes <- Reduce(intersect, lapply(ref_profile_list, rownames))
ref_profile_list_corrected <- list()
for(ref_profile_name in names(ref_profile_list)) {
  ref_profile_list_corrected[[ref_profile_name]] <- ref_profile_list[[ref_profile_name]][shared_genes,] %>% magrittr::divide_by(quantile(., 0.99)) %>% magrittr::multiply_by(1000)
}

# # Omit immune cells from profiles.
# colname_list <- lapply(ref_profile_list_corrected, colnames)
# omit_from <- list()
# omit_from[["Spleen_HCA"]] <- c("innate.lymphoid.cells", "CD34+.progenitor", "platelet")
# omit_from[["ImmuneTumor_safeTME"]] <- c("macrophages", "B.naive", "B.memory", "plasma", "T.CD4.naive", "T.CD4.memory", "T.CD8.memory", "T.CD8.naive", "NK", "pDCs", "mDCs", "monocytes.C", "monocytes.NC.I", "Treg", "endothelial.cells", "fibroblasts")
# 
# # Merge profiles.
# ref <- c()
# for(ref_profile_name in names(ref_profile_list_corrected)) {
#   ref <- cbind(ref,
#                ref_profile_list_corrected[[ref_profile_name]][, setdiff(colname_list[[ref_profile_name]], omit_from[[ref_profile_name]])]
#   )
# }
# UPDATE 2025-06-04: We will no longer remove immune cells from non-IO profiles nor merge profiles.
# Instead, we will run cell typing separately for each profile set
# and make final decisions about cell types after that. 
for(ref_profile_name in ref_profile_names) {
  ref <- ref_profile_list_corrected[[ref_profile_name]]
  
  # Update profiles for CosMx.
  system.time({
    astats <- get_anchor_stats(
      counts = counts,
      neg = negmean,
      profiles = ref
    )
  })
  # user  system elapsed
  # 36.752   1.591  35.903
  
  # Estimate per-cell bg as a fraction of total counts.
  negmean_per_totcount <- mean(negmean) / mean(rowSums(counts))
  per_cell_bg <- rowSums(counts) * negmean_per_totcount
  
  # Now choose anchors.
  system.time({
    anchors <- choose_anchors_from_stats(
      counts = counts,
      neg = negmean,
      bg = per_cell_bg,
      anchorstats = astats,
      n_cells = 500,
      min_cosine = 0.3,
      min_scaled_llr = 0.01,
      insufficient_anchors_thresh = 20
    )
  })
  # The following cell types had too few anchors and so are being removed from consideration: mantle.B.cell
  # user  system elapsed
  # 2.527   0.617   3.163
  
  # Update profile
  system.time({
    updatedprofiles <- updateReferenceProfiles(
      reference_profiles = ref,
      counts = counts,
      neg = neg,
      bg = per_cell_bg,
      anchors = anchors
    )
  })
  # user  system elapsed
  # 0.207   0.004   0.211
  final_ref_profiles[[ref_profile_name]] <- updatedprofiles$updated_profiles
}

### 5) Cohorts ----
# Per code given to us by Claire Williams, we will cohort by fluorescence.
# Get the mean fluorescence.
biomarkers <- seu.obj.umap@meta.data %>% .[,base::grep('Mean\\.', colnames(.))]
# Perform automatic cohorting:
set.seed(random_seed)
system.time({
  cohort <- fastCohorting(biomarkers,
                          gaussian_transform = TRUE)
  # ("Gaussian_transform = TRUE" maps variables to gaussians in order to 
  #  place dramatically different variables on the same scale.)
})
# user  system elapsed
# 174.689   1.172 176.121
table(cohort)
seu.obj.umap@meta.data$fastCohorting <- paste0("cohort_", cohort) %>% as.factor

## Semi-supervised clustering ----------------------------------------------------------------------------------------------------
### Initial cell typing ----
message("Splitting data.")

# Split up data into chunks
split_matrix_rows <- function(mat, n) {
  rows <- nrow(mat)
  split_indices <- split(seq_len(rows), cut(seq_len(rows), n, labels = FALSE))
  lapply(split_indices, function(idx) mat[idx, , drop = FALSE])
}
system.time({count_chunks <- split_matrix_rows(counts, n_chunks)})
# Save chunks to scratch
purrr::walk2(count_chunks, seq_along(count_chunks), function(chunk, i) {
  saveRDS(chunk, file = glue::glue("{path_to_scratch}/counts_chunk_{i}.rds"))
})
# # Check chunk size (development only)
# chunk_sizes <- sapply(count_chunks, pryr::object_size)
# # Convert to a readable format
# chunk_sizes_mb <- as.numeric(chunk_sizes) / 1024^2
# # Print summary
# summary(chunk_sizes_mb)

# Set options for parallelization
plan(multisession, workers = 4)  # Adjust number of workers to match cores
options(future.globals.maxSize = 50 * 1024^3)  # Increase to 50 GiB

res_list <- list()
system.time({
for(ref_profile_name in c(ref_profile_names, "IO")) {
  # Parallelized
  insitutype_results <- future_map(1:n_chunks, function(i) {
    set.seed(random_seed)
    path_to_scratch <- path_to_scratch
    
    # Capture all output
    log_file_err <- here::here(output_dir_logs_ct, paste0("cell-typing_", ref_profile_name, "_", i, "_err.log"))
    log_file_out <- here::here(output_dir_logs_ct, paste0("cell-typing_", ref_profile_name, "_", i, "_out.log"))
    done_flag <- here::here(output_dir_logs_ct, paste0("done_", ref_profile_name, "_", i, ".txt"))
    # Write each worker's log to a separate file (i.e., serialize logging) to prevent deadlock from multiple workers writing to the same file.
    
    log_message <- paste0(Sys.time(), " - ", glue::glue("\nCell typing with reference profile {ref_profile_name}, chunk {i}\n"))
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    
    # Get the appropriate subsets of the current chunk
    counts_sub <- readRDS(glue::glue("{path_to_scratch}/counts_chunk_{i}.rds"))
    cells_sub <- rownames(counts_sub)
    negmean_sub <- negmean[names(negmean) %in% cells_sub]
    cohort_sub <- cohort[names(cohort) %in% cells_sub]
    
    result <- NULL
    if (identical(names(negmean_sub), rownames(counts_sub)) && 
        identical(names(cohort_sub), rownames(counts_sub))) {
      log_message <- glue::glue("\nPerforming InSituType clustering.")
      cat(log_message, file = log_file_out, append = TRUE)
      flush.console()
      
      insitutype_res <- tryCatch({
        # sink(log_file_out, append = TRUE) # Capture stdout (progress bar, messages)
        insitutype(
          x = counts_sub,
          neg = negmean_sub,
          cohort = cohort_sub,
          bg = NULL,
          n_clusts = n_clusts_lower_bound:n_clusts_upper_bound,
          reference_profiles = final_ref_profiles[[ref_profile_name]],
          update_reference_profiles = FALSE,
          n_phase1 = subsamp_size_1,
          n_phase2 = subsamp_size_2,
          n_phase3 = subsamp_size_3,
          n_starts = n_iterations,
          max_iters = max_iterations
        )
        # sink()
      }, # End insitutype()
      error = function(e) {
        log_message <- glue::glue("Error running InSituType: {e$message}")
        cat(log_message, file = log_file_err, append = TRUE)
        flush.console()
        warning(log_message)
        return(e)
      }
      ) # End tryCatch()
      
      rm(counts_sub, cells_sub, negmean_sub, cohort_sub)
      gc()
      
      # Return
      result <- insitutype_res$clust
      return(result)
      
    } else {
      return(result)  # Handle mismatches gracefully
    }
    
  }, .options = furrr_options(seed = TRUE)   # ensures reproducibility
  )
  # user   system  elapsed 
  # 166.834    9.225 4164.282 
  
  # Serialized
  # for(i in 1:length(count_chunks)) {
  #   counts_sub <- count_chunks[[i]]
  #   cells_sub <- rownames(counts_sub)
  #   negmean_sub <- negmean %>% .[names(.) %in% cells_sub]
  #   cohort_sub <- cohort %>% .[names(.) %in% cells_sub]
  #   
  #   if(identical(names(negmean_sub), rownames(counts_sub)) & identical(names(cohort_sub), rownames(counts_sub))) {
  #     message("Performing InSituType clustering.")
  #     # system.time({
  #     insitutype_res <- insitutype(
  #       x = counts_sub,
  #       neg = negmean_sub,
  #       cohort = cohort_sub,
  #       bg = NULL,
  #       n_clusts = n_clusts_lower_bound:n_clusts_upper_bound,
  #       reference_profiles = final_ref_profiles[[ref_profile_name]],
  #       update_reference_profiles = FALSE,
  #       n_phase1 = subsamp_size_1,   # 200 10000
  #       n_phase2 = subsamp_size_2,   # 500 20000
  #       n_phase3 = subsamp_size_3,   # 2000 1e+05
  #       n_starts = n_iterations,   # 1 10
  #       max_iters = max_iterations   # 5 40
  #     )
  #     # })
  #     # user   system  elapsed 
  #     # 8613.198   42.265 3141.775 
  #   }
  # }
  
  res_list[[ref_profile_name]] <- insitutype_results
  gc()
}
})

# For each different reference, bind the results together.
res_df_list <- list()
for(ref_profile_name in c(ref_profile_names, "IO")) {
   res_df_list[[ref_profile_name]] <- do.call(c, res_list[[ref_profile_name]]) %>% 
     as.data.frame %>% 
     dplyr::rename(!!paste0("InSituType_", ref_profile_name) := names(.)[1]) %>% 
     tibble::rownames_to_column(var = "updated_cell_id")
}

# Add to the metadata. 
seu.obj.celltyping <- seu.obj.umap
rm(seu.obj.umap); gc()
for(ref_profile_name in names(res_df_list)) {
  seu.obj.celltyping@meta.data <- seu.obj.celltyping@meta.data %>% 
    dplyr::left_join(res_df_list[[ref_profile_name]])
  # Restore row names (`left_join` removes row names)
  rownames(seu.obj.celltyping@meta.data) <- seu.obj.celltyping@meta.data$updated_cell_id # Because left_join removes row names.
  
}

# Add broad categories.
# Rules: 
# 1) Anything called a stromal or mast cell (blood.vessel.cell, mast.cell, lymph.vessel.cell) by Lung HCA will be that cell type.
# 2) Anything called a fibroblast (Mesenchymal.Fibroblast, Secretory.reticular.Fibroblast, Pro.inflammatory.Fibroblast, Secretory.papillary.Fibroblast, Pericytes.1, Pericytes.2) by Skin HCA will be that cell type.
# 3) Anything called an immune cell by IO will be that cell type. 
# 4) Conflicts will be resolved using the following priority: Lung HCA > Skin HCA > IO.
seu.obj.celltyping@meta.data <- seu.obj.celltyping@meta.data %>% 
  # Broad categories
  dplyr::mutate(
    InSituType_broad = case_when(
      # Rule 1: Lung_HCA overrides all
      InSituType_Lung_HCA %in% c("blood.vessel.cell", "mast.cell", "lymph.vessel.cell") ~ "Stromal and mast cells",
      
      # Rule 2: Skin_HCA second priority
      InSituType_Skin_HCA %in% c(
        "Mesenchymal.Fibroblast", "Secretory.reticular.Fibroblast", 
        "Pro.inflammatory.Fibroblast", "Secretory.papillary.Fibroblast", 
        "Pericytes.1", "Pericytes.2"
      ) ~ "Fibroblasts and pericytes",
      
      # Rule 3: IO fallback
      !(InSituType_IO %in% c("Endothelial", "Fibroblast")) & !base::grepl("^[a-z]$", InSituType_IO) ~ "Immune cells",
      
      # Rule 5: Everything else
      TRUE ~ "Parenchymal or unknown"
    )
  ) %>%
  # Final categories for non-parenchymal/unknown
  dplyr::mutate(
    final_cell_type = case_when(
      # Rule 1: Use Lung_HCA label if it matches one of the specified categories
      InSituType_Lung_HCA %in% c("blood.vessel.cell", "mast.cell", "lymph.vessel.cell") ~ InSituType_Lung_HCA,
      
      # Rule 2: Use Skin_HCA label if it matches one of the specified categories
      InSituType_Skin_HCA %in% c(
        "Mesenchymal.Fibroblast", "Secretory.reticular.Fibroblast", 
        "Pro.inflammatory.Fibroblast", "Secretory.papillary.Fibroblast", 
        "Pericytes.1", "Pericytes.2"
      ) ~ InSituType_Skin_HCA,
      
      # Rule 3: Use IO label if it's not Endothelial, not Fibroblast, and not a single letter
      !(InSituType_IO %in% c("Endothelial", "Fibroblast")) & !base::grepl("^[a-z]$", InSituType_IO) ~ InSituType_IO,
      
      # Rule 5: Catch-all
      TRUE ~ "Parenchymal or unknown"
    )
  )

# Count number of conflicts
df_conflict_check <- seu.obj.celltyping@meta.data %>%
  dplyr::mutate(
    rule_1 = InSituType_Lung_HCA %in% c("blood.vessel.cell", "mast.cell", "lymph.vessel.cell"),
    rule_2 = InSituType_Skin_HCA %in% c(
      "Mesenchymal.Fibroblast", "Secretory.reticular.Fibroblast", 
      "Pro.inflammatory.Fibroblast", "Secretory.papillary.Fibroblast", 
      "Pericytes.1", "Pericytes.2"
    ),
    rule_3 = !(InSituType_IO %in% c("Endothelial", "Fibroblast")) & !base::grepl("^[a-z]$", InSituType_IO),
    n_rules_matched = rule_1 + rule_2 + rule_3
  )
conflict_res <- df_conflict_check %>%
  count(n_rules_matched)
total_cells <- sum(conflict_res$n)
total_conflicts <- sum(conflict_res %>% dplyr::filter(n_rules_matched >= 2) %>% dplyr::select(n) %>% unlist)
conflict_percent <- round(total_conflicts / total_cells * 100, digits = 2)
message(glue::glue("Of {total_cells} cells, there were {total_conflicts} conflicts in cell typing results ({conflict_percent}% of total cells)"))

# Save.

### Unsupervised clustering of parenchymal/unknown cells ----
seu.obj.celltyping.subset <- subset(seu.obj.celltyping, subset = final_cell_type == "Parenchymal or unknown")
table(seu.obj.celltyping.subset$final_cell_type)

# Run Leiden clustering
system.time({seu.obj.celltyping.subset <- seu.obj.celltyping.subset %>% FindClusters(algorithm = 4, 
                                                                                     resolution = 0.001,
                                                                                     random.seed = random_seed, 
                                                                                     verbose = TRUE)})
# resolution = 0.05 results in 38 final clusters
# resolution = 0.01 results in 36 final clusters
# user   system  elapsed 
# 2540.673  194.725 2733.657 

### Combine all results ----
# 
leiden_clusters <- seu.obj.celltyping.subset@meta.data %>% dplyr::mutate(leiden_cluster = paste0("cluster_", seurat_clusters)) %>% dplyr::select(updated_cell_id, leiden_cluster)
seu.obj.celltyping@meta.data <- seu.obj.celltyping@meta.data %>% 
  dplyr::left_join(leiden_clusters) %>% 
  dplyr::mutate(leiden_cluster = replace_na(leiden_cluster, "")) %>% 
  dplyr::mutate(final_cell_type = base::gsub("Parenchymal or unknown", "", final_cell_type)) %>% 
  dplyr::mutate(InSituType_final = paste0(final_cell_type, leiden_cluster))

rm(seu.obj.celltyping.subset); gc()

## Cluster markers ----------------------------------------------------------------------------------------------------
Idents(seu.obj.celltyping) <- seu.obj.celltyping[["InSituType_clusters"]]$`InSituType_clusters`
seu.obj.celltyping.downsamp <- subset(seu.obj.celltyping, downsample = downsampling_celltyping, seed = random_seed)
system.time(insitutype_clust_markers <- FindAllMarkers(seu.obj.celltyping.downsamp, 
                                                       logfc.threshold = logfc_threshold,
                                                       only.pos = TRUE,
                                                       random.seed = rand_seed))
# user  system elapsed 
# 231.087 112.191 344.043 
insitutype_clust_markers_filt <- insitutype_clust_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup()
# dplyr::filter(avg_log2FC > 0.5) 
date <- format(Sys.time(), "%Y-%m-%d") # https://www.statology.org/r-date-format/

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Flight paths --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cols <- InSituType::colorCellTypes(
  freqs = table(insitutype_res$clust),
  palette = "brewers"
)

flightpath <- InSituType::flightpath_layout(
  logliks = insitutype_res$logliks,
  profiles = insitutype_res$profiles
)

pdf(here(output_dir_pubs, "flight-path_unsupervised.pdf"))
par(mar = c(0, 0, 0, 0))
plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[insitutype_res$clust])
text(flightpath$clustpos[, 1], flightpath$clustpos[, 2],
     rownames(flightpath$clustpos),
     cex = 0.7
)
dev.off()

# Heatmap
pdf(here(output_dir_pubs, "unsupervised-clustering_heatmap.pdf"), height = 50)
pheatmap::pheatmap(
  sweep(
    insitutype_res$profiles, 1,
    pmax(apply(insitutype_res$profiles, 1, max), 0.2), "/"
  )[apply(insitutype_res$profiles, 1, max) > 0.2, ],
  fontsize_row = 6,
  col = colorRampPalette(c("white", "darkblue"))(100)
)
dev.off()

# Plot in space
# Will need to do one plot for each slide
for(slide_name in unique(seu.obj.celltyping@meta.data[[slide_name_var]])) {
  meta_sub <- seu.obj.celltyping@meta.data %>% 
    dplyr::filter(!!as.name(slide_name_var) == slide_name)
  
  # Extract global XY coordinates from metadata
  xy <- as.matrix(meta_sub[, dimension_name_vars])
  
  pdf(here(output_dir_pubs, paste0("unsupervised-clustering_plot-in-space_", slide_name, ".pdf")))
  plot(xy, pch = 16, cex = 0.1, col = cols[meta_sub$InSituType_clusters], asp = 1)
  legend("topright", pch = 16, col = cols, legend = names(cols), cex = 0.5)
  dev.off()
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Save Seurat object.
saveRDS(seu.obj.celltyping, paste0(output_dir_rdata, "seuratObject_celltyping.rds"))

# Save cell-type markers.
write.csv(insitutype_clust_markers, paste0(output_dir_tabular, "InSituType-semisup-cluster-markers.csv"))
write.csv(insitutype_clust_markers_filt, paste0(output_dir_tabular, "InSituType-semisup-cluster-markers_top-10.csv"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# UMAP with clustering identities --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Change Identity
Idents(seu.obj.celltyping) <- "InSituType_clusters"

# Plot UMAP
pdf(here(output_dir_pubs, "UMAP_InSituType-clustering.pdf"),
    width = 16, height = 12
)
DimPlot(seu.obj.celltyping, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE)
dev.off()