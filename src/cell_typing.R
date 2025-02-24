#!/usr/bin/env Rscript
# https://training.nextflow.io/advanced/structure/#bin

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Setup -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Source the setup.R file.
# Note that we cannot use source() directly in Nextflow; see https://stackoverflow.com/questions/76067626/how-to-call-a-function-present-in-a-different-r-script-from-an-r-executable-usin
# and https://community.seqera.io/t/source-another-r-script-in-the-bin-directory/1059
# So we will have to use the workaround suggested above if the workflow system is Nextflow.
cl_args <- commandArgs(TRUE)
workflow_system <- cl_args[2]
if(workflow_system=="Nextflow") {
  path <- Sys.getenv("PATH") |> strsplit(":")
  bin_path <- tail(path[[1]], n=1)
  source(file.path(bin_path, "setup.R"))
} else {
  bin_path <- ""
  source("src/setup.R") # I guess the path it sources from is the current working directory, not the path the R script lives in.
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
neg <- seu.obj.umap@assays$negprobes$counts %>% as.matrix %>% t
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

# For CosMx reference profiles:
message("Building reference profile matrix.")
final_ref_profile <- read_csv(ref_matrix_file) %>%
  column_to_rownames("...1") %>%
  as.matrix

# For non-CosMx reference profiles:
# ref_profile_names <- c("Spleen_HCA", "ImmuneTumor_safeTME")
# ref_profile_list <- list()
# # https://rpubs.com/sediaz/createenvironment
# for(ref_profile_name in ref_profile_names) {
#   cell_ref_env <- new.env()
#   load(paste0(cell_type_ref_mats_path, ref_profile_name, ".RData"), envir = cell_ref_env)
#   attach(cell_ref_env)
#   ref_profile_list[[ref_profile_name]] <- profile_matrix %>% as.matrix()
#   mode(ref_profile_list[[ref_profile_name]]) <- "numeric"
#   rm(list=ls(cell_ref_env), envir = cell_ref_env)
#   detach(cell_ref_env)
#   cell_ref_env <- NULL
# }
# 
# # Align genes and normalize the scale.
# shared_genes <- Reduce(intersect, lapply(ref_profile_list, rownames))
# ref_profile_list_corrected <- list()
# for(ref_profile in names(ref_profile_list)) {
#   ref_profile_list_corrected[[ref_profile]] <- ref_profile_list[[ref_profile]][shared_genes,] %>% magrittr::divide_by(quantile(., 0.99)) %>% magrittr::multiply_by(1000)
# }
# 
# # Omit immune cells from profiles.
# colname_list <- lapply(ref_profile_list_corrected, colnames)
# omit_from <- list()
# omit_from[["Spleen_HCA"]] <- c("innate.lymphoid.cells", "CD34+.progenitor", "platelet")
# omit_from[["ImmuneTumor_safeTME"]] <- c("macrophages", "B.naive", "B.memory", "plasma", "T.CD4.naive", "T.CD4.memory", "T.CD8.memory", "T.CD8.naive", "NK", "pDCs", "mDCs", "monocytes.C", "monocytes.NC.I", "Treg", "endothelial.cells", "fibroblasts")
# 
# # Merge profiles.
# ref <- c()
# for(ref_profile in names(ref_profile_list_corrected)) {
#   ref <- cbind(ref,
#                ref_profile_list_corrected[[ref_profile]][, setdiff(colname_list[[ref_profile]], omit_from[[ref_profile]])]
#   )
# }
# 
# # Update profiles for CosMx.
# system.time({
#   astats <- get_anchor_stats(
#     counts = counts,
#     neg = negmean,
#     profiles = ref
#   )
# })
# # user  system elapsed 
# # 36.752   1.591  35.903 
# 
# # Estimate per-cell bg as a fraction of total counts.
# negmean_per_totcount <- mean(negmean) / mean(rowSums(counts))
# per_cell_bg <- rowSums(counts) * negmean_per_totcount
# 
# # Now choose anchors.
# system.time({
#   anchors <- choose_anchors_from_stats(
#     counts = counts,
#     neg = negmean,
#     bg = per_cell_bg,
#     anchorstats = astats,
#     n_cells = 500,
#     min_cosine = 0.3,
#     min_scaled_llr = 0.01,
#     insufficient_anchors_thresh = 20
#   )
# })
# # The following cell types had too few anchors and so are being removed from consideration: mantle.B.cell
# # user  system elapsed 
# # 2.527   0.617   3.163 
# 
# # Update profile
# system.time({
#   updatedprofiles <- updateReferenceProfiles(
#     reference_profiles = ref,
#     counts = counts,
#     neg = neg,
#     bg = per_cell_bg,
#     anchors = anchors
#   )
# })
# # user  system elapsed 
# # 0.207   0.004   0.211
# final_ref_profile <- updatedprofiles$updated_profiles

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
message("Performing InSituType clustering.")
set.seed(random_seed)
system.time({
  insitutype_res <- insitutype(
    x = counts,
    neg = negmean,
    cohort = cohort,
    bg = NULL,
    n_clusts = n_clusts_lower_bound:n_clusts_upper_bound,
    reference_profiles = final_ref_profile,
    update_reference_profiles = FALSE,
    n_phase1 = subsamp_size_1,   # 200 10000
    n_phase2 = subsamp_size_2,   # 500 20000
    n_phase3 = subsamp_size_3,   # 2000 1e+05
    n_starts = n_iterations,   # 1 10
    max_iters = max_iterations   # 5 40
  )
})
# user   system  elapsed 
# 8613.198   42.265 3141.775 

# Add to the metadata. 
seu.obj.celltyping <- seu.obj.umap
rm(seu.obj.umap)
seu.obj.celltyping@meta.data <- seu.obj.celltyping@meta.data %>% 
  dplyr::left_join(insitutype_res$clust %>% 
                     as.data.frame %>%
                     tibble::rownames_to_column(var = "updated_cell_id") %>%
                     dplyr::rename(InSituType_clusters = 2))

# Save.
rownames(seu.obj.celltyping@meta.data) <- seu.obj.celltyping@meta.data$updated_cell_id # Because left_join removes row names.

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
for(slide_name in unique(seu.obj.celltyping@meta.data$Slide)) {
  meta_sub <- seu.obj.celltyping@meta.data %>% 
    dplyr::filter(Slide == slide_name)
  
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