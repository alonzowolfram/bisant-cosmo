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
spatial_data_path <- cl_args[5]
seu.obj.spatial <- readRDS(spatial_data_path)

# Create a subdirectory within `output_dir_logs` to hold the logs for spatial analysis.
output_dir_logs_sa <- paste0(output_dir_logs, "spatial_analysis/")
if(!dir.exists(output_dir_logs_sa)) dir.create(output_dir_logs_sa)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Data preprocessing -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Data import ---- 
metadata <- seu.obj.spatial@meta.data
counts <- seu.obj.spatial@assays$RNA$counts
exprs <- seu.obj.spatial@assays$RNA$scale.data

## Spatial data preprocessing ----
# Originally this was to create a spatialTIME object, but now that we aren't using the spatialTIME package,
# we don't need to create the object. Some of the preprocessing steps will still come in handy, though.
spat_obj2_list <- list()

if(!flagVariable(coloc_cell_clust_var)) {
  # 1) Spatial list (list of data frames, each one corresponding to a different sample).
  spatial_data <- metadata %>% 
    dplyr::select(
      uniqueFOV,
      !!as.name(coloc_cell_clust_var),
      !!as.name(dimension_name_vars_sa[1]),
      !!as.name(dimension_name_vars_sa[2]) 
    ) %>%
    mutate(positive = 1) %>% # Variable that indicates whether the cell is positive for a certain feature (e.g. cell class. We will pivot this wider in a few lines.)
    arrange(!!as.name(coloc_cell_clust_var)) %>% 
    pivot_wider(names_from = coloc_cell_clust_var,
                values_from = "positive",
                values_fill = 0)
  spatial_list <- split(spatial_data, spatial_data$uniqueFOV)
  
  # 2) Summary table (each row = 1 sample (FOV, 1 entry in spatial_list)).
  cell_types <- metadata[[coloc_cell_clust_var]] %>% unique %>% replace_na("NA")
  summary_df <- lapply(spatial_list, function(spat) {
    # Collapse to the number of positive for each cell type and
    # percent of cells positive for each cell type.
    spat %>%
      group_by(uniqueFOV) %>%
      summarize(total_cells = n(),
                across(!!cell_types, sum)) %>% # Get cell numbers.
      mutate(across(!!cell_types, ~ .x / total_cells * 100, .names = "percent_{col}")) # Percent of total.
  }) %>%
    do.call(bind_rows, .)
  
  # 3) Clinical data.
  # Needs to have a column linking it to summary and spatial data. 
  # Here, uniqueFOV.
  clinical <- metadata %>% .[,which(colnames(.) %in% c("uniqueFOV", clin_data_cols))]
  
  # Create mif (multiplex immunofluorescence). 
  system.time({
    spat_obj <- create_mif(clinical_data = clinical, # Clinical data.
                           sample_data = summary_df, # Summary data.
                           spatial_list = spatial_list, # Spatial data.
                           patient_id = "uniqueFOV", # Links clinical and summary data.
                           sample_id = "uniqueFOV") # Links spatial and summary data. 
  })
  # user  system elapsed 
  # 0.238   0.026   0.265 
  # View the object.
  spat_obj
  
  # Decide a radius to measure around each cell.
  # Find range of x and y for appropriate function range. 
  ranges <- spatial_data %>% 
    group_by(uniqueFOV) %>% 
    summarize(xmin = min(!!as.name(dimension_name_vars_sa[1])),
              xmax = max(!!as.name(dimension_name_vars_sa[1])),
              ymin = min(!!as.name(dimension_name_vars_sa[2])),
              ymax = max(!!as.name(dimension_name_vars_sa[2])),
              xrange = max(!!as.name(dimension_name_vars_sa[1])) - min(!!as.name(dimension_name_vars_sa[1])),
              yrange = max(!!as.name(dimension_name_vars_sa[2])) - min(!!as.name(dimension_name_vars_sa[2])))
  
  # Create a ppp (point pattern) object for each uniqueFOV.
  ppp_list <- list()
  for(unique_fov in names(spatial_list)) {
    message(paste0("Creating ppp object for FOV ", unique_fov))
    x_coords <- spatial_list[[unique_fov]]$x_adj
    y_coords <- spatial_list[[unique_fov]]$y_adj
    xrange <- ranges %>% dplyr::filter(uniqueFOV == unique_fov) %>% .[,c("xmin", "xmax")] %>% unlist
    yrange <- ranges %>% dplyr::filter(uniqueFOV == unique_fov) %>% .[,c("ymin", "ymax")] %>% unlist
    marks <- metadata %>% dplyr::filter(uniqueFOV == unique_fov) %>% .[[coloc_cell_clust_var]] %>% unlist
    ppp_list[[unique_fov]] <- ppp(x = x_coords, y = y_coords, window = owin(xrange, yrange), marks = as.factor(marks))
  }
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Spatial analysis --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## ................................................
##
## Cellular niche identification ----
##
## ................................................
# https://www.youtube.com/watch?v=tKNcAWK1uHg&t=1809s - A good conceptual overview.
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#overview
# https://www.sc-best-practices.org/spatial/neighborhood.html
# https://divingintogeneticsandgenomics.com/post/neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data-in-seurat-and-bioconductor/ 
# https://divingintogeneticsandgenomics.com/post/how-to-do-neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data/ (uses radius to determine neighbors)
# https://www.bioconductor.org/packages/release/bioc/vignettes/hoodscanR/inst/doc/Quick_start.html (uses k nearest cells to determine neighbors)
# From AtoMx:
# Identify distinct cellular neighborhood clusters based on cell type composition across tissue. 
# This module helps define the structural composition of a tissue automatically by looking for regional differences in cell type composition. 
# Structures can be repeated structures that are frequently found within a tissue but which are not contiguous (e.g. glomeruli in the kidney, germinal centers in the lymph node) or which are physically connected across a tissue (e.g. epithelial layer in the colon).
# Parameters:
# Method (radius, # neighboring cells)
# Radius
# Number of neighborhoods
# Input: 

# Get a matrix of xy coordinates for each cell.
coords <- seu.obj.spatial@meta.data %>% dplyr::select(!!as.name(dimension_name_vars_sa[1]), !!as.name(dimension_name_vars_sa[2]), updated_cell_id, !!as.name(niche_cell_clust_var))
all.equal(Cells(seu.obj.spatial), rownames(coords)) # Sanity check.

# (For a given FOV), for each cell, find its nearest neighbors (either by k nearest cells or cells within a given radius).
# # Set up parallel workers
plan(multicore, workers = (parallel::detectCores() - 1))
# plan(multisession, workers = 4)
options(future.globals.maxSize = 10 * 1024^3)  # Increase to 10 GiB
system.time({cluster_res_list <- future_map(unique(seu.obj.spatial@meta.data$uniqueFOV)[1:100], identify_spatial_niches)})
# For 100 FOVs:
# user  system elapsed 
# 255.903   6.647 263.324

# rbind() all the data frames in cluster_res_list and left_join (by cell_ID) the resulting data frame to the seu.obj.celltype metadata.
cluster_res_df <- cluster_res_list %>%
  purrr::map(1) %>%     # extract the first element (i.e., df1) from each sublist
  dplyr::bind_rows() %>%    # bind all df1s together
  dplyr::mutate(kmeans_cluster_corrected = kmeans_cluster %>% regexPipes::gsub("[_ ]", "-") %>% regexPipes::gsub("\\:", "\\."))
seu.obj.spatial@meta.data %>% left_join(cluster_res_df, by = "updated_cell_id")
# rbind() all the `average_cell_type_abundance_tscaled` data frames in cluster_res_list.
avg_ct_abund_mat <- do.call(rbind, lapply(cluster_res_list, `[[`, 2)) # average abundance (of _neighboring cells_) for each individual niche
nn_mat_all <- do.call(rbind, lapply(cluster_res_list, `[[`, 3)) # nearest neighbors for each individual cell
# View(avg_ct_abund_mat %>% as.matrix)
# View(nn_mat_all)

# Cluster all niches based on cell composition OR
# find an "eigen-niche" OR
# cluster all individual cells into super-niches based on their neighbors
# OK, let's not cluster all individual cells. I did 30000 cells and it still took forever. 
# Maybe we could still use Seurat's built-in Leiden/Louvain clustering, though?

# # Cluster individual cells:
# If we do this, we'll probably have to use `FindClusters()`, because regular knn is going to take forever.
# system.time({super_niche_res_2 <- cluster_super_niches(nn_mat_all, max_k = 20, max_iter = 100)})

# Cluster individual FOVs' niches:
system.time({super_niche_res_1 <- kmeans_dynamic(avg_ct_abund_mat, max_k = 15, max_iter = 100)})
# For 375 niches:
# user  system elapsed 
# 19.016   0.744  19.408 
supercluster_res_df <- super_niche_res_1$cluster_labels %>%
  tibble::enframe(name = "kmeans_cluster_corrected", value = "supercluster") %>%
  dplyr::mutate(supercluster = paste0("sc", supercluster))
# Left-join supercluster_res_df onto cluster_res_df (one with more entries is on left) to attach supercluster assignment to individual cells
cell_cluster_assignments <- cluster_res_df %>% dplyr::left_join(supercluster_res_df)

## ................................................
##
## Niche characterization ----
##
## ................................................

# Visualize cell-type composition of each supercluster
sn_obj <- CreateSeuratObject(counts = t(avg_ct_abund_mat)) # genes x cells
sn_obj_meta <- data.frame(cluster = rownames(avg_ct_abund_mat)) # For right now, metadata = cluster identifiers
rownames(sn_obj_meta) <- rownames(avg_ct_abund_mat) # This is necessary so `AddMetaData` will accept it as metadata
sn_obj <- sn_obj %>% AddMetaData(metadata = sn_obj_meta)
super_niches <- super_niche_res_1$cluster_labels %>% as.data.frame %>% tibble::rownames_to_column(var = "cluster")
colnames(super_niches)[2] <- "supercluster"
sn_obj@meta.data <- sn_obj@meta.data %>% dplyr::left_join(super_niches %>% dplyr::mutate(supercluster = paste0("sc", supercluster))) # Add the supercluster information to the object metadata
rownames(sn_obj@meta.data) <- sn_obj@meta.data$cluster # Because of that weird thing where left-joining removes rownames
Idents(sn_obj) <- sn_obj@meta.data$supercluster

average_niche_abundance <- Seurat::AverageExpression(
  sn_obj,
  assays = "RNA",
  layer = "counts",
  features = rownames(sn_obj),
  return.seurat = FALSE,
  group.by = "ident",
)
avg_niche_abund <- average_niche_abundance$RNA %>% t
# View(avg_niche_abund %>% as.matrix)
ComplexHeatmap::Heatmap(t(scale(t(average_niche_abundance$RNA))),
                        show_row_dend = FALSE,
                        show_column_dend = FALSE,
                        rect_gp = grid::gpar(type = "none"),
                        cell_fun = cell_fun,
                        col = col_fun,
                        column_names_rot = 45)

# Calculate and visualize cell-cell interactions for each super-cluster
immune_cells <- c("B cell", "Conventional dendritic cell", "Macrophage", "Mast cell", "Monocyte", "Neutrophil", "NK cell", "Plasma", "Plasmablast", "Plasmacytoid dendritic cell", "T cell CD4", "T cell CD8", "T cell regulatory")
for(super_niche in cell_cluster_assignments$supercluster %>% unique) {
  # Get all the individual niches that belong to this super-niche
  niches_sc_i <- cell_cluster_assignments %>% dplyr::filter(supercluster==super_niche) %>% .$kmeans_cluster_corrected %>% unique
  cells_sc_i <- cell_cluster_assignments %>% dplyr::filter(supercluster==super_niche) %>% .$updated_cell_id
  
  # Order the nodes for readability using the `immune_cells` vector above
  nn_mat_reordered <- nn_mat_all[cells_sc_i,order_nodes_custom(colnames(nn_mat), immune_cells)]
  
  # Calculate p-values for co-occurrence with permutation test
  cooc_perm_test_res <- cooccurrence_permutation_test(nn_mat_reordered, n_perm = 1000, seed = random_seed, n_cores = 4)
  # Build co-occurrence graph
  g <- build_cooccurrence_graph(nn_mat_reordered)
  # Filter edges by significance
  g_sig <- filter_significant_edges(g, cooc_perm_test_res, alpha = 0.05)
  # Compute node abundance
  g_sig <- compute_node_abundance(g_sig, nn_mat_reordered)
  # Graph
  # plot_cooccurrence_graph(g_sig)
  label_offset <- 1.15
  # 1. Get circular layout positions
  layout_df <- igraph::layout_in_circle(g_sig) %>% as.data.frame()
  colnames(layout_df) <- c("x", "y")
  layout_df$name <- V(g_sig)$name
  layout_df$abundance <- V(g_sig)$abundance
  # 2. Compute outer label positions (slightly beyond each node's radius)
  layout_df$label_x <- layout_df$x * label_offset
  layout_df$label_y <- layout_df$y * label_offset
  
  ggraph(g_sig, layout = "linear", circular = TRUE) +
    geom_edge_arc(aes(width = weight), alpha = 0.4, color = "gray40", strength = 1) +
    coord_fixed() +
    scale_edge_width(range = c(0.2, 3)) +
    geom_node_point(aes(size = abundance, color = size)) +
    scale_size_continuous(range = c(3, 12), name = "% abundance") +
    scale_color_gradientn(colors = c("blue", "red"), name = "Degree") +
    geom_text(data = layout_df, aes(x = label_x, y = label_y, label = name),
              size = 4, hjust = 0.5) +
    theme_void() +
    guides(edge_width = guide_legend("Normalized\nco-occurrence\nweight"))
}

## ................................................
##
## Cell-cell signaling ----
##
## ................................................


## ................................................
##
## Ligand-receptor interactions ----
##
## ................................................

if(!flagVariable(lr_interact_table_file)) {
  # Load the table if it exists.
  if(file.exists(lr_interact_table_file)) {
    # Load the ligand-receptor interaction table
    lr_interact_table <- read.csv(lr_interact_table_file, header = F)
    # Rename the columns
    colnames(lr_interact_table) <- c("receptor_id", "receptor", "ligand", "receptor_cell", "ligand_cell", "cell_type_col", "r")
    
    # Load receptor-TF (transcription factor) interactions from CellPhoneDB (CellSign module)
    # A key to the transcription_factor_input.csv file:
    # "receptor_id": contains either the name of a complex (see: "receptor_id" column in "tf_input") or a gene name - corresponding to the receptor
    # "TF_symbol": the gene name - corresponding to the TF
    # "Effect": set to -1 if the receptor has an inhibitory effect on the TF, and 1 otherwise
    
    tf_input <- read.csv("ext/transcription_factor_input.csv")
    complex_input <- read.csv("ext/complex_input.csv")
    gene_input <- read.csv("ext/gene_input.csv")
    protein_input <- read.csv("ext/protein_input.csv")
    
    # Add gene names to `complex_input` (right now, it only has Uniprot names)
    complex_input <- complex_input %>%
      left_join(gene_input %>% dplyr::select(uniprot, gene_name), by = c("uniprot_1" = "uniprot")) %>%
      rename(gene_name_1 = gene_name) %>%
      left_join(gene_input %>% dplyr::select(uniprot, gene_name), by = c("uniprot_2" = "uniprot")) %>%
      rename(gene_name_2 = gene_name) %>% 
      left_join(gene_input %>% dplyr::select(uniprot, gene_name), by = c("uniprot_3" = "uniprot")) %>%
      rename(gene_name_3 = gene_name) %>%
      left_join(gene_input %>% dplyr::select(uniprot, gene_name), by = c("uniprot_4" = "uniprot")) %>%
      rename(gene_name_4 = gene_name) %>% 
      left_join(gene_input %>% dplyr::select(uniprot, gene_name), by = c("uniprot_5" = "uniprot")) %>%
      rename(gene_name_5 = gene_name) %>% 
      relocate(gene_name_1, gene_name_2, gene_name_3, gene_name_4, gene_name_5, .before = transmembrane)
    
    # Load DoRothEA data
    if(species == "Homo sapiens") {
      data(dorothea_hs, package = "dorothea") # Human
      dorothea <- dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B", "C", "D"))
    } else {
      data(dorothea_mm, package = "dorothea") # Mouse
      dorothea <- dorothea_mm %>% dplyr::filter(confidence %in% c("A", "B", "C", "D"))
    }
    
    # Create gene sets for each TF.
    # One upregulated, one downregulated. 
    # Split into nested list
    tf_list <- dorothea %>%
      group_by(tf) %>%
      summarize(
        upregulated = list(target[mor == 1]),
        downregulated = list(target[mor == -1]),
        .groups = "drop"
      ) %>%
      # Now convert to named list of lists
      tidyr::nest(data = c(upregulated, downregulated)) %>%
      mutate(
        tf_targets = purrr::map(data, ~list(
          upregulated = .x$upregulated[[1]],
          downregulated = .x$downregulated[[1]]
        ))
      ) %>%
      select(tf, tf_targets) %>%
      tibble::deframe() # converts to named list, tf names as top level
    # Example: View targets for one TF
    tf_list[["TP53"]]
    
    # For each row in the LR table supplied by the user,
    # calculate the final ligand-receptor interaction score
    # 1) Calculate ligand-receptor colocalizations
    # 2) Calculate pathway activation
    # But first, expand the LR table by FOV (i.e., multiply each row by the # of FOVs)
    unique_fovs <- metadata$uniqueFOV %>% unique
    lr_interact_table_exp <- lr_interact_table %>%
      mutate(tmp = 1) %>%
      full_join(data.frame(fov = unique_fovs, tmp = 1), by = "tmp") %>%
      select(-c(tmp)) %>% as.data.frame()
    
    # Now we can apply() (technically, use `pmap()` from the purrr package) the function `lig_rec_coloc()` to each row in `lr_interact_table_exp`
    lr_scores_list  <- purrr::pmap(lr_interact_table_exp, lig_rec_coloc)
    lr_scores <- do.call(rbind, lr_scores_list)
    
    # Filter lr_scores to include only the unique receptor-receptor cell combinations
    # (Because we'll be calculating pathway scores for receptor cells only, and in 
    # `lr_scores`, you might have two rows with the same receptor cell, but different ligand cells)
    # And then `pmap()` the function `path_activ()` to each row in the filtered `lr_scores`
    lr_scores_unique_recept <- lr_scores %>% dplyr::select(receptor_id, receptor_cell, receptor) %>% dplyr::distinct()
    path_scores_list <- purrr::pmap(lr_scores_unique_recept, path_activ)
    path_scores <- do.call(rbind, path_scores_list)
    # Split the comma-separated values in the `auc_score` and `tf` columns
    auc_split <- strsplit(path_scores$auc_score, ",")
    tf_split <- strsplit(path_scores$tf, ",")
    # Generate a list of named vectors
    named_aucs <- mapply(function(auc, tf) {
      setNames(as.numeric(auc), paste0("auc_score_", tf))
    }, auc_split, tf_split, SIMPLIFY = FALSE)
    # Bind into a data frame
    expanded_scores <- do.call(rbind, lapply(named_aucs, as.data.frame.list))
    # Combine with original `path_scores`
    path_scores_exp <- cbind(path_scores, expanded_scores)
    
    # To get the final scores, multiply the ligand-receptor colocalization scores by the pathway scores.
    # First, we will have to left-join `path_scores_exp` to `lr_scores`
    lr_scores <- lr_scores %>% dplyr::mutate(receptor_cell_receptor = paste(receptor_cell, receptor, sep = ":"))
    path_scores_exp <- path_scores_exp %>% dplyr::mutate(receptor_cell_receptor = paste(receptor_cell, receptor, sep = ":"))
    final_scores <- lr_scores %>% dplyr::left_join(path_scores_exp %>% dplyr::select(receptor_cell_receptor, contains("auc_score_")))
    
  } # End `if()`: File path given by `lr_interact_table_file` exists
} # End `if()`: `lr_interact_table_file` is set

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Export to disk --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Save colocalization results
write.table(nng_derived, paste0(output_dir_rdata, "nearest-neighbor-G.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
saveRDS(spat_obj2_uni, paste0(output_dir_rdata, "univariate-NN-G_res.rds"))
# Save spatial niche results
# Save ligand-receptor interactions

# Update latest module completed
updateLatestModule(output_dir_rdata, current_module)
