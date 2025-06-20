# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# General helper functions --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' Function to offset coordinates from different slides
#'
#' @param df A data frame of single-cell data with at least a column for the x-coordinates, a column for the y-coordinates, and a column indicating which slide a cell is on.
#' @param slide_col The column in `df` indicating which slide a cell belongs to
#' @param x_col The column in `df` with the x-coordinates of the cells
#' @param y_col The column in `df` with the y-coordinates of the cells
#' @param buffer_factor An amount, given as a decimal, indicating amount of padding to add around a slide. For example, 1.2 would give a padding of 20%.
#'
#' @return
#' @export
#'
#' @examples
adjust_spatial_coordinates <- function(df, slide_col = "Run_Tissue_name", x_col = "CenterX_global_px", y_col = "CenterY_global_px", buffer_factor = 1.2) {
  # Get unique slide IDs
  slide_ids <- unique(df[[slide_col]])
  
  # Assign an index to each slide
  slide_index <- setNames(seq_along(slide_ids), slide_ids)
  
  # Get the number of slides, columns, and rows.
  n_slides <- length(slide_ids)
  n_cols <- ceiling(sqrt(n_slides))
  n_rows <- ceiling(n_slides / n_cols)
  
  # Get max dimensions of all spatial coords
  coords <- df[,dimension_name_vars]
  x_range <- range(coords[, 1])
  y_range <- range(coords[, 2])
  max_dims <- c(width = diff(x_range), height = diff(y_range))

  # Determine buffer based on the max width/height + margin
  buffer_x <- max(max_dims[1]) * buffer_factor
  buffer_y <- max(max_dims[2]) * buffer_factor
  
  # Apply a shift to each slide's coordinates
  idx <- slide_index[df[[slide_col]]]
  row_idx <- (idx - 1) %/% n_cols
  col_idx <- (idx - 1) %% n_cols
  
  x_shift <- col_idx * buffer_x
  y_shift <- row_idx * buffer_y
  
  df$x_adj <- df[[x_col]] + x_shift
  df$y_adj <- df[[y_col]] + y_shift
  
  return(df)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Spatial niches --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' Wrapper function to cluster a matrix using kmeans
#' Previously `cluster_super_niches`
#'
#' @param mat A matrix
#'
#' @return
#' @export
#'
#' @examples
kmeans_dynamic <- function(mat, niche_k_means = NULL, max_k = 15, n_pcs = 10, seed = 1026, max_iter = 100) {
  library(FactoMineR)
  library(factoextra)
  library(cluster)
  library(ggpubr)
  
  set.seed(seed)
  
  # 1. PCA
  pca_result <- PCA(mat, graph = FALSE, ncp = n_pcs)
  nn_pca <- pca_result$ind$coord[, 1:n_pcs]
  
  if(is.null(niche_k_means)) {
    # 2. Elbow (WSS)
    elbow_plot <- fviz_nbclust(nn_pca, kmeans, method = "wss", k.max = max_k) +
      ggtitle("Elbow Method")
    
    # 3. Silhouette
    silhouette_plot <- fviz_nbclust(nn_pca, kmeans, method = "silhouette", k.max = max_k) +
      ggtitle("Silhouette Method")
    
    # 4. Gap Statistic
    gap_stat <- clusGap(nn_pca, FUN = kmeans, nstart = 25, K.max = max_k, B = 50, iter.max = max_iter)
    gap_plot <- fviz_gap_stat(gap_stat) + ggtitle("Gap Statistic")
    
    # Choose optimal k from gap statistic
    best_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
    # Other methods to choose best k:
    # Picks the k that gives the highest gap statistic period. Not super helpful if the graph just increases monotonously(?)
    # best_k <- which.max(gap_stat$Tab[, "gap"]) 
    # Picks the k following Tibshirani (2001) - Is this even real? I just got it off ChatGTP lol
    # # Extract gap values and standard errors
    # gap_vals <- gap_stat$Tab[, "gap"]
    # gap_se <- gap_stat$Tab[, "SE.sim"]
    # 
    # # Apply the 1-standard-error rule: find the smallest k such that
    # # gap(k) >= gap(k+1) - SE(k+1)
    # best_k <- which(sapply(1:(length(gap_vals) - 1), function(k) {
    #   gap_vals[k] >= (gap_vals[k + 1] - gap_se[k + 1])
    # }))[1]  # take the first such k
    
    # 5. Final k-means
    km <- kmeans(nn_pca, centers = best_k, nstart = 25, iter.max = max_iter)
  } else {
    km <- kmeans(nn_pca, centers = niche_k_means, nstart = 25, iter.max = max_iter)
    best_k <- NA
    plots = list()
  }
  
  # 6. Output
  return(list(
    cluster_labels = km$cluster,
    best_k = best_k,
    plots = list(
      elbow = elbow_plot,
      silhouette = silhouette_plot,
      gap = gap_plot
    ),
    pca_result = pca_result
  ))
}

#' Function to create nearest-neighbors profile for a given FOV
#'
#' @param unique_fov
#'
#' @return
#' @export
#'
#' @examples
create_neighbor_profile <- function(unique_fov) {
  # Capture all output
  log_file_err <- here::here(output_dir_logs_sa, paste0("nearest-neighbor-profiles_err_", unique_fov, ".log"))
  log_file_out <- here::here(output_dir_logs_sa, paste0("nearest-neighbor-profiles_out_", unique_fov, ".log"))
  done_flag <- here::here(output_dir_logs_sa, paste0("done_nearest-neighbor-profiles_", unique_fov, ".txt"))
  # Write each worker's log to a separate file (i.e., serialize logging) to prevent deadlock from multiple workers writing to the same file.
  
  log_message <- paste0(Sys.time(), " - ", glue::glue("Creating nearest-neighbor profile for FOV {unique_fov}"))
  cat(log_message, file = log_file_out, append = TRUE)
  flush.console()
  
  fov_identifiers <- unique_fov %>% str_split(":") %>% unlist
  
  # List of nulls to return in case of error
  null_list <- list(nn = NULL, nn_id = NULL, nn_df = NULL, nn_count = NULL, nn_mat = NULL)
  # Error types for `tryCatch()`
  error_types <- c("simpleError", "error", "subscriptOutOfBoundsError") 
  
  # Get all the cells belonging only to that FOV
  cells_i <- metadata %>% dplyr::filter(uniqueFOV==unique_fov) %>% .$updated_cell_id
  coords_sub <- coords %>% .[.$updated_cell_id %in% cells_i,]
  
  ### 1) Nearest neighbors ----
  if(niche_neighbor_type=="radius") {
    system.time({
      nn <- tryCatch(dbscan::frNN(x = coords_sub[,1:2], eps = niche_neighbor_cutoff),
                     error = function(e) {
                       log_message <- glue::glue("Error running frNN for FOV {unique_fov}: {e$message} \n Skipping this FOV")
                       cat(log_message, file = log_file_err, append = TRUE)
                       flush.console()
                       warning(log_message)
                       return(e)
                     })
    })
    if(sum(error_types %in% class(nn)) > 0) return(null_list)
    # For 500 cells:
    # user  system elapsed 
    # 1.565   0.011   1.582 
    # nn$id is "a list of lists containing the indices of the cells that are within [the selected] radius."
    nn_id <- nn$id
  } else {
    # Default to k nearest neighbors.
    system.time({
      nn <- tryCatch(dbscan::kNN(x = coords_sub[,1:2], k = niche_neighbor_cutoff),
                     error = function(e) {
                       log_message <- glue::glue("Error running kNN for FOV {unique_fov}: {e$message} \n Skipping this FOV")
                       cat(log_message, file = log_file_err, append = TRUE)
                       flush.console()
                       warning(log_message)
                       return(e)
                     })
    })
    if(sum(error_types %in% class(nn)) > 0) return(null_list)
    # Here, nn$id is an n x m matrix, in which n = number of cells, and m = k, the number of neighbors. 
    # Each entry in the matrix is the index for the ith nearest neighbor for the cell of that row. 
    nn_id <- matrix_to_named_row_list(nn$id)
  }
  
  ### 2) Neighbor classifications ----
  # Identify the classifications of the neighbors found in step 1).
  # Create the count matrix. 
  nn_df <- nn_id %>%
    stack()
  # nn_df$values = the indices of the neighbors of the cells specified by nn_df$ind.
  # Add the cell type (cluster ID) to nn_df and convert to factor.
  cluster_ids <- coords_sub[[niche_cell_clust_var]] %>% unname()
  nn_df$cluster_id <- cluster_ids[nn_df$values] %>% as.factor
  # Get counts of each cell type.
  system.time({
    nn_count <- nn_df %>%
      group_by(ind) %>%
      dplyr::count(cluster_id, .drop = FALSE)
  })
  # user  system elapsed 
  # 0.123   0.000   0.124 
  
  # First, pivot wide (create a cell x cluster_id matrix.)
  nn_count_wide <- nn_count %>%
    tidyr::pivot_wider(names_from = cluster_id, values_from = n)
  nn_mat <- nn_count_wide[,-1] %>% as.matrix()
  rownames(nn_mat) <- nn_count_wide$ind
  # If any of the cell types are missing, add them to `nn_mat`, with 0s as entries.
  missing_cell_types <- setdiff(metadata[[niche_cell_clust_var]] %>% unique, colnames(nn_mat))
  mct_mat <- matrix(0, nrow = nrow(nn_mat), ncol = length(missing_cell_types))
  dimnames(mct_mat) <- list(rownames(nn_mat), missing_cell_types)
  nn_mat <- cbind(nn_mat, mct_mat) %>% .[,order(colnames(.))]
  
  ### 3) Return ----
  return(list(nn = nn, nn_id = nn_id, nn_df = nn_df, nn_count = nn_count, nn_mat = nn_mat))
}
  
#' Function to identify spatial niches for a given FOV
#'
#' @param nn_count
#'
#' @return
#' @export
#'
#' @examples
identify_spatial_niches <- function(nn_count) {
  # Picks up where `create_neighbor_profile()` left off
  log_message <- paste0(Sys.time(), " - ", glue::glue("Identifying spatial niches for FOV {unique_fov} | niche_k_means: {niche_k_means}"))
  
  # Set the null list to be returned in case of error
  null_list <- list(k_means_id = NULL, average_cell_type_abundance = NULL, nn_mat = NULL)
  error_types <- c("simpleError", "error", "subscriptOutOfBoundsError") 
  
  ### 1) Clustering ----
  # Using the neighbors' classifications as inputs, perform clustering on each of the cells.
  # First, pivot wide (create a cell x cluster_id matrix.)
  nn_count <- nn_count %>%
    tidyr::pivot_wider(names_from = cluster_id, values_from = n)
  nn_mat <- nn_count[,-1] %>% as.matrix()
  rownames(nn_mat) <- nn_count$ind
  # If any of the cell types are missing, add them to `nn_mat`, with 0s as entries.
  missing_cell_types <- setdiff(metadata[[niche_cell_clust_var]] %>% unique, colnames(nn_mat))
  mct_mat <- matrix(0, nrow = nrow(nn_mat), ncol = length(missing_cell_types))
  dimnames(mct_mat) <- list(rownames(nn_mat), missing_cell_types)
  nn_mat <- cbind(nn_mat, mct_mat) %>% .[,order(colnames(.))]
  # Perform k-means clustering. 
  set.seed(random_seed)
  k_means_res <- tryCatch(kmeans_dynamic(nn_mat, niche_k_means = niche_k_means), # kmeans(nn_mat, centers = niche_k_means),
                          error = function(e) {
                            log_message <- glue::glue("Error in clustering cells: {e$message} \n Skipping this FOV")
                            cat(log_message, file = log_file_err, append = TRUE)
                            flush.console()
                            warning(log_message)
                            return(e)
                          })  
  if(sum(error_types %in% class(k_means_res)) > 0) return(null_list)
  k_means_id <- k_means_res$cluster_labels %>%
    tibble::enframe(name = "updated_cell_id", value = "kmeans_cluster")
  metadata_sub <- metadata %>% 
    dplyr::filter(uniqueFOV==unique_fov) %>% 
    left_join(k_means_id, by = "updated_cell_id") %>% 
    dplyr::mutate(kmeans_cluster = paste(uniqueFOV %>% regexPipes::gsub("\\_", "\\-") %>% regexPipes::gsub("\\:", "."), kmeans_cluster, sep = "."))
  
  ### 2) Visualization ----
  #### Heatmap: cell type composition per niche ----
  # Calculate the average abundance of cell types per niche.
  nn_obj <- CreateSeuratObject(counts = t(nn_mat))
  if(all.equal(Cells(nn_obj), metadata_sub[["updated_cell_id"]])) nn_obj <- nn_obj %>% AddMetaData(metadata = metadata_sub)
  Idents(nn_obj) <- nn_obj@meta.data$kmeans_cluster
  
  average_cell_type_abundance <- Seurat::AverageExpression(
    nn_obj,
    assays = "RNA",
    layer = "counts",
    features = rownames(nn_obj),
    return.seurat = FALSE,
    group.by = "ident",
  )
  avg_ct_abund <- average_cell_type_abundance$RNA %>% t
  if(nrow(avg_ct_abund) < 2) rownames(avg_ct_abund) <- nn_obj@meta.data$kmeans_cluster[1] # Needed because if there's only 1 cluster (and therefore only 1 row in avg_ct_abund), the rownames drop out
  # ComplexHeatmap::Heatmap(t(scale(t(average_cell_type_abundance$RNA))),
  #                         show_row_dend = FALSE,
  #                         show_column_dend = FALSE,
  #                         rect_gp = grid::gpar(type = "none"),
  #                         cell_fun = cell_fun,
  #                         col = col_fun,
  #                         column_names_rot = 45)
  # ggsave() doesn't work for ComplexHeatmap, since CH doesn't output a plot exactly?
  
  #### Clustered dot plot: cell type composition per niche ----
  # scCustomize::Clustered_DotPlot(nn_obj, features = rownames(nn_obj), plot_km_elbow = FALSE)
  # ggsave(
  #   paste0("~/PRIME-TR/test/DSP-021_neighborhood-analysis-test_dotplot.png"),
  #   plot = last_plot(),
  #   width = 12,
  #   height = 9,
  #   units = "in",
  #   dpi = 300
  # )
  
  #### Spatial plot: niches + cell types ----
  # dat <- nn_obj@meta.data %>% dplyr::select(!!as.name(dimension_name_vars_sa[1]), !!as.name(dimension_name_vars_sa[2]), kmeans_cluster) %>% dplyr::mutate(kmeans_cluster = as.factor(as.character(kmeans_cluster)))
  # spatial_plot <- dat %>% ggplot(aes(x = !!as.name(dimension_name_vars_sa[1]), y = !!as.name(dimension_name_vars_sa[2]), color = kmeans_cluster)) + 
  #   geom_point() + 
  #   labs(title = paste0("Slide ", fov_identifiers[1], " - FOV ", fov_identifiers[2])) + 
  #   guides(color=guide_legend(title="Niche"))
  # ggsave(
  #   paste0("~/PRIME-TR/test/DSP-021_neighborhood-analysis-test_scatterplot.png"),
  #   plot = last_plot(),
  #   width = 9,
  #   height = 9,
  #   units = "in",
  #   dpi = 300
  # )
  
  ### 3) Cleanup ----
  rm(metadata_sub, nn_obj)
  gc()
  
  ### 4) Return ----
  k_means_id$kmeans_cluster <- paste0(unique_fov, ":", k_means_id$kmeans_cluster)
  return(list(k_means_id = k_means_id, average_cell_type_abundance = avg_ct_abund, nn_mat = nn_mat))
}

#' Function to order nodes
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
order_nodes_custom <- function(node_names, immune_cells) {
  # node_names <- colnames(nn_mat)
  
  single_letter <- node_names[nchar(node_names) == 1]
  immune <- setdiff(intersect(node_names, immune_cells), single_letter)
  other <- setdiff(node_names, c(single_letter, immune))
  
  ordered_nodes <- c(single_letter, immune, other)
  
  return(ordered_nodes)
  
  # nn_mat_reordered <- nn_mat[,ordered_nodes]
  # return(nn_mat_reordered)
}

#' Function to generate cell-cell interaction/neighbor graphs
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
build_cooccurrence_graph <- function(nn_mat, 
                                     normalize = TRUE, 
                                     diagonal_to_zero = FALSE,
                                     edge_weight_filter = -1) {
  if(normalize) {
    nn_mat <- nn_mat / ifelse(rowSums(nn_mat) == 0, 1, rowSums(nn_mat))
  }
  
  # Build co-occurrence matrix
  co_mat <- t(nn_mat) %*% nn_mat # See physical journal entry for 2025/06/19 for why the transpose of a matrix multiplied by that matrix gives a co-occurrence matrix
  if(diagonal_to_zero) diag(co_mat) <- 0 # If diagonal of the co-occurrence matrix is 0, we do not consider neighbors of the same type as the focus cell
  
  edge_df <- as.data.frame(as.table(co_mat)) %>%
    rename(from = Var1, to = Var2, weight = Freq) %>%
    filter(weight > edge_weight_filter)
  
  # 🔥 Explicitly define node order from column names
  node_order <- colnames(nn_mat)
  node_df <- data.frame(name = node_order)
  
  g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)
  
  V(g)$degree <- degree(g)
  V(g)$size <- V(g)$degree
  
  return(list(edge_df = edge_df, g = g))
}

#' Function to calculate significance of cooccurrences (parallelized)
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
cooccurrence_permutation_test <- function(nn_mat, 
                                          normalize = TRUE, 
                                          diagonal_to_zero = FALSE,
                                          n_perm = 1000, 
                                          seed = 123, 
                                          n_cores = parallel::detectCores() - 1) {
  set.seed(seed)
  
  cell_types <- colnames(nn_mat)
  
  # Optional: normalize each row to proportions (comment out if using counts)
  if (normalize) {
    nn_mat <- nn_mat / ifelse(rowSums(nn_mat) == 0, 1, rowSums(nn_mat)) # Prevents division by 0
  }
  observed <- t(nn_mat) %*% nn_mat
  
  if(diagonal_to_zero) diag(observed) <- 0  # Remove self-co-occurrence
  
  # Parallelize over permutations
  system.time({
    perm_results <- parallel::mclapply(1:n_perm, function(p) {
      permuted_mat <- apply(nn_mat, 2, sample)
      null_dist <- t(permuted_mat) %*% permuted_mat
      if(diagonal_to_zero) diag(null_dist) <- 0
      null_dist
    }, mc.cores = n_cores)
  })

  # Convert list of matrices to 3D array
  null_distributions <- array(unlist(perm_results), dim = c(length(cell_types), length(cell_types), n_perm),
                              dimnames = list(cell_types, cell_types, NULL))
  
  # Compute p-values
  pval_mat <- matrix(NA, nrow = length(cell_types), ncol = length(cell_types),
                     dimnames = list(cell_types, cell_types))
  for (i in seq_along(cell_types)) {
    for (j in seq_along(cell_types)) {
      obs <- observed[i, j]
      null_vals <- null_distributions[i, j, ]
      pval_mat[i, j] <- mean(null_vals >= obs)
    }
  }
  
  return(list(observed = observed, pvals = pval_mat))
}

#' Function to keep only significant edges
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
filter_significant_edges <- function(g, test_result, alpha = 0.05) {
  pval_mat <- test_result$pvals
  obs_mat <- test_result$observed
  
  edge_df <- igraph::as_data_frame(g, what = "edges")
  
  # Extract all vertices in the original graph in the current order
  node_names <- igraph::V(g)$name
  node_df <- data.frame(name = node_names)
  
  # Attach p-values and weights to edges
  edge_df$pval <- mapply(function(from, to) {
    pval_mat[from, to]
  }, edge_df$from, edge_df$to)
  
  edge_df$weight <- mapply(function(from, to) {
    obs_mat[from, to]
  }, edge_df$from, edge_df$to)
  
  # Keep only significant edges
  sig_edges <- edge_df %>% filter(pval < alpha)
  
  # Build new filtered graph, preserving vertex order
  g_sig <- graph_from_data_frame(sig_edges, directed = FALSE, vertices = node_df)
  
  # Compute node degrees
  V(g_sig)$degree <- degree(g_sig)
  V(g_sig)$size <- V(g_sig)$degree
  
  return(g_sig)
}

#' Function to compute abundance of each node
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
compute_node_abundance <- function(g, nn_mat) {
  total_cells <- nrow(nn_mat)
  abundances <- colSums(nn_mat) / sum(nn_mat) * 100
  V(g)$abundance <- abundances[V(g)$name]
  return(g)
}

#' Function to color edges by node type
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
color_edges_by_node_type <- function(g, target_celltypes = NULL, target_colors = NULL, default_color = "gray80") {
  stopifnot(igraph::is.igraph(g))
  
  # 1. Get endpoints of each edge
  edge_ends <- igraph::ends(g, igraph::E(g))  # 2-column matrix: [from, to]
  
  # 2. Get all unique node names
  node_names <- igraph::V(g)$name
  
  # 3. If no cell types provided, default to all nodes
  if (is.null(target_celltypes)) {
    target_celltypes <- sort(unique(node_names))
  }
  
  # 4. If no colors provided, auto-generate one per cell type
  if (is.null(target_colors)) {
    target_colors <- scales::hue_pal()(length(target_celltypes))
  }
  
  # 5. Named vector: cell type → color
  color_map <- setNames(target_colors, target_celltypes)
  
  # 6. Assign edge color based on whether either end is a target cell type
  edge_colors <- apply(edge_ends, 1, function(pair) {
    matches <- intersect(pair, target_celltypes)
    if (length(matches) > 0) {
      return(color_map[matches[1]])  # pick first match if both ends match
    } else {
      return(default_color)
    }
  })
  
  # 7. Assign to igraph object
  igraph::E(g)$color <- edge_colors
  
  return(g)
}

#' Function to plot co-occurrence graph
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
# plot_cooccurrence_graph <- function(g) {
#   ggraph(g, layout = "circle") +
#     geom_edge_arc(aes(width = weight), alpha = 0.4, color = "gray40") +
#     scale_edge_width(range = c(0.2, 3)) +
#     geom_node_point(aes(size = abundance, color = size)) +
#     scale_size_continuous(range = c(3, 12), name = "% abundance") +
#     scale_color_gradientn(colors = c("blue", "red"), name = "Degree") +
#     geom_node_text(aes(label = name), repel = TRUE, size = 4) +
#     theme_void() +
#     guides(edge_width = guide_legend("Co-occurrence\nweight"))
# }
plot_cooccurrence_graph <- function(g, 
                                    label_offset = 1.15, 
                                    base_curvature = 0.25, 
                                    highlight_nodes = NULL,
                                    highlight_nodes_color = NULL) {
  # 1. Generate layout matrix
  layout_mat <- igraph::layout_in_circle(g)
  
  # 2. Turn into a data frame and attach node names
  layout_df <- as.data.frame(layout_mat)
  colnames(layout_df) <- c("x", "y")
  layout_df$name <- V(g)$name  # this is the key!
  
  # 3. Reorder layout_df to exactly match vertex order
  layout_df <- layout_df[match(V(g)$name, layout_df$name), ]
  
  # 4. Add any node attributes you want to map (e.g., abundance)
  layout_df$abundance <- V(g)$abundance
  
  # 5. Compute offset label coordinates
  layout_df$label_x <- layout_df$x * label_offset
  layout_df$label_y <- layout_df$y * label_offset
  
  # If `highlight_nodes` is not NULL,
  # then mark edges containing those nodes.
  if(!is.null(highlight_nodes)) {
    g <- color_edges_by_node_type(g, target_celltypes = highlight_nodes, target_colors = highlight_nodes_color, default_color = "gray80") 
  }
  
  # Plot
  # https://ggraph.data-imaginist.com/articles/Layouts.html
  ggraph(g, layout = 'linear', circular = TRUE) +
    geom_edge_arc(aes(color = I(color), width = weight), alpha = 0.4, strength = 1) +
    coord_fixed() + 
    scale_edge_width(range = c(0.2, 3)) +
    geom_node_point(aes(size = abundance, color = size)) +
    scale_size_continuous(range = c(3, 12), name = "% abundance") +
    scale_color_gradientn(colors = c("blue", "red"), name = "Degree") +
    geom_node_text(aes(label = name), repel = T) +
    theme_void() +
    labs(title = a) + 
    guides(edge_width = guide_legend("Normalized\nco-occurrence\nweight"))
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Univariate nearest neighbors --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' Helper function to suggest the number of quadrats to use when assessing inhomogeneity.
#'
#' @param ppp A multitype point pattern object
#' @param target_points_per_quad
#' @param max_nx
#' @param quiet
#'
#' @return
#' @export
#'
#' @examples
suggest_quadrat_grid <- function(ppp, target_points_per_quad = c(5, 10), max_nx = 10, quiet = FALSE) {
  if (!inherits(ppp, "ppp")) stop("Input must be a 'ppp' object.")
  
  n_total <- npoints(ppp)
  window_size <- as.rectangle(ppp)
  xrange <- diff(window_size$xrange)
  yrange <- diff(window_size$yrange)
  aspect_ratio <- xrange / yrange
  
  suggestions <- list()
  for (nx in 2:max_nx) {
    ny <- round(nx / aspect_ratio)
    ny <- max(ny, 1)
    n_quads <- nx * ny
    pts_per_quad <- n_total / n_quads
    
    if (pts_per_quad >= min(target_points_per_quad) && pts_per_quad <= max(target_points_per_quad)) {
      suggestions[[length(suggestions) + 1]] <- list(
        nx = nx,
        ny = ny,
        quadrats = n_quads,
        points_per_quadrat = round(pts_per_quad, 2),
        aspect_ratio = round(aspect_ratio, 2)
      )
    }
  }
  
  if (length(suggestions) == 0) {
    message("⚠️ No grid found that meets target points per quadrat. Try adjusting `max_nx` or target range.")
    return(NULL)
  }
  
  suggestions_df <- do.call(rbind, lapply(suggestions, as.data.frame))
  
  if (!quiet) {
    cat("📦 Suggested Quadrat Grid Sizes:\n")
    print(suggestions_df)
  }
  
  return(suggestions_df)
}

#' Function to calculate measures of inhomogeneity for a point pattern.
#'
#' @param ppp A point pattern object
#' @param nx Number of quadrats in x-direction
#' @param ny Number of quadrats in y-direction
#' @param plot_k Whether or not to plot K vs Kinhom plot
#'
#' @return
#' @export
#'
#' @examples
assess_inhomogeneity <- function(focal_ppp, context_ppp = NULL, nx = 5, ny = 5, plot_k = TRUE) {
  if(is.null(focal_ppp))  {
    warning("   ⚠️ You must provide at least one ppp object, the focal_ppp.")
    return(NULL)
  }
    
  # library(spatstat.core) # Replaced by spatstat.explore and spatstat.model
  library(spatstat.geom)
  library(spatstat.explore)
  library(spatstat.model)
  
  cat("🔍 Assessing spatial inhomogeneity...\n\n")
  
  # 1. Quadrat test
  qt <- quadrat.test(focal_ppp, nx = nx, ny = ny)
  vmr <- qt$statistic / qt$parameter
  cat("📦 Quadrat Test:\n")
  cat(sprintf("- Chi-squared statistic = %.2f, df = %d, p = %.4f\n", qt$statistic, qt$parameter, qt$p.value))
  cat(sprintf("- Variance-to-Mean Ratio (VMR) = %.2f\n\n", vmr))
  
  # 2. Estimate intensity and compute Coefficient of Variation (CV)
  if(is.null(context_ppp)) {
    # Intensity estimate using only focal points
    lambda_hat <- density(focal_ppp, sigma = bw.diggle)
  } else {
    # Intensity estimate using ALL context points
    lambda_hat <- density(context_ppp, sigma = bw.diggle)
  }
  lambda_vals <- as.vector(lambda_hat$v)
  lambda_vals <- lambda_vals[!is.na(lambda_vals)]
  cv_lambda <- sd(lambda_vals) / mean(lambda_vals)
  
  cat("🌄 Intensity CV:\n")
  cat(sprintf("- Coefficient of Variation (CV) of λ(x,y) = %.3f\n\n", cv_lambda))
  
  # 3. Plot K vs Kinhom if desired
  if (plot_k) {
    cat("📈 Plotting K(r) vs Kinhom(r)...\n\n")
    K_reg <- Kest(ppp, correction = "border")
    K_inhom <- Kinhom(ppp, lambda = lambda_hat, correction = "border")
    
    plot(K_reg$r, K_reg$border, type = "l", lwd = 2, col = "black",
         ylim = range(K_reg$border, K_inhom$border), 
         xlab = "r", ylab = "K(r)", main = "K vs. Kinhom")
    lines(K_inhom$r, K_inhom$border, col = "blue", lwd = 2, lty = 2)
    legend("bottomright", legend = c("K(r)", "Kinhom(r)"),
           col = c("black", "blue"), lty = c(1, 2), lwd = 2)
  }
  
  invisible(list(
    quadrat_test = qt,
    vmr = vmr,
    cv_lambda = cv_lambda,
    lambda_estimate = lambda_hat
  ))
}

#' Function to calculate measures of inhomogeneity for a multitype point pattern, for three different use cases:.
#' 1) For cell type `i`: only cells of the same type are neighbors
#' 2) For cell type `i`: only cells of a different type (not `i`) are neighbors
#' 3) For cell type `i`: only cells of specified types are neighbors
#'
#' @param ppp_multitype A multitype point pattern object
#' @param types
#' @param neighbor_mode description
#' @param neighbor_dict description
#' @param nx Number of quadrats in x-direction
#' @param ny Number of quadrats in y-direction
#' @param plot_k Whether or not to plot K vs Kinhom plot
#'
#' @return
#' @export
#'
#' @examples
assess_inhomogeneity_multitype_flexible <- function(ppp_multitype,
                                                    types = NULL,
                                                    # neighbor_mode = c("same", "not_same", "custom"),
                                                    neighbor_dict = NULL,
                                                    target_points_per_quad = c(5, 10),
                                                    max_nx = 10,
                                                    nx = 5,
                                                    ny = 5,
                                                    plot_k = TRUE,
                                                    save_table = FALSE,
                                                    return_table = TRUE,
                                                    table_path = "inhomogeneity_summary_all_modes.csv") {
  library(knitr)
  library(kableExtra)
  
  if (!is.marked(ppp_multitype)) stop("Input must be a marked (multitype) `ppp` object.")
  
  types <- types %||% levels(marks(ppp_multitype))
  all_modes <- c("same", "not_same", "custom")
  # neighbor_mode <- match.arg(neighbor_mode)
  results <- list()
  
  for (mode in all_modes) {
    cat("🔄 Mode:", mode, "\n")
    
    for (type in types) {
      cat("🔍 Checking type:", type, "\n")
      
      # Select focal points (of type i)
      sub_ppp <- ppp_multitype[marks(ppp_multitype) == type]
      
      if (npoints(sub_ppp) < 10) {
        cat("   ⚠️ Too few points — skipping\n")
        next
      }
      
      # Select neighbors based on the chosen strategy
      neighbors <- switch(mode,
                          same = sub_ppp,
                          not_same = ppp_multitype[marks(ppp_multitype) != type],
                          custom = {
                            if (is.null(neighbor_dict) || is.null(neighbor_dict[[type]])) {
                              warning(sprintf("   ⚠️ Skipping %s: No custom neighbor types defined.", type))
                              next
                            }
                            target_types <- strsplit(neighbor_dict[[type]], ",\\s*")[[1]]
                            ppp_multitype[marks(ppp_multitype) %in% target_types]
                          }
                        )
      
      # Combine focal and neighbor points for spatial context
      combined_ppp <- superimpose(focal = sub_ppp, neighbors = neighbors, W = Window(ppp_multitype))
      marks(combined_ppp) <- factor(c(rep("focal", npoints(sub_ppp)), rep("neighbor", npoints(neighbors))))
      
      # Assess inhomogeneity only for the focal points
      cat(sprintf("📐 Using grid: nx = %d, ny = %d\n", grid["nx"], grid["ny"]))
      result <- assess_inhomogeneity(sub_ppp, nx = grid["nx"], ny = grid["ny"], plot_k = plot_k)
      
      results[[type]] <- data.frame(
        Cell_Type = type,
        N_Focal = npoints(sub_ppp),
        N_Neighbors = npoints(neighbors),
        VMR = round(result$vmr, 3),
        CV_lambda = round(result$cv_lambda, 3),
        Quadrat_p = signif(result$quadrat_test$p.value, 3),
        stringsAsFactors = FALSE
      )
      cat("\n")
    }
  }
  
  # Combine results into one summary data frame
  summary_df <- do.call(rbind, results)
  
  cat("\n📊 Auto-Grid Inhomogeneity Summary Table:\n")
  print(
    kable(summary_df, format = "markdown", caption = "Auto-Suggested Quadrat Summary by Cell Type") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
  )
  
  if (save_table) {
    write.csv(summary_df, file = table_path, row.names = FALSE)
    message("📁 Table saved to: ", table_path)
  }
  
  if (return_table) return(summary_df) else invisible(NULL)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Bivariate nearest neighbors --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' Utility function to safely compute Gcross between two types
#'
#' @param ppp_object A vector of n observations for variable X.
#' @param type_i An n x n matrix of weights.
#' @param type_j description
#' @param correction description
#'
#' @return
#' @export
#'
#' @examples
safe_Gcross <- function(ppp_object, type_i, type_j, correction = "border") {
  if (!inherits(ppp_object, "ppp")) stop("Input must be a 'ppp' object")
  if (is.null(marks(ppp_object))) stop("No marks found in ppp object")
  
  if (!is.factor(marks(ppp_object))) {
    marks(ppp_object) <- factor(marks(ppp_object))
  }
  
  if (!is.multitype(ppp_object)) {
    warning("Not a multitype pattern — at least two types required.")
    return(NULL)
  }
  
  types_present <- levels(marks(ppp_object))
  if (!(type_i %in% types_present && type_j %in% types_present)) {
    warning(sprintf("One or both types (%s, %s) not present in marks.", type_i, type_j))
    return(NULL)
  }
  
  tryCatch({
    Gcross(ppp_object, i = type_i, j = type_j, correction = correction)
  }, error = function(e) {
    warning(sprintf("Gcross failed for %s → %s: %s", type_i, type_j, e$message))
    return(NULL)
  })
}

#' Wrapper to compute Gcross for multiple FOVs and type pairs
#' 
#' @param ppp_list A list of ppp objects
#' @param type_pairs description
#' @param correction description
#'
#' @return
#' @export
#'
#' @examples
batch_Gcross_analysis <- function(ppp_list, type_pairs, correction = "border", estimator = "rs") {
  results <- list()
  
  for (fov_name in names(ppp_list)) {
    ppp <- ppp_list[[fov_name]]
    for (pair in type_pairs) {
      type_i <- pair[1]
      type_j <- pair[2]
      Gres <- safe_Gcross(ppp, type_i, type_j, correction = correction)
      
      if (!is.null(Gres) && estimator %in% names(Gres)) {
        df <- data.frame(
          FOV = fov_name,
          Type_i = type_i,
          Type_j = type_j,
          r = Gres$r,
          G = Gres[[estimator]],
          stringsAsFactors = FALSE
        )
        results[[paste(fov_name, type_i, type_j, sep = "_")]] <- df
      }
    }
  }
  
  if (length(results) == 0) {
    warning("No valid Gcross results generated.")
    return(NULL)
  }
  
  do.call(rbind, results)
}

# Utility function to safely compute Kcross between two types
safe_Kcross <- function(ppp_object, type_i, type_j, correction = "border") {
  if (!inherits(ppp_object, "ppp")) stop("Input must be a 'ppp' object")
  if (is.null(marks(ppp_object))) stop("No marks found in ppp object")
  
  if (!is.factor(marks(ppp_object))) {
    marks(ppp_object) <- factor(marks(ppp_object))
  }
  
  if (!is.multitype(ppp_object)) {
    warning("Not a multitype pattern — at least two types required.")
    return(NULL)
  }
  
  types_present <- levels(marks(ppp_object))
  if (!(type_i %in% types_present && type_j %in% types_present)) {
    warning(sprintf("One or both types (%s, %s) not present in marks.", type_i, type_j))
    return(NULL)
  }
  
  tryCatch({
    Kcross(ppp_object, i = type_i, j = type_j, correction = correction)
  }, error = function(e) {
    warning(sprintf("Kcross failed for %s → %s: %s", type_i, type_j, e$message))
    return(NULL)
  })
}

# Wrapper to compute Kcross for multiple FOVs and type pairs
batch_Kcross_analysis <- function(ppp_list, type_pairs, correction = "border") {
  results <- list()
  
  for (fov_name in names(ppp_list)) {
    ppp <- ppp_list[[fov_name]]
    for (pair in type_pairs) {
      type_i <- pair[1]
      type_j <- pair[2]
      Kres <- safe_Kcross(ppp, type_i, type_j, correction = correction)
      
      if (!is.null(Kres)) {
        df <- data.frame(
          FOV = fov_name,
          Type_i = type_i,
          Type_j = type_j,
          r = Kres$r,
          K = Kres$border,
          stringsAsFactors = FALSE
        )
        results[[paste(fov_name, type_i, type_j, sep = "_")]] <- df
      }
    }
  }
  
  if (length(results) == 0) {
    warning("No valid Kcross results generated.")
    return(NULL)
  }
  
  do.call(rbind, results)
}

# Summarize Gcross or Kcross results by computing AUC and peak
summarize_cross_curves <- function(df, value_col = "G") {
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  
  df %>% 
    group_by(FOV, Type_i, Type_j) %>%
    summarize(
      AUC = sum(!!sym(value_col)) * mean(diff(r)),
      Max = max(!!sym(value_col)),
      .groups = "drop"
    )
}

# Plotting function for cross statistics
plot_cross_curves <- function(df, value_col = "rs") {
  library(ggplot2)
  
  ggplot(df, aes(x = r, y = .data[[value_col]], color = FOV)) +
    geom_line() +
    facet_grid(Type_i ~ Type_j) +
    labs(x = "Distance r", y = value_col, title = paste(value_col, "cross-type curves")) +
    theme_minimal()
}

# Generate simulation envelopes for Gcross
generate_Gcross_envelopes <- function(ppp, i, j, nsim = 199, use_random_labeling = TRUE, correction = "border") {
  simulate_expr <- if (use_random_labeling) expression(rlabel(ppp)) else expression(rpoispp(lambda = intensity(ppp), win = Window(ppp)))
  envelope(ppp, fun = Gcross, i = i, j = j, correction = correction, nsim = nsim, simulate = simulate_expr)
}

# Compare two cross-type G/K curves using permutation of labels (paired AUC diff)
compare_cross_pair_auc <- function(ppp, i, j, k, value_col = "rs", n_perm = 999, correction = "border") {
  if (!is.factor(marks(ppp))) marks(ppp) <- factor(marks(ppp))
  
  G_ij <- Gcross(ppp, i = i, j = j, correction = correction)
  G_ik <- Gcross(ppp, i = i, j = k, correction = correction)
  
  r_vals <- G_ij$r
  delta_obs <- sum(G_ij[[value_col]] - G_ik[[value_col]]) * mean(diff(r_vals))
  
  all_jk <- which(marks(ppp) %in% c(j, k))
  original_marks <- marks(ppp)
  deltas <- numeric(n_perm)
  
  for (p in 1:n_perm) {
    perm_marks <- original_marks
    perm_marks[all_jk] <- sample(perm_marks[all_jk])
    marks(ppp) <- perm_marks
    
    Gp_ij <- Gcross(ppp, i = i, j = j, correction = correction)
    Gp_ik <- Gcross(ppp, i = i, j = k, correction = correction)
    
    deltas[p] <- sum(Gp_ij[[value_col]] - Gp_ik[[value_col]]) * mean(diff(r_vals))
  }
  
  pval <- mean(abs(deltas) >= abs(delta_obs))
  list(delta = delta_obs, pval = pval, deltas = deltas)
}

# Global permutation test across FOVs for one (i, j, k) comparison
global_permutation_test_deltas <- function(delta_vec, n_perm = 999) {
  delta_obs <- mean(delta_vec)
  permuted_means <- replicate(n_perm, {
    signs <- sample(c(-1, 1), length(delta_vec), replace = TRUE)
    mean(delta_vec * signs)
  })
  p_val <- mean(abs(permuted_means) >= abs(delta_obs))
  
  list(mean_delta = delta_obs,
       p_value = p_val,
       null_distribution = permuted_means)
}

# Batch wrapper across FOVs for compare_cross_pair_auc
batch_compare_cross_pair_auc <- function(ppp_list, comparison_list, value_col = "rs", n_perm = 999, correction = "border") {
  results <- list()
  
  for (fov_name in names(ppp_list)) {
    ppp <- ppp_list[[fov_name]]
    for (cmp in comparison_list) {
      type_i <- cmp[1]
      type_j <- cmp[2]
      type_k <- cmp[3]
      
      if (all(c(type_i, type_j, type_k) %in% levels(marks(ppp)))) {
        res <- try(compare_cross_pair_auc(ppp, i = type_i, j = type_j, k = type_k,
                                          value_col = value_col, n_perm = n_perm, correction = correction),
                   silent = TRUE)
        if (!inherits(res, "try-error")) {
          results[[paste(fov_name, type_i, type_j, type_k, sep = "_")]] <- data.frame(
            FOV = fov_name,
            Type_i = type_i,
            Type_j = type_j,
            Type_k = type_k,
            Delta = res$delta,
            P_value = res$pval,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(results) == 0) {
    warning("No valid comparisons returned.")
    return(NULL)
  }
  
  do.call(rbind, results)
}

# Volcano plot of AUC differences
plot_auc_volcano <- function(df, delta_col = "Delta", p_col = "P_value") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  library(ggplot2)
  library(dplyr)
  
  df <- df |> 
    mutate(log10p = -log10(.data[[p_col]]),
           sig = .data[[p_col]] < 0.05)
  
  ggplot(df, aes(x = .data[[delta_col]], y = log10p, color = sig)) +
    geom_point(size = 2, alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    facet_wrap(~FOV) +
    labs(x = expression(Delta[AUC]), y = expression(-log[10](p)), title = "Gcross Pairwise Comparison Volcano Plot") +
    scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue")) +
    theme_minimal()
}

# Multiple testing correction for batch comparison results
adjust_pvals <- function(df, p_col = "P_value", method = "fdr") {
  df$Adj_P <- p.adjust(df[[p_col]], method = method)
  return(df)
}

# Group summary (e.g., per patient from FOVs)
summarize_by_group <- function(df, group_map) {
  df$Group <- group_map[df$FOV]
  aggregate(cbind(Delta, P_value) ~ Group + Type_i + Type_j + Type_k, data = df, FUN = function(x) c(mean = mean(x), sd = sd(x)))
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Ligand-receptor interactions --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' Function to create data frame of A-B cell pairs for a given FOV, where A = ligand-expressing cell, B = receptor-expressing cell
#'
#' @param df a data frame, one row per cell, with the following columns: FOV, cell_ID, x_coord, y_coord, cell_type
#' @param ligand_cell the cell type (which should be found in cell type column of `df`) that should express the ligand
#' @param receptor_cell the cell type (which should be found in cell type column of `df`) that should express the receptor
#' @param r maximum radius between two cells for them to be considered interacting
#'
#' @return list: FOV, LR score, ligand, receptor, cell type A, cell type B
#' @export
#'
#' @examples
create_ab_pairs <- function(df, ligand_cell, receptor_cell, r) {
  # Get current FOV from df
  fov <- df$FOV[1]
  
  # Make sure cells `ligand_cell` and `receptor_cell` are present in the data.
  if(!all(c(ligand_cell, receptor_cell) %in% df$cell_type)) {
    message <- "Either the specified ligand-expressing cell or the specified receptor-expressing cell is missing from the data set. Skipping this ligand x receptor combination; returning NA for LR score."
    return(NULL)
  }
  # Make sure `r` is appropriate.
  if(!is.finite(r)) {
    message <- "Specified radius is not a finite number; returning NA for LR score."
    return(NULL)
  }
  
  # Remove all but the cells belonging to either the ligand-expressing (A) or receptor-expressing (B) cell types
  df_sub <- df %>% dplyr::filter(cell_type %in% c(ligand_cell, receptor_cell))
  
  # Create data frame of all A-B pairs, then filter to include only the A-B pairs that are within radius `r` of each other
  # Convert to data.table for speed
  dt <- data.table::as.data.table(df_sub)
  # Create all combinations of cell pairs
  pairs <- CJ(i = 1:nrow(dt), j = 1:nrow(dt))[i < j]
  # Join to get coordinates and IDs
  pairs <- pairs[
    , `:=`(
      cell_ID_i = dt[i]$cell_ID,
      x_i = dt[i]$x_coord,
      y_i = dt[i]$y_coord,
      cell_type_i = dt[i]$cell_type,
      cell_ID_j = dt[j]$cell_ID,
      x_j = dt[j]$x_coord,
      y_j = dt[j]$y_coord,
      cell_type_j = dt[j]$cell_type
    )
  ]
  # Calculate Euclidean distance
  pairs[, distance := sqrt((x_i - x_j)^2 + (y_i - y_j)^2)]
  # Filter by radius r
  ab_pairs <- pairs[distance <= r & ((cell_type_i==ligand_cell & cell_type_j==receptor_cell) | (cell_type_j==ligand_cell & cell_type_i==receptor_cell))]
  
  return(ab_pairs)
}

#' Function to calculate ligand-receptor colocalizations
#'
#' @param fov a unique identifier for the FOV
#' @param ligand_cell the cell type (which should be found in column `cell_type_col` of the `metadata` data frame in the global environment) that should express the ligand
#' @param receptor_cell the cell type (which should be found in column `cell_type_col` of the `metadata` data frame in the global environment) that should express the receptor
#' @param ligand gene name in `counts` corresponding to the ligand
#' @param receptor gene name in `counts` corresponding to the receptor
#' @param r maximum radius between two cells for them to be considered interacting
#'
#' @return list: FOV, LR score, ligand, receptor, cell type A, cell type B
#' @export
#'
#' @examples
lig_rec_coloc <- function(fov, ligand_cell, receptor_cell, ligand, receptor, receptor_id, cell_type_col, r) {
  # Create a list of NAs to be returned in case of failure
  na_list <- data.frame(
    FOV = fov,
    ligand = ligand,
    receptor = receptor,
    receptor_id = receptor_id,
    ligand_cell = NA,
    receptor_cell = NA,
    ligand_cell_type = ligand_cell,
    receptor_cell_type = receptor_cell,
    radius = r,
    lr_score = NA,
    lr_score_overall = NA,
    lr_score_overall_normalized = NA
  )
  
  # Subset the metadata to include only this FOV
  df <- metadata %>% .[.$uniqueFOV == fov,]
  
  # Split `ligand` and `receptor` into component parts, if they are multimers
  ligand_units <- ligand %>% str_split("\\+") %>% unlist()
  receptor_units <- receptor %>% str_split("\\+") %>% unlist()
  
  # Make sure both the ligand and receptor genes are present in the data.
  if(!all(c(ligand_units, receptor_units) %in% rownames(counts))) {
    message <- "Either the ligand, receptor, or both is missing from the counts data, in whole or part (multimers). Skipping this ligand x receptor combination; returning NA for LR score."
    return(na_list)
  }
  # Make sure cells `ligand_cell` and `receptor_cell` are present in the data.
  if(!all(c(ligand_cell, receptor_cell) %in% df[[cell_type_col]])) {
    message <- "Either the specified ligand-expressing cell or the specified receptor-expressing cell is missing from the data set. Skipping this ligand x receptor combination; returning NA for LR score."
    return(na_list)
  }
  # Make sure `r` is appropriate.
  if(!is.finite(r)) {
    message <- "Specified radius is not a finite number; returning NA for LR score."
    return(na_list)
  }
  
  # Subset counts (genes x cells)
  cts_sub <- counts %>% .[,colnames(.) %in% df$updated_cell_id] %>% as.matrix
  # rm(counts); gc()
  
  # Get the number of cells in this FOV
  n_cells_total <- nrow(df)
  
  # Remove all but the cells belonging to either the ligand-expressing (A) or receptor-expressing (B) cell types
  df_sub <- df %>% dplyr::filter(!!as.name(cell_type_col) %in% c(ligand_cell, receptor_cell))
  # Rename the columns.
  select_cols <- c("updated_cell_id", dimension_name_vars_sa[1], dimension_name_vars_sa[2], cell_type_col)
  df_sub <- df_sub %>% dplyr::select(all_of(select_cols))
  colnames(df_sub) <- c("cell_ID", "x_coord", "y_coord", "cell_type")
  
  ab_pairs <- create_ab_pairs(df_sub, ligand_cell, receptor_cell, r)
  if(is.null(ab_pairs) || nrow(ab_pairs) < 1) {
    message <- "Could not create A-B pairs; returning NA for LR score."
    return(na_list)
  }
  
  # For each pair that passes the distance cutoff, calculate the LR score = (counts of ligand L in A) x (counts of receptor R in B)
  lr_score_df_list <- sapply(1:nrow(ab_pairs), simplify = F, FUN = function(k) {
    cell_i <- ab_pairs[k,"cell_ID_i"] %>% unlist
    cell_j <- ab_pairs[k,"cell_ID_j"] %>% unlist
    
    cell_type_i <- ab_pairs[k,"cell_type_i"] %>% unlist
    cell_type_j <- ab_pairs[k,"cell_type_j"] %>% unlist
    ligand_cell_k <- ifelse(cell_type_i == ligand_cell, cell_i, cell_j)
    receptor_cell_k <- ifelse(cell_type_i == receptor_cell, cell_i, cell_j)
    ligand_cell_type_k <- ifelse(cell_type_i == ligand_cell, cell_type_i, cell_type_j)
    receptor_cell_type_k <- ifelse(cell_type_i == receptor_cell, cell_type_i, cell_type_j)
    
    # The final counts for the ligand/receptor will be the minimum counts out of all the subunits that make up the ligand/receptor.
    ligand_counts <- cts_sub[ligand_units, ligand_cell_k] %>% min
    receptor_counts <- cts_sub[receptor_units, receptor_cell_k] %>% min
    
    lr_score <- ligand_counts * receptor_counts
    
    return(data.frame(
      FOV = fov,
      ligand = ligand,
      receptor = receptor,
      receptor_id = receptor_id,
      ligand_cell = ligand_cell_k,
      receptor_cell = receptor_cell_k,
      ligand_cell_type = ligand_cell_type_k,
      receptor_cell_type = receptor_cell_type_k,
      radius = r,
      lr_score = lr_score
    ))
  })
  lr_score_df <- do.call(rbind, lr_score_df_list) %>% 
    dplyr::mutate(lr_score_overall = mean(lr_score), # Get the mean score for the whole FOV
                  lr_score_overall_normalized = lr_score_overall / n_cells_total) # Normalize LR score by dividing by total number of cells in FOV
  
  # Clean up
  rm(df, df_sub, cts_sub); gc()
  
  # Return a data frame containing the FOV, LR score, ligand, receptor, cell type A, and cell type B
  return(lr_score_df)
}

#' Function to calculate pathway activation
#'
#' @param fov a unique identifier for the FOV
#' @param ligand_cell the cell type (which should be found in column `cell_type_col` of the `metadata` data frame in the global environment) that should express the ligand
#' @param receptor_cell the cell type (which should be found in column `cell_type_col` of the `metadata` data frame in the global environment) that should express the receptor
#' @param ligand gene name in `counts` corresponding to the ligand
#' @param receptor gene name in `counts` corresponding to the receptor
#' @param r maximum radius between two cells for them to be considered interacting
#' @param gene_sets a named list in which each element is a character vector of genes
#'
#' @return list: FOV, LR score, ligand, receptor, cell type A, cell type B
#' @export
#'
#' @examples
path_activ <- function(receptor_id, receptor_cell, receptor, gene_sets) {
  # Create a list of NAs to be returned in case of failure
  na_list <- data.frame(
    receptor_cell = receptor_cell,
    receptor = receptor,
    receptor_id = receptor_id,
    auc_score = NA,
    tf = NA
  )
  
  # We really only need each individual receptor cell for this
  # since we're not really looking at any spatially dependent quantities
  
  # For this row, test all gene sets connected to receptor _i_ for receptor cell _i_
  # receptor -> TF symbol -> DoRothEA genes
  
  # Subset `exprs` matrix to include only the current cell
  exprs_sub <- exprs[,receptor_cell]
  # Get all the TFs associated with receptor complex `receptor_id` (`receptor_id` in `complex_input` corresponds to `receptor_id` in `tf_input`)
  tfs <- tf_input %>% .[.$receptor_id == receptor_id,] %>% .$TF_symbol
  
  # If there are no TFs, return the NA list
  if(length(tfs) < 1) {
    message("No TFs found for this receptor. Skipping")
    return(na_list)
  }
  # If none of the TFs have entries in tf_list (DoRothEA), return the NA list
  if(sum(tfs %in% names(tf_list)) < 1) {
    message("No target genes found for these TFs. Skipping")
    return(na_list)
  }
  
  # For each transcription factor in tfs,
  # get the appropriate upregulated and downregulated gene sets, and run AUC
  auc_scores <- c()
  for(tf in tfs) {
    if(!(tf %in% names(tf_list))) {
      auc_scores <- c(auc_scores, 0)
      next
    }
    
    genes_up <- tf_list[[tf]]$upregulated
    genes_down <- tf_list[[tf]]$downregulated
    
    # Run AUC
    if(length(genes_up) > 0) {
      auc_up_res <- tryCatch(AUCell::AUCell_run(exprMat = as.matrix(exprs_sub), geneSets = genes_up), # exprMat: genes x cells
                             error = function(e) {
                               warning(glue::glue("Error running AUCell for cell {receptor_cell} - receptor {receptor_id} ({receptor}) - TF {tf}: {e$message} \n Returning 0 for AUC score"))
                               return(0)
                             }
      )
      auc_up <- ifelse(class(auc_up_res)=="numeric", 0, auc_up_res@assays@data@listData$AUC %>% .[1,1])
    } else {
      auc_up <- 0
    }
    
    if(length(genes_down) > 0) {
      auc_down_res <- tryCatch(AUCell::AUCell_run(exprMat = as.matrix(exprs_sub), geneSets = genes_down), # exprMat: genes x cells
                             error = function(e) {
                               warning(glue::glue("Error running AUCell for cell {receptor_cell} - receptor {receptor_id} ({receptor}) - TF {tf}: {e$message} \n Returning 0 for AUC score"))
                               return(0)
                             }
      )
      auc_down <- ifelse(class(auc_down_res)=="numeric", 0, auc_down_res@assays@data@listData$AUC %>% .[1,1])
    } else {
      auc_down <- 0
    }
    
    # Combine AUC scores to get the final score
    norm_auc_up <- ifelse(length(genes_up) > 0, (auc_up / length(genes_up)), 0)
    norm_auc_down <- ifelse(length(genes_down) > 0, (auc_down / length(genes_down)), 0)
    auc_scores <- c(auc_scores, (norm_auc_up - norm_auc_down))
    
  } # End tf `for()` loop
  
  return(data.frame(
    receptor_cell = receptor_cell,
    receptor = receptor,
    receptor_id = receptor_id,
    auc_score = auc_scores %>% paste(collapse = ","),
    tf = tfs %>% paste(collapse = ",")
  ))
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Spatial colocalization --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#' Function to calculate univariate Moran's I.
#'
#' @param x A vector of n observations for variable X.
#' @param w An n x n matrix of weights.
#'
#' @return
#' @export
#'
#' @examples
uni_moran <- function(x, w) {
  # (xi-xbar)(xj-xbar) for all pairs
  dx <- x - mean(x)
  g <- expand.grid(dx, dx)
  xixj <- g[,1] * g[,2]
  
  # Make a matrix of the multiplied pairs.
  n <- length(x)
  pm <- matrix(xixj, ncol = n)
  
  # Multiply by the weights to be zero value for non-adjacent pairs.
  pmw <- pm * w
  
  # Sum the values.
  spmw <- sum(pmw)
  
  # Divide by the sum of the weights.
  smw <- sum(w)
  sw <- spmw / smw
  
  # Compute the inverse variance of x.
  vr <- n / sum(dx^2)
  
  # Last step to get Moran's I.
  MI <- vr * sw
  
  # Return the final value, as well as the z and spatially lagged x.
  return(MI)
}

# 
#' Function to calculate bivariate Moran's I.
#'
#' @param x A vector of n observations for variable X.
#' @param y A vector of n observations for variable Y.
#' @param w An n x n matrix of weights.
#'
#' @return
#' @export
#'
#' @examples
bi_moran <- function(x, y, w) {
  # (xi-xbar)(xj-xbar) for all pairs
  dx <- x - mean(x)
  dy <- y - mean(y)
  g <- expand.grid(dx, dy)
  xiyj <- g[,1] * g[,2]
  
  # Make a matrix of the multiplied pairs.
  n <- length(x)
  pm <- matrix(xiyj, ncol = n)
  
  # Multiply by the weights to be zero value for non-adjacent pairs.
  pmw <- pm * w
  
  # Sum the values.
  spmw <- sum(pmw)
  
  # Divide by the sum of the weights.
  smw <- sum(w)
  sw <- spmw / smw
  
  # Compute the inverse standard deviation of x and y.
  sd_inv_x <- sqrt(n / sum(dx^2))
  sd_inv_y <- sqrt(n / sum(dy^2))
  
  # Last step to get Moran's I.
  MI <- sd_inv_x * sd_inv_y * sw
  
  return(MI)
}

# 
#' Function to calculate univariate Anselin's LISA.
#'
#' @param x A vector of n observations for variable X.
#' @param y A vector of n observations for variable Y.
#' @param w An n x n matrix of weights.
#'
#' @return
#' @export
#'
#' @examples
anselin_lisa <- function(x, y = NULL, w) {
  if(is.null(y)) {
    # Univariate.
    
    # According to https://www.depts.ttu.edu/geospatial/center/gist4302/documents/lectures/Fall%202013/lecture6.pdf,
    # I_i = z_i * sum(w_ij * z_j)
    # and the summation, with index j, is across each row i of the spatial weights matrix, w.
    
    # z_i is the original variable x_i in standardized form OR deviation form.
    # We'll use deviation form for consistency.
    z <- x - mean(x)
    
    # To get sum(w_ij * z_j), this is just the dot product of w_i with z.
    # So to get I (the vector of all I_i), we can just left-multiply z by w (as long as w_ii = 0 for all i).
    I <- z * (w %*% z)
  } else {
    # Bivariate.
    zx <- x - mean(x)
    zy <- y - mean(y)
    
    I <- zx * (w %*% zy)
  }
  
  return(I)
  # A named vector. 
}