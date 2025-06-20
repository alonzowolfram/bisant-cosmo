#' calcCellSignalingScoreWrapper
#'
#' @param fov_i 
#' @param anchors 
#' @param targets 
#' @param lr_pair 
#'
#' @return a data frame
#' @export
#'
#' @examples
calcCellSignalingScoreWrapper <- function(fov_i, anchors, targets, lr_pair) {
  # Split the ligand-receptor pair.
  lr_pair_split <- lr_pair %>% str_split("_") %>% unlist
  ligand <- lr_pair_split[1]
  receptor <- lr_pair_split[2]
  # Split the anchor cells.
  anchor <- anchors %>% str_split(",") %>% unlist
  # Split the target cells.
  target <- targets %>% str_split(",") %>% unlist
  
  # Subset genes to include only the current FOV.
  meta_sub <- meta %>% dplyr::filter(uniqueFOV==!!fov_i)
  cells_current_fov <- meta_sub %>% .$cell_id # !! interprets fov as the variable, not the column name.
  genes_sub <- genes[,cells_current_fov]
  # Convert to z-scores to normalize, then exponentiate to keep everything > 0.
  genes_sub_z <- apply(genes_sub, MARGIN = 1, FUN = zscore) %>% t %>% exp # t() because apply inverts it, apparently?
  
  # Get the cell IDs for the anchor and targets.
  anchor_cell_ids <- meta %>% dplyr::filter(!!as.name(coloc_cell_clust_var) %in% anchor & uniqueFOV == !!fov_i) %>% .$cell_id
  target_cell_ids <- meta %>% dplyr::filter(!!as.name(coloc_cell_clust_var) %in% target & uniqueFOV == !!fov_i) %>% .$cell_id
  if(length(anchor_cell_ids) < 1 | length(target_cell_ids) < 1) {
    message(paste0("There are either no anchor cells or no target cells for FOV ", fov_i))
    return(NULL)
  }
  
  # Find the nearest neighbors of each anchor cell and filter to include only the ones that are target cells.
  coords_sub <- meta_sub %>%
    dplyr::select(!!as.name(dimension_name_vars[1]), !!as.name(dimension_name_vars[2]), cell_id, !!as.name(coloc_cell_clust_var))
  if(niche_neighbor_type=="radius") {
    system.time({
      nn <- dbscan::frNN(x = coords_sub[,1:2], eps = niche_neighbor_radius)
    })
    # For 500 cells:
    # user  system elapsed 
    # 1.565   0.011   1.582 
    # nn$id is "a list of lists containing the indices [row numbers in the seu.obj.sub metadata] of the cells that are within [the selected] radius."
  } else {
    # Default to k nearest neighbors.
    system.time({
      nn <- dbscan::kNN(x = coords_sub[,1:2], k = niche_k_nearest)
    })
    # Here, nn$id is an n x m matrix, in which n = number of cells, and m = k, the number of neighbors. 
    # Each entry in the matrix is the index for the ith nearest neighbor for the cell of that row. 
  }
  nn_anchor <- nn$id %>% .[names(.) %in% anchor_cell_ids]
  # Get the indices from coords_sub of all target cells.
  target_indices <- which(coords_sub[[coloc_cell_clust_var]] %in% target)
  if(length(target_indices) < 1) {
    message(paste0("There are no anchor-target neighbor pairs for FOV ", fov_i))
    return(data.frame(cell_id = anchor_cell_ids, score = 0))
  }
  # Apply calcCellSignalingScore() over nn_anchor.
  scores <- lapply(names(nn_anchor), FUN = calcCellSignalingScore, genes_sub_z = genes_sub_z, target_indices = target_indices, ligand = ligand, receptor = receptor) %>% do.call(rbind, .)
  
  # Clean up.
  rm(genes_sub, genes_sub_z, fov_i, anchor, target, lr_pair, ligand, receptor, nn, nn_anchor, target_indices)
  
  # Return value.
  return(scores)
}
# user  system elapsed 
# 2.224   0.164   2.398 

#' calcCellSignalingScore
#'
#' @param anchor_id 
#' @param genes_sub_z 
#' @param target_indices 
#' @param ligand 
#' @param receptor 
#'
#' @return a data frame
#' @export
#'
#' @examples
calcCellSignalingScore <- function(anchor_id, genes_sub_z, target_indices, ligand, receptor, include_target_cells = FALSE) {
  # genes_sub_z = an n x m matrix where n = # of genes, m = # of cells in the FOV under consideration, and each entry is an exponentiated, z-score-normalized expression value.
  # anchor_id = the cell_id of anchor cell i.
  # neighbor_indices = the indices of the columns in gene_sub_z (which should be identical with the rows in coords_sub) corresponding to cells neighbor anchor cell i, given the chosen criteria for neighboring (radius or k-nearest).
  # target_indices = indices from coords_sub of all target cells.
  # ligand = the identity of the ligand in the ligand-receptor pair.
  # receptor = the identity of the receptor in the ligand-receptor pair.
  
  # Test values.
  # genes_sub_z
  # anchor_id <- "c_2_2_500"
  # neighbor_indices <- nn_anchor[[anchor_id]]
  # target_indices
  # ligand
  # receptor
  
  # The score = (exp{z-score} of receptor expression in anchor cell) * [sum(exp{z-scores} of ligand expression in target cells)] 
  # First see if any of the neighbor cells are target cells. If not, return 0 for the score. 
  neighbor_indices <- nn_anchor[[anchor_id]]
  target_neighbor_indices <- intersect(neighbor_indices, target_indices)
  target_neighbor_id <- coords_sub[target_neighbor_indices, "cell_id"]
  if(length(target_neighbor_indices) < 1) return(data.frame(cell_id = anchor_id, score = 0))
  
  # Otherwise ... 
  # Subset the vector of ligand expression to include only the values for the target cells, then sum.
  exprs_lig <- genes_sub_z %>% .[ligand,] %>% .[target_neighbor_indices] %>% sum
  exprs_rec <- genes_sub_z %>% .[receptor,] %>% .[names(.)==anchor_id]
  score <- exprs_lig * exprs_rec
  
  # Clean up.
  rm(genes_sub_z, exprs_lig, exprs_rec)
  
  # Return the score for both the anchor cell and all neighboring target cells (if include_target_cells = TRUE).
  if(include_target_cells) {
    cells <- c(anchor_id, target_neighbor_id)
  } else {
    cells <- c(anchor_id)
  }
  return(data.frame(cell_id = cells, score = rep(score, length(cells))))
}