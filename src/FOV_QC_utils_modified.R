
#### README:
# This script: functions for performing FOV QC of CosMx data. 
# See the vignette for details on how to run
# The runFOVQC function is the entry point
# Plotting functions: FOVEffectsHeatmap and FOVEffectsSpatialPlots
# Edits made by Lon Fong (2024/06/06):
# Changed gsub to base::gsub so as not to conflict with regexPipes::gsub


#' FOV QC 
#' 
#' Run the full FOV QC workflow: break each FOV into a 7x7 grid, compare grid squares 
#' to similar squares in other FOVs, and look for bardode bits whose genes are losing signal.
#' Recommended to use this function for one tissue/slide at a time, not across multiple tissues.
#' @param counts Raw counts matrix, cells in rows, genes in columns. Can be sparse or standard format.
#' @param xy 2-column matrix of cells' xy positions, aligned to rows of counts.
#' @param fov Vector of cells' FOV IDs, aligned to rows of counts.
#' @param barcodemap Data frame with two columns: "gene" and "barcode". Download the barcodemap for your panel
#' from https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/code/FOV%20QC.
#' @param max_prop_loss Maximum loss of efficiency allowed for any bit. E.g., a value of "0.3" means all FOVs with bias <log2(1 - 0.3) will be flagged.
#' @param max_totalcounts_loss Maximum loss of total expression allowed for any FOV. E.g., a value of "0.5" means all FOVs with total counts <50% of comparable spatial regions will be flagged.
#' @export
runFOVQC <- function(counts, xy, fov, barcodemap, max_prop_loss = 0.3, max_totalcounts_loss = 0.3) {
  
  if ((max_prop_loss > 1) | (max_prop_loss < 0)) {
    stop("max_prop_loss must fall in range of 0-1.")
  }
  fov <- as.character(fov)
  ## create a matrix of barcode bit expression over sub-FOV grids:
  # define grids, get per-square gene expression:
  gridinfo <- makeGrid(xy = xy, fov = fov, squares_per_fov = 49, min_cells_per_square = 10) 
  # convert to per-square barcode bit expression:
  bitcounts <- cellxgene2squarexbit(counts = counts, 
                                    grid = gridinfo$gridid, 
                                    genes = barcodemap$gene, 
                                    barcodes = barcodemap$barcode) 
  
  # normalize:
  bitcounts <- sweep(bitcounts, 1, rowSums(bitcounts), "/") * mean(rowSums(bitcounts))
  
  ## for every grid square, match it to "control" squares from other FOVs, and get its residuals from them:
  # get neighbors:
  comparators <- getNearestNeighborsByFOV(x = bitcounts, 
                                          gridfov = gridinfo$gridfov, 
                                          n_neighbors = 10)
  
  
  ## Search for FOVs with low outlier total counts:
  totalcounts <- Matrix::rowSums(counts)
  totalcountsgrid <- by(totalcounts, gridinfo$gridid, mean)[rownames(bitcounts)]
  # get expected total counts:
  totalcountshat <- sapply(1:length(totalcountsgrid), function(i) {
    mean(totalcountsgrid[comparators[i, ]], na.rm = TRUE)
  })
  # resids of total counts:
  totalcountsresids <- log2(totalcountsgrid / totalcountshat)
  # flag FOVs:
  flaggedgrids <- totalcountsresids < log2(1 - max_totalcounts_loss)
  flaggedgridsperfov <- by(as.vector(1*flaggedgrids), gridinfo$gridfov[names(flaggedgrids)], mean)
  flaggedfovs_fortotalcounts <- names(which(flaggedgridsperfov > 0.75))
  
  ## Search for failed reporter cycles:
  # get expected:
  yhat <- t(sapply(1:nrow(bitcounts), function(i) {
    colMeans(bitcounts[comparators[i, ], ], na.rm = TRUE)
  }))
  # get resids:
  resid <- log2((bitcounts + 1) / (yhat + 1))
  rownames(resid) = rownames(bitcounts)
  
  ## summarize bias per FOV * bit:
  fovstats <- summarizeFOVBias(resid = resid, gridfov = gridinfo$gridfov, max_prop_loss = max_prop_loss)
  
  # count flags per FOV * reportercycle:
  flags_per_fov_x_reportercycle <- c()
  reportercycle <- substr(colnames(fovstats$flag), 1, nchar(colnames(fovstats$flag)) - 1)
  for (rc in unique(reportercycle)) {
    flags_per_fov_x_reportercycle <- cbind(flags_per_fov_x_reportercycle, rowMeans(fovstats$flag[, reportercycle == rc]))
  }
  colnames(flags_per_fov_x_reportercycle) <- unique(reportercycle)
  
  # collate all flagged FOVs:
  flaggedfovs_forbias <- rownames(flags_per_fov_x_reportercycle)[rowSums(flags_per_fov_x_reportercycle >= 0.5) > 0]
  flaggedfovs <- union(flaggedfovs_fortotalcounts, flaggedfovs_forbias)
  
  # report on flagged FOVs:
  if (length(flaggedfovs) > 0) {
    message(paste0("The following FOVs failed QC for one or more barcode positions: ",
                   paste0(flaggedfovs, collapse = ", ")))
  } else {
    message("All FOVs passed QC")
  }
  
  # build a manifest of flagged gene/fov pairs:
  flagged_fov_x_gene <- c()
  for (f in flaggedfovs) {
    flaggedreportercycles <- names(which(flags_per_fov_x_reportercycle[f, ] >= 0.5))
    for (rc in flaggedreportercycles) {
      reporterposition <- as.numeric(substr(rc, 14, nchar(rc) - 1)) * 2
      flaggedgenes <- barcodemap$gene[substr(barcodemap$barcode, reporterposition, reporterposition) != "."]
      # add to growing list:
      tempflags <- cbind(rep(f, length(flaggedgenes)), flaggedgenes)
      colnames(tempflags) <- c("fov", "gene")
      flagged_fov_x_gene <- rbind(flagged_fov_x_gene, tempflags)
    }
  }
  
  return(list(flaggedfovs = flaggedfovs, 
              flaggedfovs_fortotalcounts = flaggedfovs_fortotalcounts, 
              flaggedfovs_forbias = flaggedfovs_forbias,
              flagged_fov_x_gene = flagged_fov_x_gene, 
              flags_per_fov_x_reportercycle = flags_per_fov_x_reportercycle ,
              fovstats = fovstats, resid = resid, totalcountsresids = totalcountsresids, 
              gridinfo = gridinfo, xy = xy, fov = fov))
}



#' Spatial plots of FOV effects:
#' 
#' @param res Results object created by runFOVQC
#' @param outdir Directory where png plots are printed
#' @param bits Which bits to plot. Defaults to "flagged_reportercycles" bits, but can also plot "all" or "flagged_bits".
#' @param plotwidth Width in inches of png plots
#' @param plotheight Height in inches of png plots
#' @return For each bit, draws a plot of estimated FOV effects
#' @export
FOVEffectsSpatialPlots <- function(res, outdir = NULL, bits = "flagged_reportercycles", plotwidth = NULL, plotheight = NULL) {
  
  if (is.null(plotwidth)) {
    plotwidth <- diff(range(res$xy[, 1])) * 1.5
  }
  if (is.null(plotheight)) {
    plotheight <- diff(range(res$xy[, 2])) * 1.5
  }
  bits_to_plot <- match(bits, colnames(res$fovstats$flag))
  if (bits == "flagged_reportercycles") {
    flaggedreportercycles <- colnames(res$flags_per_fov_x_reportercycle)[colSums(res$flags_per_fov_x_reportercycle >= 0.5) > 0]
    names_of_bits_to_plot  <- paste0(rep(flaggedreportercycles, each = 4), rep(c("B", "G", "R","Y"), length(flaggedreportercycles)))
    bits_to_plot <- match(names_of_bits_to_plot, colnames(res$resid))
  }
  if (bits == "flagged_bits") {
    bits_to_plot  <- which(colSums(res$fovstats$flag) > 0)
  }
  if (bits == "all") {
    bits_to_plot <- 1:ncol(res$resid)
  }
  temp <- sapply(bits_to_plot, function(i) {
    if (!is.null(outdir)) {
      png(paste0(outdir, "/", make.names(colnames(res$resid)[i]), ".png"), width = plotwidth, height = plotheight, units = "in", res = 300)
    }
    par(mar = c(0,0,2,0))
    plot(res$xy, cex = 0.2, asp = 1, pch = 16,
         col = colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)[
           pmax(pmin(51 + res$resid[match(res$gridinfo$gridid, rownames(res$resid)), i] * 50, 101), 1)], 
         main = paste0(colnames(res$resid)[i], ": log2(fold-change)\nfrom comparable regions elsewhere"))
    for (f in unique(res$fov)) {
      inds <- res$fov == f
      rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), border = "black")
    }
    for (f in rownames(res$fovstats$flag)[res$fovstats$flag[, i] > 0]) {
      inds <- res$fov == f
      rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), lwd = 2, border = "yellow")
    }
    legend("right", pch = 16,
           col = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
           legend = rev(c("< -1", -0.5, 0, 0.5, "> 1")))
    if (!is.null(outdir)) {
      dev.off()
    }
  })
}


#' Spatial plot of loss in signal strength compared to similar regions:
#' 
#' @param res Results object created by runFOVQC
#' @param shownames Logical for whether to display FOV names
#' @param outdir Directory where png plots are printed
#' @param plotwidth Width in inches of png plots
#' @param plotheight Height in inches of png plots
#' @return For each bit, draws a plot of estimated FOV effects
#' @export
FOVSignalLossSpatialPlot <- function(res, shownames = TRUE, outdir = NULL, plotwidth = NULL, plotheight = NULL) {
  
  if (!is.null(outdir)) {
    png(paste0(outdir, "/signal loss.png"), width = plotwidth, height = plotheight, units = "in", res = 300)
  }
  if (is.null(plotwidth)) {
    plotwidth <- diff(range(res$xy[, 1])) * 1.5
  }
  if (is.null(plotheight)) {
    plotheight <- diff(range(res$xy[, 2])) * 1.5
  }
  
  plot(res$xy, cex = 0.2, asp = 1, pch = 16,
       col = colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)[
         pmax(pmin(51 + res$totalcountsresids[match(res$gridinfo$gridid, names(res$totalcountsresids))] * 100, 101), 1)], 
       main = "Log2 fold-change in total counts compared to similar regions")
  for (f in unique(res$fov)) {
    inds <- res$fov == f
    rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), border = "black")
  }
  for (f in res$flaggedfovs_fortotalcounts) {
    inds <- res$fov == f
    rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), border = "yellow", lwd = 2)
    if (shownames) {
      text(median(range(xy[inds, 1])), median(range(xy[inds, 2])), f, col = "green")
    }
  }
  legend("right", pch = 16,
         col = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
         legend = rev(c("< -1", -0.5, 0, 0.5, "> 1")))
  if (!is.null(outdir)) {
    dev.off()
  }
}

#' Map of which FOVs were flagged:
#' 
#' @param res Results object created by runFOVQC
#' @param shownames Logical for whether to display FOV names
#' @param outdir Directory where png plots are printed
#' @param plotwidth Width in inches of png plots
#' @param plotheight Height in inches of png plots
#' @return For each bit, draws a plot of estimated FOV effects
#' @export
mapFlaggedFOVs <- function(res, shownames = TRUE, outdir = NULL, plotwidth = NULL, plotheight = NULL) {
  
  if (!is.null(outdir)) {
    png(paste0(outdir, "/flagged FOVs.png"), width = plotwidth, height = plotheight, units = "in", res = 300)
  }
  if (is.null(plotwidth)) {
    plotwidth <- diff(range(res$xy[, 1])) * 1.5
  }
  if (is.null(plotheight)) {
    plotheight <- diff(range(res$xy[, 2])) * 1.5
  }
  
  plot(res$xy, cex = 0.2, asp = 1, pch = 16,
       col = "black", 
       main = "Flagged FOVs")
  for (f in unique(res$fov)) {
    inds <- res$fov == f
    rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), col = scales::alpha("dodgerblue2", 0.25))
  }
  for (f in res$flaggedfovs) {
    inds <- res$fov == f
    rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), col = scales::alpha("red", 0.5))
    if (shownames) {
      text(median(range(xy[inds, 1])), median(range(xy[inds, 2])), f, col = "green")
    }
  }
  if (!is.null(outdir)) {
    dev.off()
  }
}


#' Heatmap of estimated bit bias across FOVs
#' 
#' @param res Results object created by runFOVQC
#' @return Draws a heatmap
#' @export
FOVEffectsHeatmap <- function(res) {
  pheatmap::pheatmap(res$fovstats$bias * (res$fovstats$flag),
                     col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
                     breaks = seq(-2,2,length.out = 101),
                     main = "FOV bias: log2(fold-change) from comparable regions in other FOVs")
}




#' Get nearest neighbors, sampling diffusely across other FOVs. 
#' 
#' @param x Data matrix, obs in rows, variables in columns
#' @param gridfov vector of FOV IDs, aligned to rows of x
#' @param n_neighbors How many neighbors to record
#' @return A matrix of dimension nrow(obs) * n_neighbors, giving the row IDs for each obs's selected neighbors.
getNearestNeighborsByFOV <- function(x, gridfov, n_neighbors = 10) {
  gridfov <- gridfov[rownames(x)]
  # get nearest:
  topneighbors <- FNN::get.knnx(data = x, 
                                query = x, 
                                k = min(n_neighbors * 50, nrow(x)))$nn.index
  # reduce neighbors: none from own FOV, max of max_per_fov from each other fov:
  reducedneighbors <- t(sapply(1:nrow(topneighbors), function(i) {
    tempneighbors <- topneighbors[i, ]
    # only accept 1 per FOV
    tempneighbors[duplicated(gridfov[topneighbors[i, ]])] <- NA
    # none from own FOV:
    tempneighbors[gridfov[topneighbors[i, ]] == gridfov[i]] <- NA
    return(tempneighbors[!is.na(tempneighbors)][1:n_neighbors])
  }))
  
  return(reducedneighbors)
}


#' Define a grid across FOVs
#' 
#' @param xy 2-column matrix of xy locations
#' @param fov Vector of FOV IDs
#' @param squares_per_fov Number of squares to break fov into
#' @param min_cells_per_square Fewest cells to keep a square
#' @return A list: a vector giving the grid square assignment of each cell,
#'             and a vector giving the FOV each grid square belongs to.   
makeGrid <- function(xy, fov, squares_per_fov = 49, min_cells_per_square = 25) {
  
  # make grids:
  ncuts <- floor(sqrt(squares_per_fov))
  grid <- rep(NA, nrow(xy))
  gridfov <- c()
  
  for (fovid in unique(fov)) {   
    inds <- fov == fovid
    grid[inds] <- paste0(fov[inds], "_", cut(xy[inds, 1], ncuts), "_", cut(xy[inds, 2], ncuts))
    # associate these grid IDs with this FOV:
    gridfov[unique(grid[inds])] <- rep(fovid, length(unique(grid[inds])))
  }
  
  # clean up low-data grids:
  toosmallgridids <- names(which(table(grid) < min_cells_per_square))
  grid[is.element(grid, toosmallgridids)] <- NA
  
  return(list(gridid = grid, gridfov = gridfov))
}


#' Convert cell x gene matrix to grid square * bit matrix
#' @param counts Counts matrix
#' @param grid Vector assigning cells to squares, aligned to rows of counts
#' @param genes Vector of gene names (omit negprobes and falsecodes)
#' @param barcodes Vector of gene barcodes, aligned to genes
#' @return A matrix of total gene expression in gridsquares x barcode bits
cellxgene2squarexbit <- function(counts, grid, genes, barcodes) {
  
  # number of bits:
  nreportercycles <- nchar(barcodes[1]) / 2
  nbits <- nreportercycles * 4
  # parse barcodes:
  bitmat = matrix(0, length(setdiff(unique(grid), NA)), nbits)
  colnames(bitmat) <- paste0("reportercycle", rep(seq_len(nreportercycles), each = 4), c("B", "G", 'Y', "R"))
  rownames(bitmat) <- setdiff(unique(grid), NA)
  for (i in seq_len(nreportercycles)) {
    barcodeposition <- i*2
    barcodehere <- substr(barcodes, barcodeposition, barcodeposition)
    for (col in c("B", "Y", "G", "R")) {
      tempgenes <- setdiff(genes[barcodehere == col], NA)
      tempgenes <- intersect(tempgenes, colnames(counts))
      temptotal <- Matrix::rowSums(counts[, tempgenes, drop = FALSE])
      tempsquaretotal <- by(temptotal, grid, mean)
      bitmat[names(tempsquaretotal), paste0("reportercycle", i, col)] <- tempsquaretotal
    }
  }
  return(bitmat)
}



#' Summarize bias in FOVs
#' 
#' @param resid Matrix of grid square x bit residuals
#' @param gridfov Vector giving the FOV ID each grid square (rows of resid) belong to.\
#' @param max_prop_change Maximum bias allowed. E.g., a value of "0.5" means all FOVs with bias >log2(1.5) or <log2(1/1.5) will be flagged.
#' For each bit, and each FOV, get 4 matrices:
#' 1. a logical matrix of flags (TRUE = flagged), 2. the estimated per-fov bias,
#' 3. p-values for those estimates, and 4. the proportion of each FOV's grid squares agreeing on the direction of bias.
summarizeFOVBias <- function(resid, gridfov, max_prop_loss) {
  gridfov = gridfov[rownames(resid)]
  fovs = unique(gridfov)
  bias <- p <- propagree <- matrix(NA, length(fovs), ncol(resid),
                                   dimnames = list(fovs, colnames(resid)))
  
  for (bit in colnames(resid)) {
    mod = summary(lm(resid[, bit] ~ as.factor(gridfov) - 1))$coef
    bias[base::gsub("as.factor\\(gridfov\\)", "", rownames(mod)), bit] = mod[, "Estimate"]
    p[base::gsub("as.factor\\(gridfov\\)", "", rownames(mod)), bit] = mod[, "Pr(>|t|)"]
    for (fov in fovs) {
      inds <- gridfov == fov
      #propagree[as.character(fov), bit] <- mean(sign(resid[inds, bit]) == median(sign(resid[inds, bit])))
      propagree[as.character(fov), bit] <- mean(resid[inds, bit] < log2(1 - max_prop_loss) / 3)
    }
  }
  
  # flagging rule:
  flag <- (bias < log2(1 - max_prop_loss)) * (p < 0.01) * (propagree >= 0.5)
  return(list(flag = flag, bias = bias, p = p, propagree = propagree))
}
