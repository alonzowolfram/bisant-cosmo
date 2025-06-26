# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Setup --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Loading CosMx QC app. Please wait while we get everything set up.")

## Required packages ----
library(shiny)
library(tidyverse)
library(data.table) # Faster than data frames.
library(DT) # Display tables.
library(bslib) # Styling of Shiny apps.
library(ggrastr)
library(viridisLite) # Color schemes.
library(scales)
library(gridExtra) # Arranging FOV QC plots into grids.
library(grid) # Displaying grids of plots.
library(patchwork) # `wrap_plots` as a replacement for `grid.arrange`
library(shinycssloaders) # Loading animations.
library(ggpointdensity)
library(fs) # File manipulation - needed to get file creation date when selecting most recent config YAML file.

## QC parameter data frame ----
config <- data.frame(
  name = character(),
  displayname = character(),
  min = numeric(),
  max = numeric(),
  default = numeric(),
  step = numeric()
)
config <- rbind(config,
                data.frame(
                  name = c("min_counts_per_cell", "min_features_per_cell", "max_proportion_neg_counts", "min_count_distribution", "max_area","min_signal_strength"),
                  displayname = c("Minimum number of counts per cell", "Minimum number of features per cell", "Maximum proportion of negative counts per cell", "Minimum complexity per cell", "Maximum area per cell", "Minimum signal strength per cell"),
                  metadata_name = c("nCount_RNA", "nFeature_RNA", "ProportionNegative", "Complexity", "Area", "SignalStrength"),
                  min = c(0, 0, 0, 0, 0, 0),
                  max = c(1000, 1000, 1, 200, 50000, 100),
                  default = c(200, 200, 0.1, 1, 30000, 4),
                  step = c(1, 1, 0.01, 1, 50, 1),
                  flag_as_minimum = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
                )
)

## Custom functions ----
#' Make a histogram.
#' 
#' @param plot ggplot2 plot
#' @export
DSPplotHistogram <- function(plot, n_bins = 100,
                             scale_x_log10 = TRUE,
                             theme_base_size = 12,
                             fill_color = "lightgray", outline_color = "black",
                             vline_xintercept = NULL, vline_color = "red", vline_type = "dashed",
                             facet_wrap_var = NULL,
                             x_lab = NULL, y_lab = NULL, title = NULL) {
  plot <- plot + 
    geom_histogram(bins = n_bins, fill = fill_color, color = outline_color) +
    labs(x = x_lab, y = y_lab) +
    theme_bw(base_size = theme_base_size) +
    theme(panel.grid = element_blank())
  
  # Add vline
  if(!is.null(vline_xintercept)) plot <- plot + geom_vline(xintercept = vline_xintercept, color = vline_color, linetype = vline_type)
  # Scale log10
  if(scale_x_log10) plot <- plot + scale_x_log10()
  # Facet wrap
  if(!is.null(facet_wrap_var)) plot <- plot + facet_wrap(~ get(facet_wrap_var))
  
  return(plot)
}

#' Map of which FOVs were flagged:
#' 
#' @param res Results object created by runFOVQC
#' @param shownames Logical for whether to display FOV names
#' @return For each bit, draws a plot of estimated FOV effects
#' @export
# mapFlaggedFOVs <- function(res, shownames = TRUE, outdir = NULL, plotwidth = NULL, plotheight = NULL) {
#   
#   if (!is.null(outdir)) {
#     png(paste0(outdir, "/flagged FOVs.png"), width = plotwidth, height = plotheight, units = "in", res = 300)
#   }
#   if (is.null(plotwidth)) {
#     plotwidth <- diff(range(res$xy[, 1])) * 1.5
#   }
#   if (is.null(plotheight)) {
#     plotheight <- diff(range(res$xy[, 2])) * 1.5
#   }
#   
#   plot(res$xy, cex = 0.2, asp = 1, pch = 16,
#        col = "black", 
#        main = "Flagged FOVs")
#   xy <- res$xy
#   for (f in unique(res$fov)) {
#     inds <- res$fov == f
#     rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), col = scales::alpha("dodgerblue2", 0.25))
#   }
#   for (f in res$flaggedfovs) {
#     inds <- res$fov == f
#     rect(min(xy[inds, 1]), min(xy[inds, 2]), max(xy[inds, 1]), max(xy[inds, 2]), col = scales::alpha("red", 0.5))
#     if (shownames) {
#       text(median(range(xy[inds, 1])), median(range(xy[inds, 2])), f, col = "green")
#     }
#   }
#   if (!is.null(outdir)) {
#     dev.off()
#   }
# }

# Rewrite for efficiency and speed:
mapFlaggedFOVs <- function(res, shownames = TRUE, slide_name = NULL) {
  xy <- res$xy
  fov_df <- res$fov %>%
    unique() %>%
    as.data.frame() %>%
    setNames("fov") %>%
    left_join(
      xy %>%
        as.data.frame() %>%
        setNames(c("x", "y")) %>%
        cbind(fov = res$fov) %>%
        group_by(fov) %>%
        summarise(
          xmin = min(x), xmax = max(x),
          ymin = min(y), ymax = max(y),
          .groups = "drop"
        ),
      by = "fov"
    )
  
  # Flagged FOVs
  flagged_df <- fov_df %>% filter(fov %in% res$flaggedfovs)
  
  # Downsample
  sampled_data <- as.data.table(xy)
  setnames(sampled_data, colnames(xy))  # Ensure column names are preserved
  sampled_data <- sampled_data[sample(.N, min(50000, .N))]  # Downsampling
  
  # Base plot
  p <- ggplot() +
    geom_point_rast(
      data = sampled_data, 
      aes(x = .data[[colnames(sampled_data)[1]]], y = .data[[colnames(sampled_data)[2]]]), 
      color = "black", size = 0.2
    ) +
    ggtitle(ifelse(is.null(slide_name), "Flagged FOVs", paste("Flagged FOVs |", slide_name))) +
    theme_bw() +
    coord_fixed()
  
  # Add all FOVs at once
  p <- p +
    geom_rect(
      data = fov_df, 
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = alpha("dodgerblue2", 0.25), color = NA
    ) +
    geom_rect(
      data = flagged_df, 
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = alpha("red", 0.5), color = NA
    )
  
  # Add names if required
  if (shownames) {
    p <- p +
      geom_text(
        data = flagged_df, 
        aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = fov),
        color = "green", size = 3
      )
  }
  
  return(p)
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
# Code rewritten for speed and efficiency:
FOVSignalLossSpatialPlot <- function(res, shownames = TRUE) {
  # Extract relevant data
  df <- as.data.frame(res$xy)
  colnames(df) <- c("x", "y")
  df$fov <- res$fov
  df$gridid <- res$gridinfo$gridid
  
  # Add residuals (log2 fold-change signal loss)
  df$resid <- res$totalcountsresids[match(df$gridid, names(res$totalcountsresids))]
  
  # Cap and rescale residuals to color scale (like original)
  df$resid_clipped <- pmax(pmin(51 + df$resid * 25, 101), 1)
  
  # Create a continuous scale manually matching the original palette
  col_scale <- colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)
  
  # Flagged FOVs
  flagged_fovs <- res$flaggedfovs_fortotalcounts
  
  # Downsample
  sampled_data <- as.data.table(df)
  setnames(sampled_data, colnames(df))  # Ensure column names are preserved
  sampled_data <- sampled_data[sample(.N, min(50000, .N))]  # Downsampling
  
  # Build plot
  p <- ggplot(sampled_data, aes(x = x, y = y)) +
    geom_point(aes(color = resid), size = 0.2) +
    scale_color_gradientn(
      colors = col_scale,
      limits = c(-2, 2),
      name = "Log2 FC"
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    ggtitle("Log2 fold-change in total counts\ncompared to similar regions")
  
  # Draw FOV boxes
  fov_boxes <- sampled_data %>%
    group_by(fov) %>%
    summarise(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y), .groups = "drop")
  
  p <- p + 
    geom_rect(data = fov_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", inherit.aes = FALSE)
  
  # Draw flagged FOVs in yellow
  if (length(flagged_fovs) > 0) {
    flagged_boxes <- fov_boxes %>% filter(fov %in% flagged_fovs)
    
    p <- p +
      geom_rect(data = flagged_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = NA, color = "yellow", size = 1.2, inherit.aes = FALSE)
    
    if (shownames) {
      flagged_boxes <- flagged_boxes %>%
        mutate(xcenter = (xmin + xmax)/2, ycenter = (ymin + ymax)/2)
      
      p <- p +
        geom_text(data = flagged_boxes, aes(x = xcenter, y = ycenter, label = fov),
                  color = "green", size = 3, inherit.aes = FALSE)
    }
  }
  
  return(p)
}
# Original code:
# FOVSignalLossSpatialPlot.orig <- function(res, shownames = TRUE, outdir = NULL, plotwidth = NULL, plotheight = NULL) {
# 
#   if (!is.null(outdir)) {
#     png(paste0(outdir, "/signal loss.png"), width = plotwidth, height = plotheight, units = "in", res = 300)
#   }
#   if (is.null(plotwidth)) {
#     plotwidth <- diff(range(res$xy[, 1])) * 1.5
#   }
#   if (is.null(plotheight)) {
#     plotheight <- diff(range(res$xy[, 2])) * 1.5
#   }
# 
#   plot(res$xy, cex = 0.2, asp = 1, pch = 16,
#        col = colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)[
#          pmax(pmin(51 + res$totalcountsresids[match(res$gridinfo$gridid, names(res$totalcountsresids))] * 25, 101), 1)],
#        main = "Log2 fold-change in total counts compared to similar regions")
#   for (f in unique(res$fov)) {
#     inds <- res$fov == f
#     rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "black")
#   }
#   for (f in res$flaggedfovs_fortotalcounts) {
#     inds <- res$fov == f
#     rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "yellow", lwd = 2)
#     if (shownames) {
#       text(median(range(res$xy[inds, 1])), median(range(res$xy[inds, 2])), f, col = "green")
#     }
#   }
#   legend("right", pch = 16,
#          col = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
#          legend = rev(c("< -2", -1, 0, 1, "> 2")))
#   if (!is.null(outdir)) {
#     dev.off()
#   }
# }

#' Spatial plots of FOV effects:
#' 
#' @param res Results object created by runFOVQC
#' @param outdir Directory where png plots are printed
#' @param bits Which bits to plot. Defaults to "flagged_reportercycles" bits, but can also plot "all" or "flagged_bits".
#' @param plotwidth Width in inches of png plots
#' @param plotheight Height in inches of png plots
#' @return For each bit, draws a plot of estimated FOV effects
#' @export
FOVEffectsSpatialPlots <- function(res, bits = "flagged_reportercycles") {
  # Determine which bits to plot
  bits_to_plot <- NULL
  bit_names <- NULL
  
  if (bits == "flagged_reportercycles") {
    flaggedreportercycles <- colnames(res$flags_per_fov_x_reportercycle)[colSums(res$flags_per_fov_x_reportercycle >= 0.5) > 0]
    bit_names <- paste0(rep(flaggedreportercycles, each = 4), rep(c("B", "G", "R", "Y"), length(flaggedreportercycles)))
    bits_to_plot <- match(bit_names, colnames(res$resid))
  } else if (bits == "flagged_bits") {
    bits_to_plot <- which(colSums(res$fovstats$flag) > 0)
    bit_names <- colnames(res$resid)[bits_to_plot]
  } else if (bits == "all") {
    bits_to_plot <- seq_len(ncol(res$resid))
    bit_names <- colnames(res$resid)
  } else {
    bits_to_plot <- match(bits, colnames(res$resid))
    bit_names <- bits
  }
  
  plots <- list()
  col_scale <- colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)
  
  for (i in seq_along(bits_to_plot)) {
    bit_index <- bits_to_plot[i]
    bit_name <- bit_names[i]
    
    df <- as.data.frame(res$xy)
    colnames(df) <- c("x", "y")
    df$fov <- res$fov
    df$gridid <- res$gridinfo$gridid
    df$resid <- res$resid[match(df$gridid, rownames(res$resid)), bit_index]
    df$resid_clipped <- pmax(pmin(51 + df$resid * 50, 101), 1)
    
    # Compute FOV boxes
    fov_boxes <- df %>%
      group_by(fov) %>%
      summarise(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y), .groups = "drop")
    
    # Flagged FOVs for this bit
    flagged_fovs <- rownames(res$fovstats$flag)[res$fovstats$flag[, bit_index] > 0]
    
    # Downsample
    sampled_data <- as.data.table(df)
    setnames(sampled_data, colnames(df))  # Ensure column names are preserved
    sampled_data <- sampled_data[sample(.N, min(50000, .N))]  # Downsampling
    
    p <- ggplot(sampled_data, aes(x = x, y = y)) +
      geom_point(aes(color = resid), size = 0.2) +
      scale_color_gradientn(
        colors = col_scale,
        limits = c(-1, 1),
        name = "Log2 FC"
      ) +
      coord_fixed() +
      ggtitle(paste0(bit_name, ": log2(fold-change)\nfrom comparable regions elsewhere")) +
      theme_minimal() +
      theme(
        legend.position = "right",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      geom_rect(data = fov_boxes,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = NA, color = "black", inherit.aes = FALSE)
    
    # Highlight flagged FOVs
    if (length(flagged_fovs) > 0) {
      flagged_boxes <- fov_boxes %>% filter(fov %in% flagged_fovs)
      p <- p +
        geom_rect(data = flagged_boxes,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = NA, color = "yellow", size = 1.2, inherit.aes = FALSE)
    }
    
    plots[[bit_name]] <- p
  }
  
  return(plots)
}
# Original code:
# FOVEffectsSpatialPlots.orig <- function(res, outdir = NULL, bits = "flagged_reportercycles", plotwidth = NULL, plotheight = NULL) {
#   
#   if (is.null(plotwidth)) {
#     plotwidth <- diff(range(res$xy[, 1])) * 1.5
#   }
#   if (is.null(plotheight)) {
#     plotheight <- diff(range(res$xy[, 2])) * 1.5
#   }
#   bits_to_plot <- match(bits, colnames(res$fovstats$flag))
#   if (bits == "flagged_reportercycles") {
#     flaggedreportercycles <- colnames(res$flags_per_fov_x_reportercycle)[colSums(res$flags_per_fov_x_reportercycle >= 0.5) > 0]
#     names_of_bits_to_plot  <- paste0(rep(flaggedreportercycles, each = 4), rep(c("B", "G", "R","Y"), length(flaggedreportercycles)))
#     bits_to_plot <- match(names_of_bits_to_plot, colnames(res$resid))
#   }
#   if (bits == "flagged_bits") {
#     bits_to_plot  <- which(colSums(res$fovstats$flag) > 0)
#   }
#   if (bits == "all") {
#     bits_to_plot <- 1:ncol(res$resid)
#   }
#   temp <- sapply(bits_to_plot, function(i) {
#     if (!is.null(outdir)) {
#       png(paste0(outdir, "/", make.names(colnames(res$resid)[i]), ".png"), width = plotwidth, height = plotheight, units = "in", res = 300)
#     }
#     par(mar = c(0,0,2,0))
#     plot(res$xy, cex = 0.2, asp = 1, pch = 16,
#          col = colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)[
#            pmax(pmin(51 + res$resid[match(res$gridinfo$gridid, rownames(res$resid)), i] * 50, 101), 1)], 
#          main = paste0(colnames(res$resid)[i], ": log2(fold-change)\nfrom comparable regions elsewhere"))
#     for (f in unique(res$fov)) {
#       inds <- res$fov == f
#       rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), border = "black")
#     }
#     for (f in rownames(res$fovstats$flag)[res$fovstats$flag[, i] > 0]) {
#       inds <- res$fov == f
#       rect(min(res$xy[inds, 1]), min(res$xy[inds, 2]), max(res$xy[inds, 1]), max(res$xy[inds, 2]), lwd = 2, border = "yellow")
#     }
#     legend("right", pch = 16,
#            col = rev(c("darkblue", "blue", "grey80", "red", "darkred")),
#            legend = rev(c("< -1", -0.5, 0, 0.5, "> 1")))
#     if (!is.null(outdir)) {
#       dev.off()
#     }
#   })
# }

## Function to recursively assign variables.
assignVarsEnv <- function(yaml_list, env = .GlobalEnv) { # env = new.env()
  assignRecursive <- function(lst, parent_keys = NULL) {
    for (key in names(lst)) {
      full_key <- paste(c(parent_keys, key), collapse = "_")
      value <- lst[[key]]
      
      if (is.list(value)) {
        assignRecursive(value, c(parent_keys, key))
      } else {
        assign(key, value, envir = env)
      }
    }
  }
  
  assignRecursive(yaml_list)
  # return(env)  # Return the environment so you can use it
}

# Function to validate and process config variables.
validateProcessConfig <- function(config_var_metadata) { 
  config_vars <- read.csv(config_var_metadata, stringsAsFactors = FALSE)
  error_msg_list <- list()
  
  for (i in seq_len(nrow(config_vars))) {
    var_name <- config_vars$variable[i]
    required <- as.logical(config_vars$required[i])
    default_value <- config_vars$default_value[i]
    split_string <- as.logical(config_vars$split_string[i])
    
    # Debugging prints
    message(paste("Processing:", var_name))
    message(paste(" - Required:", required))
    message(paste(" - Default Value (raw):", default_value))
    message(paste(" - Exists in .GlobalEnv?", exists(var_name, envir = .GlobalEnv, inherits = FALSE)))
    
    # Convert default_value to the appropriate type
    if (!is.na(default_value) && default_value != "") {
      if (grepl("^\\d+$", default_value)) {  
        # Integer check (only digits, no decimal)
        default_value <- as.integer(default_value)
      } else if (grepl("^\\d+\\.\\d+$", default_value)) {  
        # Decimal check (digits + decimal point)
        default_value <- as.numeric(default_value)
      }
    }
    
    message(paste(" - Default Value (converted):", default_value, "Type:", typeof(default_value)))
    
    if (exists(var_name, envir = .GlobalEnv, inherits = FALSE)) { 
      value <- get(var_name, envir = .GlobalEnv)
      message(paste(" - Current Value:", value))
      
      if (required && flagVariable(value)) {
        error_msg_list[[var_name]] <- paste("Input missing for", var_name)
      }
      
      if (!required && flagVariable(value) && default_value != "") {
        assign(var_name, default_value, envir = .GlobalEnv)
        message(paste(" - Assigned default:", var_name, "=", get(var_name, envir = .GlobalEnv)))
      }
      
      if (!flagVariable(value) && split_string) {
        assign(var_name, strsplit(value, ",")[[1]] %>% trimws(), envir = .GlobalEnv)
        message(paste(" - Split and assigned:", get(var_name, envir = .GlobalEnv)))
      }
      
    } else {
      if (!required && default_value != "") {
        assign(var_name, default_value, envir = .GlobalEnv)
        message(paste(" - Assigned default:", var_name, "=", get(var_name, envir = .GlobalEnv)))
      }
    }
  }
  
  return(error_msg_list)
}

## Load configuration variables ----
# # Get list of files in `config` directory.
message("Importing configuration variables.")
file_info <- dir_info("config")
# Sort by creation time and filter to include only config files
config_expt <- file_info %>%
  dplyr::filter(type=="file" & base::grepl("config_", path)) %>%
  dplyr::arrange((desc(change_time))) %>% # change_time instead of birth_time because ... why is the birth_time off??
  dplyr::select(path) %>% unlist
config_expt <- yaml::read_yaml(config_expt)
# config_env <- assignVarsEnv(config_expt)
slide_name_var <- config_expt$experiment$annotation$slide_name_var

## Load data ----
# Load FOV QC results
message("Loading FOV QC results.")
fov_qc_res_list <- readRDS("Rdata/fov_qc_res_list.rds")
# Invert the list to make it easier to make grids. 
new_top_levels <- fov_qc_res_list[[1]] %>% names
new_lower_levels <- fov_qc_res_list %>% names
fov_qc_res_list_inv <- list()
for(top_level in new_top_levels) {
  fov_qc_res_list_inv[[top_level]] <- list()
  for(lower_level in new_lower_levels) {
    fov_qc_res_list_inv[[top_level]][[lower_level]] <- fov_qc_res_list[[lower_level]][[top_level]]
  }
}
rm(fov_qc_res_list)
fov_qc_res_list <- fov_qc_res_list_inv
# Extract unique max_prop_loss values (assuming all slides have the same structure)
max_prop_loss_values <- names(fov_qc_res_list)
message("FOV QC results loaded.")

## Pre-render FOV QC graphs ----
# The following commented-out block of code has been removed as we're switching from pre-rendering to on-demand rendering.
# message("Pre-rendering FOV QC graphs. This may take a while depending on how many slides you have.")
# message("Sit back and hold tight, and you'll be alerted when pre-rendering has finished.")
# fov_qc_grids <- list()
# fov_qc_grids[["flagged_fovs"]] <- list()
# fov_qc_grids[["signal_loss"]] <- list()
# fov_qc_grids[["fov_effects"]] <- list()
# for(max_prop_loss_val in names(fov_qc_res_list)) {
#   # Generate FOV plots for all slides
#   # Flagged FOV plots
#   withProgress(message = "Rendering flagged-FOV plots", value = 0, {
#     n_slides <- length(names(fov_qc_res_list[[1]]))
#     
#     flagged_fov_plots <- lapply(seq_along(names(fov_qc_res_list[[1]])), function(i) {
#       slide <- names(fov_qc_res_list[[1]])[i]
#       
#       incProgress(1 / n_slides, detail = paste("Slide", slide))
#       
#       mapFlaggedFOVs(fov_qc_res_list[[max_prop_loss_val]][[slide]], shownames = FALSE, slide_name = slide)
#     })
#   })
#   # flagged_fov_plots <- lapply(names(fov_qc_res_list[[1]]), function(slide) {
#   #   mapFlaggedFOVs(fov_qc_res_list[[max_prop_loss_val]][[slide]], shownames = F, slide_name = slide)
#   # })
#   # Signal loss spatial plots
#   withProgress(message = "Rendering signal-loss spatial plots", value = 0, {
#     n_slides <- length(names(fov_qc_res_list[[1]]))
#     
#     signal_loss_plots <- lapply(seq_along(names(fov_qc_res_list[[1]])), function(i) {
#       slide <- names(fov_qc_res_list[[1]])[i]
#       
#       incProgress(1 / n_slides, detail = paste("Slide", slide))
#       
#       FOVSignalLossSpatialPlot(fov_qc_res_list[[max_prop_loss_val]][[slide]], shownames = F)
#     })
#   })
#   # signal_loss_plots <- lapply(names(fov_qc_res_list[[1]]), function(slide) {
#   #   FOVSignalLossSpatialPlot(fov_qc_res_list[[max_prop_loss_val]][[slide]], shownames = F)
#   # })
#   # Per-bit FOV effects spatial plots
#   withProgress(message = "Rendering FOV-effects spatial plots", value = 0, {
#     n_slides <- length(names(fov_qc_res_list[[1]]))
#     
#     fov_effect_plots <- lapply(seq_along(names(fov_qc_res_list[[1]])), function(i) {
#       slide <- names(fov_qc_res_list[[1]])[i]
#       
#       incProgress(1 / n_slides, detail = paste("Slide", slide))
#       
#       FOVEffectsSpatialPlots(fov_qc_res_list[[max_prop_loss_val]][[slide]], bits = "flagged_reportercycles")
#     })
#   })
#   # fov_effect_plots <- lapply(names(fov_qc_res_list[[1]]), function(slide) {
#   #   FOVEffectsSpatialPlots(fov_qc_res_list[[max_prop_loss_val]][[slide]], bits = "flagged_reportercycles")
#   # })
#   
#   # Arrange plots in two columns.
#   # Flagged FOV plots
#   plot_grid <- grid.arrange(grobs = flagged_fov_plots, ncol = 2)
#   fov_qc_grids[["flagged_fovs"]][[max_prop_loss_val]] <- plot_grid
#   rm(plot_grid)
#   # Signal loss spatial plots
#   plot_grid <- grid.arrange(grobs = signal_loss_plots, ncol = 2)
#   fov_qc_grids[["signal_loss"]][[max_prop_loss_val]] <- plot_grid
#   rm(plot_grid)
#   # Per-bit FOV effects spatial plots
#   plot_grid <- grid.arrange(grobs = unlist(fov_effect_plots, recursive = FALSE), ncol = 2)
#   fov_qc_grids[["fov_effects"]][[max_prop_loss_val]] <- plot_grid
#   rm(plot_grid)
# }
# message("Pre-rendering complete!")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# UI layout --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Layout
ui <- fluidPage(
  # theme = bs_theme(version = 4, bootswatch = "cerulean"), # Use a theme
  titlePanel("CosMx QC dashboard"),
  
  ## Cell QC ----
  h1("Cell QC"),
  sidebarLayout(
    sidebarPanel(
      # Dynamically generated inputs
      uiOutput("dynamicInputs"),
      
      # Input value display
      verbatimTextOutput("inputValues"),
      
      # Download button for cutoff values
      downloadButton("download_yaml", "Download QC config (YAML)")
    ),
    mainPanel(
      # Dynamically generated histogram tabs
      withSpinner(uiOutput("histogramTabs")),
      
      # QC summary table
      withSpinner(DTOutput("flaggedTable")),
      
      # Download button for the QC summary table
      downloadButton("download_flagged_table", "Download flagged cell counts (CSV)")
    )
  )
  ,
  
  ## FOV QC ----
  ### Counts in space ----
  h1("FOV QC"),
  h2("Counts in space"),
  numericInput(
    inputId = "min_counts_per_cell",
    label = "Minimum counts-per-cell cutoff",
    value = 50,
    min = 0,
    max = 6000,
    step = 50
  ),
  "Cells that fall below the cutoff are highlighted in green.",
  br(),
  withSpinner(plotOutput(outputId = "counts_in_space", height = "1000px", width = "1400px"))
  ,
  
  ### Flagged FOVs ----
  fluidPage(
    titlePanel("Flagged FOVs"),
    
    sidebarLayout(
      sidebarPanel(
        selectInput("max_prop_loss", "Max Prop Loss", choices = names(fov_qc_res_list)),
        selectInput("selected_slide", "Slide", choices = names(fov_qc_res_list[[1]]))
      ),
      
      mainPanel(
        fluidRow(
          column(width = 12,
                 shinycssloaders::withSpinner(
                   plotOutput("flaggedFOVPlot", height = "1000px")
                 )
          )
        ),
        fluidRow(
          column(width = 12,
                 shinycssloaders::withSpinner(
                   plotOutput("signalLossPlot", height = "1000px")
                 )
          )
        ),
        fluidRow(
          column(width = 12,
                 shinycssloaders::withSpinner(
                   plotOutput("fovEffectPlot", height = "1800px", width = "100%")  # bigger for multi-plot grid
                 )
          )
        )
      )
    )
  )
  # h2("Flagged FOVs"),
  # selectInput("max_prop_loss", "max_prop_loss", choices = max_prop_loss_values),
  # selectInput("selected_slide", "Slide", choices = names(fov_qc_res_list[[1]])),
  # "Flagged FOVs are highlighted in red.",
  # withSpinner(plotOutput("flaggedFOVPlot", height = "1000px")),
  # withSpinner(plotOutput("signalLossPlot", height = "1000px")),
  # withSpinner(plotOutput("fovEffectPlot", height = "1000px"))
  # # The following code block is commented out as we're switching from pre-rendering to on-demand rendering
  # # withSpinner(uiOutput("flagged_fov_plots", width = "1400px", height = "4000px")) # , height = "3000px", width = "1400px"
  
)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Server-side functions --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
server <- function(input, output, session) {
  ## Load required data ----
  # Do it inside server function to prevent them from being loaded multiple times, I guess.
  # Load cell-QC stats.
  cell_stats <- readRDS("Rdata/flagged_metadata_raw.rds") %>% as.data.table
  # Change NaNs (which result from 0/0 division) in `ProportionNegative` to 0
  cell_stats <- cell_stats %>%
    mutate(ProportionNegative = ifelse(is.nan(ProportionNegative), 0, ProportionNegative))
  
  ## Input generation ----
  output$dynamicInputs <- renderUI({
    lapply(1:nrow(config), function(i) {
      # fluidRow(
      # Slider inputs
      # column(
      #   width = 8,
      #   sliderInput(
      #     inputId = paste0(config$name[i], "_slider"),
      #     label = config$displayname[i],
      #     min = config$min[i],
      #     max = config$max[i],
      #     value = config$default[i],
      #     step = config$step[i]
      #   )
      # )
      # ,
      # Numeric inputs
      # column(
      # width = 4,
      numericInput(
        inputId = paste0(config$name[i], "_numeric"),
        label = config$displayname[i], # NULL
        value = config$default[i],
        min = config$min[i],
        max = config$max[i],
        step = config$step[i]
      )
      # )
      # )
    })
  })
  
  ## Debouncing ----
  # # Store latest value.
  # latest_values_sliders <- lapply(1:nrow(config), function(i) {
  #   reactiveVal(reactive(input[[paste0(config$name[i], "_slider")]]))
  # })
  # # latest_values_numeric <- lapply(1:nrow(config), function(i) {
  # #   reactiveVal(reactive(input[[paste0(config$name[i], "_numeric")]]))
  # # })
  # 
  # # Debounce.
  # debounced_values <- lapply(1:length(latest_values_sliders), function(i) {
  #   debounce(reactive(latest_values_sliders[[i]]), 500) # 500ms debounce
  # })
  # # debounced_values_sliders <- lapply(1:length(latest_values_sliders), function(i) {
  # #   debounce(reactive(latest_values_sliders[[i]]), 500) # 500ms debounce
  # # })
  # # debounced_values_numeric <- lapply(1:length(latest_values_numeric), function(i) {
  # #   debounce(reactive(latest_values_numeric[[i]]), 500) # 500ms debounce
  # # })
  
  ## Slider and numeric input synchronization ----
  # # Update inputs only after debounce.
  # observe({
  #   for (i in 1:nrow(config)) {
  #     local({
  #       name <- config$name[i]
  # 
  #       # Update numeric inputs.
  #       observeEvent(debounced_values[[i]](), {
  #         updateNumericInput(session, paste0(name, "_numeric"), value = debounced_values[[i]]())
  #       }) # End `observeEvent()`.
  # 
  #       # Update slider inputs.
  #       observeEvent(input[[paste0(name, "_numeric")]], {
  #         updateSliderInput(session, paste0(name, "_slider"), value = input[[paste0(name, "_numeric")]])
  #       })
  # 
  #     }) # End `local()`.
  # 
  #   } # End `for()` loop.
  # 
  # }) # End `observe()`.
  
  # # Original code (no debounce):
  # observe({
  #   for (i in 1:nrow(config)) {
  #     local({
  #       name <- config$name[i]
  # 
  #       # Update numeric input.
  #       observeEvent(input[[paste0(name, "_slider")]], {
  #         updateNumericInput(session, paste0(name, "_numeric"), value = input[[paste0(name, "_slider")]])
  #       }) # End `observeEvent()`.
  # 
  #       # Update slider input.
  #       observeEvent(input[[paste0(name, "_numeric")]], {
  #         updateSliderInput(session, paste0(name, "_slider"), value = input[[paste0(name, "_numeric")]])
  #       }) # End `observeEvent()`.
  # 
  #     }) # End `local()`.
  #   } # End `for()` loop.
  # }) # End `observe()`.
  
  ## Flagged cell calculation ----
  flagged_table_data <- reactive({
    flagged_matrix <- sapply(1:nrow(config), function(i) {
      stat <- config$name[i]
      metadata_name <- config$metadata_name[i]
      cutoff <- input[[paste0(stat, "_numeric")]]
      is_minimum <- config$flag_as_minimum[i]
      if (is_minimum) {
        cell_stats[[metadata_name]] < cutoff
      } else {
        cell_stats[[metadata_name]] > cutoff
      }
    })
    flagged_counts <- colSums(flagged_matrix)
    total_flagged <- sum(rowSums(flagged_matrix) > 0)
    total_passed <- nrow(cell_stats) - total_flagged
    results <- data.frame(
      Stat = c(config$displayname, "Total Flagged", "Total Passed", "All cells"),
      FlaggedCount = c(flagged_counts, total_flagged, total_passed, nrow(cell_stats))
    )
    return(results)
  })
  
  ## Render flagged cell data table ----
  # (Previously Flagged cell calculation)
  output$flaggedTable <- renderDT({
    datatable(flagged_table_data(), options = list(pageLength = 10))
  })
  # output$flaggedTable <- renderDT({
  #   flagged_matrix <- sapply(1:nrow(config), function(i) {
  #     stat <- config$name[i]
  #     metadata_name <- config$metadata_name[i]
  #     cutoff <- input[[paste0(stat, "_numeric")]]
  #     is_minimum <- config$flag_as_minimum[i]
  #     if (is_minimum) {
  #       cell_stats[[metadata_name]] < cutoff
  #     } else {
  #       cell_stats[[metadata_name]] > cutoff
  #     }
  #   })
  #   flagged_counts <- colSums(flagged_matrix)
  #   total_flagged <- sum(rowSums(flagged_matrix) > 0)
  #   total_passed <- nrow(cell_stats) - total_flagged
  #   results <- data.frame(
  #     Stat = c(config$displayname, "Total Flagged", "Total Passed", "All cells"),
  #     `Number of flagged cells` = c(flagged_counts, total_flagged, total_passed, nrow(cell_stats))
  #   )
  #   datatable(results, options = list(pageLength = 10))
  # })
  
  ## Download handler for flagged cell data table ----
  output$download_flagged_table <- downloadHandler(
    filename = function() {
      paste0("flagged_cell_counts_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(flagged_table_data(), file, row.names = FALSE)
    }
  )
  
  ## Histogram tab generation ----
  output$histogramTabs <- renderUI({
    # Create a tab for each QC metric.
    nav_items <- lapply(1:nrow(config), function(i) {
      nav_panel(
        title = config$displayname[i],
        plotOutput(outputId = paste0("histogram_", config$name[i]))
      )
    })
    
    # Return the tabs.
    do.call(navset_pill, nav_items)
  })
  
  ## Histogram rendering ----
  observe({
    for (i in 1:nrow(config)) {
      local({
        # Capture the current iteration's variables.
        stat <- config$name[i]
        metadata_name <- config$metadata_name[i]
        display_name <- config$displayname[i]
        cutoff_input <- paste0(stat, "_numeric")
        step <- config$step[i]
        
        # Render histogram dynamically.
        output[[paste0("histogram_", stat)]] <- renderPlot({
          cutoff <- input[[cutoff_input]] # Get the current cutoff value.
          
          plot <- ggplot(cell_stats, aes_string(x = metadata_name))
          DSPplotHistogram(plot, n_bins = 100,
                           scale_x_log10 = FALSE,
                           theme_base_size = 12,
                           fill_color = "lightgray", outline_color = "black",
                           vline_xintercept = cutoff, vline_color = "red", vline_type = "dashed",
                           facet_wrap_var = NULL,
                           x_lab = display_name, y_lab = "Number of cells", title = NULL)
        })
      }) # End of local()
    }
  })
  
  ## QC configuration download ----
  output$download_yaml <- downloadHandler(
    filename = function() {
      paste0("cosmx_qc_inputs_", Sys.Date(), ".yaml")
    },
    content = function(file) {
      # Create a named list of all the current input values for Cell QC
      qc_inputs <- setNames(
        lapply(config$name, function(stat) {
          input[[paste0(stat, "_numeric")]]
        }),
        config$name
      )
      
      # Write to YAML
      yaml::write_yaml(qc_inputs, file)
    }
  )
  
  ## Counts-in-space scatter plot rendering ----
  # Use reactive expressions to store filtered data instead of recalculating everything inside renderPlot():
  # Precompute filtered data. (The filtering now happens only when inputs change, rather than in every plot re-render.)
  # filtered_data <- reactive({
  #   cell_stats[cell_stats$nCount_RNA > input$min_counts_per_cell,] %>% as.data.table
  # })
  observe({
    counts_in_space_cutoff <- log10(input$min_counts_per_cell + 1)  # Define your cutoff value (e.g., 100 raw counts)
    output$counts_in_space <- renderPlot({
      # sampled_data <- filtered_data()[sample(.N, min(100000, .N))]  # Downsample number of points to 100k to speed up rendering.
      sampled_data <- cell_stats[sample(.N, min(100000, .N))]
      
      ggplot(data = sampled_data, aes(x = CenterX_global_px, y = CenterY_global_px)) +
        # First layer: Full gradient for all points
        geom_point_rast(aes(color = log10(nCount_RNA + 1)), size = 0.1, alpha = 0.7) +
        scale_color_gradientn(
          colors = viridis(n = 100, option = "magma"),  # Full gradient applied here
          name = "Log10(Total Counts + 1)",
          breaks = log10(c(200, 1000, 2000) + 1),
          labels = c("200", "1000", "2000")
        ) +
        # Second layer: Green points for values below cutoff
        geom_point_rast(
          data = sampled_data[log10(nCount_RNA + 1) <= counts_in_space_cutoff],  # Subset for cutoff
          color = "green",  # Set color outside aes() to avoid scale conflicts
          size = 0.1, alpha = 0.7
        ) +
        labs(x = "x", y = "y") +
        theme_bw(base_size = 22) +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          legend.key.height = unit(2, "lines")
        ) +
        coord_fixed() +  # Keep aspect ratio 1:1
        facet_wrap(~get(slide_name_var)) # slide_ID
      #  scales = 'free' prevents unnecessary plotting, according to ChatGPT. But it can't be used with `coord_fixed()`. 
      
    })
  })
  
  
  # Flagged FOV plot rendering ----
  # Dynamically update slide choices when max_prop_loss changes
  observeEvent(input$max_prop_loss, {
    updateSelectInput(session, "selected_slide",
                      choices = names(fov_qc_res_list[[input$max_prop_loss]]),
                      selected = names(fov_qc_res_list[[input$max_prop_loss]])[1])
  })
  
  # --- Flagged FOVs ---
  output$flaggedFOVPlot <- renderPlot({
    req(input$max_prop_loss, input$selected_slide)
    
    withProgress(message = "Rendering flagged-FOV plots", value = 0, {
      res <- fov_qc_res_list[[input$max_prop_loss]][[input$selected_slide]]
      Sys.sleep(0.1)  # optional to simulate progress
      mapFlaggedFOVs(res, shownames = FALSE, slide_name = input$selected_slide)
    })
  })
  
  # --- Signal loss ---
  output$signalLossPlot <- renderPlot({
    req(input$max_prop_loss, input$selected_slide)
    
    withProgress(message = "Rendering signal-loss plots", value = 0, {
      res <- fov_qc_res_list[[input$max_prop_loss]][[input$selected_slide]]
      Sys.sleep(0.1)  # optional to simulate progress
      FOVSignalLossSpatialPlot(res, shownames = FALSE)
    })
  })
  
  # --- Per-bit FOV effects ---
  output$fovEffectPlot <- renderPlot({
    req(input$max_prop_loss, input$selected_slide)
    
    res <- fov_qc_res_list[[input$max_prop_loss]][[input$selected_slide]]
    plots <- FOVEffectsSpatialPlots(res, bits = "flagged_reportercycles")
    
    # Keep only valid ggplots
    plots <- Filter(function(p) inherits(p, "gg"), plots)
    if (length(plots) == 0) {
      grid::grid.newpage()
      grid::grid.text("No reporter-cycle QC plots available for this slide.", x = 0.05, just = "left")
      return()
    }
    
    # Tidy up margins
    plots <- lapply(plots, function(p) p + theme(plot.margin = margin(5, 5, 5, 5)))
    
    # Combine using patchwork
    p_combined <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect") +
      patchwork::plot_annotation(title = "FOV Reporter Cycle Effects")
    
    # ✅ Key: use print() to render patchwork in Shiny
    print(p_combined)
  }, width = 4800, height = 10000)
  
  # The following commented-out code block has been superceded since we switched from pre-rendering to on-demand rendering.
  # observe({
  #   # output$flagged_fov_plots <- renderPlot({
  #   #   req(input$max_prop_loss)
  #   # 
  #   #   grid.draw(fov_qc_grids[[input$max_prop_loss]])
  #   # })
  #   output$flagged_fov_plots <- renderUI({
  #     req(input$max_prop_loss)  # Ensure a valid selection
  #     plotOutput(outputId = "flaggedFOVPlot")
  #     plotOutput(outputId = "signalLossPlot")
  #     plotOutput(outputId = "fovEffectPlot")
  #   })
  #   
  #   # Flagged FOV plots
  #   output$flaggedFOVPlot <- renderPlot({
  #     req(fov_qc_grids[["flagged_fovs"]][[input$max_prop_loss]])  # Ensure a valid plot exists
  #     grid.newpage()
  #     grid.draw(fov_qc_grids[["flagged_fovs"]][[input$max_prop_loss]])
  #   }, height = 2000, width = 1400)
  #   # Signal loss plots
  #   output$signalLossPlot <- renderPlot({
  #     req(fov_qc_grids[["signal_loss"]][[input$max_prop_loss]])  # Ensure a valid plot exists
  #     grid.newpage()
  #     grid.draw(fov_qc_grids[["signal_loss"]][[input$max_prop_loss]])
  #   }, height = 2000, width = 1400)
  #   # FOV effect plots
  #   output$fovEffectPlot <- renderPlot({
  #     req(fov_qc_grids[["fov_effects"]][[input$max_prop_loss]])  # Ensure a valid plot exists
  #     grid.newpage()
  #     grid.draw(fov_qc_grids[["fov_effects"]][[input$max_prop_loss]])
  #   }, height = 2000, width = 1400)
  # })
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Application --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
shinyApp(ui, server)