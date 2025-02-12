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
library(shinycssloaders) # Loading animations.
library(ggpointdensity)

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
message("Pre-rendering FOV QC graphs. This may take a while depending on how many slides you have.")
message("Sit back and hold tight, and you'll be alerted when pre-rendering has finished.")
fov_qc_grids <- list()
for(max_prop_loss_val in names(fov_qc_res_list)) {
  # Generate FOV plots for all slides
  flagged_fov_plots <- lapply(names(fov_qc_res_list[[1]]), function(slide) {
    mapFlaggedFOVs(fov_qc_res_list[[max_prop_loss_val]][[slide]], shownames = F, slide_name = slide)
  })
  
  # Arrange plots in two columns.
  plot_grid <- grid.arrange(grobs = flagged_fov_plots, ncol = 2)
  fov_qc_grids[[max_prop_loss_val]] <- plot_grid
  
  rm(plot_grid)
}
message("Pre-rendering complete!")

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
      verbatimTextOutput("inputValues")
    ),
    mainPanel(
      # Dynamically generated histogram tabs
      withSpinner(uiOutput("histogramTabs")),
      
      # QC summary table
      withSpinner(DTOutput("flaggedTable"))
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
  "Cells that fall below the cutoff will be highlighted in green.",
  br(),
  withSpinner(plotOutput(outputId = "counts_in_space", height = "1000px", width = "1400px"))
  ,
  
  ### Flagged FOVs ----
  h2("Flagged FOVs"),
  selectInput("max_prop_loss", "max_prop_loss", choices = max_prop_loss_values),
  "Flagged FOVs will be highlighted in red.",
  withSpinner(uiOutput("flagged_fov_plots", width = "1400px", height = "4000px")) # , height = "3000px", width = "1400px"
  
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
  output$flaggedTable <- renderDT({
    # cell_stats <- readRDS("Rdata/flagged_metadata_raw.rds")
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
      `Flagged count` = c(flagged_counts, total_flagged, total_passed, nrow(cell_stats))
    )
    datatable(results, options = list(pageLength = 10))
  })
  
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
        facet_wrap(~slide_ID)
      #  scales = 'free' prevents unnecessary plotting, according to ChatGPT. But it can't be used with `coord_fixed()`. 
      
    })
  })
  
  
  # Flagged FOV plot rendering ----
  observe({
    # output$flagged_fov_plots <- renderPlot({
    #   req(input$max_prop_loss)
    # 
    #   grid.draw(fov_qc_grids[[input$max_prop_loss]])
    # })
    output$flagged_fov_plots <- renderUI({
      req(input$max_prop_loss)  # Ensure a valid selection
      plotOutput(outputId = "flaggedFOVPlot")
    })
    
    output$flaggedFOVPlot <- renderPlot({
      req(fov_qc_grids[[input$max_prop_loss]])  # Ensure a valid plot exists
      grid.newpage()
      grid.draw(fov_qc_grids[[input$max_prop_loss]])
    }, height = 2000, width = 1400)
  })
}


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Application --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
shinyApp(ui, server)