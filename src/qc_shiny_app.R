# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Setup --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Required packages ----
library(shiny)
library(tidyverse)
library(DT) # For the display tables.
library(bslib) # Styling of Shiny apps.

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
  sidebarLayout(
    sidebarPanel(
      # Dynamically generated inputs
      uiOutput("dynamicInputs"),
      
      # Input value display
      verbatimTextOutput("inputValues")
    ),
    mainPanel(
      # Dynamically generated histogram tabs
      uiOutput("histogramTabs"),
      
      # QC summary table
      DTOutput("flaggedTable")
    )
  )
  # ,
  # 
  # ## FOV QC ----
  # card(
  #   card_header("FOV QC"),
  #   layout_sidebar()
  # )
)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Server-side functions --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
server <- function(input, output, session) {
  ## Input generation ----
  output$dynamicInputs <- renderUI({
    lapply(1:nrow(config), function(i) {
      fluidRow(
        column(
          width = 8,
          sliderInput(
            inputId = paste0(config$name[i], "_slider"),
            label = config$displayname[i],
            min = config$min[i],
            max = config$max[i],
            value = config$default[i],
            step = config$step[i]
          )
        ),
        column(
          width = 4,
          numericInput(
            inputId = paste0(config$name[i], "_numeric"),
            label = NULL,
            value = config$default[i],
            min = config$min[i],
            max = config$max[i],
            step = config$step[i]
          )
        )
      )
    })
  })
  
  ## Slider and numeric input synchronization ----
  observe({
    for (i in 1:nrow(config)) {
      local({
        name <- config$name[i]
        observeEvent(input[[paste0(name, "_slider")]], {
          updateNumericInput(session, paste0(name, "_numeric"), value = input[[paste0(name, "_slider")]])
        })
        observeEvent(input[[paste0(name, "_numeric")]], {
          updateSliderInput(session, paste0(name, "_slider"), value = input[[paste0(name, "_numeric")]])
        })
      })
    }
  })
  
  ## Flagged cell calculation ----
  output$flaggedTable <- renderDT({
    cell_stats <- readRDS("Rdata/flagged_metadata_raw.rds")
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
    # Load data once.
    cell_stats <- readRDS("Rdata/flagged_metadata_raw.rds")
    
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
    cell_stats <- readRDS("Rdata/flagged_metadata_raw.rds") # Load data once.
    
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
  
  
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Application --------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
shinyApp(ui, server)