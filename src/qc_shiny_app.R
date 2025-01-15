library(shiny)
library(tidyverse)
library(DT) # For the display tables.

# Data frame containing all the QC parameters.
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
                  name = c("min_counts_per_cell","min_features_per_cell","max_proportion_neg_counts","min_count_distribution","max_area","min_signal_strength"),
                  displayname = c("Minimum number of counts per cell", "Minimum number of features per cell", "Maximum proportion of negative counts per cell", "Minimum complexity per cell", "Maximum area per cell", "Minimum signal strength per cell"),
                  metadata_name = c("nCount_RNA", "nFeature_RNA", "ProportionNegative", "Complexity", "Area", "SignalStrength"),
                  min = c(0, 0, 0, 0, 0, 0),
                  max = c(1000, 1000, 1, 200, 50000, 100),
                  default = c(200, 200, 0.1, 1, 30000, 4),
                  step = c(1, 1, 0.01, 1, 50, 1),
                  flag_as_minimum = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
                )
)

ui <- fluidPage(
  # Include custom CSS for table styling
  tags$head(
    tags$style(HTML("
      .dt-center { text-align: center; }
      .dt-left { text-align: left; }
    "))
  ),
  
  # Application title
  titlePanel("QC"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      # Dynamically generated inputs
      uiOutput("dynamicInputs"),
      
      # Display the values of the inputs
      verbatimTextOutput("inputValues")
    ),
    
    # Display the table of flagged cells
    mainPanel(
      DTOutput("flaggedTable")
    )
  )
)

server <- function(input, output, session) {
  # Dynamically generate the inputs
  output$dynamicInputs <- renderUI({
    # Create a list of UI elements
    uiList <- lapply(1:nrow(config), function(i) {
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
    
    # # Return the list of UI elements
    # do.call(tagList, uiList)
  })
  
  # Synchronize slider and numeric inputs
  observe({
    for (i in 1:nrow(config)) {
      local({
        name <- config$name[i]
        
        observeEvent(input[[paste0(name, "_slider")]], {
          updateNumericInput(
            session,
            paste0(name, "_numeric"),
            value = input[[paste0(name, "_slider")]]
          )
        })
        
        observeEvent(input[[paste0(name, "_numeric")]], {
          updateSliderInput(
            session,
            paste0(name, "_slider"),
            value = input[[paste0(name, "_numeric")]]
          )
        })
      })
    }
  })
  
  # Calculate flagged cells dynamically
  output$flaggedTable <- renderDT({
    # Load in the data object.
    cell_stats <- readRDS("Rdata/flagged_metadata_raw.rds")
    
    # Prepare an empty list to store flagged counts
    flagged_matrix <- sapply(1:nrow(config), function(i) {
      stat <- config$name[i]
      metadata_name <- config$metadata_name[i]
      cutoff <- input[[paste0(stat, "_numeric")]]
      is_minimum <- config$flag_as_minimum[i]
      
      # Determine flagged cells based on min/max logic
      if (is_minimum) {
        # sum(cell_stats[[metadata_name]] < cutoff) # Count cells below the minimum
        cell_stats[[metadata_name]] < cutoff
      } else {
        # sum(cell_stats[[metadata_name]] > cutoff) # Count cells above the maximum
        cell_stats[[metadata_name]] > cutoff
      }
    })
    
    # Calculate flagged counts for individual stats
    flagged_counts <- colSums(flagged_matrix)
    
    # Calculate total flagged cells (any stat)
    total_flagged <- sum(rowSums(flagged_matrix) > 0)
    total_passed <- nrow(cell_stats) - total_flagged
    
    # Create a data frame for display
    results <- data.frame(
      Stat = c(config$displayname, "Total Flagged", "Total Passed", "All cells"),
      `Flagged count` = c(flagged_counts, total_flagged, total_passed, nrow(cell_stats))
    )
    
    datatable(results, options = list(
      pageLength = 10,
      columnDefs = list(
        list(width = '40%', targets = 0), # Adjust width for Stat column
        list(width = '60%', targets = 1)
      )
    ))
  })
}

# Run the application.
shinyApp(ui, server)