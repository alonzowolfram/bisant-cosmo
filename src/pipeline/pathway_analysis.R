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

# Read in the DE genes table.
message("Reading in DE genes table.")
results2 <- read.csv(cl_args[5], row.names = 1, check.names = FALSE)
print(paste0("results2 dimensions: ", dim(results2)))

# If pathway_table_file is provided, check that it's a valid file.
# If not, default to hallmark, BioCarta, and Reactome pathways.
# In future iterations, we'll probably include a default pathway_table_file instead of hard-coding the pathways here.
set.default.pathways <- function() {
  cats <- c("H") # , rep("C2", 2)
  subcats <- c(NA) # , "CP:BIOCARTA", "CP:REACTOME"
  return(list(cats = cats, subcats = subcats))
}

if(!flagVariable(pathway_table_file)) { 
  if(!file.exists(pathway_table_file)) {
    msigdb_list <- set.default.pathways()
  } else {
    # Read in the file and check that it has at least one entry. 
    message("Checking provided pathway file.")
    if(base::grepl("\\.csv$", pathway_table_file)) { pathway_table <- read.csv(pathway_table_file, header = F)}
    else if(base::grepl("\\.tsv$", pathway_table_file)) { pathway_table <- read.table(pathway_table_file, header = F, sep = "\t") }
    else if(base::grepl("\\.xls.*$", pathway_table_file)) { pathway_table <- read_excel(pathway_table_file, col_names = F)} 
    else {
      warning("The pathway table provided is not in CSV, TSV, or Excel format. Check the file extension. Using default pathways instead.")
      pathway_table <- NULL
      msigdb_list <- set.default.pathways()}
  }
} else {
  msigdb_list <- set.default.pathways()
}
if(exists("pathway_table")) if(!is.null(pathway_table)) msigdb_list <- as.list(pathway_table); names(msigdb_list) <- c("cats", "subcats")

# Check that species is valid.
if(!(species %in% c("Homo sapiens", "Mus musculus"))) stop("Please provide a valid species: either 'Homo sapiens' or 'Mus musculus.'")

# Set graphical parameters.
nes_palette <- colorRampPalette(c("blue", "white", "red"))(100)

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Pathway analysis (FGSEA) ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## ................................................
##
### Data preparation ----
##
## ................................................

# Load the gene sets (from MSigDB).
# https://rpubs.com/LiYumei/806213
message("Loading gene sets.")
gene_sets <- data.frame()
for(i in 1:length(msigdb_list[["subcats"]])) {
  subcat <- msigdb_list[["subcats"]][i]
  cat <- msigdb_list[["cats"]][i]
  if(is.na(subcat)) {
    subcat <- NULL
  } else {
    subcat <- subcat
  }
  
  gene_sets <- rbind(gene_sets, 
                     msigdbr(species = species, category = cat, subcategory = subcat) %>% 
                       dplyr::distinct(gs_name, gene_symbol) %>% 
                       as.data.frame()
  )
}

# Subset to include only the pathways of interest if individual_pathways is provided.
message("Subsetting gene sets.")
if(sum(!is.null(individual_pathways)) > 0 & sum(individual_pathways != "") > 0) gene_sets <- gene_sets %>% dplyr::filter(gs_name %in% individual_pathways)

# Create an MSigDB pathway list.
message("Creating MSigDB gene set list.")
msigdbr_pathway_list <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

## ................................................
##
### Calculation and graphing ----
##
## ................................................
# Loop through DE gene table and perform FGSEA
# https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
plot_list_pathway_analysis <- list()
pathway_df <- data.frame()

message("Performing FGSEA.")
for(fixed_effect_var in unique(results2$fixed_effect)) {

  # Get all the contrasts for that fixed effect
  contrasts <- results2 %>% dplyr::filter(fixed_effect==fixed_effect_var) %>% .$contrast %>% unique  

  for(contrast_level in contrasts) {
    
    # Get all the cell types for that fixed effect × contrast
    cell_types <- results2 %>% dplyr::filter(fixed_effect==fixed_effect_var & contrast==contrast_level) %>% .$cell_type %>% unique
    
    for(cell_type_i in cell_types) {
      message(paste0("Fixed effect ", fixed_effect_var, " | contrast ", contrast_level, " | cell type ", cell_type_i))
      # Subset. 
      results2_sub <- results2 %>% dplyr::filter(
        fixed_effect==fixed_effect_var &
          contrast==contrast_level & 
          cell_type==cell_type_i
      )
      
      # Extract the DEG data--these will be the ranks we use as input to FGSEA. 
      ranks <- results2_sub$estimate
      names(ranks) <- results2_sub$gene
      ranks <- ranks %>% .[order(.)]
      
      #### Calculation ----
      message("Running FGSEA.")
      fgseaRes <- fgsea(pathways = msigdbr_pathway_list, 
                        stats    = ranks,
                        minSize  = 1, # 15
                        maxSize  = Inf) # 500
      df_sub <- fgseaRes %>% dplyr::filter(is.finite(NES))
      # Make the pathway names more readable.
      message("Cleaning FGSEA results.")
      df_sub$PathwayCleaned <- df_sub$pathway %>% 
        stringr::str_split("_") %>% 
        lapply(FUN = function(x) c(paste0(x[1], ":"), x[2:length(x)]) %>% paste(collapse=" ")) %>%  
        unlist() %>% 
        stringr::str_split(": ") %>% 
        lapply(FUN = function(x) stringr::str_to_sentence(x) %>% paste(collapse=": ")) %>% 
        unlist()
      
      # Add pathway ranking score (-log10(padj) * |NES|)
      message("Adding pathway ranking scores.")
      
      df_sub$PathwayScore <- -log10(df_sub$padj) * abs(df_sub$NES)
      # Before calculating percentiles for pathway, cull pathway list down to its final form (i.e., what will actually be graphed.)
      # If individual_pathways is set, subset to include only the pathways of interest.
      if(!is.null(individual_pathways) & sum(individual_pathways=="") < length(individual_pathways)) df_sub <- df_sub %>% dplyr::filter(pathway %in% individual_pathways)
      # If the number of pathways > limits for graphing, perform cutoff by PathwayScore.
      if(!is.null(n_max_pathways) & n_max_pathways != "") if(nrow(df_sub) > n_max_pathways) df_sub <- df_sub %>% dplyr::top_n(n = n_max_pathways, wt = PathwayScore)
      
      # Add information about the model (model number, contrast variable, current contrast.)
      message("Adding information about model.")
      
      df_sub$fixed_effect <- fixed_effect_var
      df_sub$contrast <- contrast_level
      
      # Add column indicating percentile (how high up the list a given pathway is.)
      df_sub <- df_sub %>% 
        dplyr::arrange(-PathwayScore) %>% 
        dplyr::mutate(percentile = 100 - (1:nrow(.) / nrow(.) * 100))
      df_sub$percentile <- format(round(df_sub$percentile, 1), nsmall = 1)
      
      ### Graphing ----
      message("Graphing.")
      
      upper_limit <- ifelse(max(df_sub$NES) > 0, max(df_sub$NES) %>% ceiling(), max(df_sub$NES) %>% abs %>% ceiling())
      lower_limit <- ifelse(min(df_sub$NES) < 0, min(df_sub$NES) %>% floor(), min(df_sub$NES) %>% -. %>% floor())
      upper_limit_axis <- ifelse(max(df_sub$NES) < 0, 0, upper_limit)
      lower_limit_axis <- ifelse(min(df_sub$NES) > 0, 0, lower_limit)
      dodgewidth <- position_dodge(width=0.9)
      
      # Get the contrast elements for the plot label.
      contrast_element_1 <- contrast %>% strsplit(" - ") %>% unlist %>% .[1]
      contrast_element_2 <- contrast %>% strsplit(" - ") %>% unlist %>% .[2]
      
      # Graphing parameters.
      # Bar width needs to _increase_ with the # of rows in df_sub. 0.75 is a good width for 15 rows, and 0.25 is a good width for <= 3 rows, so ... 
      # sqrt function? Or sigmoidal?
      # We'll go with a piecewise function.
      # Code base stolen from https://stackoverflow.com/a/8788595/23532435
      # Second piece of the function will be linear (y = mx + b).
      # Slope for second piece of the function: rise / run = (0.75 - 0.25) / (5 - 1) = 0.042
      # y-intercept: b = y - mx
      slope <- (0.75 - 0.25) / (5 - 1)
      y_int <- 0.75 - slope * 5
      ratio <- nrow(df_sub) / 3
      bar_width <- (ratio <= 1) * 0.25 + 
        (ratio > 1 & ratio <= 5) * (slope * ratio + y_int) + 
        (ratio > 5) * ((sqrt(ratio) - sqrt(5)) + 0.75)
      
      # Get the data frame into the proper format for graphing using geom_rect.
      df_sub_graphing <- df_sub
      df_sub_graphing$PathwayCleaned <- df_sub_graphing$PathwayCleaned %>% as.factor %>% factor(levels = df_sub %>% dplyr::arrange(NES) %>% .$PathwayCleaned) # Arrange factor levels of pathway by NES.
      df_sub_graphing$ymin <- ifelse(df_sub$NES >= 0, 0, df_sub$NES)
      df_sub_graphing$ymax <- ifelse(df_sub$NES >= 0, df_sub$NES, 0)
      df_sub_graphing$xmin <- (df_sub_graphing$PathwayCleaned %>% as.numeric) - (bar_width/2)
      df_sub_graphing$xmax <- (df_sub_graphing$PathwayCleaned %>% as.numeric) + (bar_width/2)
      
      # Graph.
      plot <- ggplot(df_sub_graphing, aes(y = NES, x = reorder(PathwayCleaned, NES), fill = NES)) +
        # Bar graph manually, using geom_rect(angle).
        geom_rect(xmin = df_sub_graphing$xmin, xmax = df_sub_graphing$xmax, ymin = df_sub_graphing$ymin, ymax = df_sub_graphing$ymax) +
        # Set the color scale for the bars.
        scale_fill_gradientn(colors = nes_palette, limits = c(lower_limit, upper_limit)) + 
        # Flip horizontal.
        coord_flip() + 
        # B&W theme.
        theme_bw() + 
        # Remove grid lines and set aspect ratio.
        theme(aspect.ratio = 1/1, # Aspect ratio needs to _decrease_ (limit 0) with the # of rows in df_sub_graphing.
              legend.position = "right",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              # axis.text.x=element_blank(), 
              # axis.ticks.x=element_blank(), 
              # axis.text.y=element_blank(), 
              axis.ticks.y=element_blank()) + 
        # Set the labels.
        labs(x = "", y = paste0(contrast_element_2, " <-> ", contrast_element_1))
      
      # Add to the list.
      plot_list_pathway_analysis[[fixed_effect_var]][[contrast_level]][[cell_type_i]] <- plot
      
      # Add df_sub to pathway_df
      pathway_df <- rbind(pathway_df, df_sub %>% dplyr::select(-leadingEdge))
      
      # Graph enrichment plots.
      message("Graphing enrichment plots.")
      for(pathway in names(msigdbr_pathway_list)) {
        plotEnrichment(msigdbr_pathway_list[[pathway]], ranks, gseaParam = 1, ticksSize = 0.2) +
          ggtitle(paste0(pathway, " | fixed effect: ", fixed_effect_var, " | contrast: ", contrast_level))
        ggsave(filename = paste0(output_dir_pubs, 
                                 paste0("fixed-effect-var_", fixed_effect_var, "_contrast_", contrast_level, "cell-type_", cell_type_i, "_pathway_", pathway, ".png") %>% regexPipes::gsub("\\/", "_")), 
               height = 9, width = 10, units = "in")
        
      } # End enrichment plots for loop
      
    } # End cell type for loop
        
  } # End contrasts for loop
  
} # End fixed_effect_var_levels for loop

## ................................................
##
### Arranging plots into grids ----
##
## ................................................

# plot_list_pathway_analysis_grid <- list()
# message("Arranging plots into grid.")
# for(fev in names(plot_list_pathway_analysis)) { # fixed-effect variable
#   for(cl in names(plot_list_pathway_analysis[[fev]])) { # contrast level
#     for(cti in names(plot_list_pathway_analysis[[fev]][[cl]])) { # cell type i
#       # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
#       # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
#       # https://stackoverflow.com/questions/78163631/r-get-legend-from-cowplot-package-no-longer-work-for-ggplot2-version-3-5-0
#       
#       p_list <- plot_list_pathway_analysis[[fev]][[cl]]
#       
#       n <- length(p_list)
#       if(n > 16) {
#         # Stop right there.
#         warning(paste0("The combination of ", fev, " - ", cl, " - ", cti, " has ", n, " plots, too many for graphing. Please graph these manually. Skipping to the next list of graphs."))
#         rm(p_list)
#         gc()
#         next
#       }
#       nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1.
#       
#       # Set the scaling factors for label and legend size.
#       sqrt_n_col <- sqrt(nCol)
#       scaling_factor <- ifelse(nCol > 1, (sqrt_n_col * nCol / 3), 0.35) # Number of rows in current grid / 3 (base number)
#       res_scaling_factor <- max(scaling_factor, 1)
#       
#       # Each plot will have a different scale, so we will not include a common legend.
#       # for(item in names(p_list)) {p_list[[item]] <- p_list[[item]] + theme(legend.position = "none")}
#       
#       # Get the model information (test/contrast variable, etc.)
#       fixed_effect_var <- fev
#       contrast_level <- cl
#       cell_type_i <- cti
#       
#       # The following code deals with subset variables, and, as for now, is not used.
#       # This may change in the future, though.
#       # if(subset_var_var=="DummySubsetVar") {
#       #   subset_by <- ""
#       #   filename_subset_by <- ""
#       # } else {
#       #   subset_by <- paste0("| Subset variable: ", subset_var_var, ", level: ", subset_var_var_level)
#       #   filename_subset_by <- paste0("subset-variable-", subset_var_var, "_level-", subset_var_var_level, "_")
#       # }
#       
#       # Arrange plots in p_list onto a grid.
#       plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
#       plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Pathway", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
#                                                          bottom = grid::textGrob("Normalized Enrichment Score", gp = grid::gpar(cex = 1.3)),
#                                                          top = grid::textGrob(paste0("Pathway analysis ", 
#                                                                                      "\nFixed-effect variable: ", fixed_effect_var, " | contrast: ", contrast_level,
#                                                                                      "\nCell type: ", cell_type_i)))
#       
#       # Save to list.
#       plot_list_pathway_analysis_grid[[fev]][[cl]] <- plot_grid#2 %>% ggpubr::as_ggplot()
#       
#     }
#   }
# }

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Exporting to disk.")

# Export FGSEA results as table.
pathway_df %>% write.csv(paste0(output_dir_tabular, "pathway-analysis_results.csv"))
# Export the raw plots as RDS file.
plot_list_pathway_analysis %>% saveRDS(paste0(output_dir_rdata, "pathway-analysis_raw-plots-list.rds"))
# Export the grid-arranged plots as RDS file.
# plot_list_pathway_analysis_grid %>% saveRDS(paste0(output_dir_rdata, "pathway-analysis_grid-arranged-plots-list.rds"))