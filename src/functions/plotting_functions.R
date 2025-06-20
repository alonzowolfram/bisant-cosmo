#' Export plot to PNG and EPS.
#' 
#' @param plot ggplot2 plot
#' @export
DSPexportPlot <- function(plot, path_out = NULL, filename = NULL, 
                          dpi = 600, png_width = NULL, png_height = NULL, png_unit = "mm") {
  # Export to PNG and EPS.
  if(is.null(path_out)) path_out <- "./"
  if(is.null(filename)) filename <- "plot"
  # EPS
  ggsave(filename = paste0(filename, ".eps"), path = path_out)
  
  # PNG
  if(is.null(png_width) || is.null(png_height)) {
    ggsave(filename = paste0(filename, ".png"),
           path = path_out,
           dpi = dpi)
  } else {
    ggsave(filename = paste0(filename, ".png"), 
           path = path_out,
           dpi = dpi, 
           width = png_width, height = png_height, unit = png_unit)
  }
}

#' Make a scatterplot.
#' 
#' @param plot ggplot2 plot
#' @param point_size Size of points (in pt) in scatterplot. 
#' @export
DSPplotScatter <- function(plot, 
                           point_size = 3, point_alpha = 1, geom_jitter = FALSE, 
                           panel_border_color = "black", panel_border_linewidth = 1,
                           facet_wrap_var = NULL,
                           coord_equal = FALSE,
                           legend_position = "right",
                           x_lab = NULL, y_lab = NULL, title = NULL) {
  # Extract the current labels from the ggplot object.
  current_labels <- plot$labels
  # See if the x and y labels are supplied.
  if(is.null(x_lab)) x_lab <- current_labels$x
  if(is.null(y_lab)) y_lab <- current_labels$y
  
  # Draw the plot.
  plot <- plot + 
    geom_point(size = point_size, alpha = point_alpha) + # Make it a scatter plot.
    theme_bw() + # Set the background theme.                                                 
    theme(panel.grid.major = element_blank(), # Remove the major lines.
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = panel_border_color, fill = NA, linewidth = panel_border_linewidth)) + # Set the panel border. https://www.statology.org/ggplot2-panel-border/
    labs(x = x_lab, y = y_lab, title = title) # Remove the minor lines.
  
  # Jitter
  if(geom_jitter) plot <- plot + geom_jitter()
  # Facet wrap
  if(!is.null(facet_wrap_var)) plot <- plot + facet_wrap(~ get(facet_wrap_var))
  # coord_equal
  if(coord_equal) plot <- plot + coord_equal()
  
  # Return plot.
  return(plot)

}

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

#' Make a volcano plot.
#' 
#' @param plot ggplot2 plot
#' @param point_size Size of points (in pt). 
#' @export