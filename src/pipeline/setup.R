## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##   License ----
##
#     bisant-cosmo is a Nextflow pipeline for processing, running QC on, and analyzing NanoString CosMx spatial single-cell data.
#     Copyright (C) 2025 Lon Fong.

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
##                                                                
## R settings ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

options(scipen = 6, digits = 4) # View outputs in non-scientific notation
# memory.limit(30000000)     # this is needed on some PCs to increase memory allowance, but has no impact on Macs

# Prevent regexPipes functions from masking base functions
# https://stackoverflow.com/a/5564654
# It's crazy how many of our problems stem from that lol 8/
grep <- base::grep
grepl <- base::grepl
gsub <- base::gsub

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Functions ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load required functions
library(tidyverse) # Data-science focused data "grammar".
## Function to add a slash to a directory path if it doesn't have one already
appendSlashToPath <- function(x) {
  ifelse(base::grepl("\\/$", x), x, paste0(x, "/"))
}
## Function to check if variable is NULL or empty.
flagVariable <- function(x) {
  return(is.null(x) || sum(x=="") >= length(x))
}

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

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Parameters from configuration YAML file ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Set the parameters passed from the configuration YAML file.
## Read in the config file. 
cl_args <- commandArgs(TRUE)
library(yaml) # For reading in YAML documents.
config <- yaml::read_yaml(cl_args[1])

# Check the required parameters passed from the configuration YAML file based on which module we're running.
current_module <- cl_args[3]

## Load the raw variables.
config_env <- assignVarsEnv(config)
## Process the variables.
config_metadata_path <- "config_variable_metadata.csv"
error_msg_list <- validateProcessConfig(config_metadata_path)

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Environments ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# # Python environment with leidenalg needed for clustering.
# [2024/09/24] Apparently this part needs to go first, otherwise another version of python will be used?
Sys.setenv(RETICULATE_PYTHON = path_to_python)
reticulate::use_condaenv("cosmx", required = TRUE)
reticulate::py_config()
system("which python")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Required libraries and functions ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# if(Sys.info()['sysname']=="Linux" & base::grep("dragon", Sys.info()['nodename'])) .libPaths("/home/lwfong/R/ubuntu/4.3.1")
library(furrr) # Parallelization. Must be loaded before `Seurat` to prevent conflicts.
library(Seurat) # [conda] Handling single-cell data
library(rmarkdown) # [conda] Rendering R Markdown reports
library(here) # [conda] 
library(ggraph) # Co-occurrence graphs
library(harmony) # [conda] Batch-effect correction
library(igraph) # Co-occurrence graphs
library(InSituType) # [manual] NanoString's in-house cell-typing package
library(pheatmap) # [conda] Heatmaps (for post-cell-typing flight maps)
library(mclust) # Clustering of PCs in spatial niche module.
library(ggrastr) # Rasterization of ggplot layers for smaller memory sizes.
library(glmmTMB) # Negative binomial models.
library(emmeans) # Pairwise comparisons after fitting non-spatial negative binomial models.
library(smiDE) # [manual] Spatial differential expression pre-processing.
library(mgcv) # Generalized Additive Modeling for spatial effects; replaces INLA because of computational costs.
library(readr)  # Writing logs
library(metafor) # Meta-analysis for combining differential expression results from different runs. 
library(msigdbr) # Offline access to mSigDB's pathways.
library(msigdbdf) # [manual] Helper package for msigdbr.
library(fgsea) # Fast GSEA for pathway analysis.
library(gridExtra) # Arranging plots into grids via `grid.arrange`
library(spatialTIME) # Spatial analysis.
# library(CellChat) # Cell-cell interaction analysis

workflow_system <- cl_args[2]
if(workflow_system=="Nextflow") {
  path <- Sys.getenv("PATH") |> strsplit(":")
  bin_path <- tail(path[[1]], n=1)
  source(file.path(bin_path, "helper_functions.R"))
  source(file.path(bin_path, "plotting_functions.R"))
  source(file.path(bin_path, "FOV_QC_utils_modified.R")) # FOV QC functions from https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/FOV-QC/code/FOV%20QC
  source(file.path(bin_path, "InSituCor_NeighborhoodCalculations.R")) # Functions from InSituCor that for some reason are not exported from the package ...
  source(file.path(bin_path, "spatial_functions.R"))
  source(file.path(bin_path, "differential-expression_functions.R"))
  source(file.path(bin_path, "ligand-receptor_interactions.R"))
} else {
  bin_path <- ""
  source("src/functions/helper_functions.R")
  source("src/functions/plotting_functions.R")
  source("src/functions/FOV_QC_utils_modified.R") # FOV QC functions from https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/FOV-QC/code/FOV%20QC
  source("src/functions/InSituCor_NeighborhoodCalculations.R") # Functions from InSituCor that for some reason are not exported from the package ...
  source("src/functions/spatial_functions.R")
  source("src/functions/differential-expression_functions.R")
  source("src/functions/ligand-receptor-interaction_functions.R")
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
##                                                                
## File/path settings ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if(workflow_system == "Nextflow") {
  output_dir <- ""
  output_dir_config <- ""
  output_dir_logs <- ""
  output_dir_logs_ct <- ""
  output_dir_rdata <- ""
  output_dir_tabular <- ""
  output_dir_pubs <- ""
} else {
  ## Output
  output_dir <- paste0(appendSlashToPath(cl_args[4]))
  ### Create the directory if it doesn't already exist
  if(!dir.exists(output_dir)) dir.create(output_dir)
  ### Create the folder structure within the output_dir
  for(subdir in c("config", "logs", "pubs", "Rdata", "tabular")) {
    subdir_path <- file.path(output_dir, subdir)
    if(!dir.exists(subdir_path)) dir.create(subdir_path)
  }
  output_dir_config <- paste0(output_dir, "config/")
  output_dir_logs <- paste0(output_dir, "logs/")
  output_dir_rdata <- paste0(output_dir, "Rdata/")
  output_dir_tabular <- paste0(output_dir, "tabular/")
  output_dir_pubs <- paste0(output_dir, "pubs/")
  
  ### Create the subdir in `logs/` for cell typing
  output_dir_logs_ct <- paste0(output_dir_logs, "cell_typing/")
  if(!dir.exists(output_dir_logs_ct)) dir.create(output_dir_logs_ct)
}

rdata_folder <- ifelse(workflow_system=="Nextflow", "", "Rdata/")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 
## Parameter check ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# If any of the required parameters are missing, 
# print a message for each one, then stop the pipeline
if(length(error_msg_list) > 0) {
    message("Error: you are missing one or more required parameters. Please see the error messages below.")
  
  for(msg in error_msg_list) {
    message(msg)
  }
  
  stop("Shutting down pipeline due to missing parameters.")
}