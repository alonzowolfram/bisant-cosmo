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

options(scipen = 6, digits = 4) # View outputs in non-scientific notation.
# memory.limit(30000000)     # this is needed on some PCs to increase memory allowance, but has no impact on Macs.

# Prevent regexPipes functions from masking base functions.
# https://stackoverflow.com/a/5564654
# It's crazy how many of our problems stem from that lol. 8/
grep <- base::grep
grepl <- base::grepl
gsub <- base::gsub

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Functions ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load required functions.
library(tidyverse) # For a data-science focused data "grammar".
## Function to add a slash to a directory path if it doesn't have one already.
appendSlashToPath <- function(x) {
  ifelse(base::grepl("\\/$", x), x, paste0(x, "/"))
}
## Function to check if variable is NULL or empty.
flagVariable <- function(x) {
  return(is.null(x) || x=="")
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

## Load the raw variables.
## Project
### Meta
meta <- config$project$meta
project_name <- meta$project_name
run_name <- meta$run_name
### Technical
path_to_python <- config$project$technical$path_to_python
## Data
data <- config$data
flat_file_dir <- data$flat_file_dir
run_summaries_dir <- data$run_summaries_dir
bacterial_probes <- data$bacterial_probes
rmd_template_file <- data$rmd_template_file

## Outputs
# (output_dir read in from command line arguments)

## Experiment
experiment <- config$experiment
### General
general <- experiment$general
species <- general$species
random_seed <- general$random_seed
### Annotation
annotation <- experiment$annotation
slide_name_var <- annotation$slide_name_var
fov_name_var <- annotation$fov_name_var
dimension_name_vars <- annotation$dimension_name_vars
total_counts_var <- annotation$total_counts_var
total_features_var <- annotation$total_features_var
area_var <- annotation$area_var
tissue_var <- annotation$tissue_var
neovariables <- annotation$neovariables
filter_vars <- annotation$filter_vars
### QC: general
general_qc <- experiment$general_qc
filter_cells <- general_qc$filter_cells
filter_fovs <- general_qc$filter_fovs
filter_neg_probes <- general_qc$filter_neg_probes
filter_targets <- general_qc$filter_targets
### QC: cells
cell_qc <- experiment$cell_qc
min_counts_per_cell <- cell_qc$min_counts_per_cell
min_features_per_cell <- cell_qc$min_features_per_cell
proportion_neg_counts <- cell_qc$proportion_neg_counts
count_dist <- cell_qc$count_dist
area_outlier_pval <- cell_qc$area_outlier_pval
max_area <- cell_qc$max_area
min_signal_strength <- cell_qc$min_signal_strength
### QC: FOVs
fov_qc <- experiment$fov_qc
min_cells_per_fov <- fov_qc$min_cells_per_fov
path_to_barcodes <- fov_qc$path_to_barcodes
probe_panel <- fov_qc$probe_panel
### QC: probes
probe_qc <- experiment$probe_qc
neg_probe_outlier_pval <- probe_qc$neg_probe_outlier_pval
### QC: targets
target_qc <- experiment$target_qc
neg_control_probe_quantile_cutoff <- target_qc$neg_control_probe_quantile_cutoff
detection_over_bg_p_value <- target_qc$detection_over_bg_p_value
filter_targets_by_neg_control_quantile <- target_qc$filter_targets_by_neg_control_quantile
filter_targets_by_detection_p_value <- target_qc$filter_targets_by_detection_p_value
### Normalization
normalization <- experiment$normalization
log_transform <- normalization$log_transform
n_variable_features <- normalization$n_variable_features
### Batch-effect correction
batch_effect_correction <- experiment$batch_effect_correction
correct_batch_effects <- batch_effect_correction$correct_batch_effects
batch_covariates <- batch_effect_correction$batch_covariates

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Required parameters ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Check the required parameters passed from the configuration YAML file based on which module we're running.
current_module <- cl_args[3]

# Set up a list to hold all the error messages. 
error_msg_list <- list()
yaml_part <- " In the configuration YAML file, please provide "

### All modules ----
if(flagVariable(path_to_python)) error_msg_list[["path_to_python"]] <- paste0("[Project settings: Technical]", yaml_part, "a path to the Python binary used by the current conda environment (`cosmx`).")

### Data ----
required_annotation_elements <- c()
if(current_module=="data_import_cleaning") {
  if(flagVariable(flat_file_dir)) error_msg_list[["flat_file_dir"]] <- paste0("[Data]", yaml_part, "the absolute path to the directory containing the CosMx flat files.")
  if(flagVariable(run_summaries_dir)) error_msg_list[["run_summaries_dir"]] <- paste0("[Data]", yaml_part, "the absolute path to the directory containing the CosMx run summaries.")
  # if(list(NULL) %in% config$experiment$annotation[required_annotation_elements]) error_msg_list[["path_to_python"]] <- "In the configuration YAML file, please provide values for all experimental annotation column name settings.") # https://stackoverflow.com/questions/12119019/select-multiple-elements-from-a-list
}

### Experiment ----
#### General ----
if(flagVariable(species)) error_msg_list[["species"]] <- paste0("[Experiment: General]", yaml_part, "the species of the experimental organism. The currently allowed values are `Homo sapiens` and `Mus musculus`.")
#### Annotation ----
if(flagVariable(slide_name_var)) error_msg_list[["slide_name_var"]] <- paste0("[Experiment: Annotation]", yaml_part, "the name of the column in the metadata containing the slide name.")
if(flagVariable(fov_name_var)) error_msg_list[["fov_name_var"]] <- paste0("[Experiment: Annotation]", yaml_part, "the name of the column in the metadata containing the FOV name.")
if(flagVariable(dimension_name_vars)) error_msg_list[["dimension_name_vars"]] <- paste0("[Experiment: Annotation]", yaml_part, "the names of the column in the metadata containing the x- and y-coordinates.")
if(flagVariable(total_counts_var)) error_msg_list[["total_counts_var"]] <- paste0("[Experiment: Annotation]", yaml_part, "the names of the column in the metadata containing the total counts.")
if(flagVariable(total_features_var)) error_msg_list[["total_features_var"]] <- paste0("[Experiment: Annotation]", yaml_part, "the names of the column in the metadata containing the total features.")
if(flagVariable(area_var)) error_msg_list[["area_var"]] <- paste0("[Experiment: Annotation]", yaml_part, "the names of the column in the metadata containing the area.")
#### FOV QC ----
if(current_module=="qc") {
  if(flagVariable(path_to_barcodes)) error_msg_list[["path_to_barcodes"]] <- paste0("[Experiment: FOV QC]", yaml_part, "a path to the RDS file containing the barcodes for the FOVs.")
  if(flagVariable(probe_panel)) error_msg_list[["probe_panel"]] <- paste0("[Experiment: FOV QC]", yaml_part, "the panel of probes used in the current experiment. Can be one of `Hs_IO`, `Hs_UCC`, `Hs_6k`, `Mm_Neuro`, `Mm_UCC`")
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Optional parameters ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Here we provide default values for parameters that are required by the pipeline
# but do not necessarily have to be provided by the user in the YAML config file.

### Data ---
if(flagVariable(rmd_template_file)) rmd_template_file <- "ext/report_template.Rmd"
### Experiment ----
#### General ----
if(flagVariable(random_seed)) random_seed <- 1026
#### Annotation ----
if(flagVariable(tissue_var)) tissue_var <- slide_name_var
#### General QC ----
if(flagVariable(filter_cells)) filter_cells <- TRUE
if(flagVariable(filter_fovs)) filter_fovs <- TRUE
if(flagVariable(filter_neg_probes)) filter_neg_probes <- FALSE
if(flagVariable(filter_targets)) filter_targets <- FALSE
#### Cell QC ----
if(flagVariable(min_counts_per_cell)) min_counts_per_cell <- 200
if(flagVariable(min_features_per_cell)) min_features_per_cell <- 200
if(flagVariable(proportion_neg_counts)) proportion_neg_counts <- 0.1
if(flagVariable(count_dist)) count_dist <- 1
if(flagVariable(area_outlier_pval)) area_outlier_pval <- 0.01
if(flagVariable(max_area)) max_area <- 30000
if(flagVariable(min_signal_strength)) min_signal_strength <- 4
#### FOV QC ----
if(flagVariable(min_cells_per_fov)) min_cells_per_fov <- 50
#### Probe QC ----
if(flagVariable(neg_probe_outlier_pval)) neg_probe_outlier_pval <- 0.01
#### Target QC ----
if(flagVariable(neg_control_probe_quantile_cutoff)) neg_control_probe_quantile_cutoff <- 0.5
if(flagVariable(detection_over_bg_p_value)) detection_over_bg_p_value <- 0.01
if(flagVariable(filter_targets_by_neg_control_quantile)) filter_targets_by_neg_control_quantile <- FALSE
if(flagVariable(filter_targets_by_detection_p_value)) filter_targets_by_detection_p_value <- FALSE
#### Normalization ----
if(flagVariable(log_transform)) log_transform <- FALSE
if(flagVariable(n_variable_features)) n_variable_features <- 1000
#### Batch-effect correction ----
if(flagVariable(correct_batch_effects) || flagVariable(batch_covariates)) correct_batch_effects <- FALSE

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Parameter cleaning and formatting ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

### Data ----
if(!flagVariable(bacterial_probes)) bacterial_probes <- bacterial_probes %>% strsplit(",") %>% unlist

### Experiment ---- 
#### Annotation ----
if(!flagVariable(dimension_name_vars)) dimension_name_vars <- dimension_name_vars %>% strsplit(",") %>% unlist
#### Batch-effect correction ----
if(!flagVariable(batch_covariates)) batch_covariates <- batch_covariates %>% strsplit(",") %>% unlist

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
library(Seurat) # [conda] Handling single-cell data.
library(rmarkdown) # [conda] Rendering R Markdown reports. 

workflow_system <- cl_args[2]
if(workflow_system=="Nextflow") {
  path <- Sys.getenv("PATH") |> strsplit(":")
  bin_path <- tail(path[[1]], n=1)
  source(file.path(bin_path, "FOV_QC_utils_modified.R")) # FOV QC functions from https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/FOV-QC/code/FOV%20QC
  source(file.path(bin_path, "helper_functions.R"))
  source(file.path(bin_path, "plotting_functions.R"))
  source(file.path(bin_path, "ligand-receptor_interactions.R"))
} else {
  bin_path <- ""
  source("src/FOV_QC_utils_modified.R") # FOV QC functions from https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/FOV-QC/code/FOV%20QC
  source("src/helper_functions.R")
  source("src/plotting_functions.R")
  source("src/ligand-receptor_interactions.R")
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
  output_dir_rdata <- ""
  output_dir_tabular <- ""
  output_dir_pubs <- ""
} else {
  ## Output
  output_dir <- paste0(appendSlashToPath(cl_args[4]))
  ### Create the directory if it doesn't already exist. 
  if(!dir.exists(output_dir)) dir.create(output_dir)
  ### Create the folder structure within the output_dir.
  for(subdir in c("config", "logs", "pubs", "Rdata", "tabular")) {
    subdir_path <- file.path(output_dir, subdir)
    if(!dir.exists(subdir_path)) dir.create(subdir_path)
  }
  output_dir_config <- paste0(output_dir, "config/")
  output_dir_logs <- paste0(output_dir, "logs/")
  output_dir_rdata <- paste0(output_dir, "Rdata/")
  output_dir_tabular <- paste0(output_dir, "tabular/")
  output_dir_pubs <- paste0(output_dir, "pubs/")
}

rdata_folder <- ifelse(workflow_system=="Nextflow", "", "Rdata/")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Parameter check ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# If any of the required parameters are missing, 
# print a message for each one, then stop the pipeline.
if(length(error_msg_list) > 0) {
    message("Error: you are missing one or more required parameters. Please see the error messages below.")
  
  for(msg in error_msg_list) {
    message(msg)
  }
  
  stop("Shutting down pipeline due to missing parameters.")
}