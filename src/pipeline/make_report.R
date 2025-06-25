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

# Load the RDS objects with the necessary data (which is just the latest module completed). 
rds_path <- cl_args[5]
latest_module <- readRDS(rds_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Make the report -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message(paste0("Rendering HTML report using the template found at ", rmd_template_file, "."))
# saveRDS(iris, paste0(output_dir_pubs, "test.rds"))
# Set the paths.
qc_metrics_path <- paste0(output_dir_rdata, "qc_metrics.rds")
bac_probe_stats_path <- paste0(output_dir_rdata, "bac-probe_stats.rds")
# Render.
rmarkdown::render(rmd_template_file,
                  output_file=paste0(output_dir_pubs, "report.html"),
                  params=list(qc_metrics_file=qc_metrics_path,
                              bac_probe_stats_file=bac_probe_stats_path))