#!/usr/bin/env Rscript

# Load necessary libraries
library(yaml)

# Get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
master_yaml <- args[1] # Template YAML file
project_yaml_list <- args[2] # Vector of YAML files to be configured according to the master YAML file

trimTrailingWhitespace <- function(x) {
  sub("[ \t]+$", "", x)
}
startsWithHash <- function(x) {
  grepl("^[ \t]*#", x)
}
emptyOrWhitespace <- function(x) {
  grepl("^\\s*$", x)
  
  # # Example usage:
  # lines <- c("   ", "", "This has content", "\t\t")
  # sapply(lines, emptyOrWhitespace)
}

mergeConfig <- function(master_file, project_file) {
  # Read config files as raw text.
  master_lines <- readLines(master_file, warn = FALSE) %>% sapply(trimTrailingWhitespace)
  project_lines <- if (file.exists(project_file)) readLines(project_file, warn = FALSE) %>% sapply(trimTrailingWhitespace) else NULL
  project_lines_non_comment <- project_lines %>% .[!startsWithHash(.)]
  
  # Process project file line by line.
  updated_lines <- c()
  for (line in master_lines) {
    # Check if it's a comment or empty/whitespace line. 
    if(startsWithHash(line) || emptyOrWhitespace(line)) {
      updated_lines <- c(updated_lines, line)
    } else {
      # Check if the variable already exists in the copy (project_file).
      # First extract the key only, not the value.
      key <- line %>% regexPipes::gsub(":.+$", "")
      # Now look for this key in project_lines.
      # Make sure to check only the non-comment lines. 
      index <- grep(key, project_lines_non_comment)
      
      if(length(index) < 1) {
        # It doesn't exist. 
        # So add the line from the master file to `updated_lines`.
        updated_lines <- c(updated_lines, line)
      } else {
        # It does.
        # So add the line from the project file to `updated_lines`.
        updated_lines <- c(updated_lines, project_lines_non_comment[index[1]])
      }
      
    }
  }
  
  # Write the updated file.
  writeLines(updated_lines, project_file)
}

for (project_yaml in project_yaml_list) {
  mergeConfig(master_yaml, project_yaml)
}