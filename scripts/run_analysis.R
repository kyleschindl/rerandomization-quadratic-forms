#!/usr/bin/env Rscript
# scripts/run_analysis.R
# Master script to reproduce all figures and tables for the paper

# 1. Load here for relative paths
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

# 2. Activate renv environment for exact package versions
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::activate()

# 3. Source helper functions
source(here("R", "replication_functions.R"))

# 4. Create output directories if missing
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"),  recursive = TRUE, showWarnings = FALSE)

# 5. Generate all figures
source(here("scripts", "generate_figures.R"))

# 6. Generate all tables
source(here("scripts", "generate_tables.R"))