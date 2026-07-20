# SurveilAMR — R environment bootstrap
# Installs every CRAN package required by r/surveilamr_analysis.R
# Run once with: Rscript r/install_packages.R

required_packages <- c(
  "readr",       # fast/chunked CSV reading for the 387 MB ATLAS extract
  "readxl",      # KEYSTONE / DREAM / GASAR .xlsx readers
  "dplyr",       # data wrangling
  "tidyr",       # reshaping (pivot_longer/wider)
  "stringr",     # string cleaning (MIC parsing, subtype normalization)
  "ggplot2",     # figures
  "scales",      # axis/percentage formatting
  "forcats",     # factor releveling for ordered plots
  "jsonlite"     # writing dataset_summary.json equivalents
)

installed <- rownames(installed.packages())
missing <- setdiff(required_packages, installed)

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}
