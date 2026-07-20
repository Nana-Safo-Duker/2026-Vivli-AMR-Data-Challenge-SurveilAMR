#!/usr/bin/env Rscript
# =============================================================================
# SurveilAMR — R analysis pipeline
#
# Reproduces the core Python pipeline (scripts/run_analysis.py,
# scripts/analyze_supplementary.py, scripts/analyze_spidaar.py) in R for
# analysts who prefer the tidyverse / ggplot2 ecosystem.
#
# Covers all 5 approved Vivli AMR Register datasets:
#   1. ATLAS_Antibiotics (Pfizer)            — data/raw/atlas_vivli_2004_2024.csv
#   2. KEYSTONE / Omadacycline (Paratek)      — data/raw/Omadacycline_2015 to 2025_Surveillance_data.xlsx
#   3. Bedaquiline DREAM (Johnson & Johnson)  — data/raw/BEDAQUILINE DREAM DATASET FOR VIVLI - 06-06-2022.xlsx
#   4. GASAR Study III (Venus Remedies)       — data/raw/GASAR Study III (n=494)_updated.xlsx
#   5. SPIDAAR RWE Study (Pfizer)             — data/raw/spidaar_isolatedata.xls, spidaar_patientdata.xls
#
# Outputs are written to data/processed/ (CSV, prefixed "r_") and
# outputs/figures/ (PNG, prefixed "r_") so they never collide with the
# Python-generated artifacts already tracked in the repository.
#
# Usage:
#   Rscript r/install_packages.R      # once, to install dependencies
#   Rscript r/surveilamr_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  NULL
}

script_dir <- get_script_dir()
root <- if (!is.null(script_dir)) normalizePath(file.path(script_dir, "..")) else getwd()

# Fall back to working directory if the resolved path looks wrong
# (e.g. when pasted into an interactive R session run from the repo root).
if (is.na(root) || !dir.exists(file.path(root, "data"))) {
  root <- getwd()
}

raw_dir  <- file.path(root, "data", "raw")
proc_dir <- file.path(root, "data", "processed")
fig_dir  <- file.path(root, "outputs", "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 11))

parse_mic <- function(x) {
  x <- str_replace_all(as.character(x), "[\u2264\u2265<>]", "")
  suppressWarnings(as.numeric(str_trim(x)))
}

# -----------------------------------------------------------------------------
# 1. ATLAS_Antibiotics — chunked read of the ~387 MB extract
# -----------------------------------------------------------------------------
run_atlas <- function() {
  atlas_path <- file.path(raw_dir, "atlas_vivli_2004_2024.csv")
  if (!file.exists(atlas_path)) {
    message("ATLAS raw file not found at ", atlas_path, " — skipping ATLAS section.")
    return(invisible(NULL))
  }

  priority_species <- c(
    "Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa",
    "Staphylococcus aureus", "Acinetobacter baumannii", "Enterococcus faecalis"
  )
  key_abx <- list(
    "Escherichia coli" = "Meropenem",
    "Klebsiella pneumoniae" = "Meropenem",
    "Pseudomonas aeruginosa" = "Meropenem",
    "Staphylococcus aureus" = "Oxacillin",
    "Acinetobacter baumannii" = "Meropenem",
    "Enterococcus faecalis" = "Vancomycin"
  )

  species_counts  <- list()
  country_counts  <- list()
  yearly_acc      <- list()
  n_rows          <- 0L

  bump_counts <- function(counts, keys) {
    tbl <- table(keys)
    for (nm in names(tbl)) {
      counts[[nm]] <- (counts[[nm]] %||% 0L) + as.integer(tbl[[nm]])
    }
    counts
  }

  process_chunk <- function(chunk, pos) {
    n_rows <<- n_rows + nrow(chunk)
    species_counts <<- bump_counts(species_counts, chunk$Species)
    country_counts <<- bump_counts(country_counts, chunk$Country)

    for (sp in priority_species) {
      ab_col <- paste0(key_abx[[sp]], "_I")
      if (!ab_col %in% names(chunk)) next
      sub <- chunk %>% filter(Species == sp, !is.na(.data[[ab_col]]))
      if (nrow(sub) == 0) next
      agg <- sub %>%
        group_by(Year) %>%
        summarise(n = n(), res = sum(.data[[ab_col]] == "Resistant"), .groups = "drop") %>%
        mutate(Species = sp, Antibiotic = key_abx[[sp]])
      key <- sp
      yearly_acc[[key]] <<- bind_rows(yearly_acc[[key]], agg)
    }
  }

  read_csv_chunked(
    atlas_path,
    callback = SideEffectChunkCallback$new(process_chunk),
    chunk_size = 100000,
    col_types = cols(.default = col_character())
  )

  yearly <- bind_rows(yearly_acc) %>%
    group_by(Species, Antibiotic, Year) %>%
    summarise(Total_Isolates = sum(n), Resistant = sum(res), .groups = "drop") %>%
    filter(Total_Isolates >= 30) %>%
    mutate(Resistance_Pct = round(100 * Resistant / Total_Isolates, 2)) %>%
    mutate(Year = as.integer(Year)) %>%
    arrange(Species, Antibiotic, Year)

  write_csv(yearly, file.path(proc_dir, "r_surveilamr_resistance_by_year.csv"))

  summary_list <- list(
    total_rows = n_rows,
    n_countries = length(country_counts),
    n_species = length(species_counts)
  )
  write_json(summary_list, file.path(proc_dir, "r_dataset_summary.json"), pretty = TRUE, auto_unbox = TRUE)

  p <- yearly %>%
    filter(Species %in% c("Escherichia coli", "Klebsiella pneumoniae", "Acinetobacter baumannii")) %>%
    ggplot(aes(Year, Resistance_Pct, color = Species)) +
    geom_line(linewidth = 1) +
    geom_point() +
    labs(
      title = "ATLAS: Primary-Antibiotic Resistance Trend (R pipeline)",
      y = "Resistance (%)", x = "Year"
    ) +
    scale_y_continuous(labels = label_percent(scale = 1))
  ggsave(file.path(fig_dir, "r_fig_atlas_resistance_trend.png"), p, width = 9, height = 5, dpi = 150)

  message(sprintf("ATLAS: processed %s isolates across %d countries.", format(n_rows, big.mark = ","), length(country_counts)))
}

# -----------------------------------------------------------------------------
# 2. KEYSTONE — Omadacycline surveillance (Paratek)
# -----------------------------------------------------------------------------
run_keystone <- function() {
  path <- file.path(raw_dir, "Omadacycline_2015 to 2025_Surveillance_data.xlsx")
  if (!file.exists(path)) return(invisible(NULL))

  df <- read_excel(path, sheet = "Line List")
  names(df) <- str_squish(str_replace_all(names(df), "\n", " "))
  df <- df %>% mutate(Oma_MIC = parse_mic(Omadacycline))

  priority <- c(
    "Staphylococcus aureus", "Escherichia coli", "Klebsiella pneumoniae",
    "Streptococcus pneumoniae", "Enterococcus faecalis", "Pseudomonas aeruginosa"
  )

  mic_summary <- df %>%
    filter(Organism %in% priority, !is.na(Oma_MIC)) %>%
    group_by(Organism) %>%
    summarise(
      n = n(),
      MIC50 = quantile(Oma_MIC, 0.5, na.rm = TRUE),
      MIC90 = quantile(Oma_MIC, 0.9, na.rm = TRUE),
      mean = round(mean(Oma_MIC, na.rm = TRUE), 4),
      .groups = "drop"
    )
  write_csv(mic_summary, file.path(proc_dir, "r_keystone_omadacycline_mic.csv"))

  p <- mic_summary %>%
    ggplot(aes(x = reorder(Organism, MIC90), y = MIC90)) +
    geom_col(fill = "#d35400") +
    coord_flip() +
    labs(title = "KEYSTONE: Omadacycline MIC90 by Pathogen (R pipeline)", x = NULL, y = "MIC90 (\u00b5g/mL)")
  ggsave(file.path(fig_dir, "r_fig_keystone_mic90.png"), p, width = 8, height = 5, dpi = 150)

  message(sprintf("KEYSTONE: n = %d isolates, %d species, %d countries.", nrow(df), n_distinct(df$Organism), n_distinct(df$Country)))
}

# -----------------------------------------------------------------------------
# 3. DREAM — Bedaquiline MDR-TB surveillance (Johnson & Johnson)
# -----------------------------------------------------------------------------
run_dream <- function() {
  path <- file.path(raw_dir, "BEDAQUILINE DREAM DATASET FOR VIVLI - 06-06-2022.xlsx")
  if (!file.exists(path)) return(invisible(NULL))

  df <- read_excel(path, sheet = "DREAM Dataset") %>%
    mutate(
      Continent = str_to_title(str_trim(as.character(Continent))),
      BDQ_MIC = parse_mic(`BDQ Broth`)
    )

  by_continent <- df %>%
    group_by(Continent) %>%
    summarise(n = n(), median_bdq = median(BDQ_MIC, na.rm = TRUE), mean_bdq = mean(BDQ_MIC, na.rm = TRUE), .groups = "drop")
  write_csv(by_continent, file.path(proc_dir, "r_dream_bdq_by_continent.csv"))

  by_year <- df %>%
    group_by(`Year Collected`) %>%
    summarise(n = n(), median_bdq = median(BDQ_MIC, na.rm = TRUE), .groups = "drop")
  write_csv(by_year, file.path(proc_dir, "r_dream_bdq_by_year.csv"))

  p <- ggplot(by_year, aes(`Year Collected`, median_bdq)) +
    geom_line(color = "#c0392b", linewidth = 1) +
    geom_point() +
    labs(title = "DREAM: Median Bedaquiline MIC Over Time (R pipeline)", y = "Median BDQ MIC (\u00b5g/mL)", x = "Year")
  ggsave(file.path(fig_dir, "r_fig_dream_bdq_trend.png"), p, width = 8, height = 5, dpi = 150)

  message(sprintf("DREAM: n = %d isolates, %d countries, median BDQ MIC = %.3f.", nrow(df), n_distinct(df$Country), median(df$BDQ_MIC, na.rm = TRUE)))
}

# -----------------------------------------------------------------------------
# 4. GASAR Study III — Gram-negative mechanisms (Venus Remedies)
# -----------------------------------------------------------------------------
run_gasar <- function() {
  path <- file.path(raw_dir, "GASAR Study III (n=494)_updated.xlsx")
  if (!file.exists(path)) return(invisible(NULL))

  df <- read_excel(path, sheet = "Sheet1") %>%
    mutate(Poly_MIC = suppressWarnings(as.numeric(`Polymyxin B MIC (mcg/ml)`)))

  phenotype_bucket <- function(x) {
    s <- str_to_upper(as.character(x))
    case_when(
      str_detect(s, "NON ESBL") & str_detect(s, "NON MBL") ~ "Non-ESBL/Non-MBL",
      str_detect(s, "ESBL") & str_detect(s, "MBL") ~ "ESBL+MBL",
      str_detect(s, "MBL") ~ "MBL",
      str_detect(s, "ESBL") ~ "ESBL",
      str_detect(s, "CARBAPENEMASE") ~ "Carbapenemase",
      TRUE ~ "Other"
    )
  }
  df <- df %>% mutate(Phenotype_Bucket = phenotype_bucket(`Phenotypic Combination`))

  pheno_counts <- df %>% count(Phenotype_Bucket, name = "Count") %>% arrange(desc(Count))
  write_csv(pheno_counts, file.path(proc_dir, "r_gasar_phenotype_counts.csv"))

  p <- ggplot(pheno_counts, aes(x = reorder(Phenotype_Bucket, Count), y = Count)) +
    geom_col(fill = "#8e44ad") +
    coord_flip() +
    labs(title = "GASAR: Phenotypic Resistance Classes (R pipeline)", x = NULL, y = "Isolates")
  ggsave(file.path(fig_dir, "r_fig_gasar_phenotypes.png"), p, width = 7, height = 5, dpi = 150)

  message(sprintf("GASAR: n = %d isolates, MBL/carbapenemase = %.1f%%.", nrow(df), 100 * mean(df$Phenotype_Bucket %in% c("MBL", "ESBL+MBL", "Carbapenemase"))))
}

# -----------------------------------------------------------------------------
# 5. SPIDAAR RWE Study — sub-Saharan Africa isolate + patient data (Pfizer)
# -----------------------------------------------------------------------------
run_spidaar <- function() {
  isolate_path <- file.path(raw_dir, "spidaar_isolatedata.xls")
  patient_path <- file.path(raw_dir, "spidaar_patientdata.xls")
  if (!file.exists(isolate_path) || !file.exists(patient_path)) return(invisible(NULL))

  group_labels <- c(
    "1" = "Staphylococcus aureus", "2" = "Streptococcus", "3" = "Enterobacterales",
    "4" = "Pseudomonads", "5" = "Acinetobacters", "6" = "Enterococcus",
    "7" = "Other Staphylococci", "8" = "Other Gram-negative", "9" = "Other Gram-positive"
  )

  isolates <- read_excel(isolate_path, sheet = "data") %>%
    mutate(Group_Label = recode(as.character(group), !!!group_labels, .default = "Unclassified"))

  group_by_country <- isolates %>% count(ctry, Group_Label, name = "Count")
  write_csv(group_by_country, file.path(proc_dir, "r_spidaar_isolate_groups_by_country.csv"))

  mdr_pct <- isolates %>% filter(mdr != 99) %>% summarise(pct = round(100 * mean(mdr == 1), 2)) %>% pull(pct)

  p <- group_by_country %>%
    ggplot(aes(x = ctry, y = Count, fill = Group_Label)) +
    geom_col(position = "stack") +
    labs(title = "SPIDAAR: Isolate Organism Groups by Country (R pipeline)", x = "Country", y = "Isolates", fill = NULL)
  ggsave(file.path(fig_dir, "r_fig_spidaar_isolates.png"), p, width = 9, height = 5, dpi = 150)

  ward_labels <- c("1" = "ICU", "2" = "Internal medicine", "3" = "Surgery wards", "4" = "Other wards")
  patients <- read_excel(patient_path, sheet = "data") %>%
    mutate(Ward_Label = recode(as.character(ward), !!!ward_labels, .default = "Unknown"))

  ward_counts <- patients %>% count(Ward_Label, name = "Count")
  write_csv(ward_counts, file.path(proc_dir, "r_spidaar_patients_by_ward.csv"))

  message(sprintf("SPIDAAR: %d isolates (MDR = %.1f%% of tested), %d patients across 4 countries.", nrow(isolates), mdr_pct, nrow(patients)))
}

# -----------------------------------------------------------------------------
main <- function() {
  message("== SurveilAMR R pipeline ==")
  run_atlas()
  run_keystone()
  run_dream()
  run_gasar()
  run_spidaar()
  message("Done. Outputs written to data/processed/ and outputs/figures/ (r_ prefix).")
}

if (identical(environment(), globalenv())) {
  main()
}
