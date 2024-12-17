
# MSV000093184 - Baseline Metabolomics Results ----------------------------

#title: "MSV000093184 - Baseline Metabolomics Results"
#author: "Mauricio Caraballo"
#purpose: clean outputs from fbmn-stats and metadata for downstream analysis using MixOmics package
#input sources: metadata and fbmn_stats_norm
#outputs: filtered tables for downstream analysis

#Install libraries (if not installed)
library(tidyverse)

## NOTE: if using the same script for analysis using mixOmics, make sure to clean the data before loading mixOmics libraries. Some functions are not compatible with Tidyverse.

## Read metadata
metadata <- read_tsv("merged_metadata.tsv")

#print(colnames(metadata))

#Filter metadata and remove ".mzML" suffix from filenames and remove columbs corresponding to blanks and blankQC filenames
# Define rows to remove
rows_to_remove <- c('blank_extraction', 'blank_analysis', 'blank_QC')

metadata_filtered <- metadata |>
  select(filename, starts_with("ATTRIBUTE_"), SampleType, mouse_genotype) |>
  mutate(filename = str_replace(filename, "\\.mzML", "")) |>
  filter(!SampleType %in% rows_to_remove)

## Read fbmn-stats output FILTERED (NO BLANKS AND QCs), imputation and normalized (TIC norm)
fbmn_stats_norm <- read_csv("fbmn_stats_data_clean_export.csv")

## Rename first column
# Rename specific columns using dplyr
fbmn_stats_norm <- fbmn_stats_norm |> rename_with(~ "filename", .cols = 1)

# Select columns "filename" and numerical names
fbmn_stats_norm <- fbmn_stats_norm %>%
  select(filename, matches("^[0-9]"))

# CHECK DATA AND METADATA FILES CONSISTENCY -------------------------------

#Check consistency between data and metadata (output => TRUE if both dataframes are consistent)
fbmn_stats_norm <- arrange(fbmn_stats_norm, filename)
metadata_filtered <- arrange(metadata_filtered, filename)

all(fbmn_stats_norm$filename == metadata_filtered$filename)

## Read data ###
