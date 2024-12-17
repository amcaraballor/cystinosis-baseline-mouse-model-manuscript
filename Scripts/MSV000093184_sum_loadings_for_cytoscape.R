
# MSV000093184 - Baseline Metabolomics Results ----------------------------

---
  #title: MSV000093184 - post mixOmics stats - sum loadings for Cytoscape visualization
  #author: Mauricio Caraballo
  #purpose: sum importance score from PLS-DA to create input file for network in Cytoscape
  #input sources: Loading tables from downstream analysis using mixOmics
  #outputs: Combine all tables and create new column "sum_importance
  ---
  

# SUM IMPORTANCE SCORE AND CREATE CYTOSCAPE OUTPUT --------------
#Install libraries (if not installed)
# Library
library(tidyverse)

# Baseline urine ----------------------------------------------------------

# List of dataframes (replace with dataframes to combine)
Urine_Loadings_list <- list(Loadings_PLSDA_baseline_urine_female_age12, Loadings_PLSDA_baseline_urine_female_age9, Loadings_PLSDA_baseline_urine_female_age6, Loadings_PLSDA_baseline_urine_female_age4, Loadings_PLSDA_baseline_urine_female_age2, Loadings_PLSDA_baseline_urine_male_age12, Loadings_PLSDA_baseline_urine_male_age9, Loadings_PLSDA_baseline_urine_male_age6, Loadings_PLSDA_baseline_urine_male_age4, Loadings_PLSDA_baseline_urine_male_age2)

# Function to add "absol_importance" column and select only "rowname" and "absol_importance"
add_absolute_importance <- function(df) {
  df %>%
    mutate(absol_importance = abs(importance)) %>%
    select(rowname, absol_importance)  # Select only necessary columns
}

# Apply the function to all dataframes in the list
Urine_Loadings_list_abs <- Urine_Loadings_list %>%
  map(add_absolute_importance)

# Combine the dataframes by "rowname" using a pivot to stack "absol_importance"
combined_Urine_Loadings <- Urine_Loadings_list_abs %>%
  reduce(full_join, by = "rowname") %>%
  pivot_longer(cols = starts_with("absol_importance"), values_to = "absol_importance", values_drop_na = TRUE) %>%
  select(rowname, absol_importance)

# Summing the "absol_importance" values by "rowname"
result_Urine_Loadings <- combined_Urine_Loadings %>%
  group_by(rowname) %>%
  summarise(sum_importance = sum(absol_importance, na.rm = TRUE)) %>%
  ungroup()

# View the final result
print(result_Urine_Loadings)

write_csv(result_Urine_Loadings, "result_Urine_sum_Loadings_baseline.csv")

## Adjust Sum Loading table to use as input in Cytoscape
result_Urine_Loadings_cytoscape <- result_Urine_Loadings %>%
  dplyr::mutate(featID = stringr::str_split_fixed(rowname, "_", 2)[, 1])

# Save results

write_csv(result_Urine_Loadings_cytoscape, "result_Urine_sum_Loadings_baseline_cytoscape.csv")


# Baseline feces ----------------------------------------------------------

# List of dataframes (replace with dataframes to combine)
baseline_feces_Loadings_list <- list(Loadings_PLSDA_baseline_feces_female_age12, Loadings_PLSDA_baseline_feces_female_age9, Loadings_PLSDA_baseline_feces_female_age6, Loadings_PLSDA_baseline_feces_female_age4, Loadings_PLSDA_baseline_feces_female_age2, Loadings_PLSDA_baseline_feces_male_age12, Loadings_PLSDA_baseline_feces_male_age9, Loadings_PLSDA_baseline_feces_male_age6, Loadings_PLSDA_baseline_feces_male_age4, Loadings_PLSDA_baseline_feces_male_age2)  # Assuming df1, df2, df3 are your dataframes

# Function to add "absol_importance" column and select only "rowname" and "absol_importance"
add_absolute_importance <- function(df) {
  df %>%
    mutate(absol_importance = abs(importance)) %>%
    select(rowname, absol_importance)  # Select only necessary columns
}

# Apply the function to all dataframes in the list
baseline_feces_Loadings_list_abs <- baseline_feces_Loadings_list %>%
  map(add_absolute_importance)

# Combine the dataframes by "rowname" using a pivot to stack "absol_importance"
combined_baseline_feces_Loadings <- baseline_feces_Loadings_list_abs %>%
  reduce(full_join, by = "rowname") %>%
  pivot_longer(cols = starts_with("absol_importance"), values_to = "absol_importance", values_drop_na = TRUE) %>%
  select(rowname, absol_importance)

# Summing the "absol_importance" values by "rowname"
result_baseline_feces_Loadings <- combined_baseline_feces_Loadings %>%
  group_by(rowname) %>%
  summarise(sum_importance = sum(absol_importance, na.rm = TRUE)) %>%
  ungroup()

# View the final result
print(result_baseline_feces_Loadings)

write_csv(result_baseline_feces_Loadings, "result_baseline_feces_sum_Loadings.csv")

## Adjust Sum Loading table to use as input in Cytoscape
result_baseline_feces_Loadings_cytoscape <- result_baseline_feces_Loadings %>%
  dplyr::mutate(featID = stringr::str_split_fixed(rowname, "_", 2)[, 1])

# Save results

write_csv(result_baseline_feces_Loadings_cytoscape, "result_baseline_feces_sum_Loadings_cytoscape.csv")


# Baseline serum ----------------------------------------------------------

# List of dataframes (replace with dataframes to combine)
baseline_serum_Loadings_list <- list(Loadings_PLSDA_baseline_serum_male_age12, Loadings_PLSDA_baseline_serum_female_age12)  # selected Loading dataframes

# Function to add "absol_importance" column and select only "rowname" and "absol_importance"
add_absolute_importance <- function(df) {
  df %>%
    mutate(absol_importance = abs(importance)) %>%
    select(rowname, absol_importance)  # Select only necessary columns
}

# Apply the function to all dataframes in the list
baseline_serum_Loadings_list_abs <- baseline_serum_Loadings_list %>%
  map(add_absolute_importance)

# Combine the dataframes by "rowname" using a pivot to stack "absol_importance"
combined_baseline_serum_Loadings <- baseline_serum_Loadings_list_abs %>%
  reduce(full_join, by = "rowname") %>%
  pivot_longer(cols = starts_with("absol_importance"), values_to = "absol_importance", values_drop_na = TRUE) %>%
  select(rowname, absol_importance)

# Summing the "absol_importance" values by "rowname"
result_baseline_serum_Loadings <- combined_baseline_serum_Loadings %>%
  group_by(rowname) %>%
  summarise(sum_importance = sum(absol_importance, na.rm = TRUE)) %>%
  ungroup()

# View the final result
print(result_baseline_serum_Loadings)

write_csv(result_baseline_serum_Loadings, "result_baseline_serum_sum_Loadings.csv")

## Adjust Sum Loading table to use as input in Cytoscape
result_baseline_serum_Loadings_cytoscape <- result_baseline_serum_Loadings %>%
  dplyr::mutate(featID = stringr::str_split_fixed(rowname, "_", 2)[, 1])

# Save results

write_csv(result_baseline_serum_Loadings_cytoscape, "result_baseline_serum_sum_Loadings_cytoscape.csv")
