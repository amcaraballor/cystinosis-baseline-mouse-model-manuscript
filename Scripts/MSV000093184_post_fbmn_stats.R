# MSV000093184 - Baseline Metabolomics Results ----------------------------

---
  #title: "MSV000093184 - post FBMN stats"
  #author: "Mauricio Caraballo"
  #purpose: "downstream analysis using MixOmics package"
  #input sources: "metadata and fbmn_stats_norm"
  #outputs: "filtered tables and plots"
  ---
  

# Data cleanup ------------------------------------------------------------

# If data has not been cleaned up before hand, perform the following steps to generate the metadata and fbmn_stats_norm files. 
# Otherwise, proceed to "Downstream analysis with mixOmics"
  
#Install libraries (if not installed)
library(tidyverse)

## Cleaning dataframe before MixOmics
## NOTE: MAKE SURE TO CLEAN THE DATA BEFORE LOADING mixOmics libraries. Some functions are not compatible with Tidyverse.

## Read metadata
metadata <- read_tsv("merged_metadata.tsv")

#print(colnames(metadata))

#FILTERED METADATA and remove 'mzML from filenames and remove blanks and blankQC

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

# Check data and metadata files for consistency -------------------------------

#Check consistency between data and metadata (output => TRUE if both dataframes are consistent)
fbmn_stats_norm <- arrange(fbmn_stats_norm, filename)
metadata_filtered <- arrange(metadata_filtered, filename)

all(fbmn_stats_norm$filename == metadata_filtered$filename)

## Read data ###


# Downstream analysis with mixOmics ---------------------------------------


# LOAD ADDITIONAL LIBRARIES -----------------------------------------------

#Only after data iss cleaned for MixOmics analysis. Avoid incompatibilities with Tidyverse
library(mixOmics)
library(ggpubr)
library(vegan)
library(pheatmap)

### Generate PCA whole and scores ###

PCA_whole <- mixOmics::pca(column_to_rownames(fbmn_stats_norm, var = "filename"), ncomp = 2, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X, metadata_filtered)

# PCA bar plot ------------------------------------------------------------

PCA_whole_plot = mixOmics::pca(column_to_rownames(fbmn_stats_norm, var = "filename"), ncomp = 10, center = TRUE, scale = TRUE)
plot(PCA_whole_plot)

# PCA raw ATTRIBUTE_diet_intervention -------------------------------------

#Filter "not applicable" rows from metadata and data by ATTRIBUTE
fbmn_stats_norm_filt_diet_intervention <- fbmn_stats_norm %>% dplyr::filter(metadata_filtered$ATTRIBUTE_diet_intervention != "not applicable")
metadata_filtered_diet_intervention <- metadata_filtered %>% dplyr::filter(metadata_filtered$ATTRIBUTE_diet_intervention != "not applicable")

# PCA whole and scores
PCA_whole_diet_intervention <- mixOmics::pca(column_to_rownames(fbmn_stats_norm_filt_diet_intervention, var = "filename"), ncomp = 2, scale = TRUE)
PCA_whole_diet_intervention_scores <- data.frame(PCA_whole_diet_intervention$variates$X, metadata_filtered_diet_intervention)

# PCA plot 
PCA_diet_intervention_plot <- PCA_whole_diet_intervention_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "ATTRIBUTE_diet_intervention", alpha = 0.6, 
            title = paste("PCA -", "ATTRIBUTE_diet_intervention", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_diet_intervention$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_diet_intervention$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_diet_intervention_scores %>% group_by(ATTRIBUTE_diet_intervention) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = ATTRIBUTE_diet_intervention), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Generate plot
plot(PCA_diet_intervention_plot)

# PCA raw ATTRIBUTE_mouse_age -------------------------------------

#Filter "not applicable" rows from metadata and data by ATTRIBUTE
fbmn_stats_norm_filt_mouse_age <- fbmn_stats_norm %>% dplyr::filter(metadata_filtered$ATTRIBUTE_mouse_age != "not applicable")
metadata_filtered_mouse_age <- metadata_filtered %>% dplyr::filter(metadata_filtered$ATTRIBUTE_mouse_age != "not applicable")

# PCA whole and scores
PCA_whole_mouse_age <- mixOmics::pca(column_to_rownames(fbmn_stats_norm_filt_mouse_age, var = "filename"), ncomp = 2, scale = TRUE)
PCA_whole_mouse_age_scores <- data.frame(PCA_whole_mouse_age$variates$X, metadata_filtered_mouse_age)

# PCA plot 
PCA_mouse_age_plot <- PCA_whole_mouse_age_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "ATTRIBUTE_mouse_age", alpha = 0.6, 
            title = paste("PCA -", "ATTRIBUTE_mouse_age", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_mouse_age$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_mouse_age$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_mouse_age_scores %>% group_by(ATTRIBUTE_mouse_age) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = ATTRIBUTE_mouse_age), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Generate plot
plot(PCA_mouse_age_plot)


# # PCA raw ATTRIBUTE_mouse_gender -----------------------------------

#Filter "not applicable" rows from metadata and data by ATTRIBUTE
fbmn_stats_norm_filt_mouse_gender <- fbmn_stats_norm %>% dplyr::filter(metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable")
metadata_filtered_mouse_gender <- metadata_filtered %>% dplyr::filter(metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable")

# PCA whole and scores
PCA_whole_mouse_gender <- mixOmics::pca(column_to_rownames(fbmn_stats_norm_filt_mouse_gender, var = "filename"), ncomp = 2, scale = TRUE)
PCA_whole_mouse_gender_scores <- data.frame(PCA_whole_mouse_gender$variates$X, metadata_filtered_mouse_gender)

# PCA plot 
PCA_mouse_gender_plot <- PCA_whole_mouse_gender_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "ATTRIBUTE_mouse_gender", alpha = 0.6, 
            title = paste("PCA -", "ATTRIBUTE_mouse_gender", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_mouse_gender$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_mouse_gender$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_mouse_gender_scores %>% group_by(ATTRIBUTE_mouse_gender) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = ATTRIBUTE_mouse_gender), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Generate plot
plot(PCA_mouse_gender_plot)


# PCA raw ATTRIBUTE_sample_type ------------------------------------------

#Filter "not applicable" rows from metadata and data by ATTRIBUTE
fbmn_stats_norm_filt_sample_type <- fbmn_stats_norm %>% dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable")
metadata_filtered_sample_type <- metadata_filtered %>% dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable")

# PCA whole and scores
PCA_whole_sample_type <- mixOmics::pca(column_to_rownames(fbmn_stats_norm_filt_sample_type, var = "filename"), ncomp = 2, scale = TRUE)
PCA_whole_sample_type_scores <- data.frame(PCA_whole_sample_type$variates$X, metadata_filtered_sample_type)

# PCA plot 
PCA_sample_type_plot <- PCA_whole_sample_type_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "ATTRIBUTE_sample_type", alpha = 0.6, 
            title = paste("PCA -", "ATTRIBUTE_sample_type", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_sample_type$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_sample_type$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_sample_type_scores %>% group_by(ATTRIBUTE_sample_type) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = ATTRIBUTE_sample_type), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Generate plot
plot(PCA_sample_type_plot)


# PCA raw ATTRIBUTE_mouse_treatment ---------------------------------------

#Filter "not applicable" rows from metadata and data by ATTRIBUTE
fbmn_stats_norm_filt_mouse_treatment <- fbmn_stats_norm %>% dplyr::filter(metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable")
metadata_filtered_mouse_treatment <- metadata_filtered %>% dplyr::filter(metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable")

# PCA whole and scores
PCA_whole_mouse_treatment <- mixOmics::pca(column_to_rownames(fbmn_stats_norm_filt_mouse_treatment, var = "filename"), ncomp = 2, scale = TRUE)
PCA_whole_mouse_treatment_scores <- data.frame(PCA_whole_mouse_treatment$variates$X, metadata_filtered_mouse_treatment)

# PCA plot 
PCA_mouse_treatment_plot <- PCA_whole_mouse_treatment_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "ATTRIBUTE_mouse_treatment", alpha = 0.6, 
            title = paste("PCA -", "ATTRIBUTE_mouse_treatment", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_mouse_treatment$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_mouse_treatment$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_mouse_treatment_scores %>% group_by(ATTRIBUTE_mouse_treatment) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = ATTRIBUTE_mouse_treatment), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Generate plot
plot(PCA_mouse_treatment_plot)


# PCA raw by mouse_genotype ----------------------------------------

#Filter "not applicable" rows from metadata and data by ATTRIBUTE
fbmn_stats_norm_filt_mouse_genotype <- fbmn_stats_norm %>% dplyr::filter(metadata_filtered$mouse_genotype != "not applicable")
metadata_filtered_mouse_genotype <- metadata_filtered %>% dplyr::filter(metadata_filtered$mouse_genotype != "not applicable")

# PCA whole and scores
PCA_whole_mouse_genotype <- mixOmics::pca(column_to_rownames(fbmn_stats_norm_filt_mouse_genotype, var = "filename"), ncomp = 2, scale = TRUE)
PCA_whole_mouse_genotype_scores <- data.frame(PCA_whole_mouse_genotype$variates$X, metadata_filtered_mouse_genotype)

# PCA plot 
PCA_mouse_genotype_plot <- PCA_whole_mouse_genotype_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "mouse_genotype", alpha = 0.6, 
            title = paste("PCA -", "mouse_genotype", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_mouse_genotype$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_mouse_genotype$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_mouse_genotype_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = mouse_genotype), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Generate plot
plot(PCA_mouse_genotype_plot)


# PERMANOVA -----------------------------------------------------------

# PERMANOVA raw
#dist_metabolites <- vegdist(column_to_rownames(fbmn_stats_norm, var = "filename"), method = "euclidean")
#permanova <- adonis2(dist_metabolites ~ PCA_whole_scores$ATTRIBUTE_week, PCA_whole_scores, na.action = na.omit)


# BASELINE STUDY - 09/10/2024 ---------------------------------------------


# Ctns-/- BASELINE STUDY => ATTRIBUTE_diet_intervention = ATTRIBUTE_mouse_treatment = "normal diet"


# BASELINE URINE ----------------------------------------------------------


# BASELINE URINE GENDERS (M & F) AND AGE ALL --------------------------------------

## BASELINE, Sample type = "urine", Shape = Genders (M/F), Age2-12 - ALL under ATTRIBUTE_mouse_treatment = normal diet
fbmn_stats_norm_filt_baseline_genotype_urine_age_all_genders_mf <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
                
  )

metadata_filtered_baseline_genotype_urine_age_all_genders_mf <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_genotype_urine_age_all_genders_mf$combined_factor <- factor(
  paste(metadata_filtered_baseline_genotype_urine_age_all_genders_mf$ATTRIBUTE_mouse_gender,
        metadata_filtered_baseline_genotype_urine_age_all_genders_mf$mouse_genotype, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_genotype_urine_age_all_genders_mf <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_genotype_urine_age_all_genders_mf, var = "filename"),
  Y = metadata_filtered_baseline_genotype_urine_age_all_genders_mf$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_genotype_urine_age_all_genders_mf_scores <- data.frame(PLSDA_baseline_genotype_urine_age_all_genders_mf$variates$X, metadata_filtered_baseline_genotype_urine_age_all_genders_mf)

# Plot PLS-DA results
PLSDA_baseline_genotype_urine_age_all_genders_mf_plot <- PLSDA_baseline_genotype_urine_age_all_genders_mf_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    shape = "ATTRIBUTE_mouse_gender",
    alpha = 0.6,
    title = "PLS-DA: baseline study by genotype, urine, ages all and gender M/F",
    xlab = paste("Component 1 (", round(PLSDA_baseline_genotype_urine_age_all_genders_mf$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_genotype_urine_age_all_genders_mf$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  stat_ellipse(aes(comp1, comp2,  linetype = ATTRIBUTE_mouse_gender), level = 0.95, show.legend = TRUE, color = "black") + # Confidence ellipses
  geom_point(data = PLSDA_baseline_genotype_urine_age_all_genders_mf_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 5, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  coord_fixed()

plot(PLSDA_baseline_genotype_urine_age_all_genders_mf_plot)


# BASELINE GENOTYPE URINE AGE GENDER - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_genotype_urine_age_all_genders_mf <- plotLoadings(PLSDA_baseline_genotype_urine_age_all_genders_mf, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_genotype_urine_age_all_genders_mf <- perf(PLSDA_baseline_genotype_urine_age_all_genders_mf, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_genotype_urine_age_all_genders_mf, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf <- as.data.frame(mixOmics::vip(PLSDA_baseline_genotype_urine_age_all_genders_mf))
VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_filter <- VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_filter$ID <- rownames(VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_select <- VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load <- VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_select %>% 
  left_join(Loadings_PLSDA_baseline_genotype_urine_age_all_genders_mf, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load_Top100 <- VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_genotype_urine_age_all_genders_mf_plot <- plotLoadings(PLSDA_baseline_genotype_urine_age_all_genders_mf, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_genotype_urine_age_all_genders_mf") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_genotype_urine_age_all_genders_mf_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load_cytoscape <- VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load, "PLSDA_baseline_genotype_urine_age_all_genders_mf_Load.csv")
write_csv(VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load_cytoscape, "VIPs_PLSDA_baseline_genotype_urine_age_all_genders_mf_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_genotype_urine_age_all_genders_mf_plot, filename = "Loadings_PLSDA_baseline_genotype_urine_age_all_genders_mf.png", device = "png", dpi = "retina")

# END OF BASELINE URINE AGEs ALL AND GENDER M/F - PLS-DA & VIP PLOTS ----------------------



# BASELINE, URINE, MALE, AGE2 - ALL under normal diet --------

fbmn_stats_norm_filt_baseline_urine_male_age2 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_male_age2 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_male_age2$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_male_age2$mouse_genotype,
        metadata_filtered_baseline_urine_male_age2$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_male_age2 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_male_age2, var = "filename"),
  Y = metadata_filtered_baseline_urine_male_age2$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_male_age2_scores <- data.frame(PLSDA_baseline_urine_male_age2$variates$X, metadata_filtered_baseline_urine_male_age2)

# Plot PLS-DA results
PLSDA_baseline_urine_male_age2_plot <- PLSDA_baseline_urine_male_age2_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, male, age 2 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_male_age2$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_male_age2$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_male_age2_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_male_age2_plot)
ggsave(plot = PLSDA_baseline_urine_male_age2_plot, filename = "PLSDA_baseline_urine_male_age2_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE URINE MALE AGE2 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_male_age2 <- plotLoadings(PLSDA_baseline_urine_male_age2, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_male_age2 <- perf(PLSDA_baseline_urine_male_age2, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_male_age2, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_male_age2 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_male_age2))
VIPs_PLSDA_baseline_urine_male_age2_filter <- VIPs_PLSDA_baseline_urine_male_age2 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_male_age2_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_male_age2_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_male_age2_select <- VIPs_PLSDA_baseline_urine_male_age2_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_male_age2_Load <- VIPs_PLSDA_baseline_urine_male_age2_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_male_age2, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_male_age2_Load_Top100 <- VIPs_PLSDA_baseline_urine_male_age2_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_male_age2_plot <- plotLoadings(PLSDA_baseline_urine_male_age2, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_male_age2_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_male_age2") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_male_age2_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_male_age2_Load_cytoscape <- VIPs_PLSDA_baseline_urine_male_age2_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_male_age2_Load, "PLSDA_baseline_urine_male_age2_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_male_age2_Load_cytoscape, "VIPs_PLSDA_baseline_urine_male_age2_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_male_age2_plot, filename = "Loadings_PLSDA_baseline_urine_male_age2.png", device = "png", dpi = "retina")

# END OF BASELINE URINE AGE2 - PLS-DA & VIP PLOTS



# START BASELINE, URINE, MALE, AGE4 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_male_age4 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_male_age4 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_male_age4$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_male_age4$mouse_genotype,
        metadata_filtered_baseline_urine_male_age4$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_male_age4 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_male_age4, var = "filename"),
  Y = metadata_filtered_baseline_urine_male_age4$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_male_age4_scores <- data.frame(PLSDA_baseline_urine_male_age4$variates$X, metadata_filtered_baseline_urine_male_age4)

# Plot PLS-DA results
PLSDA_baseline_urine_male_age4_plot <- PLSDA_baseline_urine_male_age4_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, male, age 4 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_male_age4$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_male_age4$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_male_age4_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_male_age4_plot)
ggsave(plot = PLSDA_baseline_urine_male_age4_plot, filename = "PLSDA_baseline_urine_male_age4_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE URINE MALE AGE4 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_male_age4 <- plotLoadings(PLSDA_baseline_urine_male_age4, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_male_age4 <- perf(PLSDA_baseline_urine_male_age4, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_male_age4, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_male_age4 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_male_age4))
VIPs_PLSDA_baseline_urine_male_age4_filter <- VIPs_PLSDA_baseline_urine_male_age4 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_male_age4_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_male_age4_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_male_age4_select <- VIPs_PLSDA_baseline_urine_male_age4_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_male_age4_Load <- VIPs_PLSDA_baseline_urine_male_age4_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_male_age4, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_male_age4_Load_Top100 <- VIPs_PLSDA_baseline_urine_male_age4_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_male_age4_plot <- plotLoadings(PLSDA_baseline_urine_male_age4, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_male_age4_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_male_age4") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_male_age4_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_male_age4_Load_cytoscape <- VIPs_PLSDA_baseline_urine_male_age4_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_male_age4_Load, "PLSDA_baseline_urine_male_age4_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_male_age4_Load_cytoscape, "VIPs_PLSDA_baseline_urine_male_age4_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_male_age4_plot, filename = "Loadings_PLSDA_baseline_urine_male_age4.png", device = "png", dpi = "retina")

# END OF BASELINE URINE AGE4 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, URINE, MALE, AGE6 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_male_age6 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_male_age6 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_male_age6$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_male_age6$mouse_genotype,
        metadata_filtered_baseline_urine_male_age6$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_male_age6 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_male_age6, var = "filename"),
  Y = metadata_filtered_baseline_urine_male_age6$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_male_age6_scores <- data.frame(PLSDA_baseline_urine_male_age6$variates$X, metadata_filtered_baseline_urine_male_age6)

# Plot PLS-DA results
PLSDA_baseline_urine_male_age6_plot <- PLSDA_baseline_urine_male_age6_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, male, age 6 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_male_age6$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_male_age6$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_male_age6_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_male_age6_plot)
ggsave(plot = PLSDA_baseline_urine_male_age6_plot, filename = "PLSDA_baseline_urine_male_age6_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE URINE MALE AGE6 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_male_age6 <- plotLoadings(PLSDA_baseline_urine_male_age6, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_male_age6 <- perf(PLSDA_baseline_urine_male_age6, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_male_age6, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_male_age6 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_male_age6))
VIPs_PLSDA_baseline_urine_male_age6_filter <- VIPs_PLSDA_baseline_urine_male_age6 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_male_age6_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_male_age6_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_male_age6_select <- VIPs_PLSDA_baseline_urine_male_age6_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_male_age6_Load <- VIPs_PLSDA_baseline_urine_male_age6_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_male_age6, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_male_age6_Load_Top100 <- VIPs_PLSDA_baseline_urine_male_age6_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_male_age6_plot <- plotLoadings(PLSDA_baseline_urine_male_age6, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_male_age6_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_male_age6") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_male_age6_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_male_age6_Load_cytoscape <- VIPs_PLSDA_baseline_urine_male_age6_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_male_age6_Load, "PLSDA_baseline_urine_male_age6_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_male_age6_Load_cytoscape, "VIPs_PLSDA_baseline_urine_male_age6_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_male_age6_plot, filename = "Loadings_PLSDA_baseline_urine_male_age6.png", device = "png", dpi = "retina")

# END OF BASELINE URINE AGE6 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, URINE, MALE, AGE9 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_male_age9 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_male_age9 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_male_age9$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_male_age9$mouse_genotype,
        metadata_filtered_baseline_urine_male_age9$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_male_age9 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_male_age9, var = "filename"),
  Y = metadata_filtered_baseline_urine_male_age9$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_male_age9_scores <- data.frame(PLSDA_baseline_urine_male_age9$variates$X, metadata_filtered_baseline_urine_male_age9)

# Plot PLS-DA results
PLSDA_baseline_urine_male_age9_plot <- PLSDA_baseline_urine_male_age9_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, male, age 9 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_male_age9$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_male_age9$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_male_age9_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_male_age9_plot)
ggsave(plot = PLSDA_baseline_urine_male_age9_plot, filename = "PLSDA_baseline_urine_male_age9_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE URINE MALE AGE9 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_male_age9 <- plotLoadings(PLSDA_baseline_urine_male_age9, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_male_age9 <- perf(PLSDA_baseline_urine_male_age9, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_male_age9, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_male_age9 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_male_age9))
VIPs_PLSDA_baseline_urine_male_age9_filter <- VIPs_PLSDA_baseline_urine_male_age9 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_male_age9_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_male_age9_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_male_age9_select <- VIPs_PLSDA_baseline_urine_male_age9_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_male_age9_Load <- VIPs_PLSDA_baseline_urine_male_age9_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_male_age9, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_male_age9_Load_Top100 <- VIPs_PLSDA_baseline_urine_male_age9_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_male_age9_plot <- plotLoadings(PLSDA_baseline_urine_male_age9, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_male_age9_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_male_age9") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_male_age9_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_male_age9_Load_cytoscape <- VIPs_PLSDA_baseline_urine_male_age9_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_male_age9_Load, "PLSDA_baseline_urine_male_age9_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_male_age9_Load_cytoscape, "VIPs_PLSDA_baseline_urine_male_age9_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_male_age9_plot, filename = "Loadings_PLSDA_baseline_urine_male_age9.png", device = "png", dpi = "retina")

# END OF BASELINE URINE AGE9 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, URINE, MALE, AGE12 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_male_age12 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_male_age12 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_male_age12$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_male_age12$mouse_genotype,
        metadata_filtered_baseline_urine_male_age12$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_male_age12 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_male_age12, var = "filename"),
  Y = metadata_filtered_baseline_urine_male_age12$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_male_age12_scores <- data.frame(PLSDA_baseline_urine_male_age12$variates$X, metadata_filtered_baseline_urine_male_age12)

# Plot PLS-DA results
PLSDA_baseline_urine_male_age12_plot <- PLSDA_baseline_urine_male_age12_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, male, age 12 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_male_age12$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_male_age12$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_male_age12_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_male_age12_plot)
ggsave(plot = PLSDA_baseline_urine_male_age12_plot, filename = "PLSDA_baseline_urine_male_age12_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE URINE MALE AGE12 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_male_age12 <- plotLoadings(PLSDA_baseline_urine_male_age12, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_male_age12 <- perf(PLSDA_baseline_urine_male_age12, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_male_age12, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_male_age12 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_male_age12))
VIPs_PLSDA_baseline_urine_male_age12_filter <- VIPs_PLSDA_baseline_urine_male_age12 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_male_age12_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_male_age12_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_male_age12_select <- VIPs_PLSDA_baseline_urine_male_age12_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_male_age12_Load <- VIPs_PLSDA_baseline_urine_male_age12_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_male_age12, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_male_age12_Load_Top100 <- VIPs_PLSDA_baseline_urine_male_age12_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_male_age12_plot <- plotLoadings(PLSDA_baseline_urine_male_age12, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_male_age12_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_male_age12") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_male_age12_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_male_age12_Load_cytoscape <- VIPs_PLSDA_baseline_urine_male_age12_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_male_age12_Load, "PLSDA_baseline_urine_male_age12_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_male_age12_Load_cytoscape, "VIPs_PLSDA_baseline_urine_male_age12_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_male_age12_plot, filename = "Loadings_PLSDA_baseline_urine_male_age12.png", device = "png", dpi = "retina")

# END OF BASELINE URINE MALE AGE12 - PLS-DA & VIP PLOTS ----------------------


# URINE FEMALE ------------------------------------------------------------

# BASELINE, URINE, FEMALE, AGE2 - ALL under ATTRIBUTE_mouse_treatment --------

fbmn_stats_norm_filt_baseline_urine_female_age2 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_female_age2 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_female_age2$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_female_age2$mouse_genotype,
        metadata_filtered_baseline_urine_female_age2$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_female_age2 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_female_age2, var = "filename"),
  Y = metadata_filtered_baseline_urine_female_age2$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_female_age2_scores <- data.frame(PLSDA_baseline_urine_female_age2$variates$X, metadata_filtered_baseline_urine_female_age2)

# Plot PLS-DA results
PLSDA_baseline_urine_female_age2_plot <- PLSDA_baseline_urine_female_age2_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, female, age 2 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_female_age2$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_female_age2$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_female_age2_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_female_age2_plot)
ggsave(plot = PLSDA_baseline_urine_female_age2_plot, filename = "PLSDA_baseline_urine_female_age2_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE URINE FEMALE AGE2 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_female_age2 <- plotLoadings(PLSDA_baseline_urine_female_age2, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_female_age2 <- perf(PLSDA_baseline_urine_female_age2, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_female_age2, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_female_age2 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_female_age2))
VIPs_PLSDA_baseline_urine_female_age2_filter <- VIPs_PLSDA_baseline_urine_female_age2 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_female_age2_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_female_age2_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_female_age2_select <- VIPs_PLSDA_baseline_urine_female_age2_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_female_age2_Load <- VIPs_PLSDA_baseline_urine_female_age2_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_female_age2, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_female_age2_Load_Top100 <- VIPs_PLSDA_baseline_urine_female_age2_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_female_age2_plot <- plotLoadings(PLSDA_baseline_urine_female_age2, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_female_age2_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_female_age2") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_female_age2_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_female_age2_Load_cytoscape <- VIPs_PLSDA_baseline_urine_female_age2_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_female_age2_Load, "PLSDA_baseline_urine_female_age2_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_female_age2_Load_cytoscape, "VIPs_PLSDA_baseline_urine_female_age2_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_female_age2_plot, filename = "Loadings_PLSDA_baseline_urine_female_age2.png", device = "png", dpi = "retina")

# END OF BASELINE URINE FEMALE AGE2 - PLS-DA & VIP PLOTS



# START BASELINE, URINE, FEMALE, AGE4 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_female_age4 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_female_age4 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_female_age4$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_female_age4$mouse_genotype,
        metadata_filtered_baseline_urine_female_age4$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_female_age4 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_female_age4, var = "filename"),
  Y = metadata_filtered_baseline_urine_female_age4$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_female_age4_scores <- data.frame(PLSDA_baseline_urine_female_age4$variates$X, metadata_filtered_baseline_urine_female_age4)

# Plot PLS-DA results
PLSDA_baseline_urine_female_age4_plot <- PLSDA_baseline_urine_female_age4_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, female, age 4 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_female_age4$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_female_age4$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_female_age4_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_female_age4_plot)
ggsave(plot = PLSDA_baseline_urine_female_age4_plot, filename = "PLSDA_baseline_urine_female_age4_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE URINE MALE AGE4 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_female_age4 <- plotLoadings(PLSDA_baseline_urine_female_age4, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_female_age4 <- perf(PLSDA_baseline_urine_female_age4, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_female_age4, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_female_age4 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_female_age4))
VIPs_PLSDA_baseline_urine_female_age4_filter <- VIPs_PLSDA_baseline_urine_female_age4 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_female_age4_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_female_age4_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_female_age4_select <- VIPs_PLSDA_baseline_urine_female_age4_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_female_age4_Load <- VIPs_PLSDA_baseline_urine_female_age4_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_female_age4, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_female_age4_Load_Top100 <- VIPs_PLSDA_baseline_urine_female_age4_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_female_age4_plot <- plotLoadings(PLSDA_baseline_urine_female_age4, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_female_age4_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_female_age4") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_female_age4_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_female_age4_Load_cytoscape <- VIPs_PLSDA_baseline_urine_female_age4_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_female_age4_Load, "PLSDA_baseline_urine_female_age4_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_female_age4_Load_cytoscape, "VIPs_PLSDA_baseline_urine_female_age4_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_female_age4_plot, filename = "Loadings_PLSDA_baseline_urine_female_age4.png", device = "png", dpi = "retina")

# END OF BASELINE URINE FEMALE AGE4 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, URINE, FEMALE, AGE6 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_female_age6 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_female_age6 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_female_age6$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_female_age6$mouse_genotype,
        metadata_filtered_baseline_urine_female_age6$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_female_age6 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_female_age6, var = "filename"),
  Y = metadata_filtered_baseline_urine_female_age6$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_female_age6_scores <- data.frame(PLSDA_baseline_urine_female_age6$variates$X, metadata_filtered_baseline_urine_female_age6)

# Plot PLS-DA results
PLSDA_baseline_urine_female_age6_plot <- PLSDA_baseline_urine_female_age6_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, female, age 6 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_female_age6$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_female_age6$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_female_age6_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_female_age6_plot)
ggsave(plot = PLSDA_baseline_urine_female_age6_plot, filename = "PLSDA_baseline_urine_female_age6_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE URINE FEMALE AGE6 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_female_age6 <- plotLoadings(PLSDA_baseline_urine_female_age6, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_female_age6 <- perf(PLSDA_baseline_urine_female_age6, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_female_age6, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_female_age6 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_female_age6))
VIPs_PLSDA_baseline_urine_female_age6_filter <- VIPs_PLSDA_baseline_urine_female_age6 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_female_age6_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_female_age6_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_female_age6_select <- VIPs_PLSDA_baseline_urine_female_age6_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_female_age6_Load <- VIPs_PLSDA_baseline_urine_female_age6_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_female_age6, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_female_age6_Load_Top100 <- VIPs_PLSDA_baseline_urine_female_age6_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_female_age6_plot <- plotLoadings(PLSDA_baseline_urine_female_age6, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_female_age6_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_female_age6") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_female_age6_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_female_age6_Load_cytoscape <- VIPs_PLSDA_baseline_urine_female_age6_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_female_age6_Load, "PLSDA_baseline_urine_female_age6_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_female_age6_Load_cytoscape, "VIPs_PLSDA_baseline_urine_female_age6_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_female_age6_plot, filename = "Loadings_PLSDA_baseline_urine_female_age6.png", device = "png", dpi = "retina")

# END OF BASELINE URINE FEMALE AGE6 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, URINE, FEMALE, AGE9 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_female_age9 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_female_age9 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_female_age9$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_female_age9$mouse_genotype,
        metadata_filtered_baseline_urine_female_age9$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_female_age9 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_female_age9, var = "filename"),
  Y = metadata_filtered_baseline_urine_female_age9$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_female_age9_scores <- data.frame(PLSDA_baseline_urine_female_age9$variates$X, metadata_filtered_baseline_urine_female_age9)

# Plot PLS-DA results
PLSDA_baseline_urine_female_age9_plot <- PLSDA_baseline_urine_female_age9_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, female, age 9 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_female_age9$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_female_age9$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_female_age9_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_female_age9_plot)
ggsave(plot = PLSDA_baseline_urine_female_age9_plot, filename = "PLSDA_baseline_urine_female_age9_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE URINE FEMALE AGE9 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_female_age9 <- plotLoadings(PLSDA_baseline_urine_female_age9, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_female_age9 <- perf(PLSDA_baseline_urine_female_age9, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_female_age9, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_female_age9 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_female_age9))
VIPs_PLSDA_baseline_urine_female_age9_filter <- VIPs_PLSDA_baseline_urine_female_age9 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_female_age9_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_female_age9_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_female_age9_select <- VIPs_PLSDA_baseline_urine_female_age9_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_female_age9_Load <- VIPs_PLSDA_baseline_urine_female_age9_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_female_age9, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_female_age9_Load_Top100 <- VIPs_PLSDA_baseline_urine_female_age9_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_female_age9_plot <- plotLoadings(PLSDA_baseline_urine_female_age9, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_female_age9_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_female_age9") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_female_age9_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_female_age9_Load_cytoscape <- VIPs_PLSDA_baseline_urine_female_age9_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_female_age9_Load, "PLSDA_baseline_urine_female_age9_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_female_age9_Load_cytoscape, "VIPs_PLSDA_baseline_urine_female_age9_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_female_age9_plot, filename = "Loadings_PLSDA_baseline_urine_female_age9.png", device = "png", dpi = "retina")

# END OF BASELINE URINE AGE9 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, URINE, FEMALE, AGE12 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_urine_female_age12 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_urine_female_age12 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_urine_female_age12$combined_factor <- factor(
  paste(metadata_filtered_baseline_urine_female_age12$mouse_genotype,
        metadata_filtered_baseline_urine_female_age12$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_urine_female_age12 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_urine_female_age12, var = "filename"),
  Y = metadata_filtered_baseline_urine_female_age12$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_urine_female_age12_scores <- data.frame(PLSDA_baseline_urine_female_age12$variates$X, metadata_filtered_baseline_urine_female_age12)

# Plot PLS-DA results
PLSDA_baseline_urine_female_age12_plot <- PLSDA_baseline_urine_female_age12_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (urine, female, age 12 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_urine_female_age12$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_urine_female_age12$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_urine_female_age12_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_urine_female_age12_plot)
ggsave(plot = PLSDA_baseline_urine_female_age12_plot, filename = "PLSDA_baseline_urine_female_age12_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE URINE FEMALE AGE12 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_urine_female_age12 <- plotLoadings(PLSDA_baseline_urine_female_age12, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_urine_female_age12 <- perf(PLSDA_baseline_urine_female_age12, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_urine_female_age12, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_urine_female_age12 <- as.data.frame(mixOmics::vip(PLSDA_baseline_urine_female_age12))
VIPs_PLSDA_baseline_urine_female_age12_filter <- VIPs_PLSDA_baseline_urine_female_age12 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_urine_female_age12_filter$ID <- rownames(VIPs_PLSDA_baseline_urine_female_age12_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_urine_female_age12_select <- VIPs_PLSDA_baseline_urine_female_age12_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_urine_female_age12_Load <- VIPs_PLSDA_baseline_urine_female_age12_select %>% 
  left_join(Loadings_PLSDA_baseline_urine_female_age12, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_urine_female_age12_Load_Top100 <- VIPs_PLSDA_baseline_urine_female_age12_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_urine_female_age12_plot <- plotLoadings(PLSDA_baseline_urine_female_age12, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_urine_female_age12_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_urine_female_age12") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_urine_female_age12_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_urine_female_age12_Load_cytoscape <- VIPs_PLSDA_baseline_urine_female_age12_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_urine_female_age12_Load, "PLSDA_baseline_urine_female_age12_Load.csv")
write_csv(VIPs_PLSDA_baseline_urine_female_age12_Load_cytoscape, "VIPs_PLSDA_baseline_urine_female_age12_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_urine_female_age12_plot, filename = "Loadings_PLSDA_baseline_urine_female_age12.png", device = "png", dpi = "retina")

# END OF BASELINE URINE FEMALE AGE12 - PLS-DA & VIP PLOTS ----------------------


# BASELINE FECES ALL (GENOTYPE, GENDERS & AGE) ----------------------------

# BASELINE, FECES ALL (GENOTYPE, GENDERS & AGE) ALL under ATTRIBUTE_mouse_treatment --------

fbmn_stats_norm_filt_baseline_feces_all <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_all <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_all$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_all$mouse_genotype,
        metadata_filtered_baseline_feces_all$ATTRIBUTE_mouse_gender,
        metadata_filtered_baseline_feces_all$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_all <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_all, var = "filename"),
  Y = metadata_filtered_baseline_feces_all$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_all_scores <- data.frame(PLSDA_baseline_feces_all$variates$X, metadata_filtered_baseline_feces_all)

# Plot PLS-DA results
PLSDA_baseline_feces_all_plot <- PLSDA_baseline_feces_all_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, genders (M/F), ages (2-12 months))",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_all$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_all$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = TRUE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_all_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_all_plot)
ggsave(plot = PLSDA_baseline_feces_all_plot, filename = "PLSDA_baseline_feces_all_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE feces ALL - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_all <- plotLoadings(PLSDA_baseline_feces_all, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_all <- perf(PLSDA_baseline_feces_all, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_all, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_all <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_all))
VIPs_PLSDA_baseline_feces_all_filter <- VIPs_PLSDA_baseline_feces_all %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_all_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_all_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_all_select <- VIPs_PLSDA_baseline_feces_all_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_all_Load <- VIPs_PLSDA_baseline_feces_all_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_all, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_all_Load_Top100 <- VIPs_PLSDA_baseline_feces_all_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_all_plot <- plotLoadings(PLSDA_baseline_feces_all, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_all_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_all") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_all_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_all_Load_cytoscape <- VIPs_PLSDA_baseline_feces_all_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_all_Load, "PLSDA_baseline_feces_all_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_all_Load_cytoscape, "VIPs_PLSDA_baseline_feces_all_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_all_plot, filename = "Loadings_PLSDA_baseline_feces_all.png", device = "png", dpi = "retina")

# END OF BASELINE feces ALL - PLS-DA & VIP PLOTS



# BASELINE FECES ----------------------------------------------------------


# BASELINE, FECES, MALE, AGE2 - ALL under ATTRIBUTE_mouse_treatment --------

fbmn_stats_norm_filt_baseline_feces_male_age2 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_male_age2 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_male_age2$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_male_age2$mouse_genotype,
        metadata_filtered_baseline_feces_male_age2$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_male_age2 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_male_age2, var = "filename"),
  Y = metadata_filtered_baseline_feces_male_age2$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_male_age2_scores <- data.frame(PLSDA_baseline_feces_male_age2$variates$X, metadata_filtered_baseline_feces_male_age2)

# Plot PLS-DA results
PLSDA_baseline_feces_male_age2_plot <- PLSDA_baseline_feces_male_age2_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, male, age 2 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_male_age2$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_male_age2$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_male_age2_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_male_age2_plot)
ggsave(plot = PLSDA_baseline_feces_male_age2_plot, filename = "PLSDA_baseline_feces_male_age2_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE feces MALE AGE2 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_male_age2 <- plotLoadings(PLSDA_baseline_feces_male_age2, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_male_age2 <- perf(PLSDA_baseline_feces_male_age2, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_male_age2, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_male_age2 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_male_age2))
VIPs_PLSDA_baseline_feces_male_age2_filter <- VIPs_PLSDA_baseline_feces_male_age2 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_male_age2_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_male_age2_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_male_age2_select <- VIPs_PLSDA_baseline_feces_male_age2_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_male_age2_Load <- VIPs_PLSDA_baseline_feces_male_age2_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_male_age2, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_male_age2_Load_Top100 <- VIPs_PLSDA_baseline_feces_male_age2_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_male_age2_plot <- plotLoadings(PLSDA_baseline_feces_male_age2, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_male_age2_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_male_age2") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_male_age2_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_male_age2_Load_cytoscape <- VIPs_PLSDA_baseline_feces_male_age2_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_male_age2_Load, "PLSDA_baseline_feces_male_age2_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_male_age2_Load_cytoscape, "VIPs_PLSDA_baseline_feces_male_age2_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_male_age2_plot, filename = "Loadings_PLSDA_baseline_feces_male_age2.png", device = "png", dpi = "retina")

# END OF BASELINE feces AGE2 - PLS-DA & VIP PLOTS



# START BASELINE, feces, MALE, AGE4 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_male_age4 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_male_age4 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_male_age4$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_male_age4$mouse_genotype,
        metadata_filtered_baseline_feces_male_age4$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_male_age4 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_male_age4, var = "filename"),
  Y = metadata_filtered_baseline_feces_male_age4$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_male_age4_scores <- data.frame(PLSDA_baseline_feces_male_age4$variates$X, metadata_filtered_baseline_feces_male_age4)

# Plot PLS-DA results
PLSDA_baseline_feces_male_age4_plot <- PLSDA_baseline_feces_male_age4_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, male, age 4 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_male_age4$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_male_age4$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_male_age4_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_male_age4_plot)
ggsave(plot = PLSDA_baseline_feces_male_age4_plot, filename = "PLSDA_baseline_feces_male_age4_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE feces MALE AGE4 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_male_age4 <- plotLoadings(PLSDA_baseline_feces_male_age4, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_male_age4 <- perf(PLSDA_baseline_feces_male_age4, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_male_age4, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_male_age4 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_male_age4))
VIPs_PLSDA_baseline_feces_male_age4_filter <- VIPs_PLSDA_baseline_feces_male_age4 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_male_age4_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_male_age4_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_male_age4_select <- VIPs_PLSDA_baseline_feces_male_age4_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_male_age4_Load <- VIPs_PLSDA_baseline_feces_male_age4_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_male_age4, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_male_age4_Load_Top100 <- VIPs_PLSDA_baseline_feces_male_age4_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_male_age4_plot <- plotLoadings(PLSDA_baseline_feces_male_age4, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_male_age4_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_male_age4") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_male_age4_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_male_age4_Load_cytoscape <- VIPs_PLSDA_baseline_feces_male_age4_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_male_age4_Load, "PLSDA_baseline_feces_male_age4_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_male_age4_Load_cytoscape, "VIPs_PLSDA_baseline_feces_male_age4_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_male_age4_plot, filename = "Loadings_PLSDA_baseline_feces_male_age4.png", device = "png", dpi = "retina")

# END OF BASELINE feces AGE4 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, feces, MALE, AGE6 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_male_age6 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_male_age6 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_male_age6$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_male_age6$mouse_genotype,
        metadata_filtered_baseline_feces_male_age6$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_male_age6 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_male_age6, var = "filename"),
  Y = metadata_filtered_baseline_feces_male_age6$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_male_age6_scores <- data.frame(PLSDA_baseline_feces_male_age6$variates$X, metadata_filtered_baseline_feces_male_age6)

# Plot PLS-DA results
PLSDA_baseline_feces_male_age6_plot <- PLSDA_baseline_feces_male_age6_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, male, age 6 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_male_age6$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_male_age6$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_male_age6_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_male_age6_plot)
ggsave(plot = PLSDA_baseline_feces_male_age6_plot, filename = "PLSDA_baseline_feces_male_age6_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE feces MALE AGE6 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_male_age6 <- plotLoadings(PLSDA_baseline_feces_male_age6, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_male_age6 <- perf(PLSDA_baseline_feces_male_age6, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_male_age6, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_male_age6 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_male_age6))
VIPs_PLSDA_baseline_feces_male_age6_filter <- VIPs_PLSDA_baseline_feces_male_age6 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_male_age6_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_male_age6_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_male_age6_select <- VIPs_PLSDA_baseline_feces_male_age6_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_male_age6_Load <- VIPs_PLSDA_baseline_feces_male_age6_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_male_age6, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_male_age6_Load_Top100 <- VIPs_PLSDA_baseline_feces_male_age6_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_male_age6_plot <- plotLoadings(PLSDA_baseline_feces_male_age6, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_male_age6_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_male_age6") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_male_age6_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_male_age6_Load_cytoscape <- VIPs_PLSDA_baseline_feces_male_age6_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_male_age6_Load, "PLSDA_baseline_feces_male_age6_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_male_age6_Load_cytoscape, "VIPs_PLSDA_baseline_feces_male_age6_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_male_age6_plot, filename = "Loadings_PLSDA_baseline_feces_male_age6.png", device = "png", dpi = "retina")

# END OF BASELINE feces AGE6 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, feces, MALE, AGE9 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_male_age9 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_male_age9 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_male_age9$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_male_age9$mouse_genotype,
        metadata_filtered_baseline_feces_male_age9$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_male_age9 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_male_age9, var = "filename"),
  Y = metadata_filtered_baseline_feces_male_age9$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_male_age9_scores <- data.frame(PLSDA_baseline_feces_male_age9$variates$X, metadata_filtered_baseline_feces_male_age9)

# Plot PLS-DA results
PLSDA_baseline_feces_male_age9_plot <- PLSDA_baseline_feces_male_age9_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, male, age 9 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_male_age9$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_male_age9$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_male_age9_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_male_age9_plot)
ggsave(plot = PLSDA_baseline_feces_male_age9_plot, filename = "PLSDA_baseline_feces_male_age9_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE feces MALE AGE9 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_male_age9 <- plotLoadings(PLSDA_baseline_feces_male_age9, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_male_age9 <- perf(PLSDA_baseline_feces_male_age9, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_male_age9, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_male_age9 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_male_age9))
VIPs_PLSDA_baseline_feces_male_age9_filter <- VIPs_PLSDA_baseline_feces_male_age9 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_male_age9_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_male_age9_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_male_age9_select <- VIPs_PLSDA_baseline_feces_male_age9_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_male_age9_Load <- VIPs_PLSDA_baseline_feces_male_age9_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_male_age9, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_male_age9_Load_Top100 <- VIPs_PLSDA_baseline_feces_male_age9_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_male_age9_plot <- plotLoadings(PLSDA_baseline_feces_male_age9, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_male_age9_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_male_age9") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_male_age9_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_male_age9_Load_cytoscape <- VIPs_PLSDA_baseline_feces_male_age9_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_male_age9_Load, "PLSDA_baseline_feces_male_age9_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_male_age9_Load_cytoscape, "VIPs_PLSDA_baseline_feces_male_age9_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_male_age9_plot, filename = "Loadings_PLSDA_baseline_feces_male_age9.png", device = "png", dpi = "retina")

# END OF BASELINE feces AGE9 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, feces, MALE, AGE12 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_male_age12 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_male_age12 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_male_age12$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_male_age12$mouse_genotype,
        metadata_filtered_baseline_feces_male_age12$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_male_age12 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_male_age12, var = "filename"),
  Y = metadata_filtered_baseline_feces_male_age12$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_male_age12_scores <- data.frame(PLSDA_baseline_feces_male_age12$variates$X, metadata_filtered_baseline_feces_male_age12)

# Plot PLS-DA results
PLSDA_baseline_feces_male_age12_plot <- PLSDA_baseline_feces_male_age12_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, male, age 12 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_male_age12$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_male_age12$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_male_age12_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_male_age12_plot)
ggsave(plot = PLSDA_baseline_feces_male_age12_plot, filename = "PLSDA_baseline_feces_male_age12_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE feces MALE AGE12 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_male_age12 <- plotLoadings(PLSDA_baseline_feces_male_age12, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_male_age12 <- perf(PLSDA_baseline_feces_male_age12, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_male_age12, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_male_age12 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_male_age12))
VIPs_PLSDA_baseline_feces_male_age12_filter <- VIPs_PLSDA_baseline_feces_male_age12 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_male_age12_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_male_age12_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_male_age12_select <- VIPs_PLSDA_baseline_feces_male_age12_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_male_age12_Load <- VIPs_PLSDA_baseline_feces_male_age12_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_male_age12, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_male_age12_Load_Top100 <- VIPs_PLSDA_baseline_feces_male_age12_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_male_age12_plot <- plotLoadings(PLSDA_baseline_feces_male_age12, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_male_age12_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_male_age12") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_male_age12_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_male_age12_Load_cytoscape <- VIPs_PLSDA_baseline_feces_male_age12_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_male_age12_Load, "PLSDA_baseline_feces_male_age12_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_male_age12_Load_cytoscape, "VIPs_PLSDA_baseline_feces_male_age12_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_male_age12_plot, filename = "Loadings_PLSDA_baseline_feces_male_age12.png", device = "png", dpi = "retina")

# END OF BASELINE feces MALE AGE12 - PLS-DA & VIP PLOTS ----------------------


# feces FEMALE ------------------------------------------------------------

# BASELINE, feces, FEMALE, AGE2 - ALL under ATTRIBUTE_mouse_treatment --------

fbmn_stats_norm_filt_baseline_feces_female_age2 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_female_age2 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_female_age2$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_female_age2$mouse_genotype,
        metadata_filtered_baseline_feces_female_age2$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_female_age2 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_female_age2, var = "filename"),
  Y = metadata_filtered_baseline_feces_female_age2$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_female_age2_scores <- data.frame(PLSDA_baseline_feces_female_age2$variates$X, metadata_filtered_baseline_feces_female_age2)

# Plot PLS-DA results
PLSDA_baseline_feces_female_age2_plot <- PLSDA_baseline_feces_female_age2_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, female, age 2 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_female_age2$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_female_age2$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_female_age2_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_female_age2_plot)
ggsave(plot = PLSDA_baseline_feces_female_age2_plot, filename = "PLSDA_baseline_feces_female_age2_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE feces FEMALE AGE2 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_female_age2 <- plotLoadings(PLSDA_baseline_feces_female_age2, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_female_age2 <- perf(PLSDA_baseline_feces_female_age2, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_female_age2, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_female_age2 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_female_age2))
VIPs_PLSDA_baseline_feces_female_age2_filter <- VIPs_PLSDA_baseline_feces_female_age2 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_female_age2_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_female_age2_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_female_age2_select <- VIPs_PLSDA_baseline_feces_female_age2_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_female_age2_Load <- VIPs_PLSDA_baseline_feces_female_age2_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_female_age2, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_female_age2_Load_Top100 <- VIPs_PLSDA_baseline_feces_female_age2_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_female_age2_plot <- plotLoadings(PLSDA_baseline_feces_female_age2, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_female_age2_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_female_age2") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_female_age2_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_female_age2_Load_cytoscape <- VIPs_PLSDA_baseline_feces_female_age2_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_female_age2_Load, "PLSDA_baseline_feces_female_age2_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_female_age2_Load_cytoscape, "VIPs_PLSDA_baseline_feces_female_age2_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_female_age2_plot, filename = "Loadings_PLSDA_baseline_feces_female_age2.png", device = "png", dpi = "retina")

# END OF BASELINE feces FEMALE AGE2 - PLS-DA & VIP PLOTS



# START BASELINE, feces, FEMALE, AGE4 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_female_age4 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_female_age4 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_female_age4$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_female_age4$mouse_genotype,
        metadata_filtered_baseline_feces_female_age4$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_female_age4 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_female_age4, var = "filename"),
  Y = metadata_filtered_baseline_feces_female_age4$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_female_age4_scores <- data.frame(PLSDA_baseline_feces_female_age4$variates$X, metadata_filtered_baseline_feces_female_age4)

# Plot PLS-DA results
PLSDA_baseline_feces_female_age4_plot <- PLSDA_baseline_feces_female_age4_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, female, age 4 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_female_age4$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_female_age4$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_female_age4_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_female_age4_plot)
ggsave(plot = PLSDA_baseline_feces_female_age4_plot, filename = "PLSDA_baseline_feces_female_age4_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE feces MALE AGE4 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_female_age4 <- plotLoadings(PLSDA_baseline_feces_female_age4, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_female_age4 <- perf(PLSDA_baseline_feces_female_age4, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_female_age4, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_female_age4 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_female_age4))
VIPs_PLSDA_baseline_feces_female_age4_filter <- VIPs_PLSDA_baseline_feces_female_age4 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_female_age4_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_female_age4_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_female_age4_select <- VIPs_PLSDA_baseline_feces_female_age4_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_female_age4_Load <- VIPs_PLSDA_baseline_feces_female_age4_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_female_age4, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_female_age4_Load_Top100 <- VIPs_PLSDA_baseline_feces_female_age4_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_female_age4_plot <- plotLoadings(PLSDA_baseline_feces_female_age4, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_female_age4_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_female_age4") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_female_age4_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_female_age4_Load_cytoscape <- VIPs_PLSDA_baseline_feces_female_age4_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_female_age4_Load, "PLSDA_baseline_feces_female_age4_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_female_age4_Load_cytoscape, "VIPs_PLSDA_baseline_feces_female_age4_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_female_age4_plot, filename = "Loadings_PLSDA_baseline_feces_female_age4.png", device = "png", dpi = "retina")

# END OF BASELINE feces FEMALE AGE4 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, feces, FEMALE, AGE6 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_female_age6 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_female_age6 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_female_age6$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_female_age6$mouse_genotype,
        metadata_filtered_baseline_feces_female_age6$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_female_age6 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_female_age6, var = "filename"),
  Y = metadata_filtered_baseline_feces_female_age6$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_female_age6_scores <- data.frame(PLSDA_baseline_feces_female_age6$variates$X, metadata_filtered_baseline_feces_female_age6)

# Plot PLS-DA results
PLSDA_baseline_feces_female_age6_plot <- PLSDA_baseline_feces_female_age6_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, female, age 6 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_female_age6$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_female_age6$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_female_age6_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_female_age6_plot)
ggsave(plot = PLSDA_baseline_feces_female_age6_plot, filename = "PLSDA_baseline_feces_female_age6_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE feces FEMALE AGE6 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_female_age6 <- plotLoadings(PLSDA_baseline_feces_female_age6, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_female_age6 <- perf(PLSDA_baseline_feces_female_age6, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_female_age6, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_female_age6 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_female_age6))
VIPs_PLSDA_baseline_feces_female_age6_filter <- VIPs_PLSDA_baseline_feces_female_age6 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_female_age6_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_female_age6_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_female_age6_select <- VIPs_PLSDA_baseline_feces_female_age6_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_female_age6_Load <- VIPs_PLSDA_baseline_feces_female_age6_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_female_age6, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_female_age6_Load_Top100 <- VIPs_PLSDA_baseline_feces_female_age6_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_female_age6_plot <- plotLoadings(PLSDA_baseline_feces_female_age6, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_female_age6_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_female_age6") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_female_age6_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_female_age6_Load_cytoscape <- VIPs_PLSDA_baseline_feces_female_age6_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_female_age6_Load, "PLSDA_baseline_feces_female_age6_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_female_age6_Load_cytoscape, "VIPs_PLSDA_baseline_feces_female_age6_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_female_age6_plot, filename = "Loadings_PLSDA_baseline_feces_female_age6.png", device = "png", dpi = "retina")

# END OF BASELINE feces FEMALE AGE6 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, feces, FEMALE, AGE9 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_female_age9 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_female_age9 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_female_age9$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_female_age9$mouse_genotype,
        metadata_filtered_baseline_feces_female_age9$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_female_age9 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_female_age9, var = "filename"),
  Y = metadata_filtered_baseline_feces_female_age9$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_female_age9_scores <- data.frame(PLSDA_baseline_feces_female_age9$variates$X, metadata_filtered_baseline_feces_female_age9)

# Plot PLS-DA results
PLSDA_baseline_feces_female_age9_plot <- PLSDA_baseline_feces_female_age9_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, female, age 9 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_female_age9$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_female_age9$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_female_age9_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_female_age9_plot)
ggsave(plot = PLSDA_baseline_feces_female_age9_plot, filename = "PLSDA_baseline_feces_female_age9_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE feces FEMALE AGE9 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_female_age9 <- plotLoadings(PLSDA_baseline_feces_female_age9, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_female_age9 <- perf(PLSDA_baseline_feces_female_age9, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_female_age9, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_female_age9 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_female_age9))
VIPs_PLSDA_baseline_feces_female_age9_filter <- VIPs_PLSDA_baseline_feces_female_age9 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_female_age9_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_female_age9_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_female_age9_select <- VIPs_PLSDA_baseline_feces_female_age9_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_female_age9_Load <- VIPs_PLSDA_baseline_feces_female_age9_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_female_age9, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_female_age9_Load_Top100 <- VIPs_PLSDA_baseline_feces_female_age9_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_female_age9_plot <- plotLoadings(PLSDA_baseline_feces_female_age9, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_female_age9_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_female_age9") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_female_age9_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_female_age9_Load_cytoscape <- VIPs_PLSDA_baseline_feces_female_age9_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_female_age9_Load, "PLSDA_baseline_feces_female_age9_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_female_age9_Load_cytoscape, "VIPs_PLSDA_baseline_feces_female_age9_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_female_age9_plot, filename = "Loadings_PLSDA_baseline_feces_female_age9.png", device = "png", dpi = "retina")

# END OF BASELINE feces AGE9 - PLS-DA & VIP PLOTS ----------------------


# START BASELINE, feces, FEMALE, AGE12 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" --------

fbmn_stats_norm_filt_baseline_feces_female_age12 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_feces_female_age12 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "serum" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_feces_female_age12$combined_factor <- factor(
  paste(metadata_filtered_baseline_feces_female_age12$mouse_genotype,
        metadata_filtered_baseline_feces_female_age12$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_feces_female_age12 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_feces_female_age12, var = "filename"),
  Y = metadata_filtered_baseline_feces_female_age12$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_feces_female_age12_scores <- data.frame(PLSDA_baseline_feces_female_age12$variates$X, metadata_filtered_baseline_feces_female_age12)

# Plot PLS-DA results
PLSDA_baseline_feces_female_age12_plot <- PLSDA_baseline_feces_female_age12_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (feces, female, age 12 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_feces_female_age12$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_feces_female_age12$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_feces_female_age12_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_feces_female_age12_plot)
ggsave(plot = PLSDA_baseline_feces_female_age12_plot, filename = "PLSDA_baseline_feces_female_age12_plot.png", device = "png", dpi = "retina")


# BASELINE GENOTYPE feces FEMALE AGE12 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_feces_female_age12 <- plotLoadings(PLSDA_baseline_feces_female_age12, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_feces_female_age12 <- perf(PLSDA_baseline_feces_female_age12, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_feces_female_age12, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_feces_female_age12 <- as.data.frame(mixOmics::vip(PLSDA_baseline_feces_female_age12))
VIPs_PLSDA_baseline_feces_female_age12_filter <- VIPs_PLSDA_baseline_feces_female_age12 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_feces_female_age12_filter$ID <- rownames(VIPs_PLSDA_baseline_feces_female_age12_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_feces_female_age12_select <- VIPs_PLSDA_baseline_feces_female_age12_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_feces_female_age12_Load <- VIPs_PLSDA_baseline_feces_female_age12_select %>% 
  left_join(Loadings_PLSDA_baseline_feces_female_age12, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_feces_female_age12_Load_Top100 <- VIPs_PLSDA_baseline_feces_female_age12_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_feces_female_age12_plot <- plotLoadings(PLSDA_baseline_feces_female_age12, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_feces_female_age12_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_feces_female_age12") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_feces_female_age12_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_feces_female_age12_Load_cytoscape <- VIPs_PLSDA_baseline_feces_female_age12_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_feces_female_age12_Load, "PLSDA_baseline_feces_female_age12_Load.csv")
write_csv(VIPs_PLSDA_baseline_feces_female_age12_Load_cytoscape, "VIPs_PLSDA_baseline_feces_female_age12_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_feces_female_age12_plot, filename = "Loadings_PLSDA_baseline_feces_female_age12.png", device = "png", dpi = "retina")

# END OF BASELINE feces FEMALE AGE12 - PLS-DA & VIP PLOTS ----------------------


# END OF BASELINE FECES  --------------------------------------------------


# BASELINE SERUM - ONLY HAS (9 months not found in metadata) 12 MONTHS -------------------------------


# BASELINE SERUM GENOTYPE, GENDERS AND AGES (12 MONTHS) ALL ------------------------------------------------------

fbmn_stats_norm_filt_baseline_serum_all <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "control diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_serum_all <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "control diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  #metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_serum_all$combined_factor <- factor(
  paste(metadata_filtered_baseline_serum_all$mouse_genotype,
        metadata_filtered_baseline_serum_all$ATTRIBUTE_mouse_gender,
        metadata_filtered_baseline_serum_all$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_serum_all <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_serum_all, var = "filename"),
  Y = metadata_filtered_baseline_serum_all$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_serum_all_scores <- data.frame(PLSDA_baseline_serum_all$variates$X, metadata_filtered_baseline_serum_all)

# Plot PLS-DA results
PLSDA_baseline_serum_all_plot <- PLSDA_baseline_serum_all_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (serum, male, age 12 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_serum_all$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_serum_all$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = TRUE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_serum_all_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_serum_all_plot)
ggsave(plot = PLSDA_baseline_serum_all_plot, filename = "PLSDA_baseline_serum_all_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE serum ALL - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_serum_all <- plotLoadings(PLSDA_baseline_serum_all, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_serum_all <- perf(PLSDA_baseline_serum_all, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_serum_all, legend = TRUE, legend.position = "center")

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_serum_all <- as.data.frame(mixOmics::vip(PLSDA_baseline_serum_all))
VIPs_PLSDA_baseline_serum_all_filter <- VIPs_PLSDA_baseline_serum_all %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_serum_all_filter$ID <- rownames(VIPs_PLSDA_baseline_serum_all_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_serum_all_select <- VIPs_PLSDA_baseline_serum_all_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_serum_all_Load <- VIPs_PLSDA_baseline_serum_all_select %>% 
  left_join(Loadings_PLSDA_baseline_serum_all, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_serum_all_Load_Top100 <- VIPs_PLSDA_baseline_serum_all_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_serum_all_plot <- plotLoadings(PLSDA_baseline_serum_all, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_serum_all_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_serum_all") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_serum_all_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_serum_all_Load_cytoscape <- VIPs_PLSDA_baseline_serum_all_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_serum_all_Load, "PLSDA_baseline_serum_all_Load.csv")
write_csv(VIPs_PLSDA_baseline_serum_all_Load_cytoscape, "VIPs_PLSDA_baseline_serum_all_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_serum_all_plot, filename = "Loadings_PLSDA_baseline_serum_all.png", device = "png", dpi = "retina")

# END OF BASELINE SERUM ALL - PLS-DA & VIP PLOTS ----------------------


# BASELINE GENOTYPE SERUM MALE AGE12 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" ---------------------------

fbmn_stats_norm_filt_baseline_serum_male_age12 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "control diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_serum_male_age12 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "control diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "female" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "9" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_serum_male_age12$combined_factor <- factor(
  paste(metadata_filtered_baseline_serum_male_age12$mouse_genotype,
        metadata_filtered_baseline_serum_male_age12$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_serum_male_age12 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_serum_male_age12, var = "filename"),
  Y = metadata_filtered_baseline_serum_male_age12$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_serum_male_age12_scores <- data.frame(PLSDA_baseline_serum_male_age12$variates$X, metadata_filtered_baseline_serum_male_age12)

# Plot PLS-DA results
PLSDA_baseline_serum_male_age12_plot <- PLSDA_baseline_serum_male_age12_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (serum, male, age 12 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_serum_male_age12$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_serum_male_age12$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_serum_male_age12_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_serum_male_age12_plot)
ggsave(plot = PLSDA_baseline_serum_male_age12_plot, filename = "PLSDA_baseline_serum_male_age12_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE serum MALE AGE12 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_serum_male_age12 <- plotLoadings(PLSDA_baseline_serum_male_age12, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_serum_male_age12 <- perf(PLSDA_baseline_serum_male_age12, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_serum_male_age12, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_serum_male_age12 <- as.data.frame(mixOmics::vip(PLSDA_baseline_serum_male_age12))
VIPs_PLSDA_baseline_serum_male_age12_filter <- VIPs_PLSDA_baseline_serum_male_age12 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_serum_male_age12_filter$ID <- rownames(VIPs_PLSDA_baseline_serum_male_age12_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_serum_male_age12_select <- VIPs_PLSDA_baseline_serum_male_age12_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_serum_male_age12_Load <- VIPs_PLSDA_baseline_serum_male_age12_select %>% 
  left_join(Loadings_PLSDA_baseline_serum_male_age12, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_serum_male_age12_Load_Top100 <- VIPs_PLSDA_baseline_serum_male_age12_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_serum_male_age12_plot <- plotLoadings(PLSDA_baseline_serum_male_age12, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_serum_male_age12_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_serum_male_age12") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_serum_male_age12_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_serum_male_age12_Load_cytoscape <- VIPs_PLSDA_baseline_serum_male_age12_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_serum_male_age12_Load, "PLSDA_baseline_serum_male_age12_Load.csv")
write_csv(VIPs_PLSDA_baseline_serum_male_age12_Load_cytoscape, "VIPs_PLSDA_baseline_serum_male_age12_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_serum_male_age12_plot, filename = "Loadings_PLSDA_baseline_serum_male_age12.png", device = "png", dpi = "retina")

# END OF BASELINE SERUM AGE12 - PLS-DA & VIP PLOTS ----------------------

# BASELINE GENOTYPE SERUM FEMALE AGE12 - ALL under ATTRIBUTE_mouse_treatment = "normal diet" ---------------------------

fbmn_stats_norm_filt_baseline_serum_female_age12 <- fbmn_stats_norm %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "control diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

metadata_filtered_baseline_serum_female_age12 <- metadata_filtered %>% 
  dplyr::filter(metadata_filtered$ATTRIBUTE_sample_type != "not applicable" &
                  metadata_filtered$ATTRIBUTE_sample_type != "feces" &
                  metadata_filtered$ATTRIBUTE_sample_type != "urine" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "Dapensutrile" &
                  metadata_filtered$ATTRIBUTE_mouse_treatment != "vehicle control" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "glycine diet" &
                  #metadata_filtered&ATTRIBUTE_mouse_treatment != "control diet" &
                  #metadata_filtered$ATTRIBUTE_mouse_treatment != "normal diet" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_gender != "male" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "not applicable" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "2" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "4" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "6" &
                  metadata_filtered$ATTRIBUTE_mouse_age != "12" &
                  metadata_filtered$mouse_genotype != "not applicable" &
                  metadata_filtered$mouse_genotype != "wild type (c57bl)" &
                  metadata_filtered$mouse_genotype != "ctns knockout (cystinosis mouse model)"
  )

# Create a combined factor for the response variable
metadata_filtered_baseline_serum_female_age12$combined_factor <- factor(
  paste(metadata_filtered_baseline_serum_female_age12$mouse_genotype,
        metadata_filtered_baseline_serum_female_age12$ATTRIBUTE_mouse_age, sep = "_")
)

# PLS-DA analysis with the combined factor as the response variable
PLSDA_baseline_serum_female_age12 <- mixOmics::plsda(
  X = column_to_rownames(fbmn_stats_norm_filt_baseline_serum_female_age12, var = "filename"),
  Y = metadata_filtered_baseline_serum_female_age12$combined_factor,
  ncomp = 3,
  scale = TRUE
)

# Extract scores
PLSDA_baseline_serum_female_age12_scores <- data.frame(PLSDA_baseline_serum_female_age12$variates$X, metadata_filtered_baseline_serum_female_age12)

# Plot PLS-DA results
PLSDA_baseline_serum_female_age12_plot <- PLSDA_baseline_serum_female_age12_scores %>%
  ggscatter(
    x = "comp1",
    y = "comp2",
    color = "mouse_genotype",
    #shape = "ATTRIBUTE_mouse_age",
    alpha = 0.6,
    title = "PLS-DA: baseline study (serum, female, age 12 months)",
    xlab = paste("Component 1 (", round(PLSDA_baseline_serum_female_age12$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
    ylab = paste("Component 2 (", round(PLSDA_baseline_serum_female_age12$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
    legend.title = "Group",
    ggtheme = theme_classic(),
    palette = c("#F8766D", "#619CFF", "#00BFC4", "#C77CFF")
  ) +
  #stat_ellipse(aes(comp1, comp2, linetype = ATTRIBUTE_mouse_age), type = "norm", level = 0.95, show.legend = FALSE, color = "black") +
  stat_ellipse(aes(comp1, comp2, color = mouse_genotype), type = "norm", level = 0.95, show.legend = FALSE) + # Confidence ellipses
  geom_point(data = PLSDA_baseline_serum_female_age12_scores %>% group_by(mouse_genotype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = mouse_genotype), size = 2, shape = 8) +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 5)
  ) +
  
  coord_fixed()

plot(PLSDA_baseline_serum_female_age12_plot)
ggsave(plot = PLSDA_baseline_serum_female_age12_plot, filename = "PLSDA_baseline_serum_female_age12_plot.png", device = "png", dpi = "retina")

# BASELINE GENOTYPE serum FEMALE AGE12 - PLS-DA & VIP PLOTS ---------------------------

# Extract loadings
Loadings_PLSDA_baseline_serum_female_age12 <- plotLoadings(PLSDA_baseline_serum_female_age12, plot = TRUE, contrib = "max") %>%
  rownames_to_column(var = "rowname") %>% 
  dplyr::select(rowname, GroupContrib, importance)

# Perform performance evaluation
perf_PLSDA_baseline_serum_female_age12 <- perf(PLSDA_baseline_serum_female_age12, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_PLSDA_baseline_serum_female_age12, legend = FALSE)

# Calculate VIP scores and filter by threshold
VIPs_PLSDA_baseline_serum_female_age12 <- as.data.frame(mixOmics::vip(PLSDA_baseline_serum_female_age12))
VIPs_PLSDA_baseline_serum_female_age12_filter <- VIPs_PLSDA_baseline_serum_female_age12 %>% dplyr::filter(comp1 > 1)
VIPs_PLSDA_baseline_serum_female_age12_filter$ID <- rownames(VIPs_PLSDA_baseline_serum_female_age12_filter)

# Select relevant VIPs and merge with loadings
VIPs_PLSDA_baseline_serum_female_age12_select <- VIPs_PLSDA_baseline_serum_female_age12_filter %>% dplyr::select(ID, comp1)
VIPs_PLSDA_baseline_serum_female_age12_Load <- VIPs_PLSDA_baseline_serum_female_age12_select %>% 
  left_join(Loadings_PLSDA_baseline_serum_female_age12, by = c("ID" = "rowname")) %>% 
  arrange(desc(comp1))

# Select top 100 VIPs
VIPs_PLSDA_baseline_serum_female_age12_Load_Top100 <- VIPs_PLSDA_baseline_serum_female_age12_Load %>% head(100)

# Plot top 100 loadings
Loadings_PLSDA_baseline_serum_female_age12_plot <- plotLoadings(PLSDA_baseline_serum_female_age12, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_PLSDA_baseline_serum_female_age12_Load_Top100$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "descending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLSDA_baseline_serum_female_age12") + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot(Loadings_PLSDA_baseline_serum_female_age12_plot, legend = TRUE)

## Adjust VIP Loading table to use as input in Cytoscape
VIPs_PLSDA_baseline_serum_female_age12_Load_cytoscape <- VIPs_PLSDA_baseline_serum_female_age12_Load %>%
  dplyr::mutate(featID = stringr::str_split_fixed(ID, "_", 2)[, 1])

# Save results

write_csv(VIPs_PLSDA_baseline_serum_female_age12_Load, "PLSDA_baseline_serum_female_age12_Load.csv")
write_csv(VIPs_PLSDA_baseline_serum_female_age12_Load_cytoscape, "VIPs_PLSDA_baseline_serum_female_age12_Load_cytoscape.csv")

ggsave(plot = Loadings_PLSDA_baseline_serum_female_age12_plot, filename = "Loadings_PLSDA_baseline_serum_female_age12.png", device = "png", dpi = "retina")

# END OF BASELINE SERUM AGE12 - PLS-DA & VIP PLOTS ----------------------

