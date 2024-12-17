# cystinosis-baseline-mouse-model-manuscript
LC-MS/MS metabolomics profile phenotype from gut (fecal), urine and serum on Ctns -/- mice versus wild type controls (c57BL/6J). 

This repository contains data processing steps from a subset of LC-MS/MS metabolomics dataset using the following softwares and resources:
* A subset of LC-MS/MS data corresponding to the baseline triple metabolomics (fecal, urine and serum) profile of a mouse model Ctns -/- mice versus wild type controls (c57BL/6J).
* The baseline data is part of the public dataset [MSV000093184](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=fccbbebb0870401bbfc3199c10df6c3f).
* The baseline subset has been processed using [MZmine4](https://mzio.io/) version 4.0.3. from the MZmine4 outputs using the following workflows:
* Feature-based Molecular Networking in [GNPS2](https://gnps2.org/homepage)
* Data cleaning (blank removal, imputation and normalization) in [FBMN-stats](https://fbmn-statsguide.gnps2.org/)
* Statistical analysis using [MixOmics](https://mixomics.org/) package in [Rstudio](https://www.rstudio.com/products/rstudio/download/).

The `MZmine4_outputs/` directory contains `.mgf` and `.csv` files used as inputs for FBMN and stats. 

The `Scripts/` directory contains `.R` files used for cleaning dataframes, generate tables and graphics from the statistical analysis. 
It requires R and RStudio. 

* R: [R Download](https://cran.r-project.org/bin/).
* RStudio: [RStudio Download](https://www.rstudio.com/products/rstudio/download/).

# Table of contents
1. Data processing with [MZmine4](https://mzio.io/).
2. Feature-based Molecular Networking (FBMN)
3. Data cleaning 
4. Statistical analysis using mixOmics
5. Additional outputs to map into networks
6. Cytoscape visualization
7. Results

# 1. Data processing with MZmine4

With this software, the input files necessary for performing FBMN in GNPS2 are generated and described below: 

* Batch: the `.mzbatch` file contains the settings used to process the data (`.mzML` files can be obtained from the dataset [MSV000093184](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=fccbbebb0870401bbfc3199c10df6c3f)).
* Feature table: a `dataset_quant.csv` file corresponds to the MZmine4 output containing the features in rows, samples in columns and peak areas as cell values. This file is the "inputfeatures" used in FBMN GNPS2.
* Edges annotation: a `dataset_edges_msannotation.csv` file corresponds to the MZmine4 output containing the annotations of edges from Ion Identity Molecular Networking (some ions might correspond to adducts of the same molecule and this helps with annotating them).
* Spectral data: a `dataset.mgf` contains the MS1 (precursor ion) and MS2 (fragmentation spectra) information. This is the file used as "inputspectra" in FBMN GNPS2.

**Note:** By **features** we mean detected ion MS1, fragmentation spectrum MS2, retention time (rt) from chromatography, and peak area.

For more details and tutorials, visit [MZmine4](https://mzio.io/) and the following publication, although on MZmine3, provides a step-by-step guide for data processing which remains applicable to MZmine4: 

[Heuckeroth, S., Damiani, T., Smirnov, A. et al. Reproducible mass spectrometry data processing and compound annotation in MZmine 3. Nat Protoc 19, 2597–2641 (2024)](https://doi.org/10.1038/s41596-024-00996-y).

# 2. Feature-based Molecular Networking (FBMN)

FBMN provides the network where features (detected ion MS1, fragmentation spectrum MS2, retention time (rt) from chromatography, and peak area) are clustered based on spectral similarity. Chemically related molecules, represented as nodes in a network, are connected forming molecular families. In addition to the molecular networking, a library search is performed providing annotations to the detected features. 

The FBMN job (task=e92c0b19522946af89db73401c39d672) in GNPS2 is available through the link [FBMN job](https://gnps2.org/status?task=e92c0b19522946af89db73401c39d672)

You can visualize the FBMN results by clicking "Library Results" (to display spectral library matches) or visualize the entire molecular network by clicking "Visualize Full Network in Browser" ("Visualize Full Network w/Singletons in Browser" will allow you to visualize nodes that are not connected to any other one). 

For more details about molecular networking, feature-based molecular networking and ion identity, please read the following publications: 
[Aron, A.T., Gentry, E.C., McPhail, K.L. et al. Reproducible molecular networking of untargeted mass spectrometry data using GNPS. Nat Protoc 15, 1954–1991 (2020)](https://doi.org/10.1038/s41596-020-0317-5)

[Nothias, LF., Petras, D., Schmid, R. et al. Feature-based molecular networking in the GNPS analysis environment. Nat Methods 17, 905–908 (2020)](https://doi.org/10.1038/s41592-020-0933-6). 

[Schmid, R., Petras, D., Nothias, LF. et al. Ion identity molecular networking for mass spectrometry-based metabolomics in the GNPS environment. Nat Commun 12, 3832 (2021)](https://doi.org/10.1038/s41467-021-23953-9).

# 3. Data cleaning

By data cleaning, the following steps, which can be optional, were applied: Blank removal, imputation and normalization. This can be done using [FBMN-stats](https://fbmn-statsguide.gnps2.org/) by providing the FBMN job task (e92c0b19522946af89db73401c39d672) under GNPS task ID and by clicking "Load files from GNPS" (red bottom)

* Once the quantification table (feature table) and metadata are retrieved, **Data Cleanup** can be performed.
* **Samples** from the "attribute for sample selection" are selected by "SampleType", "sample selection" = "Animal"
* **Blanks** from the "attribute for blank selection" are selected by "SampleType", "blank selection" = "blank_analysis", "blank_extraction", "blank_QC"
* cutoff threshold for blank removal: 0.30
* Imputation: These values will be filled with random number between 1 and 171 (Limit of Detection) during imputation.
* Normalization: Total Ion Current (TIC) or sample-centric normalization
* Submit and download as .csv
* *Optional* Quick statistical analysis and overview of the data can be done in the fbmn-statsguide app

These data cleaning step provide the `fbmn_stats_data_clean_export.csv` used for downstream process. Due to his size (>20Mb), this file is available in the [MSV000093184](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=fccbbebb0870401bbfc3199c10df6c3f). You can access it and download directly or copy/paste the ftp URL to your computer browser: ftp://massive.ucsd.edu/v06/MSV000093184//updates/2024-12-16_amcaraballor_398a311b/other/fbmn_stats_and_metadata/fbmn_stats_data_clean_export.csv

**Metadata**: A metadata `merged_metadata.tsv` is provided under the `MixOmics`.  this file is available in the [MSV000093184](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=fccbbebb0870401bbfc3199c10df6c3f). You can access it and download directly or copy/paste the ftp URL to your computer browser: ftp://massive.ucsd.edu/v06/MSV000093184//updates/2024-12-16_amcaraballor_398a311b/other/fbmn_stats_and_metadata/merged_metadata.tsv

# 4. Statistical analysis using mixOmics

After data cleaning, the mixOmics package is used for statistical analysis. The following Rscripts are provided under the `Scripts/` directory:

* MSV000093184_dataframe_cleaning.R
* MSV000093184_post_fbmn_stats.R
* MSV000093184_sum_loadings_for_cytoscape.R
* MSV000093184_boxplots_selected_features.R

# 5. Additional outputs to map into networks

Additional outputs are generated (e.g., loading and VIPs scores) to map into networks and visualize in [Cytoscape](https://cytoscape.org/).

These outputs are also provided under the `MixOmics/` folder.

# 6. Cytoscape visualization

Additional `FBMN_MSV000093184.cys` output already containing all the additional information created through the downstream analysis described in the previous steps is available under the `Cytoscape` folder.

# 7. Results

After the steps described above, one can visualize features (molecules) of interest that have been selected after the downstream statistical analysis. 


This example shows the full molecular network mapped by the sum of importance scores from the PLS-DA models in the baseline urine subset (higher and purple indicate high sum score). Spectral libraries provide a match to p-cresol glucuronide (_m/z_ 302.1241).  


Boxplot of the p-cresol glucuronide (featureID 1530, _m/z_ 302.1241, rt 2.41 min) by gender, age and phenotype.
