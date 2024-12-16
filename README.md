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

The `Scripts/` directory contains `.Rmd` files that generate the tables and graphics shown below. 
It requires R, RStudio, and the rmarkdown package. 

* R: [R Download](https://cran.r-project.org/bin/)
* RStudio: [RStudio Download](https://www.rstudio.com/products/rstudio/download/)
* rmarkdown can be installed using the install packages interface in RStudio

# Table of contents
1. Data processing in [MZmine4](https://mzio.io/).
2. [FBMN job](https://gnps2.org/status?task=e92c0b19522946af89db73401c39d672)
3. Blank removal, imputation and normalization. This can be done using [FBMN-stats](https://fbmn-statsguide.gnps2.org/) by providing the FBMN job task (e92c0b19522946af89db73401c39d672)
4. Statistical analysis using MixOmics
5. Generate input tables with loading and VIPs scores to map into networks and visualize in [Cytoscape](https://cytoscape.org/).

# 1. Data processing in MZmine4

With this software, the input files necessary for performing FBMN in GNPS2 are generated and described below: 

* Batch: the `.mzbatch` file contains the settings used to process the data (`.mzML` files can be obtained from the dataset [MSV000093184](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=fccbbebb0870401bbfc3199c10df6c3f)).
* Feature table: a `dataset_quant.csv` file corresponds to the MZmine4 output containing the features in rows, samples in columns and peak areas as cell values. This file is the "inputfeatures" used in FBMN GNPS2.
* Edges annotation: a `dataset_edges_msannotation.csv` file corresponds to the MZmine4 output containing the annotations of edges from Ion Identity Molecular Networking (some ions might correspond to adducts of the same molecule and this helps with annotating them). 

Note: By *features* we mean detected ion MS1, fragmentation spectrum MS2, retention time (rt) from chromatography, and peak area.

# 2. Feature-based Molecular Networking

FBMN provides the network where features (detected ion MS1, fragmentation spectrum MS2, retention time (rt) from chromatography, and peak area) are clustered based on spectral similarity. Chemically related molecules, represented as nodes in a network, are connected forming molecular families. In addition to the molecular networking, a library search is performed providing annotations to the detected features. 

[example molecular family]

For more details about molecular networking, feature-based molecular networking and ion identity, please read the following publications: 
[Aron, A.T., Gentry, E.C., McPhail, K.L. et al. Reproducible molecular networking of untargeted mass spectrometry data using GNPS. Nat Protoc 15, 1954–1991 (2020)](https://doi.org/10.1038/s41596-020-0317-5)

[Nothias, LF., Petras, D., Schmid, R. et al. Feature-based molecular networking in the GNPS analysis environment. Nat Methods 17, 905–908 (2020)](https://doi.org/10.1038/s41592-020-0933-6). 

[Schmid, R., Petras, D., Nothias, LF. et al. Ion identity molecular networking for mass spectrometry-based metabolomics in the GNPS environment. Nat Commun 12, 3832 (2021)](https://doi.org/10.1038/s41467-021-23953-9).
