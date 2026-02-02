# 03_scRNA-seq Analysis

This directory contains the single-cell RNA-seq (scRNA-seq) analysis workflow used in this study.\
The code here covers data loading, quality control, normalization, dimensionality reduction, clustering, and downstream analyses.

Due to file size constraints, raw inputs, intermediate objects, and processed results (e.g. `.rds` files) are **not stored in this GitHub repository** and are instead deposited on Zenodo.

## Data Availability

All scRNA-seq input files and processed results required to reproduce or explore the analysis are available at:

**Zenodo DOI:** <https://doi.org/10.5281/zenodo.18390360>

This includes, but is not limited to: - Seurat objects (`.rds`) - Count matrices and metadata - Intermediate analysis objects - Processed results used for figures and downstream analyses

## Directory Structure

Typical contents of this folder include: - R scripts for scRNA-seq preprocessing and analysis - Markdown / HTML reports generated from the analysis - Visualization outputs (e.g. PCA, UMAP, clustering plots)

Each script or subdirectory is documented inline or via accompanying README files where relevant.

## Analysis Overview

The scRNA-seq workflow generally includes: 1. Loading of processed count data 2. Quality control and filtering 3. Normalization and feature selection 4. Dimensionality reduction (PCA, UMAP) 5. Cell clustering and annotation 6. Downstream analyses (e.g. differential expression, module scoring)

Exact parameters and steps are documented in the scripts and rendered reports.

## Reproducibility Notes

-   Analysis was performed in **R** using **Seurat** and related packages.
-   Package versions and session information are provided in the rendered reports when applicable.
-   Users are encouraged to download the Zenodo archive and run the scripts in this directory to fully reproduce the analysis.

------------------------------------------------------------------------

For questions or issues related to this analysis, please refer to the corresponding scripts or open an issue in the repository.
