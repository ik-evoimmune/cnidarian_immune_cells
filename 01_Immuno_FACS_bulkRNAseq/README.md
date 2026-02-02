# Bulk RNA-seq analysis of FACS-sorted RLRb populations

This directory contains the code used to analyze bulk RNA-seq data generated from FACS-sorted RLRb-high and RLRb-low cells following intracellular immunostaining.

Raw sequencing reads generated for this study are publicly available via the NCBI Sequence Read Archive (SRA) under BioProject accession number **PRJNA1207058**.

## Directory Overview

-   `01_data_processing/`
    -   Raw data handling and preprocessing steps, including input files, trimming resources, and data preparation.
    -   Includes an HTML/Markdown report documenting preprocessing.
-   `02_differential_expression_analysis/`
    -   Differential expression analysis using DESeq2.
    -   Exploratory analyses (PCA, heatmaps, volcano plots).
    -   Functional interpretation using GSEA with clusterProfiler.

Each subdirectory contains its own README describing the analysis steps and outputs in detail.
