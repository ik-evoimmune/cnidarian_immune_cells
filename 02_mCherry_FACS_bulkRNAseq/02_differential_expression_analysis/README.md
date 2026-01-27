# Differential Expression and Functional Analysis

This directory contains the bulk RNA-seq differential expression and downstream exploratory analyses comparing mCherry+ and mCherry- samples.

## Overview

The analysis workflow consists of the following steps:

1.  **Count-based differential expression analysis**
    -   Gene-level count matrix generated using FeatureCounts.
    -   Differential expression performed using **DESeq2**, including normalization and statistical testing.
    -   Output includes log2 fold changes, adjusted p-values, and significance calls.
2.  **Exploratory data analysis**
    -   **Principal Component Analysis (PCA)** to assess sample relationships and variance structure.
    -   **Heatmaps** of differentially expressed genes to visualize expression patterns across samples.
    -   **Volcano plots** to summarize effect size versus statistical significance.
3.  **Gene Set Enrichment Analysis (GSEA)**
    -   Ranked gene lists derived from DESeq2 results.
    -   Enrichment analysis performed using **clusterProfiler**.
    -   Outputs include enriched pathways and functional categories, along with summary tables and figures.

## Directory Structure

-   `data/`
    -   Raw and processed input/output files, including:
        -   FeatureCounts count matrix
        -   Differential expression tables
        -   Lists of up- and down-regulated genes
        -   GSEA results and summary figures
-   `script/`
    -   R scripts used for DESeq2 analysis, visualization, and GSEA

## Outputs

Key result files include: - `DE_table.csv` – full DESeq2 results - `mCherry_upregulated.txt` / `mCherry_downregulated.txt` – filtered gene lists - `Results_GSEA_New.csv` – GSEA enrichment results - `GSEA_metrics_table.pdf/png` – summary visualization of enriched terms - `S4_table.csv` – supplementary results table

This directory contains all analyses required to reproduce the differential expression, visualization, and functional interpretation presented in the manuscript.
