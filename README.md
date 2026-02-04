# cnidarian_immune_cells

Code and analysis workflows to reproduce the study\
*Functional characterization of specialized immune cells in a cnidarian reveals an ancestral antiviral program*.

This repository contains the complete computational pipelines and analysis scripts used in the study, integrating bulk RNA sequencing, single-cell RNA sequencing, imaging cytometry, proteomics, and cross-species comparative analyses in *Nematostella vectensis*.

All analyses are provided as documented, reproducible workflows, with rendered reports (`.md` / `.html`) accompanying the source code.

Large input files and processed objects are deposited separately on Zenodo and are linked from the relevant analysis directories.

------------------------------------------------------------------------

## Repository Structure

### 01_Immuno_FACS_bulkRNAseq

Bulk RNA-seq analysis of FACS-isolated immune cell populations, including data preprocessing, quality control, differential expression, and functional enrichment analyses.

### 02_mCherry_FACS_bulkRNAseq

Bulk RNA-seq analysis of mCherry-labeled cell populations, structured analogously to `01_Immuno_FACS_bulkRNAseq` and analyzed independently.

### 03_scRNAseq

Single-cell RNA-seq analysis workflows, including: - Metacell construction and analysis - Seurat-based clustering and cell type annotation - Pseudotime and trajectory inference - Cell Ranger quality control summaries

These analyses underpin cell type identification and functional characterization at single-cell resolution.

### 04_Orthology_go_pipeline

Orthology inference and Gene Ontology annotation pipelines used for functional interpretation and cross-species comparisons.

### 05_Imaging_cytometry_analysis

Analysis of imaging flow cytometry data, including preprocessing, dimensionality reduction, visualization, and quantitative comparison of immune cell populations.

### 06_Mass_spec_analysis

Proteomics analysis workflow based on mass spectrometry data, including statistical analysis and result visualization.

### 07_Cross_species_comparison

Comparative analyses integrating orthology, bulk RNA-seq, and single-cell results to identify conserved and lineage-specific immune programs.

### 08_phagocytosis

Analysis of phagocytosis assays across different cell fractions and experimental conditions.\
This directory contains statistical analyses and visualization scripts used to quantify phagocytic activity and generate the corresponding figures.

### 09_supplementary_files

Supplementary annotation tables and reference files used across multiple analyses.

------------------------------------------------------------------------

## Reproducibility and Data Availability

-   Analyses are primarily implemented in **R**, with additional **Python** and **shell** scripts where required.
-   Each major analysis directory contains a dedicated `README.md` describing its contents and execution logic.
-   Rendered reports provide a transparent record of parameters, intermediate steps, and outputs.
-   Large datasets, intermediate objects, and processed results are deposited on Zenodo and are referenced from the corresponding analysis directories.
-   Raw and processed **single-cell RNA sequencing data** are available via Zenodo (DOI: https://doi.org/10.5281/zenodo.18390360) and are linked from the `03_scRNAseq` directory.
-   **Cytometry data**, including imaging flow cytometry and phagocytosis assay datasets, are deposited in a separate Zenodo repository (DOI: https://doi.org/10.5281/zenodo.18212341) and are linked from the relevant analysis directories.

------------------------------------------------------------------------

## License

This repository is distributed under the terms specified in the `LICENSE` file.
