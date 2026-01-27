# Cross-species Immune Transcriptomic Analysis

## *Nematostella vectensis* and *Stylophora pistillata*

This folder contains data, code, and results for a comparative transcriptomic analysis of immune responses in the cnidarian species *Nematostella vectensis* and *Stylophora pistillata*. The analysis focuses on gene- and orthogroup-level regulation across species, with an emphasis on immune genes and lineage-specific expansions.

------------------------------------------------------------------------

## Overview

The main goals of this analysis are to:

-   Compare differential gene expression responses to immune stimulation across species
-   Integrate datasets using orthogroups and orthology relationships
-   Assess concordance of transcriptional responses using multiple aggregation strategies
-   Characterize ortholog multiplicity per gene and its impact on cross-species comparisons

------------------------------------------------------------------------

## Dataset Description

### *Nematostella vectensis*

Transcriptomic data derived from immune stimulation experiments in *N. vectensis*.\
These experiments include exposure to innate immune stimulation with dsRNA designed to activate antiviral pathways, followed by differential expression analysis.

Results are summarized at the gene level and integrated with orthogroup information for cross-species comparison.

------------------------------------------------------------------------

### *Stylophora pistillata*

Transcriptomic data from *S. pistillata* immune activation experiments described in:

> **Li et al., 2023.**\
> *cGLRs are a diverse family of pattern recognition receptors in animal innate immunity.*\
> **Cell**

In this dataset, corals were stimulated with cyclic dinucleotide immune ligands:

-   **2′3′-cGAMP - for the comparison we used this condition which was the most potent.**
-   **3′3′-cUA**

These ligands activate conserved innate immune signaling pathways downstream of cGAS-like receptors. Differential expression was computed relative to control conditions.

------------------------------------------------------------------------

## Folder Structure

```         
├── data/
│ ├── mast_results.csv # Motif / regulatory analysis results
│ ├── Orthogroups.GeneCount.tsv # Orthogroup gene counts
│ └── Spi_DE_2.3.csv # Stylophora differential expression results
│
├── script/
│ ├── cross_species_analysis.R # Main analysis script
│ ├── cross_species_analysis.md # Analysis report (Markdown)
│ ├── cross_species_analysis.html # Rendered analysis report
│ └── cross_species_analysis_files/ # Supporting files for report figures
│
├── results/
│ ├── orthogroup_classified.csv
│ ├── shared_vs_specific.csv
│ ├── shared_upregulated_orthogroups_new.txt
│ ├── shared_downregulated_orthogroups_new.txt
│ ├── orthogroup_venn_diagrams.png
│ ├── ortholog_correlation_full_labeled.png
│ └── ortholog_correlation_full_labeled.pdf
│
└── README.md
```

------------------------------------------------------------------------

## Analysis Summary

-   Differential expression results from both species are integrated using orthogroups.
-   Orthogroup-level data are collapsed to the gene level using:
    -   **Median aggregation** (typical response)
    -   **Max aggregation** (peak inducible response)
-   Correlations of log2 fold changes across species are computed under multiple filtering strategies.
-   The distribution of ortholog counts per gene reveals strong lineage-specific expansions.

------------------------------------------------------------------------

## Reproducibility

The full analysis can be reproduced by running:

\`\`\`bash Rscript script/cross_species_analysis.R
