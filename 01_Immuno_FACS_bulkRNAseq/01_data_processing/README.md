# 01_data_processing

This directory contains the code used to analyze bulk RNA-seq data generated from FACS-sorted RLRb-high and RLRb-low cells following intracellular immunostaining.

## Directory structure

``` text
01_data_processing/
├── data/
│   └── Input files and reference resources used during preprocessing
├── scripts/
│   ├── download_sra.sh
│   ├── Trimmomatic.sh
│   ├── fastqc_alt.sh
│   ├── multiqc.sh
│   ├── STAR_Index.sh
│   ├── STAR_align.sh
│   └── CountEstimation_transcript.sh
├── 01_data_processing.md
└── 01_data_processing.html
```

## Processing steps

The following steps are performed in this directory:

### Data retrieval

Raw sequencing data are downloaded from public repositories using `download_sra.sh`.

### Read preprocessing

Adapter removal and quality trimming of single-end FASTQ files using Trimmomatic (`Trimmomatic.sh`).

### Quality control

Per-sample quality assessment using FastQC (fastqc_alt.sh) and aggregation of QC metrics using MultiQC (`multiqc.sh`).

### Reference genome preparation

Generation of a STAR genome index using a Nematostella vectensis reference genome augmented with the mCherry reporter sequence (`STAR_Index.sh`).

### Read alignment

Alignment of processed reads to the indexed genome using STAR (`STAR_align.sh`).

Read counting / quantification Transcript-level read count estimation for downstream analysis (`CountEstimation_transcript.sh`).

## Usage

Scripts are designed to be executed on an HPC cluster using SLURM and should be run from the project root directory, with file paths specified relative to the repository.

Example:

``` bash
sbatch 01_data_processing/scripts/STAR_align.sh
```

## Notes

Paths and environment initialization commands may require adaptation to local computing environments.

This directory contains preprocessing outputs only. Downstream statistical analyses and visualization in R are documented in `02_differential_expression_analysis`.
