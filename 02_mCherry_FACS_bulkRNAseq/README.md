# RNA-seq Analysis Pipeline

This document outlines a complete RNA-seq workflow, from raw read download to downstream differential expression and functional enrichment analysis.\
Each section provides placeholders for commands, parameters, and notes to ensure a reproducible and well-documented analysis.

------------------------------------------------------------------------

## Table of Contents

1.  [1. SRA Toolkit]
2.  [Renaming FASTQ Files](#2-renaming-fastq-files)
3.  [FastQC](#3-fastqc)
4.  [MultiQC](#4-multiqc)
5.  [Trimmomatic](#5-trimmomatic)
6.  [STAR Alignment](#6-star-alignment)
7.  [featureCounts](#7-featurecounts)
8.  [Downstream Analysis in R](#8-downstream-analysis-in-r)
9.  [Summary and File Structure](#9-summary-and-file-structure)

------------------------------------------------------------------------

## 1. SRA Toolkit

**Purpose:**\
Download raw sequencing data from NCBI SRA as `.fastq.gz` files.

------------------------------------------------------------------------

``` bash

#!/bin/bash
#SBATCH -c 4
#SBATCH --mem-per-cpu 4000
#SBATCH --time=96:00:00
#SBATCH -J wget
#SBATCH -o Download.%A.out
#SBATCH -e Download.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<itamar.kozlovski@mail.huji.ac.il>

module load sratoolkit

while read -r SRA_ID; do
    echo "Downloading $SRA_ID"
    prefetch "$SRA_ID"
    echo "Converting $SRA_ID to Fastq"
    fasterq-dump "$SRA_ID"
done < sra_list.txt
```

## 2. Renaming FASTQ files

This command will rename the files according to their mCherry status (e.g. neg1 = mCherry negative cells replicate 1; pos2 = mCherry positive cells replicate 2 etc.)

```{bash}

```
