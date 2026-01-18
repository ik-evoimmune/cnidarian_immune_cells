# RNA-seq Analysis Pipeline

This document outlines a complete RNA-seq workflow, from raw read download to downstream differential expression and functional enrichment analysis.\
Each section provides placeholders for commands, parameters, and notes to ensure a reproducible and well-documented analysis.

------------------------------------------------------------------------

## Table of Contents

1.  [1. SRA Toolkit]
2.  [2. QC using FASTQC]
3.  [3. Summarizing the results using multiqc]
4.  [Trimmomatic](#5-trimmomatic)
5.  [STAR Alignment](#6-star-alignment)
6.  [featureCounts](#7-featurecounts)
7.  [Downstream Analysis in R](#8-downstream-analysis-in-r)
8.  [Summary and File Structure](#9-summary-and-file-structure)

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

## 2. QC using FASTQC

This script will preform quality check on the raw fastq files. The output will be save into a folder called QC.

```{bash}
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=3:00:00
#SBATCH -J fastqc
#SBATCH -o FastQC.%A.out
#SBATCH -e FastQC.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your.email@domain.com>

mkdir -p QC  # Create the QC directory if it doesn't already exist
for file in *.fastq
do
  fastqc -o QC $file
done

```

## 3. Summarizing the results using multiqc {data-link="1. SRA Toolkit"}

```{bash}
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=3:00:00
#SBATCH -J multiqc
#SBATCH -o multiQC.%A.out
#SBATCH -e multiQC.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your.email@domain.com>

eval "$(/path/to/your/miniconda3/bin/conda shell.zsh hook)"

conda activate multiqc


multiqc .

```
