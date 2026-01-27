---
title: "RNA-seq Analysis Pipeline"
author: "Itamar Kozlovski"
date: "2026-01-15"
---

# RNA-seq Analysis Pipeline

This document outlines a complete RNA-seq workflow, from raw read download to downstream differential expression and functional enrichment analysis.\
Each section provides placeholders for commands, parameters, and notes to ensure a reproducible and well-documented analysis.

------------------------------------------------------------------------

## Table of Contents

1.  [1. SRA Toolkit]
2.  [2. QC using FASTQC]
3.  [3. Quality control summary with MultiQC]
4.  [4. Read trimming with Trimmomatic]
5.  [5. Star genome index generation]
6.  [6. Star alignment]
7.  [7. Read count estimation]
8.  [8. Downstream analysis in R]

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

``` bash
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

## 3. Quality control summary with MultiQC

Sequencing quality control reports were aggregated using **MultiQC**. The script scans the current directory for supported QC outputs (e.g. FastQC, STAR logs) and generates a consolidated summary report for downstream inspection.

``` bash
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

## 4. Read trimming with Trimmomatic

This script performs preprocessing of sequencing data, including adapter removal and quality trimming, using Trimmomatic. All FASTQ files in the directory are processed in single-end mode with standard quality-control parameters to generate cleaned reads suitable for downstream analyses.

``` bash
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=24:00:00
#SBATCH -J Trimmomatic
#SBATCH -o trimmo.%A.out
#SBATCH -e trimmo.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your.email@domain.com>

module load trimmomatic
for file in *.fastq; do trimmomatic SE -phred33 $file "`basename $file .fastq`.trimmomatic_out.fastq" ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
 done
```

## 5. Star genome index generation

Prior to genome indexing, the **mCherry reporter sequence `data/mCherry.fasta`** was concatenated to the *Nematostella vectensis* reference genome FASTA, and corresponding gene annotations were added to the GTF file. This augmented reference allows accurate alignment and quantification of reads derived from the mCherry transgene.

Genome indexing was then performed using **STAR** in `genomeGenerate` mode. The index was built with 8 threads using the concatenated genome FASTA and modified GTF annotation, with splice junctions generated using a read-length–appropriate overhang (`sjdbOverhang = 84`). The resulting indexed genome was used for downstream RNA-seq read alignment.

``` bash
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 20000
#SBATCH --time=24:00:00
#SBATCH -J STAR
#SBATCH -o STAR.%A.out
#SBATCH -e STAR.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your.email@domain.com>

module load star

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /path/to/project/indexed_genome \
  --genomeFastaFiles /path/to/project/genome/genome_with_transgene.fasta \
  --sjdbGTFfile /path/to/project/genome/annotation_with_transgene.gtf \
  --sjdbOverhang 84 \
  --genomeSAindexNbases 13
```

## 6. Star alignment

This script performs read alignment using STAR on **trimmed FASTQ files** generated from the preprocessing step. Each trimmed FASTQ file is aligned independently to a pre-built STAR genome index, and the output is a **coordinate-sorted BAM file** suitable for downstream RNA-seq analyses.

The script is designed to run on a SLURM-based high-performance computing cluster and assumes that trimming and genome indexing have been completed prior to execution.

``` bash
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=24:00:00
#SBATCH -J STAR_align
#SBATCH -o STAR.%A.out
#SBATCH -e STAR.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your.email@domain.com>

# Load STAR module
module load star

# Define paths
GENOME_DIR=/path/to/STAR_index
INPUT_DIR=/path/to/fastq_files
OUTPUT_DIR=/path/to/alignment_output

# Create output directory if it does not exist
mkdir -p "${OUTPUT_DIR}"

# Run STAR alignment for each FASTQ file
for file in "${INPUT_DIR}"/*out.fastq; do
    BASENAME=$(basename "${file}" .fastq)

    STAR \
        --runMode alignReads \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${file}" \
        --runThreadN 8 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${OUTPUT_DIR}/${BASENAME}_"
done
```

## 7. Read count estimation

This script quantifies aligned RNA-seq reads using **featureCounts**. The input consists of **coordinate-sorted BAM files** generated from the STAR alignment step. Reads are assigned to transcript features based on a provided GTF annotation file, and the output is a **count matrix** that can be used for downstream differential expression analysis.

The script is designed to run on a SLURM-based high-performance computing cluster and assumes that read alignment has been completed prior to execution.

``` bash
#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=24:00:00
#SBATCH -J featureCounts
#SBATCH -o featureCounts.%A.out
#SBATCH -e featureCounts.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your.email@domain.com>

# Load subread / featureCounts module
module load subread

# Define paths
ANNOTATION_GTF=/path/to/annotation.gtf
BAM_DIR=/path/to/bam_files
OUTPUT_FILE=counts_matrix.txt

# Run featureCounts (transcript-level counting)
featureCounts \
    -a "${ANNOTATION_GTF}" \
    -o "${OUTPUT_FILE}" \
    -g transcript_id \
    "${BAM_DIR}"/*.bam
```

## 8. Downstream analysis in R

All downstream analyses were performed in **R**, including differential expression analysis and data visualization. Scripts, processed results, and additional documentation are provided in the `02_differential_expression_analysis` directory.
