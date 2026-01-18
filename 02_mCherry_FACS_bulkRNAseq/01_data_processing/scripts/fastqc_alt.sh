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

