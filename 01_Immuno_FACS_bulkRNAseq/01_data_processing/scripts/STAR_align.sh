#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=24:00:00
#SBATCH -J STAR
#SBATCH -o STAR.%A.out
#SBATCH -e STAR.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<itamar.kozlovski@mail.huji.ac.il>

module load star
for file in *out.fastq; do STAR --runMode alignReads --genomeDir /sci/labs/yehum79/itamar273/immune_cells_RNAseq_New_Aanalysis/mCherry_RLRb_FACS/reanalyzed_07102024/indexed_genome --outSAMtype BAM SortedByCoordinate --readFilesIn ${file} --runThreadN 8 --outFileNamePrefix /sci/home/itamar273/immune_cells/mCherry_RLRb_FACS/reanalyzed_07102024/alignment/mapped${file}; done
