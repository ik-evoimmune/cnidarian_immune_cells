#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=24:00:00
#SBATCH -J Trimmomatic
#SBATCH -o trimmo.%A.out
#SBATCH -e trimmo.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<itamar.kozlovski@mail.huji.ac.il>

module load trimmomatic
for file in *.fastq; do trimmomatic SE -phred33 $file "`basename $file .fastq`.trimmomatic_out.fastq" ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
 done
