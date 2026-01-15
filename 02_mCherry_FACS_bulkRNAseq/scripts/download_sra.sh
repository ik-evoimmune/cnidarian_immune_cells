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
