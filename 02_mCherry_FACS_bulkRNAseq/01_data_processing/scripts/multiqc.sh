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
