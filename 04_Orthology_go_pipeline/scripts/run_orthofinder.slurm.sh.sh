#!/bin/bash
#SBATCH -c 12
#SBATCH --mem-per-cpu 2000
#SBATCH --time=96:00:00
#SBATCH -J orthofinder
#SBATCH -o orthofinder.%A.out
#SBATCH -e orthofinder.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<ajjb80@gmail.com>

eval "$(/sci/labs/yehum79/adrianjjb/miniconda3/bin/conda shell.zsh hook)"
conda activate orthofinder

orthofinder -f  /sci/labs/yehum79/adrianjjb/Itamar_Orthology/
