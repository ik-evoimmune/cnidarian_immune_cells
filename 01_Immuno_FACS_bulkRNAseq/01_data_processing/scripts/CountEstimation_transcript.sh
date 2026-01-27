#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 2000
#SBATCH --time=24:00:00
#SBATCH -J FeatureCounts
#SBATCH -o FeatureCount.%A.out
#SBATCH -e FeatureCounts.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<itamar.kozlovski@mail.huji.ac.il>




/sci/labs/yehum79/itamar273/subread/subread-2.0.7-source/bin/./featureCounts \
-a /sci/home/itamar273/immune_cells/scRNAseq_analysis/gtf_file/Nvec_v4_merged_annotation_sort_mCherry.gtf \
-o count_matrix_updated.txt \
-g transcript_id \
 *.bam 
