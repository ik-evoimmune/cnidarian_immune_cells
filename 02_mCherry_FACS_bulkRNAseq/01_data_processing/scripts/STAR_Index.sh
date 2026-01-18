#!/bin/bash
#SBATCH -c 8
#SBATCH --mem-per-cpu 20000
#SBATCH --time=24:00:00
#SBATCH -J STAR
#SBATCH -o STAR.%A.out
#SBATCH -e STAR.%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<itamar.kozlovski@mail.huji.ac.il>

module load star

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /sci/labs/yehum79/itamar273/immune_cells_RNAseq_New_Aanalysis/mCherry_RLRb_FACS/reanalyzed_07102024/indexed_genome --genomeFastaFiles /sci/labs/yehum79/itamar273/immune_cells_RNAseq_New_Aanalysis/mCherry_RLRb_FACS/reanalyzed_07102024/mCherry_genome/Nvec_vc1.1_gDNA.mCherry.fasta --sjdbGTFfile /sci/labs/yehum79/itamar273/immune_cells_RNAseq_New_Aanalysis/mCherry_RLRb_FACS/reanalyzed_07102024/mCherry_genome/Nvec_vc1.1_long.annot_mCherry.gtf --sjdbOverhang 84 --genomeSAindexNbases 13
