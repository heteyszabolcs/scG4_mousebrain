#!/bin/bash -l
#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M snowy
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools
module load R_packages

cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

# Rscript bw_heatmap_PQS.R -bw "../data/bw/canonical_PQS_mm10_binary.bw" -l "ext.\ PQS" -o "../results/deeptools/unsorted_cluster0-1_canPQS.pdf" \
    # -ymax 1 -yzmax 1
# Rscript bw_heatmap_PQS.R -bw "../data/bw/extended_canonical_PQS_mm10_binary.bw" -l "ext.\ PQS" -o "../results/deeptools/GSM3003547_ext_canPQS.png" \
    # -ymax 1 -yzmax 1
	
Rscript bw_heatmap_overlap.R -o "../results/deeptools/mESC-MEF-overlap.pdf" -ymax 100 -yzmax 100


