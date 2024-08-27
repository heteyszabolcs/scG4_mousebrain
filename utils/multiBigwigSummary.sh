#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 16:00:00
#SBATCH -M snowy
#SBATCH -J bigwigSummary

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

# modules
module load bioinfo-tools
module load deepTools

multiBigwigSummary bins -b ../data/Marek_data/mESC-MEF/cluster_0.bw \
 ../data/Marek_data/mESC-MEF/cluster_1.bw \
 ../data/Marek_data/mESC-MEF/G4_H33WT_SL_CnT_R1.mm10.bw \
 ../data/Marek_data/mESC-MEF/G4_MEF_CnT_R1.mm10.bw \
 --BED ../data/Marek_data/mESC-MEF/regions.bed \
 -o ../results/deeptools/multibwsumm_results.npz

plotCorrelation -in ../results/deeptools/multibwsumm_results.npz --corMethod pearson --skipZeros  --whatToPlot heatmap -o ../results/deeptools/heatmap_PearsonCorr_bigwigScores.png --outFileCorMatrix ../results/deeptools/PearsonCorr_bigwigScores.tab

