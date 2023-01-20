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

# computeMatrix reference-point -o ../results/deeptools/mESC_MEF_PQS.mat.gz -S \
 # ../data/bw/mm10_canPQS-regex_binary.bw \
 # -R  ../results/Seurat/callpeaks_mESC-MEF/peak_sets/0_peaks_res0.1.bed \
 # ../results/Seurat/callpeaks_mESC-MEF/peak_sets/1_peaks_res0.1.bed \
  # --referencePoint center \
 # -b 3000 -a 3000 --samplesLabel "ext. PQS" --skipZeros --missingDataAsZero --scale 100

 plotHeatmap -m ../results/deeptools/mESC_MEF_PQS.mat.gz \
 -out "../results/deeptools/mESC-MEF_clusters_PQS.pdf" \
 --refPointLabel "G4" \
 --whatToShow "heatmap and colorbar" \
 --heatmapHeight 14 \
 --colorMap "Greys" \
 --yMin 1 \
 --yMax 1 \
 --zMax 1 \
 -z "0" \
 --yAxisLabel "" \
 --xAxisLabel "" \


plotHeatmap -m ../results/deeptools/mESC_MEF_PQS.mat.gz \
 -out "../results/deeptools/mESC-MEF_clusters_PQS.png" \
 --refPointLabel "G4" \
 --whatToShow "heatmap and colorbar" \
 --heatmapHeight 14 \
 --colorMap "Greys" \
 --yMin 1 \
 --yMax 1 \
 --zMax 1 \
 -z "0" \
 --yAxisLabel "" \
 --xAxisLabel "" \
