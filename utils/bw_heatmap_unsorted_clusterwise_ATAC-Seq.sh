#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools
module load R_packages

cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

computeMatrix reference-point -o ../results/deeptools/matrix_unsorted_clusterwise_ATAC.mat.gz -S \
 ../data/bw/ATAC-H33WT.mm10.bw \
 -R ../results/Seurat/callpeaks_unsorted/peak_sets/0_peaks.bed \
 ../results/Seurat/callpeaks_unsorted/peak_sets/1_peaks.bed \
 ../results/Seurat/callpeaks_unsorted/peak_sets/2_peaks.bed \
 ../results/Seurat/callpeaks_unsorted/peak_sets/3_peaks.bed \
 ../results/Seurat/callpeaks_unsorted/peak_sets/4_peaks.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "mESC ATAC-Seq" --skipZeros --missingDataAsZero
 
plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise_ATAC.mat.gz \
 -out "../results/deeptools/unsorted_clusterwise_ATAC.pdf" \
 --refPointLabel "G4" \
 --whatToShow "heatmap and colorbar" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "0" "1" "2" "3" "4" \
 --yAxisLabel "" \
 --xAxisLabel "" \


plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise_ATAC.mat.gz \
 -out "../results/deeptools/unsorted_clusterwise_ATAC.png" \
 --refPointLabel "G4" \
 --whatToShow "heatmap and colorbar" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "0" "1" "2" "3" "4" \
 --yAxisLabel "" \
 --xAxisLabel "" \
