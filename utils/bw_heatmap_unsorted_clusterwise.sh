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

# computeMatrix reference-point -o ../results/deeptools/matrix_unsorted_clusterwise.mat.gz -S \
 # ../data/bw/GSM4303808_WT-H3K4me3-1_mm10.bigwig \
  # ../data/bw/GSM4303796_WT-H3K4me1-1_mm10.bigwig \
 # ../data/bw/GSM4377779_WT-H3K27me3-1_mm10.bigwig \
 # ../data/bw/GSM4205678_WT-H3K27ac-1.10M_depth_mm10.bigwig \
 # ../data/bw/ATAC-H33WT.mm10.bw \
 # -R ../results/Seurat/callpeaks_unsorted/peak_sets/0_peaks.bed \
 # ../results/Seurat/callpeaks_unsorted/peak_sets/1_peaks.bed \
 # ../results/Seurat/callpeaks_unsorted/peak_sets/2_peaks.bed \
 # ../results/Seurat/callpeaks_unsorted/peak_sets/3_peaks.bed \
 # ../results/Seurat/callpeaks_unsorted/peak_sets/4_peaks.bed \
 # ../results/Seurat/callpeaks_unsorted/peak_sets/5_peaks.bed \
 # ../results/Seurat/callpeaks_unsorted/peak_sets/6_peaks.bed \
 # -b 3000 -a 3000 --samplesLabel "H3K4me3" "H3K4me1" "H3K27me3" "H3K27ac" "ATAC-Seq" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise.mat.gz \
 -out "../results/deeptools/unsorted_clusterwise.pdf" \
 --refPointLabel "G4" \
 --whatToShow "heatmap and colorbar" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "0" "1" "2" "3" "4" "5" "6" \
 --yAxisLabel "" \
 --xAxisLabel "" \


plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise.mat.gz \
 -out "../results/deeptools/unsorted_clusterwise.png" \
 --refPointLabel "G4" \
 --whatToShow "heatmap and colorbar" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "0" "1" "2" "3" "4" "5" "6" \
 --yAxisLabel "" \
 --xAxisLabel "" \
