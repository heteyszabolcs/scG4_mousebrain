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

computeMatrix reference-point -o ../results/deeptools/matrix_unsorted_ATAC.mat.gz \
 -S ../data/bw/ATAC-H33WT.mm10.bw \
 -R ../results/Seurat/callpeaks_unsorted/peak_sets/0_peaks_res0.1.bed \
	../results/Seurat/callpeaks_unsorted/peak_sets/1_peaks_res0.1.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "mESC ATAC-Seq" \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_unsorted_ATAC.mat.gz \
 -out "../results/deeptools/unsorted_ATAC.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 250 \
 --zMax 250 \
 -z "cluster 0" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_unsorted_ATAC.mat.gz \
 -out "../results/deeptools/unsorted_ATAC.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 250 \
 --zMax 250 \
 -z "cluster 0" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \
