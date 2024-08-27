#!/bin/bash -l
#SBATCH -A naiss2024-22-108
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools
module load R_packages

cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

# computeMatrix reference-point -o ../results/deeptools/matrix_unsorted_PQS.mat.gz -S \
 # ../data/bw/all_PQS_mm10_binary.bw \
 # ../data/bw/ATAC-H33WT.mm10.bw \
 # -R ../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak \
 # ../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak \
 # --referencePoint center \
 # -b 3000 -a 3000 --samplesLabel "PQS" "mESC ATAC-Seq" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_unsorted_PQS.mat.gz \
 -out "../results/deeptools/unsorted_clusters_PQS.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --yMin 0 0 \
 --yMax 1 250 \
 -z "cluster 0" "cluster 1" \
 --whatToShow "heatmap only" \
 --colorList "white, black" "white, #addd8e" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --legendLocation "none"

plotHeatmap -m ../results/deeptools/matrix_unsorted_PQS.mat.gz \
 -out "../results/deeptools/unsorted_clusters_PQS.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --yMin 0 0 \
 --yMax 1 250 \
 --sortUsing "max" \
 -z "cluster 0" "cluster 1" \
 --whatToShow "heatmap only" \
 --colorList "white, black" "white, #addd8e" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --legendLocation "none"