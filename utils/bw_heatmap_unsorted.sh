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

# computeMatrix reference-point -o ../results/deeptools/matrix_unsorted.mat.gz -S \
 # ../data/CellRanger/unsorted/possorted_RPGC.bigwig \
 # -R ../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak \
 # ../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak \
 # --referencePoint center \
 # -b 3000 -a 3000 --samplesLabel "G4" --skipZeros --missingDataAsZero

# plotHeatmap -m ../results/deeptools/matrix_unsorted.mat.gz \
 # -out "../results/deeptools/unsorted_clusters.pdf" \
 # --refPointLabel "G4" \
 # --heatmapHeight 14 \
 # --yMin 0 \
 # --yMax 15 \
 # -z "cluster 0" "cluster 1" \
 # --colorMap "Blues" \
 # --yAxisLabel "" \
 # --xAxisLabel "" \
 # --legendLocation "none"
 #--whatToShow "heatmap only" \

# plotHeatmap -m ../results/deeptools/matrix_unsorted.mat.gz \
 # -out "../results/deeptools/unsorted_clusters.png" \
 # --refPointLabel "G4" \
 # --heatmapHeight 14 \
 # --yMin 0 \
 # --yMax 15 \
 # --sortUsing "max" \
 # -z "cluster 0" "cluster 1" \
 # --colorMap "Blues" \
 # --yAxisLabel "" \
 # --xAxisLabel "" \
 # --legendLocation "none"
 #--whatToShow "heatmap only" \
 
 plotProfile -m ../results/deeptools/matrix_unsorted.mat.gz \
              -out ../results/deeptools/unsorted_clusters_profilePlot.pdf \
              --plotTitle "" \
			  --yMin 0 \
			  --yMax 40