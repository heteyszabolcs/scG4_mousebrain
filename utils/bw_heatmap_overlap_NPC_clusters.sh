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

computeMatrix reference-point -o ../results/deeptools/matrix_mESC_NPC_ol_clusters.mat.gz -S \
../results/Seurat/callpeaks_unsorted/cluster_bigwigs/0.bam.bw \
../results/Seurat/callpeaks_unsorted/cluster_bigwigs/1.bam.bw \
../results/Seurat/callpeaks_unsorted/cluster_bigwigs/2.bam.bw \
../results/Seurat/callpeaks_unsorted/cluster_bigwigs/3.bam.bw \
../results/Seurat/callpeaks_unsorted/cluster_bigwigs/4.bam.bw \
../results/Seurat/callpeaks_unsorted/cluster_bigwigs/5.bam.bw \
../results/Seurat/callpeaks_unsorted/cluster_bigwigs/6.bam.bw \
-R ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/only_mESC_bulk.bed \
../results/GenomicRanges/unsorted_outputs/bulk_overlaps/intersection.bed \
../results/GenomicRanges/unsorted_outputs/bulk_overlaps/only_NPC_bulk.bed \
-b 3000 -a 3000 --samplesLabel "cluster 0" "cluster 1" "cluster 2" "cluster 3" "cluster 4" "cluster 5" "cluster 6" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_ol_clusters.mat.gz \
 -out "../results/deeptools/mESC_NPC_ol-clusters.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "mESC" "both" "NPC" \
 --whatToShow "heatmap only" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --legendLocation "none"

plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_ol_clusters.mat.gz \
 -out "../results/deeptools/mESC_NPC_ol-clusters.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "mESC" "both" "NPC" \
 --whatToShow "heatmap only" \
 --yAxisLabel "" \
 --xAxisLabel "" \
 --legendLocation "none"