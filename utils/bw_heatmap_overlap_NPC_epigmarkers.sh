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

computeMatrix reference-point -o ../results/deeptools/matrix_mESC_NPC_epigmarkers.mat.gz -S \
 ../data/bw/GSM4303808_WT-H3K4me3-1_mm10.bigwig \
  ../data/bw/GSM4303796_WT-H3K4me1-1_mm10.bigwig \
 ../data/bw/GSM4377779_WT-H3K27me3-1_mm10.bigwig \
 ../data/bw/GSM4205678_WT-H3K27ac-1.10M_depth_mm10.bigwig \
 -R ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_cl0_only.bed \
 ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_intersection.bed \
 ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_cl1_only.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "H3K4me3" "H3K4me1" "H3K27me3" "H3K27ac" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_epigmarkers.mat.gz \
 -out "../results/deeptools/mESC_NPC_ol_epigmarkers.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 15 \
 --zMax 15 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \


plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_epigmarkers.mat.gz \
 -out "../results/deeptools/mESC_NPC_ol_epigmarkers.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Greens" \
 --yMin 0 \
 --yMax 15 \
 --zMax 15 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \
