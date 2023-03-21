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

computeMatrix reference-point -o ../results/deeptools/matrix_mESC_MEF_ol_ATAC.mat.gz \
 -S ../data/bw/ATAC-H33WT.mm10.bw \
 -R ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/MACS2_cl0_only.bed \
 ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/MACS2_intersection.bed \
 ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/MACS2_cl1_only.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "mESC ATAC-Seq" \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_mESC_MEF_ol_ATAC.mat.gz \
 -out "../results/deeptools/mESC_MEF_ol_ATAC_Marek_MACS2.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Blues" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_mESC_MEF_ol_ATAC.mat.gz \
 -out "../results/deeptools/mESC_MEF_ol_ATAC_Marek_MACS2.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Blues" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \
