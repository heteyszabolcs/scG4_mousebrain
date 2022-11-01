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

computeMatrix reference-point -o \
 ../results/deeptools/matrix_mESC_MEF_epigmarkers.mat.gz \
 -S ../data/bw/bulk_G4_CnT_mESC_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
 ../data/bw/bulk_G4_CnT_MEF_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
 ../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/1_res0.1.bigwig \
 ../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/0_res0.1.bigwig \
 -R ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/only_mESC_bulk.bed \
 ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/intersection.bed \
 ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/only_MEF_bulk.bed \
 -b 3000 -a 3000 \
 --samplesLabel "H3K4me3" "H3K4me1" "H3K27me3" "H3K27ac" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_mESC_MEF_epigmarkers.mat.gz \
 -out "../results/deeptools/mESC_MEF_ol_epigmarkers.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Blues" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "mESC" "both" "MEF" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_mESC_MEF_epigmarkers.mat.gz \
 -out "../results/deeptools/mESC_MEF_ol_epigmarkers.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Blues" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "mESC" "both" "MEF" \
 --yAxisLabel "" \
 --xAxisLabel "" \