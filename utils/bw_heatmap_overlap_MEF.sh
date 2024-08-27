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

computeMatrix reference-point -o ../results/deeptools/matrix_mESC_MEF_overlap.mat.gz \
 -S ../data/Marek_data/mESC-MEF/G4_H33WT_SL_CnT_R1.mm10.bw \
 ../data/Marek_data/mESC-MEF/G4_MEF_CnT_R1.mm10.bw \
 ../results/Seurat/final/mesc_mef/cluster_spec_bigwigs/1.bam_RPGC.bigwig \
 ../results/Seurat/final/mesc_mef/cluster_spec_bigwigs/0.bam_RPGC.bigwig \
 ../data/bw/ATAC-H33WT.mm10.bw \
 -R ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/MACS2_cl0_only.bed \
 ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/MACS2_intersection.bed \
 ../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/MACS2_cl1_only.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "bulk mESC" "bulk MEF" "mESC cluster 1" "MEF cluster 0" "mESC ATAC-Seq" \
 --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_mESC_MEF_overlap.mat.gz \
 -out "../results/deeptools/mESC_MEF_overlap_Marek_MACS2.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorList "white, #fee0d2" "white, #deebf7" "white, #fc9272" "white, #9ecae1" "white, #bdbdbd"  \
 --yMin 0 0 0 0 0 \
 --yMax 60 60 100 100 100 \
 --zMin 0 0 0 0 0 \
 --zMax 60 60 100 100 100 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \

plotHeatmap -m ../results/deeptools/matrix_mESC_MEF_overlap.mat.gz \
 -out "../results/deeptools/mESC_MEF_overlap_Marek_MACS2.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorList "white, #fee0d2" "white, #deebf7" "white, #fc9272" "white, #9ecae1" "white, #bdbdbd"  \
 --yMin 0 0 0 0 0 \
 --yMax 60 60 100 100 100 \
 --zMin 0 0 0 0 0 \
 --zMax 60 60 100 100 100 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \
