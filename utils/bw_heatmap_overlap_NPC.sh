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

#computeMatrix reference-point -o ../results/deeptools/matrix_mESC_NPC_overlap.mat.gz -S \
# ../data/bw/bulk_G4_CnT_mESC_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
# ../data/bw/bulk_G4_CnT_neuron_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
# ../data/bw/bulk_G4_CnT_NPC_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
# ../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/1.bam_RPGC.bigwig \
# ../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/0.bam_RPGC.bigwig \
# -R ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/unique_cluster0-median_filtered-GR.bed \
# ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_intersection.bed \
# ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/unique_cluster1-median_filtered-GR.bed \
# --referencePoint center \
# -b 3000 -a 3000 --samplesLabel "bulk mESC" "bulk neuron" "bulk NPC" "cluster 1" "cluster 0" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_overlap.mat.gz \
 -out "../results/deeptools/mESC_NPC_overlap_MACS2.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorList "white, #addd8e" "white, #addd8e" "white, #addd8e" "white, #addd8e" "white, #addd8e"  \
 --yMin 0 0 0 0 0 \
 --yMax 40 40 40 80 40 \
 --zMin 0 0 0 0 0 \
 --zMax 40 40 40 80 40 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \


plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_overlap.mat.gz \
 -out "../results/deeptools/mESC_NPC_overlap_MACS2.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorList "white, #addd8e" "white, #addd8e" "white, #addd8e" "white, #addd8e" "white, #addd8e"  \
 --yMin 0 0 0 0 0 \
 --yMax 40 40 40 80 40 \
 --zMin 0 0 0 0 0 \
 --zMax 40 40 40 80 40 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \
