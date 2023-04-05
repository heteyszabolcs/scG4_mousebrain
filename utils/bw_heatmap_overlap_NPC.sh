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

# computeMatrix reference-point -o ../results/deeptools/matrix_mESC_NPC_overlap.mat.gz -S \
 # ../data/bw/bulk_G4_CnT_mESC_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
 # ../data/bw/bulk_G4_CnT_neuron_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
 # ../data/bw/bulk_G4_CnT_NPC_rep1_R1.mLb.clN_mm10_RPGC.bigwig \
 # ../results/Seurat/callpeaks_unsorted/cluster_bigwigs/1_res0.1_brain_cells.bw \
 # ../results/Seurat/callpeaks_unsorted/cluster_bigwigs/0_res0.1_oligo.bw \
 # -R ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_cl0_only.bed \
 # ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_intersection.bed \
 # ../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_cl1_only.bed \
 # --referencePoint center \
 # -b 3000 -a 3000 --samplesLabel "bulk mESC" "bulk neuron" "bulk NPC" "brain cl 1" "oligodendr cl 0" --skipZeros --missingDataAsZero

plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_overlap.mat.gz \
 -out "../results/deeptools/mESC_NPC_overlap_MACS2.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Blues" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \


plotHeatmap -m ../results/deeptools/matrix_mESC_NPC_overlap.mat.gz \
 -out "../results/deeptools/mESC_NPC_overlap_MACS2.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorMap "Blues" \
 --yMin 0 \
 --yMax 100 \
 --zMax 100 \
 -z "cluster 0" "both" "cluster 1" \
 --yAxisLabel "" \
 --xAxisLabel "" \
