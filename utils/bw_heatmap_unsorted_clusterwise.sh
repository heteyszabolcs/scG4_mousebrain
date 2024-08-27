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

computeMatrix reference-point -o ../results/deeptools/matrix_unsorted_clusterwise_max.mat.gz -S \
 ../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/0.bam_RPGC.bigwig \
 ../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/1.bam_RPGC.bigwig \
 ../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/2.bam_RPGC.bigwig \
 ../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/3.bam_RPGC.bigwig \
 ../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/4.bam_RPGC.bigwig \
 ../data/bw/all_PQS_mm10_binary.bw \
 -R /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl0_peaks.bed \
 /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl1_peaks.bed \
 /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl2_peaks.bed \
 /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl3_peaks.bed \
 /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl4_peaks.bed \
 --referencePoint center \
 -b 3000 -a 3000 --samplesLabel "cluster 0" "cluster 1" "cluster 2" "cluster 3" "cluster 4" "PQS" \
 --averageTypeBins "max" \
 --skipZeros --missingDataAsZero --outFileNameMatrix ../results/deeptools/matrix_unsorted_clusterwise_max.tab
 
plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise_max.mat.gz \
 -out "../results/deeptools/unsorted_clusterwise_max.pdf" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorList "white, #3182bd" "white, #3182bd" "white, #3182bd" "white, #3182bd" "white, #3182bd" "white, black" \
 --yMin 0 0 0 0 0 0 \
 --yMax 40 40 40 40 40 1 \
 --zMin 0 0 0 0 0 0 \
 --zMax 40 40 40 40 40 1 \
 --yAxisLabel "" \
 --xAxisLabel "" \
 -z "cluster 0" "cluster 1" "cluster 2" "cluster 3" "cluster 4"

plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise_max.mat.gz \
 -out "../results/deeptools/unsorted_clusterwise_max.png" \
 --refPointLabel "G4" \
 --heatmapHeight 14 \
 --colorList "white, #3182bd" "white, #3182bd" "white, #3182bd" "white, #3182bd" "white, #3182bd" "white, black" \
 --yMin 0 0 0 0 0 0 \
 --yMax 40 40 40 40 40 1 \
 --zMin 0 0 0 0 0 0 \
 --zMax 40 40 40 40 40 1 \
 --yAxisLabel "" \
 --xAxisLabel "" \
 -z "cluster 0" "cluster 1" "cluster 2" "cluster 3" "cluster 4"

# computeMatrix reference-point -o ../results/deeptools/matrix_unsorted_clusterwise_pqs.mat.gz -S \
  # ../data/bw/all_PQS_mm10_binary.bw \
 # -R /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl0_peaks.bed \
 # /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl1_peaks.bed \
 # /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl2_peaks.bed \
 # /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl3_peaks.bed \
 # /crex/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/GenomicRanges/unsorted_outputs/unsorted-unique_cl4_peaks.bed \
 # --referencePoint center \
 # -b 3000 -a 3000 --samplesLabel "PQS" \
 # --skipZeros --missingDataAsZero
 
 # plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise_pqs.mat.gz \
 # -out "../results/deeptools/unsorted_clusterwise_pqs.png" \
 # --refPointLabel "G4" \
 # --heatmapHeight 14 \
 # --whatToShow "heatmap only" \
 # --colorMap "binary" \
 # --yMin 0 \
 # --yMax 0.01 \
 # --zMin 0 \
 # --zMax 0.01 \
 # --yAxisLabel "" \
 # --xAxisLabel "" \
 # -z "cluster 0" "cluster 1" "cluster 2" "cluster 3" "cluster 4"
 
  # plotHeatmap -m ../results/deeptools/matrix_unsorted_clusterwise_pqs.mat.gz \
 # -out "../results/deeptools/unsorted_clusterwise_pqs.pdf" \
 # --refPointLabel "G4" \
 # --heatmapHeight 14 \
 # --whatToShow "heatmap only" \
 # --colorMap "binary" \
 # --yMin 0 \
 # --yMax 0.01 \
 # --zMin 0 \
 # --zMax 0.01 \
 # --yAxisLabel "" \
 # --xAxisLabel "" \
 # -z "cluster 0" "cluster 1" "cluster 2" "cluster 3" "cluster 4"