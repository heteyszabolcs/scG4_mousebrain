#!/bin/bash -l 
#SBATCH -A naiss2023-22-84 
#SBATCH -p core 
#SBATCH -n 8 
#SBATCH -t 20:00:00 
#SBATCH -M rackham
#SBATCH -J macs2_deeptools

# load modules
module load bioinfo-tools
module load deepTools
module load R_packages
module load MACS/2.2.6
module load samtools/1.17

# change to dir containing the R script
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/final_codes

# run analysis on unsorted dataset
Rscript bam_process.R -s "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.8/outputs/Seurat_object.Rds" \
 -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.8"