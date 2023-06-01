#!/bin/bash -l 
#SBATCH -A naiss2023-22-84 
#SBATCH -p core 
#SBATCH -n 8 
#SBATCH -t 20:00:00 
#SBATCH -M rackham
#SBATCH -J sinto

# load modules
module load bioinfo-tools
module load R_packages

# change to utils dir
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/final_codes

# source pyenv environment
source /home/szabolcs/.pyenv/versions/3.8.1/envs/sinto/bin/activate

# run analysis on unsorted dataset
#Rscript sinto.R -s "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.8/outputs/Seurat_object.Rds" \
# -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.8" \
# -b "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/CellRanger/unsorted/possorted_bam.bam"

#Rscript sinto.R -s "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.1/outputs/Seurat_object.Rds" \
# -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.1" \
# -b "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/CellRanger/unsorted/possorted_bam.bam"

# run analysis on GFP sorted dataset
Rscript sinto.R -s "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/sorted_brain/outputs/Seurat_object.Rds" \
 -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/sorted_brain" \
 -b "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/CellRanger/GFP_sorted/possorted_bam.bam"
