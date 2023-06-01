#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -M rackham
#SBATCH -J seurat

# load modules
module load bioinfo-tools
module load R_packages

# change to dir where the analysis codes are
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/final_codes

# run analysis on unsorted brain dataset
#Rscript seurat_workflow.R -c "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/CellRanger/unsorted/" \
# -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.8" \
# -r 0.8

#Rscript seurat_workflow.R -c "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/CellRanger/unsorted/" \
# -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/unsorted_brain/res0.1" \
# -r 0.1

# run analysis on sorted brain dataset
Rscript seurat_workflow.R -c "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/CellRanger/GFP_sorted/" \
 -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/sorted_brain" \
 -r 0.8
