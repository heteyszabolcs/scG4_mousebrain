#!/bin/bash -l 
#SBATCH -A naiss2023-22-84 
#SBATCH -p core 
#SBATCH -n 8 
#SBATCH -t 20:00:00 
#SBATCH -M rackham
#SBATCH -J integration

# load modules
module load bioinfo-tools
module load R_packages

# change to dir containing the R script
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/final_codes

# run integration on GFP sorted dataset
Rscript seurat_integration.R -s "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/sorted_brain/outputs/Seurat_object.Rds" \
 -w "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/final/sorted_brain/integration"