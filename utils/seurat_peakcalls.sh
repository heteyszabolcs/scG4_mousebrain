#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J seurat_peakcall

module load bioinfo-tools
module load R_packages
module load MACS/2.2.6

# working directory
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/

Rscript /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/seurat_peakcalls.R
