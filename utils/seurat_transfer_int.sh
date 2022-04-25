#!/bin/bash -l
#SBATCH -M snowy

module load bioinfo-tools
module load R_packages

# working directory
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/

Rscript /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/seurat_transfer_int.R
