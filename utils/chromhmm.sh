#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J wigglescout

module load bioinfo-tools
module load R_packages

# working directory
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/

Rscript /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/chromhmm.R
