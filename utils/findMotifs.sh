#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH -M snowy
#SBATCH -J motif_analysis

module load bioinfo-tools
module load HOMER

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

# make output dir
mkdir -p ../results/Seurat/callpeaks_GFPsorted/findMotif_output/

# HOMER motif analysis
findMotifs.pl ../results/Seurat/callpeaks_GFPsorted/motif_analysis_input.txt \
 mouse ../results/Seurat/callpeaks_GFPsorted/findMotif_output/ \
 -start -500 -end 100 -len 8,10,12