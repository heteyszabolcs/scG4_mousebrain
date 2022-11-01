#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 16:00:00
#SBATCH -M snowy
#SBATCH -J cellranger

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

# modules
module load bioinfo-tools
module load cellranger/6.1.2

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/GSE163484

# call cellranger count
cellranger count --fastqs=/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/GSE163484 \
--id=rep2 \
--sample=rep2 \
--transcriptome=/sw/data/Chromium/cellranger-data/2020-A/refdata-gex-mm10-2020-A
