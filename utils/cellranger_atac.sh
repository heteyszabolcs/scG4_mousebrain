#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 48:00:00
#SBATCH -M snowy

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/GEOPIPE/scCutnTag/Marek_scCnT/gfp_negative_k4me3

# modules
module load bioinfo-tools
module load cellranger-ATAC/2.0.0

# call cellranger-atac count
cellranger-atac count --id=Marek_gfp_neg_k4me3 --fastqs=/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/Marek_data/histone_scCnT/gfp_negative_k4me3/fastq --reference=/sw/data/Chromium/cellranger-ATAC-data/2.0.0/rackham/refdata-cellranger-arc-mm10-2020-A-2.0.0
