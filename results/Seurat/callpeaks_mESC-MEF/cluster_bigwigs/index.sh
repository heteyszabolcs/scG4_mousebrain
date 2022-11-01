#!/bin/bash -l 
#SBATCH -A snic2020-15-9 
#SBATCH -p core 
#SBATCH -n 4 
#SBATCH -t 16:00:00 
#SBATCH -M snowy
#SBATCH -J samtools

module load bioinfo-tools
module load samtools

cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs

samtools index 0.bam
