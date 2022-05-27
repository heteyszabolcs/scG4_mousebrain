#!/bin/bash -l
#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M snowy

# working directory
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils/

module load bioinfo-tools
module load deepTools

sinto filterbarcodes -b ../data/CellRanger/GFP_sorted/possorted_bam.bam -c $1 --outdir ../results/Seurat/callpeaks_GFPsorted/ 

bamCoverage -b ../results/Seurat/callpeaks_GFPsorted/1.bam -o ../results/Seurat/callpeaks_GFPsorted/$1.bw

rm 1.bam 
