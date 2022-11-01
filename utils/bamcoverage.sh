#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 16:00:00
#SBATCH -M snowy
#SBATCH -J bamcoverage

# go to workdir
cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

# modules
module load bioinfo-tools
module load deepTools

bamCoverage --bam $1 -o $2 \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2652783500