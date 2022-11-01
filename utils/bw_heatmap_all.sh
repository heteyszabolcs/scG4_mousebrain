#!/bin/bash -l
#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M snowy
#SBATCH -J deeptools

module load bioinfo-tools
module load deepTools
module load R_packages

cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

for i in H3K27me3 H3K4me3 H3K27ac H3K36me3
do
	for j in Astrocytes mOL OPC OEC VLMC
	do
		for k in 0 1 2 3 4 
		do
			  Rscript bw_heatmap_all.R -bw "../data/GSE157637/${i}_${j}.bw" \
			 -b "../data/bed/${k}_peaks_lanceotron.bed" \
			 -r "G4_cluster_${k}" -l "${i}_${j}" \
			 -o "../results/deeptools/cluster${k}_${i}_${j}_lanceotron.png" \
			 -ymax 10 -yzmax 10
		done			
	done
done

for i in H3K27me3 H3K4me3 H3K27ac H3K36me3
do
	for j in Astrocytes mOL OPC OEC VLMC
	do
		for k in 0 1 2 3 4 
		do
			  Rscript bw_heatmap_all.R -bw "../data/GSE157637/${i}_${j}.bw" \
			 -b "../data/bed/${k}_peaks_robust_peaks.bed" \
			 -r "G4_cluster_${k}" -l "${i}_${j}" \
			 -o "../results/deeptools/cluster${k}_${i}_${j}_seurat_macs2.png" \
			 -ymax 10 -yzmax 10
		done			
	done
done


for cluster in 0 1 2 3 4 
do
	Rscript bw_heatmap_all.R -bw "/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/data/GSE163484/rep1/outs/GSE163484_rep1.bw" \
	-b "../data/bed/${cluster}_peaks_lanceotron.bed" \
	-r "G4_cluster_${cluster}" -l "RNA_Seq_rep1" \
	-o "../results/deeptools/cluster${cluster}_RNA_Seq_lanceotron.png" \
	-ymax 50 -yzmax 50
done			


