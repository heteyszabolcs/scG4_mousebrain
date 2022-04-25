#!/bin/bash -l
#SBATCH -M snowy

module load bioinfo-tools
module load deepTools
module load R_packages

cd /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/utils

# H3K27Ac ChIP-Seq
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ESC_H33WT_H3K27ac.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3_H3K27ac" -o "../results/deeptools/H3.3_H3K27ac_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ESC_H33KO_H3K27ac.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3KO_H3K27ac" -o "../results/deeptools/H3.3KO_H3K27ac_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/bw/Chronis2017_Klf4.bs5.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Klf4" -o "../results/deeptools/Klf4_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/bw/Chronis2017_Oct4.bs5.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Oct4" -o "../results/deeptools/Oct4_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/bw/Chronis2017_Sox2.bs5.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Sox2" -o "../results/deeptools/Sox2_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/bw/King2017_BRG1fl_BRG1_1.bs5.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Brg1fl_Brg1" -o "../results/deeptools/Brg1fl_Brg1_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/bw/King2017_BRG1KO_BRG1_1.bs5.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Brg1KO" -o "../results/deeptools/Brg1KO_ltr.png" \
# -ymax 10 -yzmax 10

# ATAC-Seq
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ATAC-H33WT.mm9.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3" -o "../results/deeptools/H3.3_ltr.png" \
# -ymax 150 -yzmax 150
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ATAC-H33KO.mm9.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3KO" -o "../results/deeptools/H3.3KO_ltr.png" \
# -ymax 150 -yzmax 150
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ATAC-H33KO-H32rescue-802.mm9.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3KO-H3.2rescue" -o "../results/deeptools/H3.3KO-H3.2rescue_ltr.png" \
# -ymax 150 -yzmax 150
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ATAC-H33KO-H33mut-rescue-803.mm9.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3KO-H3.3mut_rescue" -o "../results/deeptools/HH3.3KO-H3.3mut_rescue_ltr.png" \
# -ymax 150 -yzmax 150
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ATAC-H33KO-H33rescue-801.mm9.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3KO-H3.3_rescue" -o "../results/deeptools/HH3.3KO-H3.3_rescue_ltr.png" \
# -ymax 150 -yzmax 150
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ATAC-H33KO-smarcad1KD.mm9.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3KO-Smarcad1KD" -o "../results/deeptools/H3.3KO-Smarcad1KD_ltr.png" \
# -ymax 150 -yzmax 150
# Rscript bw_heatmap_ltr.R -bw "../data/bw/ATAC-H33WT-smarcad1KD.mm9.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3-Smarcad1KD" -o "../results/deeptools/H3.3-Smarcad1KD_ltr.png" \
# -ymax 150 -yzmax 150Rscript bw_heatmap_ltr.R -bw "../data/bw/ESC_H33WT_H3K27ac.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3.3_H3K27ac" -o "../results/deeptools/H3.3_H3K27ac_ltr.png" \
#-ymax 10 -yzmax 10

# Bartosovic H3K27ac scCut&Tag - H33_dependent_G4_JL.v2
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27ac_Astrocytes.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27ac-Astrocytes" -o "../results/deeptools/Bartosovic_H3K27ac-Astrocytes" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27ac_mOL.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27ac-mOL" -o "../results/deeptools/Bartosovic_H3K27ac-mOL" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27ac_VLMC.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27ac-VLMC" -o "../results/deeptools/Bartosovic_H3K27ac-VLMC" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27ac_OEC.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27ac-OEC" -o "../results/deeptools/Bartosovic_H3K27ac-OEC" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27ac_OPC.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27ac-OPC" -o "../results/deeptools/Bartosovic_H3K27ac-OPC" \
# -ymax 10 -yzmax 10

# Bartosovic H3K27ac scCut&Tag - RepMasker_lt200bp.LTRIS2.bed
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27ac_Astrocytes.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27ac-Astrocytes" -o "../results/deeptools/Bartosovic_H3K27ac-Astrocytes_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27ac_mOL.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27ac-mOL" -o "../results/deeptools/Bartosovic_H3K27ac-mOL_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27ac_VLMC.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27ac-VLMC" -o "../results/deeptools/Bartosovic_H3K27ac-VLMC_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27ac_OEC.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27ac-OEC" -o "../results/deeptools/Bartosovic_H3K27ac-OEC_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27ac_OPC.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27ac-OPC" -o "../results/deeptools/Bartosovic_H3K27ac-OPC_ltr.png" \
# -ymax 10 -yzmax 10

# Bartosovic Rad21 - H33_dependent_G4_JL.v2
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/Rad21_OEC.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "Rad21-OEC" -o "../results/deeptools/Bartosovic_Rad21-OEC.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/Rad21_Astrocytes.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "Rad21-Astrocytes" -o "../results/deeptools/Bartosovic_Rad21-Astrocytes.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/Rad21_mOL.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "Rad21-mOL" -o "../results/deeptools/Bartosovic_Rad21-mOL.png" \
# -ymax 10 -yzmax 10

# # Bartosovic Rad21 - RepMasker_lt200bp.LTRIS2.bed
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/Rad21_OEC.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Rad21-OEC" -o "../results/deeptools/Bartosovic_Rad21-OEC_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/Rad21_Astrocytes.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Rad21-Astrocytes" -o "../results/deeptools/Bartosovic_Rad21-Astrocytes_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/Rad21_mOL.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "Rad21-mOL" -o "../results/deeptools/Bartosovic_Rad21-mOL_ltr.png" \
# -ymax 10 -yzmax 10

# # Bartosovic H3K27me3 scCut&Tag - H33_dependent_G4_JL.v2
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27me3_Astrocytes.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27me3-Astrocytes" -o "../results/deeptools/Bartosovic_H3K27me3-Astrocytes" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27me3_mOL.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27me3-mOL" -o "../results/deeptools/Bartosovic_H3K27me3-mOL" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27me3_VLMC.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27me3-VLMC" -o "../results/deeptools/Bartosovic_H3K27me3-VLMC" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27me3_OEC.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27me3-OEC" -o "../results/deeptools/Bartosovic_H3K27me3-OEC" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_H33_g4.R -bw "../data/GSE157637/H3K27me3_OPC.bw" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "H3K27me3-OPC" -o "../results/deeptools/Bartosovic_H3K27me3-OPC" \
# -ymax 10 -yzmax 10

# # Pleth H3K27ac pre-iPSC ChIP-Seq - H33_dependent_G4_JL.v2
# Rscript bw_heatmap_H33_g4.R -bw "../data/bw/Plath_et_al-pre_IPS.K27AC.gaussian_10.wig" -b "../data/bed/H33_dependent_G4_JL.v2.bed" -l "pre-iPSC_H3K27ac" -o "../results/deeptools/pre-iPSC_H3K27ac.png" \
# -ymax 10 -yzmax 10

# # Bartosovic H3K27me3 scCut&Tag - RepMasker_lt200bp.LTRIS2.bed
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27me3_Astrocytes.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27me3-Astrocytes" -o "../results/deeptools/Bartosovic_H3K27me3-Astrocytes_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27me3_mOL.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27me3-mOL" -o "../results/deeptools/Bartosovic_H3K27me3-mOL_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27me3_VLMC.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27me3-VLMC" -o "../results/deeptools/Bartosovic_H3K27me3-VLMC_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27me3_OEC.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27me3-OEC" -o "../results/deeptools/Bartosovic_H3K27me3-OEC_ltr.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_ltr.R -bw "../data/GSE157637/H3K27me3_OPC.bw" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "H3K27me3-OPC" -o "../results/deeptools/Bartosovic_H3K27me3-OPC_ltr.png" \
# -ymax 10 -yzmax 10

# # Pleth H3K27ac pre-iPSC ChIP-Seq - RepMasker_lt200bp.LTRIS2.bed
# Rscript bw_heatmap_ltr.R -bw "../data/bw/Plath_et_al-pre_IPS.K27AC.gaussian_10.wig" -b "../data/bed/RepMasker_lt200bp.LTRIS2.bed" -l "pre-iPSC_H3K27ac" -o "../results/deeptools/pre-iPSC_H3K27ac_ltr.png" \
 # -ymax 10 -yzmax 10

# Bartosovic H3K27ac scCut&Tag - ESC_Enhancer_CruzMolina.active.bed (active enhancers, CruzMolina et al.)
# Rscript bw_heatmap_activeenh_CruzMolina.R -bw "../data/GSE157637/H3K27ac_Astrocytes.bw" -b "../data/bed/ESC_Enhancer_CruzMolina.active.bed" -l "H3K27ac-Astrocytes" -o "../results/deeptools/Bartosovic_H3K27ac-Astrocytes_actenh.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_CruzMolina.R -bw "../data/GSE157637/H3K27ac_mOL.bw" -b "../data/bed/ESC_Enhancer_CruzMolina.active.bed" -l "H3K27ac-mOL" -o "../results/deeptools/Bartosovic_H3K27ac-mOL_actenh.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_CruzMolina.R -bw "../data/GSE157637/H3K27ac_VLMC.bw" -b "../data/bed/ESC_Enhancer_CruzMolina.active.bed" -l "H3K27ac-VLMC" -o "../results/deeptools/Bartosovic_H3K27ac-VLMC_actenh.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_CruzMolina.R -bw "../data/GSE157637/H3K27ac_OEC.bw" -b "../data/bed/ESC_Enhancer_CruzMolina.active.bed" -l "H3K27ac-OEC" -o "../results/deeptools/Bartosovic_H3K27ac-OEC_actenh.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_CruzMolina.R -bw "../data/GSE157637/H3K27ac_OPC.bw" -b "../data/bed/ESC_Enhancer_CruzMolina.active.bed" -l "H3K27ac-OPC" -o "../results/deeptools/Bartosovic_H3K27ac-OPC_actenh.png" \
# -ymax 10 -yzmax 10

## enhancer subsets
# # G4 scCut&Tag - ESC_Enhancer_CruzMolina.active.bed (active enhancers, CruzMolina et al.)
# Rscript bw_heatmap_activeenh_CruzMolina.R -bw "../data/bw/G4_mES_H33WT_NC_R1_S33_L001.mm9.bw" -b "../data/bed/ESC_Enhancer_CruzMolina.active.bed" -l "G4-H3.3WT" -o "../results/deeptools/G4-H3.3WT_actenh.png" \
# -ymax 50 -yzmax 50
# Rscript bw_heatmap_activeenh_CruzMolina.R -bw "../data/bw/G4_mES_H33KO_NC_R1_S37_L001.mm9.bw" -b "../data/bed/ESC_Enhancer_CruzMolina.active.bed" -l "G4-H3.3KO" -o "../results/deeptools/G4-H3.3KO_actenh.png" \
# -ymax 50 -yzmax 50

# Bartosovic H3K27ac scCut&Tag - ESC_Enhancer_CruzMolina.active.bed (active enhancers, Glaser et al.)
# Rscript bw_heatmap_activeenh_Glaser.R -bw "../data/GSE157637/H3K27ac_Astrocytes.bw" -b "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed" -l "H3K27ac-Astrocytes" -o "../results/deeptools/Bartosovic_H3K27ac-Astrocytes_actenh_Glaser.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_Glaser.R -bw "../data/GSE157637/H3K27ac_mOL.bw" -b "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed" -l "H3K27ac-mOL" -o "../results/deeptools/Bartosovic_H3K27ac-mOL_actenh_Glaser.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_Glaser.R -bw "../data/GSE157637/H3K27ac_VLMC.bw" -b "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed" -l "H3K27ac-VLMC" -o "../results/deeptools/Bartosovic_H3K27ac-VLMC_actenh_Glaser.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_Glaser.R -bw "../data/GSE157637/H3K27ac_OEC.bw" -b "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed" -l "H3K27ac-OEC" -o "../results/deeptools/Bartosovic_H3K27ac-OEC_actenh_Glaser.png" \
# -ymax 10 -yzmax 10
# Rscript bw_heatmap_activeenh_Glaser.R -bw "../data/GSE157637/H3K27ac_OPC.bw" -b "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed" -l "H3K27ac-OPC" -o "../results/deeptools/Bartosovic_H3K27ac-OPC_actenh_Glaser.png" \
# -ymax 10 -yzmax 10

# # G4 scCut&Tag - ESC_Enhancer_CruzMolina.active.bed (active enhancers, Glaser et al.)
# Rscript bw_heatmap_activeenh_Glaser.R -bw "../data/bw/G4_mES_H33WT_NC_R1_S33_L001.mm9.bw" -b "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed" -l "G4-H3.3WT" -o "../results/deeptools/G4-H3.3WT_actenh_Glaser.png" \
# -ymax 70 -yzmax 70
# Rscript bw_heatmap_activeenh_Glaser.R -bw "../data/bw/G4_mES_H33KO_NC_R1_S37_L001.mm9.bw" -b "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed" -l "G4-H3.3KO" -o "../results/deeptools/G4-H3.3KO_actenh_Glaser.png" \
# -ymax 70 -yzmax 70

# H3.3K27Ac ChIP-Seq - unique_G4_NPC.bed
Rscript bw_heatmap_NPC_G4s.R -bw "../data/bw/ESC_H33WT_H3K27ac.bw" -b "../data/bed/unique_G4_NPC.bed" -l "H3K27ac-H3.3WT" -o "../results/deeptools/H3.3WT_NPC_G4s.png" \
 -ymax 70 -yzmax 70
Rscript bw_heatmap_NPC_G4s.R -bw "../data/bw/ESC_H33KO_H3K27ac.bw" -b "../data/bed/unique_G4_NPC.bed" -l "H3K27ac-H3.3KO" -o "../results/deeptools/H3.3KO_NPC_G4s.png" \
 -ymax 70 -yzmax 70

