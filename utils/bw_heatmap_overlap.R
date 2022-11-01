#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("glue"))

# create parser object
parser <- ArgumentParser()

# add arguments and specify our desired options 
parser$add_argument("-bw", "--bigwig", type = "character",
                    help = "bigwig to be visualized")
parser$add_argument("-b", "--bed", type = "character",
                    help = "bed to subset")
parser$add_argument("-l", "--label", type = "character",
                    help = "y axis label")
parser$add_argument("-bl", "--bed_label", type = "character",
                    help = "bed file label")
parser$add_argument("-o", "--output", type = "character",
                    help = "output filename")
parser$add_argument("-ymax", "--ymax", type = "integer",
                    help = "ymax of profile plot", default = 150)
parser$add_argument("-yzmax", "--yzmax", type = "integer",
                    help = "z max of heatmap", default = 150)

args <- parser$parse_args()

#bigwig = args$bigwig
#bed = args$bed
#label = args$label
#bed_label = args$bed_label
output = args$output
ymax = args$ymax
zmax = args$yzmax

# deeptools heatmap
system(glue("computeMatrix reference-point -S \"../data/bw/bulk_G4_CnT_mESC_rep1_R1.mLb.clN_mm10.bigWig\" \\
 \"../data/bw/bulk_G4_CnT_MEF_rep1_R1.mLb.clN_mm10.bigWig\" \\
 \"../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/1_res0.1.bigwig\" \\
 \"../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/0_res0.1.bigwig\" \\
 -R  \"../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/only_mESC_bulk.bed\" \\ 
 \"../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/intersection.bed\" \\ 
 \"../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/only_MEF_bulk.bed\" \\
 -b 3000 \\
 -a 3000 \\
 --samplesLabel \"bulk\ mESC\" \"bulk\ MEF\" \"mESC\ cluster\ 1\" \"MEF\ cluster\ 0\" \\
 --skipZeros \\
 --missingDataAsZero \\
 -o ../results/deeptools/matrix.mat.gz"))

system(glue("plotHeatmap -m ../results/deeptools/matrix.mat.gz \\
-out {output} \\
--refPointLabel \"G4\" \\
--heatmapHeight 14 \\
--whatToShow \'heatmap only\' \\
--colorMap \"Blues\" \\
--yMin 0 \\
--yMax {ymax} \\
--zMax {zmax} \\
-z \"mESC\" \"both\" \"MEF\" \\
--yAxisLabel \"\" \\
--xAxisLabel \"\" \\
--legendLocation \"none\""))



