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

bigwig = args$bigwig
#bed = args$bed
label = args$label
#bed_label = args$bed_label
output = args$output
ymax = args$ymax
zmax = args$yzmax

# deeptools heatmap
system(glue("computeMatrix reference-point -S {bigwig} \\
-R  ../results/Seurat/callpeaks_unsorted/peak_sets/0_peaks.bed \\
      ../results/Seurat/callpeaks_unsorted/peak_sets/1_peaks.bed \\
      ../results/Seurat/callpeaks_unsorted/peak_sets/2_peaks.bed \\
      ../results/Seurat/callpeaks_unsorted/peak_sets/3_peaks.bed \\
      ../results/Seurat/callpeaks_unsorted/peak_sets/4_peaks.bed \\
      ../results/Seurat/callpeaks_unsorted/peak_sets/5_peaks.bed \\
      ../results/Seurat/callpeaks_unsorted/peak_sets/6_peaks.bed \\
-b 3000 \\
-a 3000 \\
--skipZeros \\
--missingDataAsZero \\
--samplesLabel {label} \\
-o ../results/deeptools/matrix.mat.gz"))

system(glue("plotHeatmap -m ../results/deeptools/matrix.mat.gz \\
-out {output} \\
--refPointLabel \"G4\" \\
--heatmapHeight 14 \\
--whatToShow \'heatmap and colorbar\' \\
--colorMap \"Reds\" \\
--yMin 0 \\
--yMax {ymax} \\
--zMax {zmax} \\
-z 0 1 2 3 4 5 6 \\
--yAxisLabel \"\" \\
--xAxisLabel \"\" \\
--legendLocation \"none\""))



