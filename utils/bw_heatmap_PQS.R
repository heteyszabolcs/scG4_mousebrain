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
-R  \"../results/Seurat/callpeaks_unsorted/peak_sets/0_peaks_res0.1.bed\" \"../results/Seurat/callpeaks_unsorted/peak_sets/1_peaks_res0.1.bed\" \\
-b 500 \\
-a 500 \\
--skipZeros \\
--missingDataAsZero \\
--samplesLabel {label} \\
-o ../results/deeptools/matrix.mat.gz"))

system(glue("plotHeatmap -m ../results/deeptools/matrix.mat.gz \\
-out {output} \\
--refPointLabel \"G4\" \\
--heatmapHeight 14 \\
--whatToShow \'heatmap only\' \\
--legendLocation \"none\" \\
--colorMap \"binary\" \\
--yMin 0 \\
--yMax {ymax} \\
--zMax {zmax} \\
-z \"cluster0\" \"cluster1\" \\
--yAxisLabel \"\" \\
--xAxisLabel \"\" \\
--legendLocation \"none\""))



