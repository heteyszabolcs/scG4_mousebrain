#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("glue"))

# bw_path = "../data/LTRIS2_BRG1_H33_and_G4s/bw/"
# bed_path = "../data/LTRIS2_BRG1_H33_and_G4s/bed/"

# bigwigs = list.files(bw_path)
# beds = list.files(bed_path)

# create parser object
parser <- ArgumentParser()

# add arguments and specify our desired options 
parser$add_argument("-bw", "--bigwig", type = "character",
                    help = "bigwig to be visualized")
parser$add_argument("-b", "--bed", type = "character",
                    help = "bed to subset")
parser$add_argument("-l", "--label", type = "character",
                    help = "y axis label")
parser$add_argument("-o", "--output", type = "character",
                    help = "output filename")
parser$add_argument("-ymax", "--ymax", type = "integer",
                    help = "ymax of profile plot", default = 150)
parser$add_argument("-yzmax", "--yzmax", type = "integer",
                    help = "z max of heatmap", default = 150)

args <- parser$parse_args()

# args$bigwig = "../data/bw/ATAC-H33KO-H32rescue-802.mm9.bw"
# args$bed = "../data/bed/H33_dependent_G4_JL.v2.bed"
# args$label = "H3.3KO-H3.2rescue"
# args$output = "../results/deeptools/H3.3KO-H3.2rescue_G4_hm.png"

bigwig = args$bigwig
bed = args$bed
label = args$label
output = args$output
ymax = args$ymax
zmax = args$yzmax

# deeptools heatmap
system(glue("computeMatrix reference-point -S {bigwig} \\
-R  {bed} \\
-b 3000 \\
-a 3000 \\
--samplesLabel {label} \\
--skipZeros -o /proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/deeptools/matrix.mat.gz"))

system(glue("plotHeatmap -m ../results/deeptools/matrix.mat.gz \\
-out {output} \\
--refPointLabel \"act.enh.\" \\
--heatmapHeight 14 \\
--colorMap \"Reds\" \\
--yMin 0 \\
--yMax {ymax} \\
--zMax {zmax} \\
--yAxisLabel \"\" \\
--xAxisLabel \"\" \\
--regionsLabel \"active enh (Glaser et al.)\" \\
--legendLocation \"none\""))


