print("Open virtual environment")
# load virtual environment
system("source /home/szabolcs/.pyenv/versions/3.8.1/envs/sinto/bin/activate")

print("Load R packages")
# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("argparse")
  library("tidyverse")
})

# create parser object
parser = ArgumentParser()

parser$add_argument("-s", "--seurat_object", type = "character",
                    help = "path to Seurat object with clusters")
parser$add_argument("-w", "--workdir", type = "character",
                    help = "path of working dir")
parser$add_argument("-b", "--bam", type = "character",
                    help = "path to sorted bam file of cellranger")
args = parser$parse_args()

# add Seurat path
seurat = args$seurat_object
# open Seurat object
seurat = readRDS(seurat)
# add working directory
workdir = args$workdir
# bam
bam = args$bam

# mm10 effective genome size
mm10 = 2652783500

# create barcodes
print("Retrieve cluster specific barcodes")

system(paste0("mkdir -p ", workdir, "/cluster_spec_barcodes"))

barcodes = seurat@meta.data %>%
  rownames_to_column("barcodes") %>%
  dplyr::select(barcodes, seurat_clusters)
write_tsv(barcodes, glue("{workdir}/cluster_spec_barcodes/barcodes_per_cluster.tsv"))

for(cluster in unique(barcodes$seurat_clusters)) {
  subset = barcodes %>% dplyr::filter(seurat_clusters == cluster)
  write_tsv(subset, glue("{workdir}/cluster_spec_barcodes/barcodes_cluster_{as.character(cluster)}.tsv"), col_names = FALSE)
}

# sinto
print("Generate cluster specific bam files")

barcode_folder = list.files(glue("{workdir}/cluster_spec_barcodes/"), full.names = TRUE, pattern = "*barcodes_cluster_*")

for(barcodes in barcode_folder) {
  print(paste0("Working on: ", barcodes))
  sinto = glue("sinto filterbarcodes -b {bam} -c {barcodes} --outdir {workdir}/cluster_spec_barcodes")
  system(sinto)
}

system(paste0("mkdir -p ", workdir, "/cluster_spec_bams"))
system(paste0("mv ", workdir, "/cluster_spec_barcodes/", "*.bam ", workdir, "/cluster_spec_bams"))

# # peak calling
# system(paste0("mkdir -p ", workdir, "/cluster_spec_peaks"))
# 
# peaks = CallPeaks(
#   object = seurat,
#   group.by = "seurat_clusters",
#   cleanup = FALSE,
#   outdir = paste0(workdir, "/cluster_spec_peaks"),
#   effective.genome.size = mm10
# )
# 
# write.table(
#   as.data.frame(peaks),
#   file = glue("{workdir}/cluster_spec_peaks/peaks_per_clusters.bed"),
#   quote = F,
#   sep = "\t",
#   row.names = F,
#   col.names = F
# )
# 
# 
# # convert cluster specific bams to bigwigs
# setwd(dir = paste0(workdir, "/cluster_spec_bams"))
# 
# print("Convert bam files to RPGC bigwigs")
# system(paste0("for i in *.bam; do samtools index ${i}; done"))
# system(paste0("for i in *.bam; do bamCoverage --bam ${i} -o ${i}_RPGC.bigwig --binSize 10 --normalizeUsing RPGC ",
#        "--effectiveGenomeSize ", mm10, "; done"))
# 
# system(paste0("mkdir -p ", workdir, "/cluster_spec_bigwigs"))
# system(paste0("mv *.bigwig ", workdir, "/cluster_spec_bigwigs"))





