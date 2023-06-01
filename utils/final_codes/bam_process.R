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
args = parser$parse_args()

# add Seurat path
seurat = args$seurat_object
# open Seurat object
seurat = readRDS(seurat)
# add working directory
workdir = args$workdir
# mm10 effective genome size
mm10 = 2652783500

# peak calling
print("Cluster-specific peak calling (MACS2)")

system(paste0("mkdir -p ", workdir, "/cluster_spec_peaks"))

peaks = CallPeaks(
  object = seurat,
  group.by = "seurat_clusters",
  cleanup = FALSE,
  outdir = paste0(workdir, "/cluster_spec_peaks"),
  effective.genome.size = mm10
)

write.table(
  as.data.frame(peaks),
  file = glue("{workdir}/cluster_spec_peaks/peaks_per_clusters.bed"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)


# convert cluster specific bams to bigwigs
setwd(dir = paste0(workdir, "/cluster_spec_bams"))

print("Convert bam files to RPGC bigwigs (deeptools)")
system(paste0("for i in *.bam; do samtools index ${i}; done"))
system(paste0("for i in *.bam; do bamCoverage --bam ${i} -o ${i}_RPGC.bigwig --binSize 10 --normalizeUsing RPGC ",
       "--effectiveGenomeSize ", mm10, "; done"))

system(paste0("mkdir -p ", workdir, "/cluster_spec_bigwigs"))
system(paste0("mv *.bigwig ", workdir, "/cluster_spec_bigwigs"))


