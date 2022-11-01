suppressPackageStartupMessages({
  library("tidyverse")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
})

# result / peak folder
result_folder = "../results/Seurat/callpeaks_unsorted/"
bigwig_folder = "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/"

cluster_bigwigs = list.files(bigwig_folder, pattern = "[0-9].bw", full.names = TRUE)
#cluster_bigwigs = grep('_res', cluster_bigwigs, value = TRUE, invert = TRUE)
chromhmm = plot_bw_loci_summary_heatmap(cluster_bigwigs, loci = "../data/bed/chromhmm_hindbrain_0_mm10_15_posterior.bed")

ggsave(
  glue("{result_folder}chromhmm_mESC-unsorted_cls-hindbrain_postnatal_mm10.pdf"),
  plot = chromhmm,
  width = 10,
  height = 10,
  device = "pdf"
)

