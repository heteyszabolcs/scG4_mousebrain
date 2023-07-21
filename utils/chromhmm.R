suppressPackageStartupMessages({
  library("tidyverse")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
})

# result / peak folder
result_folder = "../results/Seurat/callpeaks_unsorted/"
bigwig_folder = "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"

cluster_bigwigs = list.files(bigwig_folder, pattern = "*.bigwig", full.names = TRUE)
#cluster_bigwigs = grep('_res', cluster_bigwigs, value = TRUE, invert = TRUE)
chromhmm = plot_bw_loci_summary_heatmap(cluster_bigwigs, loci = "../data/bed/chromhmm_mESC_E14_12.bed")

ggsave(
  glue("{result_folder}chromhmm_mESC-unsorted_cls_mm10.pdf"),
  plot = chromhmm,
  width = 10,
  height = 10,
  device = "pdf"
)

