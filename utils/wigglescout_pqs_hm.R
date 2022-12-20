suppressPackageStartupMessages({
  library("tidyverse")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
})


# result / peak folder
result_folder = "../results/Seurat/callpeaks_unsorted/"
peak_folder = "../results/Seurat/callpeaks_unsorted/peak_sets/"
bw_folder = "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/"

# GSM3003547 G4-Seq peaks
g4_seq = fread("../data/bed/GSM3003547_G4seq_mm10_all_w15_th-1_minus.hits.max.K.w50.25_sorted.bed")
g4_seq = g4_seq %>% arrange(desc(V4)) %>% top_n(100, wt = V4)
write_tsv(g4_seq, glue("{result_folder}GSM3003547_G4seq_mm10_all_w15_th-1_minus.hits.max.K.w50.25_top100.bed"),
          col_names = FALSE)
g4_seq = glue("{result_folder}GSM3003547_G4seq_mm10_all_w15_th-1_minus.hits.max.K.w50.25_top100.bed")

# PQS mm10 bigwig
pqs = "../data/bw/mm10_canPQS-regex_binary.bw"

# scG4 peak set
peaks = fread(glue("{peak_folder}cluster_1_res0.1_Lanc_peaks.tsv"))
peaks = peaks %>% arrange(desc(`Peak Score`)) %>% top_n(100, wt = `Peak Score`) %>% 
  select("Chr", "Start", "End")
write_tsv(peaks, glue("{result_folder}cluster_1_res0.1_Lanc_top_100.bed"), col_names = FALSE)

bed = glue("{result_folder}cluster_1_res0.1_Lanc_top_100.bed")

# scG4 bigwig
signal = glue("{bw_folder}1_res0.1_brain_cells.bw")


hm = plot_bw_heatmap(
  pqs,
  loci = bed,
  mode = "center",
  upstream = 10000,
  downstream = 10000,
  default_na = 0,
  bin_size = 100,
  cmap = "Greens"
)
hm

ggsave(
  glue("{result_folder}PQS_hm_top_Lanc_peaks_cl1.pdf"),
  width = 7,
  height = 7,
  device = "pdf"
)

hm2 = plot_bw_heatmap(
  pqs,
  loci = g4_seq,
  mode = "center",
  upstream = 10000,
  downstream = 10000,
  default_na = 100,
  cmap = "Greens"
)
hm2
