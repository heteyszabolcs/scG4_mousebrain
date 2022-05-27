# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("GenomicFeatures")
  library("GenomicRanges")
  library("wigglescout")
  library("plyranges")
  library("ArchR")
})

# export folder
result_folder = "../results/GenomicRanges/"
peak_folder = "../results/Seurat/callpeaks_GFPsorted/"

# data
lanceotron_peaks = "../results/Seurat/callpeaks_GFPsorted/"
cm_enh = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
gl_enh = "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed"
ltrs = "../data/bed/RepMasker_lt200bp.LTRIS2.bed"

list_peaks = list.files(lanceotron_peaks, pattern = "_lanceotron*")

read_peaks = function(peak_set) {
  cluster = str_split(peak_set, "_")[[1]][1]
  peak_set = fread(glue("{lanceotron_peaks}{peak_set}"))
  peak_set = as_tibble(peak_set)
  peak_set = peak_set %>% mutate(Seurat_cluster = cluster)
  return(peak_set)
}

list_peaks = lapply(list_peaks, read_peaks)
lanceotron_peaks = bind_rows(list_peaks)
lanceotron_peaks_filt = lanceotron_peaks %>% filter(`Peak Score` >= 0.10) %>% 
  filter(str_detect(Chr, "chr"))

# export filtered Lanceotron peak set
write_tsv(lanceotron_peaks_filt, glue("{peak_folder}lanceotron_0.10filt_peak_set.tsv"))


hist = lanceotron_peaks_filt %>%
  ggplot(., aes(x = `Peak Score`, fill = Seurat_cluster)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "Blues") +
  labs(
    title = "Peak Score >= 0.25",
    x = "LanceOtron peak score",
    y = "Density",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
hist

ggsave(
  glue("{result_folder}Lanceotron_peak_score_distr.png"),
  plot = hist,
  width = 10,
  height = 10,
  dpi = 300,
)







