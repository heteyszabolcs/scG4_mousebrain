# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("ggrepel")
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
archr_peaks = "../data/ArchR_peaks/mouse_brain/"
cm_enh = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
gl_enh = "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed"
ltrs = "../data/bed/RepMasker_lt200bp.LTRIS2.bed"

# create gr objects
cm_enh = fread(cm_enh)
cm_enh$V6 = "Cruz-Molina_active_enh"
cm_enh = GRanges(
  seqnames = cm_enh$V1,
  ranges = IRanges(
    start = cm_enh$V2,
    end = cm_enh$V3,
    names = cm_enh$V6
  )
)

gl_enh = fread(gl_enh)
gl_enh$V6 = "Glaser_active_enh"
gl_enh = GRanges(
  seqnames = gl_enh$V1,
  ranges = IRanges(
    start = gl_enh$V2,
    end = gl_enh$V3,
    names = gl_enh$V6
  )
)

list_peaks = list.files(archr_peaks, pattern = "*.rds")
read_peaks = function(peak_set) {
  peak_set = readRDS(glue("{archr_peaks}{peak_set}"))
  peak_set = as_tibble(peak_set)
  return(peak_set)
}
                 
list_peaks = lapply(list_peaks, read_peaks)
archr_peaks = bind_rows(list_peaks)
archr_peaks = archr_peaks %>% filter(score >= quantile(score, .90))


bar = archr_peaks %>% 
  group_by(GroupReplicate) %>% 
  summarise(count = n()) %>% 
  filter(str_detect(GroupReplicate, "GFP")) %>% 
  ggplot(data = ., aes(
    x = reorder(GroupReplicate, -count),
    y = count,
    fill = GroupReplicate
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "Greens") +
  labs(
    title = 
      "ArchR peak calling per Seurat cluster (without filters)",
    x = "Seurat cluster",
    y = "# of G4 peaks",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
bar

ggsave(
  glue("{result_folder}ArchR_peak_bars.png"),
  plot = bar,
  width = 10,
  height = 10,
  dpi = 300,
)

bar_regulatory = archr_peaks %>% 
  group_by(peakType) %>% 
  summarise(count = n()) %>% 
  ggplot(data = ., aes(
    x = reorder(peakType, -count),
    y = count,
    fill = peakType
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "Greens") +
  labs(
    title = 
      "ArchR peak types",
    x = "Seurat cluster",
    y = "# of G4 peaks",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
bar_regulatory

ggsave(
  glue("{result_folder}ArchR_peak_type_bars.png"),
  plot = bar_regulatory,
  width = 10,
  height = 10,
  dpi = 300,
)

hist = archr_peaks %>%
  filter(str_detect(GroupReplicate, "GFP")) %>% 
  ggplot(., aes(x = score, fill = GroupReplicate)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "Greens") +
  labs(
    title = "",
    x = "ArchR peak score",
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
  glue("{result_folder}ArchR_peak_score_distr.png"),
  plot = hist,
  width = 10,
  height = 10,
  dpi = 300,
)

## complete ArchR peaks with enhancer signatures
archr_peaks = archr_peaks %>% mutate(name = "ArchR peak")
archr_peaks_gr = GRanges(
  seqnames = archr_peaks$seqnames,
  ranges = IRanges(
    start = archr_peaks$start,
    end = archr_peaks$end,
    names = cm_enh$name
  )
)

# intersect with Cruz-Molina enhancers
cm_int = suppressWarnings(findOverlaps(archr_peaks_gr, cm_enh, type = "any", minoverlap = 1))
archr_peaks = archr_peaks %>% mutate(Cruz_Molina_enh = "0")
for (hit in cm_int@from) {
  archr_peaks[hit, which(colnames(archr_peaks) == "Cruz_Molina_enh")] = "1"
}
# intersect with Glaser enhancers
gl_int = suppressWarnings(findOverlaps(archr_peaks_gr, gl_enh, type = "any", minoverlap = 1))
archr_peaks = archr_peaks %>% mutate(Glaser_enh = "0")
for (hit in gl_int@from) {
  archr_peaks[hit, which(colnames(archr_peaks) == "Glaser_enh")] = "1"
}
  
peak_rank = archr_peaks %>%
  filter(str_detect(GroupReplicate, "GFP")) %>% 
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-score), y = score)) +
  geom_point(
    stat = 'identity',
    aes(col = GroupReplicate),
    size = 3,
    alpha = 0.4
  ) +
  scale_color_brewer(palette = "Greens") +
  labs(
    title = "ArchR peak scores",
    x = "ArchR peaks",
    y = "ArchR peak score",
    color = "Seurat cluster"
  ) +
  geom_label_repel(
    aes(label = ifelse(
      Cruz_Molina_enh == 1 |
        Glaser_enh == 1,
      as.character(nearestGene),
      ''
    )),
    box.padding   = 0.9,
    max.overlaps = Inf,
    point.padding = 0.5,
    size = 2,
    segment.color = 'black'
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 20),
    plot.title = element_text(size = 15)
  )
peak_rank

ggsave(
  glue("{result_folder}ArchR_peak_rank.png"),
  plot = peak_rank,
  width = 10,
  height = 10,
  dpi = 300,
)

# filter ArchR peaks based on X percentile
archr_peaks_filt = archr_peaks %>% filter(score >= quantile(score, .50))
write_tsv(archr_peaks_filt, glue("{peak_folder}archr_0.50perc_peak_set.tsv"))
  
  

