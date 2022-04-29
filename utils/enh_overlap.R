# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("GenomicFeatures")
  library("GenomicRanges")
})

result_folder = "../results/GenomicRanges/"

g4_peaks = fread("../results/Seurat/callpeaks_GFPsorted/peaks_per_clusters.bed")
g4_peaks$V7 = "G4"
# G4 length distributions
hist = g4_peaks %>% mutate(diff = V3 - V2) %>%
  filter(V6 %in% c("0", "1", "2", "3", "4")) %>%
  ggplot(., aes(x = diff, fill = V6)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "G4 length distributions / Seurat cluster", 
       x = "length (bp)", 
       y = "Density",
       fill = "Seurat cluster") +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size = 13, color = "black")) 
hist

ggsave(
  glue("{result_folder}G4_length_distr_Seurat_clst.png"),
  plot = hist,
  width = 10,
  height = 10,
  dpi = 300,
)

# convert to GRanges
g4_peaks =
  lapply(split(g4_peaks, g4_peaks$V6), function(i) {
    GRanges(seqnames = i$V1,
            ranges = IRanges(
              start = i$V2,
              end = i$V3,
              names = i$V7
            ))
  })

# enhancers
cm_enh = fread("../data/bed/ESC_Enhancer_CruzMolina.active.bed")
cm_enh$V6 = "Cruz-Molina_active_enh"
cm_enh = GRanges(
  seqnames = cm_enh$V1,
  ranges = IRanges(
    start = cm_enh$V2,
    end = cm_enh$V3,
    names = cm_enh$V6
  )
)

# query: g4_peaks
# subject: enhancers

g4s_no_ol = g4_peaks[c("0", "1", "2", "3", "4")]
g4_enh_quant = numeric()

for(cluster in names(g4s_no_ol)) {
  ol = suppressWarnings(findOverlaps(g4_peaks[[cluster]], cm_enh, type = "any", minoverlap = 1))
  g4_enh_quant = c(g4_enh_quant, length(ol))
}

bar = tibble(overlap = g4_enh_quant, cluster = names(g4s_no_ol), fill_col = names(g4s_no_ol)) %>% 
  ggplot(data = ., aes(x = reorder(cluster, -overlap), y = overlap, fill = fill_col)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = expression(paste("G4 overlaps with active enhancers of ", italic("Cruz-Molina et al."))), 
       x = "Seurat cluster", 
       y = "# of G4 - active enhancer overlaps",
       fill = "Seurat cluster") +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(size = 13, color = "black")) 
bar  

ggsave(
  glue("{result_folder}G4_overlaps_w_CruzM_active_enh.png"),
  plot = bar,
  width = 10,
  height = 10,
  dpi = 300,
)






