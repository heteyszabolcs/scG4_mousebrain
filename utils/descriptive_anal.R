# packages
suppressPackageStartupMessages({
  library("data.table")
  library("glue")
  library("tidyverse")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("GenomicFeatures")
  library("GenomeInfoDb")
  library("cowplot")
  #library("DoubletFinder")
})

# path to result folder
result_folder = "../results/Seurat/callpeaks_mESC-MEF/"

peaks = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/0_peaks_res0.1.bed")

txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

annot1 = annotatePeak("../results/Seurat/callpeaks_mESC-MEF/peak_sets/0_peaks_res0.1.bed", 
                     tssRegion = c(-1000, 1000), TxDb = txdb, verbose = FALSE)
annot2 = annotatePeak("../results/Seurat/callpeaks_mESC-MEF/peak_sets/1_peaks_res0.1.bed", 
                     tssRegion = c(-1000, 1000), TxDb = txdb, verbose = FALSE)

df1 = annot1@annoStat
df1 = df1 %>% mutate(cluster = "mESC−MEF G4 scCnT (cluster 0)")
df2 = annot2@annoStat
df2 = df2 %>% mutate(cluster = "mESC−MEF G4 scCnT (cluster 1)")

df = rbind(df1, df2)

ggplot(df, aes(fill = Feature, y = Frequency, x = cluster)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette = "Spectral") +
  labs(
    title = "",
    x = "",
    y = "Frequency",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  ) 

ggsave(
  glue("{result_folder}mESC-MEF_plotAnnot-stacked.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 300,
)

ggsave(
  glue("{result_folder}mESC-MEF_plotAnnot-stacked.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8,
  device = "pdf"
)

annot1_df = as.data.frame(annot1@anno)

options(scipen=10000)
hist = ggplot(annot1_df, aes(x = distanceToTSS)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  xlim(-100000, 100000) +
  #scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "Distance to TSS",
       x = "",
       y = "Density",
       fill = "") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black"))
hist
    
