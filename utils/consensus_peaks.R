# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("GenomicFeatures")
  library("ChIPseeker")
  library("GenomicRanges")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
})

# result / peak folder
peak_folder = "../results/Seurat/callpeaks_GFPsorted/"
bed_folder = "../data/bed/"
result_folder = "../results/GenomicRanges/"
# peak sets
signac_macs2 = "enhancer_analysis_output.tsv"
archr = "archr_0.50perc_peak_set.tsv"
lanceotron = "lanceotron_0.10filt_peak_set.tsv"

# filtered MACS2 peak calls
signac_macs2 = fread(glue("{peak_folder}{signac_macs2}"))
signac_macs2$calling_strategy = "Signac MACS2"
signac_macs2 = GRanges(
  seqnames = signac_macs2$seqnames,
  ranges = IRanges(
    start = signac_macs2$start,
    end = signac_macs2$end,
    names = signac_macs2$calling_strategy,
    score = signac_macs2$signalValue
  )
)

# coverage plot
pdf(file = glue("{result_folder}Signac_MACS2_coverage_plot.pdf"),   # The directory you want to save the file in
    width = 5, 
    height = 5)
covplot(signac_macs2, weightCol = "score", title = "Signac's MACS2")
dev.off()

# filtered Archr peak calls
archr = fread(glue("{peak_folder}{archr}"))
archr$calling_strategy = "ArchR"
archr = GRanges(
  seqnames = archr$seqnames,
  ranges = IRanges(
    start = archr$start,
    end = archr$end,
    names = archr$calling_strategy,
    score = archr$score
  )
)

# coverage plot
pdf(file = glue("{result_folder}ArchR_MACS2_coverage_plot.pdf"),   # The directory you want to save the file in
    width = 5, 
    height = 5)
covplot(archr, weightCol = "score", title = "ArchR's MACS2")
dev.off()

# filtered lanceotron peak calls
lanceotron = fread(glue("{peak_folder}{lanceotron}"))
lanceotron$calling_strategy = "LanceOtron"
lanceotron = lanceotron %>% mutate(score = `Peak Score` * 10)
lanceotron = GRanges(
  seqnames = lanceotron$Chr,
  ranges = IRanges(
    start = lanceotron$Start,
    end = lanceotron$End,
    names = lanceotron$calling_strategy,
    score = lanceotron$score
  )
)

# coverage plot
pdf(file = glue("{result_folder}LanceOtron_coverage_plot.pdf"),   # The directory you want to save the file in
    width = 5, 
    height = 5)
covplot(lanceotron, weightCol = "score", title = "LanceOtron")
dev.off()

# Venn diagram of peak sets
peak_list = list(lanceotron, signac_macs2, archr)
names(peak_list) = c("LanceOtron", "Signac MACS2", "ArchR MACS2")
pdf(file = glue("{result_folder}peak_calling_VennPlot.pdf"),   # The directory you want to save the file in
    width = 8, 
    height = 8)
vennplot(peak_list)
dev.off()

# feature bars
# Signac MACS2
signac_macs2 = signac_macs2 %>% 
  as.tibble() %>% 
  distinct() %>% 
  mutate(peak_id = row_number())
signac_macs2 = GRanges(
  seqnames = signac_macs2$seqnames,
  ranges = IRanges(
    start = signac_macs2$start,
    end = signac_macs2$end,
    names = signac_macs2$peak_id,
    score = signac_macs2$score
  )
)
signac_macs2_anot = annotatePeak(signac_macs2, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
pdf(file = glue("{result_folder}Signac_MACS2_feature_distr.pdf"),   # The directory you want to save the file in
    width = 8, 
    height = 3)
plotAnnoBar(signac_macs2_anot, title = "Signac's MACS2")
dev.off()

# LanceOtron
lanceotron = lanceotron %>% 
  as.tibble() %>% 
  distinct() %>% 
  mutate(peak_id = row_number())
lanceotron = GRanges(
  seqnames = lanceotron$seqnames,
  ranges = IRanges(
    start = lanceotron$start,
    end = lanceotron$end,
    names = lanceotron$peak_id,
    score = lanceotron$score
  )
)
lanceotron_anot = annotatePeak(lanceotron, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
pdf(file = glue("{result_folder}lanceotron_feature_distr.pdf"),   # The directory you want to save the file in
    width = 8, 
    height = 3)
plotAnnoBar(lanceotron_anot, title = "LanceOtron")
dev.off()

# ArchR MACS2
archr = archr %>% 
  as.tibble() %>% 
  distinct() %>% 
  mutate(peak_id = row_number())
archr = GRanges(
  seqnames = archr$seqnames,
  ranges = IRanges(
    start = archr$start,
    end = archr$end,
    names = archr$peak_id,
    score = archr$score
  )
)
archr_anot = annotatePeak(archr, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
pdf(file = glue("{result_folder}archr_MACS2_feature_distr.pdf"),   # The directory you want to save the file in
    width = 8, 
    height = 3)
plotAnnoBar(archr_anot, title = "ArchR's MACS2")
dev.off()

# create intersection by GenomicRanges
consensus = as.tibble(Reduce(subsetByOverlaps, list(signac_macs2, archr, lanceotron)))

signac_macs2 = "enhancer_analysis_output.tsv"
signac_macs2 = fread(glue("{peak_folder}{signac_macs2}"))

consensus = consensus %>% inner_join(.,
                                     signac_macs2,
                                     by = c(
                                       "seqnames" = "seqnames",
                                       "start" = "start",
                                       "end" = "end"
                                     ))

bed = consensus %>% dplyr::select(seqnames, start, end, score)
# consensus bed file
write_tsv(bed, glue("{bed_folder}consensus_G4s.bed"), col_names = FALSE)

# visualize scCut&Tag epigenetic marks over consensus G4 peak regions
h3k27ac_hist = consensus %>% pivot_longer(
  .,
  cols = starts_with("H3K27ac"),
  names_to = "scCut&Tag data",
  values_to = "score"
) %>% dplyr::select(`scCut&Tag data`, score) %>% 
  ggplot(., aes(x = score, fill = `scCut&Tag data`)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "Greens") +
  xlim(0, 60) +
  labs(
    title = "H3K27ac over common G4 peaks",
    x = "H3K27ac signalValue (MACS2)",
    y = "Density",
    fill = "scCut&Tag data"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
h3k27ac_hist

ggsave(
  glue("{result_folder}Consensus_peaks-H3K27ac_distr.png"),
  plot = h3k27ac_hist,
  width = 10,
  height = 10,
  dpi = 300,
)

h3k4me3_hist = consensus %>% pivot_longer(
  .,
  cols = starts_with("H3K4me3"),
  names_to = "scCut&Tag data",
  values_to = "score"
) %>% dplyr::select(`scCut&Tag data`, score) %>% 
  ggplot(., aes(x = score, fill = `scCut&Tag data`)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "Greens") +
  labs(
    title = "H3K4me3 over common G4 peaks",
    x = "H3K4me3 signalValue (MACS2)",
    y = "Density",
    fill = "scCut&Tag data"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
h3k4me3_hist

ggsave(
  glue("{result_folder}Consensus_peaks-H3K4me3_distr.png"),
  plot = h3k4me3_hist,
  width = 10,
  height = 10,
  dpi = 300,
)


