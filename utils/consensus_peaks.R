# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("GenomicFeatures")
  library("ChIPseeker")
  library("ChIPpeakAnno")
  library("Vennerable")
  library("GenomicRanges")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("biomaRt")
})

# result / peak folder
peak_folder = "../results/Seurat/callpeaks_GFPsorted/peak_sets/"
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
    score = signac_macs2$signalValue,
    cluster = signac_macs2$Seurat_cluster
  )
)

# coverage plot
pdf(
  file = glue("{result_folder}Signac_MACS2_coverage_plot.pdf"),
  # The directory you want to save the file in
  width = 5,
  height = 5
)
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
    score = archr$score,
    replicate = archr$GroupReplicate
  )
)

# coverage plot
pdf(
  file = glue("{result_folder}ArchR_MACS2_coverage_plot.pdf"),
  # The directory you want to save the file in
  width = 5,
  height = 5
)
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
    score = lanceotron$score,
    cluster = lanceotron$Seurat_cluster
  )
)

# coverage plot
pdf(
  file = glue("{result_folder}LanceOtron_coverage_plot.pdf"),
  # The directory you want to save the file in
  width = 5,
  height = 5
)
covplot(lanceotron, weightCol = "score", title = "LanceOtron")
dev.off()

# Venn diagram of peak sets
# peak_list = list(lanceotron, signac_macs2, archr)
# names(peak_list) = c("LanceOtron", "Signac MACS2", "ArchR MACS2")
# pdf(
#   file = glue("{result_folder}peak_calling_VennPlot.pdf"),
#   # The directory you want to save the file in
#   width = 8,
#   height = 8
# )
# vennplot(peak_list)
# dev.off()

# retrieve unique peaks 
l_unique = lanceotron[-queryHits(findOverlaps(lanceotron, signac_macs2, type = "any", minoverlap = 1)),] 
l_unique = l_unique[-queryHits(findOverlaps(l_unique, archr, type = "any", minoverlap = 1)),] 
l_unique = as_tibble(l_unique)
l_unique = l_unique %>% mutate(peak_calling = "LanceOtron")
write_tsv(l_unique,
          glue("{peak_folder}lanceotron_unique_peaks.tsv"),
          col_names = FALSE)

a_unique = archr[-queryHits(findOverlaps(archr, signac_macs2, type = "any", minoverlap = 1)),] 
a_unique = a_unique[-queryHits(findOverlaps(a_unique, lanceotron, type = "any", minoverlap = 1)),] 
a_unique = as_tibble(a_unique)
a_unique = a_unique %>% mutate(peak_calling = "ArchR MACS2")
write_tsv(a_unique,
          glue("{peak_folder}archr_unique_peaks.tsv"),
          col_names = FALSE)

s_unique = signac_macs2[-queryHits(findOverlaps(signac_macs2, lanceotron, type = "any", minoverlap = 1)),] 
s_unique = s_unique[-queryHits(findOverlaps(s_unique, archr, type = "any", minoverlap = 1)),] 
s_unique = as_tibble(s_unique)
s_unique = s_unique %>% mutate(peak_calling = "Signac MACS2")
write_tsv(s_unique,
          glue("{peak_folder}signac_macs2_unique_peaks.tsv"),
          col_names = FALSE)

# create Venn diagram
signac_macs2$type = "Signac MACS2"
lanceotron$type = "LanceOtron"
archr$type = "ArchR"
gr = c(lanceotron, signac_macs2, archr)
grl = splitAsList(gr, gr$type)
grl = unique(grl)

res = makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl))

venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Venn(SetNames=SetNames, Weight=Weight)
}

pdf(
  file = glue("{result_folder}peak_calling_VennPlot.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 8
)
v <- venn_cnt2venn(res$vennCounts)
plot(v, doWeights = FALSE)
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
pdf(
  file = glue("{result_folder}Signac_MACS2_feature_distr.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 3
)
plotAnnoBar(signac_macs2_anot, title = "Signac MACS2")
dev.off()

signac_macs2_tibble = as.tibble(signac_macs2_anot@anno)
signac_macs2_distal_igenic = signac_macs2_tibble %>% filter(annotation == "Distal Intergenic") %>%
  mutate(norm_score = score / max(score)) %>%
  mutate(peak_calling = "Signac MACS2")
signac_macs2_promoter = signac_macs2_tibble %>% filter(str_detect(annotation, "Promoter")) %>%
  mutate(norm_score = score / max(score)) %>%
  mutate(peak_calling = "Signac MACS2")


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
pdf(
  file = glue("{result_folder}lanceotron_feature_distr.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 3
)
plotAnnoBar(lanceotron_anot, title = "LanceOtron")
dev.off()

lanc_tibble = as.tibble(lanceotron_anot@anno)
lanc_distal_igenic = lanc_tibble %>% filter(annotation == "Distal Intergenic") %>%
  mutate(norm_score = score / max(score)) %>%
  mutate(peak_calling = "LanceOtron")
lanc_promoter = lanc_tibble %>% filter(str_detect(annotation, "Promoter")) %>%
  mutate(norm_score = score / max(score)) %>%
  mutate(peak_calling = "LanceOtron")

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
pdf(
  file = glue("{result_folder}archr_MACS2_feature_distr.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 3
)
plotAnnoBar(archr_anot, title = "ArchR MACS2")
dev.off()

archr_tibble = as.tibble(archr_anot@anno)
archr_distal_igenic = archr_tibble %>%
  filter(annotation == "Distal Intergenic") %>%
  mutate(norm_score = score / max(score)) %>%
  mutate(peak_calling = "ArchR MACS2")
archr_promoter = archr_tibble %>% filter(str_detect(annotation, "Promoter")) %>%
  mutate(norm_score = score / max(score)) %>%
  mutate(peak_calling = "ArchR MACS2")

# create merged table about Distal region peaks
distal_igenic_all = rbind(archr_distal_igenic,
                          lanc_distal_igenic,
                          signac_macs2_distal_igenic)
# create merged table about Promoter peaks
promoter_all = rbind(archr_promoter, lanc_promoter, signac_macs2_promoter)

# export bed files
bed_distal = distal_igenic_all %>% dplyr::select(seqnames, start, end, score)
bed_prom = promoter_all %>% dplyr::select(seqnames, start, end, score)
write_tsv(
  bed_distal,
  glue("{bed_folder}ArchRSignacLanc_distal_igenic_G4_peaks.bed"),
  col_names = FALSE
)
write_tsv(bed_prom,
          glue("{bed_folder}ArchRSignacLanc_promoter_G4_peaks.bed"),
          col_names = FALSE)

# annotation of entrez ids of promoter peaks
ensembl = useEnsembl(biomart = "genes")
datasets = listDatasets(ensembl)
# filters = listFilters(ensembl)
# attributes = listAttributes(ensembl)
ensembl = useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
annot = getBM(attributes = "mgi_symbol",
      filters = 'entrezgene_id',
      values = promoter_all$geneId, 
      mart = ensembl)

write_tsv(promoter_all,
          glue("{peak_folder}ArchRSignacLanc_promoter_G4_peaks.tsv"),
          col_names = FALSE)
write_tsv(annot,
          glue("{peak_folder}ArchRSignacLanc_promoter_G4_peaks_genes.tsv"),
          col_names = FALSE)


# consensus bed file
hist = distal_igenic_all %>%
  ggplot(., aes(x = norm_score, fill = peak_calling)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "Set3") +
  xlab(c(0, 1)) +
  labs(title = "",
       x = "norm. peak score",
       y = "Density",
       fill = "peak calling method") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
hist

ggsave(
  glue("{result_folder}Norm_peak_score_distr.png"),
  plot = hist,
  width = 10,
  height = 10,
  dpi = 300,
)

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
                                     )) %>%
  dplyr::select(
    "seqnames",
    "start",
    "end",
    "score",
    "Distance_to_TSS",
    "Gene_name",
    starts_with("H3K27ac"),
    starts_with("H3K4me3")
  )

bed = consensus %>% dplyr::select(seqnames, start, end, score)
# consensus bed file
write_tsv(bed, glue("{bed_folder}consensus_G4s.bed"), col_names = FALSE)
write_tsv(consensus, glue("{peak_folder}consensus_G4s.tsv"))
gene_names = consensus %>% dplyr::select(Gene_name)
write_tsv(gene_names, glue("{peak_folder}consensus_G4s_genes.tsv"))

# visualize scCut&Tag epigenetic marks over consensus G4 peak regions
h3k27ac_hist = consensus %>% pivot_longer(
  .,
  cols = starts_with("H3K27ac"),
  names_to = "scCut&Tag data",
  values_to = "H3K27ac"
) %>%
  dplyr::select(`scCut&Tag data`, H3K27ac) %>%
  ggplot(., aes(x = H3K27ac, fill = `scCut&Tag data`)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "Greens") +
  xlim(0, 60) +
  labs(title = "H3K27ac over common G4 peaks",
       x = "H3K27ac signalValue (MACS2)",
       y = "Density",
       fill = "scCut&Tag data") +
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
  values_to = "H3K4me3"
) %>% dplyr::select(`scCut&Tag data`, H3K4me3) %>%
  ggplot(., aes(x = H3K4me3, fill = `scCut&Tag data`)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "Greens") +
  labs(title = "H3K4me3 over common G4 peaks",
       x = "H3K4me3 signalValue (MACS2)",
       y = "Density",
       fill = "scCut&Tag data") +
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
