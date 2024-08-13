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
  library("cotools")
  library("ComplexHeatmap")
  library("circlize")
  library("wigglescout")
})

## Venn-diagrams
# Venn with 3 sets - bulk mESC and cluster0/1 of mESC-MEF scCnT
result_folder = "../results/GenomicRanges/mESC-MEF_outputs/"

cluster0 = fread("../results/Seurat/final/mesc_mef/cluster_spec_peaks/0_peaks.narrowPeak")
cluster0$type = "cluster0"
cluster0 = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    names = cluster0$type,
  )
)

cluster1 = fread("../results/Seurat/final/mesc_mef/cluster_spec_peaks/1_peaks.narrowPeak")
cluster1$type = "cluster1"
cluster1 = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    names = cluster1$type,
  )
)

bulk_mesc_peak = fread("../data/bed/bulk_G4_CnT_mESC_rep1_R1.mLb.clN_peaks.broadPeak")
bulk_mesc_peak$type = "mESC_bulk_CnT"
bulk_mesc_peak = GRanges(
  seqnames = bulk_mesc_peak$V1,
  ranges = IRanges(
    start = bulk_mesc_peak$V2,
    end = bulk_mesc_peak$V3,
    names = bulk_mesc_peak$type,
  )
)

cluster0$type = "cluster_0"
cluster1$type = "cluster_1"
bulk_mesc_peak$type = "mESC_bulk_CnT"
gr = c(cluster0, cluster1, bulk_mesc_peak)
grl = splitAsList(gr, gr$type)
grl = unique(grl)

res = makeVennDiagram(Peaks = grl, NameOfPeaks = names(grl))

venn_cnt2venn <- function(venn_cnt) {
  n <- which(colnames(venn_cnt) == "Counts") - 1
  SetNames = colnames(venn_cnt)[1:n]
  Weight = venn_cnt[, "Counts"]
  names(Weight) <- apply(venn_cnt[, 1:n], 1, paste, collapse = "")
  Venn(SetNames = SetNames, Weight = Weight)
}

pdf(
  file = glue("{result_folder}bulk_mESC_MEF-MESCcluster0_1.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 8
)
v <- venn_cnt2venn(res$vennCounts)
print(plot(v, doWeights = FALSE))
dev.off()

# Venn with 3 sets - bulk MEF and cluster0/1 of mESC-MEF scCnT
cluster0 = fread("../results/Seurat/final/mesc_mef/cluster_spec_peaks/0_peaks.narrowPeak")
cluster0$type = "cluster0"
cluster0 = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    names = cluster0$type,
  )
)

cluster1 = fread("../results/Seurat/final/mesc_mef/cluster_spec_peaks/1_peaks.narrowPeak")
cluster1$type = "cluster1"
cluster1 = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    names = cluster1$type,
  )
)

bulk_mef_peak = fread("../data/bed/bulk_G4_CnT_MEF_rep1_R1.mLb.clN_peaks.broadPeak")
bulk_mef_peak$type = "MEF_bulk_CnT"
bulk_mef_peak = GRanges(
  seqnames = bulk_mef_peak$V1,
  ranges = IRanges(
    start = bulk_mef_peak$V2,
    end = bulk_mef_peak$V3,
    names = bulk_mef_peak$type,
  )
)

cluster0$type = "cluster_0"
cluster1$type = "cluster_1"
bulk_mef_peak$type = "MEF_bulk_CnT"
gr = c(cluster0, cluster1, bulk_mef_peak)
grl = splitAsList(gr, gr$type)
grl = unique(grl)

res = makeVennDiagram(Peaks = grl, NameOfPeaks = names(grl))

venn_cnt2venn <- function(venn_cnt) {
  n <- which(colnames(venn_cnt) == "Counts") - 1
  SetNames = colnames(venn_cnt)[1:n]
  Weight = venn_cnt[, "Counts"]
  names(Weight) <- apply(venn_cnt[, 1:n], 1, paste, collapse = "")
  Venn(SetNames = SetNames, Weight = Weight)
}

pdf(
  file = glue("{result_folder}bulk_MEFC_MEF-MEFcluster0_1.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 8
)
v <- venn_cnt2venn(res$vennCounts)
print(plot(v, doWeights = FALSE))
dev.off()

## mESC - MEF - G4 scCut&Tag Jaccard analysis
# jaccard indices (package cotools) and GenometriCorrelation analysis (package GenometriCorr)
peaks = list(cluster0, cluster1, bulk_mesc_peak, bulk_mef_peak)
jaccards = matrix(NA_real_, length(peaks), length(peaks))
colnames(jaccards) = c("cluster 0", "cluster 1", "bulk mESC CnT", "bulk MEF CnT")
rownames(jaccards) = c("cluster 0", "cluster 1", "bulk mESC CnT", "bulk MEF CnT")
for(i in seq(1, ncol(jaccards))) {
  for(j in seq(1, nrow(jaccards))) {
    jaccard = genomicCorr.jaccard(peaks[[i]], peaks[[j]])
    jaccards[i, j] = jaccard
  }
}
is.matrix(jaccards)

pdf(
  file = "../results/GenomicRanges/mESC-MEF_outputs/mESC_MEF_res0.1_Jaccard_hm.pdf",
  width = 6,
  height = 6
)
col_fun = colorRamp2(c(0, 0.25, 0.5), c("#9ecae1", "white", "#fc9272"))
Heatmap(
  jaccards,
  column_title = "",
  row_title = "",
  name = "Jaccard index",
  # row_km = 2,
  # column_km = 1,
  clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(8, "cm"),
  heatmap_height = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", jaccards[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

# Venn with 3 sets - bulk NPC and cluster0/1 of unsorted scCnT
result_folder = "../results/GenomicRanges/unsorted_outputs/"

cluster0_unsorted = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")
cluster0_unsorted$type = "cluster0"
cluster0_unsorted = GRanges(
  seqnames = cluster0_unsorted$V1,
  ranges = IRanges(
    start = cluster0_unsorted$V2,
    end = cluster0_unsorted$V3,
    names = cluster0_unsorted$type,
  )
)

cluster1_unsorted = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
cluster1_unsorted$type = "cluster1"
cluster1_unsorted = GRanges(
  seqnames = cluster1_unsorted$V1,
  ranges = IRanges(
    start = cluster1_unsorted$V2,
    end = cluster1_unsorted$V3,
    names = cluster1_unsorted$type,
  )
)

bulk_npc_peak = fread("../data/bed/bulk_CnT_G4_NPC_mm10.bed")
bulk_npc_peak$type = "NPC_bulk_CnT"
bulk_npc_peak = GRanges(
  seqnames = bulk_npc_peak$V1,
  ranges = IRanges(
    start = bulk_npc_peak$V2,
    end = bulk_npc_peak$V3,
    names = bulk_npc_peak$type,
  )
)

cluster0_unsorted$type = "cluster_0"
cluster1_unsorted$type = "cluster_1"
bulk_npc_peak$type = "NPC_bulk_CnT"
gr = c(cluster0_unsorted, cluster1_unsorted, bulk_npc_peak)
grl = splitAsList(gr, gr$type)
grl = unique(grl)

res = makeVennDiagram(Peaks = grl, NameOfPeaks = names(grl))

venn_cnt2venn <- function(venn_cnt) {
  n <- which(colnames(venn_cnt) == "Counts") - 1
  SetNames = colnames(venn_cnt)[1:n]
  Weight = venn_cnt[, "Counts"]
  names(Weight) <- apply(venn_cnt[, 1:n], 1, paste, collapse = "")
  Venn(SetNames = SetNames, Weight = Weight)
}

pdf(
  file = glue("{result_folder}bulk_NPC_unsorted_cluster0_1.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 8
)
v <- venn_cnt2venn(res$vennCounts)
print(plot(v, doWeights = FALSE))
dev.off()

# Venn with 3 sets - bulk neuron and cluster0/1 of unsorted scCnT
cluster0_unsorted = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")
cluster0_unsorted$type = "cluster0"
cluster0_unsorted = GRanges(
  seqnames = cluster0_unsorted$V1,
  ranges = IRanges(
    start = cluster0_unsorted$V2,
    end = cluster0_unsorted$V3,
    names = cluster0_unsorted$type,
  )
)

cluster1_unsorted = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
cluster1_unsorted$type = "cluster1"
cluster1_unsorted = GRanges(
  seqnames = cluster1_unsorted$V1,
  ranges = IRanges(
    start = cluster1_unsorted$V2,
    end = cluster1_unsorted$V3,
    names = cluster1_unsorted$type,
  )
)

bulk_neuron_peak = fread("../data/bed/bulk_CnT_G4_neuron_mm10.bed")
bulk_neuron_peak$type = "neuron_bulk_CnT"
bulk_neuron_peak = GRanges(
  seqnames = bulk_neuron_peak$V1,
  ranges = IRanges(
    start = bulk_neuron_peak$V2,
    end = bulk_neuron_peak$V3,
    names = bulk_neuron_peak$type,
  )
)

cluster0_unsorted$type = "cluster_0"
cluster1_unsorted$type = "cluster_1"
bulk_neuron_peak$type = "neuron_bulk_CnT"
gr = c(cluster0_unsorted, cluster1_unsorted, bulk_neuron_peak)
grl = splitAsList(gr, gr$type)
grl = unique(grl)

res = makeVennDiagram(Peaks = grl, NameOfPeaks = names(grl))

venn_cnt2venn <- function(venn_cnt) {
  n <- which(colnames(venn_cnt) == "Counts") - 1
  SetNames = colnames(venn_cnt)[1:n]
  Weight = venn_cnt[, "Counts"]
  names(Weight) <- apply(venn_cnt[, 1:n], 1, paste, collapse = "")
  Venn(SetNames = SetNames, Weight = Weight)
}

pdf(
  file = glue("{result_folder}bulk_neuron_unsorted_cluster0_1.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 8
)
v <- venn_cnt2venn(res$vennCounts)
print(plot(v, doWeights = FALSE))
dev.off()

gr = c(cluster0_unsorted, cluster1_unsorted)
grl = splitAsList(gr, gr$type)
grl = unique(grl)

res = makeVennDiagram(Peaks = grl, NameOfPeaks = names(grl))

# cluster 1 & 0 overlap
venn_cnt2venn <- function(venn_cnt) {
  n <- which(colnames(venn_cnt) == "Counts") - 1
  SetNames = colnames(venn_cnt)[1:n]
  Weight = venn_cnt[, "Counts"]
  names(Weight) <- apply(venn_cnt[, 1:n], 1, paste, collapse = "")
  Venn(SetNames = SetNames, Weight = Weight)
}

## NPC - neuron - brain G4 scCut&Tag Jaccard analysis
# jaccard indices (package cotools) and GenometriCorrelation analysis (package GenometriCorr)
peaks = list(cluster0_unsorted, cluster1_unsorted, bulk_neuron_peak, bulk_npc_peak)
jaccards = matrix(NA_real_, length(peaks), length(peaks))
colnames(jaccards) = c("cluster 0", "cluster 1", "bulk neuron CnT", "bulk NPC CnT")
rownames(jaccards) = c("cluster 0", "cluster 1", "bulk neuron CnT", "bulk NPC CnT")
for(i in seq(1, ncol(jaccards))) {
  for(j in seq(1, nrow(jaccards))) {
    jaccard = genomicCorr.jaccard(peaks[[i]], peaks[[j]])
    jaccards[i, j] = jaccard
  }
}
is.matrix(jaccards)

pdf(
  file = glue("{result_folder}unsorted_res0.1_Jaccard_hm.pdf"),
  width = 6,
  height = 6
)
col_fun = colorRamp2(c(0, 0.25, 0.5), c("#9ecae1", "white", "#fc9272"))
Heatmap(
  jaccards,
  column_title = "",
  row_title = "",
  name = "Jaccard index",
  # row_km = 2,
  # column_km = 1,
  clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(8, "cm"),
  heatmap_height = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", jaccards[i, j]), x, y, gp = gpar(fontsize = 10))}
)
dev.off()

## unsorted scCut&Tag
## generate bed files from overlap analysis for deeptools heatmaps
bed_output = "../results/GenomicRanges/unsorted_outputs/bulk_overlaps/"
peak_set = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/"
bigwigs = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/"
#list.files(bigwigs)

# scCnT data
cluster0 = fread(glue("{peak_set}0_peaks.narrowPeak"))

cluster0_gr = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    signalValue = cluster0$V7
  )
)

cluster1 = fread(glue("{peak_set}1_peaks.narrowPeak"))

cluster1_gr = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    signalValue = cluster1$V7
  )
)

all = rbind(cluster0, cluster1)
all = all %>% mutate(type = "unsorted")

all = GRanges(
  seqnames = all$V1,
  ranges = IRanges(
    start = all$V2,
    end = all$V3,
    signalValue = all$V7
  )
)

# unique cluster 0
read_cov = bw_loci(glue(bigwigs, "0.bam_RPGC.bigwig"), glue(bigwigs, "1.bam_RPGC.bigwig"), 
                   loci = all, )
read_cov = as.data.frame(read_cov)
read_cov = inner_join(read_cov, as.data.frame(all), by = c("seqnames", "start", "end"))

only_cluster0 = read_cov %>% dplyr::filter(!X0.bam_RPGC == Inf) %>% 
  dplyr::filter(X0.bam_RPGC >= 4) %>% 
  dplyr::select(seqnames, start, end, signalValue)

write_tsv(only_cluster0, glue("{bed_output}unsorted_cluster0_unique_peaks_res0.1.tsv"))

only_cluster0 = only_cluster0 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster0,
  glue("{bed_output}MACS2_cl0_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# unique cluster 1
read_cov = bw_loci(glue(bigwigs, "1.bam_RPGC.bigwig"), glue(bigwigs, "0.bam_RPGC.bigwig") , loci = all)
read_cov = as.data.frame(read_cov)
read_cov = inner_join(read_cov, as.data.frame(all), by = c("seqnames", "start", "end"))

only_cluster1 = read_cov %>% dplyr::filter(!X1.bam_RPGC == Inf) %>% 
  dplyr::filter(X1.bam_RPGC >= 4) %>% 
  dplyr::select(seqnames, start, end, signalValue)

write_tsv(only_cluster1, glue("{bed_output}unsorted_cluster1_unique_peaks_res0.1.tsv"))

only_cluster1 = only_cluster1 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster1,
  glue("{bed_output}MACS2_cl1_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# both cluster 0,1
ol = findOverlaps(cluster0_gr, cluster1_gr, type = "any", ignore.strand = FALSE)

# intersection
both = cluster0_gr[queryHits(ol)]
both = data.frame(
  seqnames = seqnames(both),
  starts = start(both) - 1,
  ends = end(both)
)

write.table(
  both,
  glue("{bed_output}MACS2_intersection.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# with GenomicRanges
cluster0_unique = cluster0_gr[-queryHits(findOverlaps(cluster0_gr, cluster1_gr, type="any")),]
cluster0_unique = as.data.frame(cluster0_unique)
cluster0_unique = cluster0_unique[which
                                  (cluster0_unique$signalValue > median(cluster0_unique$signalValue)),]
cluster0_unique %>% dplyr::select(seqnames, start, end) %>% 
  write_tsv(., glue("{bed_output}unique_cluster0-median_filtered-GR.bed"), col_names = FALSE)

cluster1_unique = cluster1_gr[-queryHits(findOverlaps(cluster1_gr, cluster0_gr, type="any")),]
cluster1_unique = as.data.frame(cluster1_unique)
cluster1_unique = cluster1_unique[which
                                  (cluster1_unique$signalValue > median(cluster1_unique$signalValue)),]
cluster1_unique %>% dplyr::select(seqnames, start, end) %>% 
  write_tsv(., glue("{bed_output}unique_cluster1-median_filtered-GR.bed"), col_names = FALSE)


# unsorted scCut&Tag reso 0.8
peak_set = "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_peaks/"
result_folder = "../results/GenomicRanges/unsorted_outputs/"
# scCnT data
cluster0 = fread(glue("{peak_set}0_peaks.narrowPeak"))
cluster0$type = "cluster0"
cluster0_gr = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    names = cluster0$type,
  )
)
cluster0_gr$type = "cluster0"


cluster1 = fread(glue("{peak_set}1_peaks.narrowPeak"))
cluster1$type = "cluster1"
cluster1_gr = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    names = cluster1$type,
  )
)
cluster1_gr$type = "cluster1"

cluster2 = fread(glue("{peak_set}2_peaks.narrowPeak"))
cluster2$type = "cluster2"
cluster2_gr = GRanges(
  seqnames = cluster2$V1,
  ranges = IRanges(
    start = cluster2$V2,
    end = cluster2$V3,
    names = cluster2$type,
  )
)
cluster2_gr$type = "cluster2"

cluster3 = fread(glue("{peak_set}3_peaks.narrowPeak"))
cluster3$type = "cluster3"
cluster3_gr = GRanges(
  seqnames = cluster3$V1,
  ranges = IRanges(
    start = cluster3$V2,
    end = cluster3$V3,
    names = cluster3$type,
  )
)
cluster3_gr$type = "cluster3"

cluster4 = fread(glue("{peak_set}4_peaks.narrowPeak"))
cluster4_gr = GRanges(
  seqnames = cluster4$V1,
  ranges = IRanges(
    start = cluster4$V2,
    end = cluster4$V3,
    names = cluster4$type,
  )
)
cluster4_gr$type = "cluster4"

gr = c(cluster0_gr, cluster1_gr, cluster2_gr, cluster3_gr, cluster4_gr)
grl = splitAsList(gr, gr$type)
grl = unique(grl)

pdf(
  file = glue("{result_folder}unsorted_cluster0_1_2_3_4.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 8
)
res = makeVennDiagram(Peaks = grl, NameOfPeaks = names(grl))
print(res)
dev.off()

venn_cnt2venn <- function(venn_cnt) {
  n <- which(colnames(venn_cnt) == "Counts") - 1
  SetNames = colnames(venn_cnt)[1:n]
  Weight = venn_cnt[, "Counts"]
  names(Weight) <- apply(venn_cnt[, 1:n], 1, paste, collapse = "")
  Venn(SetNames = SetNames, Weight = Weight)
}


v <- venn_cnt2venn(res$vennCounts)
print(plot(v, doWeights = FALSE))


## mESC-MEF scCut&Tag
## generate bed files from overlap analysis for deeptools heatmaps
bed_output = "../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/"
peak_set = "../results/Seurat/final/mesc_mef/cluster_spec_peaks/"
bigwigs = "../results/Seurat/final/mesc_mef/cluster_spec_bigwigs/"
#list.files(bigwigs)

# scCnT data
cluster0 = fread(glue("{peak_set}0_peaks.narrowPeak"))

cluster0_gr = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    signalValue = cluster0$V7
  )
)

cluster1 = fread(glue("{peak_set}1_peaks.narrowPeak"))

cluster1_gr = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    signalValue = cluster1$V7
  )
)

all = rbind(cluster0, cluster1)
all = all %>% mutate(type = "unsorted")

all = GRanges(
  seqnames = all$V1,
  ranges = IRanges(
    start = all$V2,
    end = all$V3,
    signalValue = all$V7
  )
)

# unique cluster 0
read_cov = bw_loci(glue(bigwigs, "0.bam_RPGC.bigwig"), glue(bigwigs, "1.bam_RPGC.bigwig") , loci = all)
read_cov = as.data.frame(read_cov)
read_cov = inner_join(read_cov, as.data.frame(all), by = c("seqnames", "start", "end"))

only_cluster0 = read_cov %>% dplyr::filter(!X0.bam_RPGC == Inf) %>% 
  dplyr::filter(X0.bam_RPGC >= 4) %>% 
  dplyr::select(seqnames, start, end, signalValue)

write_tsv(only_cluster0, glue("{bed_output}unsorted_cluster0_unique_peaks_res0.1.tsv"))

only_cluster0 = only_cluster0 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster0,
  glue("{bed_output}MACS2_cl0_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# unique cluster 1
read_cov = bw_loci(glue(bigwigs, "1.bam_RPGC.bigwig"), glue(bigwigs, "0.bam_RPGC.bigwig") , loci = all)
read_cov = as.data.frame(read_cov)
read_cov = inner_join(read_cov, as.data.frame(all), by = c("seqnames", "start", "end"))

only_cluster1 = read_cov %>% dplyr::filter(!X1.bam_RPGC == Inf) %>% 
  dplyr::filter(X1.bam_RPGC >= 4) %>% 
  dplyr::select(seqnames, start, end, signalValue)

write_tsv(only_cluster1, glue("{bed_output}unsorted_cluster1_unique_peaks_res0.1.tsv"))

only_cluster1 = only_cluster1 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster1,
  glue("{bed_output}MACS2_cl1_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# both cluster 0,1
ol = findOverlaps(cluster0_gr, cluster1_gr, type = "any", ignore.strand = FALSE)

# intersection
both = cluster0_gr[queryHits(ol)]
both = data.frame(
  seqnames = seqnames(both),
  starts = start(both) - 1,
  ends = end(both)
)

write.table(
  both,
  glue("{bed_output}MACS2_intersection.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

