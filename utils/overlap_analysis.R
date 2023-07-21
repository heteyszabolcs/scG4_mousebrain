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
})

# result / peak folder
peak_folder = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/"
bed_folder = "../data/bed/"
result_folder = "../results/GenomicRanges/unsorted_outputs/"

create_venn = function(peak1,
                       peak2,
                       type1,
                       type2,
                       filename) {
  # filtered MACS2 peak calls
  peak = fread(glue("{peak_folder}{peak1}"))
  peak$type = type1
  peak = GRanges(
    seqnames = peak$V1,
    ranges = IRanges(
      start = peak$V2,
      end = peak$V3,
      names = peak$type,
    )
  )
  
  bulk_peak = fread(glue("{bed_folder}{peak2}"))
  bulk_peak$type = type2
  bulk_peak = GRanges(
    seqnames = bulk_peak$V1,
    ranges = IRanges(
      start = bulk_peak$V2,
      end = bulk_peak$V3,
      names = bulk_peak$type,
    )
  )
  
  peak$type = type1
  bulk_peak$type = type2
  gr = c(peak, bulk_peak)
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
    file = glue("{result_folder}{filename}"),
    # The directory you want to save the file in
    width = 8,
    height = 8
  )
  v <- venn_cnt2venn(res$vennCounts)
  print(plot(v, doWeights = FALSE))
  dev.off()
  
}

# mESC-MEF scCut&Tag data
# cluster 0 and cluster 1 of Seurat clustering with resolution 0.1
create_venn(
  peak1 = "0_peaks_res0.1.narrowPeak",
  peak2 = "bulk_CnT_G4_mES_mm10.bed",
  type1 = "scCnT mESC-MEF",
  type2 = "bulk CnT mESC",
  filename = "0_peaks_res0.1_vs_mESC_CnT-Venn.pdf"
)

create_venn(
  peak1 = "1_peaks_res0.1.narrowPeak",
  peak2 = "bulk_CnT_G4_mES_mm10.bed",
  type1 = "scCnT mESC-MEF",
  type2 = "bulk CnT mESC",
  filename = "1_peaks_res0.1_vs_mESC_CnT-Venn.pdf"
)

# Venn with 3 sets - bulk mESC and cluster0/1 of mESC-MEF scCnT
cluster0 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/0_peaks_res0.1.narrowPeak")
cluster0$type = "cluster0"
cluster0 = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    names = cluster0$type,
  )
)

cluster1 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/1_peaks_res0.1.narrowPeak")
cluster1$type = "cluster1"
cluster1 = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    names = cluster1$type,
  )
)

bulk_mesc_peak = fread("../data/bed/bulk_CnT_G4_mES_mm10.bed")
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
cluster0 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/0_peaks_res0.1.narrowPeak")
cluster0$type = "cluster0"
cluster0 = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    names = cluster0$type,
  )
)

cluster1 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/1_peaks_res0.1.narrowPeak")
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
bulk_peak$type = "MEF_bulk_CnT"
gr = c(cluster0, cluster1, bulk_peak)
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
  file = glue("{result_folder}bulk_mESC_MEF-MEFcluster0_1.pdf"),
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
bulk_peak$type = "NPC_bulk_CnT"
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
bulk_peak$type = "neuron_bulk_CnT"
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


### create bed files for bulk mESC and sc mESC-MEF
mesc_mef = fread("../data/CellRanger/mES-mEF/peaks_noheader.bed")
mesc_mef$type = "mESC-MEF"
mesc_mef = GRanges(
  seqnames = mesc_mef$V1,
  ranges = IRanges(
    start = mesc_mef$V2,
    end = mesc_mef$V3,
    names = mesc_mef$type,
  )
)

mesc_cnt = fread("../data/bed/bulk_CnT_G4_mES_mm10.bed")
mesc_cnt$type = "bulk_mESC"
mesc_cnt = GRanges(
  seqnames = mesc_cnt$V1,
  ranges = IRanges(
    start = mesc_cnt$V2,
    end = mesc_cnt$V3,
    names = mesc_cnt$type,
  )
)

npc_cnt = fread("../data/bed/bulk_CnT_G4_NPC_mm10.bed")
npc_cnt$type = "bulk_NPC"
npc_cnt = GRanges(
  seqnames = npc_cnt$V1,
  ranges = IRanges(
    start = npc_cnt$V2,
    end = npc_cnt$V3,
    names = npc_cnt$type,
  )
)

# overlaps between bulk mESC CnT and cluster 0 of mESC-MEF scCnT
mesc_cnt = fread("../data/bed/bulk_CnT_G4_mES_mm10.bed")
mesc_cnt$type = "sc_mESC_MEF"
mesc_cnt = GRanges(
  seqnames = mesc_cnt$V1,
  ranges = IRanges(
    start = mesc_cnt$V2,
    end = mesc_cnt$V3,
    names = mesc_cnt$type,
  )
)

cluster0 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/0_peaks_res0.1.narrowPeak")
cluster0$type = "cluster0"
cluster0 = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    names = cluster0$type,
  )
)

ol = findOverlaps(mesc_cnt, cluster0, type = "any", ignore.strand = FALSE)

only_bulk = mesc_cnt[-queryHits(ol)]
only_bulk = data.frame(
  seqnames = seqnames(only_bulk),
  starts = start(only_bulk) - 1,
  ends = end(only_bulk)
)

write.table(
  only_bulk,
  "../data/bed/mESC_cl0_ol-bulk_mESC-only.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

only_sc = cluster0[-subjectHits(ol)]
only_sc = data.frame(
  seqnames = seqnames(only_sc),
  starts = start(only_sc) - 1,
  ends = end(only_sc)
)

write.table(
  only_sc,
  "../data/bed/mESC_cl0_ol-sc_mESC_MEF-only.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

both_bulk = mesc_cnt[queryHits(ol)]
both_bulk = data.frame(
  seqnames = seqnames(both_bulk),
  starts = start(both_bulk) - 1,
  ends = end(both_bulk)
)

write.table(
  both_bulk,
  "../data/bed/mESC_cl0_ol-bulk_mESC-intersect.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

both_sc = cluster0[subjectHits(ol)]
both_sc = data.frame(
  seqnames = seqnames(both_sc),
  starts = start(both_sc) - 1,
  ends = end(both_sc)
)

write.table(
  both_sc,
  "../data/bed/mESC_cl0_ol-sc_mESC_MEF-intersect.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

cluster1 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/1_peaks_res0.1.narrowPeak")
cluster1$type = "cluster1"
cluster1 = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    names = cluster1$type,
  )
)

only_sc = cluster0[-subjectHits(ol)]
ol_cl1_onlysc = findOverlaps(cluster1, only_sc, type = "any", ignore.strand = FALSE)
only_cl0 = cluster0[-subjectHits(ol)]
only_cl0 = data.frame(
  seqnames = seqnames(only_cl0),
  starts = start(only_cl0) - 1,
  ends = end(only_cl0)
)

# overlaps between bulk mESC CnT and mESC-MEF scCnT
ol = findOverlaps(mesc_cnt, mesc_mef, type = "any", ignore.strand = FALSE)

only_bulk = mesc_cnt[-queryHits(ol)]

only_bulk = data.frame(
  seqnames = seqnames(only_bulk),
  starts = start(only_bulk) - 1,
  ends = end(only_bulk)
)

write.table(
  only_bulk,
  "../data/bed/mESC_overlap-bulk_mESC-only.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

only_sc = mesc_mef[-subjectHits(ol)]

only_sc = data.frame(
  seqnames = seqnames(only_sc),
  starts = start(only_sc) - 1,
  ends = end(only_sc)
)

write.table(
  only_sc,
  "../data/bed/mESC_overlap-sc_mESC_MEF-only.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

ol_bulk = mesc_cnt[queryHits(ol)]

ol_bulk = data.frame(
  seqnames = seqnames(ol_bulk),
  starts = start(ol_bulk) - 1,
  ends = end(ol_bulk)
)

write.table(
  ol_bulk,
  "../data/bed/mESC_overlap-bulk_mESC-intersect.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

ol_sc = mesc_mef[subjectHits(ol)]

ol_sc = data.frame(
  seqnames = seqnames(ol_sc),
  starts = start(ol_sc) - 1,
  ends = end(ol_sc)
)

write.table(
  ol_sc,
  "../data/bed/mESC_overlap-sc_mESC_MEF-intersect.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# overlaps between bulk NPC CnT and mESC-MEF scCnT
ol = findOverlaps(npc_cnt, mesc_mef, type = "any", ignore.strand = FALSE)

only_bulk = npc_cnt[-queryHits(ol)]

only_bulk = data.frame(
  seqnames = seqnames(only_bulk),
  starts = start(only_bulk) - 1,
  ends = end(only_bulk)
)

write.table(
  only_bulk,
  "../data/bed/NPC_overlap-bulk_NPC-only.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

only_sc = mesc_mef[-subjectHits(ol)]

only_sc = data.frame(
  seqnames = seqnames(only_sc),
  starts = start(only_sc) - 1,
  ends = end(only_sc)
)

write.table(
  only_sc,
  "../data/bed/NPC_overlap-sc_mESC_MEF-only.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

ol_bulk = npc_cnt[queryHits(ol)]

ol_bulk = data.frame(
  seqnames = seqnames(ol_bulk),
  starts = start(ol_bulk) - 1,
  ends = end(ol_bulk)
)

write.table(
  ol_bulk,
  "../data/bed/NPC_overlap-bulk_NPC-intersect.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

ol_sc = mesc_mef[subjectHits(ol)]

ol_sc = data.frame(
  seqnames = seqnames(ol_sc),
  starts = start(ol_sc) - 1,
  ends = end(ol_sc)
)

write.table(
  ol_sc,
  "../data/bed/NPC_overlap-sc_mESC_MEF-intersect.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)
