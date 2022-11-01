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
  library("wigglescout")
})

bed_output = "../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/"
bigwigs = "../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/"
list.files(bigwigs)

# scCnT data
cluster0 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/0_lanceotron_peaks_res0.1.tsv")
cluster0 = cluster0 %>% filter(`Peak Score` > 0) %>% dplyr::select(Chr, Start, End) %>% mutate(type = "cluster 0")
cluster0 = GRanges(
  seqnames = cluster0$Chr,
  ranges = IRanges(
    start = cluster0$Start,
    end = cluster0$End,
    names = cluster0$type,
  )
)

cluster1 = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/1_lanceotron_peaks_res0.1.tsv")
cluster1 = cluster1 %>% filter(`Peak Score` > 0) %>% dplyr::select(Chr, Start, End) %>% mutate(type = "cluster 1")
cluster1 = GRanges(
  seqnames = cluster1$Chr,
  ranges = IRanges(
    start = cluster1$Start,
    end = cluster1$End,
    names = cluster1$type,
  )
)

read_cov = bw_loci(glue(bigwigs, "0_res0.1.bigwig"), glue(bigwigs, "1_res0.1.bigwig") , loci = cluster0)
read_cov = as.data.frame(read_cov)
only_cluster0 = read_cov %>% filter(!X0_res0.1 == Inf) %>% filter(X0_res0.1 >= 2)
only_cluster0 = only_cluster0 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster0,
  glue("{bed_output}Lanceotron_cl0_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

read_cov = bw_loci(glue(bigwigs, "1_res0.1.bigwig"), glue(bigwigs, "0_res0.1.bigwig"), loci = cluster1)
read_cov = as.data.frame(read_cov)
only_cluster1 = read_cov %>% filter(!X1_res0.1 == Inf) %>% filter(X1_res0.1 >= 2)
only_cluster1 = only_cluster1 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster1,
  glue("{bed_output}Lanceotron_cl1_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# bulk NPC CnT data
npc = fread("../data/bed/bulk_G4_CnT_NPC_rep1_R1.mLb.clN_peaks.broadPeak")
npc$type = "bulk_mESC"
npc = GRanges(
  seqnames = npc$V1,
  ranges = IRanges(
    start = npc$V2,
    end = npc$V3,
    names = npc$type,
  )
)

# bulk MEF
mef = fread("../data/bed/bulk_G4_CnT_MEF_rep1_R1.mLb.clN_peaks.broadPeak")
mef$type = "mESC-MEF_cluster0"
mef = GRanges(
  seqnames = mef$V1,
  ranges = IRanges(
    start = mef$V2,
    end = mef$V3,
    names = mef$type,
  )
)

# bulk mESC CnT data
mesc_cnt = fread("../data/bed/bulk_G4_CnT_mESC_rep1_R1.mLb.clN_peaks.broadPeak")
mesc_cnt$type = "bulk_mESC"
mesc_cnt = GRanges(
  seqnames = mesc_cnt$V1,
  ranges = IRanges(
    start = mesc_cnt$V2,
    end = mesc_cnt$V3,
    names = mesc_cnt$type,
  )
)

## overlap between cluster 0,1
ol = findOverlaps(cluster0, cluster1, type = "any", ignore.strand = FALSE)

# intersection
both = cluster0[queryHits(ol)]
both = data.frame(
  seqnames = seqnames(both),
  starts = start(both) - 1,
  ends = end(both)
)

write.table(
  both,
  glue("{bed_output}Lanceotron_intersection.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)





