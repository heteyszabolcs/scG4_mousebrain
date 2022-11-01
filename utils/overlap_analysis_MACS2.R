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

bed_output = "../results/GenomicRanges/unsorted_outputs/bulk_overlaps/"
peak_set = "../results/Seurat/callpeaks_unsorted/peak_sets/"
bigwigs = "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/"
list.files(bigwigs)

# scCnT data
cluster0 = fread("../results/Seurat/callpeaks_unsorted/peak_sets/0_peaks_res0.1.narrowPeak")

cluster0_gr = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    signalValue = cluster0$V7
  )
)

cluster1 = fread("../results/Seurat/callpeaks_unsorted/peak_sets/1_peaks_res0.1.narrowPeak")

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


read_cov = bw_loci(glue(bigwigs, "0_res0.1_oligo.bw"), glue(bigwigs, "1_res0.1_brain_cells.bw") , loci = all)
read_cov = as.data.frame(read_cov)
read_cov = inner_join(read_cov, as.data.frame(all), by = c("seqnames", "start", "end"))

only_cluster0 = read_cov %>% dplyr::filter(!X0_res0.1_oligo == Inf) %>% 
  dplyr::filter(X0_res0.1_oligo >= 4) %>% 
  dplyr::select(seqnames, start, end, signalValue)

write_tsv(only_cluster0, glue("{peak_set}0_unique_peaks_res0.1.tsv"))

only_cluster0 = only_cluster0 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster0,
  glue("{bed_output}MACS2_cl0_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

read_cov = bw_loci(glue(bigwigs, "1_res0.1_brain_cells.bw"), glue(bigwigs, "0_res0.1_oligo.bw") , loci = all)
read_cov = as.data.frame(read_cov)
read_cov = inner_join(read_cov, as.data.frame(all), by = c("seqnames", "start", "end"))

only_cluster1 = read_cov %>% dplyr::filter(!X1_res0.1_brain_cells == Inf) %>% 
  dplyr::filter(X1_res0.1_brain_cells >= 4) %>% 
  dplyr::select(seqnames, start, end, signalValue)

write_tsv(only_cluster1, glue("{peak_set}1_unique_peaks_res0.1.tsv"))

only_cluster1 = only_cluster1 %>% dplyr::select(seqnames, start, end)

write.table(
  only_cluster1,
  glue("{bed_output}MACS2_cl1_only.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

## overlap between cluster 0,1
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





