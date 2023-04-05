library("wigglescout")
library("data.table")
library("tidyverse")

k4me3 = "../data/bw/GSM4303808_WT-H3K4me3-1_mm10.bigwig"
bed1 = "../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MACS2_intersection.bed"

values = bw_loci("../data/bw/GSM4303808_WT-H3K4me3-1_mm10.bigwig", loci = bed1)
values = as.data.frame(values)
values = values %>% arrange(desc(GSM4303808_WT.H3K4me3.1_mm10))

new_bed = values[,1:3]

write_tsv(new_bed, "../results/GenomicRanges/unsorted_outputs/bulk_overlaps/MMACS2_intersection_K4me3_signal_sorted.bed", 
          col_names = FALSE)
                            