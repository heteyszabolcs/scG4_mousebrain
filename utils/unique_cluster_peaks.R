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
  library("ArchR")
  library("GenomicRanges")
  library("cowplot")
})

# result / peak folder
peak_folder = "../results/Seurat/callpeaks_GFPsorted/peak_sets/"

peak_list = list.files(glue("{peak_folder}"), pattern = "^[0-9]_peaks.bed$")

# generate genomicrange object
create_gr = function(bed, name) {
  bed = fread(glue("{peak_folder}{bed}"))
  rownumber = nrow(bed)
  gr = GRanges(seqnames = bed$V1,
               ranges = IRanges(
                 start = bed$V2,
                 end = bed$V3,
                 names = rep(name, rownumber)
               ))
  
  return(gr)
}

# Signac MACS2 peaks
peak_list
peak0 = create_gr(bed = "0_peaks.bed", name = "0")
peak1 = create_gr(bed = "1_peaks.bed", name = "1")
peak2 = create_gr(bed = "2_peaks.bed", name = "2")
peak3 = create_gr(bed = "3_peaks.bed", name = "3")
peak4 = create_gr(bed = "4_peaks.bed", name = "4")

# create Venn diagram - overlap between clusters
peak0$type = "cluster 0"
peak1$type = "cluster 1"
peak2$type = "cluster 2"
peak3$type = "cluster 3"
peak4$type = "cluster 4"

gr = c(peak0, peak1, peak2, peak3, peak4)
grl = splitAsList(gr, gr$type)
#grl = unique(grl)
res = makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl))

venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Venn(SetNames=SetNames, Weight=Weight)
}

pdf(
  file = glue("{peak_folder}peak_calling_per_cluster_VennPlot.pdf"),
  # The directory you want to save the file in
  width = 8,
  height = 8
)
v <- venn_cnt2venn(res$vennCounts)
plot(v, doWeights = TRUE)
dev.off()

# cluster 0
ol_1 = peak0[queryHits(findOverlaps(
  peak0,
  peak1,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)),]
ol_2 = peak0[queryHits(findOverlaps(
  peak0,
  peak2,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)),]
ol_3 = peak0[queryHits(findOverlaps(
  peak0,
  peak3,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)),]
ol_4 = peak0[queryHits(findOverlaps(
  peak0,
  peak4,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)),]
ol = c(ol_1, ol_2, ol_3, ol_4)

unique_0 = GenomicRanges::setdiff(peak0, ol)
unique_0 = as_tibble(unique_0)
bed0 = unique_0[, 1:3]
write_tsv(bed0,
          glue("{peak_folder}unique_0_peaks_lanceotron.bed"),
          col_names = FALSE)

# cluster 1
ol_1 = peak1[queryHits(findOverlaps(
  peak1,
  peak0,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_2 = peak1[queryHits(findOverlaps(
  peak1,
  peak2,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_3 = peak1[queryHits(findOverlaps(
  peak1,
  peak3,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_4 = peak1[queryHits(findOverlaps(
  peak1,
  peak4,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol = c(ol_1, ol_2, ol_3, ol_4)

unique_1 = GenomicRanges::setdiff(peak1, ol)
unique_1 = as_tibble(unique_1)
bed1 = unique_1[, 1:3]
write_tsv(bed1,
          glue("{peak_folder}unique_1_peaks_lanceotron.bed"),
          col_names = FALSE)

# cluster 2
ol_1 = peak2[queryHits(findOverlaps(
  peak2,
  peak1,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_2 = peak2[queryHits(findOverlaps(
  peak2,
  peak0,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_3 = peak2[queryHits(findOverlaps(
  peak2,
  peak3,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_4 = peak2[queryHits(findOverlaps(
  peak2,
  peak4,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol = c(ol_1, ol_2, ol_3, ol_4)

unique_2 = GenomicRanges::setdiff(peak2, ol)
unique_2 = as_tibble(unique_2)
bed2 = unique_2[, 1:3]
write_tsv(bed2,
          glue("{peak_folder}unique_2_peaks_lanceotron.bed"),
          col_names = FALSE)

# cluster 3
ol_1 = peak3[queryHits(findOverlaps(
  peak3,
  peak1,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_2 = peak3[queryHits(findOverlaps(
  peak3,
  peak2,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_3 = peak3[queryHits(findOverlaps(
  peak3,
  peak0,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_4 = peak3[queryHits(findOverlaps(
  peak3,
  peak4,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol = c(ol_1, ol_2, ol_3, ol_4)

unique_3 = GenomicRanges::setdiff(peak3, ol)
unique_3 = as_tibble(unique_3)
bed3 = unique_3[, 1:3]
write_tsv(bed3,
          glue("{peak_folder}unique_3_peaks_lanceotron.bed"),
          col_names = FALSE)

# cluster 4
ol_1 = peak4[queryHits(findOverlaps(
  peak4,
  peak1,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_2 = peak4[queryHits(findOverlaps(
  peak4,
  peak2,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_3 = peak4[queryHits(findOverlaps(
  peak4,
  peak3,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol_4 = peak4[queryHits(findOverlaps(
  peak4,
  peak0,
  maxgap = -1L,
  minoverlap = 1,
  type = c("any")
)), ]
ol = c(ol_1, ol_2, ol_3, ol_4)

unique_4 = GenomicRanges::setdiff(peak4, ol)
unique_4 = as_tibble(unique_4)
bed4 = unique_4[, 1:3]
write_tsv(bed4,
          glue("{peak_folder}unique_4_peaks_lanceotron.bed"),
          col_names = FALSE)
