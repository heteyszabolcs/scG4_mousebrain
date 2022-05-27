library("devtools")
library("GenomicRanges")
library("genomation")
library("gplots")

# result folder
result_folder = "../results/GenomicRanges/"

# ChromHMM bed, peak bed
chromhmm = "../data/bed/ENCFF850ZVL_mm10_chromhmm_brain.bed"
peaks = "../data/bed/consensus_G4s.bed"
chrHMM = readBed(chromhmm)
chrHMM.list = GenomicRanges::split(chrHMM, chrHMM$name,drop=TRUE)

peaks = fread(peaks)
peaks$V5 = "consensus_G4_peak"
scores = peaks$V4
peaks = GRanges(
  seqnames = peaks$V1,
  ranges = IRanges(
    start = peaks$V2,
    end = peaks$V3,
    names = archr$V5
  )
)
values(peaks) <- DataFrame(score = scores)

# genomation annotation
peak2ann = annotateWithFeatures(peaks, chrHMM.list)
peak.list=GenomicRanges::GRangesList("cons. G4 peaks" = peaks)
peak2ann.l=annotateWithFeatures(peak.list,chrHMM.list)

# export
pdf(file = glue("{result_folder}genomation_ChromHMM.pdf"),
    width = 8, 
    height = 2)
heatTargetAnnotation(peak2ann.l, plot = TRUE, col = c("white", "#2ca25f"))
dev.off()

# mat = heatTargetAnnotation(peak2ann.l, plot=FALSE)
# library(gplots)
# heatmap.2(mat,col="topo.colors",cellnote=round(mat),
#           notecol="black",trace="none",cexCol=0.5,cexRow=0.8)

