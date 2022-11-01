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
  library("GenometriCorr")
})


bed_output = "../results/GenomicRanges/mESC-MEF_outputs/bulk_overlaps/"

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

# mESC bulk 
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

## overlap between mESC CnT and MEF CnT
ol = findOverlaps(mesc_cnt, mef, type = "any", ignore.strand = FALSE)

# peaks only in bulk mESC
only_mesc = mesc_cnt[-queryHits(ol)]
only_mesc = data.frame(
  seqnames = seqnames(only_mesc),
  starts = start(only_mesc) - 1,
  ends = end(only_mesc)
)

write.table(
  only_mesc,
  glue("{bed_output}only_mESC_bulk.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# peaks only in bulk MEF
only_mef = mef[-subjectHits(ol)]
only_mef = data.frame(
  seqnames = seqnames(only_mef),
  starts = start(only_mef) - 1,
  ends = end(only_mef)
)

write.table(
  only_mef,
  glue("{bed_output}only_MEF_bulk.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# intersection
both_bulk = mesc_cnt[queryHits(ol)]
both_bulk = data.frame(
  seqnames = seqnames(both_bulk),
  starts = start(both_bulk) - 1,
  ends = end(both_bulk)
)

write.table(
  both_bulk,
  glue("{bed_output}intersection.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

## overlap between mESC CnT and NPC CnT
bed_output = "../results/GenomicRanges/unsorted_outputs/bulk_overlaps/"

ol = findOverlaps(mesc_cnt, npc, type = "any", ignore.strand = FALSE)

# peaks only in bulk mESC
only_mesc = mesc_cnt[-queryHits(ol)]
only_mesc = data.frame(
  seqnames = seqnames(only_mesc),
  starts = start(only_mesc) - 1,
  ends = end(only_mesc)
)

write.table(
  only_mesc,
  glue("{bed_output}only_mESC_bulk.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# peaks only in bulk MEF
only_npc = npc[-subjectHits(ol)]
only_npc = data.frame(
  seqnames = seqnames(only_npc),
  starts = start(only_npc) - 1,
  ends = end(only_npc)
)

write.table(
  only_npc,
  glue("{bed_output}only_NPC_bulk.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# intersection
both_bulk = mesc_cnt[queryHits(ol)]
both_bulk = data.frame(
  seqnames = seqnames(both_bulk),
  starts = start(both_bulk) - 1,
  ends = end(both_bulk)
)

write.table(
  both_bulk,
  glue("{bed_output}intersection.bed"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)









