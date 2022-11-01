# packages
suppressPackageStartupMessages({
  library("rtracklayer")
  library("glue")
  library("data.table")
  library("tidyverse")
})

# liftover function
mm9_to_mm10_liftover = function(bw = "../data/bw/GSM5625015_NPC_G4_CnT_Rep2_Batch2.bw",
                                chrom_sizes = "../data/mm10.chrom.sizes.txt",
                                chain = "../data/bw/mm9ToMm10.over.chain",
                                export = "../data/bw/GSM5625015_NPC_G4_CnT_Rep2_Batch2_mm10.bw") {
  mm10_chrom = fread(chrom_sizes, header = FALSE)
  bigwig = import(bw)
  chain = import.chain(chain)
  lo = liftOver(bigwig, chain)
  lo = unlist(lo)
  genome(lo) = "mm10"
  sl = tibble(chrom = names(seqlengths(lo)))
  sl = sl %>% left_join(., mm10_chrom, by = c("chrom" = "V1"))
  seqlengths(lo) = sl$V2
  hits = findOverlaps(lo, drop.self = TRUE)
  lo = lo[-queryHits(hits)]
  export(object = lo, con = export)

}

mm9_to_mm10_liftover()











