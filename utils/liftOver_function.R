library("rtracklayer")
library("glue")

help("liftOver")

# liftover function
mm9_to_mm10_liftover = function(bw = "../data/bw/ATAC-H33WT.mm9.bw",
                                chrom_sizes = "../data/mm10.chrom.sizes.txt",
                                chain = "../data/bw/mm9ToMm10.over.chain",
                                export = "../data/bw/ATAC-H33WT.mm10.bw") {
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











