# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("ggpubr")
  library("cowplot")
  library("Seurat")
  library("wigglescout")
  library("GenomicRanges")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
})

# get unique peaks by wigglescout
get_unique_ws = function(bw,
                         bw_backgr,
                         subset,
                         thr = 4) {
  label = "fold_change"
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  
  read_cov = bw_loci(
    bwfiles = bw,
    bg_bwfiles = bw_backgr,
    labels = label,
    loci = subset
  )
  
  annot = annotatePeak(
    read_cov,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db"
  )
  annot = as.data.frame(annot)
  annot = annot %>% dplyr::select(seqnames,
                                  start,
                                  end,
                                  starts_with("fold_change"),
                                  distanceToTSS,
                                  gene_symbol = SYMBOL) %>% dplyr::filter(abs(fold_change) > thr) %>%
    dplyr::filter(!fold_change == Inf)
  
  return(annot)
  
}

# MEF-mESC scCnT
peaks = fread("../results/Seurat/callpeaks_mESC-MEF/peak_sets/peaks_per_clusters_res0.1.bed")
peaks = GRanges(
  seqnames = peaks$V1,
  ranges = IRanges(
    start = peaks$V2,
    end = peaks$V3
  )
)

cl1 = get_unique_ws(bw = "../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/1_res0.1.bigwig", 
                      bw_backgr = "../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/0_res0.1.bigwig",
                      subset = peaks)
cl0 = get_unique_ws(bw = "../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/0_res0.1.bigwig", 
                    bw_backgr = "../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/1_res0.1.bigwig",
                    subset = peaks)

# highly predicted MOL cells
peaks = fread("../results/Seurat/callpeaks_GFPsorted/high_pred_MOL_peaks.bed")
peaks = GRanges(
  seqnames = peaks$V1,
  ranges = IRanges(
    start = peaks$V2,
    end = peaks$V3
  )
)

mol = get_unique_ws(bw = "../results/Seurat/callpeaks_GFPsorted/high_pred_MOL.bw", 
                    bw_backgr = "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/1_res0.1_brain_cells.bw",
                    subset = peaks)

mol_prom = mol %>% dplyr::filter(abs(distanceToTSS) < 3000)





