suppressPackageStartupMessages({
  library("tidyverse")
  library("Signac")
  library("data.table")
  library("GenomicRanges")
  library("wigglescout")
})


make_coverageplot = function(object, gene, group, bulk = FALSE) {
  print(gene)
  
  if (!bulk) {
    plot = CoveragePlot(
      object = object,
      region = gene,
      show.bulk = FALSE,
      annotation = TRUE,
      group.by = group,
      peaks = TRUE
    ) + scale_fill_brewer(type = "seq", palette = "Set3")
    return(print(plot))
  } else {
    plot = CoveragePlot(
      object = object,
      region = gene,
      show.bulk = TRUE,
      annotation = TRUE,
      group.by = group,
      peaks = TRUE
    ) + scale_fill_brewer(type = "seq", palette = "Set3")
    return(print(plot))
  }
}

# G4 data
# Seurat objects
sorted = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")

make_coverageplot(object = sorted, gene = "Cspg5", group = "seurat_clusters", bulk = TRUE)


markers = fread("../data/GSE75330/marker_genes.txt", header = FALSE)
opc = markers %>% filter(V2 == "OPC") %>% pull(V1)
opc = opc[which(!opc %in% c("Igfb2", "Lhfpl3"))]

sapply(opc, make_coverageplot, object = sorted, group = "seurat_clusters", bulk = TRUE)

unsorted = readRDS(file = "../results/Seurat/callpeaks_unsorted/unsorted.Rds")
make_coverageplot(object = unsorted, gene = "chr13-5694290-5700626", group = "seurat_clusters", bulk = TRUE)

# export
ggsave(
  "../results/genome_browser/GFPsorted-marcks.pdf",
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)






