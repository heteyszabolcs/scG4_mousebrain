suppressPackageStartupMessages({
  library("tidyverse")
  library("Signac")
  library("data.table")
  library("GenomicRanges")
  library("wigglescout")
})

# make coverage plot 
# object: ChromatinAssay object made by Signac (e.g. scATAC-Seq, scCut&Tag objects)
# gene: region of interest
# bulk: visualize pseudobulk 
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

# G4 scCut&Tag data
# Seurat objects
sorted = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")
unsorted = readRDS(file = "../results/Seurat/callpeaks_unsorted/unsorted.Rds")
mesc_mef = readRDS(file = "../results/Seurat/callpeaks_mESC-MEF/mESC-MEF_res0.1.Rds")

# scRNA-Seq oligodendrocyte markers (Marques et al.)
markers = fread("../data/GSE75330/marker_genes.txt", header = FALSE)
opc = markers %>% dplyr::filter(V2 == "OPC") %>% pull(V1)
opc = opc[which(!opc %in% c("Igfb2", "Lhfpl3"))]

# runs
make_coverageplot(object = sorted, gene = "Cspg5", group = "seurat_clusters", bulk = TRUE)


sapply(opc, make_coverageplot, object = sorted, group = "seurat_clusters", bulk = TRUE)


make_coverageplot(object = unsorted, gene = "Cspg5", group = "seurat_clusters", bulk = TRUE)

make_coverageplot(object = unsorted, gene = "Myc", group = "seurat_clusters", bulk = TRUE)
ggsave(
  "../results/Seurat/callpeaks_unsorted/coverageplot_Myc_common.pdf",
  plot = last_plot(),
  width = 6,
  height = 6
)



## TilePlots
# export
ggsave(
  "../results/genome_browser/GFPsorted-marcks.pdf",
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)
 
TilePlot(object = mesc_mef, region = c("Cdhr3"), tile.cells = 100, order.by = "total") + 
  scale_fill_gradient(low = "white", high = "#3182bd")

# export
ggsave(
  "../results/Seurat/callpeaks_mESC-MEF/TilePlot_Cdhr3_cluster0_spec.pdf",
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)

TilePlot(object = unsorted, region = "Th", tile.cells = 100, order.by = "total") + 
  scale_fill_gradient(low = "white", high = "#3182bd")
# export
ggsave(
  "../results/Seurat/callpeaks_unsorted/TilePlot_Th_common.pdf",
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)

TilePlot(object = unsorted, region = "Myc", tile.cells = 100, order.by = "total") + 
  scale_fill_gradient(low = "white", high = "#3182bd")
# export
ggsave(
  "../results/Seurat/callpeaks_unsorted/TilePlot_Myc_common.pdf",
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)

TilePlot(object = mesc_mef, region = c("Lin28a"), tile.cells = 100, order.by = "total") + scale_fill_gradient(low = "white", high = "#3182bd")
# export
ggsave(
  "../results/Seurat/callpeaks_mESC-MEF/TilePlot_Lin28a_cluster1_spec.pdf",
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)


