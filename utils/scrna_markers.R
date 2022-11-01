suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("data.table")
  library("glue")
  library("tidyverse")
  library("glue")
})

# result / marker list (Marques et al., 2017)
result_folder = "../results/Seurat/"
markers = fread("../data/GSE75330/marker_genes.txt", header = FALSE)

# G4 Seurat object
g4 = readRDS("../results/Seurat/callpeaks_unsorted/unsorted.Rds")
DefaultAssay(g4) = "GA"

# featurePlot visualization. How do G4 levels of oligodendrocyte scRNA-Seq marker genes cluster in on the UMAP plot? 
g4_fplot = function(markers) {
  FeaturePlot(
    g4,
    features = markers,
    pt.size = 1,
    blend = FALSE,
    cols = c("lightgrey", "#31a354")
  )
}

# oligodendrocyte markers
cell_types = unique(markers$V2)

for (cell_type in cell_types) {
  print(cell_type)
  list()
  
  filt = markers %>% filter(V2 == cell_type)
  filt = filt$V1
  
  if (length(filt) <= 16) {
    featureplot = g4_fplot(markers = filt)
    ggsave(
      glue(
        "{result_folder}scRNA_markers-{cell_type}_G4_featureplot.png"
      ),
      plot = featureplot,
      width = 10,
      height = 10,
      dpi = 500,
    )
  } else {
    filt1 = filt[1:16]
    featureplot = g4_fplot(markers = filt1)
    ggsave(
      glue(
        "{result_folder}scRNA_markers-{cell_type}_G4_featureplot_1.png"
      ),
      plot = featureplot,
      width = 10,
      height = 10,
      dpi = 500,
    )
    filt2 = filt[17:length(filt)]
    featureplot = g4_fplot(markers = filt2)
    ggsave(
      glue(
        "{result_folder}scRNA_markers-{cell_type}_G4_featureplot_2.png"
      ),
      plot = featureplot,
      width = 10,
      height = 10,
      dpi = 500,
    )
  }
}
