library("Seurat")
library("Signac")
library("data.table")
library("glue")
library("tidyverse")

# path to result folder
result_path = "../results/Seurat/"

featureplot = function(ptm_data,
                       marker_path = "../data/markers_brain_nbiotech/",
                       cell_type = "Astrocytes") {
  # label
  label = strsplit(ptm_data, "../data/GSE157637/")[[1]][2]
  label = strsplit(label, "_seurat_object.Rds")[[1]][1]
  
  markers = list.files(marker_path, full.names = TRUE, pattern = glue("{label}*.txt"))
  
  if(length(markers) == 0) {
    return(print(glue("No marker found for {label}")))
  }
  
  
  # single-cell objects
  seurat1 = readRDS(ptm_data)
  DefaultAssay(seurat1) = "GA"
  seurat2 = readRDS("../data/merged/Seurat_merged.Rds")
  DefaultAssay(seurat2) = "GA"
  
  # markers
  marker = fread(markers, header = TRUE)
  marker = marker %>%
    filter(cluster == cell_type) %>%
    arrange(desc(avg_logFC))
  
  features = marker$closest_gene[1:8]
  
  dimplot1 = DimPlot(seurat1, reduction = "umap", label = TRUE)
  dimplot1 = dimplot1 + ggtitle(glue("{label} scCut&Tag (Bartosovic et al.)"))
  
  ggsave(
    glue("{result_path}dimplot_{label}_scCut&Tag.png"),
    plot = dimplot1,
    width = 10,
    height = 10,
    dpi = 300,
  )

  featureplot1 = FeaturePlot(seurat1, features = features, reduction = "umap")
  
  ggsave(
    glue("{result_path}featureplot_{label}_scCut&Tag.png"),
    plot = featureplot1,
    width = 10,
    height = 10,
    dpi = 300,
  )
  
  dimplot2 = DimPlot(seurat2, reduction = "umap", label = TRUE)
  dimplot2 = dimplot2 + ggtitle("G4 scCut&Tag")
  dimplot2
  
  ggsave(
    glue("{result_path}dimplot_G4_scCut&Tag.png"),
    plot = dimplot2,
    width = 10,
    height = 10,
    dpi = 300,
  )
  
  featureplot2 = FeaturePlot(seurat2, features = features, reduction = "umap")
  
  ggsave(
    glue("{result_path}featureplot_{label}_G4_scCut&Tag.png"),
    plot = featureplot2,
    width = 10,
    height = 10,
    dpi = 300,
  )
  
  return(print(glue("{label} is done.")))
  
}

marek_data = list.files("../data/GSE157637/", pattern = "_seurat_object.Rds", full.names = TRUE)
lapply(marek_data, featureplot)

















