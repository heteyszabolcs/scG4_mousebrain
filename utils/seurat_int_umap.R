library("Seurat")
library("Signac")
library("tidyverse")
library("data.table")
library("glue")
library("devtools")

# list of Seurat integrated objects
integrated_rds = list.files("../results/Seurat/integration/GA/", pattern = "_integrated.rds", full.names = TRUE)

# create a UMAP plot about the integrated data grouped by cell types (annotated by Marek Bartosovic)
umap_group_by_celltypes = function(integrated_seurat = "../results/Seurat/integration/GA/H3K4me3_G4_integrated.rds") {
  
  label = strsplit(integrated_seurat, 
                   "../results/Seurat/integration/GA/")[[1]][2]
  label = strsplit(label, ".rds")[[1]][1]
  
  integrated_seurat = readRDS(integrated_seurat)
  DefaultAssay(integrated_seurat) = "GA"
  integrated_seurat = ScaleData(integrated_seurat, verbose = FALSE)
  integrated_seurat = FindVariableFeatures(object = integrated_seurat)
  integrated_seurat = RunPCA(integrated_seurat, npcs = 30, verbose = FALSE, 
                             features = VariableFeatures(object = integrated_seurat))
  integrated_seurat = RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)
  integrated_seurat = FindNeighbors(integrated_seurat, reduction = "pca", dims = 1:30)
  integrated_seurat = FindClusters(integrated_seurat, resolution = 0.5)
  
  umap = DimPlot(integrated_seurat, reduction = "umap", group.by = "cell_type")
  umap + ggtitle(label) 
  
  ggsave(
    glue("../results/Seurat/Seurat_int_umap-{label}.png"),
    plot = umap,
    width = 10,
    height = 10,
    dpi = 300,
  )
  
}

lapply(integrated_rds, umap_group_by_celltypes)



