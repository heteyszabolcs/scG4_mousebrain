library("Seurat")
library("Signac")
library("tidyverse")
library("data.table")
library("glue")
library("devtools")

# path to result folder
result_path = "../results/Seurat/"

# function for Seurat integration
# source: https://satijalab.org/seurat/articles/integration_introduction.html
integration = function(seurat1 = "../data/GSE157637/Olig2_seurat_object.Rds",
                       seurat2 = "../data/merged/Seurat_merged.Rds",
                       label1 = "Olig2",
                       label2 = "G4",
                       assay = "bins_5000",
                       export_rds = TRUE,
                       export_image = TRUE) {
  print(glue("{label1} with {label2}"))
  
  seurat1 = readRDS(seurat1)
  DefaultAssay(seurat1) <- assay
  seurat1 = AddMetaData(object = seurat1,
                        metadata = label1,
                        col.name = "group")
  for (i in Assays(seurat1)) {
    if (i != assay) {
      seurat1[[i]] = NULL
    }
  }
  
  seurat2 = readRDS(seurat2)
  DefaultAssay(seurat2) <- assay
  seurat2 = AddMetaData(object = seurat2,
                        metadata = label2,
                        col.name = "group")
  for (i in Assays(seurat2)) {
    if (i != assay) {
      seurat2[[i]] = NULL
    }
  }
  
  comb = merge(
    seurat1,
    y = seurat2,
    add.cell.ids = c(label1, label2),
    project = "int"
  )
  
  rm(seurat1)
  rm(seurat2)
  
  # split combined object keeping atomic feature
  split = SplitObject(comb, split.by = "group")
  
  # normalize and identify variable features for each dataset independently
  split <- lapply(
    X = split,
    FUN = function(x) {
      x <- NormalizeData(x)
      x <-
        FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    }
  )
  
  # select features that are repeatedly variable across datasets for integration
  features = SelectIntegrationFeatures(object.list = split)
  # find anchors
  anchors = FindIntegrationAnchors(object.list = split, anchor.features = features)
  # this command creates an 'integrated' data assay
  integrated = IntegrateData(anchorset = anchors)
  # export rds
  if (export_rds) {
    saveRDS(integrated,
            glue("{result_path}{label1}_{label2}_integrated.rds"))
  }
  
  # integration basted on Stuart et al. (Seurat protocol)
  # run the standard workflow for visualization and clustering
  standard_wf <- ScaleData(integrated, verbose = FALSE)
  standard_wf <- RunPCA(standard_wf, npcs = 30, verbose = FALSE)
  standard_wf <-
    RunUMAP(standard_wf, reduction = "pca", dims = 1:30)
  standard_wf <-
    FindNeighbors(standard_wf, reduction = "pca", dims = 1:30)
  standard_wf <- FindClusters(standard_wf, resolution = 0.5)
  
  # Visualization
  p1 <-
    DimPlot(standard_wf, reduction = "umap", group.by = "group")
  p1 + ggtitle(glue("scC&T integration, {label1} - {label2}")) + scale_color_manual(labels = c(label2, glue("Bartosovic et al. {label1}")),
                                                                                    values = c("#3182bd", "#e34a33"))
  
  if (export_image) {
    ggsave(
      glue("{result_path}Seurat_int-{label1}_{label2}.png"),
      plot = last_plot(),
      width = 10,
      height = 10,
      dpi = 300,
    )
  }
  
  return(print(p1))
}

# run
integration(seurat1 = "../data/GSE157637/H3K27ac_seurat_object.Rds",
            label1 = "H3K27ac",
            assay = "GA")
integration(seurat1 = "../data/GSE157637/H3K4me3_seurat_object.Rds",
            label1 = "H3K4me3",
            assay = "GA")
integration(seurat1 = "../data/GSE157637/H3K27me3_seurat_object.Rds",
            label1 = "H3K27me3",
            assay = "GA")
integration(seurat1 = "../data/GSE157637/Rad21_seurat_object.Rds",
            label1 = "Rad21",
            assay = "GA")
integration(seurat1 = "../data/GSE157637/Olig2_seurat_object.Rds",
            label1 = "Olig2",
            assay = "GA")
integration(seurat1 = "../data/GSE157637/H3K36me3_seurat_object.Rds",
            label1 = "H3K36me3",
            assay = "GA")
