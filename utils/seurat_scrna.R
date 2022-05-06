# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("GenomicFeatures")
  library("Matrix")
  library("cowplot")
})

# scRNA-Seq folder (Marek's postnatal mouse brain scRNA-Seq data)
cellranger_folder = "../data/GSE163484/"
# path to result folder
result_folder = "../results/Seurat/"

### replication 1 - GSM4979874 ### 
# read and process data
counts1 = Read10X_h5(
  filename = glue(
    "{cellranger_folder}scRNASeq_GSM4979874_filtered_feature_bc_matrix_rep1.h5"
  )
)
rna1 = CreateSeuratObject(
  counts = counts1,
  project = "mouse_brain",
  min.cells = 3,
  min.features = 200
)

rna1[["percent.mt"]] <- PercentageFeatureSet(rna1, pattern = "^MT-")
VlnPlot(
  rna1,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
rna1 = FindVariableFeatures(rna1, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(rna1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(rna1)
plot2 <- LabelPoints(plot = plot1,
                     points = top10,
                     repel = TRUE)
varfeatures = plot1 + plot2

ggsave(
  glue("{result_folder}VariableFeatures_scRNASeq_rep1_Marek.png"),
  plot = varfeatures,
  width = 14,
  height = 8,
  dpi = 300,
)

all.genes = rownames(rna1)
rna1 = ScaleData(rna1, features = all.genes)
rna1[["replicate"]] <- "replicate 1"


# Perform linear dimensional reduction
rna1 = RunPCA(rna1, features = VariableFeatures(object = rna1))
rna1 = FindNeighbors(rna1, dims = 1:10)
rna1 = FindClusters(rna1, resolution = 0.5)
# Run non-linear dimensional reduction
rna1 = RunUMAP(
  rna1,
  dims = 1:10,
  n.neighbors = 40,
  n.components = 2
)
umap1 = DimPlot(rna1, reduction = "umap")
umap1


ggsave(
  glue("{result_folder}UMAP_nn40_ncomp2_scRNASeq_rep1_Marek.png"),
  plot = umap1,
  width = 10,
  height = 10,
  dpi = 300,
)

### replication 2 - GSM4979875 ###
# read and process data
counts2 = Read10X_h5(
  filename = glue(
    "{cellranger_folder}scRNASeq_GSM4979875_filtered_feature_bc_matrix_rep2.h5"
  )
)
rna2 = CreateSeuratObject(
  counts = counts2,
  project = "mouse_brain",
  min.cells = 3,
  min.features = 200
)

rna2[["percent.mt"]] <- PercentageFeatureSet(rna2, pattern = "^MT-")
VlnPlot(
  rna2,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
rna2 = FindVariableFeatures(rna2, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(rna2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(rna2)
plot2 <- LabelPoints(plot = plot1,
                     points = top10,
                     repel = TRUE)
varfeatures = plot1 + plot2

ggsave(
  glue("{result_folder}VariableFeatures_scRNASeq_rep2_Marek.png"),
  plot = varfeatures,
  width = 14,
  height = 8,
  dpi = 300,
)

all.genes = rownames(rna2)
rna2 = ScaleData(rna2, features = all.genes)
rna2[["replicate"]] <- "replicate 2"

# Perform linear dimensional reduction
rna2 = RunPCA(rna2, features = VariableFeatures(object = rna2))
rna2 = FindNeighbors(rna2, dims = 1:10)
rna2 = FindClusters(rna2, resolution = 0.5)
# Run non-linear dimensional reduction
rna2 = RunUMAP(
  rna2,
  dims = 1:10,
  n.neighbors = 40,
  n.components = 2
)
umap2 = DimPlot(rna2, reduction = "umap")
umap2

# replication 2
ggsave(
  glue("{result_folder}UMAP_nn40_ncomp2_scRNASeq_rep2_Marek.png"),
  plot = umap2,
  width = 10,
  height = 10,
  dpi = 300,
)

## integrate the two replicates by CCA-based Seurat integration protocol
rna = list(rna1, rna2)

rm(rna1)
rm(rna2)

# anchor finding (searching best "buddies")
anchors = FindIntegrationAnchors(object.list = rna,
                                 anchor.features = 2000,
                                 dims = 1:30)
integrate = IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(integrate) = "integrated"
integrate = ScaleData(integrate, do.center = T, do.scale = F)
integrate = RunPCA(
  integrate,
  npcs = 40,
  ndims.print = 1:5,
  nfeatures.print = 5
)
Idents(integrate) = "replicate"

# PCA
pca_int = DimPlot(
  integrate,
  dims = c(1, 2),
  reduction = "pca",
  split.by = "replicate"
)

ggsave(
  glue("{result_folder}PCA_integrated_scRNA-Seq_Marek.png"),
  plot = pca_int,
  width = 10,
  height = 10,
  dpi = 300,
)

# UMAP
integrate <-
  RunUMAP(
    integrate,
    dims = 1:30,
    reduction = "pca",
    n.neighbors = 15,
    min.dist = 0.5,
    spread = 1,
    metric = "euclidean",
    seed.use = 1
  )
umap_int1 <-
  DimPlot(integrate, reduction = "umap", group.by = "seurat_clusters")
umap_int1
umap_int2 <-
  DimPlot(integrate, reduction = "umap", group.by = "replicate")
umap_int2
umap_grid = plot_grid(umap_int1,
                      umap_int2)

ggsave(
  glue("{result_folder}UMAPs_integrated_scRNA-Seq_Marek.png"),
  plot = umap_grid,
  width = 10,
  height = 10,
  dpi = 300,
)

# exporting integrated Seurat object
saveRDS(integrate, glue("{result_folder}integration/mneural_scRNA-Seq_Marek.Rds"))

# fetching gene level data
# example:
# gene = FetchData(integrate, vars = "Plp1")



