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
ptm_objects = "../data/GSE157637/"

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

saveRDS(rna1, glue("{result_folder}scRNASeq_GSM4979874_rep1.rds"))

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

# create normalized expression data
rna2 = NormalizeData(rna2, normalization.method = "LogNormalize")
lognorm = rna2[["RNA"]]@data
writeMM(lognorm, glue("{result_folder}scRNASeq_GSM4979875_rep1_norm.mtx"))

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

## integrate the replicate1 with H3K4me3 scC&T data by CCA-based Seurat integration protocol
# load Bartosovic H3K4me3 scC&T data and process it:
k4me3 = readRDS(glue("{ptm_objects}H3K4me3_seurat_object.Rds"))
k4me3 = AddMetaData(object = k4me3,
                      metadata = "scCut&Tag H3K4me3",
                      col.name = "group")
DefaultAssay(k4me3) = "GA"
k4me3[["bins_5000"]] = NULL
k4me3[["peaks"]] = NULL
k4me3[["PA"]] = NULL


k4me3 = FindVariableFeatures(k4me3, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(k4me3)
k4me3 = ScaleData(k4me3, features = all.genes)
k4me3[["replicate"]] <- "H3K4me3"

# Perform linear dimensional reduction
k4me3 = RunPCA(k4me3, features = VariableFeatures(object = k4me3))
k4me3 = FindNeighbors(k4me3, dims = 1:10)
k4me3 = FindClusters(k4me3, resolution = 0.5)
# Run non-linear dimensional reduction
k4me3 = RunUMAP(
  k4me3,
  dims = 1:10,
  n.neighbors = 40,
  n.components = 2
)
k4me3_umap = DimPlot(k4me3, reduction = "umap")
k4me3_umap

ggsave(
  glue("{result_folder}UMAP_H3K4me3_scCT_Marek.png"),
  plot = k4me3_umap,
  width = 10,
  height = 10,
  dpi = 300,
)

# load scRNA rep1 data and integrate:
rna1 = readRDS(glue("{result_folder}scRNASeq_GSM4979874_rep1.rds"))

# integrate the H3K4me3 data with scRNA-Seq rep1 by CCA-based Seurat integration protocol
int = list(rna1, k4me3)
rm(rna1)
rm(k4me3)

# anchor finding (searching best "buddies")
rna1_k4me3_anchors = FindIntegrationAnchors(object.list = int,
                                 anchor.features = 2000,
                                 dims = 1:30)
rm(int)

rna1_k4me3_integrate = IntegrateData(anchorset = rna1_k4me3_anchors, dims = 1:30)
DefaultAssay(rna1_k4me3_integrate) = "integrated"
rm(rna1_k4me3_anchors)

rna1_k4me3_integrate = ScaleData(rna1_k4me3_integrate, do.center = T, do.scale = F)
rna1_k4me3_integrate = RunPCA(
  rna1_k4me3_integrate,
  npcs = 40,
  ndims.print = 1:5,
  nfeatures.print = 5
)
Idents(rna1_k4me3_integrate) = "replicate"

# PCA
rna1_k4me3_pca_int = DimPlot(
  rna1_k4me3_integrate,
  dims = c(1, 2),
  reduction = "pca",
  split.by = "replicate"
)
rna1_k4me3_pca_int

ggsave(
  glue("{result_folder}PCA_integrated_H3K4me3-scRNASeq_Marek.png"),
  plot = rna1_k4me3_pca_int,
  width = 10,
  height = 10,
  dpi = 300,
)

# UMAP
rna1_k4me3_integrate <-
  RunUMAP(
    rna1_k4me3_integrate,
    dims = 1:30,
    reduction = "pca",
    n.neighbors = 15,
    min.dist = 0.5,
    spread = 1,
    metric = "euclidean",
    seed.use = 1
  )
rna1_k4me3_umap_int1 <-
  DimPlot(rna1_k4me3_integrate, reduction = "umap", group.by = "seurat_clusters")
rna1_k4me3_umap_int1
rna1_k4me3_umap_int2 <-
  DimPlot(rna1_k4me3_integrate, reduction = "umap", group.by = "replicate")
rna1_k4me3_umap_int2
rna1_k4me3_umap_int3 <-
DimPlot(rna1_k4me3_integrate, reduction = "umap", group.by = "cell_type")
rna1_k4me3_umap_int3
rna1_k4me3_umap_grid = plot_grid(rna1_k4me3_umap_int1,
                                 rna1_k4me3_umap_int2,
                                 rna1_k4me3_umap_int3)

ggsave(
  glue("{result_folder}UMAP_integrated_H3K4me3-scRNASeq_Marek.png"),
  plot = rna1_k4me3_umap_grid,
  width = 10,
  height = 10,
  dpi = 300,
)

# exporting integrated Seurat object
saveRDS(rna1_k4me3_integrate, glue("{result_folder}integration/H3K4me3_scRNA-Seq_Marek.Rds"))


