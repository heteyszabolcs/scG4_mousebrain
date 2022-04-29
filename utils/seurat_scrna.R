# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("GenomicFeatures")
})

# scRNA-Seq folder (Marek's mouse brain scRNA-Seq data)
cellranger_folder = "../data/GSE163484/"
# path to result folder
result_folder = "../results/Seurat/"

# read and process data
counts1 = Read10X_h5(filename = glue("{cellranger_folder}scRNASeq_GSM4979874_filtered_feature_bc_matrix_rep1.h5"))
rna1 = CreateSeuratObject(counts = counts1, project = "mouse_brain", min.cells = 3, min.features = 200)
rna1
                          
counts2 = Read10X_h5(filename = glue("{cellranger_folder}scRNASeq_GSM4979875_filtered_feature_bc_matrix_rep2.h5"))
rna2 = CreateSeuratObject(counts = counts2, project = "mouse_brain", min.cells = 3, min.features = 200)
rna2

rna = merge(rna1, y = rna2, add.cell.ids = c("rep1", "rep2"), project = "mouse_brain")

rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(rna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(rna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
varfeatures = plot1 + plot2

ggsave(
  glue("{result_folder}VariableFeatures_scRNASeq_Marek.png"),
  plot = varfeatures,
  width = 14,
  height = 8,
  dpi = 300,
)

all.genes = rownames(rna)
rna = ScaleData(rna, features = all.genes)

# Perform linear dimensional reduction
rna = RunPCA(rna, features = VariableFeatures(object = rna))
rna = FindNeighbors(rna, dims = 1:10)
rna = FindClusters(rna, resolution = 0.5)
# Run non-linear dimensional reduction
rna = RunUMAP(rna, dims = 1:10, n.neighbors = 40, n.components = 2)
umap = DimPlot(rna, reduction = "umap")

ggsave(
  glue("{result_folder}UMAP_nn40_ncomp2_scRNASeq_Marek.png"),
  plot = umap,
  width = 10,
  height = 10,
  dpi = 300,
)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
rna.markers = FindAllMarkers(rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.0)
rna.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# exporting
saveRDS(rna, "{result_folder}mneural_scRNA-Seq_Marek.Rds")
write_tsv(rna.markers, "{result_folder}mneural_scRNA-Seq_Marek-markers.tsv")


