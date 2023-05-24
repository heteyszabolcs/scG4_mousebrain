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
  library("data.table")
})

# scRNA-Seq folder (Marek's postnatal mouse brain scRNA-Seq data)
cellranger_folder = "../data/GSE163484/"
marques_folder = "../data/GSE75330/"
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
  n.neighbors = 20,
  min.dist = 0.010,
  n.components = 2
)
umap1 = DimPlot(rna1, reduction = "umap", pt.size = 3, label = TRUE) + NoLegend()
umap1

FeaturePlot(rna1, features = "Pdgfra", reduction = "umap", pt.size = 3, cols = c("#f0f0f0", "#de2d26")) 

ggsave(
  glue("{result_folder}UMAP_nn10_ncomp2_scRNASeq_rep1_Marek.png"),
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
umap2 = DimPlot(rna2, reduction = "umap", pt.size = 3, label = TRUE) + NoLegend()
umap2

ggsave(
  glue("{result_folder}UMAP_nn40_ncomp2_scRNASeq_rep2_Marek.png"),
  plot = umap2,
  width = 10,
  height = 10,
  dpi = 300,
)

# create normalized expression data
rna2 = NormalizeData(rna2, normalization.method = "LogNormalize")
lognorm = rna2[["RNA"]]@data
#writeMM(lognorm, glue("{result_folder}scRNASeq_GSM4979875_rep1_norm.mtx"))

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

# scRNA of Marques et al. - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5221728/
marques = fread(glue("{marques_folder}GSE75330_Marques_et_al_mol_counts2.tab"))
rownames(marques) = marques$cellid
marques = marques %>% dplyr::select(-cellid)
marques = CreateSeuratObject(counts = marques)

marques[["percent.mt"]] <- PercentageFeatureSet(marques, pattern = "^MT-")
marques = FindVariableFeatures(marques, selection.method = "vst", nfeatures = 2000)

# Perform linear dimensional reduction
all.genes = rownames(marques)
marques = ScaleData(marques, features = all.genes)
marques = RunPCA(marques, features = VariableFeatures(object = marques))
marques = FindNeighbors(marques, dims = 1:10)
marques = FindClusters(marques, resolution = 0.5)
# Run non-linear dimensional reduction
marques = RunUMAP(
  marques,
  dims = 1:10,
  n.neighbors = 40,
  n.components = 2
)


#marques_markers = FindAllMarkers(marques, logfc.threshold = 0.20)

## marker genes are from Marques et al paper
# OPC markers (progenitor oligodendrocytes)
opc = c("Pdgfra", "Ptprz1", "Cspg4", "Vcan") 
opc_featureplot = FeaturePlot(marques, features = opc, reduction = "umap", pt.size = 1, cols = c("#f0f0f0", "red")) 

ggsave(
  glue("{result_folder}Feature_OPC_scRNASeq_GSE75330.png"),
  plot = opc_featureplot,
  width = 7,
  height = 7,
  dpi = 500,
)

# COP markers (commited oligodendrocytes)
cop = c("Vcan", "Sox6", "Gpr17", "Neu4")
cop_featureplot = FeaturePlot(marques, features = cop, reduction = "umap", pt.size = 1, cols = c("#f0f0f0", "red")) 

ggsave(
  glue("{result_folder}Feature_COP_scRNASeq_GSE75330.png"),
  plot = cop_featureplot,
  width = 7,
  height = 7,
  dpi = 500,
)

# MOL markers (mature oligodendrocytes)
mol = c("Fosb", "Btg2", "Apod", "Klk6")
mol_featureplot = FeaturePlot(marques, features = mol, reduction = "umap", pt.size = 1, cols = c("#f0f0f0", "red"))

ggsave(
  glue("{result_folder}Feature_MOL_scRNASeq_GSE75330.png"),
  plot = mol_featureplot,
  width = 7,
  height = 7,
  dpi = 500,
)

# VLMC markers (mature oligodendrocytes)
vlmc = c("Tbx18", "Vtn", "Lum", "Col1a2")
vlmc_featureplot = FeaturePlot(marques, features = vlmc, reduction = "umap", pt.size = 1, cols = c("#f0f0f0", "red")) 

# MFOL markers (myelin-forming oligodendrocytes)
mfol = c("Ctps", "Opalin", "Serinc5")
mfol_featureplot = FeaturePlot(marques, features = mfol, reduction = "umap", pt.size = 1, cols = c("#f0f0f0", "red")) 

ggsave(
  glue("{result_folder}Feature_MFOL_scRNASeq_GSE75330.png"),
  plot = mfol_featureplot,
  width = 7,
  height = 7,
  dpi = 500,
)

# NFOL markers (newly formed oligodendrocytes)
nfol = c("Gpr17", "Nkx2-2", "Bmp4", "Neu4")
nfol_featureplot = FeaturePlot(marques, features = nfol, reduction = "umap", pt.size = 1, cols = c("#f0f0f0", "red")) 

ggsave(
  glue("{result_folder}Feature_VLMC_scRNASeq_GSE75330.png"),
  plot = vlmc_featureplot,
  width = 7,
  height = 7,
  dpi = 500,
)

new.cluster.ids = c("MOL", "MFOL/MOL", "MOL", "MOL", "MOL", "MFOL", "COP", "OPC", "MOL", "MOL", "VLMC")
names(new.cluster.ids) = levels(marques)
marques = RenameIdents(marques, new.cluster.ids)

umap3 = DimPlot(marques, reduction = "umap", pt.size = 3, label = TRUE) + NoLegend()
umap3

ggsave(
  glue("{result_folder}scRNASeq_GSE75330_UMAP.png"),
  plot = umap3,
  width = 10,
  height = 10,
  dpi = 300,
)

marques[["cell_type"]] = marques@active.ident

# export
saveRDS(marques, glue("{result_folder}scRNASeq_GSE75330.rds"))

# fetching gene level data
# example:
# gene = FetchData(integrate, vars = "Plp1")

# log norm data of Marques et al. scRNA
marques_lognorm = as.data.frame(marques[["RNA"]]@data)
write.table(marques_lognorm, glue("{result_folder}scRNASeq_GSE75330_lognorm.tsv"), row.names = TRUE, quote = FALSE, col.names = TRUE)




