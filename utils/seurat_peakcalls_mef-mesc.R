# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("EnsDb.Mmusculus.v79")
  library("GenomicFeatures")
  library("GenomeInfoDb")
  library("cowplot")
  #library("DoubletFinder")
})

set.seed(5)

# Cellranger-ATAC folder
cellranger_folder = "../data/CellRanger/mES-mEF/"
# path to result folder
result_folder = "../results/Seurat/callpeaks_mESC-MEF/"

# read and process data
counts = Read10X_h5(filename = glue("{cellranger_folder}filtered_peak_bc_matrix.h5"))
metadata = read.csv(
  file = glue("{cellranger_folder}singlecell.csv"),
  header = TRUE,
  row.names = 1
)

chrom_assay = CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = glue("{cellranger_folder}fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)

### clustering with resolution 0.1
g4_res0.1 = CreateSeuratObject(counts = chrom_assay,
                        assay = "peaks",
                        meta.data = metadata)

# check peaks in GRanges format
granges(g4_res0.1)

# give EnsDB mm10 annotation
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) = 'UCSC'
Annotation(g4_res0.1) = annotations

# normalization
g4_res0.1 = RunTFIDF(g4_res0.1)
g4_res0.1 = FindTopFeatures(g4_res0.1, min.cutoff = 'q0')
g4_res0.1 = RunSVD(g4_res0.1)
DepthCor(g4_res0.1) # 1st component highly correlate with sequencing depth

# Non-linear dimension reduction and clustering
g4_res0.1 = RunUMAP(object = g4_res0.1,
             reduction = 'lsi',
             dims = 2:30)
g4_res0.1 = FindNeighbors(object = g4_res0.1,
                   reduction = 'lsi',
                   dims = 2:30)
g4_res0.1 = FindClusters(object = g4_res0.1,
                  verbose = FALSE,
                  resolution = 0.1,
                  algorithm = 3)

# QC violin plots
nF_violin = VlnPlot(g4_res0.1, group.by = "seurat_clusters", features = "nFeature_peaks", pt.size = 0.1) +
  scale_fill_manual(values = c("#9ecae1", "#fc9272")) +
  ggtitle("nFeature (peaks)") +
  xlab("cluster") + 
  ylim(0, 15000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
TSS_violin = VlnPlot(g4_res0.1, group.by = "seurat_clusters", features = "TSS_fragments", pt.size = 0.1) +
  scale_fill_manual(values = c("#9ecae1", "#fc9272")) +
  ggtitle("TSS fragments") +
  xlab("cluster") + 
  ylim(0, 15000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
mito_violin = VlnPlot(g4_res0.1, group.by = "seurat_clusters", features = "mitochondrial", pt.size = 0.1) +
  scale_fill_manual(values = c("#9ecae1", "#fc9272")) +
  ggtitle("mitochondrial fragments") +
  xlab("cluster") + 
  ylim(0, 15000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
qc_violins = ggarrange(nF_violin, TSS_violin, mito_violin)

ggsave(
  glue("{result_folder}QC_violins.png"),
  plot = qc_violins,
  width = 12,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}QC_violins.pdf"),
  plot = qc_violins,
  device = "pdf",
  width = 12,
  height = 10,
  dpi = 300,
)

cols = c(
  "0" = "#9ecae1",
  "1" = "#fc9272"
)

dim_res0.1 = DimPlot(object = g4_res0.1, label = TRUE, pt.size = 1.5, label.size = 7, repel = TRUE) + 
  #NoLegend() +
  #scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  scale_colour_manual(values = cols, breaks = c("0", "1"), labels = c("MEF", "mESC")) +
  ggtitle("mESC-MEF G4 scCut&Tag, res = 0.1") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim_res0.1

### clustering without resolution contraints
g4 = CreateSeuratObject(counts = chrom_assay,
                        assay = "peaks",
                        meta.data = metadata)

# check peaks in GRanges format
granges(g4)

# give EnsDB mm10 annotation
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) = 'UCSC'
Annotation(g4) = annotations

# normalization
g4 = RunTFIDF(g4)
g4 = FindTopFeatures(g4, min.cutoff = 'q0')
g4 = RunSVD(g4)
DepthCor(g4) # 1st component highly correlate with sequencing depth

# Non-linear dimension reduction and clustering
g4 = RunUMAP(object = g4,
             reduction = 'lsi',
             dims = 2:30)

# find doublets using DoubletFinder 
# g4_doublet_test = Seurat::Read10X_h5(filename = glue("{cellranger_folder}filtered_peak_bc_matrix.h5"), use.names = T)
# g4_doublet_test = CreateSeuratObject(g4_doublet_test, project = "MEF-mESC")
# g4_doublet_test = NormalizeData(g4_doublet_test)
# g4_doublet_test = FindVariableFeatures(g4_doublet_test, verbose = F)
# g4_doublet_test = ScaleData(g4_doublet_test, vars.to.regress = c("nFeature_RNA", "percent_mito"),
#                       verbose = F)
# g4_doublet_test = RunPCA(g4_doublet_test, verbose = F, npcs = 20)
# g4_doublet_test = RunUMAP(g4_doublet_test, dims = 1:10, verbose = F)
# 
# nExp <- round(ncol(g4_doublet_test) * 0.04)  # expect 4% doublets
# g4_doublet_test <- doubletFinder_v3(g4_doublet_test, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
# DF.name = colnames(g4_doublet_test@meta.data)[grepl("DF.classification", colnames(g4_doublet_test@meta.data))]
# umap = DimPlot(g4_doublet_test, group.by = "orig.ident", pt.size = 0.3) +
#   xlim(-10, 10) +
#   ylim(-10, 10) +
#   ggtitle("MEF-mESC G4 scCnT")
# 
# # indicate doublets
# doublets = DimPlot(g4_doublet_test, group.by = DF.name, pt.size = 0.3) +
#   xlim(-10, 10) +
#   ylim(-10, 10) +
#   ggtitle("Doublet finder")
# 
# # compare UMAPs
# plot_grid(umap, doublets, nrow = 1, ncol = 2)
# 
# ggsave(
#   glue("{result_folder}Seurat_mESC-MEF_doublets.png"),
#   plot = last_plot(),
#   width = 10,
#   height = 5,
#   dpi = 300,
# )

## carry on with dimension plots
g4 = FindNeighbors(object = g4,
                   reduction = 'lsi',
                   dims = 2:30)
g4 = FindClusters(object = g4,
                  verbose = FALSE,
                  algorithm = 3)


dim = DimPlot(object = g4, label = TRUE, pt.size = 2, label.size = 7, repel = TRUE) + 
  NoLegend() +
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("mESC-MEF G4 scCut&Tag") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim

dim_blank = DimPlot(object = g4, label = TRUE, pt.size = 2, label.size = 7, repel = TRUE) + 
  NoLegend() +
  NoAxes() +
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  #ggtitle("mESC-MEF G4 scCut&Tag") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 0, color = "black"),
    axis.text.y = element_text(size = 0, color = "black")
  )
dim_blank

# tuning
# for(i in seq(5, 50, by = 5)) {
#   tune = RunUMAP(object = g4,
#                reduction = 'lsi',
#                n.neighbors = i,
#                dims = 2:30)
#   print(DimPlot(object = tune, label = TRUE, pt.size = 6) + 
#     NoLegend() +
#     scale_color_brewer(palette = "Set3") +
#     xlim(-10, 10) + 
#     ylim(-10, 10) + 
#     ggtitle(as.character(i)))
#   
# }
  
ggsave(
  glue("{result_folder}Seurat_mESC-MEF_UMAP.png"),
  plot = dim,
  width = 10,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_mESC-MEF_UMAP.pdf"),
  plot = dim,
  width = 10,
  height = 10,
  device = "pdf"
)

ggsave(
  glue("{result_folder}Seurat_mESC-MEF_UMAP_blankedplot.pdf"),
  plot = dim_blank,
  width = 10,
  height = 10,
  device = "pdf"
)

ggsave(
  glue("{result_folder}Seurat_mESC-MEF_UMAP_blankedplot.png"),
  plot = dim_blank,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  glue("{result_folder}Seurat_mESC-MEF_UMAP_res0.1.png"),
  plot = dim_res0.1,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  glue("{result_folder}Seurat_mESC-MEF_UMAP_res0.1.pdf"),
  plot = dim_res0.1,
  width = 10,
  height = 10,
  device = "pdf"
)

mesc_mef_umaps = plot_grid(dim, dim_res0.1, ncol = 2, nrow = 1)

ggsave(
  glue("{result_folder}Seurat_mESC-MEF_UMAPs.png"),
  plot = mesc_mef_umaps,
  width = 10,
  height = 5,
  dpi = 300,
)

### Create a gene activity matrix
gene.activities = GeneActivity(g4_res0.1)

g4_res0.1[['GA']] = CreateAssayObject(counts = gene.activities)
g4_res0.1 = NormalizeData(
  object = g4_res0.1,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(g4$nCount_peaks)
)

# extract cell barcodes per cluster
barcodes = g4_res0.1@meta.data %>% 
  rownames_to_column("barcodes") %>%
  dplyr::select(barcodes, seurat_clusters)
write_tsv(barcodes, glue("{result_folder}barcodes_per_cluster_res0.1.tsv"))

for(cluster in unique(barcodes$seurat_clusters)) {
  subset = barcodes %>% dplyr::filter(seurat_clusters == cluster)
  write_tsv(subset, glue("{result_folder}barcodes_cluster_{as.character(cluster)}_res0.1.tsv"), col_names = FALSE)
}

# export Rds
saveRDS(g4_res0.1, glue("{result_folder}mESC-MEF.Rds"))

# peak calling per Seurat clusters
peaks = CallPeaks(
  object = g4_res0.1,
  group.by = "seurat_clusters",
  cleanup = FALSE,
  outdir = result_folder,
  effective.genome.size = 2652783500
)

write.table(
  as.data.frame(peaks),
  file = glue("{result_folder}peaks_per_clusters_res0.1.bed"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)



