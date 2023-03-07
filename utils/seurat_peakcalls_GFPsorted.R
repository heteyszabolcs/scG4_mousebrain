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
  library("Matrix")
  library("ggpubr")
  library("DoubletFinder")
})

set.seed(5)

# Cellranger-ATAC folder
cellranger_folder = "../data/CellRanger/GFP_sorted/"
# path to result folder
result_folder = "../results/Seurat/callpeaks_GFPsorted/"

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
g4 = FindNeighbors(object = g4,
                   reduction = 'lsi',
                   dims = 2:30)
g4 = FindClusters(object = g4,
                  verbose = FALSE,
                  algorithm = 3)

# QC violin plots
nF_violin = VlnPlot(g4, group.by = "seurat_clusters", features = "nFeature_peaks", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
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
TSS_violin = VlnPlot(g4, group.by = "seurat_clusters", features = "TSS_fragments", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
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
mito_violin = VlnPlot(g4, group.by = "seurat_clusters", features = "mitochondrial", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
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

#find doublets using DoubletFinder
g4_doublet_test = Seurat::Read10X_h5(filename = glue("{cellranger_folder}filtered_peak_bc_matrix.h5"), use.names = T)
g4_doublet_test = CreateSeuratObject(g4_doublet_test, project = "unsorted")
g4_doublet_test = NormalizeData(g4_doublet_test)
g4_doublet_test = FindVariableFeatures(g4_doublet_test, verbose = F)
g4_doublet_test = ScaleData(g4_doublet_test, vars.to.regress = c("nFeature_RNA", "percent_mito"),
                            verbose = F)
g4_doublet_test = RunPCA(g4_doublet_test, verbose = F, npcs = 20)
g4_doublet_test = RunUMAP(g4_doublet_test, dims = 1:10, verbose = F)

nExp <- round(ncol(g4_doublet_test) * 0.04)  # expect 4% doublets
g4_doublet_test <- doubletFinder_v3(g4_doublet_test, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
singlets = rownames(g4_doublet_test@meta.data[which(g4_doublet_test@meta.data[,5] == "Singlet"),])

# keep only singlets!
meta = g4@meta.data
g4@meta.data = meta %>% mutate(doublet_test = ifelse(rownames(meta) %in% singlets, "Singlet", NA_character_))
g4 = subset(x = g4, subset = doublet_test == "Singlet")
# removing cluster 2 - cluster 2 can be omitted
g4 = subset(x = g4, idents = 4, invert = TRUE)

dim = DimPlot(object = g4, label = TRUE, pt.size = 6) + 
  NoLegend() +
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("GFP-sorted scCutnTag")

ggsave(
  glue("{result_folder}Seurat_GFPsorted_UMAP.png"),
  plot = dim,
  width = 10,
  height = 10,
  dpi = 300,
)

# QC violin plots
nF_violin = VlnPlot(g4, group.by = "seurat_clusters", features = "nFeature_peaks", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("nFeature (peaks)") +
  xlab("cluster") + 
  ylim(0, 2000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
TSS_violin = VlnPlot(g4, group.by = "seurat_clusters", features = "TSS_fragments", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("TSS fragments") +
  xlab("cluster") + 
  ylim(0, 2000) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
mito_violin = VlnPlot(g4, group.by = "seurat_clusters", features = "mitochondrial", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("mitochondrial fragments") +
  xlab("cluster") + 
  ylim(0, 5000) +
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

# Create a gene activity matrix
# Signac function: https://stuartlab.org/signac/reference/geneactivity
# Compute counts per cell in gene body and promoter region
gene.activities = GeneActivity(g4)

g4[['GA']] = CreateAssayObject(counts = gene.activities)
g4 = NormalizeData(
  object = g4,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(g4$nCount_peaks)
)

# find marker genes
g4_markers = FindAllMarkers(g4, only.pos = TRUE, assay = "GA")
write_tsv(g4_markers, glue("{result_folder}FindAllMarkers_GA_output.tsv"))

g4_top.markers = g4_markers %>% group_by(cluster) %>% top_n(wt = avg_log2FC,n = 5)
DoHeatmap(object = g4, features = g4_top.markers$gene, slot = 'data', raster = TRUE, assay = "GA") + 
  scale_fill_gradient2(low="white", mid="white", high="red") +
  labs(title = " ", fill = "gene act. score") +
  theme(
    plot.title = element_text(size = 8),
    axis.text.y = element_text(size = 13, color = "black")
  )

ggsave(
  glue("{result_folder}GA_heatmap.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 10,
  height = 3
)

# extract cell barcodes per cluster
barcodes = g4@meta.data %>% 
  rownames_to_column("barcodes") %>%
  dplyr::select(barcodes, seurat_clusters)
write_tsv(barcodes, glue("{result_folder}barcodes_per_cluster.tsv"))

for(cluster in unique(barcodes$seurat_clusters)) {
  subset = barcodes %>% dplyr::filter(seurat_clusters == cluster)
  write_tsv(subset, glue("{result_folder}barcodes_cluster_{as.character(cluster)}.tsv"), col_names = FALSE)
}

# export Rds
saveRDS(g4, glue("{result_folder}GFPsorted.Rds"))

# peak calling per Seurat clusters
peaks = CallPeaks(
  object = g4,
  group.by = "seurat_clusters",
  cleanup = FALSE,
  outdir = result_folder,
  effective.genome.size = 2652783500
)

write.table(
  as.data.frame(peaks),
  file = glue("{result_folder}peaks_per_clusters.bed"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)



