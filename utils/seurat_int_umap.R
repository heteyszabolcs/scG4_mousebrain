# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("ggplot2")
  library("EnsDb.Mmusculus.v79")
  library("ensembldb")
  library("GenomicRanges")
  library("dplyr")
  library("glue")
  library("tidyverse")
  library("data.table")
  library("gridExtra")
  library("cowplot")
  library("ggrastr")
  library("RColorBrewer")
  library("scclusteval")
})

set.seed(5)

# path to result folder
result_folder = "../results/Seurat/callpeaks_unsorted/"

# Seurat object (brain unsorted)
unsorted = readRDS(file = "../results/Seurat/callpeaks_unsorted/unsorted.Rds")

# Marques et al. oligo scRNA data
rna = read.table(
  "../data/GSE75330/GSE75330_Marques_et_al_mol_counts2.tab",
  stringsAsFactors = FALSE,
  header = FALSE
)
annot = readRDS(file = "../data/GSE75330/Marques2016annotation.rds")

colnames(rna) = gsub("-", "_", sub('-', '', gsub("C1-", "", rna[1,])))
rna = rna[-1,]
genes = rna[, 1]
rna = rna[,-1]
rna = as.matrix(rna)
rna = apply(rna, 2, as.numeric)
rownames(rna) = genes

rna = CreateSeuratObject(counts = rna,
                         meta.data = annot,
                         assay = 'RNA')

# Seurat workflow
all.genes = rownames(rna)
rna = NormalizeData(rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna = ScaleData(rna, features = all.genes)
rna = RunPCA(rna, features = VariableFeatures(object = rna))
rna = RunUMAP(rna, dims = 1:20)

new_ids = as.character(rna@meta.data$cell_class)
new_ids[new_ids == 'NFOL2'] = 'NFOL'
new_ids[new_ids == 'NFOL1'] = 'NFOL'
new_ids[new_ids == 'MFOL2'] = 'MFOL'
new_ids[new_ids == 'MFOL1'] = 'MFOL'
new_ids[new_ids == 'MOL1'] = 'MOL'
new_ids[new_ids == 'MOL2'] = 'MOL'
new_ids[new_ids == 'MOL3'] = 'MOL'
new_ids[new_ids == 'MOL4'] = 'MOL'
new_ids[new_ids == 'MOL5'] = 'MOL'
new_ids[new_ids == 'MOL6'] = 'MOL'
rna@meta.data$merged_cell_class = new_ids

p1 = DimPlot(unsorted,
             pt.size = 2,
             label.size = 7,
             repel = TRUE,
             raster = TRUE) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
p2 = DimPlot(
  rna,
  pt.size = 2,
  label.size = 7,
  group.by = 'merged_cell_class',
  repel = TRUE,
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

# find all markers in scRNA-Seq data
rna@active.ident = rna$cell_class
markers = FindAllMarkers(rna)
markers.pos = markers[markers$p_val < 0.05 &
                        markers$avg_logFC > 0.5,]

write_tsv(markers,
          "../results/Seurat/scRNA-Seq_Marques_et_al-FindAllMarkers_output.tsv")

DefaultAssay(unsorted) = "GA"
DefaultAssay(rna) = "RNA"

common.genes = intersect(rownames(rna), rownames(unsorted))

# anchor identification between G4 scCnT and scRNA-Seq data sets
transfer.anchors = FindTransferAnchors(
  reference = rna,
  query = unsorted,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

genes.use = VariableFeatures(rna)
refdata = GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

# imputation - data transfer
imputation = TransferData(
  anchorset = transfer.anchors,
  refdata = refdata,
  weight.reduction = unsorted[["lsi"]],
  dims = 1:50
)

unsorted[['RNA']] = imputation
saveRDS(unsorted, glue("{result_folder}unsorted_int_Marques.Rds"))
unsorted = readRDS(glue("{result_folder}unsorted_int_Marques.Rds"))

rna@meta.data$data_type = "scRNA-Seq (Marques et al.)"
unsorted@meta.data$data_type = "unsorted G4 scCut&Tag"

coembed = merge(x = rna, y = unsorted)

coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed = RunUMAP(coembed, dims = 1:30)

new_ids = as.character(coembed@meta.data$cell_class)
new_ids[new_ids == 'NFOL2'] = 'NFOL'
new_ids[new_ids == 'NFOL1'] = 'NFOL'
new_ids[new_ids == 'MFOL2'] = 'MFOL'
new_ids[new_ids == 'MFOL1'] = 'MFOL'
new_ids[new_ids == 'MOL1'] = 'MOL'
new_ids[new_ids == 'MOL2'] = 'MOL'
new_ids[new_ids == 'MOL3'] = 'MOL'
new_ids[new_ids == 'MOL4'] = 'MOL'
new_ids[new_ids == 'MOL5'] = 'MOL'
new_ids[new_ids == 'MOL6'] = 'MOL'
coembed@meta.data$cell_class = new_ids

# Ridge p
# RidgePlot(object = dataB, features.plot = rownames(cluster1.markers)[1:6])
cols = brewer.pal(n = 6, name = "Greens")
RidgePlot(object = coembed, features = "Th" , assay = "GA", group.by = "cell_class", cols = cols)

## UMAP dimension plots
DimPlot(rna)
coembed_cells = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'cell_class',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed_cells2 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'cell_class',
  repel = TRUE,
  label = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend() + NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_celltypes.pdf"),
  plot = coembed_cells2,
  width = 10,
  height = 10,
  device = "pdf"
)

coembed_gfp = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'GFP',
  repel = TRUE,
  na.value = "grey0",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_color_manual(values = c("#bdbdbd", "#addd8e")) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed_clusters = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed_experiments = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'data_type',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_colour_manual(values = c("#fc9272", "#9ecae1")) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

coembed@meta.data$seurat_clusters[is.na(coembed@meta.data$seurat_clusters)] = "scRNA-Seq"

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#de2d26",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd",
  "6" = "#bdbdbd"
)

coembed_cluster0 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "0",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols,
                      breaks = "0",
                      labels = "cluster 0") +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl0.pdf"),
  plot = coembed_cluster0,
  width = 10,
  height = 10,
  device = "pdf"
)

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#de2d26",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd",
  "6" = "#bdbdbd"
)

coembed_cluster1 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "1",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols,
                      breaks = "1",
                      labels = "cluster 1") +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl1.pdf"),
  plot = coembed_cluster1,
  width = 10,
  height = 10,
  device = "pdf"
)

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#de2d26",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd",
  "6" = "#bdbdbd"
)

coembed_cluster2 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "2",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols,
                      breaks = "2",
                      labels = "cluster 2") +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl2.pdf"),
  plot = coembed_cluster2,
  width = 10,
  height = 10,
  device = "pdf"
)

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#de2d26",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd",
  "6" = "#bdbdbd"
)

coembed_cluster3 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "3",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols,
                      breaks = "3",
                      labels = "cluster 3") +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl3.pdf"),
  plot = coembed_cluster3,
  width = 10,
  height = 10,
  device = "pdf"
)

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd",
  "4" = "#de2d26",
  "5" = "#bdbdbd",
  "6" = "#bdbdbd"
)

coembed_cluster4 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "4",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols,
                      breaks = "4",
                      labels = "cluster 4") +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl4.pdf"),
  plot = coembed_cluster4,
  width = 10,
  height = 10,
  device = "pdf"
)

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd",
  "5" = "#de2d26",
  "6" = "#bdbdbd"
)

coembed_cluster5 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "5",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols,
                      breaks = "5",
                      labels = "cluster 5") +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl5.pdf"),
  plot = coembed_cluster5,
  width = 10,
  height = 10,
  device = "pdf"
)

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd",
  "6" = "#de2d26"
)

coembed_cluster6 = DimPlot(
  coembed,
  pt.size = 1,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "6",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols,
                      breaks = "6",
                      labels = "cluster 6") +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl6.pdf"),
  plot = coembed_cluster6,
  width = 10,
  height = 10,
  device = "pdf"
)

clusters = plot_grid(
  coembed_cluster0,
  coembed_cluster1,
  coembed_cluster2,
  coembed_cluster3,
  coembed_cluster4,
  coembed_cluster5,
  coembed_cluster6,
  ncol = 2,
  nrow = 4
)

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl_all.pdf"),
  plot = clusters,
  width = 10,
  height = 10,
  device = "pdf"
)

ggsave(
  glue("{result_folder}Seurat_Marques_integration_cl_all.png"),
  plot = clusters,
  width = 10,
  height = 10,
  dpi = 500
)

coembed_ps = coembed_cells + coembed_gfp + coembed_clusters + coembed_experiments
coembed_ps

ggsave(
  glue("{result_folder}Seurat_Marques_coembeds.pdf"),
  plot = coembed_ps,
  width = 16,
  height = 6,
  device = "pdf"
)

coembed_ps = coembed_cells + coembed_clusters
coembed_ps

ggsave(
  glue("{result_folder}Seurat_Marques_coembeds.pdf"),
  plot = coembed_ps,
  width = 12,
  height = 6,
  device = "pdf"
)

umaps = plot_grid(p1,
                  p2,
                  coembed_clusters,
                  coembed_cells,
                  ncol = 2,
                  nrow = 2)

ggsave(
  glue("{result_folder}Seurat_unsorted_Marques_UMAPs.png"),
  plot = umaps,
  width = 12,
  height = 12,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_unsorted_Marques_UMAPs.pdf"),
  plot = umaps,
  width = 12,
  height = 12,
  device = "pdf"
)
