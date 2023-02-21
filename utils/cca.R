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
  library("ggpubr")
  library("matrixStats")
  library("ComplexHeatmap")
  library("circlize")
})

set.seed(5)

# path to result folder
cellranger_folder = "../data/CellRanger/GFP_sorted/"
# path to result folder
result_folder = "../results/Seurat/callpeaks_GFPsorted/"

# G4 data
g4 = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")

# scRNA data
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


rna@meta.data$data_type = "scRNA-Seq (Marques et al.)"
g4@meta.data$data_type = "sorted G4 scCut&Tag"

merged_seurat = merge(
  rna,
  y = g4,
  add.cell.ids = c("scRNA-Seq (Marques et al.)", "sorted G4 scCut&Tag"),
  project = "int"
)

merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat <- NormalizeData(object = merged_seurat)
merged_seurat <- FindVariableFeatures(object = merged_seurat)
merged_seurat <- ScaleData(object = merged_seurat)
merged_seurat <- RunPCA(object = merged_seurat)
ElbowPlot(merged_seurat)
merged_seurat <- FindNeighbors(object = merged_seurat, dims = 1:10)
merged_seurat <- FindClusters(object = merged_seurat)
merged_seurat <- RunUMAP(object = merged_seurat, dims = 1:10)

x = merged_seurat@meta.data

# plot
p1 <- DimPlot(merged_seurat, reduction = 'umap', group.by = 'data_type')
p1
p2 <- DimPlot(merged_seurat, reduction = 'umap', group.by = 'cell_class')
p2



obj.list <- SplitObject(merged_seurat, split.by = 'data_type')
# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 500)
rm(g4); rm(rna)

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

# for(i in 1:length(obj.list)){
#   print(obj.list[[i]])
#   obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
#   obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
#                                         selection.method = "vst", nfeatures = 100)
# }
# rm(coembed)

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 500)
# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)


# Feature plots
create_expr_feature_plot = function(marker_gene) {
  # keep G4 scCut&Tag part of coembedding
  coembed.scrna = coembed[,coembed$data_type == "scRNA-Seq (Marques et al.)"]
  
  # feature plot
  plot = FeaturePlot(
    object = coembed.scrna,
    features = marker_gene,
    min.cutoff = min(coembed.scrna@assays$RNA@data[marker_gene,]),
    max.cutoff = max(coembed.scrna@assays$RNA@data[marker_gene,]),
    raster = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_gradient2(low = "#0d0a1e", mid = "#c44a46", high = "#eef07a",
                          midpoint = mean(c(min(coembed.scrna@assays$RNA@data[marker_gene,]),
                                            max(coembed.scrna@assays$RNA@data[marker_gene,])))) +
    theme(
      legend.position = 'bottom',
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoAxes()
  return(print(plot))
}

create_g4_feature_plot = function(marker_gene) {
  # keep G4 scCut&Tag part of coembedding
  coembed.g4 = coembed[,coembed$data_type == "sorted G4 scCut&Tag"]
  
  # feature plot
  plot = FeaturePlot(
    object = coembed.g4,
    features = marker_gene,
    min.cutoff = min(coembed.g4@assays$GA@data[marker_gene,]),
    max.cutoff = max(coembed.g4@assays$GA@data[marker_gene,]),
    raster = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_viridis_c() +
    theme(
      legend.position = 'bottom',
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoAxes()
  return(print(plot))
}

# pull gene symbol with the highest avg_log2FC score per cluster 
get_best_marker = function(cell_type) {
  coembed.g4 = coembed[,coembed$data_type == "sorted G4 scCut&Tag"]
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = 1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data)) {
        return(i)
      }
    }
  }
}

get_underexpr_marker = function(cell_type) {
  coembed.g4 = coembed[,coembed$data_type == "sorted G4 scCut&Tag"]
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = -1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data)) {
        return(i)
      }
    }
  }
}

cell_types = unique(coembed@meta.data$cell_class)
cell_types = cell_types[!is.na(cell_types)]

underexpr_marques_markers = sapply(cell_types, get_underexpr_marker)
best_marques_markers = sapply(cell_types, get_best_marker)
g4_feature_plots = lapply(best_marques_markers, create_g4_feature_plot)
g4_feature_plots = ggarrange(plotlist = g4_feature_plots)
expr_feature_plots = lapply(best_marques_markers, create_expr_feature_plot)
expr_feature_plots = ggarrange(plotlist = expr_feature_plots)

ggsave(
  glue("{result_folder}G4_Feature_plots_Marques_markers.pdf"),
  plot = g4_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{result_folder}expr_Feature_plots_Marques_markers.pdf"),
  plot = expr_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

g4_feature_plots = lapply(underexpr_marques_markers, create_g4_feature_plot)
g4_feature_plots = ggarrange(plotlist = g4_feature_plots)
expr_feature_plots = lapply(underexpr_marques_markers, create_expr_feature_plot)
expr_feature_plots = ggarrange(plotlist = expr_feature_plots)

ggsave(
  glue("{result_folder}G4_Feature_plots_Marques_neg_markers.pdf"),
  plot = g4_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{result_folder}expr_Feature_plots_Marques_neg_markers.pdf"),
  plot = expr_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

# Violin plots
violins = lapply(best_marques_markers, function(x)
  VlnPlot(
    object = rna,
    features = x,
    raster = TRUE,
    log = FALSE,
    split.by = "merged_cell_class",
    group.by = "merged_cell_class"
  ) +
    ylim(0, 5) +
    scale_fill_manual(values = rep("#a6bddb", 6)) +
    xlab(" ") +
    guides(legend = "none", fill = "none"))

expr_violins = ggarrange(plotlist = violins)
expr_violins

g4_violins = lapply(best_marques_markers, function(x)
  VlnPlot(
    object = unsorted_res0.1[, unsorted_res0.1$seurat_clusters != "2"],
    features = x,
    raster = TRUE,
    log = FALSE,
    split.by = "seurat_clusters",
    group.by = "seurat_clusters"
  ) +
    scale_fill_manual(values = c("#addd8e", "#bdbdbd")) +
    xlab(" ") +
    ylim(0, 0.00005) +
    guides(legend = "none", fill = "none") +
    stat_compare_means(label.y = 4e-05, label.x = 1.2, size = 3) +
    stat_compare_means(
      label = "p.signif",
      method = "t.test",
      ref.group = ".all.",
      label.y = 50
    ))


g4_violins = ggarrange(plotlist = g4_violins)
g4_violins

ggsave(
  glue("{result_folder}G4_Violin_plots_Marques_markers.pdf"),
  plot = g4_violins,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{result_folder}expr_Violin_plots_Marques_markers.pdf"),
  plot = expr_violins,
  width = 8,
  height = 6,
  device = "pdf"
)

# heatmap of expression markers
hm = DoHeatmap(rna, features = best_marques_markers,
               group.by = "merged_cell_class", raster = FALSE) +
  scale_fill_gradientn(colors = c("#0d0a1e", "#c44a46", "#eef07a")) +
  theme(
    axis.text.y = element_text(size = 18, color = "black")
  ) + NoLegend()

ggsave(
  glue("{result_folder}expr_heatmap_Marques_markers.pdf"),
  plot = hm,
  width = 10,
  height = 7
)

sorted_cluster2 = subset(sorted, subset = seurat_clusters == 2)
sorted_cluster2 = as.matrix(sorted_cluster2@assays$GA@counts)
most_var_cluster2 = rownames(sorted_cluster2)[order(rowMeans(sorted_cluster2), decreasing = TRUE)][1:200]
most_var_cluster2 = markers %>% dplyr::filter(gene %in% most_var_cluster2) %>% 
  dplyr::filter(avg_log2FC > 0) %>% mutate(cluster_aggr = ifelse(str_detect(cluster, "NFOL"), "NFOL", cluster)) %>% 
  mutate(cluster_aggr = ifelse(str_detect(cluster, "MFOL"), "MFOL", cluster_aggr)) %>% 
  mutate(cluster_aggr = ifelse(str_detect(cluster, "MOL"), "MOL", cluster_aggr)) %>% 
  group_by(cluster_aggr) %>% summarise(cluster_n = n())

bars = ggplot(data = most_var_cluster2, aes(y = cluster_n, x = reorder(cluster_aggr, -cluster_n))) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "GFP+ Seurat cluster 2",
       x = " ",
       y = "# of positive expr. marker") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 15, color = "black"))
bars

ggsave(
  glue("{result_folder}Seurat_cl2-quant_of_pos_markers.pdf"),
  plot = bars,
  width = 7,
  height = 5
)

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
  na.value = "white",
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
  "4" = "#bdbdbd"
  # "5" = "#bdbdbd",
  # "6" = "#bdbdbd"
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
  "4" = "#bdbdbd"
  # "5" = "#bdbdbd",
  # "6" = "#bdbdbd"
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
  "4" = "#bdbdbd"
  # "5" = "#bdbdbd",
  # "6" = "#bdbdbd"
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
  "4" = "#bdbdbd")
# "5" = "#bdbdbd",
# "6" = "#bdbdbd"


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
  "4" = "#de2d26")
# "5" = "#bdbdbd",
# "6" = "#bdbdbd"

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
  "4" = "#bdbdbd")
# "5" = "#de2d26",
# "6" = "#bdbdbd"


# coembed_cluster5 = DimPlot(
#   coembed,
#   pt.size = 1,
#   label.size = 7,
#   group.by = 'seurat_clusters',
#   repel = TRUE,
#   order = "5",
#   raster = TRUE
# ) +
#   xlim(-15, 15) +
#   ylim(-15, 15) +
#   scale_colour_manual(values = cols,
#                       breaks = "5",
#                       labels = "cluster 5") +
#   ggtitle(" ") +
#   theme(
#     text = element_text(size = 25),
#     plot.title = element_text(size = 20),
#     axis.text.x = element_text(size = 25, color = "black"),
#     axis.text.y = element_text(size = 25, color = "black")
#   ) +
#   NoAxes()
# 
# ggsave(
#   glue("{result_folder}Seurat_Marques_integration_cl5.pdf"),
#   plot = coembed_cluster5,
#   width = 10,
#   height = 10,
#   device = "pdf"
# )

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd"
  #   "5" = "#bdbdbd",
  #   "6" = "#de2d26"
)

# coembed_cluster6 = DimPlot(
#   coembed,
#   pt.size = 1,
#   label.size = 7,
#   group.by = 'seurat_clusters',
#   repel = TRUE,
#   order = "6",
#   raster = TRUE
# ) +
#   xlim(-15, 15) +
#   ylim(-15, 15) +
#   scale_colour_manual(values = cols,
#                       breaks = "6",
#                       labels = "cluster 6") +
#   ggtitle(" ") +
#   theme(
#     text = element_text(size = 25),
#     plot.title = element_text(size = 20),
#     axis.text.x = element_text(size = 25, color = "black"),
#     axis.text.y = element_text(size = 25, color = "black")
#   ) +
#   NoAxes()

# ggsave(
#   glue("{result_folder}Seurat_Marques_integration_cl6.pdf"),
#   plot = coembed_cluster6,
#   width = 10,
#   height = 10,
#   device = "pdf"
# )

clusters = plot_grid(
  coembed_cluster0,
  coembed_cluster1,
  coembed_cluster2,
  coembed_cluster3,
  coembed_cluster4,
  #coembed_cluster6,
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
  glue("{result_folder}Seurat_sorted_Marques_UMAPs.png"),
  plot = umaps,
  width = 12,
  height = 12,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_sorted_Marques_UMAPs.pdf"),
  plot = umaps,
  width = 12,
  height = 12,
  device = "pdf"
)