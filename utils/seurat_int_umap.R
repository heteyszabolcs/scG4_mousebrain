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

# G4 - scRNA data integration
# main source: 
# https://github.com/Castelo-Branco-lab/scCut-Tag_2020/blob/master/notebooks/integration/integration_H3K4me3_marques.Rmd

# path to result folder
result_folder = "../results/Seurat/callpeaks_GFPsorted/"

# Seurat object (brain unsorted / brain GFP sorted)
unsorted = readRDS(file = "../results/Seurat/callpeaks_unsorted/unsorted.Rds")
unsorted_res0.1 = readRDS(file = "../results/Seurat/callpeaks_unsorted/unsorted_res0.1.Rds")
sorted = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")

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

p1 = DimPlot(sorted,
             pt.size = 2,
             label.size = 7,
             repel = TRUE,
             raster = TRUE) +
  scale_color_brewer(palette = "Pastel1") +
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

# export scRNA Seurat object
saveRDS(rna, glue("../results/Seurat/scRNASeq_GSE75330.rds"))

#markers = FindAllMarkers(rna)
markers = read_tsv("../results/Seurat/scRNA-Seq_Marques_et_al-FindAllMarkers_output.tsv")
markers.pos = markers[markers$p_val < 0.05 &
                        markers$avg_log2FC > 0.5,]

# write_tsv(markers,
#           "../results/Seurat/scRNA-Seq_Marques_et_al-FindAllMarkers_output.tsv")

DefaultAssay(sorted) = "GA"
DefaultAssay(rna) = "RNA"

common.genes = intersect(rownames(rna), rownames(sorted))

# anchor identification between G4 scCnT and scRNA-Seq data sets
transfer.anchors = FindTransferAnchors(
  reference = rna,
  query = sorted,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

# extract anchor matrix of AnchorSet object
g4_anchor_matrix = as_tibble(transfer.anchors@anchors)
g4_anchor_matrix = g4_anchor_matrix %>%
  mutate(anchor1_barcode = colnames(rna@assays$RNA@counts)[g4_anchor_matrix$cell1]) %>%
  mutate(anchor2_barcode = colnames(sorted@assays$GA@counts)[g4_anchor_matrix$cell2])
write_tsv(g4_anchor_matrix, glue("{result_folder}sorted_anchor_matrix.tsv"))

genes.use = VariableFeatures(rna)
refdata = GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

## Transfer data (Seurat)
# imputation - data transfer
imputation = TransferData(
  anchorset = transfer.anchors,
  refdata = refdata,
  weight.reduction = sorted[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

# predict cell labels and add to metadata table
cell_class_pred = TransferData(
  anchorset = transfer.anchors,
  refdata = rna@meta.data$merged_cell_class,
  weight.reduction = sorted[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

predicted_labels = cell_class_pred@data
cell_types = rownames(predicted_labels)
predicted_labels = as.tibble(predicted_labels) %>% mutate(pred_cell_type = cell_types)
predicted_labels = predicted_labels %>% pivot_longer(., cols = "AAACGAAAGAAGCCGT-1":"TTTGTGTTCTCGCGTT-1", 
                                                     names_to = "cell_id", values_to = "pred_max_score")
predicted_labels = predicted_labels %>% group_by(cell_id) %>% dplyr::slice(which.max(pred_max_score))

sorted_meta = sorted@meta.data
sorted_meta = sorted_meta %>% mutate(rownames_to_column(., var = "barcode")) %>% 
  inner_join(., predicted_labels, by = c("barcode" = "cell_id")) 
barcodes = sorted_meta %>% pull(barcode)
sorted_meta = sorted_meta %>% dplyr::select(-barcode)
rownames(sorted_meta) = barcodes
sorted@meta.data = sorted_meta

# export G4 object with predicted labels
saveRDS(cell_class_pred, glue("{result_folder}sorted_cell_label_preds.Rds"))

# visualizing predctions
# UMAP
pred_umap = DimPlot(sorted,
             pt.size = 2,
             label.size = 7,
             repel = TRUE,
             raster = TRUE,
             group.by = "pred_cell_type") +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(" ") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) + NoAxes()
pred_umap

ggsave(
  glue("{result_folder}predicted_labels_UMAP.pdf"),
  plot = pred_umap,
  width = 13,
  height = 10,
  device = "pdf"
)

# prediction score heatmap
predictions = as.matrix(cell_class_pred@data)
predictions = predictions[rownames(predictions) != "max",]
col_fun = colorRamp2(c(0, 0.5, 1), c("#421654", "#458f8a", "#f0e527"))
pdf(
  file = glue("{result_folder}predicton_score_heatmap.pdf"),
  width = 8,
  height = 6
)
Heatmap(
  predictions,
  column_title = " ",
  row_title = " ",
  name = "pred. score",
  # row_km = 2,
  # column_km = 2,
  clustering_method_rows = "complete",
  col = col_fun,
  #rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  heatmap_width = unit(17, "cm"),
  heatmap_height = unit(6, "cm"),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 90
)
dev.off()

# prediction score boxplot
types = rownames(predictions)
predictions = as.data.frame(predictions)
predictions$types = types
predictions_long = predictions %>% pivot_longer(., cols = c("AAACGAAAGAAGCCGT-1":"TTTGTGTTCTCGCGTT-1"),
                                                names_to = "cell_id", values_to = "pred_score")
meta = as_tibble(sorted@meta.data)
meta = meta %>% mutate(cell_id = rownames(sorted@meta.data))
meta = predictions_long %>% inner_join(., meta, by = "cell_id")

pred_boxplots = lapply(types, function(x) {
  meta = meta %>% dplyr::filter(types == x)
  plot = ggplot(meta,
                aes(x = seurat_clusters, y = pred_score, fill = seurat_clusters)) +
    geom_boxplot(color = "black") +
    scale_fill_brewer(palette = "Set3") +
    ylim(0, 1) +
    labs(
      title = x,
      x = "",
      y = "",
      fill = ""
    ) +
    theme_classic() +
    guides(fill = "none") +
    theme(
      text = element_text(size = 9),
      plot.title = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black")
    )
  return(print(plot))
})

pred_boxplots = ggarrange(plotlist = pred_boxplots)

ggsave(
  glue("{result_folder}predicton_score_boxplots.pdf"),
  plot = pred_boxplots,
  width = 10,
  height = 7,
  device = "pdf"
)

sorted[['RNA']] = imputation
saveRDS(sorted, glue("{result_folder}sorted_int_Marques.Rds"))

# coembedding G4 and Marques et al. scRNA-Seq
sorted = readRDS(glue("{result_folder}sorted_int_Marques.Rds"))

rna@meta.data$data_type = "scRNA-Seq (Marques et al.)"
sorted@meta.data$data_type = "sorted G4 scCut&Tag"

coembed = merge(x = rna, y = sorted)

coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
ElbowPlot(coembed)
coembed = RunUMAP(coembed, dims = 1:15)

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
  print(cell_type)
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
  scale_color_brewer(palette = "Pastel1") +
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
  ) + NoAxes()

coembed@meta.data$seurat_clusters[is.na(coembed@meta.data$seurat_clusters)] = "scRNA-Seq"

ggsave(
  glue("{result_folder}Seurat_Marques_integration_datatypes.pdf"),
  plot = coembed_experiments,
  width = 13,
  height = 10,
  device = "pdf"
)

cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#de2d26",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd"
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

clusters = plot_grid(
  coembed_cluster0,
  coembed_cluster1,
  coembed_cluster2,
  coembed_cluster3,
  ncol = 2,
  nrow = 2
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
