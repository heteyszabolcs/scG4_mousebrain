print("Load R packages")
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
  library("viridis")
  library("argparse")
  library("ggpubr")
  library("matrixStats")
  library("ComplexHeatmap")
  library("circlize")
})

set.seed(5)

# create parser object
parser = ArgumentParser()

parser$add_argument("-w", "--workdir", type = "character",
                    help = "path of working dir")
parser$add_argument("-s", "--seurat_object", type = "character",
                    help = "path to processed scG4 Seurat object with clusters")
parser$add_argument("-r", "--reference", type = "character",
                    help = "select reference scRNA-Seq Seurat object")

args = parser$parse_args()

# add Seurat path
g4 = args$seurat_object
# open Seurat object
g4 = readRDS(g4)
# add working directory
workdir = args$workdir
# reference scRNA-Seq Seurat object
reference = args$reference

# make folders for outputs
system(paste0("mkdir -p ", workdir, "/plots"))
system(paste0("mkdir -p ", workdir, "/outputs"))

# G4 - scRNA data integration
# Seurat workflow
print(paste0("Seurat workflow on scRNA-Seq data"))

rna = readRDS(reference)
all.genes = rownames(rna)
rna = NormalizeData(rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna = ScaleData(rna, features = all.genes)
rna = RunPCA(rna, features = VariableFeatures(object = rna))
rna = RunUMAP(rna, dims = 1:20)

p1 = DimPlot(g4,
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
  group.by = 'cell_type',
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
rna@active.ident = factor(rna$cell_type, levels = unique(rna$cell_type))

markers = FindAllMarkers(rna)
write_tsv(markers, glue("{workdir}/outputs/scRNA-Seq-FindAllMarkers_output.tsv"))
markers.pos = markers[markers$p_val < 0.05 &
                        markers$avg_log2FC > 0.5,]                    

# export scRNA Seurat object
saveRDS(rna, glue("{workdir}/outputs/scRNA_Seq_Seurat_object.Rds"))

# set assays
DefaultAssay(g4) = "GA"
DefaultAssay(rna) = "RNA"

common.genes = intersect(rownames(rna), rownames(g4))

# anchor identification between G4 scCnT and scRNA-Seq data sets
transfer.anchors = FindTransferAnchors(
  reference = rna,
  query = g4,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

# extract anchor matrix of AnchorSet object
anchor_matrix = as_tibble(transfer.anchors@anchors)
anchor_matrix = anchor_matrix %>%
  mutate(anchor1_barcode = colnames(rna@assays$RNA@counts)[anchor_matrix$cell1]) %>%
  mutate(anchor2_barcode = colnames(g4@assays$GA@counts)[anchor_matrix$cell2])
write_tsv(anchor_matrix, glue("{workdir}/outputs/anchor_matrix.tsv"))

genes.use = VariableFeatures(rna)
refdata = GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

## Transfer data (Seurat)
# imputation - data transfer
imputation = TransferData(
  anchorset = transfer.anchors,
  refdata = refdata,
  weight.reduction = g4[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

# predict cell labels and add to metadata table
cell_class_pred = TransferData(
  anchorset = transfer.anchors,
  refdata = rna@meta.data$cell_type,
  weight.reduction = g4[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

predicted_labels = cell_class_pred@data
cell_types = rownames(predicted_labels)
predicted_labels = as.tibble(predicted_labels) %>% mutate(pred_cell_type = cell_types)
predicted_labels = predicted_labels %>% pivot_longer(., cols = "AAACGAAAGAAGCCGT-1":"TTTGTGTTCTCGCGTT-1", 
                                                     names_to = "cell_id", values_to = "pred_max_score")
predicted_labels = predicted_labels %>% group_by(cell_id) %>% dplyr::slice(which.max(pred_max_score))

g4_meta = g4@meta.data
g4_meta = g4_meta %>% tibble::rownames_to_column(., var = "barcode") %>% 
  inner_join(., predicted_labels, by = c("barcode" = "cell_id")) 
barcodes = g4_meta %>% pull(barcode)
g4_meta = g4_meta %>% dplyr::select(-barcode)
rownames(g4_meta) = barcodes
g4@meta.data = g4_meta

# export G4 object with predicted labels
saveRDS(cell_class_pred, glue("{workdir}/outputs/g4_cell_label_preds.Rds"))

# visualizing predctions
# UMAP
pred_umap = DimPlot(g4,
                    pt.size = 2,
                    label.size = 7,
                    repel = TRUE,
                    raster = TRUE,
                    group.by = "pred_cell_type") +
  scale_color_brewer(palette = "Pastel1") +
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
  glue("{workdir}/plots/predicted_labels_UMAP.pdf"),
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
  file = glue("{workdir}/plots/predicton_score_heatmap.pdf"),
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
meta = as_tibble(g4@meta.data)
meta = meta %>% mutate(cell_id = rownames(g4@meta.data))
meta = predictions_long %>% inner_join(., meta, by = "cell_id")

pred_boxplots = lapply(types, function(x) {
  meta = meta %>% dplyr::filter(types == x)
  plot = ggplot(meta,
                aes(x = seurat_clusters, y = pred_score, fill = seurat_clusters)) +
    geom_boxplot(color = "black") +
    scale_fill_brewer(palette = "Pastel1") +
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
    ) + stat_compare_means(label = "p.signif", label.y = 0.9)
  return(print(plot))
})

pred_boxplots = ggarrange(plotlist = pred_boxplots)

ggsave(
  glue("{workdir}/plots/predicton_score_boxplots.pdf"),
  plot = pred_boxplots,
  width = 10,
  height = 7,
  device = "pdf"
)

g4[['RNA']] = imputation
saveRDS(g4, glue("{workdir}/integration/outputs/G4_scRNA_integration.Rds"))

# coembedding G4 and scRNA-Seq
g4 = readRDS(glue("{workdir}/integration/outputs/G4_scRNA_integration.Rds"))

rna@meta.data$data_type = "scRNA-Seq"
g4@meta.data$data_type = "G4 scCut&Tag"

coembed = merge(x = rna, y = g4)

coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
ElbowPlot(coembed)
coembed = RunUMAP(coembed, dims = 1:15)

# functions for feature plots
create_expr_feature_plot = function(marker_gene) {
  # keep scRNA-Seq part of coembedding
  coembed.scrna = coembed[,coembed$data_type == "scRNA-Seq"]
  
  # feature plot
  plot = FeaturePlot(
    object = coembed.scrna,
    features = marker_gene,
    min.cutoff = min(coembed.scrna@assays$RNA@data[marker_gene,]),
    max.cutoff = max(coembed.scrna@assays$RNA@data[marker_gene,]),
    pt.size = 4,
    raster = TRUE,
    order = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_gradient2(low = "#edf8b1", mid = "#7fcdbb", high = "#225ea8",
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
  coembed.g4 = coembed[,coembed$data_type == "G4 scCut&Tag"]
  
  plot = FeaturePlot(
    object = coembed.g4,
    features = marker_gene,
    min.cutoff = min(coembed.g4@assays$GA@data[marker_gene,]),
    max.cutoff = max(coembed.g4@assays$GA@data[marker_gene,]),
    pt.size = 4,
    raster = TRUE,
    order = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_gradient2(low = "#fee5d9", mid = "#fcae91", high = "#cb181d",
                          midpoint = mean(c(min(coembed.g4@assays$GA@data[marker_gene,]),
                                            max(coembed.g4@assays$GA@data[marker_gene,])))) +
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

check_featureplot_warnings = function(gene) {
  if (!gene %in% rownames(coembed.g4@assays$GA@data))
  {
    return(FALSE)
  }
  tt <- tryCatch(
    create_g4_feature_plot(gene),
    error = function(e)
      e,
    warning = function(w)
      w
  )
  tt2 <- tryCatch(
    create_g4_feature_plot(gene),
    error = function(e)
      e,
    warning = function(w)
      w
  )
  
  if (is(tt, "warning")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
  
}

get_best_marker = function(cell_type) {
  print(cell_type)
  coembed.g4 = coembed[,coembed$data_type == "G4 scCut&Tag"]
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = 1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data) &
      check_featureplot_warnings(marker)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = 20, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data) &
          check_featureplot_warnings(i)) {
        return(i)
      }
    }
  }
}

get_underexpr_marker = function(cell_type) {
  print(cell_type)
  coembed.g4 = coembed[,coembed$data_type == "G4 scCut&Tag"]
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = -1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data) &
      check_featureplot_warnings(marker)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = -20, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data) &
          check_featureplot_warnings(i)) {
        return(i)
      }
    }
  }
}

cell_types = unique(coembed@meta.data$cell_type)
cell_types = cell_types[!is.na(cell_types)]

underexpr_scrna_markers = sapply(cell_types, get_underexpr_marker)
best_scrna_markers = sapply(cell_types, get_best_marker)
g4_feature_plots = lapply(best_scrna_markers, create_g4_feature_plot)
g4_feature_plots = ggarrange(plotlist = g4_feature_plots)
expr_feature_plots = lapply(best_scrna_markers, create_expr_feature_plot)
expr_feature_plots = ggarrange(plotlist = expr_feature_plots)

ggsave(
  glue("{workdir}/plots/G4_Feature_plots_scRNA-Seq_markers.pdf"),
  plot = g4_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{workdir}/plots/expr_Feature_plots_scRNA-Seq_markers.pdf"),
  plot = expr_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

g4_feature_plots = lapply(underexpr_scrna_markers, create_g4_feature_plot)
g4_feature_plots = ggarrange(plotlist = g4_feature_plots)
expr_feature_plots = lapply(underexpr_scrna_markers, create_expr_feature_plot)
expr_feature_plots = ggarrange(plotlist = expr_feature_plots)

ggsave(
  glue("{workdir}/plots/G4_Feature_plots_scRNA-Seq_neg_markers.pdf"),
  plot = g4_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{workdir}/plots/expr_Feature_plots_scRNA-Seq_neg_markers.pdf"),
  plot = expr_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

# Violin plots
violins = lapply(best_scrna_markers, function(x) {
  print(x)
  plot = VlnPlot(
    object = rna,
    features = x,
    raster = FALSE,
    log = FALSE,
    split.by = "cell_type",
    group.by = "cell_type"
  ) +
    ylab("Norm. expression level") +
    ylim(0, 6) +
    scale_fill_continuous(guide=FALSE) +
    scale_fill_manual(values = rep("#a6bddb", 8)) +
    xlab(" ") +
    guides(legend = "none", fill = "none") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 11, color = "black")
    )
  print(plot)
  return(plot)
})

expr_violins = ggarrange(plotlist = violins)
expr_violins

g4_violins = lapply(best_scrna_markers, function(x)
{
  print(x)
  plot = VlnPlot(
    object = g4,
    features = x,
    raster = TRUE,
    log = FALSE,
    split.by = "seurat_clusters",
    group.by = "seurat_clusters"
  ) +
    ylab("Norm. G4 score") +
    ylim(0, 2) +
    scale_fill_continuous(guide=FALSE) +
    scale_fill_manual(values = rep("#fcae91", 8)) +
    xlab(" ") +
    guides(legend = "none", fill = "none") +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 11, color = "black")
    )
  print(plot)
  return(plot)
})

g4_violins = ggarrange(plotlist = g4_violins)
g4_violins

ggsave(
  glue("{workdir}/plots/G4_Violin_plots_scRNA-Seq_markers.pdf"),
  plot = g4_violins,
  width = 8,
  height = 11,
  device = "pdf"
)

ggsave(
  glue("{workdir}/plots/expr_Violin_plots_scRNA-Seq_markers.pdf"),
  plot = expr_violins,
  width = 8,
  height = 11,
  device = "pdf"
)

# heatmap of expression markers
hm = DoHeatmap(rna, features = best_scrna_markers,
               group.by = "cell_type", raster = FALSE, 
               group.colors = brewer.pal(8,"Pastel1"), angle = 90) +
  scale_fill_gradientn(colors = c("#421654", "#458f8a", "#f0e527")) +
  theme(
    axis.text.y = element_text(size = 18, color = "black")
  ) + NoLegend()

ggsave(
  glue("{workdir}/plots/expr_heatmap_scRNA-Seq_markers.pdf"),
  plot = hm,
  width = 11,
  height = 10
)

## UMAP dimension plots
DimPlot(rna)
coembed_cells = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'cell_type',
  repel = TRUE,
  na.value = "grey50",
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

coembed_cells2 = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'cell_type',
  repel = TRUE,
  label = TRUE,
  na.value = "grey50",
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
  ) +
  NoLegend() + NoAxes()

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_celltypes.pdf"),
  plot = coembed_cells2,
  width = 6,
  height = 6,
  device = "pdf"
)

coembed_clusters = DimPlot(
  coembed,
  pt.size = 2,
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
  pt.size = 2,
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

coembed_predcelltype = DimPlot(
  coembed,
  pt.size = 2,
  label.size = 7,
  group.by = 'pred_cell_type',
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

cells = coembed@meta.data %>% dplyr::filter(data_type == "G4 scCut&Tag") %>% rownames
coembed_predscore = FeaturePlot(
  coembed,
  pt.size = 2,
  cells = cells,
  label.size = 7,features = 'pred_max_score',
  raster = TRUE
) +
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
coembed_ps = coembed_cells + coembed_experiments + coembed_predcelltype + coembed_predscore
coembed_ps

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_coembeds.pdf"),
  plot = coembed_ps,
  width = 12,
  height = 8,
  device = "pdf"
)

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_integration_datatypes.pdf"),
  plot = coembed_experiments,
  width = 6,
  height = 6,
  device = "pdf"
)

umaps = plot_grid(p1,
                  p2,
                  coembed_experiments,
                  coembed_cells,
                  ncol = 2,
                  nrow = 2)

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_UMAPs.png"),
  plot = umaps,
  width = 12,
  height = 12,
  dpi = 300,
)

ggsave(
  glue("{workdir}/plots/Seurat_scRNA-Seq_UMAPs.pdf"),
  plot = umaps,
  width = 12,
  height = 12,
  device = "pdf"
)

# coenrichment analysis of G4 and scRNA-Seq signals
# marker gene statistics
# print("Coenrichment analysis on positive scRNA-Seq markers")
# 
# types = rownames(cell_class_pred@data)
# cell_class_pred = as_tibble(cell_class_pred@data)
# cell_class_pred = cell_class_pred %>% mutate(type = types) %>% dplyr::select(type, everything())
# 
# coembed.g4 = coembed[,coembed$data_type == "G4 scCut&Tag"]
# coembed.scrna = coembed[,coembed$data_type == "scRNA-Seq (Marques et al.)"]
# 
# max_ga = sapply(best_scrna_markers, function(x) max(coembed.g4@assays$GA[x,]))
# max_expr = sapply(best_scrna_markers, function(x) max(coembed.scrna@assays$RNA[x,]))
# 
# # set minimum GA threshold above that the cells will be labeled
# ga_thr = c(0.5, 0.5, 0.5, 0.0, 0.0, 0.0)
# marker_stat = tibble(marker = best_scrna_markers, max_GA = max_ga, 
#                      max_expression = max_expr, GA_thr = ga_thr)

# visualizations
# vis_dim_plot = function(gene) {
#   print(gene)
#   ga_thr = unname(marker_stat[which(marker_stat$marker == gene),]$GA_thr)
# 
#   get_cell_ids = coembed.g4@assays$GA@data[gene,]
#   get_cell_ids = names(get_cell_ids[get_cell_ids > ga_thr])
# 
#   predictions = cell_class_pred %>%
#     pivot_longer(
#       .,
#       cols = c(colnames(cell_class_pred)[2]:colnames(cell_class_pred)[dim(cell_class_pred)[2]]),
#       names_to = "cell_id",
#       values_to = "pred_score"
#     ) %>%
#     dplyr::filter(!type == "max")
# 
#   plot = DimPlot(
#     object = coembed.g4,
#     pt.size = 2,
#     cells.highlight = get_cell_ids,
#     cols.highlight = "red",
#     cols = "gray",
#     #rder = TRUE,
#     #raster = TRUE
#   ) + NoAxes() + ggtitle(gene) + NoLegend()
#   return(print(plot))
# }
# 
# dims = lapply(best_scrna_markers, vis_dim_plot)
# 
# for(i in seq(length(marker_stat$marker))) {
#   ggsave(
#     glue(
#       "{workdir}/plots/coenrich_cluster_anal_dimplot_{marker_stat$marker[i]}.pdf"
#     ),
#     plot = dims[i][[1]],
#     width = 4,
#     height = 3,
#     device = "pdf"
#   )
# }
# 
# vis_pred_boxplot = function(gene) {
#   print(gene)
# 
#   get_cell_ids = coembed.g4@assays$GA@data[gene,]
#   get_cell_ids = names(get_cell_ids[get_cell_ids > ga_thr])
# 
#   predictions = cell_class_pred %>%
#     pivot_longer(
#       .,
#       cols = c(colnames(cell_class_pred)[2]:colnames(cell_class_pred)[dim(cell_class_pred)[2]]),
#       names_to = "cell_id",
#       values_to = "pred_score"
#     ) %>%
#     dplyr::filter(!type == "max")
# 
#   predictions = predictions %>% dplyr::filter(cell_id %in% get_cell_ids)
# 
#   plot = ggplot(predictions,
#                 aes(x = type, y = pred_score, fill = type)) +
#     geom_boxplot(color = "black") +
#     scale_fill_brewer(palette = "Pastel1") +
#     ylim(0, 1) +
#     labs(
#       title = "",
#       x = "",
#       y = "prediction score",
#       fill = ""
#     ) +
#     theme_classic() +
#     ggtitle(gene) +
#     guides(fill = "none") +
#     theme(
#       text = element_text(size = 9),
#       plot.title = element_text(size = 15, face = "bold"),
#       axis.text.x = element_text(size = 12, color = "black"),
#       axis.text.y = element_text(size = 12, color = "black"),
#       axis.title = element_text(size = 20, color = "black")
#     ) + stat_compare_means(label = "p.signif", label.y = 0.9)
# 
#   return(print(plot))
# }
# 
# boxplots = lapply(best_scrna_markers, vis_pred_boxplot)
# boxplots[1]
# 
# for(i in seq(length(marker_stat$marker))) {
#   ggsave(
#     glue(
#       "{workdir}/plots/coenrich_cluster_anal_pred_box_{marker_stat$marker[i]}.pdf"
#     ),
#     plot = boxplots[i][[1]],
#     width = 6,
#     height = 4,
#     device = "pdf"
#   )
# }