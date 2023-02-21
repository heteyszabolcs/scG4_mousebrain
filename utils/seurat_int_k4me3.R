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
  library("ggpubr")
  library("ComplexHeatmap")
  library("circlize")
})

# scCut&Tag - scRNA data integration
# main source: 
# https://github.com/Castelo-Branco-lab/scCut-Tag_2020/blob/master/notebooks/integration/integration_H3K4me3_marques.Rmd

# path to result folder
result_folder = "../results/Seurat/"

# Seurat object of Bartosovic et al. (GSE157637)
k4me3 = readRDS(file = "../data/GSE157637/H3K4me3_seurat_object.Rds")

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

p1 = DimPlot(k4me3,
             pt.size = 2,
             label.size = 7,
             repel = TRUE,
             raster = TRUE,
             group.by = "seurat_clusters") +
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
p1

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

DefaultAssay(k4me3) = "GA"
DefaultAssay(rna) = "RNA"

common.genes = intersect(rownames(rna), rownames(k4me3))

# anchor identification between G4 scCnT and scRNA-Seq data sets
transfer.anchors = FindTransferAnchors(
  reference = rna,
  query = k4me3,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

genes.use = VariableFeatures(rna)
refdata = GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

## Transfer data (Seurat)
# imputation - data transfer
imputation = TransferData(
  anchorset = transfer.anchors,
  refdata = refdata,
  weight.reduction = k4me3[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

# predict cell labels
cell_class_pred = TransferData(
  anchorset = transfer.anchors,
  refdata = rna@meta.data$merged_cell_class,
  weight.reduction = k4me3[["lsi"]],
  dims = 1:50,
  prediction.assay = TRUE
)

saveRDS(cell_class_pred, glue("{result_folder}H3K4me3_cell_label_preds.Rds"))

predictions = as.matrix(cell_class_pred@data)
predictions = predictions[rownames(predictions) != "max",]
col_fun = colorRamp2(c(0, 0.5, 1), c("#421654", "#458f8a", "#f0e527"))
pdf(
  file = glue("{result_folder}H3K4me3_predicton_score_heatmap.pdf"),
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

types = rownames(predictions)
predictions = as.data.frame(predictions)
predictions$types = types

predictions_long = predictions %>% pivot_longer(., cols = c("H3K4me3_N1_AAACGAAAGGCTCCTG-1":"H3K4me3_N4_TTTGTGTTCGTGGCGT-1"),
                                                names_to = "cell_id", values_to = "pred_score")
meta = as_tibble(k4me3@meta.data)

meta = meta %>% mutate(cell_id = rownames(k4me3@meta.data))
meta = predictions_long %>% inner_join(., meta, by = "cell_id")
h3k4me3_high_preds = meta %>% dplyr::filter(pred_score >= 0.75)

write_tsv(h3k4me3_high_preds, glue("{result_folder}H3K4me3_highly_pred_celltypes.tsv"))

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
  glue("{result_folder}H3K4me3_predicton_score_boxplots.pdf"),
  plot = pred_boxplots,
  width = 10,
  height = 7,
  device = "pdf"
)

k4me3[['RNA']] = imputation
saveRDS(k4me3, glue("{result_folder}H3K4me3_int_Marques.Rds"))
k4me3 = readRDS(glue("{result_folder}H3K4me3_int_Marques.Rds"))

#markers = FindAllMarkers(rna)
markers = read_tsv("../results/Seurat/scRNA-Seq_Marques_et_al-FindAllMarkers_output.tsv")
markers.pos = markers[markers$p_val < 0.05 &
                        markers$avg_log2FC > 0.5,]

rna@meta.data$data_type = "scRNA-Seq (Marques et al.)"
k4me3@meta.data$data_type = "H3K4me3 scCut&Tag"

coembed = merge(x = rna, y = k4me3)

coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
ElbowPlot(coembed)
coembed = RunUMAP(coembed, dims = 1:11)

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

create_k4me3_feature_plot = function(marker_gene) {
  # keep G4 scCut&Tag part of coembedding
  coembed.g4 = coembed[,coembed$data_type == "H3K4me3 scCut&Tag"]
  
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
  coembed.g4 = coembed[,coembed$data_type == "H3K4me3 scCut&Tag"]
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
  coembed.g4 = coembed[,coembed$data_type == "H3K4me3 scCut&Tag"]
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
h3k4me3_feature_plots = lapply(best_marques_markers, create_k4me3_feature_plot)
h3k4me3_feature_plots = ggarrange(plotlist = h3k4me3_feature_plots)
expr_feature_plots = lapply(best_marques_markers, create_expr_feature_plot)
expr_feature_plots = ggarrange(plotlist = expr_feature_plots)

ggsave(
  glue("{result_folder}H3K4me3_Feature_plots_Marques_markers-K4me3_coemb.pdf"),
  plot = h3k4me3_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{result_folder}expr_Feature_plots_Marques_markers-K4me3_coemb.pdf"),
  plot = expr_feature_plots,
  width = 10,
  height = 7,
  device = "pdf"
)

## UMAP dimension plots
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
coembed_cells2

