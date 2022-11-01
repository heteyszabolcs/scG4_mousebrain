library("ggplot2")
library("glue")
library("GenomicRanges")
library("Seurat")
library("cowplot")

set.seed(5)

# path to result folder
result_folder = "../results/Seurat/"

marques_scrna <-
  readRDS(file = "../results/Seurat/scRNASeq_GSE75330.rds")
g4 = readRDS(file = "../results/Seurat/callpeaks_unsorted/unsorted.Rds")
g4_sorted = readRDS(file = "../results/Seurat/GFP_sorted.Rds")

# UMAPs
p1 =
  DimPlot(marques_scrna, label = TRUE, pt.size = 2) + 
  scale_color_brewer(palette = "Set3") +
  NoLegend() + 
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("scRNA-Seq (Marques et al.)")
p1
p2 =
  DimPlot(g4, label = TRUE, pt.size = 2) + 
  scale_color_brewer(palette = "Set3") +
  NoLegend() + 
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("unsorted G4 scCutnTag")
p2

cols = c(
  "0" = "#addd8e",
  "1" = "#bdbdbd",
  "2" = "#addd8e",
  "6" = "#bdbdbd",
  "3" = "#addd8e",
  "5" = "#addd8e",
  "4" = "#bdbdbd"
)

colored_by_sorting =
  DimPlot(
    g4,
    label = FALSE,
    pt.size = 1.5,
  ) +
  scale_colour_manual(values = cols, breaks = c("0", "4"), labels = c("GFP+", "unsorted")) +
  ggtitle("Oligodendrocyte population") +
  xlim(-10, 10) + ylim(-10, 10)
colored_by_sorting

ggsave(
  glue(
    "{result_folder}UMAP_colored_by_sorting.png"
  ),
  plot = colored_by_sorting,
  width = 10,
  height = 10,
  dpi = 500,
)

sorted = DimPlot(g4_sorted, label = TRUE, pt.size = 2) + 
  scale_color_brewer(palette = "Set3") +
  NoLegend() + 
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("sorted G4 scCutnTag")
sorted

umaps = plot_grid(p1, p2, colored_by_sorting, sorted, ncol = 2, nrow = 2)
umaps

ggsave(
  glue(
    "{result_folder}UMAP_Summary.png"
  ),
  plot = umaps,
  width = 10,
  height = 10,
  dpi = 500,
)


## transfer learning, query: unsorted G4 set
DefaultAssay(marques_scrna) <- "RNA"
DefaultAssay(g4) <- "GA"

# focus on common genes
common.genes <- intersect(rownames(marques_scrna), rownames(g4))

# g4[["bins_25000"]] = NULL
# g4[["peaks"]] = NULL
# g4[["PA"]] = NULL
# g4[["bins_5000"]] = NULL

marques_scrna = subset(x = marques_scrna, downsample = 500)
g4 = subset(x = g4, downsample = 500)

# Transfer categorical data across single-cell datasets
transfer.anchors <- FindTransferAnchors(
  reference = marques_scrna,
  query = g4,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

genes.use <- VariableFeatures(marques_scrna)
refdata <- GetAssayData(marques_scrna, assay = "RNA", slot = "data")

imputation <-
  TransferData(
    anchorset = transfer.anchors,
    refdata = marques_scrna$cell_type,
    weight.reduction = g4[["lsi"]],
    dims = 1:50
  )
g4 = AddMetaData(g4, metadata = imputation)

# UMAP of transfered data
imputed_umap = DimPlot(
  g4,
  reduction = "umap",
  group.by = "predicted.id",
  pt.size = 3,
  label = FALSE) + 
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("Label imputation - unsorted")
imputed_umap

saveRDS(g4,
        glue(
          "../results/Seurat/scRNASeq_GSE75330-unsorted_G4_integrated.Rds"
        ))

## transfer learning, query: sorted G4 set
DefaultAssay(marques_scrna) <- "RNA"
DefaultAssay(g4_sorted) <- "GA"

# focus on common genes
common.genes <- intersect(rownames(marques_scrna), rownames(g4_sorted))

marques_scrna = subset(x = marques_scrna, downsample = 500)
g4_sorted = subset(x = g4_sorted, downsample = 500)

# Transfer categorical data across single-cell datasets
transfer.anchors <- FindTransferAnchors(
  reference = marques_scrna,
  query = g4_sorted,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'RNA',
  k.filter = NA,
  features = common.genes
)

genes.use <- VariableFeatures(marques_scrna)
refdata <- GetAssayData(marques_scrna, assay = "RNA", slot = "data")

imputation <-
  TransferData(
    anchorset = transfer.anchors,
    refdata = marques_scrna$cell_type,
    weight.reduction = g4_sorted[["lsi"]],
    dims = 1:50
  )
g4_sorted = AddMetaData(g4_sorted, metadata = imputation)

# UMAP of transfered data
imputed_umap_sorted = DimPlot(
  g4_sorted,
  reduction = "umap",
  group.by = "predicted.id",
  pt.size = 3,
  label = FALSE) + 
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("Label imputation - GFP+ population")
imputed_umap_sorted

saveRDS(g4_sorted,
        glue(
          "../results/Seurat/scRNASeq_GSE75330-sorted_G4_integrated.Rds"
        ))

imputed_umaps = plot_grid(imputed_umap, imputed_umap_sorted, colored_by_sorting, sorted, ncol = 2, nrow = 1)
imputed_umaps

ggsave(
  glue(
    "{result_folder}scRNASeq_GSE75330-imputed_UMAPs.png"
  ),
  plot = imputed_umaps,
  width = 10,
  height = 5,
  dpi = 500,
)

# retrieve prediction scores
pred_scores_1 = tibble(prediction_score = g4$prediction.score.COP, cell_type = "COP")
pred_scores_2 = tibble(prediction_score = g4$prediction.score.MFOL.MOL, cell_type = "MFOL/MOL")
pred_scores_3 = tibble(prediction_score = g4$prediction.score.OPC, cell_type = "OPC")
pred_scores_4 = tibble(prediction_score = g4$prediction.score.MFOL, cell_type = "MFOL")
pred_scores_5 = tibble(prediction_score = g4$prediction.score.VLMC, cell_type = "VLMC")
pred_scores_6 = tibble(prediction_score = g4$prediction.score.MOL, cell_type = "MOL")
pred_scores = rbind(pred_scores_1, pred_scores_2, pred_scores_3, pred_scores_4, pred_scores_5, pred_scores_6)

pred_scores_1_sorted = tibble(prediction_score = g4_sorted$prediction.score.COP, cell_type = "COP")
pred_scores_2_sorted = tibble(prediction_score = g4_sorted$prediction.score.MFOL.MOL, cell_type = "MFOL/MOL")
pred_scores_3_sorted = tibble(prediction_score = g4_sorted$prediction.score.OPC, cell_type = "OPC")
pred_scores_4_sorted = tibble(prediction_score = g4_sorted$prediction.score.MFOL, cell_type = "MFOL")
pred_scores_5_sorted = tibble(prediction_score = g4_sorted$prediction.score.VLMC, cell_type = "VLMC")
pred_scores_6_sorted = tibble(prediction_score = g4_sorted$prediction.score.MOL, cell_type = "MOL")
pred_scores_sorted = rbind(pred_scores_1_sorted, pred_scores_2_sorted, pred_scores_3_sorted, pred_scores_4_sorted, 
                           pred_scores_5_sorted, pred_scores_6_sorted)

box = ggplot(pred_scores, aes(
  x = cell_type,
  y = prediction_score
)) +
  geom_boxplot(fill = "#fc9272") +
  theme_minimal() +
  ylim(0, 1) +
  labs(title = "Prediction scores - unsorted",
       fill = "") +
  xlab(label = "imputed cell type") +
  ylab(label = "prediction score") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(size = 14)
  ) 
box

box_sorted = ggplot(pred_scores_sorted, aes(
  x = cell_type,
  y = prediction_score
)) +
  geom_boxplot(fill = "#fc9272") +
  theme_minimal() +
  ylim(0, 1) +
  labs(title = "Prediction scores - GFP+ population",
       fill = "") +
  xlab(label = "imputed cell type") +
  ylab(label = "prediction score") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(size = 14)
  ) 
box_sorted

boxplots = plot_grid(box, box_sorted, ncol = 2, nrow = 1)
boxplots

# saving grid
ggsave(
  glue(
    "{result_folder}scRNASeq_GSE75330_G4_transfer_predscores.png"
  ),
  plot = boxplots,
  width = 10,
  height = 5,
  dpi = 500
)


