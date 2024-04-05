if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "Seurat",
  "glue",
  "tidyverse",
  "data.table",
  "ggplot2",
  "ggpubr",
  "anndata",
  "zellkonverter",
  "RColorBrewer",
  "viridis"
)

# result folder
result_folder = "../results/scBridge/"

# Seurat integrated object
g4 = readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/G4_Marques_scRNA_integration.Rds")

# read anndata objects (scBridge outputs)
scbr_g4 = readH5AD("../results/scBridge/output/scCutTag_gene_activity_scores-integrated.h5ad")
scbr_g4 = as.Seurat(scbr_g4, counts = "X", data = NULL)
scbr_rna = readH5AD("../results/scBridge/output/Marques_scRNA-Seq-integrated.h5ad")
scbr_rna = as.Seurat(scbr_rna, counts = "X", data = NULL)
comb = readH5AD("../results/scBridge/output/combined.h5ad")
comb = as.Seurat(comb, counts = "X", data = NULL)

# Seurat UMAPs
set3 = brewer.pal(6, "Set3")
cols = c('COP'=set3[1],'MFOL'=set3[2],'MOL'=set3[3],'NFOL'=set3[4],
         'OPC'=set3[5],'PPR'=set3[6])

# cell types
scbr_rna_pred = DimPlot(
  object = scbr_rna,
  group.by = "CellType",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  cols = cols,
  order = c('PPR', 'MFOL', 'NFOL', 'COP', 'OPC', 'MOL')
) +
  xlim(-12, 25) +
  ylim(-15, 15) +
  ggtitle("Cell type (scRNA-Seq)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

scbr_g4@meta.data = scbr_g4@meta.data %>% 
  mutate(Prediction = as.character(Prediction)) %>% 
  mutate(Prediction = ifelse(str_detect(Prediction, "Novel"), "unreliable", Prediction)) %>% 
  mutate(Prediction = as.factor(Prediction))

cols = c('COP'=set3[1],'MFOL'=set3[2],'MOL'=set3[3],'NFOL'=set3[4],
         'OPC'=set3[5],'PPR'=set3[6],
         'unreliable'='#f0f0f0')

# scBridge prediction score
scbr_g4_pred = DimPlot(
  object = scbr_g4,
  group.by = "Prediction",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  cols = cols,
  order = c('PPR', 'MFOL', 'NFOL', 'COP', 'OPC', 'MOL', 'unreliable')
) +
  xlim(-12, 25) +
  ylim(-15, 15) +
  ggtitle("Prediction") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

# scBridge reliability
rel = FeaturePlot(object = scbr_g4, features = 'Reliability') +
  scale_color_viridis() +
  xlim(-12, 25) +
  ylim(-15, 15) +
  ggtitle("Reliability") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

g4@meta.data = g4@meta.data %>% rownames_to_column(., var = "cell_id")
scbr_g4@meta.data = scbr_g4@meta.data %>% rownames_to_column(., var = "cell_id")
scbr_g4@meta.data = scbr_g4@meta.data %>% inner_join(., g4@meta.data, by = "cell_id")
rownames(scbr_g4@meta.data) = scbr_g4@meta.data$cell_id

scbr_g4@meta.data = scbr_g4@meta.data %>%
  dplyr::select(
    cell_id,
    Domain,
    Prediction,
    Reliability,
    Seurat_prediction = pred_cell_type,
    Seurat_pred_score = pred_max_score
  )

seurat_pred_score = FeaturePlot(object = scbr_g4, features = 'Seurat_pred_score') +
  scale_color_viridis() +
  xlim(-12, 25) +
  ylim(-15, 15) +
  ggtitle("Seurat prediction score") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
seurat_pred_score

cols = c('COP'=set3[1],'MFOL'=set3[2],'MOL'=set3[3],'NFOL'=set3[4],
         'OPC'=set3[5],'PPR'=set3[6])

seurat_pred =DimPlot(
  object = scbr_g4,
  group.by = "Seurat_prediction",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  cols = cols,
  order = c('PPR', 'MFOL', 'NFOL', 'COP', 'OPC', 'MOL')
) +
  xlim(-12, 25) +
  ylim(-15, 15) +
  ggtitle("Seurat prediction") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
seurat_pred

scbr_g4@meta.data = scbr_g4@meta.data %>% mutate(MOL_status = ifelse(Prediction == "MOL", "MOL", "non-MOL"))
cols = c(
  "non-MOL" = '#f0f0f0',
  "MOL" = "#de2d26")

mol = DimPlot(
  object = scbr_g4,
  group.by = "MOL_status",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  order = c("MOL", "non-MOL")
) +
  xlim(-12, 25) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols) +
  ggtitle("") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
mol

ggsave(
  glue("{result_folder}Seurat_predMOL-UMAP.pdf"),
  plot = mol,
  width = 8,
  height = 8,
  device = "pdf"
)

# experimental system (modality)
comb@meta.data = comb@meta.data %>% 
  mutate(Domain = as.character(Domain)) %>% 
  mutate(Domain = ifelse(str_detect(Domain, "scCutTag_"), "G4 scCut&Tag", 
                             Domain)) %>% 
  mutate(Domain = ifelse(str_detect(Domain, "Marques_"), "scRNA-Seq", 
                             Domain)) %>% 
  mutate(Domain = as.character(Domain)) 

domain = DimPlot(
  object = comb,
  group.by = "Domain",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
) +
  xlim(-12, 25) +
  ylim(-15, 15) +
  ggtitle("Modality") +
  scale_fill_manual(values = c("#fc9272", "#a6bddb")) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

umaps = ggarrange(
  plotlist = list(scbr_rna_pred, scbr_g4_pred, rel, domain, seurat_pred_score, seurat_pred),
  ncol = 2,
  nrow = 3
)
ggsave(
  glue("{result_folder}Seurat_UMAPs.pdf"),
  plot = umaps,
  width = 18,
  height = 18,
  device = "pdf"
)

featureplots = list()
genes = c("Otol1", "Trim34b", "Mmp20", "Serpinb3b", "Glra3", "Trim12a")
for(gene in genes) {
  plot = FeaturePlot(object = scbr_g4, features = gene, order = TRUE, reduction = "X_umap", pt.size = 2) +
    #scale_color_gradient2(low = "#f0f0f0", mid = "#f0f0f0", high = "red") +
    scale_color_viridis() +
    xlim(-12, 25) +
    ylim(-15, 15) +
    ggtitle(gene) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    )
  featureplots[[gene]] = plot
}

mol_featureplotes = ggarrange(
  plotlist = featureplots,
  ncol = 2,
  nrow = 3
)

ggsave(
  glue("{result_folder}Seurat_UMAPs-MOL_spec_G4_markers.pdf"),
  plot = mol_featureplotes,
  width = 18,
  height = 18,
  device = "pdf"
)


# scBridge output tables
rel = fread("../results/scBridge/output/scbridge_reliability.csv")
pred = fread("../results/scBridge/output/scbridge_predictions.csv", header = TRUE)

meta = g4@meta.data %>% 
  left_join(., rel, c("cell_id" = "V1")) %>% rename(scBridge_reliability = "Reliability") %>% 
  left_join(., pred, c("cell_id" = "V1")) %>% rename(scBridge_prediction = "Prediction")

rownames(meta) = rownames(g4@meta.data) 
g4@meta.data = meta

both_preds_mol = meta %>% filter(scBridge_prediction == "MOL" & pred_cell_type == "MOL")

cor(g4@meta.data$pred_max_score, g4@meta.data$scBridge_reliability, method = "spearman")
ggplot(g4@meta.data, aes(x = scBridge_reliability, y = pred_max_score, color = scBridge_prediction)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Pastel1") +
  labs(
    title = " ",
    x = "scBridge reliability",
    y = "Seurat prediction score",
    color = " "
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"))

types = unique(meta$pred_cell_type)
pred_boxplots = lapply(types, function(x) {
  meta = meta %>% dplyr::filter(scBridge_prediction == x)
  plot = ggplot(meta,
                aes(x = seurat_clusters, y = scBridge_reliability, fill = seurat_clusters)) +
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
pred_boxplots

DimPlot(
  g4,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = TRUE,
  group.by = "scBridge_prediction"
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

FeaturePlot(
  object = g4,
  features = "scBridge_reliability",
  min.cutoff = 0,
  max.cutoff = max(g4@meta.data$scBridge_reliability),
  pt.size = 2,
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-10, 10) +
  scale_color_gradient2(low = "white", mid = "orange", high = "red",
                        midpoint = mean(c(min(g4@meta.data$scBridge_reliability),
                                          max(g4@meta.data$scBridge_reliability)))) +
  theme(
    legend.position = 'bottom',
    text = element_text(size = 15),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

FeaturePlot(
  object = g4,
  features = "scBridge_reliability",
  cells = rownames(both_preds_mol),
  min.cutoff = 0,
  max.cutoff = max(g4@meta.data$scBridge_reliability),
  pt.size = 2,
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-10, 10) +
  scale_color_gradient2(low = "white", mid = "orange", high = "red",
                        midpoint = mean(c(min(g4@meta.data$scBridge_reliability),
                                          max(g4@meta.data$scBridge_reliability)))) +
  theme(
    legend.position = 'bottom',
    text = element_text(size = 15),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()
