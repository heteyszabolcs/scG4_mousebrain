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
g4 = readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/G4_scRNA_integration.Rds")

# read anndata objects (scBridge outputs)
scbr_g4 = readH5AD("../results/scBridge/output/scCutTag_gene_activity_scores-integrated.h5ad")
scbr_g4 = as.Seurat(scbr_g4, counts = "X", data = NULL)
scbr_rna = readH5AD("../results/scBridge/output/Bartosovic_scRNA-Seq-integrated.h5ad")
scbr_rna = as.Seurat(scbr_rna, counts = "X", data = NULL)
comb = readH5AD("../results/scBridge/output/combined.h5ad")
comb = as.Seurat(comb, counts = "X", data = NULL)

scbr_rna@meta.data = scbr_rna@meta.data %>% 
  mutate(CellType = str_replace_all(CellType, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(CellType = str_replace_all(CellType, pattern = "Oligodendrocytes", replacement = "MOL"))

scbr_g4@meta.data = scbr_g4@meta.data %>% 
  mutate(Prediction = str_replace_all(Prediction, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(Prediction = str_replace_all(Prediction, pattern = "Oligodendrocytes", replacement = "MOL"))

glue(
  "number of AST: {as.character(scbr_g4@meta.data %>% dplyr::filter(Prediction == 'AST') %>% rownames %>% length)}
     number of non-AST: {as.character(scbr_g4@meta.data %>% dplyr::filter(Prediction != 'AST') %>% rownames %>% length)}"
)

glue(
  "Av. reliability of AST: {as.character(scbr_g4@meta.data %>% dplyr::filter(Prediction == 'AST') %>% pull('Reliability') %>% 
  mean %>% round(2))}"
)

glue(
  "# of reliable cells: {as.character(scbr_g4@meta.data %>% dplyr::filter(Reliability > 0.9) %>% rownames %>% 
  length)} 
  # of all cells: {as.character(dim(scbr_g4@meta.data)[1])}"
)

# Seurat UMAPs
set3 = brewer.pal(8, "Set3")
cols = c('AST'=set3[1],'COP-NFOL'=set3[2],'MOL'=set3[3],'OPC'=set3[4],
         'OEC'=set3[5],'VEC'=set3[6],'VLMC'=set3[7], 'Pericytes'=set3[8])

# cell types
scbr_rna_pred = DimPlot(
  object = scbr_rna,
  group.by = "CellType",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  raster = TRUE,
  repel = TRUE,
  cols = cols,
  order = c('Pericytes', 'VLMC', 'VEC', 'COP-NFOL', 'OPC', 'OEC', 'MOL', 'AST')
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Cell type (scRNA-Seq)") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_rna_pred

scbr_g4@meta.data = scbr_g4@meta.data %>% 
  mutate(Prediction = as.character(Prediction)) %>% 
  mutate(Prediction = ifelse(str_detect(Prediction, "Novel"), "unreliable", Prediction)) %>% 
  mutate(Prediction = as.factor(Prediction))

cols = c('AST'=set3[1],'COP-NFOL'=set3[2],'MOL'=set3[3],'OPC'=set3[4],
         'OEC'=set3[5],'VEC'=set3[6],'VLMC'=set3[7], 'Pericytes'=set3[8], 'unreliable'='#f0f0f0')
         
# scBridge prediction score
scbr_g4_pred = DimPlot(
  object = scbr_g4,
  group.by = "Prediction",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = TRUE,
  cols = cols,
  order = c('Pericytes', 'VLMC', 'VEC', 'COP-NFOL', 'OPC', 'OEC', 'MOL', 'AST', 'unreliable')
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Prediction") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
scbr_g4_pred

# scBridge reliability
rel = FeaturePlot(object = scbr_g4, features = 'Reliability', raster = TRUE) +
  scale_color_viridis() +
  xlim(-12, 25) +
  ylim(-20, 20) +
  ggtitle("Reliability") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
rel

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

seurat_pred_score = FeaturePlot(object = scbr_g4, features = 'Seurat_pred_score', raster = TRUE) +
  scale_color_viridis() +
  xlim(-12, 25) +
  ylim(-20, 20) +
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

scbr_g4@meta.data = scbr_g4@meta.data %>% 
  mutate(Seurat_prediction = 
           str_replace_all(Seurat_prediction, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(Seurat_prediction = 
           str_replace_all(Seurat_prediction, pattern = "Oligodendrocytes", replacement = "MOL"))

cols = c('AST'=set3[1],'COP-NFOL'=set3[2],'MOL'=set3[3],'OPC'=set3[4],
         'OEC'=set3[5],'VEC'=set3[6],'VLMC'=set3[7], 'Pericytes'=set3[8])

seurat_pred =DimPlot(
  object = scbr_g4,
  group.by = "Seurat_prediction",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = TRUE,
  cols = cols,
  order = c('Pericytes', 'VLMC', 'VEC', 'COP-NFOL', 'OPC', 'OEC', 'MOL', 'AST')
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
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

scbr_g4@meta.data = scbr_g4@meta.data %>% mutate(AST_status = ifelse(Prediction == "AST", "AST", "non-AST"))
cols = c(
  "non-AST" = '#f0f0f0',
  "AST" = "#de2d26")

ast = DimPlot(
  object = scbr_g4,
  group.by = "AST_status",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = TRUE,
  order = c("AST", "non-AST")
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
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
ast

ggsave(
  glue("{result_folder}Seurat_predAST-UMAP.pdf"),
  plot = ast,
  width = 8,
  height = 8,
  device = "pdf"
)

# experimental system (modality)
comb@meta.data = comb@meta.data %>% 
  mutate(Domain = as.character(Domain)) %>% 
  mutate(Domain = ifelse(str_detect(Domain, "scCutTag_"), "G4 scCut&Tag", 
                             Domain)) %>% 
  mutate(Domain = ifelse(str_detect(Domain, "Bartosovic"), "scRNA-Seq", 
                             Domain)) %>% 
  mutate(Domain = as.character(Domain)) 

domain = DimPlot(
  object = comb,
  group.by = "Domain",
  label = FALSE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE,
  raster = TRUE,
  alpha = 0.1
) +
  xlim(-12, 25) +
  ylim(-20, 20) +
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
domain

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
genes = c("Pygb", "Pitpnc1", "Gpam", "Nwd1")
for(gene in genes) {
  plot = FeaturePlot(object = scbr_g4, features = gene, order = TRUE, reduction = "X_umap", pt.size = 2) +
    #scale_color_gradient2(low = "#f0f0f0", mid = "#f0f0f0", high = "red") +
    scale_color_viridis() +
    xlim(-12, 25) +
    ylim(-20, 20) +
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

ast_featureplotes = ggarrange(
  plotlist = featureplots,
  ncol = 2,
  nrow = 2
)

ggsave(
  glue("{result_folder}Seurat_UMAPs-AST_spec_G4_markers.pdf"),
  plot = ast_featureplotes,
  width = 18,
  height = 18,
  device = "pdf"
)

featureplots_rna = list()
genes = c("Pygb", "Pitpnc1", "Gpam", "Nwd1")
for(gene in genes) {
  plot = FeaturePlot(object = scbr_rna, features = gene, 
                     order = FALSE, raster = TRUE, reduction = "X_umap", pt.size = 2) +
    #scale_color_gradient2(low = "#f0f0f0", mid = "#f0f0f0", high = "red") +
    scale_color_viridis() +
    xlim(-12, 25) +
    ylim(-20, 20) +
    ggtitle(gene) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    )
  featureplots_rna[[gene]] = plot
}

ast_featureplotes_rna = ggarrange(
  plotlist = featureplots_rna,
  ncol = 2,
  nrow = 2
)

ggsave(
  glue("{result_folder}Seurat_UMAPs-AST_spec_G4_markers-RNA_level.pdf"),
  plot = ast_featureplotes_rna,
  width = 18,
  height = 18,
  device = "pdf"
)

# scBridge output tables
rel = fread("../results/scBridge/output/scbridge_reliability.csv")
pred = fread("../results/scBridge/output/scbridge_predictions.csv", header = TRUE)

meta = g4@meta.data %>% 
  left_join(., rel, c("cell_id" = "V1")) %>% dplyr::rename(scBridge_reliability = "Reliability") %>% 
  left_join(., pred, c("cell_id" = "V1")) %>% dplyr::rename(scBridge_prediction = "Prediction")

rownames(meta) = rownames(g4@meta.data) 
g4@meta.data = meta

both_preds_ast = meta %>% filter(scBridge_prediction == "AST" & pred_cell_type == "AST")

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



