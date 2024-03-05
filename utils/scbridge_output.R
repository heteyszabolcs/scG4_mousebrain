if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("Seurat",
               "glue",
               "tidyverse",
               "data.table",
               "ggplot2",
               "ggpubr"
)

# result folder
result_folder = "../results/scBridge/"

# integrated object
g4 = readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/G4_Marques_scRNA_integration.Rds")

# outputs
rel = fread("../results/scBridge/output/scbridge_reliability.csv")
pred = fread("../results/scBridge/output/scbridge_predictions.csv", header = TRUE)

meta = g4@meta.data %>% mutate(cell_id = rownames(g4@meta.data)) %>% 
  left_join(., rel, c("cell_id" = "V1")) %>% rename(scBridge_reliability = "Reliability") %>% 
  left_join(., pred, c("cell_id" = "V1")) %>% rename(scBridge_prediction = "Prediction")

rownames(meta) = rownames(g4@meta.data) 
g4@meta.data = meta

cor(g4@meta.data$pred_max_score, g4@meta.data$scBridge_reliability, method = "spearman")
ggplot(g4@meta.data, aes(x = scBridge_reliability, y = pred_max_score, color = scBridge_prediction)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Pastel1") +
  labs(
    title = " ",
    x = " ",
    y = " ",
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
