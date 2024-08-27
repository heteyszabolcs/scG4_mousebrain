if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse",
               "Seurat",
               "Signac",
               "anndata",
               "data.table",
               "ggplot2",
               "ggpubr",
               "glue"
)

result_folder = "../results/Seurat/"

# Seurat integrated object
g4 = readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/G4_scRNA_integration.Rds")
g4@meta.data = g4@meta.data %>% 
  mutate(pred_cell_type = str_replace_all(pred_cell_type, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(pred_cell_type = str_replace_all(pred_cell_type, pattern = "Oligodendrocytes", replacement = "MOL"))

order = g4@meta.data %>% group_by(pred_cell_type) %>% summarize(med = median(pred_max_score)) %>% 
  arrange(desc(med)) %>% pull(pred_cell_type)
order = factor(g4@meta.data$pred_cell_type, levels = order)


stats = compare_means(pred_max_score ~ pred_cell_type,  data = g4@meta.data)

ggplot(g4@meta.data, aes(x = order, y = pred_max_score)) +
  geom_boxplot(color = "black", fill = "#9ecae1") +
  ylim(0, 1) +
  labs(
    title = "Max Seurat prediction score distributions",
    x = "predicted label",
    y = "Seurat prediction score"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 13, color = "black")
  ) 

ggsave(
  glue("{result_folder}Seurat-pred_score_distr.pdf"),
  plot = last_plot(),
  width = 6,
  height = 4,
  device = "pdf"
)
