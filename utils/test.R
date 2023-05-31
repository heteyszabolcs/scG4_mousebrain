library("tidyverse")
library("ggplot2")
library("ggpubr")
library("Seurat")

x = tibble(type = c(rep("A", 10), c(rep("B", 10))), value = runif(20, 1, 99))
           

plot = ggplot(x, aes(x = type, y = value, fill = type)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Blues") +
  ylim(0, 150) +
  labs(
    title = "",
    x = "something",
    y = "log something expr.",
    fill = "type"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  ) +
  stat_compare_means(label.y = 135, label.x = 1.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 100)
plot

ggsave(
  "../results/test.png",
  plot = plot,
  width = 10,
  height = 10,
  dpi = 300,
)

seurat = readRDS("../results/Seurat/callpeaks_unsorted/unsorted.Rds")

dim = DimPlot(
  object = seurat,
  label = TRUE,
  pt.size = 2,
  label.size = 7,
  repel = TRUE
) +
  NoLegend() +
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) +
  ylim(-10, 10) +
  ggtitle("") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

ggsave(
  "../results/umap_test.png",
  plot = dim,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  "../results/umap_test.pdf",
  plot = dim,
  width = 10,
  height = 10
)