# packages
suppressPackageStartupMessages({
  library("cluster")
  library("Seurat")
  library("Signac")
  library("ggplot2")
  library("dplyr")
  library("glue")
  library("tidyverse")
  library("data.table")
  library("ggpubr")
})

result_folder = "../results/Seurat/"

# silhouette score (cluster package)
# a measure of how similar a cell is to its own cluster compared to other clusters
# source: Tim Stuart git repo 
# https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36

# function
make_silhouette_bp = function(seurat) {
  # runPCA
  seurat = ScaleData(seurat, verbose = FALSE)
  seurat = RunPCA(seurat, verbose = FALSE)
  dist.matrix = dist(Embeddings(object = seurat[["pca"]])[, 1:30])
  clusters = seurat$seurat_clusters
  
  # calculate scores
  sil = silhouette(x = as.numeric(x = clusters), dist = dist.matrix)
  sil = as_tibble(sil)
  
  # reset Seurat cluster names
  sil_clusters = as.character(sort(unique(sil$cluster)))
  orig_clusters = sort(levels(clusters))
  for (i in seq(1, length(orig_clusters))) {
    sil = sil %>% mutate(cluster = as.character(cluster)) %>%
      mutate(cluster = ifelse(cluster == sil_clusters[i], orig_clusters[i], cluster))
  }
  
  # boxplot
  bp = ggplot(sil, aes(x = cluster, y = sil_width, fill = cluster)) +
    geom_boxplot(color = "black") +
    scale_fill_brewer(palette = "Pastel1") +
    ylim(-1, 1.2) +
    labs(
      title = "Silhouette scores",
      x = "",
      y = "silhouette coefficient",
      fill = ""
    ) +
    theme_classic() +
    guides(fill = "none") +
    theme(
      text = element_text(size = 9),
      plot.title = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black")
    ) + stat_compare_means(label = "p.signif", label.y = 1)
  bp
  
}

# run Seurat objects
unsorted_res0.8 = readRDS("../results/Seurat/final/unsorted_brain/res0.8/outputs/Seurat_object.Rds")
unsorted_res0.8_bp = make_silhouette_bp(seurat = unsorted_res0.8)
sorted_res0.8 = readRDS("../results/Seurat/final/sorted_brain/res0.8/outputs/Seurat_object.Rds")
sorted_res0.8_bp = make_silhouette_bp(seurat = sorted_res0.8)

ggsave(
  glue("{result_folder}sorted_res0.8_silhouette_bp.pdf"),
  plot = sorted_res0.8_bp,
  width = 6,
  height = 6,
  device = "pdf"
)
