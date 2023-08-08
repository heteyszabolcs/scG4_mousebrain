# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("tidyverse")
  library("data.table")
  library("glue")
  library("devtools")
  library("pracma")
})

set.seed(5)

# read into Seurat
seurat = readRDS("../results/Seurat/final/unsorted_brain/res0.8/outputs/Seurat_object.Rds")
umap_folder = "../results/Seurat/callpeaks_unsorted/UMAPs/"

# testing of number of n of neighbours and min distance parameters
umap_test = function(n, dist) {
  seurat = RunUMAP(
    object = seurat,
    reduction = 'lsi',
    dims = 2:30,
    n.neighbors = n,
    min.dist = dist
  )
  seurat = FindNeighbors(object = seurat,
                         reduction = 'lsi',
                         dims = 2:30)
  seurat = FindClusters(
    object = seurat,
    verbose = FALSE,
    resolution = 0.1,
    algorithm = 3
  )
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
    ggtitle(glue("n.neigh = {as.character(n)}, min.dist = {as.character(dist)}, res = 0.1")) +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    )
  
  ggsave(
    glue("{umap_folder}neigh{as.character(n)}_mindist{as.character(dist)}_res0.1.pdf"),
    plot = dim,
    device = "pdf",
    width = 8,
    height = 8,
    dpi = 300,
  )
  
  return(print(dim))
}

ns = seq(5, 50, 2)
dists = round(linspace(0.001, 0.5, n = 25), 4)
lapply(ns, function(x) sapply(dists, function(y) umap_test(n = x, dist = y)))

# testing of clustering algorithms
# 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm, Leiden requires the leidenalg python.
# umap_algo_test = function(algo) {
#   seurat = RunUMAP(
#     object = seurat,
#     reduction = 'lsi',
#     dims = 2:30
#   )
#   seurat = FindNeighbors(object = seurat,
#                          reduction = 'lsi',
#                          dims = 2:30)
#   seurat = FindClusters(
#     object = seurat,
#     verbose = FALSE,
#     resolution = 0.8,
#     algorithm = algo
#   )
#   dim = DimPlot(
#     object = seurat,
#     label = TRUE,
#     pt.size = 2,
#     label.size = 7,
#     repel = TRUE
#   ) +
#     NoLegend() +
#     scale_color_brewer(palette = "Set3") +
#     xlim(-10, 10) +
#     ylim(-10, 10) +
#     ggtitle(glue("clust algo = {as.character(algo)}, res = 0.1")) +
#     theme(
#       text = element_text(size = 25),
#       plot.title = element_text(size = 20),
#       axis.text.x = element_text(size = 25, color = "black"),
#       axis.text.y = element_text(size = 25, color = "black")
#     )
#   
#   if(algo == 1) {
#     name_of_algo = "original_Louvain"
#   } else if(algo == 2) {
#     name_of_algo = "multilevel_Louvain"
#   } else if(algo == 3) {
#     name_of_algo = "SLM"
#   } else if(algo == 4) {
#     name_of_algo = "Leiden"
#   }
#   
#   ggsave(
#     glue("{umap_folder}algo_{name_of_algo}_res0.8.pdf"),
#     plot = dim,
#     device = "pdf",
#     width = 8,
#     height = 8,
#     dpi = 300,
#   )
#   
#   return(print(dim))
# }
# 
# algos = c(1, 2, 3)
# lapply(algos, umap_algo_test)





