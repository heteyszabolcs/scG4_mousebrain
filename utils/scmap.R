library("SingleCellExperiment")
library("scmap")
library("glue")

result_folder = "../results/Seurat/callpeaks_GFPsorted/"
seurat_folder = "../results/Seurat/"

sorted = readRDS(glue("{result_folder}sorted_int_Marques.Rds"))
rna = readRDS(glue("{seurat_folder}scRNASeq_GSE75330.rds"))

# prepare scRNA-Seq data
rna_se = as.SingleCellExperiment(rna)
rowData(rna_se)$feature_symbol = rownames(rna_se)
rna_se = rna_se[!duplicated(rownames(rna_se)), ]
rna_se = selectFeatures(rna_se, suppress_plot = FALSE)
rna_se = indexCluster(rna_se, cluster_col = "merged_cell_class")

heatmap(as.matrix(metadata(rna_se)$scmap_cluster_index))

# prepare G4 data
g4 = as.SingleCellExperiment(sorted)
rowData(g4)$feature_symbol = rownames(g4)
g4 = g4[!duplicated(rownames(g4)), ]
g4 = selectFeatures(g4, suppress_plot = FALSE)
g4 = indexCluster(g4, cluster_col = "seurat_clusters")

# project sc G4 data to scRNA-Seq
scmapCluster_results = scmapCluster(
  projection = g4, 
  index_list = list(
    yan = metadata(rna_se)$scmap_cluster_index
  ), threshold = 0.1
)

# Sankey plot
plot(
  getSankey(
    colData(g4)$seurat_clusters, 
    scmapCluster_results$scmap_cluster_labs[,'yan'],
    plot_height = 400, colors = c("#8dd3c7", "#ffffb3", "#bebada", "fb8072", "#80b1d3")
  )
)

similarities = tibble(similarity = scmapCluster_results$scmap_cluster_siml)



