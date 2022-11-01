# packages
suppressPackageStartupMessages({
  library("glue")
  library("Seurat")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("ggplot2")
})

# result / peak folder
peak_folder = "../results/Seurat/callpeaks_GFPsorted/peak_sets/"
result_folder = "../results/Seurat/"
cellranger_folder = "../data/GSE163484/"

# Signac MACS2 peak set
signac_macs2 = "enhancer_analysis_output.tsv"

# Lanceotron peak set (annotated)
lanceotron_annot = "lanceotron_0.10filt_peak_set_annot.tsv"
lanceotron_annot = fread(glue("{peak_folder}{lanceotron_annot}"))
lanceotron_annot = lanceotron_annot %>%
  separate(
    `PeakID (cmd=annotatePeaks.pl ../data/bed/lanceotron_0.10filt_peak_set.bed mm10)`,
    "NA-",
    into = c("temp", "ID")
  ) %>% select(-temp) %>% mutate(as.character(ID))
lanceotron = "lanceotron_0.10filt_peak_set.tsv"
lanceotron = fread(glue("{peak_folder}{lanceotron}"))
lanceotron = lanceotron %>% rowid_to_column(., "ID") %>% mutate(ID = as.character(ID)) %>%
  inner_join(., lanceotron_annot, by = c("ID" = "ID")) %>%
  select(
    Chr = Chr.x,
    Start = Start.x,
    End = End.x,
    Peak_score = `Peak Score.x`,
    Seurat_cluster,
    Annotation,
    `Distance to TSS`,
    `Gene Name`
  )
rm(lanceotron_annot)

### replication 1 - GSM4979874 ### 
# read and process data
counts1 = Read10X_h5(
  filename = glue(
    "{cellranger_folder}scRNASeq_GSM4979874_filtered_feature_bc_matrix_rep1.h5"
  )
)
rna1 = CreateSeuratObject(
  counts = counts1,
  project = "mouse_brain",
  min.cells = 3,
  min.features = 200
)

all.genes = rownames(rna1)
rna1 = ScaleData(rna1, features = all.genes)
rna1[["replicate"]] <- "replicate 1"
# Perform linear dimensional reduction and clustering
rna1 = FindVariableFeatures(rna1, selection.method = "vst", nfeatures = 2000)
rna1 = RunPCA(rna1, features = VariableFeatures(object = rna1))
rna1 = FindNeighbors(rna1, dims = 1:10)
rna1 = FindClusters(rna1, resolution = 0.5)

# Run non-linear dimensional reduction
rna1 = RunUMAP(
  rna1,
  dims = 1:10,
  n.neighbors = 40,
  n.components = 2
)
umap1 = DimPlot(rna1, reduction = "umap")
umap1

#  global-scaling normalization by LogNormalize:
# Feature counts for each cell are divided by the total counts for that cell
# and multiplied by the scale.factor. This is then natural-log transformed using log1p
rna1 = NormalizeData(rna1,
                     normalization.method = "LogNormalize")

norm1 = rna1[["RNA"]]@data

### replication 2 - GSM4979874 ###
# read and process data
counts2 = Read10X_h5(
  filename = glue(
    "{cellranger_folder}scRNASeq_GSM4979875_filtered_feature_bc_matrix_rep2.h5"
  )
)
rna2 = CreateSeuratObject(
  counts = counts2,
  project = "mouse_brain",
  min.cells = 3,
  min.features = 200
)

all.genes = rownames(rna2)
rna2 = ScaleData(rna2, features = all.genes)
rna2[["replicate"]] <- "replicate 2"
# Perform linear dimensional reduction and clustering
rna2 = FindVariableFeatures(rna2, selection.method = "vst", nfeatures = 2000)
rna2 = RunPCA(rna2, features = VariableFeatures(object = rna2))
rna2 = FindNeighbors(rna2, dims = 1:10)
rna2 = FindClusters(rna2, resolution = 0.5)
# Run non-linear dimensional reduction
rna2 = RunUMAP(
  rna2,
  dims = 1:10,
  n.neighbors = 40,
  n.components = 2
)
umap2 = DimPlot(rna2, reduction = "umap")
umap2

rna = merge(rna1, y = rna2, add.cell.ids = c("replicate1", "replicate2"), project = "scRNA_brain")

#  global-scaling normalization by LogNormalize:
# Feature counts for each cell are divided by the total counts for that cell
# and multiplied by the scale.factor. This is then natural-log transformed using log1p
rna = NormalizeData(rna,
                     normalization.method = "LogNormalize")

saveRDS(rna, glue("{result_folder}scRNASeq_GSM4979874-75.rds"))

rm(rna1); rm(rna2)
rm(norm1); rm(norm2)

# filtered MACS2 peak calls
signac_macs2 = fread(glue("{peak_folder}{signac_macs2}"))

create_heatmap = function(genes = cluster4,
                          data = rna,
                          scrna_cluster_n = 15,
                          title) {
  # fetch normalized expression
  norm = data[["RNA"]]@data
  
  existing_gene_symbols = character()
  for (gene in genes) {
    if (gene %in% norm@Dimnames[[1]]) {
      existing_gene_symbols = c(existing_gene_symbols, gene)
    }
  }
  
  existing_gene_symbols = existing_gene_symbols[1:10]
  
  ms = list()
  for (cluster in seq(0, scrna_cluster_n)) {
    m = t(as.matrix(norm[existing_gene_symbols,]))
    cluster_barcodes = WhichCells(data, idents = cluster)
    m = m[cluster_barcodes,]
    t = tibble(means = colMeans(m))
    t = t %>%
      mutate(gene_symbol = existing_gene_symbols) %>%
      mutate(Seurat_cluster = as.character(cluster))
    ms[[as.character(cluster)]] = t
  }
  ms = bind_rows(ms)
  
  y_order = ms %>% group_by(gene_symbol) %>% summarise(mean = mean(means)) %>% arrange(mean) %>%
    pull(gene_symbol)
  y_order = factor(ms$gene_symbol, levels = y_order)
  
  x_order = ms %>% group_by(Seurat_cluster) %>% summarise(mean = mean(means)) %>% arrange(desc(mean)) %>%
    pull(Seurat_cluster)
  x_order = factor(ms$Seurat_cluster, levels = x_order)
  
  hm = ggplot(ms, aes(x = x_order, y = y_order, fill = means)) +
    geom_tile(color = "white",
              lwd = 1.0,
              linetype = 1) +
    scale_fill_gradient2(
      low = "#075AFF",
      mid = "#FFFFCC",
      high = "#FF0000",
      limits = c(0, 3)
    ) +
    xlab(label = "scRNA-Seq Seurat cluster") +
    ylab(label = "nearest gene to top G4 peak") +
    labs(fill = "log norm", title = title) +
    theme(
      axis.text.x = element_text(
        color = "black",
        size = 15,
        angle = 0,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title = element_text(size = 14)
    ) +
    coord_fixed()
  print(hm)
  
  return(hm)
  
}

create_lanc_heatmap = function(genes = lanc_top_genes_igenic,
                               data = rna2,
                               scrna_cluster_n = 15,
                               title) {
  # fetch normalized expression
  norm = data[["RNA"]]@data
  
  existing_gene_symbols = character()
  for (gene in genes) {
    if (gene %in% norm@Dimnames[[1]]) {
      existing_gene_symbols = c(existing_gene_symbols, gene)
    }
  }
  
  existing_gene_symbols = existing_gene_symbols[1:20]
  
  ms = list()
  for (cluster in seq(0, scrna_cluster_n)) {
    m = t(as.matrix(norm[existing_gene_symbols,]))
    cluster_barcodes = WhichCells(data, idents = cluster)
    m = m[cluster_barcodes,]
    t = tibble(means = colMeans(m))
    t = t %>%
      mutate(gene_symbol = existing_gene_symbols) %>%
      mutate(Seurat_cluster = as.character(cluster))
    ms[[as.character(cluster)]] = t
  }
  ms = bind_rows(ms)
  
  y_order = ms %>% group_by(gene_symbol) %>% summarise(mean = mean(means)) %>% arrange(mean) %>%
    pull(gene_symbol)
  y_order = factor(ms$gene_symbol, levels = y_order)
  
  x_order = ms %>% group_by(Seurat_cluster) %>% summarise(mean = mean(means)) %>% arrange(desc(mean)) %>%
    pull(Seurat_cluster)
  x_order = factor(ms$Seurat_cluster, levels = x_order)
  
  hm = ggplot(ms, aes(x = x_order, y = y_order, fill = means)) +
    geom_tile(color = "white",
              lwd = 1.0,
              linetype = 1) +
    scale_fill_gradient2(
      low = "#075AFF",
      mid = "#FFFFCC",
      high = "#FF0000",
      limits = c(0, 3)
    ) +
    xlab(label = "scRNA-Seq Seurat cluster") +
    ylab(label = "nearest gene to top G4 peak") +
    labs(fill = "log norm", title = title) +
    theme(
      axis.text.x = element_text(
        color = "black",
        size = 15,
        angle = 0,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title = element_text(size = 14)
    ) +
    coord_fixed()
  print(hm)
  
  return(hm)
  
}


# filtered MACS2 peak calls
get_top_genes = function(peak_set = signac_macs2, cluster) {
  cluster = peak_set %>% filter(Seurat_cluster == cluster) %>%
    filter(Distance_to_TSS > -3000) %>%
    filter(Distance_to_TSS < 3000) %>%
    arrange(desc(signalValue)) %>%
    top_n(15, signalValue) %>%
    pull(Gene_name)
  
  return(cluster)
  
}

cluster0 = get_top_genes(cluster = 0)
cluster0 = create_heatmap(genes = cluster0, title = "Seurat cluster 0")
ggsave(
  glue("{result_folder}heatmap_cluster0_signacMACS2_topgenes.png"),
  plot = cluster0,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster1 = get_top_genes(cluster = 1)
cluster1 = create_heatmap(genes = cluster1, title = "Seurat cluster 1")
ggsave(
  glue("{result_folder}heatmap_cluster1_signacMACS2_topgenes.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster2 = get_top_genes(cluster = 2)
cluster2 = create_heatmap(genes = cluster2, title = "Seurat cluster 2")
ggsave(
  glue("{result_folder}heatmap_cluster2_signacMACS2_topgenes.png"),
  plot = cluster2,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster3 = get_top_genes(cluster = 3)
cluster3 = create_heatmap(genes = cluster3, title = "Seurat cluster 3")
ggsave(
  glue("{result_folder}heatmap_cluster3_signacMACS2_topgenes.png"),
  plot = cluster3,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster4 = get_top_genes(cluster = 4)
cluster4 = create_heatmap(genes = cluster4, title = "Seurat cluster 4")
ggsave(
  glue("{result_folder}heatmap_cluster4_signacMACS2_topgenes.png"),
  plot = cluster4,
  width = 10,
  height = 10,
  dpi = 300,
)


# filtered Lanceotron peak calls
get_top_lanc_genes = function(peak_set = lanceotron, cluster) {
  cluster = peak_set %>% filter(Seurat_cluster == cluster) %>%
    filter(str_detect(Annotation, "TSS")) %>%
    arrange(desc(Peak_score)) %>%
    top_n(15, Peak_score) %>%
    pull("Gene Name")
  return(cluster)
}

cluster0 = get_top_lanc_genes(cluster = 0)
cluster0 = create_heatmap(genes = cluster0, title = "Seurat cluster 0")
ggsave(
  glue("{result_folder}heatmap_cluster0_LanceOtron_topgenes.png"),
  plot = cluster0,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster1 = get_top_lanc_genes(cluster = 1)
cluster1 = create_heatmap(genes = cluster1, title = "Seurat cluster 1")
ggsave(
  glue("{result_folder}heatmap_cluster1_LanceOtron_topgenes.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster2 = get_top_lanc_genes(cluster = 2)
cluster2 = create_heatmap(genes = cluster2, title = "Seurat cluster 2")
ggsave(
  glue("{result_folder}heatmap_cluster2_LanceOtron_topgenes.png"),
  plot = cluster2,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster3 = get_top_lanc_genes(cluster = 3)
cluster3 = create_heatmap(genes = cluster3, title = "Seurat cluster 3")
ggsave(
  glue("{result_folder}heatmap_cluster3_LanceOtron_topgenes.png"),
  plot = cluster3,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster4 = get_top_lanc_genes(cluster = 4)
cluster4 = create_heatmap(genes = cluster4, title = "Seurat cluster 4")
ggsave(
  glue("{result_folder}heatmap_cluster4_LanceOtron_topgenes.png"),
  plot = cluster4,
  width = 10,
  height = 10,
  dpi = 300,
)

# lanc_top_genes = get_top_lanc_genes()
# lanc_hm = create_lanc_heatmap(genes = lanc_top_genes, title = "Lanceotron top TSS proximal peaks")
# ggsave(
#   glue("{result_folder}heatmap_Lanceotron_topgenes.png"),
#   plot = lanc_hm,
#   width = 10,
#   height = 10,
#   dpi = 300,
# )
