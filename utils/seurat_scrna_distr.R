# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("GenomicFeatures")
  library("Matrix")
  library("cowplot")
  library("ggpubr")
})

# scRNA-Seq folder (Marek's postnatal mouse brain scRNA-Seq data)
cellranger_folder = "../data/GSE163484/"
# path to result folder
result_folder = "../results/Seurat/"

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

rm(rna1); rm(rna2)

# function for visualization of expression of genes near to top peaks 
expr_distr = function(genes,
                      data = rna,
                      title, 
                      scrna_cluster_n = 15) {
  
  # fetch normalized expression
  norm = data[["RNA"]]@data
  
  # exclude missing genes
  existing_gene_symbols = character()
  for(gene in genes) {
    if(gene %in% norm@Dimnames[[1]]) {
      existing_gene_symbols = c(existing_gene_symbols, gene)
    }
  }
  
  ms = list()
  for(cluster in seq(0, scrna_cluster_n)) {
    m = t(as.matrix(norm[existing_gene_symbols, ]))
    cluster_barcodes = WhichCells(data, idents = cluster)
    m = m[cluster_barcodes, ]
    cluster_barcodes = rownames(m)
    m = as_tibble(m)
    m = pivot_longer(m,
                     existing_gene_symbols,
                     names_to = "gene_symbol",
                     values_to = "log_norm_expr") 
    m = m %>% mutate(cluster = as.character(cluster))
    ms[[as.character(cluster)]] = m
  }
  ms = bind_rows(ms)
  
  ms = ms %>% arrange(cluster)
  
  # boxplot
  p = ggplot(ms, aes(
    x = reorder(cluster,-log_norm_expr),
    y = log_norm_expr
  )) +
    geom_boxplot(fill = "#2ca25f") +
    theme_minimal() +
    ylim(0, 10) +
    labs(title = title,
         fill = "") +
    xlab(label = "scRNA-Seq Seurat cluster") +
    ylab(label = "log2 normalized expr.") +
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
    ) 
  print(p)
  
  return(p)
  
}

expr_distr_genelevel = function(genes = c("Tpmt", "Ccne1"),
                                data = rna,
                                title = "ArchR cluster 6") {
  
  # fetch normalized expression
  norm = data[["RNA"]]@data
  
  # exclude missing genes
  existing_gene_symbols = character()
  for(gene in genes) {
    if(gene %in% norm@Dimnames[[1]]) {
      existing_gene_symbols = c(existing_gene_symbols, gene)
    }
  }
  
  
  m = t(as.matrix(norm[existing_gene_symbols,]))
  cluster_barcodes = rownames(m)
  m = as_tibble(m)
  m = pivot_longer(m,
                   existing_gene_symbols,
                   names_to = "gene_symbol",
                   values_to = "log_norm_expr")
  
  m = m %>% arrange(gene_symbol)
  
  # boxplot
  p = ggplot(m, aes(
    x = reorder(gene_symbol,-log_norm_expr),
    y = log_norm_expr
  )) +
    geom_boxplot(fill = "#2ca25f") +
    theme_minimal() +
    ylim(0, 10) +
    labs(title = title,
         fill = "") +
    xlab(label = "nearest gene to top G4 peak") +
    ylab(label = "log2 normalized expr.") +
    theme(
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
  print(p)
  
  return(p)
  
}

# get gene lists for top Signac MACS2 peaks
peak_folder = "../results/Seurat/callpeaks_GFPsorted/peak_sets/"

get_top_genes = function(peak_set = glue("{peak_folder}enhancer_analysis_output.tsv"), cluster) {
  peak_set = read_tsv(peak_set, show_col_types = FALSE)
  cluster = peak_set %>% filter(Seurat_cluster == cluster) %>% 
    filter(Distance_to_TSS > -3000) %>% 
    filter(Distance_to_TSS < 3000) %>% 
    arrange(desc(signalValue)) %>% 
    top_n(10, signalValue) %>% 
    pull(Gene_name)
  
  return(cluster)

}

# make boxplots for each Seurat clusters
cluster0 = expr_distr(genes = get_top_genes(cluster = 0), title = "Seurat cluster 0", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster0_signacMACS2_topgenes.png"),
  plot = cluster0,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster1 = expr_distr(genes = get_top_genes(cluster = 1), title = "Seurat cluster 1", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster1_signacMACS2_topgenes.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster2 = expr_distr(genes = get_top_genes(cluster = 2), title = "Seurat cluster 2", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster2_signacMACS2_topgenes.png"),
  plot = cluster2,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster3 = expr_distr(genes = get_top_genes(cluster = 3), title = "Seurat cluster 3", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster3_signacMACS2_topgenes.png"),
  plot = cluster3,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster4 = expr_distr(genes = get_top_genes(cluster = 4), title = "Seurat cluster 4", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster4_signacMACS2_topgenes.png"),
  plot = cluster4,
  width = 10,
  height = 10,
  dpi = 300,
)

# Perform linear dimensional reduction and clustering
rna = ScaleData(rna, features = all.genes)
rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna = RunPCA(rna, features = VariableFeatures(object = rna))
rna = FindNeighbors(rna, dims = 1:10)
rna = FindClusters(rna, resolution = 0.5)
# Run non-linear dimensional reduction
rna = RunUMAP(
  rna,
  dims = 1:10,
  n.neighbors = 40,
  n.components = 2
)
umap = DimPlot(rna, reduction = "umap")
umap

ggsave(
  glue("{result_folder}UMAP_GSM4979874_scRNASeq.png"),
  plot = umap2,
  width = 10,
  height = 10,
  dpi = 300,
)

# at G4 related gene level
cluster0 = expr_distr_genelevel(genes = get_top_genes(cluster = 0), title = "Seurat cluster 0")
ggsave(
  glue("{result_folder}boxplot_cluster0_signacMACS2_topgenes_genelevel.png"),
  plot = cluster0,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster1 = expr_distr_genelevel(genes = get_top_genes(cluster = 1), title = "Seurat cluster 1")
ggsave(
  glue("{result_folder}boxplot_cluster1_signacMACS2_topgenes_genelevel.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster2 = expr_distr_genelevel(genes = get_top_genes(cluster = 2), title = "Seurat cluster 2")
ggsave(
  glue("{result_folder}boxplot_cluster2_signacMACS2_topgenes_genelevel.png"),
  plot = cluster2,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster3 = expr_distr_genelevel(genes = get_top_genes(cluster = 3), title = "Seurat cluster 3")
ggsave(
  glue("{result_folder}boxplot_cluster3_signacMACS2_topgenes_genelevel.png"),
  plot = cluster3,
  width = 10,
  height = 10,
  dpi = 300,
)

cluster4 = expr_distr_genelevel(genes = get_top_genes(cluster = 4), title = "Seurat cluster 4")
ggsave(
  glue("{result_folder}boxplot_cluster4_signacMACS2_topgenes_genelevel.png"),
  plot = cluster4,
  width = 10,
  height = 10,
  dpi = 300,
)


# get gene lists for top ArchR MACS2 peaks 
get_archr_top_genes = function(peak_set = glue("{peak_folder}archr_0.50perc_peak_set.tsv"), 
                               cluster) {
  peak_set = read_tsv(peak_set, show_col_types = FALSE)
  cluster = peak_set %>% filter(peakType == "Promoter") %>% 
    filter(distToGeneStart > -3000) %>% 
    filter(distToGeneStart < 3000) %>% 
    arrange(desc(score)) %>% 
    top_n(10, score) %>% 
    pull(nearestGene)
  
  return(cluster)
  
}

# make boxplots for each ArchR GFP sorted clusters
archr_peaks = read_tsv(glue("{peak_folder}archr_0.50perc_peak_set.tsv", show_col_types = FALSE))
clusters = archr_peaks %>% pull(GroupReplicate) %>% unique
clusters = clusters[str_detect(clusters, "GFP")]

archr_cluster3 = expr_distr(genes = get_archr_top_genes(cluster = "C3._.GFP"), 
                            title = "ArchR cluster 3", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster3_ArchRMACS2_topgenes.png"),
  plot = cluster0,
  width = 10,
  height = 10,
  dpi = 300,
)

archr_cluster4 = expr_distr(genes = get_archr_top_genes(cluster = "C4._.GFP"), 
                            title = "ArchR cluster 4", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster4_ArchRMACS2_topgenes.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)

archr_cluster6 = expr_distr(genes = get_archr_top_genes(cluster = "C6._.GFP"), 
                            title = "ArchR cluster 6", scrna_cluster_n = 15)
ggsave(
  glue("{result_folder}boxplot_cluster6_ArchRMACS2_topgenes.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)

archr_cluster3 = expr_distr_genelevel(genes = get_archr_top_genes(cluster = "C3._.GFP"), 
                            title = "ArchR cluster 3")
ggsave(
  glue("{result_folder}boxplot_cluster3_ArchRMACS2_topgenes_genelevel.png"),
  plot = cluster0,
  width = 10,
  height = 10,
  dpi = 300,
)

archr_cluster4 = expr_distr_genelevel(genes = get_archr_top_genes(cluster = "C4._.GFP"), 
                            title = "ArchR cluster 4")
ggsave(
  glue("{result_folder}boxplot_cluster4_ArchRMACS2_topgenes_genelevel.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)

archr_cluster6 = expr_distr_genelevel(genes = get_archr_top_genes(cluster = "C6._.GFP"), 
                            title = "ArchR cluster 6")
ggsave(
  glue("{result_folder}boxplot_cluster6_ArchRMACS2_topgenes_genelevel.png"),
  plot = cluster1,
  width = 10,
  height = 10,
  dpi = 300,
)
