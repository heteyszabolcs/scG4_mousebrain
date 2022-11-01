# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("ggpubr")
  library("cowplot")
  library("Seurat")
  library("wigglescout")
  library("GenomicRanges")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
})

# result / peak folder
result_folder = "../results/Seurat/callpeaks_unsorted/"
peaks = "../results/Seurat/callpeaks_unsorted/Grubbs_test-unique_G4_peaks_0.01_joined.tsv"
peak_folder = "../results/Seurat/callpeaks_unsorted/peak_sets/"
marques_scrna = "../results/Seurat/scRNASeq_GSE75330.rds"
marques_scrna_lognorm = "../results/Seurat/scRNASeq_GSE75330_lognorm.tsv"

unique_peaks = fread(peaks)
prom_peaks = unique_peaks %>% filter(abs(`Distance to TSS`) < 1500)
genes = prom_peaks %>% pull(`Gene Name`)

marques_scrna =readRDS(marques_scrna)
marques_scrna_lognorm = fread("../results/Seurat/scRNASeq_GSE75330_lognorm.tsv")
marques_scrna_lognorm = marques_scrna_lognorm %>% rename(gene_name = V1)

clusters = tibble("cluster" = marques_scrna@active.ident, "cell_id" = names(marques_scrna@active.ident))

  # visualizations
prom_uniques = prom_peaks %>% group_by(unique) %>% summarise(count = n()) %>% arrange(desc(count)) %>% 
  ggplot(data = ., aes(
    x = reorder(unique ,-count),
    y = count,
    fill = unique
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "Reds") +
  labs(
    title = "# of unique G4 locations (-1.5/+1.5 kb TSS)",
    x = "Seurat cluster",
    y = "# of G4 structures",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
prom_uniques

# ggsave(
#   glue("{result_folder}Grubbs_test-unique_G4_peaks_bar.png"),
#   plot = n_uniques,
#   width = 10,
#   height = 10,
#   dpi = 300,
# )

get_expression = function(cell_type, g4_cluster) {
  
  genes = prom_peaks %>% filter(unique == g4_cluster) %>% pull(`Gene Name`)
  
  cell_ids = clusters %>% filter(cluster == cell_type) %>% pull(cluster) %>% names
  
  expr = marques_scrna_lognorm %>% 
    filter(gene_name %in% genes) %>% 
    dplyr::select(cell_ids) %>% 
    pivot_longer(starts_with("C"), names_to = "cell_id", values_to = "expression") %>% 
    mutate(cell_type = cell_type, "G4_cluster" = g4_cluster) %>% 
    dplyr::select(G4_cluster, cell_type, expression)
  
  return(expr)
  
}

# gene = marques_scrna_lognorm %>% filter(gene_name == "Gusb") %>% t

cell_types = levels(unique(unname(marques_scrna@active.ident)))

## cluster 0
cluster_0 = lapply(cell_types, get_expression, g4_cluster = "cluster_0")
cluster_0 = do.call(rbind, cluster_0)

cluster_0_bp = ggplot(cluster_0, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Reds") +
  ylim(0, 100) +
  labs(
    title = "Unique promoter peaks of G4 cluster 0 (-1.5/+1.5 kb TSS)",
    x = "scRNA-Seq (Marques et al.) Seurat cluster",
    y = "log normalized expr.",
    fill = "cell type"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  ) +
  stat_compare_means(label.y = 75, label.x = 3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 50) 
cluster_0_bp

## cluster 3
cluster_3 = lapply(cell_types, get_expression, g4_cluster = "cluster_3")
cluster_3 = do.call(rbind, cluster_3)

cluster_3_bp = ggplot(cluster_3, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Reds") +
  ylim(0, 100) +
  labs(
    title = "Unique promoter peaks of G4 cluster 3 (-1.5/+1.5 kb TSS)",
    x = "scRNA-Seq (Marques et al.) Seurat cluster",
    y = "log normalized expr.",
    fill = "cell type"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  ) +
  stat_compare_means(label.y = 75, label.x = 3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 50) 
cluster_3_bp

## cluster 4
cluster_4 = lapply(cell_types, get_expression, g4_cluster = "cluster_4")
cluster_4 = do.call(rbind, cluster_4)

cluster_4_bp = ggplot(cluster_4, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Reds") +
  ylim(0, 100) +
  labs( 
    title = "Unique promoter peaks of G4 cluster 4 (-1.5/+1.5 kb TSS)",
    x = "scRNA-Seq (Marques et al.) Seurat cluster",
    y = "log normalized expr.",
    fill = "cell type"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  ) +
  stat_compare_means(label.y = 75, label.x = 3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 50) 
cluster_4_bp

## cluster 5
cluster_5 = lapply(cell_types, get_expression, g4_cluster = "cluster_5")
cluster_5 = do.call(rbind, cluster_5)

cluster_5_bp = ggplot(cluster_5, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Reds") +
  ylim(0, 100) +
  labs(
    title = "Unique promoter peaks of G4 cluster 5 (-1.5/+1.5 kb TSS)",
    x = "scRNA-Seq (Marques et al.) Seurat cluster",
    y = "log normalized expr.",
    fill = "cell type"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  ) +
  stat_compare_means(label.y = 75, label.x = 3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 50) 
cluster_5_bp

bp_grid = plot_grid(cluster_0_bp, cluster_3_bp, cluster_4_bp, cluster_5_bp)

ggsave(
  glue("{result_folder}Unique_prom_G4s-expr_boxplot.png"),
  plot = bp_grid,
  width = 10,
  height = 10,
  dpi = 300,
)


