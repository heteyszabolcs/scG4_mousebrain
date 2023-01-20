suppressPackageStartupMessages({
  library("tidyverse")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("Seurat")
  library("matrixStats")
  library("ComplexHeatmap")
  library("circlize")
  library("outliers")
})

set.seed(42)

# result / peak folder
result_folder = "../results/Seurat/callpeaks_unsorted/"
bigwig_folder = "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/"
peak_sets = "../results/Seurat/callpeaks_unsorted/peak_sets/"
markers = "../data/GSE75330/marker_genes.txt"

scrna = "../results/Seurat/scRNASeq_GSM4979874-75.rds"
marques_scrna = "../results/Seurat/scRNASeq_GSE75330.rds"

cluster_bigwigs = list.files(bigwig_folder, pattern = "*[0-9].bw", full.names = TRUE)
peak_set = "enhancer_analysis_output.bed"
bed = glue("{peak_sets}{peak_set}")

# create read coverage matrix
read_cov = bw_loci(cluster_bigwigs, loci = bed)
read_cov = as.data.frame(read_cov)
mat = read_cov %>% dplyr::select(starts_with("X")) %>% dplyr::select(
  "0" = X0,
  "1" = X1,
  "2" = X2,
  "3" = X3,
  "4" = X4,
  "5" = X5,
  "6" = X6
)
mat = as.matrix(mat)

# Grubbs test to find outliers
# The Grubbs test allows to detect whether the highest or lowest value in a dataset is an outlier.
grubbs.ps = numeric()
for (i in 1:nrow(mat)) {
  p = unname(grubbs.test(mat[i,])$p.value)
  grubbs.ps = c(grubbs.ps, p)
}
mat_grubbs = mat[which(grubbs.ps < 0.001, arr.ind = TRUE), ]
read_cov_grubbs = read_cov[which(grubbs.ps < 0.001, arr.ind = TRUE), ]
bed_grubbs = read_cov_grubbs[, 1:3]
write_tsv(
  bed_grubbs,
  glue("{result_folder}Grubbs_test-unique_G4_peaks_0.001.bed"),
  col_names = FALSE
)

mat_log = log(mat_grubbs + 1)
colnames(mat_log) = c(
  "0",
  "1",
  "2",
  "3",
  "4",
  "5",
  "6"
)

pdf(
  file = glue(
    "{result_folder}read_density_hm-Signac_MACS2_peaks-Grubbs0.001.pdf"
  ),
  width = 4,
  height = 6
)
col_fun = colorRamp2(c(0, 60, 120), c("white","#fec44f", "#f03b20"))
unique_hm = Heatmap(
  mat_grubbs,
  column_title = " ",
  row_title = "Signac MACS2 peak (Grubb's test, p < 0.001)",
  name = "read cov",
  row_km = 2,
  #clustering_method_rows = "complete",
  col = col_fun,
  # rect_gp = gpar(col = "black", lwd = 0.0),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(13, "cm"),
  row_names_gp = gpar(fontsize = 3),
  column_names_rot = 0
)
unique_hm
dev.off()

png(
  file = glue("{result_folder}read_density_hm-Signac_MACS2_peaks-Grubbs0.001.png"),
  width = 15,
  height = 15,
  units = "cm",
  res = 300
)
print(unique_hm)
dev.off()

tibble = as_tibble(mat_log)
tibble = tibble %>% mutate(max = pmax(rownames(tibble)), type = colnames(tibble)[max.col(tibble[, 1:7])])
mat_grubbs = as_tibble(mat_grubbs)
mat_grubbs_max = mat_grubbs %>% mutate(max = pmax(rownames(mat_grubbs)), type = colnames(mat_grubbs)[max.col(mat_grubbs[, 1:7])]) %>%
  mutate(
    seqnames = read_cov_grubbs$seqnames,
    start = read_cov_grubbs$start,
    end = read_cov_grubbs$end,
    width = read_cov_grubbs$width
  ) %>%
  dplyr::select(seqnames, start, end, width, starts_with("cluster"), unique = type) 


n_uniques = tibble %>% group_by(type) %>% summarise(count = dplyr::n()) %>% arrange(desc(count)) %>%
  ggplot(data = ., aes(
    x = reorder(type, -count),
    y = count,
    fill = type
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "Reds") +
  labs(title = "# of unique G4 locations (Grubb's test, p < 0.001)",
       x = "Seurat cluster",
       y = "# of G4 structures",
       fill = "Seurat cluster") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
n_uniques

ggsave(
  glue("{result_folder}Grubbs_test-unique_G4_peaks_bar.png"),
  plot = n_uniques,
  width = 10,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Grubbs_test-unique_G4_peaks_bar.pdf"),
  plot = n_uniques,
  width = 10,
  height = 5,
  device = "pdf"
)

# if HOMER annotation is available
read_cov_grubbs_annot = "../results/Seurat/callpeaks_unsorted/Grubbs_test-unique_G4_peaks_0.001_annot.tsv"
read_cov_grubbs_annot = fread(read_cov_grubbs_annot)
read_cov_grubbs_annot = read_cov_grubbs_annot %>% dplyr::select(rowid = starts_with("PeakID"), everything()) %>% 
  mutate(rowid = row_number())

read_cov_grubbs_annot = read_cov_grubbs %>% mutate(rowid = row_number()) %>% 
  inner_join(., read_cov_grubbs_annot, by = c("rowid" = "rowid")) %>%
  dplyr::select(
    "cluster 0" = X0,
    "cluster 1" = X1,
    "cluster 2" = X2,
    "cluster 3" = X3,
    "cluster 4" = X4,
    "cluster 5" = X5,
    "cluster 6" = X6,
    `Distance to TSS`,
    `Gene Name`
  ) %>% 
  mutate(rowid = row_number())

cluster6 = mat_grubbs_max %>% mutate(rowid = row_number()) %>%
  inner_join(., read_cov_grubbs_annot, by = c("rowid" = "rowid")) %>%
  dplyr::filter(unique == "6") %>%
  dplyr::select(Chr = "seqnames", Start = "start", End = "end")

write_tsv(cluster6,
          glue("{result_folder}Grubbs_test-unique_G4_cluster6.bed"),
          col_names = FALSE)

# inner join with annotation - be careful, here some genes will be omitted as they have not been annotated
mat_grubbs_max = mat_grubbs_max %>% 
  mutate(rowid = row_number()) %>% inner_join(., read_cov_grubbs_annot, by = c("rowid" = "rowid")) %>%
  dplyr::select(
    "cluster 0",
    "cluster 1",
    "cluster 2",
    "cluster 3",
    "cluster 4",
    "cluster 5",
    "cluster 6",
    unique,
    `Distance to TSS`,
    `Gene Name`
  )

write_tsv(
  mat_grubbs_max,
  glue(
    "{result_folder}Grubbs_test-unique_G4_peaks_0.001_joined.tsv"
  ),
  col_names = TRUE
)

# integrate with scRNA-Seq data
# GSE163484, GSE163485 (two replicates)
scrna = readRDS(scrna)
norm = scrna[["RNA"]]@data

existing_gene_symbols = character()
for (gene in read_cov_grubbs_annot$`Gene Name`) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# heatmap
ms = list()
for (cluster in seq(0, 15)) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(scrna, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}
ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name"))
rows = ms$gene_symbol
ms = ms %>% dplyr::select("0":"15")
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_peaks-scRNA_GSE163484-85.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
bartosovic_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Bartosovic et al.) Seurat clusters",
  row_title = "unique G4 location",
  name = "norm. expr.",
  row_km = 3,
  column_km = 1,
  #clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  heatmap_width = unit(6, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 0
)
bartosovic_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_peaks-scRNA_GSE163484-85.pdf"
  ),
  width = 5,
  height = 5,
)
print(bartosovic_hm)
dev.off()

# integrate with oligodendrocyte markers (Marques et al.)
markers = fread(markers, header = FALSE) # marker genes from Zhang et al. and Marques et al.

existing_gene_symbols = character()
for (gene in markers$V1) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# heatmap
hm_markers = list()
for (cluster in seq(0, 15)) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(scrna, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  hm_markers[[as.character(cluster)]] = t
}
hm_markers = bind_rows(hm_markers)
hm_markers = hm_markers %>% distinct_all(., .keep_all = TRUE)
hm_markers = hm_markers %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

rows = hm_markers$gene_symbol
hm_markers = hm_markers %>% dplyr::select("0":"15")
hm_markers = as.matrix(hm_markers)
rownames(hm_markers) = rows

png(
  file = glue(
    "{result_folder}heatmap_oligo_markers-scRNA_GSE163484-85.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
bartosovic_marker_hm = Heatmap(
  hm_markers,
  column_title = "scRNA-Seq (Bartosovic et al.) Seurat clusters",
  row_title = "oligodendrocyte marker gene (Marques et al.)",
  name = "norm. expr.",
  row_km = 2,
  column_km = 2,
  rect_gp = gpar(col = "black", lwd = 0.3),
  #clustering_method_rows = "complete",
  col = col_fun,
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  heatmap_width = unit(12, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 5),
  column_names_rot = 0
)
bartosovic_marker_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_oligo_markers-scRNA_GSE163484-85.pdf"
  ),
  width = 6,
  height = 5,
)
print(bartosovic_marker_hm)
dev.off()

# integrate with scRNA-Seq data
# GSE75330 (Marques et al.)
marques = readRDS(marques_scrna)
# log normalized expression matrix
norm = marques[["RNA"]]@data

existing_gene_symbols = character()
for (gene in read_cov_grubbs_annot$`Gene Name`) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# heatmap
ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}
ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name"))
rows = ms$gene_symbol
ms = ms %>% dplyr::select("MOL":"VLMC")
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_peaks-scRNA_GSE75330.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
marques_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location",
  name = "norm. expr.",
  row_km = 3,
  column_km = 1,
  #clustering_method_rows = "complete",
  col = col_fun,
  #top_annotation = ha,
  show_column_dend = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.2),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  # heatmap_width = unit(12, "cm"),
  # heatmap_height = unit(12, "cm"),
  width = unit(2, "cm"),
  height = unit(14, "cm"),
  row_names_gp = gpar(fontsize = 1.7),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
marques_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_peaks-scRNA_GSE75330.pdf"
  ),
  width = 4.5,
  height = 7,
)
print(marques_hm)
dev.off()


# Marques et al. expressions near TSS
read_cov_grubbs_tss = read_cov_grubbs_annot %>% dplyr::filter(abs(`Distance to TSS`) < 3000)

existing_gene_symbols = character()
for (gene in read_cov_grubbs_tss$`Gene Name`) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# heatmap
ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}
ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name"))
rows = ms$gene_symbol
ms = ms %>% dplyr::select("MOL":"VLMC")
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_promoter_peaks-scRNA_GSE75330.png"
  ),
  width = 18,
  height = 19,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
marques_proms_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location (+3/-3 kb TSS)",
  name = "norm. expr.",
  row_km = 3,
  column_km = 1,
  #clustering_method_rows = "complete",
  col = col_fun,
  #top_annotation = ha,
  show_column_dend = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.3),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  # heatmap_width = unit(12, "cm"),
  # heatmap_height = unit(12, "cm"),
  width = unit(8, "cm"),
  height = unit(19, "cm"),
  row_names_gp = gpar(fontsize = 3.5),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 90
)
marques_proms_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_promoter_peaks-scRNA_GSE75330.pdf"
  ),
  width = 8,
  height = 10
)
print(marques_proms_hm)
dev.off()

# integrate with oligodendrocyte markers (Marques et al.)
# GSE75330 (Marques et al.)
existing_gene_symbols = character()
for (gene in markers$V1) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# heatmap
hm_markers = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  hm_markers[[as.character(cluster)]] = t
}
hm_markers = bind_rows(hm_markers)
hm_markers = hm_markers %>% distinct_all(., .keep_all = TRUE)
hm_markers = hm_markers %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

rows = hm_markers$gene_symbol
hm_markers = hm_markers %>% select(levels(unique(Idents(marques))))
hm_markers = as.matrix(hm_markers)
rownames(hm_markers) = rows

png(
  file = glue("{result_folder}heatmap_oligo_markers-scRNA_GSE75330.png"),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
marques_markers_hm = Heatmap(
  hm_markers,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "oligodendrocyte marker gene (Marques et al.)",
  name = "norm. expr.",
  row_km = 3,
  column_km = 0,
  rect_gp = gpar(col = "black", lwd = 0.2),
  # clustering_method_rows = "complete",
  # clustering_method_columns = "complete",
  col = col_fun,
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(2, "cm"),
  height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
marques_markers_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_oligo_markers-scRNA_GSE75330.pdf"
  ),
  width = 6,
  height = 6
)
print(marques_markers_hm)
dev.off()

# integrate with scRNA-Seq data
# GSE75330 (Marques et al.)
existing_gene_symbols = character()
for (gene in read_cov_grubbs_annot$`Gene Name`) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# heatmap
ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}
ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name")) %>% 
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% select(levels(unique(Idents(marques))))
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_peaks-scRNA_GSE75330.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
unique_marques_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location",
  name = "norm. expr.",
  row_km = 3,
  column_km = 1,
  #clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.2),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(1.5, "cm"),
  height = unit(16, "cm"),
  row_names_gp = gpar(fontsize = 2.5),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
unique_marques_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_peaks-scRNA_GSE75330_.pdf"
  ),
  width = 6,
  height = 7.5
)
print(unique_marques_hm)
dev.off()

unique_marques_hm2 = Heatmap(
  ms,
  # column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location",
  name = "norm. expr.",
  row_km = 3,
  column_km = 1,
  show_row_names = FALSE,
  #clustering_method_rows = "complete",
  col = col_fun,
  # rect_gp = gpar(col = "black", lwd = 0.2),
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(1.5, "cm"),
  height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 2.5),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
unique_marques_hm2

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_peaks-scRNA_GSE75330_2.pdf"
  ),
  width = 4,
  height = 6
)
print(unique_marques_hm2)
dev.off()

# heatmaps for G4 Seurat clusters indicating Marques et al. expression levels
# cluster 2
# cluster2_unique = mat_grubbs_max %>% filter(unique == "cluster_2") %>% pull("Gene Name")
#
# existing_gene_symbols = character()
# for (gene in cluster2_unique) {
#   if (gene %in% norm@Dimnames[[1]]) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
#
# ms = list()
# for (cluster in levels(unique(Idents(marques)))) {
#   m = t(as.matrix(norm[existing_gene_symbols,]))
#   cluster_barcodes = WhichCells(marques, idents = cluster)
#   m = m[cluster_barcodes,]
#   t = tibble(means = colMeans(m))
#   t = t %>%
#     mutate(gene_symbol = existing_gene_symbols) %>%
#     mutate(Seurat_cluster = as.character(cluster))
#   ms[[as.character(cluster)]] = t
# }
#
# ms = bind_rows(ms)
# ms = ms %>% distinct_all(., .keep_all = TRUE)
# ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")
#
# ms = ms %>% inner_join(., read_cov_grubbs, by = c("gene_symbol" = "Gene Name"))
# rows = ms$gene_symbol
# ms = ms %>% select(levels(unique(Idents(marques))))
# ms = as.matrix(ms)
# rownames(ms) = rows
#
# png(
#   file = glue("{result_folder}heatmap_uniqueMACS2_cl2-scRNA_GSE75330.png"),
#   width = 15,
#   height = 15,
#   units = 'cm',
#   res = 500
# )
# col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
# Heatmap(
#   ms,
#   column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
#   row_title = "unique G4 location - cluster 2",
#   name = "norm. expr.",
#   row_km = 1,
#   column_km = 1,
#   rect_gp = gpar(col = "black", lwd = 0.2),
#   #clustering_method_rows = "complete",
#   col = col_fun,
#   #top_annotation = ha,
#   show_column_dend = TRUE,
#   cluster_columns = FALSE,
#   cluster_rows = TRUE,
#   show_row_dend = FALSE,
#   width = unit(10, "cm"),
#   height = unit(5, "cm"),
#   row_names_gp = gpar(fontsize = 8),
#   column_names_gp = gpar(fontsize = 12),
#   column_names_rot = 90
# )
# dev.off()

# cluster 3
cluster3_unique = mat_grubbs_max %>% filter(unique == "3") %>% pull("Gene Name")

existing_gene_symbols = character()
for (gene in cluster3_unique) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}

ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name")) %>% 
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% select(levels(unique(Idents(marques))))
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl3-scRNA_GSE75330.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
cl_3_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location - cluster 3",
  name = "norm. expr.",
  row_km = 2,
  #column_km = 1,
  #clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.2),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(3, "cm"),
  height = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 90
)
cl_3_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl3-scRNA_GSE75330.pdf"
  ),
  width = 5.5,
  height = 5,
)
print(cl_3_hm)
dev.off()

# cluster 4
cluster4_unique = mat_grubbs_max %>% filter(unique == "4") %>% pull("Gene Name")

existing_gene_symbols = character()
for (gene in cluster4_unique) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}

ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name")) %>% 
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% select(levels(unique(Idents(marques))))
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl4-scRNA_GSE75330.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
cl_4_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location - cluster 4",
  name = "norm. expr.",
  row_km = 1,
  #column_km = 2,
  #clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.2),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(3, "cm"),
  height = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
cl_4_hm
dev.off()

# cluster 1
# cluster1_unique = mat_grubbs_max %>% filter(unique == "cluster_1") %>% pull("Gene Name")
#
# existing_gene_symbols = character()
# for (gene in cluster1_unique) {
#   if (gene %in% norm@Dimnames[[1]]) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
#
# ms = list()
# for (cluster in levels(unique(Idents(marques)))) {
#   m = t(as.matrix(norm[existing_gene_symbols,]))
#   cluster_barcodes = WhichCells(marques, idents = cluster)
#   m = m[cluster_barcodes,]
#   t = tibble(means = colMeans(m))
#   t = t %>%
#     mutate(gene_symbol = existing_gene_symbols) %>%
#     mutate(Seurat_cluster = as.character(cluster))
#   ms[[as.character(cluster)]] = t
# }
#
# ms = bind_rows(ms)
# ms = ms %>% distinct_all(., .keep_all = TRUE)
# ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")
#
# ms = ms %>% inner_join(., read_cov_grubbs, by = c("gene_symbol" = "Gene Name"))
# rows = ms$gene_symbol
# ms = ms %>% select(levels(unique(Idents(marques))))
# ms = as.matrix(ms)
# rownames(ms) = rows
#
# png(
#   file = glue("{result_folder}heatmap_uniqueMACS2_cl1-scRNA_GSE75330.png"),
#   width = 15,
#   height = 15,
#   units = 'cm',
#   res = 500
# )
# col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
# Heatmap(
#   ms,
#   column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
#   row_title = "unique G4 location - cluster 1",
#   name = "norm. expr.",
#   row_km = 3,
#   #column_km = 2,
#   #clustering_method_rows = "complete",
#   col = col_fun,
#   rect_gp = gpar(col = "black", lwd = 0.2),
#   #top_annotation = ha,
#   show_column_dend = TRUE,
#   cluster_columns = FALSE,
#   cluster_rows = TRUE,
#   show_row_dend = FALSE,
#   width = unit(1.5, "cm"),
#   height = unit(10, "cm"),
#   row_names_gp = gpar(fontsize = 8),
#   column_names_gp = gpar(fontsize = 6),
#   column_names_rot = 90
# )
# dev.off()

# cluster 0
cluster0_unique = mat_grubbs_max %>% filter(unique == "0") %>% pull("Gene Name")

existing_gene_symbols = character()
for (gene in cluster0_unique) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}

ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name")) %>% 
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% select(levels(unique(Idents(marques))))
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl0-scRNA_GSE75330.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
cl_0_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location - cluster 0",
  name = "norm. expr.",
  row_km = 1,
  #column_km = 2,
  #clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.2),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(3, "cm"),
  height = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
cl_0_hm
dev.off()

# cluster 5
cluster5_unique = mat_grubbs_max %>% filter(unique == "5") %>% pull("Gene Name")

existing_gene_symbols = character()
for (gene in cluster5_unique) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols, ]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes, ]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}

ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name")) %>%
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% select(levels(unique(Idents(marques))))
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl5-scRNA_GSE75330.png"
  ),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
cl_5_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location - cluster 5",
  name = "norm. expr.",
  row_km = 4,
  #column_km = 2,
  #clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.2),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(1, "cm"),
  height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 5),
  column_names_rot = 90
)
cl_5_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl5-scRNA_GSE75330.pdf"
  ),
  width = 4.5,
  height = 5.5,
)
print(cl_5_hm)
dev.off()

cl_5_hm2 = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location - cluster 5",
  name = "norm. expr.",
  row_km = 4,
  #column_km = 2,
  #clustering_method_rows = "complete",
  col = col_fun,
  #rect_gp = gpar(col = "black", lwd = 0.2),
  #top_annotation = ha,
  show_column_dend = TRUE,
  show_row_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(1.5, "cm"),
  height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 5),
  column_names_rot = 90
)
cl_5_hm2

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl5-scRNA_GSE75330.pdf"
  ),
  width = 4.5,
  height = 5.5,
)
print(cl_5_hm2)
dev.off()


# cluster 6
cluster6_unique = mat_grubbs_max %>% filter(unique == "6") %>% pull("Gene Name")

existing_gene_symbols = character()
for (gene in cluster6_unique) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

ms = list()
for (cluster in levels(unique(Idents(marques)))) {
  m = t(as.matrix(norm[existing_gene_symbols,]))
  cluster_barcodes = WhichCells(marques, idents = cluster)
  m = m[cluster_barcodes,]
  t = tibble(means = colMeans(m))
  t = t %>%
    mutate(gene_symbol = existing_gene_symbols) %>%
    mutate(Seurat_cluster = as.character(cluster))
  ms[[as.character(cluster)]] = t
}

ms = bind_rows(ms)
ms = ms %>% distinct_all(., .keep_all = TRUE)
ms = ms %>% pivot_wider(., names_from = "Seurat_cluster", values_from = "means")

ms = ms %>% inner_join(., read_cov_grubbs_annot, by = c("gene_symbol" = "Gene Name")) %>%
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% select(levels(unique(Idents(marques))))
ms = as.matrix(ms)
rownames(ms) = rows

png(
  file = glue("{result_folder}heatmap_uniqueMACS2_cl6-scRNA_GSE75330.png"),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0, 3, 6), c("#9ecae1", "white", "#fc9272"))
cl_6_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "unique G4 location - cluster 6",
  name = "norm. expr.",
  row_km = 3,
  #column_km = 2,
  #clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.2),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  width = unit(1, "cm"),
  height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 5),
  column_names_rot = 90
)
cl_6_hm
dev.off()

pdf(
  file = glue(
    "{result_folder}heatmap_uniqueMACS2_cl6-scRNA_GSE75330.pdf"
  ),
  width = 5.5,
  height = 5.5,
)
print(cl_6_hm)
dev.off()
