# packages
suppressPackageStartupMessages({
  library("data.table")
  library("ComplexHeatmap")
  library("glue")
  library("tidyverse")
  library("ComplexHeatmap")
  library("EnhancedVolcano")
  library("patchwork")
  library("circlize")
  library("Seurat")
  library("Matrix")
})

# folders
result_folder = "../results/Seurat/callpeaks_unsorted/"
marker_genes = "../data/GSE75330/marker_genes.txt" # Marques et al. scRNA-Seq marker genes
marques_marker_genes = fread(marker_genes, header = FALSE)

# scRNA-Seq data
scrna = "../results/Seurat/scRNASeq_GSM4979874-75.rds"
marques_scrna = "../results/Seurat/scRNASeq_GSE75330.rds"

markers = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_output_annot.tsv") # it comes from "seurat_peakcalls*" script following HOMER annot
markers[markers == ""] <- NA
markers = markers %>% na.omit()

# Grubb's test output
grubbs = fread(
  "../results/Seurat/callpeaks_unsorted/Grubbs_test-unique_G4_peaks_0.001_joined.tsv"
)
print(glue("Grubbs test found {as.character(length(unique(grubbs$`Gene Name`)))} unique peak"))
grubbs = grubbs %>% dplyr::filter(abs(`Distance to TSS`) < 3000)

# filter and export to bed
clusterwise = markers %>% dplyr::filter(abs(avg_log2FC) > 0.2 &
                                          p_val_adj < 0.05)
clusters = unique(clusterwise$cluster)
for (i in clusters) {
  subset = clusterwise %>% dplyr::filter(cluster == i) %>% dplyr::select(chr, start, end)
  write_tsv(
    subset,
    glue(
      "{result_folder}FindAllMarker_padj0.05_fc0.2_cluster_{i}.bed"
    ),
    col_names = FALSE
  )
}

# export clusters to tsv
for (i in clusters) {
  subset = markers %>% dplyr::filter(cluster == i) %>%
    dplyr::select(`Gene Name`, avg_log2FC, p_val_adj, `Distance to TSS`) %>%
    dplyr::filter(abs(`Distance to TSS`) < 3000)
  write_tsv(
    subset,
    glue("{result_folder}FindAllMarker_cluster_{i}_promoters.tsv"),
    col_names = TRUE
  )
}

# keep significant marker regions proximal to TSS
markers = markers %>% dplyr::filter(abs(`Distance to TSS`) < 3000) %>%
  dplyr::filter(p_val_adj < 0.001)
print(glue("FindAllMarkers found {as.character(length(unique(markers$`Gene Name`)))} promoter proximal G4 markers"))
markers_grubbs_int = intersect(unique(markers$`Gene Name`), unique(grubbs$`Gene Name`))
print(
  glue(
    "{as.character(length(markers_grubbs_int))} of {as.character(dim(grubbs)[1])} outliers have been found by Seurat"
  )
)

markers_all = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_output_annot.tsv")
setdiff(unique(grubbs[which(unique == "4")]$`Gene Name`), unique(markers[which(cluster == "4")]$`Gene Name`) 
          )


markers_wide = pivot_wider(
  markers,
  names_from = cluster,
  id_cols = "Gene Name",
  values_from = avg_log2FC,
  values_fn = mean
)
row_names = markers_wide %>% pull("Gene Name")
markers_wide = markers_wide %>% dplyr::select(-"Gene Name") %>% as.matrix
rownames(markers_wide) = row_names

markers_wide[is.na(markers_wide)] = 0

G4_prox_marker_genes = marques_marker_genes$V1[which(marques_marker_genes$V1 %in% rownames(markers_wide))]
G4_prox_marker_genes = append(G4_prox_marker_genes, "Tcn2")
G4_prox_marker_genes_pos = which(rownames(markers_wide) %in% G4_prox_marker_genes)

pdf(
  file = glue("{result_folder}FindAllMarkers_hm_p0.001.pdf"),
  width = 4,
  height = 6
)
ha = rowAnnotation(foo = anno_mark(at = G4_prox_marker_genes_pos, labels = rownames(markers_wide)[G4_prox_marker_genes_pos]))
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "#fec44f", "red"))
hm = Heatmap(
  markers_wide,
  column_title = "marker regions, adj. p < 0.001",
  row_title = "genes close to G4 (+/- 3 kb to TSS)",
  name = "Avg logFC",
  right_annotation = ha,
  col = col_fun,
  heatmap_legend_param = list(
    at = c(0, 0.5, 1),
    labels = c("< 0.25", "0.5", "1"),
    legend_height = unit(2, "cm")
  ),
  #rect_gp = gpar(col = "black", lwd = 0.01),
  show_column_dend = FALSE,
  show_row_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(13, "cm"),
  row_names_gp = gpar(fontsize = 1),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 0
)
hm
dev.off()

png(
  file = glue("{result_folder}FindAllMarkers_hm_p0.001.png"),
  width = 10,
  height = 12,
  units = "cm",
  res = 300
)
print(hm)
dev.off()

# integrate with scRNA-Seq data
# GSE75330 (Marques et al.)
marques = readRDS(marques_scrna)
# log normalized expression matrix
norm = marques[["RNA"]]@data

existing_gene_symbols = character()
for (gene in markers$`Gene Name`) {
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

ms = ms %>% inner_join(., markers, by = c("gene_symbol" = "Gene Name")) %>% dplyr::filter(abs(avg_log2FC) > 0.5) %>%
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% dplyr::select("MOL":"VLMC")
ms = as.matrix(ms)
rownames(ms) = rows

pdf(
  file = glue(
    "{result_folder}heatmap_unsorted-marker_genes-scRNA_GSE75330.pdf"
  ),
  width = 5,
  height = 5
)
col_fun = colorRamp2(c(0, 1, 2), c("#9ecae1", "white", "#fc9272"))
marques_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Marques et al.) Seurat clusters",
  row_title = "G4 marker genes (avg log2FC > 0.5)",
  name = "norm. expr.",
  clustering_method_rows = "complete",
  col = col_fun,
  #top_annotation = ha,
  show_column_dend = TRUE,
  #rect_gp = gpar(col = "black", lwd = 0.2),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  # heatmap_width = unit(12, "cm"),
  # heatmap_height = unit(12, "cm"),
  width = unit(3, "cm"),
  height = unit(6, "cm"),
  show_row_names = FALSE,
  #row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 13),
  column_names_rot = 90
)
marques_hm
dev.off()

### Volcano and scatterplot visualization

# keep only those G4 markers that markers in one given Seurat cluster
cluster_spec_markers = markers %>% distinct(`Gene Name`, cluster, .keep_all = TRUE) %>%
  group_by(`Gene Name`) %>% count(`Gene Name`) %>% filter(n == 1) %>% pull(`Gene Name`) %>% unique
cluster_spec_markers = markers %>% filter(`Gene Name` %in% cluster_spec_markers)
sign_g4_markers = cluster_spec_markers %>% pull(`Gene Name`) %>% unique

## run FindAllMarkers on Marques scRNA-Seq object (default is wilcoxon test)
# marques_markers = FindAllMarkers(marques,
#                                  logfc.threshold = 0.2)
# write_tsv(marques_markers,
#           "../data/GSE75330/FindAllMarker_fc0.2_wilcox.tsv")
marques_markers = read_tsv("../data/GSE75330/FindAllMarker_fc0.2_wilcox.tsv")
marques_markers_filt = marques_markers %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 2)

# combine with our G4 marker genes
all_markers = marques_markers_filt %>% inner_join(., cluster_spec_markers, by = c("gene" = "Gene Name")) %>%
  dplyr::select(
    gene_name = gene,
    scRNA_Seq_fc = avg_log2FC.x,
    G4_fc = avg_log2FC.y,
    scRNA_Seq_padj = p_val_adj.x,
    G4_padj = p_val_adj.y,
    scRNA_Seq_cluster = cluster.x,
    G4_cluster = cluster.y
  )

cluster_quant = all_markers %>% group_by(scRNA_Seq_cluster) %>% count(scRNA_Seq_cluster)

# volcano plots
volcano = function(cluster, color = "black", title) {
  # keyvalues = ifelse(
  #   cluster$G4_cluster == 4, '#8bb0d1',
  #   ifelse(cluster$G4_cluster == 5, '#f0b265',
  #          '#b5db68'))
  # keyvalues[is.na(keyvalues)] <- 'black'
  # names(keyvalues)[keyvalues == '#8bb0d1'] <- '4'
  # names(keyvalues)[keyvalues == '#f0b265'] <- '5'
  # names(keyvalues)[keyvalues == '#b5db68'] <- '6'
  
  plot = EnhancedVolcano(
    cluster,
    lab = cluster$gene_name,
    labSize = 4,
    title = title,
    pCutoff = 10e-25,
    FCcutoff = 0.5,
    colAlpha = 1,
    legendLabels = c(" ", " ", " ", " "),
    legendLabSize = 14,
    legendIconSize = 0,
    legendDropLevels = TRUE,
    legendPosition = "right",
    drawConnectors = TRUE,
    max.overlaps = 100,
    #colCustom = keyvalues,
    x = 'G4_fc',
    y = 'G4_padj',
    subtitle = "",
    col = c('grey', 'grey', 'grey', color)
  )
  return(print(plot))
  
}

print(unique(all_markers$scRNA_Seq_cluster))

# volcano plots of each Marques et al. scRNA-Seq cluster
pdf(
  file = glue("{result_folder}volcano_MOL-G4_fc.pdf"),
  width = 6,
  height = 6
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "MOL" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "MOL"
)
dev.off()
png(
  file = glue("{result_folder}volcano_MOL-G4_fc.png"),
  width = 12,
  height = 12,
  units = "cm",
  res = 300
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "MOL" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "MOL"
)
dev.off()
pdf(
  file = glue("{result_folder}volcano_MFOL-MOL-G4_fc.pdf"),
  width = 6,
  height = 6
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "MFOL/MOL" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "MFOL/MOL"
)
dev.off()
png(
  file = glue("{result_folder}volcano_MFOL-MOL-G4_fc.png"),
  width = 12,
  height = 12,
  units = "cm",
  res = 300
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "MFOL/MOL" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "MFOL/MOL"
)
dev.off()
pdf(
  file = glue("{result_folder}volcano_COP-G4_fc.pdf"),
  width = 6,
  height = 6
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "COP" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "COP"
)
dev.off()
png(
  file = glue("{result_folder}volcano_COP-G4_fc.png"),
  width = 12,
  height = 12,
  units = "cm",
  res = 300
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "COP" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "COP"
)
dev.off()
pdf(
  file = glue("{result_folder}volcano_VLMC-G4_fc.pdf"),
  width = 6,
  height = 6
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "VLMC" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "VLMC"
)
dev.off()
png(
  file = glue("{result_folder}volcano_VLMC-G4_fc.png"),
  width = 12,
  height = 12,
  units = "cm",
  res = 300
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "VLMC" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "VLMC"
)
dev.off()
pdf(
  file = glue("{result_folder}volcano_OPC-G4_fc.pdf"),
  width = 6,
  height = 6
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "OPC" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "OPC"
)
dev.off()
png(
  file = glue("{result_folder}volcano_OPC-G4_fc.png"),
  width = 12,
  height = 12,
  units = "cm",
  res = 300
)
volcano(
  all_markers %>% dplyr::filter(scRNA_Seq_cluster == "OPC" &
                                  G4_cluster == 5),
  color = "#f0b265",
  title = "OPC"
)
dev.off()

# lower scRNA-Seq fold change thr and color by G4 Seurat clusters
marques_markers_filt = marques_markers %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>%
  dplyr::filter(gene %in% sign_g4_markers)

all_markers = marques_markers_filt %>% inner_join(., markers, by = c("gene" = "Gene Name")) %>%
  dplyr::select(
    gene_name = gene,
    scRNA_Seq_fc = avg_log2FC.x,
    G4_fc = avg_log2FC.y,
    scRNA_Seq_padj = p_val_adj.x,
    G4_padj = p_val_adj.y,
    scRNA_Seq_cluster = cluster.x,
    G4_cluster = cluster.y
  )

keyvalues = ifelse(
  all_markers$G4_cluster == 4,
  '#8bb0d1',
  ifelse(
    all_markers$G4_cluster == 5,
    '#f0b265',
    ifelse(
      all_markers$G4_cluster == 6,
      '#b5db68',
      ifelse(
        all_markers$G4_cluster == 3,
        '#EB8073',
        ifelse(
          all_markers$G4_cluster == 2,
          '#BDB8D8',
          ifelse(
            all_markers$G4_cluster == 1,
            '#FDFFB2',
            ifelse(all_markers$G4_cluster == 0, '#99D1C4',
                   'grey')
          )
        )
      )
    )
  )
)
keyvalues[is.na(keyvalues)] <- 'grey'
names(keyvalues)[keyvalues == '#99D1C4'] <- '0'
names(keyvalues)[keyvalues == '#FDFFB2'] <- '1'
names(keyvalues)[keyvalues == '#BDB8D8'] <- '2'
names(keyvalues)[keyvalues == '#EB8073'] <- '3'
names(keyvalues)[keyvalues == '#8bb0d1'] <- '4'
names(keyvalues)[keyvalues == '#f0b265'] <- '5'
names(keyvalues)[keyvalues == '#b5db68'] <- '6'

# G4 fold change volcano
pdf(
  file = glue("{result_folder}volcano_full-G4_fc.pdf"),
  width = 6,
  height = 6
)
full_volc_g4 = EnhancedVolcano(
  all_markers,
  lab = all_markers$gene_name,
  labSize = 4,
  title = "|scRNA-Seq log2FC| > 1, adj. p. < 0.05" ,
  pCutoff = 10e-45,
  FCcutoff = 0.5,
  colAlpha = 1,
  selectLab = " ",
  legendLabels = c(" ", " ", " ", " "),
  legendLabSize = 14,
  legendIconSize = 5,
  legendDropLevels = FALSE,
  legendPosition = "right",
  drawConnectors = TRUE,
  max.overlaps = 100,
  
  colCustom = keyvalues,
  x = 'G4_fc',
  y = 'G4_padj',
  subtitle = "",
  col = c('grey', 'grey', 'grey', color)
)
full_volc_g4
dev.off()

png(
  file = glue("{result_folder}volcano_full-G4_fc.png"),
  width = 16,
  height = 12,
  units = "cm",
  res = 300
)
print(full_volc_g4)
dev.off()

# scRNA-Seq fold change volcano
pdf(
  file = glue("{result_folder}volcano_full-scRNA-Seq_fc.pdf"),
  width = 6,
  height = 6
)
full_volc_scrna = EnhancedVolcano(
  all_markers,
  lab = all_markers$gene_name,
  labSize = 4,
  title = "|scRNA-Seq log2FC| > 1, adj. p. < 0.05" ,
  pCutoff = 10e-45,
  FCcutoff = 1,
  colAlpha = 1,
  selectLab = " ",
  legendLabels = c(" ", " ", " ", " "),
  legendLabSize = 14,
  legendIconSize = 5,
  legendDropLevels = FALSE,
  legendPosition = "right",
  drawConnectors = TRUE,
  max.overlaps = 100,
  colCustom = keyvalues,
  x = 'scRNA_Seq_fc',
  y = 'scRNA_Seq_padj',
  subtitle = "",
  col = c('grey', 'grey', 'grey', color)
)
full_volc_scrna
dev.off()

png(
  file = glue("{result_folder}volcano_full-scRNA-Seq_fc.png"),
  width = 16,
  height = 12,
  units = "cm",
  res = 300
)
print(full_volc_scrna)
dev.off()

# scatterplot
all_markers = marques_markers %>% inner_join(., markers, by = c("gene" = "Gene Name")) %>%
  dplyr::select(
    gene_name = gene,
    scRNA_Seq_fc = avg_log2FC.x,
    G4_fc = avg_log2FC.y,
    scRNA_Seq_padj = p_val_adj.x,
    G4_padj = p_val_adj.y,
    scRNA_Seq_cluster = cluster.x,
    G4_cluster = cluster.y
  )

ggplot(all_markers, aes(
  x = scRNA_Seq_fc,
  y = G4_fc,
  color = factor(G4_cluster)
)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "scRNA-Seq vs. G4 markers",
    x = "Marques et al. scRNA-Seq fold change",
    y = "G4 fold change",
    color = " "
  ) +
  theme_classic() +
  labs(color = "Seurat cluster") +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )

ggsave(
  plot = last_plot(),
  glue("{result_folder}FindMarker_anal-fc_scatterplot.pdf"),
  width = 6,
  height = 6
)
ggsave(
  plot = last_plot(),
  glue("{result_folder}FindMarker_anal-fc_scatterplot.png"),
  width = 6,
  height = 6,
  dpi = 300
)

ggplot(all_markers, aes(
  x = scRNA_Seq_fc,
  y = G4_fc,
  color = factor(scRNA_Seq_cluster)
)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "scRNA-Seq vs. G4 markers",
    x = "Marques et al. scRNA-Seq fold change",
    y = "G4 fold change",
    color = " "
  ) +
  theme_classic() +
  labs(color = "cell type") +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )

ggsave(
  plot = last_plot(),
  glue("{result_folder}FindMarker_anal-fc_scatterplot2.pdf"),
  width = 6,
  height = 6
)

# GSE163484, GSE163485 (two replicates)
scrna = readRDS(scrna)
norm = scrna[["RNA"]]@data

existing_gene_symbols = character()
for (gene in markers$`Gene Name`) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# heatmap
ms = list()
for (cluster in levels(unique(Idents(scrna)))) {
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

ms = ms %>% inner_join(., markers, by = c("gene_symbol" = "Gene Name")) %>% dplyr::filter(abs(avg_log2FC) > 0.5) %>%
  distinct(gene_symbol, .keep_all = TRUE)
rows = ms$gene_symbol
ms = ms %>% relocate(.) %>% dplyr::select("0":"9")
ms = as.matrix(ms)
rownames(ms) = rows

pdf(
  file = glue(
    "{result_folder}heatmap_unsorted-marker_genes-scRNA_GSE163484-85.pdf"
  ),
  width = 5,
  height = 5
)
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
bartosovic_hm = Heatmap(
  ms,
  column_title = "scRNA-Seq (Bartosovic et al.) Seurat clusters",
  row_title = "G4 marker genes (avg log2FC > 0.5)",
  name = "norm. expr.",
  clustering_method_rows = "complete",
  col = col_fun,
  #top_annotation = ha,
  show_column_dend = FALSE,
  #rect_gp = gpar(col = "black", lwd = 0.2),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  # heatmap_width = unit(12, "cm"),
  # heatmap_height = unit(12, "cm"),
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  show_row_names = FALSE,
  #row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 13),
  column_names_rot = 90
)
bartosovic_hm
dev.off()

## run FindAllMarkers on Bartosovic scRNA-Seq object (default is wilcoxon test)
#bartosovic_markers = FindAllMarkers(scrna, logfc.threshold = 0.2)
#write_tsv(bartosovic_markers, "../data/GSE163484/FindAllMarker_fc0.2_wilcox.tsv")
