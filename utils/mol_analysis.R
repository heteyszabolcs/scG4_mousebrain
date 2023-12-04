suppressPackageStartupMessages({
  library("tidyverse")
  library("Signac")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("Seurat")
  library("matrixStats")
  library("ComplexHeatmap")
  library("circlize")
  library("outliers")
  library("topGO")
})

# path to result folder
result_folder = "../results/Seurat/"

# helper functions
source("C:/Szabolcs/Karolinska/Data/scripts/get_unique_peaks_ws.R")
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# MOL markers
markers = read_tsv("../results/Seurat/final/sorted_brain/res0.1/integration/outputs/scRNA-Seq_Marques_et_al-FindAllMarkers_output.tsv")
markers = markers[markers$p_val < 0.05 &
                        markers$avg_log2FC > 0.5, ]

zeisel_markers = read_tsv("../data/Zeisel_et_al_2015/Zeisel_et_al_markers.txt", col_names = FALSE)
marques_markers = read_tsv("../data/GSE75330/marker_genes.txt", col_names = FALSE)
zeisel_markers = zeisel_markers %>% dplyr::filter(X2 == "Oligodendrocyte")
marques_markers = marques_markers %>% dplyr::filter(X2 == "MOL")
literature_markers = unique(c(marques_markers$X1, zeisel_markers$X1))

# anchor matrix of GFP+ - scRNA integration
anchors = read_tsv("../results/Seurat/final/sorted_brain/res0.1/integration/outputs/anchor_matrix.tsv")

# prediction scores
pred_score = readRDS("../results/Seurat/final/sorted_brain/res0.1/integration/outputs/g4_cell_label_preds.Rds")
pred_score = t(pred_score@data)
ids = rownames(pred_score)
pred_score = as_tibble(pred_score)
pred_score = pred_score %>% dplyr::select(-max) %>% mutate(cell_id = ids) %>% dplyr::filter(MOL > 0.75)
mol_ids = pred_score %>% pull(cell_id)

mol_markers = markers %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "MOL"), "MOL", cluster)) %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "NFOL"), "NFOL", cluster)) %>% 
  dplyr::filter(str_detect(cluster, "MOL")) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(500, wt = avg_log2FC)
all_mol_markers = mol_markers$gene
mol_markers = unique(mol_markers$gene)[1:20]

opc_markers = markers %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "OPC"), "OPC", cluster)) %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "NFOL"), "NFOL", cluster)) %>% 
  dplyr::filter(str_detect(cluster, "OPC")) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(100, wt = avg_log2FC)
opc_markers = unique(opc_markers$gene)[1:20]

# Seurat objects
sorted = readRDS(file = "../results/Seurat/final/sorted_brain/res0.1/outputs/Seurat_object.Rds")
g4_counts = as.matrix(sorted@assays$GA@data)
peaks = as.matrix(sorted@assays$peaks@data)

# highly predicted MOL cell ids
barcodes = pred_score$cell_id[which(pred_score$MOL > 0.75)]
barcodes_non_mol = rownames(sorted@meta.data)[which(!rownames(sorted@meta.data) %in% barcodes)]

diff_peaks_MOL_vs_nonMOL = FindMarkers(
  sorted,
  ident.1 = barcodes,
  ident.2 = barcodes_non_mol,
  only.pos = TRUE,
  assay = "peaks",
  logfc.threshold = 0,
  test.use = "LR",
  latent.vars = "peak_region_fragments"
)

# compute G4 fold changes
# condition: predicted MOL cells
# background: non-oligodendrocyte cells (sorted)
peaks = peaks[,barcodes]
peaks = tibble(peak = rownames(peaks))
peaks = peaks %>% separate(peak, sep = "-", into = c("V1", "V2", "V3")) %>% 
  mutate(V2 = as.numeric(V2), V3 = as.numeric(V3))
peaks$type = "MOL"
peaks = GRanges(
  seqnames = peaks$V1,
  ranges = IRanges(
    start = peaks$V2,
    end = peaks$V3,
    names = peaks$type,
  )
)

mol = get_unique_ws(bw = "../results/Seurat/callpeaks_GFPsorted/high_pred_MOL.bw", 
                    bw_backgr = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/1.bam_RPGC.bigwig",
                    subset = peaks)

# G4 signals around genes listed as MOL markers (literature or marker analysis)
mol_filt = mol %>% dplyr::filter(gene_symbol %in% all_mol_markers | gene_symbol %in% literature_markers)

# all sorted peakset
sorted_peaks = fread("../data/CellRanger/GFP_sorted/peaks.bed")
colnames(sorted_peaks) = c("V1", "V2", "V3")
# annotation
sorted_peaks = mm10_annotation(regions = sorted_peaks, seqname_col = "V1", start_col = "V2", end_col = "V3", 
                                 feature_1 = NULL)

# genes with the highest G4 signal over promoters
genes = mol_filt %>% dplyr::filter(abs(distanceToTSS) < 3000) %>% 
  dplyr::filter(fold_change > 10) %>% 
  pull(gene_symbol)

## topGO analysis (source: https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf)
# create input for GO analysis
input = rep(0, length(unique(sorted_peaks$SYMBOL)))
names(input) = unique(sorted_peaks$SYMBOL)
input[names(input) %in% genes] = 1

GOdata <- new(
  "topGOdata",
  ontology = "BP",
  # use biological process ontology
  allGenes = input,
  geneSelectionFun = function(x)
    (x == 1),
  annot = annFUN.org,
  mapping = "org.Mm.eg.db",
  ID = "symbol"
)

# Fisher test
resultFisher <-
  runTest(GOdata, algorithm = "elim", statistic = "fisher")
out <-
  GenTable(GOdata,
           Fisher = resultFisher,
           topNodes = 20,
           numChar = 60)

out = out %>% dplyr::select(Term, Fisher)

## expression heatmaps
# Marques et al. oligo scRNA data
rna = read.table(
  "../data/GSE75330/GSE75330_Marques_et_al_mol_counts2.tab",
  stringsAsFactors = FALSE,
  header = FALSE
)
annot = readRDS(file = "../data/GSE75330/Marques2016annotation.rds")

colnames(rna) = gsub("-", "_", sub('-', '', gsub("C1-", "", rna[1,])))
rna = rna[-1,]
genes_rna = rna[, 1]
rna = rna[,-1]
rna = as.matrix(rna)
rna = apply(rna, 2, as.numeric)
rownames(rna) = genes_rna

# create Seurat object and log normalize
rna = CreateSeuratObject(counts = rna,
                         meta.data = annot,
                         assay = 'RNA')
rna = NormalizeData(rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
norm = rna[["RNA"]]@data # normalized expression levels

# retrieve the expression levels of genes with highest G4 signal (100-fold higher than background!)
top_mol_genes = mol_filt %>% dplyr::filter(abs(distanceToTSS) < 3000) %>% 
  dplyr::filter(fold_change > 10) %>% arrange(desc(fold_change)) %>% 
  top_n(100, wt = fold_change) %>% pull(gene_symbol)

existing_gene_symbols = character()
for (gene in top_mol_genes) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}
#existing_gene_symbols = existing_gene_symbols[1:10]

# simplify oligo subtypes
new_ids = as.character(rna@meta.data$cell_class)
new_ids[new_ids == 'NFOL2'] = 'NFOL'
new_ids[new_ids == 'NFOL1'] = 'NFOL'
new_ids[new_ids == 'MFOL2'] = 'MFOL'
new_ids[new_ids == 'MFOL1'] = 'MFOL'
new_ids[new_ids == 'MOL1'] = 'MOL'
new_ids[new_ids == 'MOL2'] = 'MOL'
new_ids[new_ids == 'MOL3'] = 'MOL'
new_ids[new_ids == 'MOL4'] = 'MOL'
new_ids[new_ids == 'MOL5'] = 'MOL'
new_ids[new_ids == 'MOL6'] = 'MOL'
new_ids[new_ids == 'PPR'] = 'VLMC'
rna@meta.data$merged_cell_class = new_ids

cell_types = as.character(unique(rna@meta.data$merged_cell_class))

# helper function - collect aggregated mean expression by gene
collect = function(x) {
  m = as.matrix(norm[which(rownames(norm) %in% existing_gene_symbols), rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
  t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
  return(t)
}
expr = lapply(cell_types, collect)
expr = bind_rows(expr)
mat = pivot_wider(expr, id_cols = gene_symbol, names_from = type, values_from = means) 

# visualize by heatmap - genes w/ top MOL specific G4 promoters
y_order = mat %>% mutate(mean = rowMeans(dplyr::select(mat, -gene_symbol), na.rm = TRUE)) %>% arrange(mean) %>%
  pull(gene_symbol)
y_order = factor(expr$gene_symbol, levels = y_order)

x_order = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
x_order = colMeans(x_order)
x_order = names(x_order[order(x_order, decreasing = TRUE)])
x_order = factor(expr$type, levels = x_order)

hm_expr_G4s = ggplot(expr, aes(x = x_order, y = y_order, fill = means)) +
  geom_tile(color = "black",
            lwd = 1.0,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = max(expr$means) / 2,
    limits = c(0, 5)
  ) +
  xlab(label = "oligodendrocyte subtype") +
  ylab(label = "gene with MOL specific promoter G4") +
  labs(fill = "log norm", title = "") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(hm_expr_G4s)

ggsave(
  glue("{result_folder}top_MOL_G4_prom_signals-norm_expr_hm.pdf"),
  plot = hm_expr_G4s,
  width = 8,
  height = 7,
  device = "pdf"
)

norm_mol = norm[existing_gene_symbols, 
                 rownames(rna@meta.data[which(rna@meta.data$merged_cell_class == "MOL"),])]
norm_mol = as_tibble(norm_mol)
norm_mol = norm_mol %>% mutate(gene = existing_gene_symbols) %>% 
  dplyr::select(gene, everything())
norm_mol = pivot_longer(norm_mol, cols =
                          colnames(norm_mol)[2]:colnames(norm_mol)[ncol(norm_mol)],
                        names_to = "barcode", values_to = "expr")
norm_mol = norm_mol %>% dplyr::filter(expr > 0) %>% 
  mutate(prediction = "MOL", med = round(median(expr), 2))

norm_non_mol = norm[existing_gene_symbols, 
                    rownames(rna@meta.data[which(rna@meta.data$merged_cell_class != "MOL"),])]
norm_non_mol = as_tibble(norm_non_mol)
norm_non_mol = norm_non_mol %>% mutate(gene = existing_gene_symbols) %>% 
  dplyr::select(gene, everything())
norm_non_mol = pivot_longer(norm_non_mol, cols =
                          colnames(norm_non_mol)[2]:colnames(norm_non_mol)[ncol(norm_non_mol)],
                        names_to = "barcode", values_to = "expr")
norm_non_mol = norm_non_mol %>% dplyr::filter(expr > 0) %>% 
  mutate(prediction = "non-MOL", med = round(median(expr), 2))
norm_expr_mol_markers = bind_rows(list(norm_mol, norm_non_mol))

expr_violin = ggplot(norm_expr_mol_markers, aes(x = prediction, y = expr)) +
  geom_violin(color = "black", fill = "#4682b4") +
  ylim(0, 10) +
  labs(
    title = "genes with MOL specific promoter G4",
    x = "scRNA-Seq cluster",
    y = "norm. expression"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 7.5, label.x = 1.25, size = 6) +
  geom_text(data = norm_expr_mol_markers, aes(x = prediction, y = med, label = med),
            size = 5, vjust = -0.5, color = "white")
expr_violin

ggsave(
  glue("{result_folder}top_MOL_G4_prom_signals-norm_expr_violin.pdf"),
  plot = expr_violin,
  width = 8,
  height = 7,
  device = "pdf"
)

# retrieve the expression levels of MOL marker genes of FindMarker analysis
existing_gene_symbols = character()
for (gene in mol_markers) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}
existing_gene_symbols = existing_gene_symbols[1:10]

# helper function - collect aggregated mean expression by gene
collect = function(x) {
  m = as.matrix(norm[which(rownames(norm) %in% existing_gene_symbols), 
                     rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
  t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
  return(t)
}

expr = lapply(cell_types, collect)
expr = bind_rows(expr)
mat = pivot_wider(expr, id_cols = gene_symbol, names_from = type, values_from = means)

# visualize by heatmap - MOL markers
y_order = mat %>% mutate(mean = rowMeans(dplyr::select(mat, -gene_symbol), na.rm = TRUE)) %>% arrange(mean) %>%
  pull(gene_symbol)
y_order = factor(expr$gene_symbol, levels = y_order)

x_order = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
x_order = colMeans(x_order)
x_order = names(x_order[order(x_order, decreasing = TRUE)])
x_order = factor(expr$type, levels = x_order)

hm_expr_mol_markers = ggplot(expr, aes(x = x_order, y = y_order, fill = means)) +
  geom_tile(color = "black",
            lwd = 1.0,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = max(expr$means) / 2,
    limits = c(0, 6)
  ) +
  xlab(label = "oligodendrocyte subtype") +
  ylab(label = "top MOL expression markers") +
  labs(fill = "log norm. expr.", title = "") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(hm_expr_mol_markers)

ggsave(
  glue("{result_folder}top_MOL_markers-norm_expr_hm.pdf"),
  plot = hm_expr_mol_markers,
  width = 8,
  height = 7,
  device = "pdf"
)

# retrieve the GA scores of highly predicted MOL cells
g4_high_pred_counts = rowMeans(g4_counts[existing_gene_symbols, barcodes])
g4_high_pred_counts = tibble(gene_symbol = names(g4_high_pred_counts), G4_count = g4_high_pred_counts, predictions = "MOL")
g4_nonhigh_pred_counts = rowMeans(g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)])
g4_nonhigh_pred_counts = tibble(gene_symbol = names(g4_nonhigh_pred_counts), G4_count = g4_nonhigh_pred_counts, predictions = "non-MOL")
g4_counts_MOL_markers = bind_rows(list(g4_nonhigh_pred_counts, g4_high_pred_counts))

y_order2 = factor(g4_counts_MOL_markers$gene_symbol, levels = levels(y_order))
x_order2 = factor(g4_counts_MOL_markers$predictions, levels = c("MOL", "non-MOL"))

prediction_hm = ggplot(g4_counts_MOL_markers, aes(x = x_order2, y = y_order2, fill = G4_count)) +
  geom_tile(color = "black",
            lwd = 1.0,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = max(g4_counts_MOL_markers$G4_count) / 2,
    limits = c(0, max(g4_counts_MOL_markers$G4_count))
  ) +
  xlab(label = "prediction") +
  ylab(label = "top MOL expression markers") +
  labs(fill = "log norm. G4", title = "") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 45,
      hjust = 1,
      vjust = 1.2
    ),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(prediction_hm)

ggsave(
  glue("{result_folder}top_MOL_markers-G4_count_hm.pdf"),
  plot = prediction_hm,
  width = 6,
  height = 7,
  device = "pdf"
)

# retrieve the GA scores of highly predicted MOL cells (top 100)
existing_gene_symbols = character()
for (gene in all_mol_markers) {
  if (gene %in% rownames(g4_counts)) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

# make input for violins (1)
x = g4_counts[existing_gene_symbols, barcodes]
rows = rownames(x)
x_long = pivot_longer(
  as_tibble(x),
  cols =
    colnames(x)[1]:colnames(x)[ncol(x)],
  names_to = "barcode",
  values_to = "G4_count"
)
x_long = x_long %>% 
  dplyr::filter(G4_count > 0) %>% 
  mutate(predictions = "MOL", med = round(median(G4_count), 3)) 

y = g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)]
rows = rownames(y)
y_long =
  pivot_longer(
    as_tibble(y),
    cols = colnames(y)[1]:colnames(y)[ncol(y)],
    names_to = "barcode",
    values_to = "G4_count"
  )
y_long = y_long %>% 
  dplyr::filter(G4_count > 0) %>% 
  mutate(predictions = "non-MOL", med = round(median(G4_count), 3)) 

x_y = bind_rows(list(x_long, y_long))

prediction_boxplot = ggplot(x_y, aes(x = predictions, y = G4_count)) +
geom_boxplot(color = "black", fill = "#4682b4") +
  ylim(0, 2.5) +
  labs(
    title = "top MOL expression markers",
    x = "prediction",
    y = "log norm. G4"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 2, label.x = 1.25, size = 6) +
  geom_text(data = x_y, aes(x = predictions, y = med, label = med),
            size = 5, vjust = -0.5, color = "white")
prediction_boxplot

ggsave(
  glue("{result_folder}MOL-nonMOL_G4_count_boxplot.pdf"),
  plot = prediction_boxplot,
  width = 5,
  height = 5,
  device = "pdf"
)

# anchor score (coming from scRNA-Seq integration) distribution of highly predicted MOL cells 
anchors = anchors %>% mutate(predictions = ifelse(anchor2_barcode %in% barcodes, "MOL", "non-MOL"))
anchor_score_boxplot = ggplot(anchors, aes(x = predictions, y = score)) +
  geom_boxplot(color = "black", fill = "#4682b4") +
  ylim(0, 1) +
  labs(title = "",
       x = "prediction",
       y = "anchor score") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 1.0,
                     label.x = 1.25,
                     size = 4)
anchor_score_boxplot

ggsave(
  glue("{result_folder}MOL-nonMOL_anchorscore_boxplot.pdf"),
  plot = anchor_score_boxplot,
  width = 5,
  height = 4,
  device = "pdf"
)

# retrieve the GA scores of highly predicted MOL cells (at Zeisel & Marques oligodendrocyte marker genes)
existing_gene_symbols = character()
for (gene in literature_markers) {
  if (gene %in% rownames(g4_counts)) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}

g4_lit_counts = rowMeans(g4_counts[existing_gene_symbols, barcodes])
g4_lit_counts = tibble(gene_symbol = names(g4_lit_counts), G4_count = g4_lit_counts, predictions = "MOL")
g4_nonhigh_pred_counts = rowMeans(g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)])
g4_nonhigh_pred_counts = tibble(gene_symbol = names(g4_nonhigh_pred_counts), 
                                G4_count = g4_nonhigh_pred_counts, predictions = "non-MOL")
g4_counts_MOL_markers = bind_rows(list(g4_nonhigh_pred_counts, g4_lit_counts))
y_order = g4_counts_MOL_markers %>% arrange(G4_count) %>% 
  pull(gene_symbol) %>% unique 

y_order2 = factor(g4_counts_MOL_markers$gene_symbol, levels = y_order)
x_order2 = factor(g4_counts_MOL_markers$predictions, levels = c("MOL", "non-MOL"))

prediction_hm2 = ggplot(g4_counts_MOL_markers, aes(x = x_order2, y = y_order2, fill = G4_count)) +
  geom_tile(color = "black",
            lwd = 1.0,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = max(g4_counts_MOL_markers$G4_count) / 2,
    limits = c(0, max(g4_counts_MOL_markers$G4_count))
  ) +
  xlab(label = "prediction") +
  ylab(label = "top MOL expression markers") +
  labs(fill = "log norm. G4", title = "Zeisel et al. & Marques et al. oligodendr. markers") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 45,
      hjust = 1,
      vjust = 1.2
    ),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(prediction_hm2)

ggsave(
  glue("{result_folder}top_MOL_markers-Zeisel.Marques-G4_count_hm.pdf"),
  plot = prediction_hm2,
  width = 6,
  height = 9,
  device = "pdf"
)

# make input for violins (2)
x = g4_counts[existing_gene_symbols, barcodes]
rows = rownames(x)
x_long = pivot_longer(
  as_tibble(x),
  cols =
    colnames(x)[1]:colnames(x)[ncol(x)],
  names_to = "barcode",
  values_to = "G4_count"
)
x_long = x_long %>% mutate(predictions = "MOL")

y = g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)]
rows = rownames(y)
y_long =
  pivot_longer(
    as_tibble(y),
    cols = colnames(y)[1]:colnames(y)[ncol(y)],
    names_to = "barcode",
    values_to = "G4_count"
  )
y_long = y_long %>% mutate(predictions = "non-MOL")

x_y = bind_rows(list(x_long, y_long))
x_y = x_y %>% dplyr::filter(G4_count > 0)

prediction_violin2 = ggplot(x_y, aes(x = predictions, y = G4_count)) +
  geom_violin(color = "black", fill = "#4682b4") +
  ylim(0, 4) +
  labs(
    title = "markers from Zeisel et al and Marques et al",
    x = "prediction",
    y = "log norm. G4"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 3, label.x = 1.25, size = 6)
prediction_violin2

ggsave(
  glue("{result_folder}MOL-nonMOL_Zeisel_markers.pdf"),
  plot = prediction_violin2,
  width = 5,
  height = 4,
  device = "pdf"
)


## retrieve the GA scores of highly predicted MOL cells (at top MOL specific G4 genes)
existing_gene_symbols = character()
for (gene in top_mol_genes) {
  if (gene %in% rownames(g4_counts)) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}
existing_gene_symbols = existing_gene_symbols[1:10]
existing_gene_symbols = existing_gene_symbols[!is.na(existing_gene_symbols)]

g4_top_mol_counts = rowMeans(g4_counts[existing_gene_symbols, barcodes])
g4_top_mol_counts = tibble(gene_symbol = names(g4_top_mol_counts), G4_count = g4_top_mol_counts, predictions = "MOL")
g4_nonhigh_pred_counts = rowMeans(g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)])
g4_nonhigh_pred_counts = tibble(gene_symbol = names(g4_nonhigh_pred_counts), 
                                G4_count = g4_nonhigh_pred_counts, predictions = "non-MOL")
g4_counts_top_mol_genes = bind_rows(list(g4_nonhigh_pred_counts, g4_top_mol_counts))
y_order = g4_counts_top_mol_genes %>% arrange(G4_count) %>% 
  pull(gene_symbol) %>% unique 

y_order3 = factor(g4_counts_top_mol_genes$gene_symbol, levels = y_order)
x_order3 = factor(g4_counts_top_mol_genes$predictions, levels = c("MOL", "non-MOL"))

prediction_hm3 = ggplot(g4_counts_top_mol_genes, aes(x = x_order3, y = y_order3, fill = G4_count)) +
  geom_tile(color = "black",
            lwd = 1.0,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = max(g4_counts_top_mol_genes$G4_count) / 2,
    limits = c(0, max(g4_counts_top_mol_genes$G4_count))
  ) +
  xlab(label = "prediction") +
  ylab(label = "gene with MOL specific promoter G4") +
  labs(fill = "log norm. G4", title = "") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 45,
      hjust = 1,
      vjust = 1.2
    ),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(prediction_hm3)

ggsave(
  glue("{result_folder}top_MOL_G4_prom_signals-G4_count_hm.pdf"),
  plot = prediction_hm3,
  width = 6,
  height = 9,
  device = "pdf"
)

# make input for violins (3)
x = g4_counts[existing_gene_symbols, barcodes]
rows = rownames(x)
x_long = pivot_longer(
  as_tibble(x),
  cols =
    colnames(x)[1]:colnames(x)[ncol(x)],
  names_to = "barcode",
  values_to = "G4_count"
)
x_long = x_long %>% mutate(predictions = "MOL")

y = g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)]
rows = rownames(y)
y_long =
  pivot_longer(
    as_tibble(y),
    cols = colnames(y)[1]:colnames(y)[ncol(y)],
    names_to = "barcode",
    values_to = "G4_count"
  )
y_long = y_long %>% mutate(predictions = "non-MOL")

x_y = bind_rows(list(x_long, y_long))
x_y = x_y %>% dplyr::filter(G4_count > 0)

prediction_violin3 = ggplot(x_y, aes(x = predictions, y = G4_count)) +
  geom_violin(color = "black", fill = "#4682b4") +
  ylim(0, 3) +
  labs(
    title = "gene with MOL specific promoter G4",
    x = "prediction",
    y = "log norm. G4"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 2.5, label.x = 1.25, size = 6)
prediction_violin3

ggsave(
  glue("{result_folder}MOL-nonMOL_MOLspec_promG4s.pdf"),
  plot = prediction_violin3,
  width = 5,
  height = 4,
  device = "pdf"
)

# other markers
# helper function - collect aggregated mean expression by gene
collect = function(x) {
  m = as.matrix(norm[which(rownames(norm) %in% c("Olfml1", "Spock1")), 
                     rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
  t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
  return(t)
}
misc = lapply(cell_types, collect)
misc = bind_rows(misc)
mat = pivot_wider(misc, id_cols = gene_symbol, names_from = type, values_from = means)

# visualize by heatmap - MOL markers
y_order3 = mat %>% mutate(mean = rowMeans(dplyr::select(mat, -gene_symbol), na.rm = TRUE)) %>% arrange(mean) %>%
  pull(gene_symbol)
y_order3 = factor(misc$gene_symbol, levels = y_order3)

x_order3 = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
x_order3 = colMeans(x_order3)
x_order3 = names(x_order3[order(x_order3, decreasing = TRUE)])
x_order3 = factor(misc$type, levels = x_order3)

misc_hm = ggplot(misc, aes(x = x_order3, y = y_order3, fill = means)) +
  geom_tile(color = "black",
            lwd = 1.0,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = max(misc$means) / 2,
    limits = c(0, 0.5)
  ) +
  xlab(label = "oligodendrocyte subtype") +
  ylab(label = "") +
  labs(fill = "log norm. expr.", title = "") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 14,
      angle = 0,
      hjust = 0.5,
      vjust = 1.2
    ),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(misc_hm)

ggsave(
  glue("{result_folder}top_MOL_markers-G4_count_hm-example.pdf"),
  plot = misc_hm,
  width = 6,
  height = 4,
  device = "pdf"
)

# retrieve genes showing highest GA score in predicted MOL cells
top_G4_count_of_MOLs = rowMeans(g4_counts[, barcodes])
top_G4_count_of_MOLs = top_G4_count_of_MOLs[order(top_G4_count_of_MOLs, decreasing = TRUE)][1:100]
top_G4_count_of_MOLs = names(top_G4_count_of_MOLs)

existing_gene_symbols = character()
for (gene in top_G4_count_of_MOLs) {
  if (gene %in% norm@Dimnames[[1]]) {
    existing_gene_symbols = c(existing_gene_symbols, gene)
  }
}
existing_gene_symbols = existing_gene_symbols[1:10]

collect = function(x) {
  m = as.matrix(norm[which(rownames(norm) %in% existing_gene_symbols), rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
  t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
  return(t)
}
expr = lapply(cell_types, collect)
expr = bind_rows(expr)
mat = pivot_wider(expr, id_cols = gene_symbol, names_from = type, values_from = means)

# visualize by heatmap - genes w/ top MOL specific G4 promoters
y_order = mat %>% mutate(mean = rowMeans(dplyr::select(mat, -gene_symbol), na.rm = TRUE)) %>% arrange(mean) %>%
  pull(gene_symbol)
y_order = factor(expr$gene_symbol, levels = y_order)

x_order = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
x_order = colMeans(x_order)
x_order = names(x_order[order(x_order, decreasing = TRUE)])
x_order = factor(expr$type, levels = x_order)

# visualize by heatmap - genes with highest G4 count in predicted MOLs
hm_top_G4_count_of_MOLs = ggplot(expr, aes(x = x_order, y = y_order, fill = means)) +
  geom_tile(color = "black",
            lwd = 1.0,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#a6bddb",
    mid = "#FFFFCC",
    high = "#fc9272",
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  xlab(label = "oligodendrocyte subtype") +
  ylab(label = "genes with highest G4 count in MOLs") +
  labs(fill = "log norm", title = "") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(color = "black", size = 20),
    axis.title = element_text(size = 20)
  ) +
  coord_fixed()
print(hm_top_G4_count_of_MOLs)

ggsave(
  glue("{result_folder}highest_G4_in_MOLs_hm.pdf"),
  plot = hm_top_G4_count_of_MOLs,
  width = 12,
  height = 10,
  device = "pdf"
)

# genome browser examples
# Seurat objects
sorted = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")

meta = sorted@meta.data
meta = meta %>% mutate(MOL_status = ifelse(rownames(meta) %in% barcodes, "predicted MOL", "predicted non-MOL"))
sorted@meta.data = meta

# for(gene in genes) {
#   plot = CoveragePlot(
#     object = sorted,
#     region = gene,
#     annotation = TRUE,
#     show.bulk = TRUE,
#     ymax = 5,
#     group.by = "MOL_status",
#     peaks = TRUE
#   ) + scale_fill_brewer(type = "seq", palette = "Set3")
#   print(plot)
# }

tmem88b = CoveragePlot(
  object = sorted,
  region = "Tmem88b",
  annotation = TRUE,
  show.bulk = TRUE,
  ymax = 10,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")

ndrg1 = CoveragePlot(
  object = sorted,
  region = "Ndrg1",
  annotation = TRUE,
  show.bulk = TRUE,
  ymax = 10,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")

serpinb1a = CoveragePlot(
  object = sorted,
  region = "Serpinb1a",
  annotation = TRUE,
  show.bulk = TRUE,
  ymax = 10,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")

oxr1 = CoveragePlot(
  object = sorted,
  region = "Oxr1",
  annotation = TRUE,
  show.bulk = TRUE,
  ymax = 5,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3") 

anln = CoveragePlot(
  object = sorted,
  region = "Anln",
  annotation = TRUE,
  show.bulk = TRUE,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")

sepp1 = CoveragePlot(
  object = sorted,
  region = "Sepp1",
  annotation = TRUE,
  show.bulk = TRUE,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")

examples = ggarrange(ndrg1, anln, serpinb1a, ncol = 1, nrow = 3)

ggsave(
  glue("{result_folder}MOL_markers-G4_score_examples.pdf"),
  plot = examples,
  width = 6,
  height = 9,
  device = "pdf"
)

olfml1 = CoveragePlot(
  object = sorted,
  region = "Olfml1",
  annotation = TRUE,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")

ggsave(
  glue("{result_folder}MOL_markers-olfml1.pdf"),
  plot = olfml1,
  width = 6,
  height = 6,
  device = "pdf"
)

cspg5 = CoveragePlot(
  object = sorted,
  region = "Cspg5",
  annotation = TRUE,
  group.by = "seurat_clusters",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")

zfp407 = CoveragePlot(
  object = sorted,
  region = "Zfp407",
  annotation = TRUE,
  show.bulk = TRUE,
  group.by = "MOL_status",
  peaks = TRUE
) + scale_fill_brewer(type = "seq", palette = "Set3")
