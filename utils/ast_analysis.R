if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse",
               "Seurat",
               "Signac",
               "wigglescout",
               "data.table",
               "ggplot2",
               "ggpubr",
               "glue",
               "matrixStats",
               "ComplexHeatmap",
               "circlize",
               "topGO",
               "outliers",
               "enrichR",
               "ggrepel",
               "ggrastr"
)

# path to result folder
result_folder = "../results/Seurat/"

# helper functions
source("C:/Szabolcs/Karolinska/Data/scripts/get_unique_peaks_ws.R")
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# MOL markers
markers = read_tsv("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/scRNA-Seq-FindAllMarkers_output.tsv")
markers = markers[markers$p_val < 0.05 &
                        markers$avg_log2FC > 0.5, ]
markers = markers %>% 
  mutate(cluster = str_replace_all(cluster, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(cluster = str_replace_all(cluster, pattern = "Oligodendrocytes", replacement = "MOL"))

# scBridge outputs
rel = fread("../results/scBridge/output/scbridge_reliability.csv")
pred = fread("../results/scBridge/output/scbridge_predictions.csv", header = TRUE)
pred = pred %>% 
  mutate(Prediction = str_replace_all(Prediction, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(Prediction = str_replace_all(Prediction, pattern = "Oligodendrocytes", replacement = "MOL"))

zeisel_markers = read_tsv("../data/Zeisel_et_al_2015/Zeisel_et_al_markers.txt", col_names = FALSE)
marques_markers = read_tsv("../data/GSE75330/marker_genes.txt", col_names = FALSE)
zeisel_markers = zeisel_markers %>% dplyr::filter(X2 == "Oligodendrocyte")
marques_markers = marques_markers %>% dplyr::filter(X2 == "MOL")
literature_markers = unique(c(marques_markers$X1, zeisel_markers$X1))

# anchor matrix of GFP+ - scRNA integration
anchors = read_tsv("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/anchor_matrix.tsv")

# prediction scores
pred_score = readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/g4_cell_label_preds.Rds")

pred_score = t(pred_score@data)
ids = rownames(pred_score)
pred_score = as_tibble(pred_score)
colnames(pred_score) = c("OEC","AST","MOL","Pericytes","VEC","VLMC","COP-NFOL","OPC","max")
pred_score = pred_score %>% 
  dplyr::select(-max) %>% mutate(cell_id = ids) %>% dplyr::filter(MOL > 0.50)
ast_ids = pred_score %>% pull(cell_id)

ast_markers = markers %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "AST"), "AST", cluster)) %>% 
  dplyr::filter(str_detect(cluster, "AST")) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(500, wt = avg_log2FC)
all_ast_markers = ast_markers$gene
ast_markers = unique(ast_markers$gene)[1:20]

opc_markers = markers %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "OPC"), "OPC", cluster)) %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "NFOL"), "NFOL", cluster)) %>% 
  dplyr::filter(str_detect(cluster, "OPC")) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(100, wt = avg_log2FC)
opc_markers = unique(opc_markers$gene)[1:20]

# Seurat objects
sorted = readRDS(file = "../results/Seurat/final/sorted_brain/res0.1/outputs/Seurat_object.Rds")
rna = readRDS(file = "../data/GSE163484/brain_object.Rds")
rna@meta.data = rna@meta.data %>% 
  mutate(cell_type = str_replace_all(cell_type, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(cell_type = str_replace_all(cell_type, pattern = "Oligodendrocytes", replacement = "MOL"))
norm = rna[['RNA']]@data

g4_counts = as.matrix(sorted@assays$GA@data)
peaks = as.matrix(sorted@assays$peaks@data)

# highly predicted MOL cell ids
barcodes_seurat = pred_score$cell_id[which(pred_score$AST > 0.50)]
barcodes_scbridge = pred %>% filter(Prediction == "AST") %>% pull(V1)
barcodes_scbridge_others = pred %>% filter(Prediction != "AST") %>% 
  filter(Prediction != "Novel (Most Unreliable)") %>% 
  pull(V1)

barcodes = unique(c(barcodes_seurat, barcodes_scbridge))
barcodes_df = tibble(AST = barcodes, type = "AST")
write_tsv(barcodes_df, glue("{result_folder}Pred_AST-barcodes.tsv"))

barcodes_non_ast = rownames(sorted@meta.data)[which(!rownames(sorted@meta.data) %in% barcodes)]
barcodes_non_ast_df = tibble(non_AST = barcodes_non_ast, type = "non_AST")
write_tsv(barcodes_non_ast_df, glue("{result_folder}Pred_non_AST-barcodes.tsv"))

# differential GA scores
# diff_AST_vs_nonAST_seurat = FindMarkers(
#   sorted,
#   ident.1 = barcodes,
#   ident.2 = barcodes_non_ast,
#   only.pos = FALSE,
#   assay = "GA",
#   logfc.threshold = 0,
#   test.use = "LR",
#   latent.vars = "peak_region_fragments"
# )

# write.table(diff_AST_vs_nonAST_seurat, glue("{result_folder}diff_AST_vs_nonAST_seurat.tsv"),
#             quote = FALSE)
diff_AST_vs_nonAST_seurat = read.table(glue("{result_folder}diff_AST_vs_nonAST_seurat.tsv"))

sign_AST_vs_nonAST_seurat = diff_AST_vs_nonAST_seurat %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames

# diff_AST_vs_nonAST_scBr = FindMarkers(
#   sorted,
#   ident.1 = barcodes_scbridge,
#   ident.2 = barcodes_scbridge_others,
#   only.pos = FALSE,
#   assay = "GA",
#   logfc.threshold = 0,
#   test.use = "LR",
#   latent.vars = "peak_region_fragments"
# )
# 
# write.table(diff_AST_vs_nonAST_scBr, glue("{result_folder}diff_AST_vs_nonAST_scBridge.tsv"),
#             quote = FALSE)

diff_AST_vs_nonAST_scBr = read.table(glue("{result_folder}diff_AST_vs_nonAST_scBridge.tsv"))

sign_AST_vs_nonAST_scBr = diff_AST_vs_nonAST_scBr %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames

# compute G4 fold changes
# condition: predicted MOL cells
# background: non-oligodendrocyte cells (sorted)
# peaks = peaks[,barcodes]
# peaks = tibble(peak = rownames(peaks))
# peaks = peaks %>% separate(peak, sep = "-", into = c("V1", "V2", "V3")) %>% 
#   mutate(V2 = as.numeric(V2), V3 = as.numeric(V3))
# peaks$type = "MOL"
# peaks = GRanges(
#   seqnames = peaks$V1,
#   ranges = IRanges(
#     start = peaks$V2,
#     end = peaks$V3,
#     names = peaks$type,
#   )
# )
# 
# mol = get_unique_ws(bw = "../results/Seurat/callpeaks_GFPsorted/high_pred_MOL.bw", 
#                     bw_backgr = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/1.bam_RPGC.bigwig",
#                     subset = peaks)
# 
# # G4 signals around genes listed as MOL markers (literature or marker analysis)
# mol_filt = mol %>% dplyr::filter(gene_symbol %in% all_mol_markers | gene_symbol %in% literature_markers)
# 
# # all sorted peakset
# sorted_peaks = fread("../data/CellRanger/GFP_sorted/peaks.bed")
# colnames(sorted_peaks) = c("V1", "V2", "V3")
# # annotation
# sorted_peaks = mm10_annotation(regions = sorted_peaks, seqname_col = "V1", start_col = "V2", end_col = "V3", 
#                                  feature_1 = NULL)
# 
# # genes with the highest G4 signal over promoters
# genes = mol_filt %>% dplyr::filter(abs(distanceToTSS) < 3000) %>% 
#   dplyr::filter(fold_change > 10) %>% 
#   pull(gene_symbol)
# 
# ## expression heatmaps
# # Marques et al. oligo scRNA data
# rna = read.table(
#   "../data/GSE75330/GSE75330_Marques_et_al_mol_counts2.tab",
#   stringsAsFactors = FALSE,
#   header = FALSE
# )
# annot = readRDS(file = "../data/GSE75330/Marques2016annotation.rds")
# 
# colnames(rna) = gsub("-", "_", sub('-', '', gsub("C1-", "", rna[1,])))
# rna = rna[-1,]
# genes_rna = rna[, 1]
# rna = rna[,-1]
# rna = as.matrix(rna)
# rna = apply(rna, 2, as.numeric)
# rownames(rna) = genes_rna

# create Seurat object and log normalize
# rna = CreateSeuratObject(counts = rna,
#                          meta.data = annot,
#                          assay = 'RNA')
# rna = NormalizeData(rna,
#                     normalization.method = "LogNormalize",
#                     scale.factor = 10000)
# norm = rna[["RNA"]]$data # normalized expression levels
# 
# # retrieve the expression levels of genes with highest G4 signal (100-fold higher than background!)
# top_mol_genes = mol_filt %>% dplyr::filter(abs(distanceToTSS) < 3000) %>% 
#   dplyr::filter(fold_change > 10) %>% arrange(desc(fold_change)) %>% 
#   top_n(100, wt = fold_change) %>% pull(gene_symbol)
# 
# existing_gene_symbols = character()
# for (gene in top_mol_genes) {
#   if (gene %in% norm@Dimnames[[1]]) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
# #existing_gene_symbols = existing_gene_symbols[1:10]
# 
# # simplify oligo subtypes
# new_ids = as.character(rna@meta.data$cell_class)
# new_ids[new_ids == 'NFOL2'] = 'NFOL'
# new_ids[new_ids == 'NFOL1'] = 'NFOL'
# new_ids[new_ids == 'MFOL2'] = 'MFOL'
# new_ids[new_ids == 'MFOL1'] = 'MFOL'
# new_ids[new_ids == 'MOL1'] = 'MOL'
# new_ids[new_ids == 'MOL2'] = 'MOL'
# new_ids[new_ids == 'MOL3'] = 'MOL'
# new_ids[new_ids == 'MOL4'] = 'MOL'
# new_ids[new_ids == 'MOL5'] = 'MOL'
# new_ids[new_ids == 'MOL6'] = 'MOL'
# new_ids[new_ids == 'PPR'] = 'VLMC'
# rna@meta.data$merged_cell_class = new_ids
# 
# cell_types = as.character(unique(rna@meta.data$merged_cell_class))
# 
# ## visualizations ##
# # helper function - collect aggregated mean expression by gene
# collect = function(x) {
#   m = as.matrix(norm[which(rownames(norm) %in% existing_gene_symbols), rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
#   t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
#   return(t)
# }
# expr = lapply(cell_types, collect)
# expr = bind_rows(expr)
# mat = pivot_wider(expr, id_cols = gene_symbol, names_from = type, values_from = means) 
# 
# # visualize by heatmap - genes w/ top MOL specific G4 promoters
# y_order = mat %>% mutate(mean = rowMeans(dplyr::select(mat, -gene_symbol), na.rm = TRUE)) %>% arrange(mean) %>%
#   pull(gene_symbol)
# y_order = factor(expr$gene_symbol, levels = y_order)
# 
# x_order = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
# x_order = colMeans(x_order)
# x_order = names(x_order[order(x_order, decreasing = TRUE)])
# x_order = factor(expr$type, levels = x_order)
# 
# hm_expr_G4s = ggplot(expr, aes(x = x_order, y = y_order, fill = means)) +
#   geom_tile(color = "black",
#             lwd = 1.0,
#             linetype = 1) +
#   scale_fill_gradient2(
#     low = "#a6bddb",
#     mid = "#FFFFCC",
#     high = "#fc9272",
#     midpoint = max(expr$means) / 2,
#     limits = c(0, 5)
#   ) +
#   xlab(label = "oligodendrocyte subtype") +
#   ylab(label = "gene with MOL specific promoter G4") +
#   labs(fill = "log norm", title = "") +
#   theme_classic() +
#   theme(
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(
#       color = "black",
#       size = 15,
#       angle = 90,
#       hjust = 1,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(color = "black", size = 20),
#     axis.title = element_text(size = 20)
#   ) +
#   coord_fixed()
# print(hm_expr_G4s)
# 
# ggsave(
#   glue("{result_folder}top_MOL_G4_prom_signals-norm_expr_hm.pdf"),
#   plot = hm_expr_G4s,
#   width = 8,
#   height = 7,
#   device = "pdf"
# )
# 
# norm_mol = norm[existing_gene_symbols, 
#                  rownames(rna@meta.data[which(rna@meta.data$merged_cell_class == "MOL"),])]
# norm_mol = as_tibble(norm_mol)
# norm_mol = norm_mol %>% mutate(gene = existing_gene_symbols) %>% 
#   dplyr::select(gene, everything())
# norm_mol = pivot_longer(norm_mol, cols =
#                           colnames(norm_mol)[2]:colnames(norm_mol)[ncol(norm_mol)],
#                         names_to = "barcode", values_to = "expr")
# norm_mol = norm_mol %>% dplyr::filter(expr > 0) %>% 
#   mutate(prediction = "MOL", med = round(median(expr), 2))
# 
# norm_non_mol = norm[existing_gene_symbols, 
#                     rownames(rna@meta.data[which(rna@meta.data$merged_cell_class != "MOL"),])]
# norm_non_mol = as_tibble(norm_non_mol)
# norm_non_mol = norm_non_mol %>% mutate(gene = existing_gene_symbols) %>% 
#   dplyr::select(gene, everything())
# norm_non_mol = pivot_longer(norm_non_mol, cols =
#                           colnames(norm_non_mol)[2]:colnames(norm_non_mol)[ncol(norm_non_mol)],
#                         names_to = "barcode", values_to = "expr")
# norm_non_mol = norm_non_mol %>% dplyr::filter(expr > 0) %>% 
#   mutate(prediction = "non-MOL", med = round(median(expr), 2))
# norm_expr_mol_markers = bind_rows(list(norm_mol, norm_non_mol))
# 
# expr_violin = ggplot(norm_expr_mol_markers, aes(x = prediction, y = expr)) +
#   geom_violin(color = "black", fill = "#4682b4") +
#   ylim(0, 10) +
#   labs(
#     title = "genes with MOL specific promoter G4",
#     x = "scRNA-Seq cluster",
#     y = "norm. expression"
#   ) +
#   theme_classic() +
#   theme(
#     text = element_text(size = 20),
#     plot.title = element_text(size = 15),
#     axis.text.x = element_text(size = 20, color = "black"),
#     axis.text.y = element_text(size = 20, color = "black")
#   ) +
#   stat_compare_means(label.y = 7.5, label.x = 1.25, size = 6) +
#   geom_text(data = norm_expr_mol_markers, aes(x = prediction, y = med, label = med),
#             size = 5, vjust = -0.5, color = "white")
# expr_violin
# 
# ggsave(
#   glue("{result_folder}top_MOL_G4_prom_signals-norm_expr_violin.pdf"),
#   plot = expr_violin,
#   width = 8,
#   height = 7,
#   device = "pdf"
# )
# 
# # retrieve the expression levels of MOL marker genes of FindMarker analysis
# existing_gene_symbols = character()
# for (gene in mol_markers) {
#   if (gene %in% norm@Dimnames[[1]]) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
# existing_gene_symbols = existing_gene_symbols[1:10]
# 
# # helper function - collect aggregated mean expression by gene
# collect = function(x) {
#   m = as.matrix(norm[which(rownames(norm) %in% existing_gene_symbols), 
#                      rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
#   t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
#   return(t)
# }
# 
# expr = lapply(cell_types, collect)
# expr = bind_rows(expr)
# mat = pivot_wider(expr, id_cols = gene_symbol, names_from = type, values_from = means)
# 
# # visualize by heatmap - MOL markers
# y_order = mat %>% mutate(mean = rowMeans(dplyr::select(mat, -gene_symbol), na.rm = TRUE)) %>% arrange(mean) %>%
#   pull(gene_symbol)
# y_order = factor(expr$gene_symbol, levels = y_order)
# 
# x_order = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
# x_order = colMeans(x_order)
# x_order = names(x_order[order(x_order, decreasing = TRUE)])
# x_order = factor(expr$type, levels = x_order)
# 
# hm_expr_mol_markers = ggplot(expr, aes(x = x_order, y = y_order, fill = means)) +
#   geom_tile(color = "black",
#             lwd = 1.0,
#             linetype = 1) +
#   scale_fill_gradient2(
#     low = "#a6bddb",
#     mid = "#FFFFCC",
#     high = "#fc9272",
#     midpoint = max(expr$means) / 2,
#     limits = c(0, 6)
#   ) +
#   xlab(label = "oligodendrocyte subtype") +
#   ylab(label = "top MOL expression markers") +
#   labs(fill = "log norm. expr.", title = "") +
#   theme_classic() +
#   theme(
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(
#       color = "black",
#       size = 15,
#       angle = 90,
#       hjust = 1,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(color = "black", size = 20),
#     axis.title = element_text(size = 20)
#   ) +
#   coord_fixed()
# print(hm_expr_mol_markers)
# 
# ggsave(
#   glue("{result_folder}top_MOL_markers-norm_expr_hm.pdf"),
#   plot = hm_expr_mol_markers,
#   width = 8,
#   height = 7,
#   device = "pdf"
# )
# 
# # retrieve the GA scores of highly predicted MOL cells
# g4_high_pred_counts = rowMeans(g4_counts[existing_gene_symbols, barcodes])
# g4_high_pred_counts = tibble(gene_symbol = names(g4_high_pred_counts), G4_count = g4_high_pred_counts, predictions = "MOL")
# g4_nonhigh_pred_counts = rowMeans(g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)])
# g4_nonhigh_pred_counts = tibble(gene_symbol = names(g4_nonhigh_pred_counts), G4_count = g4_nonhigh_pred_counts, predictions = "non-MOL")
# g4_counts_MOL_markers = bind_rows(list(g4_nonhigh_pred_counts, g4_high_pred_counts))
# 
# y_order2 = factor(g4_counts_MOL_markers$gene_symbol, levels = levels(y_order))
# x_order2 = factor(g4_counts_MOL_markers$predictions, levels = c("MOL", "non-MOL"))
# 
# prediction_hm = ggplot(g4_counts_MOL_markers, aes(x = x_order2, y = y_order2, fill = G4_count)) +
#   geom_tile(color = "black",
#             lwd = 1.0,
#             linetype = 1) +
#   scale_fill_gradient2(
#     low = "#a6bddb",
#     mid = "#FFFFCC",
#     high = "#fc9272",
#     midpoint = max(g4_counts_MOL_markers$G4_count) / 2,
#     limits = c(0, max(g4_counts_MOL_markers$G4_count))
#   ) +
#   xlab(label = "prediction") +
#   ylab(label = "top MOL expression markers") +
#   labs(fill = "log norm. G4", title = "") +
#   theme_classic() +
#   theme(
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(
#       color = "black",
#       size = 15,
#       angle = 45,
#       hjust = 1,
#       vjust = 1.2
#     ),
#     axis.text.y = element_text(color = "black", size = 20),
#     axis.title = element_text(size = 20)
#   ) +
#   coord_fixed()
# print(prediction_hm)
# 
# ggsave(
#   glue("{result_folder}top_MOL_markers-G4_count_hm.pdf"),
#   plot = prediction_hm,
#   width = 6,
#   height = 7,
#   device = "pdf"
# )
# 
# # retrieve the GA scores of highly predicted MOL cells (top 100)
# existing_gene_symbols = character()
# for (gene in all_mol_markers) {
#   if (gene %in% rownames(g4_counts)) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
# 
# # make input for violins (1)
# x = g4_counts[existing_gene_symbols, barcodes]
# rows = rownames(x)
# x_long = pivot_longer(
#   as_tibble(x),
#   cols =
#     colnames(x)[1]:colnames(x)[ncol(x)],
#   names_to = "barcode",
#   values_to = "G4_count"
# )
# x_long = x_long %>% 
#   dplyr::filter(G4_count > 0) %>% 
#   mutate(predictions = "MOL", med = round(median(G4_count), 3)) 
# 
# y = g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)]
# rows = rownames(y)
# y_long =
#   pivot_longer(
#     as_tibble(y),
#     cols = colnames(y)[1]:colnames(y)[ncol(y)],
#     names_to = "barcode",
#     values_to = "G4_count"
#   )
# y_long = y_long %>% 
#   dplyr::filter(G4_count > 0) %>% 
#   mutate(predictions = "non-MOL", med = round(median(G4_count), 3)) 
# 
# x_y = bind_rows(list(x_long, y_long))
# 
# prediction_boxplot = ggplot(x_y, aes(x = predictions, y = G4_count)) +
# geom_boxplot(color = "black", fill = "#4682b4") +
#   ylim(0, 2.5) +
#   labs(
#     title = "top MOL expression markers",
#     x = "prediction",
#     y = "log norm. G4"
#   ) +
#   theme_classic() +
#   theme(
#     text = element_text(size = 20),
#     plot.title = element_text(size = 15),
#     axis.text.x = element_text(size = 20, color = "black"),
#     axis.text.y = element_text(size = 20, color = "black")
#   ) +
#   stat_compare_means(label.y = 2, label.x = 1.25, size = 6) +
#   geom_text(data = x_y, aes(x = predictions, y = med, label = med),
#             size = 5, vjust = -0.5, color = "white")
# prediction_boxplot
# 
# ggsave(
#   glue("{result_folder}MOL-nonMOL_G4_count_boxplot.pdf"),
#   plot = prediction_boxplot,
#   width = 5,
#   height = 5,
#   device = "pdf"
# )
# 
# # anchor score (coming from scRNA-Seq integration) distribution of highly predicted MOL cells 
# anchors = anchors %>% mutate(predictions = ifelse(anchor2_barcode %in% barcodes, "MOL", "non-MOL"))
# anchor_score_boxplot = ggplot(anchors, aes(x = predictions, y = score)) +
#   geom_boxplot(color = "black", fill = "#4682b4") +
#   ylim(0, 1) +
#   labs(title = "",
#        x = "prediction",
#        y = "anchor score") +
#   theme_classic() +
#   theme(
#     text = element_text(size = 20),
#     plot.title = element_text(size = 10),
#     axis.text.x = element_text(size = 20, color = "black"),
#     axis.text.y = element_text(size = 20, color = "black")
#   ) +
#   stat_compare_means(label.y = 1.0,
#                      label.x = 1.25,
#                      size = 4)
# anchor_score_boxplot
# 
# ggsave(
#   glue("{result_folder}MOL-nonMOL_anchorscore_boxplot.pdf"),
#   plot = anchor_score_boxplot,
#   width = 5,
#   height = 4,
#   device = "pdf"
# )
# 
# # retrieve the GA scores of highly predicted MOL cells (at Zeisel & Marques oligodendrocyte marker genes)
# existing_gene_symbols = character()
# for (gene in literature_markers) {
#   if (gene %in% rownames(g4_counts)) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
# 
# g4_lit_counts = rowMeans(g4_counts[existing_gene_symbols, barcodes])
# g4_lit_counts = tibble(gene_symbol = names(g4_lit_counts), G4_count = g4_lit_counts, predictions = "MOL")
# g4_nonhigh_pred_counts = rowMeans(g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)])
# g4_nonhigh_pred_counts = tibble(gene_symbol = names(g4_nonhigh_pred_counts), 
#                                 G4_count = g4_nonhigh_pred_counts, predictions = "non-MOL")
# g4_counts_MOL_markers = bind_rows(list(g4_nonhigh_pred_counts, g4_lit_counts))
# y_order = g4_counts_MOL_markers %>% arrange(G4_count) %>% 
#   pull(gene_symbol) %>% unique 
# 
# y_order2 = factor(g4_counts_MOL_markers$gene_symbol, levels = y_order)
# x_order2 = factor(g4_counts_MOL_markers$predictions, levels = c("MOL", "non-MOL"))
# 
# prediction_hm2 = ggplot(g4_counts_MOL_markers, aes(x = x_order2, y = y_order2, fill = G4_count)) +
#   geom_tile(color = "black",
#             lwd = 1.0,
#             linetype = 1) +
#   scale_fill_gradient2(
#     low = "#a6bddb",
#     mid = "#FFFFCC",
#     high = "#fc9272",
#     midpoint = max(g4_counts_MOL_markers$G4_count) / 2,
#     limits = c(0, max(g4_counts_MOL_markers$G4_count))
#   ) +
#   xlab(label = "prediction") +
#   ylab(label = "top MOL expression markers") +
#   labs(fill = "log norm. G4", title = "Zeisel et al. & Marques et al. oligodendr. markers") +
#   theme_classic() +
#   theme(
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(
#       color = "black",
#       size = 15,
#       angle = 45,
#       hjust = 1,
#       vjust = 1.2
#     ),
#     plot.title = element_text(hjust = 0.5),
#     axis.text.y = element_text(color = "black", size = 12),
#     axis.title = element_text(size = 20)
#   ) +
#   coord_fixed()
# print(prediction_hm2)
# 
# ggsave(
#   glue("{result_folder}top_MOL_markers-Zeisel.Marques-G4_count_hm.pdf"),
#   plot = prediction_hm2,
#   width = 6,
#   height = 9,
#   device = "pdf"
# )
# 
# # make input for violins (2)
# x = g4_counts[existing_gene_symbols, barcodes]
# rows = rownames(x)
# x_long = pivot_longer(
#   as_tibble(x),
#   cols =
#     colnames(x)[1]:colnames(x)[ncol(x)],
#   names_to = "barcode",
#   values_to = "G4_count"
# )
# x_long = x_long %>% mutate(predictions = "MOL",  mean = round(mean(G4_count), 3))
# 
# y = g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)]
# rows = rownames(y)
# y_long =
#   pivot_longer(
#     as_tibble(y),
#     cols = colnames(y)[1]:colnames(y)[ncol(y)],
#     names_to = "barcode",
#     values_to = "G4_count"
#   )
# y_long = y_long %>% mutate(predictions = "non-MOL", mean = round(mean(G4_count), 3))
# 
# x_y = bind_rows(list(x_long, y_long))
# x_y = x_y %>% dplyr::filter(G4_count > 0)
# 
# prediction_violin2 = ggplot(x_y, aes(x = predictions, y = G4_count)) +
#   geom_violin(color = "black", fill = "#4682b4") +
#   ylim(0, 2) +
#   labs(
#     title = "markers from Zeisel et al and Marques et al",
#     x = "prediction",
#     y = "log norm. G4"
#   ) +
#   theme_classic() +
#   theme(
#     text = element_text(size = 20),
#     plot.title = element_text(size = 15),
#     axis.text.x = element_text(size = 20, color = "black"),
#     axis.text.y = element_text(size = 20, color = "black")
#   ) +
#   stat_compare_means(label.y = 1.75, label.x = 1.25, size = 6) +
#   geom_text(data = x_y, aes(x = predictions, y = 0.5, label = mean),
#             size = 5, vjust = -0.5, color = "white")
# prediction_violin2
# 
# ggsave(
#   glue("{result_folder}MOL-nonMOL_Zeisel_markers.pdf"),
#   plot = prediction_violin2,
#   width = 5,
#   height = 4,
#   device = "pdf"
# )
# 
# ## retrieve the GA scores of highly predicted MOL cells (at top MOL specific G4 genes)
# existing_gene_symbols = character()
# for (gene in top_mol_genes) {
#   if (gene %in% rownames(g4_counts)) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
# #existing_gene_symbols = existing_gene_symbols[1:10]
# existing_gene_symbols = existing_gene_symbols[!is.na(existing_gene_symbols)]
# 
# g4_top_mol_counts = rowMeans(g4_counts[existing_gene_symbols, barcodes])
# g4_top_mol_counts = tibble(gene_symbol = names(g4_top_mol_counts), G4_count = g4_top_mol_counts, predictions = "MOL")
# g4_nonhigh_pred_counts = rowMeans(g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)])
# g4_nonhigh_pred_counts = tibble(gene_symbol = names(g4_nonhigh_pred_counts), 
#                                 G4_count = g4_nonhigh_pred_counts, predictions = "non-MOL")
# g4_counts_top_mol_genes = bind_rows(list(g4_nonhigh_pred_counts, g4_top_mol_counts))
# y_order = g4_counts_top_mol_genes %>% arrange(G4_count) %>% 
#   pull(gene_symbol) %>% unique 
# 
# y_order3 = factor(g4_counts_top_mol_genes$gene_symbol, levels = y_order)
# x_order3 = factor(g4_counts_top_mol_genes$predictions, levels = c("MOL", "non-MOL"))
# 
# prediction_hm3 = ggplot(g4_counts_top_mol_genes, aes(x = x_order3, y = y_order3, fill = G4_count)) +
#   geom_tile(color = "black",
#             lwd = 1.0,
#             linetype = 1) +
#   scale_fill_gradient2(
#     low = "#a6bddb",
#     mid = "#FFFFCC",
#     high = "#fc9272",
#     midpoint = max(g4_counts_top_mol_genes$G4_count) / 2,
#     limits = c(0, max(g4_counts_top_mol_genes$G4_count))
#   ) +
#   xlab(label = "prediction") +
#   ylab(label = "gene with MOL specific promoter G4") +
#   labs(fill = "log norm. G4", title = "") +
#   theme_classic() +
#   theme(
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(
#       color = "black",
#       size = 15,
#       angle = 45,
#       hjust = 1,
#       vjust = 1.2
#     ),
#     plot.title = element_text(hjust = 0.5),
#     axis.text.y = element_text(color = "black", size = 20),
#     axis.title = element_text(size = 20)
#   ) +
#   coord_fixed()
# print(prediction_hm3)
# 
# ggsave(
#   glue("{result_folder}top_MOL_G4_prom_signals-G4_count_hm.pdf"),
#   plot = prediction_hm3,
#   width = 6,
#   height = 9,
#   device = "pdf"
# )
# 
# # make input for violins (3)
# x = g4_counts[existing_gene_symbols, barcodes]
# rows = rownames(x)
# x_long = pivot_longer(
#   as_tibble(x),
#   cols =
#     colnames(x)[1]:colnames(x)[ncol(x)],
#   names_to = "barcode",
#   values_to = "G4_count"
# )
# x_long = x_long %>% mutate(predictions = "MOL",  mean = round(mean(G4_count), 3))
# 
# y = g4_counts[existing_gene_symbols, which(colnames(g4_counts) != barcodes)]
# rows = rownames(y)
# y_long =
#   pivot_longer(
#     as_tibble(y),
#     cols = colnames(y)[1]:colnames(y)[ncol(y)],
#     names_to = "barcode",
#     values_to = "G4_count"
#   )
# y_long = y_long %>% mutate(predictions = "non-MOL", mean = round(mean(G4_count), 3))
# 
# x_y = bind_rows(list(x_long, y_long))
# x_y = x_y %>% dplyr::filter(G4_count > 0)
# 
# prediction_violin3 = ggplot(x_y, aes(x = predictions, y = G4_count)) +
#   geom_violin(color = "black", fill = "#4682b4") +
#   ylim(0, 3) +
#   labs(
#     title = "gene with MOL specific promoter G4",
#     x = "prediction",
#     y = "log norm. G4"
#   ) +
#   theme_classic() +
#   theme(
#     text = element_text(size = 20),
#     plot.title = element_text(size = 15),
#     axis.text.x = element_text(size = 20, color = "black"),
#     axis.text.y = element_text(size = 20, color = "black")
#   ) +
#   stat_compare_means(label.y = 2.5, label.x = 1.25, size = 6) +
#   geom_text(data = x_y, aes(x = predictions, y = 0.5, label = mean),
#             size = 5, vjust = -0.5, color = "white")
# prediction_violin3
# 
# ggsave(
#   glue("{result_folder}MOL-nonMOL_MOLspec_promG4s.pdf"),
#   plot = prediction_violin3,
#   width = 5,
#   height = 4,
#   device = "pdf"
# )
# 
# # other markers
# # helper function - collect aggregated mean expression by gene
# collect = function(x) {
#   m = as.matrix(norm[which(rownames(norm) %in% c("Olfml1", "Spock1")), 
#                      rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
#   t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
#   return(t)
# }
# misc = lapply(cell_types, collect)
# misc = bind_rows(misc)
# mat = pivot_wider(misc, id_cols = gene_symbol, names_from = type, values_from = means)
# 
# # visualize by heatmap - MOL markers
# y_order3 = mat %>% mutate(mean = rowMeans(dplyr::select(mat, -gene_symbol), na.rm = TRUE)) %>% arrange(mean) %>%
#   pull(gene_symbol)
# y_order3 = factor(misc$gene_symbol, levels = y_order3)
# 
# x_order3 = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
# x_order3 = colMeans(x_order3)
# x_order3 = names(x_order3[order(x_order3, decreasing = TRUE)])
# x_order3 = factor(misc$type, levels = x_order3)
# 
# misc_hm = ggplot(misc, aes(x = x_order3, y = y_order3, fill = means)) +
#   geom_tile(color = "black",
#             lwd = 1.0,
#             linetype = 1) +
#   scale_fill_gradient2(
#     low = "#a6bddb",
#     mid = "#FFFFCC",
#     high = "#fc9272",
#     midpoint = max(misc$means) / 2,
#     limits = c(0, 0.5)
#   ) +
#   xlab(label = "oligodendrocyte subtype") +
#   ylab(label = "") +
#   labs(fill = "log norm. expr.", title = "") +
#   theme_classic() +
#   theme(
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(
#       color = "black",
#       size = 14,
#       angle = 0,
#       hjust = 0.5,
#       vjust = 1.2
#     ),
#     axis.text.y = element_text(color = "black", size = 20),
#     axis.title = element_text(size = 20)
#   ) +
#   coord_fixed()
# print(misc_hm)
# 
# ggsave(
#   glue("{result_folder}top_MOL_markers-G4_count_hm-example.pdf"),
#   plot = misc_hm,
#   width = 6,
#   height = 4,
#   device = "pdf"
# )
# 
# # retrieve genes showing highest GA score in predicted MOL cells
# top_G4_count_of_MOLs = rowMeans(g4_counts[, barcodes])
# top_G4_count_of_MOLs = top_G4_count_of_MOLs[order(top_G4_count_of_MOLs, decreasing = TRUE)][1:100]
# top_G4_count_of_MOLs = names(top_G4_count_of_MOLs)
# 
# existing_gene_symbols = character()
# for (gene in top_G4_count_of_MOLs) {
#   if (gene %in% norm@Dimnames[[1]]) {
#     existing_gene_symbols = c(existing_gene_symbols, gene)
#   }
# }
# #existing_gene_symbols = existing_gene_symbols[1:10]
# 
# collect = function(x) {
#   m = as.matrix(norm[which(rownames(norm) %in% existing_gene_symbols), rna@meta.data$cellid[which(rna@meta.data$merged_cell_class == x)]])
#   t = tibble(means = rowMeans(m), gene_symbol = rownames(m), type = x)
#   return(t)
# }
# expr = lapply(cell_types, collect)
# expr = bind_rows(expr)
# mat = pivot_wider(expr, id_cols = gene_symbol, names_from = type, values_from = means)
# highest_vars = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
# rownames(highest_vars) = mat$gene_symbol
# highest_vars = mat$gene_symbol[order(rowVars(highest_vars), decreasing = TRUE)][1:20]
# expr = expr %>% dplyr::filter(gene_symbol %in% highest_vars)
# 
# # visualize by heatmap - genes w/ top MOL specific G4 promoters
# mat = mat %>% dplyr::filter(gene_symbol %in% highest_vars)
# y_order = mat %>% dplyr::filter(gene_symbol %in% highest_vars)
# y_order = factor(expr$gene_symbol, levels = rev(highest_vars))
# 
# x_order = mat %>% dplyr::select(-gene_symbol) %>% as.matrix
# x_order = colMeans(x_order)
# x_order = names(x_order[order(x_order, decreasing = TRUE)])
# x_order = factor(expr$type, levels = x_order)
# 
# # visualize by heatmap - genes with highest G4 count in predicted MOLs
# hm_top_G4_count_of_MOLs = ggplot(expr, aes(x = x_order, y = y_order, fill = means)) +
#   geom_tile(color = "black",
#             lwd = 1.0,
#             linetype = 1) +
#   scale_fill_gradient2(
#     low = "#a6bddb",
#     mid = "#FFFFCC",
#     high = "#fc9272",
#     midpoint = 2,
#     limits = c(0, 4)
#   ) +
#   xlab(label = "oligodendrocyte subtype") +
#   ylab(label = "genes with highest G4 count in MOLs") +
#   labs(fill = "log norm", title = "") +
#   theme_classic() +
#   theme(
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.x = element_text(
#       color = "black",
#       size = 15,
#       angle = 90,
#       hjust = 1,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(color = "black", size = 20),
#     axis.title = element_text(size = 20)
#   ) +
#   coord_fixed()
# print(hm_top_G4_count_of_MOLs)
# 
# ggsave(
#   glue("{result_folder}highest_G4_in_MOLs_hm.pdf"),
#   plot = hm_top_G4_count_of_MOLs,
#   width = 12,
#   height = 10,
#   device = "pdf"
# )

# genome browser examples
# Seurat object
sorted = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")

meta = sorted@meta.data

# meta = meta %>%
#   rownames_to_column("cell_id") %>%
#   left_join(., pred, by = c("cell_id" = "V1")) %>%
#   dplyr::rename(scBridge_prediction = Prediction) 
# rownames(meta) = meta$cell_id

meta = meta %>%
  mutate(AST_status = ifelse(rownames(meta) %in% barcodes, "predicted AST", "predicted non-AST")) %>%
  mutate(AST_status_scBridge = ifelse(rownames(meta) %in% barcodes_scbridge, "predicted AST", "predicted non-AST"))
sorted@meta.data = meta

genes = markers %>% dplyr::filter(cluster == "AST") %>% arrange(desc(avg_log2FC)) %>% top_n(6, wt = avg_log2FC) %>% pull(gene)
ctrl = markers %>% dplyr::filter(cluster == "OPC") %>% arrange(desc(avg_log2FC)) %>% top_n(6, wt = avg_log2FC) %>% pull(gene)

plots = list()
for(i in seq(1:length(genes))) {
  plots[[i]] = CoveragePlot(
    object = sorted,
    region = genes[i],
    annotation = TRUE,
    show.bulk = TRUE,
    ymax = 10,
    group.by = "AST_status",
    peaks = TRUE
  )
  print(plot)
}

examples = ggarrange(plotlist = plots, ncol = 2, nrow = 3)

ggsave(
  glue("{result_folder}scBr_predictions-browser_example.pdf"),
  plot = examples,
  width = 12,
  height = 12,
  device = "pdf"
)

plots = list()
for(i in seq(1:length(genes))) {
  mol = norm[genes[i], rownames(rna@meta.data[which(rna@meta.data$cell_type == "AST"),])]
  mol = tibble(type = "AST", expr = mol)
  nonmol = norm[genes[i], rownames(rna@meta.data[which(rna@meta.data$cell_type != "non-AST"),])]
  nonmol = tibble(type = "non-AST", expr = nonmol)
  both = rbind(mol, nonmol)
  plots[[i]] = ggplot(both, aes(x = type, y = expr, fill = type)) +
    geom_boxplot(color = "black") +
    ylim(0, 200) +
    labs(
      title = genes[i],
      x = "",
      y = "log norm. expr"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black")
    ) +
    stat_compare_means(label.y = 160, label.x = 1.15, size = 6)
}
  
bp_examples = ggarrange(plotlist = plots, ncol = 2, nrow = 3)

ggsave(
  glue("{result_folder}AST_markers-browser_example-expr_bps.pdf"),
  plot = bp_examples,
  width = 12,
  height = 12,
  device = "pdf"
)

plots = list()
for(i in seq(1:length(sign_AST_vs_nonAST_scBr[1:6]))) {
  plots[[i]] = CoveragePlot(
    object = sorted,
    region = sign_AST_vs_nonAST_scBr[i],
    annotation = TRUE,
    show.bulk = TRUE,
    ymax = 5,
    group.by = "AST_status_scBridge",
    peaks = TRUE
  )
  print(plots[[i]])
}

sign_scBridge_genes = ggarrange(plotlist = plots, ncol = 2, nrow = 3)

# ggsave(
#   glue("{result_folder}MOL_markers-browser_example-sign_up_scBr.png"),
#   plot = sign_scBridge_genes,
#   width = 12,
#   height = 12,
#   dpi = 1200
# )

ggsave(
  glue("{result_folder}AST_markers-browser_example-sign_up_scBr.pdf"),
  plot = sign_scBridge_genes,
  width = 12,
  height = 12,
  device = "pdf"
)

plots = list()
for(i in seq(1:length(sign_AST_vs_nonAST_scBr[1:6]))) {
  mol = norm[sign_AST_vs_nonAST_scBr[i], rownames(rna@meta.data[which(rna@meta.data$cell_type == "AST"),])]
  mol = tibble(type = "AST", expr = mol)
  nonmol = norm[sign_AST_vs_nonAST_scBr[i], rownames(rna@meta.data[which(rna@meta.data$cell_type != "AST"),])]
  nonmol = tibble(type = "non-AST", expr = nonmol)
  both = rbind(mol, nonmol)
  plots[[i]] = ggplot(both, aes(x = type, y = expr, fill = type)) +
    geom_boxplot(color = "black") +
    ylim(0, 8) +
    labs(
      title = sign_AST_vs_nonAST_scBr[i],
      x = "",
      y = "log norm. expr"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black")
    ) +
    stat_compare_means(label.y = 7, label.x = 1.15, size = 6)
}

sign_scBridge_genes_expr_bp = ggarrange(plotlist = plots, ncol = 2, nrow = 3)

ggsave(
  glue("{result_folder}scBr_MOL_markers-browser_example-expr_bps.pdf"),
  plot = sign_scBridge_genes_expr_bp,
  width = 12,
  height = 12,
  device = "pdf"
)

plots = list()
for(i in 1:length(sign_AST_vs_nonAST_scBr[1:6])) {
  mol = g4_counts[sign_AST_vs_nonAST_scBr[i], which(colnames(g4_counts) %in% barcodes_scbridge)]
  mol = tibble(type = "AST", g4 = mol)
  nonmol = g4_counts[sign_AST_vs_nonAST_scBr[i], which(colnames(g4_counts) %in% barcodes_scbridge_others)]
  nonmol = tibble(type = "non-AST", g4 = nonmol)
  both = rbind(mol, nonmol)
  both = both %>% dplyr::filter(g4 > 0)
  plots[[i]] = ggplot(both, aes(x = type, y = g4, fill = type)) +
    geom_boxplot(color = "black") +
    ylim(0, 1.5) +
    labs(
      title = sign_AST_vs_nonAST_scBr[i],
      x = "",
      y = "log norm. G4"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black")
    ) +
    stat_compare_means(label.y = 1.2, label.x = 1.15, size = 6)
}

sign_scBridge_genes_g4_bp = ggarrange(plotlist = plots, ncol = 2, nrow = 3)

ggsave(
  glue("{result_folder}scBr_AST_markers-browser_example-g4count_bps.pdf"),
  plot = sign_scBridge_genes_g4_bp,
  width = 12,
  height = 12,
  device = "pdf"
)

# volcano plot
volc_input = diff_AST_vs_nonAST_scBr %>% 
  mutate(gene_name = rownames(.)) 

volc_input = volc_input %>% mutate(group = case_when(
  avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "up",
  avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "down",
  avg_log2FC >= -0.5 & avg_log2FC <= 0.5 ~ "unaltered",
  p_val_adj > 0.05 ~ "unaltered"
)) %>%
  mutate(sign_label = case_when(avg_log2FC > 0.5 & p_val_adj < 0.05 ~ gene_name,
                                avg_log2FC < -0.5 & p_val_adj < 0.05 ~ gene_name,
                                avg_log2FC >= -0.5 & avg_log2FC <= 0.5 ~ ""))
labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

# plot
ggplot_volc = volc_input %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-6, 6, 2)),  	 
                     limits = c(-6, 6)) +
  scale_y_continuous(breaks = c(seq(0, 10, 2)),  	 
                     limits = c(0, 10)) +
  labs(
    title = "predicted AST vs other G4 differences",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 7, max.overlaps = 100) # add labels
ggplot_volc

# export in pdf and png
ggsave(
  glue("{result_folder}scBr_predAST_vs_nonAST_volc.pdf"),
  width = 10,
  height = 10,
  device = "pdf"
)
ggsave(
  glue("{result_folder}scBr_predAST_vs_nonAST_volc.png"),
  width = 10,
  height = 10,
  dpi = 300
)
