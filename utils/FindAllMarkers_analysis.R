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
  library("ggpubr")
  library("ggrastr")
})

set.seed(5)

# helper function
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# folders
result_folder = "../results/Seurat/callpeaks_unsorted/"
marker_genes = "../data/GSE75330/marker_genes.txt" # Marques et al. scRNA-Seq marker genes
marques_marker_genes = fread(marker_genes, header = FALSE)

# scRNA-Seq data
scrna = "../results/Seurat/scRNASeq_GSM4979874-75.rds"
marques_scrna = "../results/Seurat/scRNASeq_GSE75330.rds"


# G4 scCut&Tag marker regions
# Wilcoxon rank sum based differential analysis
markers = fread("../results/Seurat/final/unsorted_brain/res0.8/outputs/FindAllMarkers_logreg_output.tsv")
markers = markers %>% rownames_to_column() %>% mutate(rowname = as.character(rowname))
markers_annot = mm10_annotation(markers, seqname_col = "chr", start_col = "start", end_col = "end", feature_1 = "rowname", feature_2 = "rowname")
markers = markers %>% inner_join(., markers_annot, by = c("rowname" = "feature_1")) %>% 
  dplyr::select(chr, start = start.x, end = end.x, p_val, p_val_adj, avg_log2FC, cluster, distanceToTSS, gene = SYMBOL)
write_tsv(markers, "../results/Seurat/callpeaks_unsorted/FindAllMarkers_logreg_output_annot.tsv")
markers = markers %>% na.omit()

# logistic regression based differential analysis
markers_lr = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_logreg_output.tsv")
markers_lr = markers_lr %>% mutate(region = paste(chr, start, end, sep = "-"))
markers_lr_annot = mm10_annotation(regions = markers_lr, start_col = "start", end_col = "end", seqname_col = "chr") %>% 
  dplyr::filter(abs(distanceToTSS) < 3000) %>% dplyr::select(seqnames, start, end, gene_symbol = SYMBOL) %>% 
  mutate(region = paste(seqnames, start, end, sep = "-"))
markers_lr = markers_lr %>% inner_join(., markers_lr_annot, by = "region") %>% 
  dplyr::select(region, gene_symbol, avg_log2FC, p_val, p_val_adj, cluster) %>% distinct_all()

mean_fc_up = markers_lr %>% dplyr::filter(avg_log2FC > 0) %>% pull(avg_log2FC) %>% mean
mean_fc_down = markers_lr %>% dplyr::filter(avg_log2FC < 0) %>% pull(avg_log2FC) %>% mean
print(paste0("Mean up G4-quadruplexed fold change: ", as.character(round(mean_fc_up, 2))))
print(paste0("Mean down G4-quadruplexed fold change: ", as.character(round(mean_fc_down, 2))))

volc_input = markers_lr %>% 
  # group_by(gene_symbol) %>%
  # dplyr::filter(avg_log2FC == max(abs(avg_log2FC), na.rm=TRUE)) %>% 
  mutate(group = case_when(
  avg_log2FC > 0.1 & p_val_adj < 0.05 ~ "up",
  avg_log2FC < -0.1 & p_val_adj < 0.05 ~ "down",
  avg_log2FC >= -0.1 & avg_log2FC <= 0.1 ~ "unaltered"
))
volc_input = volc_input %>% mutate(sign_label = case_when(
  avg_log2FC > 0.5 & p_val_adj < 0.05 ~ gene_symbol,
  avg_log2FC < -0.5& p_val_adj < 0.05 ~ gene_symbol,
  avg_log2FC >= -0.5 & avg_log2FC <= 0.5 ~ ""
))
  
labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

# plot
ggplot_volc = volc_input %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             fill = as.character(cluster),    
             size = group,
             alpha = group)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-0.1, 0.1),
             linetype = "dashed") +
  #scale_fill_manual(values = cols) + # Modify point colour
  scale_fill_manual(values = c("#bcbada", "#ffffb3", "#8dd3c7", "#80b1d3")) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4.5, 4.5, 1)),  	 
                     limits = c(-4.5, 4.5)) +
  labs(
    title = "Differential G-quadruplexed regions (+/- 3 kb to TSS)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  ylim(0, 30) +
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
  geom_text_repel(label = labels, size = 6, max.overlaps = 12) # add labels
ggplot_volc

ggsave(
  glue("{result_folder}FindAllMarkers_volc_logreg.png"),
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}FindAllMarkers_volc_logreg.pdf"),
  device = "pdf",
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300,
)

# assign expression levels
scrna = readRDS(scrna)
marques_scrna = readRDS(marques_scrna)

interesting_genes = volc_input %>% arrange(desc(abs(avg_log2FC))) %>% pull(gene_symbol)
interesting_genes = interesting_genes[1:20]

# creating inputs for visualization
marques_lognorm = as.data.frame(marques_scrna[["RNA"]]@data)
#summary(rowMeans(marques_lognorm))

# plot 1 as described in question
ggplot(b, aes(x = population, y = value)) + geom_histogram(aes(fill = lambda), stat = "identity", position = "dodge") 

valid = c()
for(gene in interesting_genes) {
  if(gene %in% rownames(marques_lognorm)) {
    valid = c(valid, gene)
  }
}
interesting_genes = valid[1:10]

marques_lognorm_filt = marques_lognorm[interesting_genes,]
marques_bp_input = marques_lognorm_filt %>% t %>% as_tibble %>%
  pivot_longer(interesting_genes[1]:interesting_genes[10],
               values_to = "log2_norm_expr",
               names_to = "gene") %>% 
  mutate(source = "Marques et al.")

marques_means = rowMeans(marques_lognorm_filt)
marques_means = tibble(source = "Marques et al.", gene = names(marques_means), 
                       pseudobulked_mean_expr = round(marques_means, 3))
order = marques_means %>% arrange(desc(pseudobulked_mean_expr)) %>% pull(gene)

bartosovic_lognorm = as.data.frame(scrna[["RNA"]]@data)
summary(rowMeans(bartosovic_lognorm))
bartosovic_lognorm_filt = bartosovic_lognorm[interesting_genes,]

bartosovic_bp_input = bartosovic_lognorm_filt %>% t %>% as_tibble %>%
  pivot_longer(interesting_genes[1]:interesting_genes[10],
               values_to = "log2_norm_expr",
               names_to = "gene") %>% 
  mutate(source = "Bartosovic et al.")

# boxplot inputs
bp_input = rbind(bartosovic_bp_input, marques_bp_input)
bp_input = bp_input %>% dplyr::filter(log2_norm_expr > 0)

# heatmap input (aggr. log2 norm. expression)
bartosovic_means = rowMeans(bartosovic_lognorm_filt)
bartosovic_means = tibble(source = "Bartosovic et al.", gene = names(bartosovic_means), 
                       pseudobulked_mean_expr = round(bartosovic_means, 3))
  

means = rbind(bartosovic_means, marques_means)
order = factor(means$gene, levels = order)

# remove huge objects
rm(scrna); rm(marques_scrna)
rm(bartosovic_lognorm); rm(marques_lognorm)

# heatmap
hm = ggplot(means, aes(x = order, y = source, fill = pseudobulked_mean_expr)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#9ecae1",
    mid = "white",
    high = "#fc9272",
    midpoint = 1,
    limits = c(0, 2)
  ) +
  theme_classic() +
  xlab(label = "G4 enrichment") +
  ylab(label = "") +
  labs(fill = "log2 norm. expr.") +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 13,
      angle = 45,
      hjust = 1,
      vjust = 1
    ), 
    axis.text.y = element_text(color = "black", size = 10)
  ) +
  coord_fixed()
print(hm)

ggsave(
  glue("{result_folder}FindAllMarkers_hm_posmarkers.png"),
  plot = last_plot(),
  width = 5,
  height = 3.5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}FindAllMarkers_hm_posmarkers.pdf"),
  device = "pdf",
  plot = last_plot(),
  width = 5,
  height = 3.5,
  dpi = 300,
)

# boxplot
order = marques_means %>% arrange(desc(pseudobulked_mean_expr)) %>% pull(gene)
order = factor(bp_input$gene, levels = order)

jitter = ggplot(bp_input, aes(x = order, y = log2_norm_expr)) +
  #geom_boxplot(color = "black", outlier.shape = NA, outlier.stroke = NA, outlier.fill = NA)
  geom_jitter(aes(color = source), size = 1, alpha = 0.5) +
  scale_color_manual(values = c("#bdbdbd", "#fec44f")) +
  guides(
    alpha = FALSE,
    size = FALSE,
    color = guide_legend(override.aes = list(size = 5))
  ) +
  ylim(0, 5) +
  labs(
    title = "",
    x = "",
    y = "log2 norm. expr.",
    color = "scRNA-Seq"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(
      size = 12,
      color = "black",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(size = 12, color = "black")
  ) + stat_compare_means(aes(group = source), label = "p.signif")
jitter = rasterize(jitter, layers='Point', dpi=300)

ggsave(
  glue("{result_folder}FindAllMarkers_jitter_posmarkers.png"),
  plot = last_plot(),
  width = 5,
  height = 3.5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}FindAllMarkers_jitter_posmarkers.pdf"),
  device = "pdf",
  plot = last_plot(),
  width = 5,
  height = 3.5,
  dpi = 300,
)

# logistic regression based differential analysis
markers_lr_res0.1 = fread("../results/Seurat/final/unsorted_brain/res0.1/outputs/FindAllMarkers_logreg_output.tsv")
markers_lr_res0.1 = markers_lr_res0.1 %>% mutate(region = paste(chr, start, end, sep = "-"))
markers_lr_res0.1_annot = mm10_annotation(regions = markers_lr_res0.1, start_col = "start", end_col = "end", seqname_col = "chr") %>% 
  dplyr::filter(abs(distanceToTSS) < 3000) %>% 
  dplyr::select(seqnames, start, end, gene_symbol = SYMBOL) %>% 
  mutate(region = paste(seqnames, start, end, sep = "-"))
markers_lr_res0.1 = markers_lr_res0.1 %>% inner_join(., markers_lr_res0.1_annot, by = "region") %>% 
  dplyr::select(region, gene_symbol, avg_log2FC, p_val, p_val_adj, cluster) %>% distinct_all()


mean_fc_up = markers_lr_res0.1 %>% dplyr::filter(avg_log2FC > 0) %>% pull(avg_log2FC) %>% mean
mean_fc_down = markers_lr_res0.1 %>% dplyr::filter(avg_log2FC < 0) %>% pull(avg_log2FC) %>% mean
print(paste0("Mean up G4-quadruplexed fold change in res0.1 data: ", as.character(round(mean_fc_up, 2))))
print(paste0("Mean down G4-quadruplexed fold change in res0.1 data: ", as.character(round(mean_fc_down, 2))))

number_of_up = markers_lr_res0.1 %>% dplyr::filter(avg_log2FC > 0) %>% pull(gene_symbol) %>% unique %>% length
print(paste0("Number of up G4-quadruplexed promoters in res0.1 data: ", as.character(round(number_of_up, 2))))

volc_input_res0.1 = markers_lr_res0.1 %>% 
  group_by(gene_symbol) %>%
  dplyr::filter(avg_log2FC == max(abs(avg_log2FC), na.rm=TRUE)) %>% 
  mutate(group = case_when(
    avg_log2FC > 0.1 & p_val_adj < 0.05 ~ "up",
    avg_log2FC < -0.1 & p_val_adj < 0.05 ~ "down",
    avg_log2FC >= -0.1 & avg_log2FC <= 0.1 ~ "unaltered",
    TRUE ~ "non sign.")
  )
volc_input_res0.1 = volc_input_res0.1 %>% mutate(sign_label = case_when(
  avg_log2FC > 0.5 & p_val_adj < 0.005 ~ gene_symbol,
  avg_log2FC < -0.5& p_val_adj < 0.005 ~ gene_symbol,
  avg_log2FC >= -0.5 & avg_log2FC <= 0.5 ~ "",
  TRUE ~ "non sign."
))

labels = volc_input_res0.1 %>% pull(sign_label)

# plot
ggplot_volc_res0.1 = volc_input_res0.1 %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             fill = as.character(cluster),    
             size = group,
             alpha = group)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") +
  # geom_hline(yintercept = -log10(0.05),
  #            linetype = "dashed") +
  # geom_vline(xintercept = c(-0.1, 0.1),
  #            linetype = "dashed") +
  #scale_fill_manual(values = cols) + # Modify point colour
  scale_fill_manual(values = c("#bcbcbc", "#b1d689")) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-3, 3, 0.5)),  	 
                     limits = c(-3, 3)) +
  labs(
    title = "Differential G-quadruplexed regions (+/- 3 kb to TSS)",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  ylim(0, 30) +
  xlim(0.8, 3) +
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
  geom_text_repel(label = labels, size = 6, max.overlaps = 8) # add labels
ggplot_volc_res0.1

ggsave(
  glue("{Relb}FindAllMarkers_volc_logreg_res0.1.png"),
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}FindAllMarkers_volc_logreg_res0.1.pdf"),
  device = "pdf",
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300,
)

# analyse the whole genome not only the promoters
markers_lr = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_logreg_output.tsv")
markers_lr = markers_lr %>% mutate(region = paste(chr, start, end, sep = "-"))
markers_lr_annot = mm10_annotation(regions = markers_lr, start_col = "start", end_col = "end", seqname_col = "chr") %>% 
  dplyr::select(seqnames, start, end, gene_symbol = SYMBOL) %>% 
  mutate(region = paste(seqnames, start, end, sep = "-"))
markers_lr = markers_lr %>% inner_join(., markers_lr_annot, by = "region") %>% 
  dplyr::select(region, gene_symbol, avg_log2FC, p_val, p_val_adj, cluster) %>% distinct_all()

# testing cluster uniqueness by FindMarkers 
seurat = readRDS("../results/Seurat/final/unsorted_brain/res0.8/outputs/Seurat_object.Rds")

# in case of resolution 0.8
clusters = c(0, 1, 2, 3, 4)
for(i in clusters) {
  for (j in clusters) {
    if (i == j) {
      next
    }
    out = FindMarkers(seurat,
                      ident.1 = i,
                      ident.2 = j,
                      verbose = FALSE)
    sign = out %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(abs(avg_log2FC) > 0.1) %>%
      nrow
    if (sign == 0) {
      sign = 0
    }
    print(
      paste0(
        as.character(i),
        " vs ",
        as.character(j) ,
        " - ",
        " number of marker regions: ",
        as.character(sign)
      )
    )
  }
}

