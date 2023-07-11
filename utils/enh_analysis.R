# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("GenomicFeatures")
  library("GenomicRanges")
  library("wigglescout")
  library("plyranges")
  library("wigglescout")
  library("ggpubr")
  library("ggrastr")
  #library("ArchR")
})

# export folder
result_folder = "../results/GenomicRanges/unsorted_outputs/"
# cluster-spec bigwigs
bigwigs = "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"

# Marek K27ac
k27ac = "../data/Marek_data/histone_scCnT/k27ac/possorted_RPGC.bigwig"
# mESC K27ac
mesc_k27ac = "../data/bw/Martire2019_ESC_H33WT_H3K27ac.bw"
mesc_k27ac2 = "../data/bw/Creyghton_et_al_mESC_H3K27ac.bw"

# data
g4_peaks = "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_peaks/peaks_per_clusters.bed"
cm_enh = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
gl_enh = "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed"
ltrs = "../data/bed/RepMasker_lt200bp.LTRIS2.bed"

# Cruz-Molina et al. mESC p300 ChIP-Seq peaks
cm_p300 = "../data/GSE89211_CruzMolina/P300_results/mESC_p300_bothR1_R2.narrowPeak"

# cicero co-G4 networks extended with +/- 1000 bp
load(file = "../results/Seurat/callpeaks_unsorted/unsorted_cicero_coG4networks.Rds")
CCAN_assigns = CCAN_assigns %>% separate(Peak, c("chr", "start", "end"), sep = "-") %>% 
  mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
  mutate(start = start - 1000, end = end + 1000)

CCAN_assigns$type = "Cicero, co-G4 network"
CCAN_assigns = GRanges(
  seqnames = CCAN_assigns$chr,
  ranges = IRanges(
    start = CCAN_assigns$start,
    end = CCAN_assigns$end,
    names = CCAN_assigns$type
  )
)

## data exploration on G4 peak set
# G4 peak set
g4_peaks = fread(g4_peaks)
g4_peaks$V7 = "G4"

# enhancer sets (mESC Cruz-Molina and Glaser)
cm_enh = fread(cm_enh)
cm_enh$V6 = "Cruz-Molina_active_enh"
cm_enh = GRanges(
  seqnames = cm_enh$V1,
  ranges = IRanges(
    start = cm_enh$V2,
    end = cm_enh$V3,
    names = cm_enh$V6
  )
)

# p300
cm_p300 = fread(cm_p300)
cm_p300$V7 = "Cruz-Molina_p300"
cm_p300 = GRanges(
  seqnames = cm_p300$V1,
  ranges = IRanges(
    start = cm_p300$V2,
    end = cm_p300$V3,
    names = cm_p300$V7
  )
)

# both p300 and enhancer
p300_enh_ol = findOverlaps(cm_p300, cm_enh, type = "any", ignore.strand = FALSE)
cm_enh_p300 = cm_p300[queryHits(p300_enh_ol)]


# Marek's brain K27ac at mESC Cruz-Molina enhancers
k27ac_enh = bw_loci(k27ac, loci = cm_enh_p300)
k27ac_enh = as_tibble(k27ac_enh)
k27ac_enh = k27ac_enh %>% dplyr::filter(possorted_RPGC > quantile(k27ac_enh$possorted_RPGC, 0.75))
k27ac_enh$type = "Cruz-Molina_active_enh"
k27ac_enh = GRanges(
  seqnames = k27ac_enh$seqnames,
  ranges = IRanges(
    start = k27ac_enh$start,
    end = k27ac_enh$end,
    names = k27ac_enh$type
  )
)

# Highly acetylated active enhancers colocalized with co-G4 sites
k27ac_coG4 = findOverlaps(k27ac_enh, CCAN_assigns, type = "any", ignore.strand = FALSE)
k27ac_coG4 = k27ac_enh[queryHits(k27ac_coG4)]

load(file = "../results/Seurat/callpeaks_unsorted/unsorted_cicero.Rds")

# visualization
gene_anno = rtracklayer::readGFF("../data/gencode.vM10.annotation.gff3")

gene_anno$chromosome = paste0("chr", gene_anno$seqid)
gene_anno$gene = gene_anno$gene_id
gene_anno$transcript = gene_anno$transcript_id
gene_anno$symbol = gene_anno$gene_name

pdf(file = "../results/genome_browser/Figure_3/unsorted_brain/cicero/Cicero_CCAN_chr3_9485472-9555738.pdf",   
    width = 4, 
    height = 4) 
plot_connections(conns, "chr5", 64662854, 64743058,
                 gene_model = as.data.frame(gene_anno),
                 coaccess_cutoff = .0,
                 connection_width = .5)
dev.off()


# mESC K27ac
mesc_k27ac_enh = bw_loci(mesc_k27ac, loci = cm_enh_p300)
mesc_k27ac_enh = as_tibble(mesc_k27ac_enh)
mesc_k27ac_enh = mesc_k27ac_enh %>%
  na.omit()
mesc_k27ac_enh = mesc_k27ac_enh %>%
  dplyr::filter(
    Martire2019_ESC_H33WT_H3K27ac > quantile(mesc_k27ac_enh$Martire2019_ESC_H33WT_H3K27ac, 0.75)
  )
mesc_k27ac_enh$type = "Martire_active_enh"
mesc_k27ac_enh = GRanges(
  seqnames = mesc_k27ac_enh$seqnames,
  ranges = IRanges(
    start = mesc_k27ac_enh$start,
    end = mesc_k27ac_enh$end,
    names = mesc_k27ac_enh$type
  )
)

mesc_k27ac_enh2 = bw_loci(mesc_k27ac2, loci = cm_enh_p300)
mesc_k27ac_enh2 = as_tibble(mesc_k27ac_enh2)
mesc_k27ac_enh2 = mesc_k27ac_enh2 %>%
  na.omit()
mesc_k27ac_enh2 = mesc_k27ac_enh2 %>% dplyr::filter(
  Creyghton_et_al_mESC_H3K27ac > quantile(mesc_k27ac_enh2$Creyghton_et_al_mESC_H3K27ac, 0.75)
)
mesc_k27ac_enh2$type = "Creyghton_active_enh"
mesc_k27ac_enh2 = GRanges(
  seqnames = mesc_k27ac_enh2$seqnames,
  ranges = IRanges(
    start = mesc_k27ac_enh2$start,
    end = mesc_k27ac_enh2$end,
    names = mesc_k27ac_enh2$type
  )
)

# Glaser et al enhancer set
gl_enh = fread(gl_enh)
gl_enh$V6 = "Glaser_active_enh"
gl_enh = GRanges(
  seqnames = gl_enh$V1,
  ranges = IRanges(
    start = gl_enh$V2,
    end = gl_enh$V3,
    names = gl_enh$V6
  )
)

### wigglescout analysis
scatter_function = function(bigwig) {
  p = plot_bw_bins_scatter(
    x = glue("{bigwigs}{bigwig}"),
    y = k27ac,
    genome = "mm10",
    selection = cm_enh_p300,
    verbose = FALSE
  ) +
    ggrastr::geom_point_rast(color = "#636363") +
    labs(
      title = bigwig,
      x = "G4",
      y = "H3K27ac (Bartosovic et al.)",
      color = ""
    ) +
    xlim(0,30) +
    ylim(0,30) +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 15, color = "black"),
      axis.text.y = element_text(size = 15, color = "black")
    ) + stat_cor(
      method = "pearson",
      label.x = 15,
      label.y = 20,
      size = 5,
      p.accuracy = 0.001
    )
  p = rasterize(p, layers='Point', dpi=100)
  return(p)
}

bigwig_samples = list.files(bigwigs, full.names = TRUE)
plot_input = bw_loci(bigwig_samples, loci = k27ac_enh)
plot_input = as_tibble(plot_input)
plot_input_res0.8 = plot_input %>% 
  dplyr::select(seqnames, start, end, "0" = X0.bam_RPGC, "1" = X1.bam_RPGC, "2" = X2.bam_RPGC, "3" = X3.bam_RPGC, "4" = X4.bam_RPGC) %>% 
  pivot_longer("0":"4", names_to = "cluster", values_to = "G4") 
  # dplyr::filter(cluster != "0 - oligodendr") %>% 
  # dplyr::filter(cluster != "1 - misc. brain") 

plot_input_res0.8 = plot_input_res0.8 %>% mutate(k27ac_source = "Bartosovic et al.")

mesc_plot_input = bw_loci(bigwig_samples, loci = mesc_k27ac_enh)
mesc_plot_input = as_tibble(mesc_plot_input)
mesc_plot_input_res0.8 = mesc_plot_input %>% 
  dplyr::select(seqnames, start, end, "0" = X0.bam_RPGC, "1" = X1.bam_RPGC, "2" = X2.bam_RPGC, "3" = X3.bam_RPGC, "4" = X4.bam_RPGC) %>%  
  pivot_longer("0":"4", names_to = "cluster", values_to = "G4")  
  # dplyr::filter(cluster != "0 - oligodendr") %>% 
  # dplyr::filter(cluster != "1 - misc. brain") 

mesc_plot_input_res0.8 = mesc_plot_input_res0.8 %>% mutate(k27ac_source = "Martier et al.")

mesc_plot_input2 = bw_loci(bigwig_samples, loci = mesc_k27ac_enh2)
mesc_plot_input2 = as_tibble(mesc_plot_input)
mesc_plot_input_res0.8_2 = mesc_plot_input2 %>% 
  dplyr::select(seqnames, start, end, "0" = X0.bam_RPGC, "1" = X1.bam_RPGC, "2" = X2.bam_RPGC, "3" = X3.bam_RPGC, "4" = X4.bam_RPGC) %>%  
  pivot_longer("0":"4", names_to = "cluster", values_to = "G4")
  # dplyr::filter(cluster != "0 - oligodendr") %>% 
  # dplyr::filter(cluster != "1 - misc. brain") 

mesc_plot_input_res0.8_2 = mesc_plot_input_res0.8_2 %>% mutate(k27ac_source = "Creyghton et al.")

plot_input = rbind(mesc_plot_input_res0.8, mesc_plot_input_res0.8_2, plot_input_res0.8) 
plot_input = plot_input %>% dplyr::filter(G4 > 0)

order = plot_input_res0.8 %>% group_by(cluster) %>% summarise(mean = mean(G4)) %>% arrange(desc(mean)) %>% pull(cluster) 
order = factor(plot_input$cluster, levels = order)

grouped_jitter1 = ggplot(plot_input, aes(x = order, y = G4, color = k27ac_source)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5), alpha = 1, pch = 19) +
  scale_color_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  ylim(0, 200) +
  labs(
    title = "",
    x = "",
    y = "G4 signal",
    color = "H3K27ac"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) + stat_compare_means(aes(group = k27ac_source), label = "p.signif")
grouped_jitter1
grouped_jitter1 = rasterize(grouped_jitter1, layers='Point', dpi=300)

# make facet at level of clusters
# comparison table for pairwise statistics
comparisons = list(c("Martier et al.", "Creyghton et al."), 
                       c("Creyghton et al.", "Bartosovic et al."), 
                       c("Martier et al.", "Bartosovic et al."))

order = plot_input %>% group_by(k27ac_source) %>% summarise(mean = mean(G4)) %>% 
  arrange(desc(mean)) %>% pull(k27ac_source) 
order = factor(plot_input$k27ac_source, levels = order)

facet = ggplot(plot_input, aes(x = order, y = G4, color = k27ac_source)) +
  geom_boxplot() +
  facet_wrap(~ cluster) +
  scale_color_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  ylim(0, 200) +
  labs(
    title = "",
    x = "",
    y = "G4 signal",
    color = "H3K27ac"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) + stat_compare_means(comparisons = comparisons, label.y = c(150, 165, 180),
                         label = "p.signif")
facet
facet = rasterize(facet, layers='Point', dpi=300)

ggsave(
  glue("{result_folder}CM_enhancers-G4s_over_high_K27ac.png"),
  plot = grouped_jitter1,
  width = 8,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}CM_enhancers-G4s_over_high_K27ac.pdf"),
  plot = grouped_jitter1,
  device = "pdf",
  width = 8,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}CM_enhancers-G4s_over_high_K27ac-facet.png"),
  plot = facet,
  width = 8,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}CM_enhancers-G4s_over_high_K27ac-facet.pdf"),
  plot = facet,
  device = "pdf",
  width = 8,
  height = 5,
  dpi = 300,
)

order = plot_input_res0.8 %>% group_by(cluster) %>% summarise(mean = mean(G4)) %>% arrange(desc(mean)) %>% pull(cluster) 
order = factor(plot_input_res0.8$cluster, levels = order)

jitter1 = ggplot(plot_input_res0.8, aes(x = order, y = G4)) +
 geom_jitter(color = "#bdbdbd") +
  ylim(0, 200) +
  labs(
    title = "Cruz-Molina et al. enhancers with high K27ac in mouse brain",
    x = "",
    y = "G4 signal"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 5),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) + stat_compare_means(label = "p.signif")
jitter1
jitter1 = rasterize(jitter1, layers='Point', dpi=300)

ggsave(
  glue("{result_folder}CM_enhancers-G4s_over_high_Bartosovic_K27ac.png"),
  plot = jitter1,
  width = 8,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}CM_enhancers-G4s_over_high_Bartosovic_K27ac.pdf"),
  plot = jitter1,
  device = "pdf",
  width = 8,
  height = 5,
  dpi = 300,
)


bigwig_samples = c("0.bam_RPGC.bigwig", "1.bam_RPGC.bigwig", "2.bam_RPGC.bigwig", "3.bam_RPGC.bigwig", "4.bam_RPGC.bigwig")
scatters = lapply(bigwig_samples, scatter_function)
scatters = ggarrange(plotlist = scatters)
scatters

ggsave(
  glue("{result_folder}CM_enhancers-K27ac_scatters.png"),
  plot = scatters,
  width = 12,
  height = 8,
  dpi = 300,
)

ggsave(
  glue("{result_folder}CM_enhancers-K27ac_scatters.pdf"),
  plot = scatters,
  device = "pdf",
  width = 12,
  height = 8,
  dpi = 300,
)

# G4 length distributions
hist = g4_peaks %>% mutate(diff = V3 - V2) %>%
  dplyr::filter(V6 %in% c("0", "1", "2", "3", "4")) %>%
  ggplot(., aes(x = diff, fill = V6)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  xlim(200, 750) +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "G4 length distributions / Seurat cluster",
       x = "length (bp)",
       y = "Density",
       fill = "Seurat cluster") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
hist

ggsave(
  glue("{result_folder}G4_length_distr_Seurat_clst.pdf"),
  plot = hist,
  width = 10,
  height = 7,
  device = "pdf"
)

# convert G4 peak ranges to GRanges
g4_peaks =
  lapply(split(g4_peaks, g4_peaks$V6), function(i) {
    GRanges(seqnames = i$V1,
            ranges = IRanges(
              start = i$V2,
              end = i$V3,
              names = i$V7
            ))
  })


## find overlaps between G4 peaks and active enhancers
# query: g4_peaks
# subject: enhancers
g4s_no_ol = g4_peaks[c("0", "1", "2", "3", "4")]
g4_cm_enh_quant = numeric()
for (cluster in names(g4s_no_ol)) {
  ol = suppressWarnings(findOverlaps(g4_peaks[[cluster]], cm_enh, type = "any", minoverlap = 1))
  g4_cm_enh_quant = c(g4_cm_enh_quant, length(ol))
}

# query: g4_peaks
# subject: enhancers with Marek K27ac signal
g4_k27ac_cm_enh_quant = numeric()
for (cluster in names(g4s_no_ol)) {
  ol = suppressWarnings(findOverlaps(g4_peaks[[cluster]], k27ac_enh, type = "any", minoverlap = 1))
  g4_k27ac_cm_enh_quant = c(g4_k27ac_cm_enh_quant, length(ol))
}

get_signals_over_k27ac_enh = function(selected_cluster = "4") {
  require("wigglescout")
  cl_cm_ol = suppressWarnings(subsetByOverlaps(g4_peaks[[selected_cluster]], k27ac_enh, type = "any", minoverlap = 1))
  cl_cm_ol = as.tibble(cl_cm_ol)
  
  bigwigs = c(
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/0.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/1.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/2.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/3.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/4.bam_RPGC.bigwig",
    k27ac
  )
  
  bed = cl_cm_ol[, 1:3]
  write_tsv(
    bed,
    glue(
      "../data/bed/Cruz_Molina_enh-unsorted_cluster{selected_cluster}.bed"
    ),
    col_names = FALSE
  )
  bed = glue("../data/bed/Cruz_Molina_enh-unsorted_cluster{selected_cluster}.bed")
  selected_bw = glue("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/{selected_cluster}.bam_RPGC.bigwig")
  bigwigs = c(bigwigs[which(bigwigs != selected_bw)], selected_bw)
  
  read_cov = bw_loci(bigwigs, loci = bed)
  read_cov = as.data.frame(read_cov)
  
  columns = colnames(read_cov)[grepl("X", colnames(read_cov)[which(colnames(read_cov) != glue("X{selected_cluster}.bam_RPGC"))])]
  read_cov_filt = read_cov %>% 
    mutate(average = rowMeans(dplyr::select(., columns))) %>% 
    mutate(enrichment = get(glue("X{selected_cluster}.bam_RPGC")) / average) %>% 
    dplyr::filter(enrichment >= 10)
  
  return(read_cov_filt)
  
}

get_signals_over_CCAN_enh = function(selected_cluster = "4") {
  require("wigglescout")
  cl_cicero_ol = suppressWarnings(subsetByOverlaps(g4_peaks[[selected_cluster]], k27ac_coG4, type = "any", minoverlap = 1))
  cl_cicero_ol = as.tibble(cl_cicero_ol)
  
  bigwigs = c(
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/0.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/1.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/2.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/3.bam_RPGC.bigwig",
    "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/4.bam_RPGC.bigwig",
    k27ac
  )
  
  bed = cl_cicero_ol[, 1:3]
  write_tsv(
    bed,
    glue(
      "../data/bed/cicero_CCAN_enh-unsorted_cluster{selected_cluster}.bed"
    ),
    col_names = FALSE
  )
  bed = glue("../data/bed/cicero_CCAN_enh-unsorted_cluster{selected_cluster}.bed")
  selected_bw = glue("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/{selected_cluster}.bam_RPGC.bigwig")
  bigwigs = c(bigwigs[which(bigwigs != selected_bw)], selected_bw)
  
  read_cov = bw_loci(bigwigs, loci = bed)
  read_cov = as.data.frame(read_cov)
  
  columns = colnames(read_cov)[grepl("X", colnames(read_cov)[which(colnames(read_cov) != glue("X{selected_cluster}.bam_RPGC"))])]
  read_cov_filt = read_cov %>% 
    mutate(average = rowMeans(dplyr::select(., columns))) %>% 
    mutate(enrichment = get(glue("X{selected_cluster}.bam_RPGC")) / average) %>% 
    dplyr::filter(enrichment >= 10)
  
  return(read_cov_filt)
  
}


# get signals test
cl_4_spec = get_signals_over_k27ac_enh()
cl_3_spec = get_signals_over_k27ac_enh(selected_cluster = "3")
cl_0_spec_CCAN = get_signals_over_CCAN_enh(selected_cluster = "0")
cl_1_spec_CCAN = get_signals_over_CCAN_enh(selected_cluster = "1")
cl_2_spec_CCAN = get_signals_over_CCAN_enh(selected_cluster = "2")
cl_3_spec_CCAN = get_signals_over_CCAN_enh(selected_cluster = "3")
cl_4_spec_CCAN = get_signals_over_CCAN_enh(selected_cluster = "4")

g4_gl_enh_quant = numeric()
for (cluster in names(g4s_no_ol)) {
  ol = suppressWarnings(findOverlaps(g4_peaks[[cluster]], gl_enh, type = "any", minoverlap = 1))
  g4_gl_enh_quant = c(g4_gl_enh_quant, length(ol))
}

clusters = c("0", "1", "2", "3", "4")
cm_bar = tibble(
  overlap = g4_k27ac_cm_enh_quant,
  peak_count = sapply(clusters, function(x) length(g4_peaks[[which(names(g4_peaks) == x)]])),
  enhancer_ratio = (overlap / peak_count) * 100,
  cluster = names(g4s_no_ol),
  fill_col = names(g4s_no_ol)
) %>%
  ggplot(data = ., aes(
    x = reorder(cluster,-enhancer_ratio),
    y = enhancer_ratio,
    fill = fill_col
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, 1)) +
  labs(
    title = expression(paste(
      "G4 overlaps with active enhancers of ",
      italic("Cruz-Molina et al.")
    )),
    x = "Seurat cluster",
    y = "active enhancer %",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black")
  )
cm_bar

cm_bar_input = tibble(
  overlap = g4_k27ac_cm_enh_quant,
  peak_count = sapply(clusters, function(x) length(g4_peaks[[which(names(g4_peaks) == x)]])),
  enhancer_ratio = (overlap / peak_count) * 100,
  cluster = names(g4s_no_ol),
  fill_col = names(g4s_no_ol)
)

print(glue("Average enhancer percentage across clusters: {round(mean(cm_bar_input$enhancer_ratio), 2)} % (SD: {round(sd(cm_bar_input$enhancer_ratio), 2)})"))
print(glue("Average peak number over active enhancer: {round(mean(cm_bar_input$overlap), 0)} (SD: {round(sd(cm_bar_input$overlap), 0)})"))

ggsave(
  glue("{result_folder}G4_overlaps_w_CruzM_active_enh.pdf"),
  plot = cm_bar,
  width = 10,
  height = 7,
  device = "pdf"
)

clusters = c("0", "1", "2", "3", "4")
gl_bar = tibble(
  overlap = g4_gl_enh_quant,
  peak_count = sapply(clusters, function(x) length(g4_peaks[[which(names(g4_peaks) == x)]])),
  enhancer_ratio = (overlap / peak_count) * 100,
  cluster = names(g4s_no_ol),
  fill_col = names(g4s_no_ol)
) %>%
  ggplot(data = ., aes(
    x = reorder(cluster,-enhancer_ratio),
    y = enhancer_ratio,
    fill = fill_col
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title = expression(paste(
      "G4 overlaps with active enhancers of ",
      italic("Glaser et al.")
    )),
    x = "Seurat cluster",
    y = "active enhancer %",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
gl_bar

ggsave(
  glue("{result_folder}G4_overlaps_w_Glaser_active_enh.pdf"),
  plot = gl_bar,
  width = 10,
  height = 7,
  device = "pdf"
)

