# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ggfortify")
  library("ggrepel")
  library("Seurat")
})

# result / peak folder
result_folder = "../results/GenomicRanges/"
bigwig_folder = "../data/bw/"
bed_folder = "../data/bed/"

# get variable features of mESC-MEF data
mesc_mef = readRDS("../results/Seurat/callpeaks_mESC-MEF/mESC-MEF_res0.1.Rds")
x = mesc_mef@meta.data %>% select(starts_with("nFeature_bins"))

mesc_mef = FindVariableFeatures(
  mesc_mef,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  nfeatures = 150,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE
)

top_features = mesc_mef@assays$peaks@var.features

bed = tibble(V1 = sapply(top_features, function(x) strsplit(x, "-")[[1]][1]),
             V2 = sapply(top_features, function(x) strsplit(x, "-")[[1]][2]),
             V3 = sapply(top_features, function(x) strsplit(x, "-")[[1]][3]))

write_tsv(bed, glue("{bed_folder}mESC-MEF_scCnT-150_top_features.bed"), col_names = FALSE)
mesc_mef_bed = glue("{bed_folder}mESC-MEF_scCnT-150_top_features.bed")

unsorted = readRDS("../results/Seurat/callpeaks_unsorted/unsorted.Rds")
unsorted = FindVariableFeatures(
  unsorted,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  nfeatures = 150,
  num.bin = 20,
  binning.method = "equal_width",
  verbose = TRUE
)

top_features = unsorted@assays$peaks@var.features

bed = tibble(V1 = sapply(top_features, function(x) strsplit(x, "-")[[1]][1]),
             V2 = sapply(top_features, function(x) strsplit(x, "-")[[1]][2]),
             V3 = sapply(top_features, function(x) strsplit(x, "-")[[1]][3]))


write_tsv(bed, glue("{bed_folder}unsorted_scCnT-150_top_features.bed"), col_names = FALSE)
unsorted_bed = glue("{bed_folder}unsorted_scCnT-150_top_features.bed")
             
# bigwig files
bigwigs = c(
  glue(
    "{bigwig_folder}bulk_G4_CnT_MEF_rep1_R1.mLb.clN_mm10.bigWig"
  ),
  glue(
    "{bigwig_folder}bulk_G4_CnT_MEF_rep2_R1.mLb.clN_mm10.bigWig"
  ),
  glue(
    "{bigwig_folder}bulk_G4_CnT_mESC_rep1_R1.mLb.clN_mm10.bigWig"
  ),
  glue(
    "{bigwig_folder}bulk_G4_CnT_mESC_rep2_R1.mLb.clN_mm10.bigWig"
  ),
  glue(
    "../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/0_res0.1.bigwig"
  ),
  glue(
    "../results/Seurat/callpeaks_mESC-MEF/cluster_bigwigs/1_res0.1.bigwig"
  )
)

# create read coverage matrix by 10 kb resolution
read_cov = bw_bins(bigwigs, bin_size = 10000, genome = "mm10")
read_cov = as.data.frame(read_cov)
mat = read_cov %>% select(starts_with("X"), starts_with("bulk")) %>%
  select(
    "mESC-MEF cluster 0" = "X0_res0.1",
    "mESC-MEF cluster 1" = "X1_res0.1",
    "MEF bulk, rep 1" = "bulk_G4_CnT_MEF_rep1_R1.mLb.clN_mm10",
    "MEF bulk, rep 2" = "bulk_G4_CnT_MEF_rep2_R1.mLb.clN_mm10",
    "mESC bulk, rep 1" = "bulk_G4_CnT_mESC_rep1_R1.mLb.clN_mm10",
    "mESC bulk, rep 2" = "bulk_G4_CnT_mESC_rep2_R1.mLb.clN_mm10",
  )
mat = mat %>% na.omit()
mat = t(as.matrix(mat))

# generate PCA
pca = prcomp(mat, scale. = TRUE)

# add colors
values = c("black", "#636363", "#bdbdbd", "#f0f0f0", "#9ecae1", "#fc9272")

# plot PCA
p = autoplot(pca,
             data = mat,
             loadings = FALSE,
             label = FALSE) +
  geom_point(size = 6, aes(color = factor(rownames(pca$x), levels = c("MEF bulk, rep 1", 
                                                                      "MEF bulk, rep 2",
                                                                      "mESC bulk, rep 1",
                                                                      "mESC bulk, rep 2",
                                                                      "mESC-MEF cluster 0", 
                                                                      "mESC-MEF cluster 1"
                                                                      )))) +
  labs(title = "PCA on 10 kb bins", color = " ") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 14, color = "black")
  ) +
  scale_color_manual(values = values)
p

ggsave(
  glue("{result_folder}PCA_10kb_bins_with_reps.png"),
  width = 7,
  height = 7,
  dpi = 500,
)
ggsave(
  glue("{result_folder}PCA_10kb_bins_with_reps.pdf"),
  width = 7,
  height = 7,
  device = "pdf"
)

# PCA on 150 top variable feature of mESC-MEF dataset
top_features = bw_loci(bigwigs, loci = "../data/bed/mESC-MEF_scCnT-150_top_features.bed")
top_features = as.data.frame(top_features)

top_mat = top_features %>% select(starts_with("X"), starts_with("bulk")) %>%
  select(
    "mESC-MEF cluster 0" = "X0_res0.1",
    "mESC-MEF cluster 1" = "X1_res0.1",
    "MEF bulk, rep 1" = "bulk_G4_CnT_MEF_rep1_R1.mLb.clN_mm10",
    "MEF bulk, rep 2" = "bulk_G4_CnT_MEF_rep2_R1.mLb.clN_mm10",
    "mESC bulk, rep 1" = "bulk_G4_CnT_mESC_rep1_R1.mLb.clN_mm10",
    "mESC bulk, rep 2" = "bulk_G4_CnT_mESC_rep2_R1.mLb.clN_mm10",
  )
top_mat = top_mat %>% na.omit()
top_mat = t(as.matrix(top_mat))

# generate PCA
pca = prcomp(top_mat, scale. = TRUE)

# add colors
values = c("black", "#636363", "#bdbdbd", "#f0f0f0", "#9ecae1", "#fc9272")

# plot PCA
top_p = autoplot(pca,
             data = top_mat,
             loadings = FALSE,
             label = FALSE) +
  geom_point(size = 6, aes(color = factor(rownames(pca$x), levels = c("MEF bulk, rep 1", 
                                                                      "MEF bulk, rep 2",
                                                                      "mESC bulk, rep 1",
                                                                      "mESC bulk, rep 2",
                                                                      "mESC-MEF cluster 0", 
                                                                      "mESC-MEF cluster 1"
  )))) +
  labs(title = "PCA on top 150 most variable features", color = " ") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 14, color = "black")
  ) +
  scale_color_manual(values = values)
top_p

ggsave(
  glue("{result_folder}PCA_top150features_with_reps.png"),
  width = 7,
  height = 7,
  dpi = 500,
)
ggsave(
  glue("{result_folder}PCA_top150features_with_reps.pdf"),
  width = 7,
  height = 7,
  device = "pdf"
)

# PCA on marker regions (Marek's method)
# bigwig files
bigwig_folder = "../data/Marek_data/mESC-MEF/"

bigwigs = c(
  glue(
    "{bigwig_folder}cluster_0.bw"
  ),
  glue(
    "{bigwig_folder}cluster_1.bw"
  ),
  glue(
    "{bigwig_folder}cluster_2.bw"
  ),
  glue(
    "{bigwig_folder}G4_H33WT_SL_CnT_R1.mm10.bw"
  ),
  glue(
    "{bigwig_folder}G4_H33WT_SL_CnT_R2.mm10.bw"
  ),
  glue(
    "{bigwig_folder}G4_MEF_CnT_R1.mm10.bw"
  ),
  glue(
    "{bigwig_folder}G4_MEF_CnT_R2.mm10.bw"
  )
)

bed = glue(
  "{bigwig_folder}regions.bed"
)

# create read coverage matrix over marker regions
read_cov = bw_loci(bigwigs, loci = bed)
read_cov = as.data.frame(read_cov)
markers_mat = read_cov %>% select(starts_with("cluster"), starts_with("G4")) %>%
  select(
    "mESC-MEF cluster 0" = "cluster_0",
    "mESC-MEF cluster 1" = "cluster_1",
    "mESC-MEF cluster 2" = "cluster_2",
    "MEF bulk, rep 1" = "G4_MEF_CnT_R1.mm10",
    "MEF bulk, rep 2" = "G4_MEF_CnT_R2.mm10",
    "mESC bulk, rep 1" = "G4_H33WT_SL_CnT_R1.mm10",
    "mESC bulk, rep 2" = "G4_H33WT_SL_CnT_R2.mm10",
  )


markers_mat = markers_mat %>% na.omit()
markers_mat = t(as.matrix(markers_mat))

# generate PCA
pca = prcomp(markers_mat, scale. = TRUE)

# add colors
values = c("black", "#636363", "#bdbdbd", "#f0f0f0", "#9ecae1", "#fc9272", "red")

# plot PCA
markers_p = autoplot(pca,
                 data = markers_mat,
                 loadings = FALSE,
                 label = FALSE) +
  geom_point(size = 6, aes(color = factor(rownames(pca$x), levels = c("MEF bulk, rep 1", 
                                                                      "MEF bulk, rep 2",
                                                                      "mESC bulk, rep 1",
                                                                      "mESC bulk, rep 2",
                                                                      "mESC-MEF cluster 0", 
                                                                      "mESC-MEF cluster 1",
                                                                      "mESC-MEF cluster 2"
  )))) +
  labs(title = "PCA on feature markers", color = " ") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 14, color = "black")
  ) +
  scale_color_manual(values = values)
markers_p

ggsave(
  glue("{result_folder}PCA_on_markers_with_reps.png"),
  width = 7,
  height = 7,
  dpi = 500,
)
ggsave(
  glue("{result_folder}PCA_on_markers_with_reps.pdf"),
  width = 7,
  height = 7,
  device = "pdf"
)
