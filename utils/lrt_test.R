suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("wigglescout")
  library("outliers")
  library("DESeq2")
  library("ComplexHeatmap")
  library("circlize")
})

# result / peak folder
result_folder = "../results/Seurat/callpeaks_unsorted/"
bigwig_folder = "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/"
peak_sets = "../results/Seurat/callpeaks_unsorted/peak_sets/"

# make metadata for DESeq2 LRT calculation
coldata = "../results/Seurat/callpeaks_unsorted/DESeq2_coldata.txt"
coldata = fread(coldata)

cluster_bigwigs = list.files(bigwig_folder, pattern = "*.bw", full.names = TRUE)
peak_set = "enhancer_analysis_output.bed"
bed = glue("{peak_sets}{peak_set}")

ranges = fread(bed)
ranges = ranges %>% mutate(row_id = row.names(.)) %>% mutate(name = paste0(V1, "_", V2, "_", V3, "_", row_id))

# create read coverage matrix
read_cov = bw_loci(cluster_bigwigs, loci = bed)
read_cov = as.data.frame(read_cov)
mat = read_cov %>% select(starts_with("X")) %>% select(
  "0" = X0.bam,
  "1" = X1.bam,
  "2" = X2.bam,
  "3" = X3.bam,
  "4" = X4.bam,
  "5" = X5.bam,
  "6" = X6.bam
)
rownames(mat) = ranges$name

# DESeq2 workflow
mat_rounded = round(mat)
dds = DESeqDataSetFromMatrix(countData = mat_rounded,
                             colData = coldata,
                             design = ~ scCnT)

# using LRT to find unique peaks
# source: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
dds_lrt = DESeq(dds, test = "LRT", reduced = ~ 1)
res_LRT = results(dds_lrt)
res_LRT = as.data.frame(res_LRT)

# keep regions below adj.p = 0.1
sig_res_LRT <- res_LRT %>%
  rownames_to_column(var = "locus") %>%
  as_tibble() %>%
  filter(padj < 0.05)

lrt_uniques = mat[sig_res_LRT$locus, ]


col_fun = colorRamp2(c(0, 100, 200), c("#f0f0f0", "#ffeda0", "#feb24c"))

# heatmap
png(
  file = glue("{result_folder}DESeq_LRT-unique_G4_heatmap.png"),
  width = 7,
  height = 7,
  units = 'in',
  res = 500
)
hm = Heatmap(
  lrt_uniques,
  name = "read density",
  clustering_distance_rows = "pearson",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.5),
  column_title = "DESeq2 LRT test, adj. p < 0.05",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  width = unit(8, "cm"),
  height = unit(14, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 90
)
hm
dev.off()

# make bed and export
lrt_uniques_bed = lrt_uniques %>% mutate(ranges = rownames(lrt_uniques)) %>%
  select(ranges) %>% separate(ranges, sep = "_", into = c("seq", "start", "end"))

write_tsv(
  lrt_uniques_bed,
  glue(
    "{result_folder}LRT_test-unique_G4-adjp0.1.bed"),
    col_names = FALSE
)
