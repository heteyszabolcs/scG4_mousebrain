# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("topGO")
  library("EnsDb.Mmusculus.v79")
  library("ComplexHeatmap")
  library("circlize")
  library("enrichR")
})

### run on Linux (e.g. Uppmax). Windows returns a mysterious bug for topGO...
# folders
result_folder = "../results/Seurat/callpeaks_unsorted/"

# source annotation script for mm10
source("/crex/proj/snic2020-6-3/SZABOLCS/scripts/annotation.R")

# create topGO input. Consider promotor proximal peaks + peaks higher than median signalValue (MACS2)
# input: output of HOMER
create_input = function(peaks) {
  annotations = mm10_annotation(regions = peaks, seqname_col = "V1", start_col = "V2", end_col = "V3", feature_1 = "V5", feature_2 = "V5", feature_3 = "V5")
  
  genes = annotations %>% dplyr::filter(abs(`distanceToTSS`) < 3000) %>%
    dplyr::filter(feature_1 > median(feature_1)) %>% 
    pull(SYMBOL) 

  input = rep(0, dim(annotations)[1])
  names(input) = annotations$SYMBOL
  input[names(input) %in% genes] = 1
  
  return(input)
}

## topGO analysis
# create topGO output
# with top 10 enriched GO terms
# reference: "org.Mm.eg.db"
create_go_matrix = function(genes, colname) {
  # find biological process ontology
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    # use biological process ontology
    allGenes = genes,
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
             topNodes = 10,
             numChar = 60)
  
  out = out %>% dplyr::select(Term, Fisher)
  colnames(out) = c("Term", colname)
  
  return(out)
  
}

# unsorted brain scCutnTag Seurat clusters (res 0.1)
peaks_1 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
peaks_0 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")

# cluster res 0.1
cl0 = create_go_matrix(
  create_input(peaks_0),
  colname = "cluster 0"
)

cl1 = create_go_matrix(
  create_input(peaks_1),
  colname = "cluster 1"
)

go_outputs = list(cl0, cl1)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
full = full %>% distinct(Term, .keep_all = TRUE) %>% 
  mutate(`cluster 0` = str_replace(`cluster 0`, "< ", ""), `cluster 1` = str_replace(`cluster 1`, "< ", ""))
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)

# visualize enrichments by heatmap
generate_heatmap = function(matrix) {
  col_fun = colorRamp2(c(1e-10, 1e-12, 1e-16, 1e-20, 1e-30), c("grey", "#ece7f2", "#e5f5f9", "#99d8c9", "#2ca25f"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05, -0.001),
      labels = c("not enriched", "p < 1e-20", "p < 1e-10"),
      legend_height = unit(2, "cm")
    ),
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    heatmap_width = unit(7, "cm"),
    heatmap_height = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90
  )
  hm
}

pdf(
  file = glue("{result_folder}topGO_clusterwise_res0.1.pdf"),
  width = 5,
  height = 5
)
print(generate_heatmap(full))
dev.off()

png(
  file = glue("{result_folder}topGO_clusterwise_res0.1.png"),
  width = 11,
  height = 11,
  unit = "cm",
  res = 500
)
print(generate_heatmap(full))
dev.off()

# unsorted brain scCutnTag Seurat clusters (res 0.1)
peaks_0 = fread("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_peaks/0_peaks.narrowPeak")
peaks_1 = fread("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_peaks/1_peaks.narrowPeak")
peaks_2 = fread("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_peaks/2_peaks.narrowPeak")
peaks_3 = fread("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_peaks/3_peaks.narrowPeak")
peaks_4 = fread("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_peaks/4_peaks.narrowPeak")


cl0 = create_go_matrix(
  create_input(peaks_0),
  colname = "cluster 0"
)
cl1 = create_go_matrix(
  create_input(peaks_1),
  colname = "cluster 1"
)
cl2 = create_go_matrix(
  create_input(peaks_2),
  colname = "cluster 2"
)
cl3 = create_go_matrix(
  create_input(peaks_3),
  colname = "cluster 3"
)
cl4 = create_go_matrix(
  create_input(peaks_4),
  colname = "cluster 4"
)

go_outputs = list(cl0, cl1, cl2, cl3, cl4)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)
full = full[which(rownames(full) != "biological_process"),]

generate_heatmap = function(matrix) {
  col_fun = colorRamp2(c(0.05, 1e-10, 1e-14, 1e-28, 1e-22), c("grey", "#ece7f2", "#e5f5f9", "#99d8c9", "#2ca25f"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05, -0.001),
      labels = c("not enriched", "p < 1e-20", "p < 0.05"),
      legend_height = unit(2, "cm")
    ),
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    heatmap_width = unit(7, "cm"),
    heatmap_height = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90
  )
  hm
}

pdf(
  file = glue("{result_folder}topGO_clusterwise.pdf"),
  width = 5,
  height = 5
)
print(generate_heatmap(full))
dev.off()

png(
  file = glue("{result_folder}topGO_clusterwise.png"),
  width = 11,
  height = 11,
  unit = "cm",
  res = 500
)
print(generate_heatmap(full))
dev.off()

# visualize enrichments by heatmap
generate_heatmap = function(matrix) {
  col_fun = colorRamp2(c(0.05,  0.025, 0.010, 0.005, 0.001), c("grey", "#ece7f2", "#e5f5f9", "#99d8c9", "#2ca25f"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05, -0.001),
      labels = c("not enriched", "p < 0.002", "p < 0.05"),
      legend_height = unit(2, "cm")
    ),
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    heatmap_width = unit(7, "cm"),
    heatmap_height = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90
  )
  hm
}

# FindAllMarkers output
cl1 = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_logreg_res0.1-cluster1.tsv")
cl0 = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_logreg_res0.1-cluster0.tsv")

peaks_1 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
peaks_0 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")

create_input_fam = function(fc_table, background = peaks_1) {
  annotations = mm10_annotation(regions = background, seqname_col = "V1", start_col = "V2", end_col = "V3", feature_1 = "V5", feature_2 = "V5", feature_3 = "V5")
  background = annotations %>% dplyr::filter(abs(`distanceToTSS`) < 3000) %>%
    dplyr::filter(feature_1 > median(feature_1)) %>% 
    pull(SYMBOL) 
  
  background = annotations %>% dplyr::filter(abs(`distanceToTSS`) < 3000) %>%
    dplyr::filter(feature_1 > median(feature_1)) %>% 
    pull(SYMBOL) 
  
  genes = fc_table %>% pull(gene_symbol) 
  
  input = rep(0, length(background))
  names(input) = background
  input[names(input) %in% genes] = 1
  
  return(input)
}

cl1 = create_go_matrix(
  create_input_fam(fc_table = cl1, background = peaks_1),
  colname = "cluster 1"
)
cl0 = create_go_matrix(
  create_input_fam(fc_table = cl0, background = peaks_0),
  colname = "cluster 0"
)

go_outputs = list(cl0, cl1)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)
full = full[which(rownames(full) != "biological_process"),]

pdf(
  file = glue("{result_folder}topGO_findmarkers_res0.1.pdf"),
  width = 5,
  height = 5
)
print(generate_heatmap(full))
dev.off()

png(
  file = glue("{result_folder}topGO_findmarkers_res0.1.png"),
  width = 11,
  height = 11,
  unit = "cm",
  res = 500
)
print(generate_heatmap(full))
dev.off()
