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
})

# folders
result_folder = "../results/Seurat/callpeaks_unsorted/"

# create topGO input. Consider promotor proximal peaks + peaks higher than median signalValue (MACS2)
# input: output of HOMER
create_input = function(path_to_peaks) {
  annotations = fread(path_to_peaks)
  
  genes = annotations %>% dplyr::filter(abs(`Distance to TSS`) < 3000) %>%
    dplyr::filter(signalValue > median(signalValue)) %>% 
    pull(`Gene Name`) 

  input = rep(0, dim(annotations)[1])
  names(input) = annotations$`Gene Name`
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

# unsorted brain scCutnTag Seurat clusters

# cluster res 0.1
cl0 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/0_unique_peaks_res0.1_annot.tsv"
  ),
  colname = "cluster 0"
)

cl1 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/1_unique_peaks_res0.1_annot.tsv"
  ),
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
  col_fun = colorRamp2(c(1, 0.0001, 0.000001), c("#bdbdbd", "#ffeda0", "#de2d26"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05, -0.001),
      labels = c("not enriched", "1e-6", "1e-3"),
      legend_height = unit(1, "cm")
    ),
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    heatmap_width = unit(7, "cm"),
    heatmap_height = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 8),
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

# cluster res default
cl0 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/0_peaks_annot.tsv"
  ),
  colname = "cluster 0"
)
cl1 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/1_peaks_annot.tsv"
  ),
  colname = "cluster 1"
)
cl2 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/2_peaks_annot.tsv"
  ),
  colname = "cluster 2"
)
cl3 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/3_peaks_annot.tsv"
  ),
  colname = "cluster 3"
)
cl4 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/4_peaks_annot.tsv"
  ),
  colname = "cluster 4"
)
cl5 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/5_peaks_annot.tsv"
  ),
  colname = "cluster 5"
)
cl6 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/peak_sets/6_peaks_annot.tsv"
  ),
  colname = "cluster 6"
)

go_outputs = list(cl0, cl1, cl2, cl3, cl4, cl5, cl6)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)

# visualize enrichments by heatmap
generate_heatmap = function(matrix) {
  col_fun = colorRamp2(c(1, 0.05, 0.001), c("#bdbdbd", "#ffeda0", "#de2d26"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05, -0.001),
      labels = c("not enriched", "0.001", "0.05"),
      legend_height = unit(1, "cm")
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

# visualize enriched GO terms in FindAllMarkers output
create_input = function(path_to_markers) {
  annotations = fread(path_to_markers)
  
  genes = annotations %>% dplyr::filter(abs(`Distance to TSS`) < 3000) %>%
    dplyr::filter(avg_log2FC > 0.2 & p_val_adj < 0.05) %>% 
    pull(`Gene Name`) 
  
  input = rep(0, dim(annotations)[1])
  names(input) = annotations$`Gene Name`
  input[names(input) %in% genes] = 1
  
  return(input)
}

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

cl0 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/FindAllMarker_cluster_0_promoters.tsv"
  ),
  colname = "cluster 0"
)

cl1 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/FindAllMarker_cluster_1_promoters.tsv"
  ),
  colname = "cluster 1"
)

cl2 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/FindAllMarker_cluster_2_promoters.tsv"
  ),
  colname = "cluster 2"
)

cl3 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/FindAllMarker_cluster_3_promoters.tsv"
  ),
  colname = "cluster 3"
)

cl4 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/FindAllMarker_cluster_4_promoters.tsv"
  ),
  colname = "cluster 4"
)

cl5 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/FindAllMarker_cluster_5_promoters.tsv"
  ),
  colname = "cluster 5"
)

cl6 = create_go_matrix(
  create_input(
    "../results/Seurat/callpeaks_unsorted/FindAllMarker_cluster_6_promoters.tsv"
  ),
  colname = "cluster 6"
)

go_outputs = list(cl0, cl1, cl2, cl3, cl4, cl5, cl6)

full <- Reduce(function(x, y, ...)
  merge(x, y, all = TRUE, ...),
  go_outputs)

full[is.na(full)] = 1
full = full %>% distinct(Term, .keep_all = TRUE) %>% 
  mutate(`cluster 0` = str_replace(`cluster 0`, "< ", ""), `cluster 1` = str_replace(`cluster 1`, "< ", ""))
rownames(full) = full$Term
full = full %>% dplyr::select(-Term) %>% mutate_if(is.character, as.numeric)
full = as.matrix(full)
full[is.na(full)] = 0

full = full[-which(rowMeans(full) == 1, arr.ind = FALSE),]


generate_heatmap = function(matrix) {
  col_fun = colorRamp2(c(0.05, 1), c("#de2d26", "#bdbdbd"))
  hm = Heatmap(
    matrix,
    column_title = "",
    row_title = "",
    name = "Fisher test",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Fisher test",
      at = c(1, -0.05),
      labels = c("not enriched", "p < 0.05"),
      legend_height = unit(1, "cm")
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
  file = glue("{result_folder}topGO_clusterwise_G4markers.pdf"),
  width = 5,
  height = 5
)
print(generate_heatmap(full))
dev.off()

