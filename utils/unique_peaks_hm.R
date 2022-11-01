suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ComplexHeatmap")
  library("circlize")
})

# result / peak folder
result_folder = "../results/Seurat/callpeaks_unsorted/"

uniques = fread("../results/Seurat/callpeaks_unsorted/Grubbs_test-unique_G4_peaks_0.01_joined.tsv")

rows = uniques %>% filter(!`Gene Name` %in% c("Tecrl", "Letm1", "Slit2", "Rab3d")) %>% pull(`Gene Name`)
uniques = uniques %>% filter(!`Gene Name` %in% c("Tecrl", "Letm1", "Slit2", "Rab3d")) %>% dplyr::select(-`Gene Name`, -`Distance to TSS`, -unique)
mat = as.matrix(uniques)
rownames(mat) = rows

col_fun = colorRamp2(c(0, 100, 200), c("#f0f0f0", "#ffeda0", "#feb24c"))
# ha = HeatmapAnnotation(
#   sex = c("401", "975", "980", "983a", "LT2e", "WA14"),
#   col = list(sex = c("LT2e" = "#636363", "980" = "#636363", "WA14" = "#f0f0f0", "983a" = "#f0f0f0", "975" = "#636363", "401" = "#f0f0f0")
#   ),
#   gp = gpar(col = "black"),
#   show_legend = FALSE
# )

png(
  file = glue("{result_folder}Grubbs_test-unique_G4_heatmap.png"),
  width = 5,
  height = 7,
  units = 'in',
  res = 500
)
top_hm = Heatmap(
  mat,
  name = "read density",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "closest gene symbol",
  row_title_side = "right",
  #rect_gp = gpar(col = "black", lwd = 0.0),
  # top_annotation = ha,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  width = unit(5, "cm"),
  height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 13),
  column_names_rot = 90
)

#lgd = Legend(at = c("XX", "XY"), title = "sex", legend_gp = gpar(fill = c("#636363", "#f0f0f0")))
#top_hm = draw(top_hm, heatmap_legend_list = lgd)
top_hm
dev.off()

pdf(
  file = glue("{result_folder}Grubbs_test-unique_G4_heatmap.pdf"),
  width = 5,
  height = 7
)
top_hm = Heatmap(
  mat,
  name = "read density",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "closest gene symbol",
  row_title_side = "right",
  #rect_gp = gpar(col = "black", lwd = 0.0),
  # top_annotation = ha,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  width = unit(5, "cm"),
  height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 2),
  column_names_gp = gpar(fontsize = 13),
  column_names_rot = 90
)

#lgd = Legend(at = c("XX", "XY"), title = "sex", legend_gp = gpar(fill = c("#636363", "#f0f0f0")))
#top_hm = draw(top_hm, heatmap_legend_list = lgd)
top_hm
dev.off()
