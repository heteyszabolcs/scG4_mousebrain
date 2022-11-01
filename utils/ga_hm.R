# packages
suppressPackageStartupMessages({
  library("data.table")
  library("ComplexHeatmap")
  library("glue")
  library("tidyverse")
  library("ComplexHeatmap")
  library("circlize")
  library("matrixStats")
})

# result folder
result_folder = "../results/Seurat/callpeaks_unsorted/"

# create gene activity object
unsorted = readRDS(glue("{result_folder}unsorted_int_Marques.Rds"))
ga = unsorted@assays$GA@data
ident = unsorted@active.ident

# retrieve genes
get_genes = function(cell_type) {
  markers = fread("../data/GSE75330/marker_genes.txt", header = FALSE)
  markers = markers %>% dplyr::filter(V2 == cell_type) %>% pull(V1)
  markers = intersect(markers, rownames(unsorted@assays$GA@data))
  return(markers)
}

# retrieve identifications
get_ga = function(cluster) {
  idents = names(ident[which(ident == cluster)])
  return(idents)
}

# gene activity heatmap
ga_heatmap = function(cell_type) {
  clusters = c("0", "1", "2", "3", "4", "5", "6")
  mat = lapply(clusters, function(x)
    colMeans(t(as.matrix(ga[get_genes(cell_type), get_ga(x)]))))
  mat = do.call(rbind, mat)
  mat = t(mat)
  mat = mat[order(rowVars(mat), decreasing = TRUE),]
  
  if (!nrow(mat) < 10) {
    mat = mat[1:10,]
  }
  
  colnames(mat) = c("0", "1", "2", "3", "4", "5", "6")
  mat = as.matrix(mat)
  col_fun = colorRamp2(c(0, 0.05, 0.1), c("#9ecae1", "#deebf7", "#fc9272"))
  
  hm = Heatmap(
    mat,
    column_title = " ",
    row_title = paste0(cell_type, " marker"),
    name = "G4 gene act.",
    col = col_fun,
    # heatmap_legend_param = list(
    #   at = c(0, 0.5, 1),
    #   labels = c("< 0.25", "0.5", "1"),
    #   legend_height = unit(2, "cm")
    # ),
    rect_gp = gpar(col = "black", lwd = 1),
    show_column_dend = FALSE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    heatmap_width = unit(8, "cm"),
    heatmap_height = unit(12, "cm"),
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 0
  )
  return(print(hm))
}

unique(markers$V2)
markers = fread("../data/GSE75330/marker_genes.txt", header = FALSE)
lapply(unique(markers$V2), ga_heatmap)

pdf(
  file = glue("{result_folder}G4_gene_activity_OPC.pdf"),
  width = 5,
  height = 5
)
ga_heatmap("OPC")
dev.off()

pdf(
  file = glue("{result_folder}G4_gene_activity_newOL.pdf"),
  width = 5,
  height = 5
)
ga_heatmap("newly_formed_OL")
dev.off()







