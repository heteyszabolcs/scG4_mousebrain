# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("ggplot2")
  library("glue")
  library("tidyverse")
  library("data.table")
  library("RColorBrewer")
})

# path to result folder
result_folder = "../results/Seurat/callpeaks_GFPsorted/"

# cell type predictions coming from Seurat's label transferring (seurat_int_umap.R)
cell_class_pred = readRDS("../results/Seurat/final/sorted_brain/integration/outputs/g4_cell_label_preds.Rds")
types = rownames(cell_class_pred@data)
cell_class_pred = as_tibble(cell_class_pred@data)
cell_class_pred = cell_class_pred %>% mutate(type = types) %>% dplyr::select(type, everything())

# coembedding G4 and Marques et al. scRNA-Seq
# main source: 
# https://github.com/Castelo-Branco-lab/scCut-Tag_2020/blob/master/notebooks/integration/integration_H3K4me3_marques.Rmd
# Seurat objects
sorted = readRDS("../results/Seurat/final/sorted_brain/integration/outputs/G4_Marques_scRNA_integration.Rds")
rna = readRDS("../results/Seurat/final/sorted_brain/integration/outputs/scRNA_Seq_Marques_et_el.Rds")
genes.use = VariableFeatures(rna)

rna@meta.data$data_type = "scRNA-Seq (Marques et al.)"
sorted@meta.data$data_type = "G4 scCut&Tag"

coembed = merge(x = rna, y = sorted)

coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
ElbowPlot(coembed)
coembed = RunUMAP(coembed, dims = 1:15)

# merge together Marques et al. cell labels
new_ids = as.character(coembed@meta.data$cell_class)
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
coembed@meta.data$cell_class = new_ids

coembed.g4 = coembed[,coembed$data_type == "G4 scCut&Tag"]
coembed.scrna = coembed[,coembed$data_type == "scRNA-Seq (Marques et al.)"]

# marker gene statistics
markers = c("Ptprz1", "Gpr17", "Slc1a1", "Ctps", "Sepp1", "Igf2")
max_ga = sapply(markers, function(x) max(coembed.g4@assays$GA[x,]))
max_expr = sapply(markers, function(x) max(coembed.scrna@assays$RNA[x,]))
ga_thr = c(0.5, 0.5, 0.5, 0.0, 0.0, 0.0)
#quant_99 = sapply(markers, function(x) quantile(as.vector(coembed.g4@assays$GA[x, ]), 0.9999))

marker_stat = tibble(marker = markers, max_GA = max_ga, max_expression = max_expr, GA_thr = ga_thr)

# visualizations
vis_dim_plot = function(gene) {
  print(gene)
  ga_thr = unname(marker_stat[which(marker_stat$marker == gene),]$GA_thr)
  
  get_cell_ids = coembed.g4@assays$GA@data[gene,]
  get_cell_ids = names(get_cell_ids[get_cell_ids > ga_thr])
  
  predictions = cell_class_pred %>%
    pivot_longer(
      .,
      cols = c(colnames(cell_class_pred)[2]:colnames(cell_class_pred)[dim(cell_class_pred)[2]]),
      names_to = "cell_id",
      values_to = "pred_score"
    ) %>%
    dplyr::filter(!type == "max")
  
  plot = DimPlot(
    object = coembed.g4,
    cells.highlight = get_cell_ids,
    cols.highlight = "red",
    cols = "gray",
    #order = TRUE,
    #raster = TRUE
  ) + NoAxes() + ggtitle(gene) + NoLegend()
  return(print(plot))
}

dims = lapply(markers, vis_dim_plot)
dims[[5]]

FeaturePlot(
  object = coembed.g4,
  features = "Sepp1",
  min.cutoff = min(coembed.g4@assays$GA@data["Sepp1",]),
  max.cutoff = max(coembed.g4@assays$GA@data["Sepp1",]),
  pt.size = 2,
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-10, 10) +
  scale_color_viridis_c() +
  theme(
    legend.position = 'bottom',
    text = element_text(size = 15),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

for(i in seq(length(marker_stat$marker))) {
  ggsave(
    glue(
      "{result_folder}coenrich_cluster_anal_dimplot_{marker_stat$marker[i]}.pdf"
    ),
    plot = dims[i][[1]],
    width = 4,
    height = 3,
    device = "pdf"
  )
}

vis_pred_boxplot = function(gene) {
  print(gene)
  
  get_cell_ids = coembed.g4@assays$GA@data[gene,]
  get_cell_ids = names(get_cell_ids[get_cell_ids > ga_thr])
  
  predictions = cell_class_pred %>%
    pivot_longer(
      .,
      cols = c(colnames(cell_class_pred)[2]:colnames(cell_class_pred)[dim(cell_class_pred)[2]]),
      names_to = "cell_id",
      values_to = "pred_score"
    ) %>%
    dplyr::filter(!type == "max")
  
  predictions = predictions %>% dplyr::filter(cell_id %in% get_cell_ids)
  
  plot = ggplot(predictions,
                aes(x = type, y = pred_score, fill = type)) +
    geom_boxplot(color = "black") +
    scale_fill_brewer(palette = "Set3") +
    ylim(0, 1) +
    labs(
      title = "",
      x = "",
      y = "prediction score",
      fill = ""
    ) +
    theme_classic() +
    ggtitle(gene) +
    guides(fill = "none") +
    theme(
      text = element_text(size = 9),
      plot.title = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 20, color = "black")
    )
  
  return(print(plot))
}

boxplots = lapply(markers, vis_pred_boxplot)
boxplots[1]

for(i in seq(length(marker_stat$marker))) {
  ggsave(
    glue(
      "{result_folder}coenrich_cluster_anal_pred_box_{marker_stat$marker[i]}.pdf"
    ),
    plot = boxplots[i][[1]],
    width = 6,
    height = 4,
    device = "pdf"
  )
}

# sorted_chroma = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")
# 
# x = sorted@meta.data
# x = x %>% rownames_to_column(var = "cell_id")
# x = x %>% dplyr::select(cell_id, pred_cell_type)
# 
# y = sorted_chroma@meta.data
# y = y %>% rownames_to_column(var = "cell_id")
# 
# z = inner_join(x, y, by = c("cell_id" = "cell_id"))
# 
# sorted_chroma@meta.data = z
# 
# TilePlot(object = sorted_chroma, region = "Lin28a", 
#          group.by = "pred_cell_type", tile.cells = 100, order.by = "total") + 
#   scale_fill_gradient(low = "white", high = "black")





