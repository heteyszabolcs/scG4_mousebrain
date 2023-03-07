# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("tidyverse")
  library("data.table")
  library("glue")
  library("devtools")
})

# path to result folder
result_path = "../results/Seurat/"

# function for Seurat integration
# source: https://satijalab.org/seurat/articles/integration_introduction.html
integration = function(seurat1 = "../data/GSE157637/Olig2_seurat_object.Rds",
                       seurat2 = "../data/merged/Seurat_merged.Rds",
                       label1 = "Olig2",
                       label2 = "G4",
                       assay = "bins_5000",
                       export_rds = TRUE,
                       export_image = TRUE) {
  print(glue("{label1} with {label2}"))
  
  seurat1 = readRDS(seurat1)
  
  if(assay == "GA") {
    DefaultAssay(seurat1) <- assay
  } else {
    DefaultAssay(seurat1) <- "RNA"
  }
  
  seurat1 = AddMetaData(object = seurat1,
                        metadata = label1,
                        col.name = "group")
  for (i in Assays(seurat1)) {
    if (i != assay) {
      seurat1[[i]] = NULL
    }
  }
  
  seurat2 = readRDS(seurat2)
  DefaultAssay(seurat2) <- "GA"
  seurat2 = AddMetaData(object = seurat2,
                        metadata = label2,
                        col.name = "group")
  for (i in Assays(seurat2)) {
    if (i != "GA") {
      seurat2[[i]] = NULL
    }
  }
  
  comb = merge(
    seurat1,
    y = seurat2,
    add.cell.ids = c(label1, label2),
    project = "int"
  )
  
  rm(seurat1)
  rm(seurat2)
  
  # split combined object keeping atomic feature
  split = SplitObject(comb, split.by = "group")
  
  # normalize and identify variable features for each dataset independently
  split <- lapply(
    X = split,
    FUN = function(x) {
      x <- NormalizeData(x)
      x <-
        FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    }
  )
  
  # select features that are repeatedly variable across datasets for integration
  features = SelectIntegrationFeatures(object.list = split)
  # find anchors
  anchors = FindIntegrationAnchors(object.list = split, anchor.features = features)
  # this command creates an 'integrated' data assay
  integrated = IntegrateData(anchorset = anchors)
  # export rds
  if (export_rds) {
    saveRDS(integrated,
            glue("{result_path}{label1}_{label2}_integrated.rds"))
  }
  
  # integration based on Stuart et al. (Seurat protocol)
  # run the standard workflow for visualization and clustering
  standard_wf <- ScaleData(integrated, verbose = FALSE)
  standard_wf <- RunPCA(standard_wf, npcs = 30, verbose = FALSE)
  standard_wf <-
    RunUMAP(standard_wf, reduction = "pca", dims = 1:30)
  standard_wf <-
    FindNeighbors(standard_wf, reduction = "pca", dims = 1:30)
  standard_wf <- FindClusters(standard_wf, resolution = 0.5)
  
  #exporting
  saveRDS(standard_wf, glue("{result_path}{label1}__{label2}_integrated.rds"))
  
  # Visualization
  p1 <- DimPlot(
    standard_wf,
    reduction = "umap",
    pt.size = 1,
    label.size = 7,
    group.by = 'group',
    repel = TRUE,
    na.value = "grey50",
    raster = TRUE
  ) +
    scale_fill_manual(labels = c(label2, glue("Bartosovic et al. {label1}"), values = c("#fc9272", "#9ecae1"))) +
    xlim(-15, 15) +
    ylim(-15, 15) +
    ggtitle(glue("scC&T integration, {label1} - {label2}")) +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    )
  
  if (export_image) {
    ggsave(
      glue("{result_path}Seurat_int-{label1}_{label2}.png"),
      plot = last_plot(),
      width = 10,
      height = 10,
      dpi = 300,
    )
    ggsave(
      glue("{result_path}Seurat_int-{label1}_{label2}.pdf"),
      plot = last_plot(),
      width = 10,
      height = 7,
      device = "pdf"
    )
  }
  
  return(print(p1))
}

# run
# integration(seurat1 = "../data/GSE157637/H3K27ac_seurat_object.Rds",
#             label1 = "H3K27ac",
#             assay = "GA")
# integration(seurat1 = "../data/GSE157637/H3K4me3_seurat_object.Rds",
#             label1 = "H3K4me3",
#             assay = "GA")
# integration(seurat1 = "../data/GSE157637/H3K27me3_seurat_object.Rds",
#             label1 = "H3K27me3",
#             assay = "GA")
# integration(seurat1 = "../data/GSE157637/Rad21_seurat_object.Rds",
#             label1 = "Rad21",
#             assay = "GA")
# integration(seurat1 = "../data/GSE157637/Olig2_seurat_object.Rds",
#             label1 = "Olig2",
#             assay = "GA")
# integration(seurat1 = "../data/GSE157637/H3K36me3_seurat_object.Rds",
#             label1 = "H3K36me3",
#             assay = "GA")

# integration of GFP sorted G4 data with Bartosovic et al. H3K4me3
integration(seurat1 = "../data/GSE157637/H3K4me3_seurat_object.Rds",
            seurat2 = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds",
            label1 = "H3K4me3",
            label2 = "G4",
            assay = "GA")

# H3K4me3 - G4 UMAPs 
int = readRDS("../results/Seurat/H3K4me3_G4_integrated.rds")

int <- ScaleData(int, verbose = FALSE)
int <- RunPCA(int, npcs = 30, verbose = FALSE)
int <-
  RunUMAP(int, reduction = "pca", dims = 1:30)
int <-
  FindNeighbors(int, reduction = "pca", dims = 1:30)
int <- FindClusters(int, resolution = 0.5)

DimPlot(
  int,
  reduction = "umap",
  pt.size = 1,
  label.size = 7,
  group.by = 'group',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_colour_manual(values = c("#fc9272", "#9ecae1")) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle(glue("scC&T integration, H3K4me3 - G4")) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  "../results/Seurat/callpeaks_GFPsorted/H3K4me3_G4_integrated_UMAP.pdf",
  plot = last_plot(),
  width = 10,
  height = 7,
  device = "pdf"
)

FeaturePlot(
  object = int,
  features = "marker_Oligodendrocytes",
  min.cutoff = min(int@meta.data$marker_Oligodendrocytes),
  max.cutoff = max(int@meta.data$marker_Oligodendrocytes),
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_color_viridis_c() +
  ggtitle("scC&T integration, H3K4me3 - G4", subtitle = "Oligodendrocyte markers") +
  theme(
    legend.position = 'bottom',
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  "../results/Seurat/callpeaks_GFPsorted/H3K4me3_G4_int_UMAP_oligo_mark.pdf",
  plot = last_plot(),
  width = 10,
  height = 7,
  device = "pdf"
)

g4_meta = readRDS("../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")
g4_preds = readRDS("../results/Seurat/callpeaks_GFPsorted/sorted_cell_label_preds.Rds")

g4_meta = g4_meta@meta.data
g4_ids = rownames(g4_meta)
g4_meta = g4_meta %>% mutate(cell_id = paste0("G4_", g4_ids))

int_meta = int@meta.data
int_ids = rownames(int_meta)
int_meta = int_meta %>% mutate(cell_id = int_ids)

comb_meta = inner_join(int_meta, g4_meta, by = "cell_id") 
comb_meta = comb_meta %>%
  rename(G4_seurat_clusters = seurat_clusters.y) %>% 
  rename(int_seurat_clusters = seurat_clusters.x) 
comb_meta = as.data.frame(comb_meta)
rownames(comb_meta) = comb_meta$cell_id
  
int@meta.data = comb_meta

DimPlot(
  int,
  reduction = "umap",
  pt.size = 1,
  label.size = 7,
  group.by = 'G4_seurat_clusters',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle("scC&T integration, H3K4me3 - G4", subtitle = "GFP+ G4 Seurat clusters") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  "../results/Seurat/callpeaks_GFPsorted/H3K4me3_G4_int_UMAP_G4_clusters.pdf",
  plot = last_plot(),
  width = 10,
  height = 7,
  device = "pdf"
)

g4_preds_meta = g4_preds@data
cell_type = rownames(g4_preds_meta)
g4_preds_meta = as_tibble(g4_preds_meta)
g4_preds_meta = g4_preds_meta %>% mutate(cell_type = cell_type)
g4_preds_meta = pivot_longer(g4_preds_meta, cols = "AAACGAAAGAAGCCGT-1":"TTTGTGTTCTCGCGTT-1",
                            names_to = "cell_id", values_to = "prediction_score")
g4_preds_meta = g4_preds_meta %>% mutate(cell_id = paste0("G4_", cell_id)) %>% 
  group_by(cell_id) %>% 
  filter(prediction_score == max(prediction_score)) %>% filter(!cell_type == "max")

predscore_hist = ggplot(g4_preds_meta, aes(x = prediction_score)) +
  geom_histogram(position = "identity", alpha = 1, aes(fill = after_stat(x > 0.75))) +
  scale_fill_manual(values = c("#f0f0f0", "#fc9272")) +
  labs(title = "Prediction score distributions",
       x = "prediction score",
       y = "") +
  xlim(0, 1) +
  theme_classic() +
  theme(text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black"))
predscore_hist
    
high_preds = g4_preds_meta %>% filter(prediction_score >= 0.75)
high_preds_cellids = high_preds %>% pull(cell_id)
high_preds_mol = high_preds %>% filter(cell_type == "MOL") %>% 
  mutate(barcode = str_split(cell_id, "_")[[1]][2])
write_tsv(high_preds_mol, "../results/Seurat/callpeaks_GFPsorted/high_pred_MOL.tsv")

int_meta = int@meta.data
comb_meta2 = inner_join(int_meta, g4_preds_meta, by = "cell_id") 
comb_meta2 = comb_meta2 %>% distinct(cell_id, .keep_all = TRUE) %>% 
  rename(prediction_score = prediction_score) %>% 
  rename(pred_cell_type = cell_type.y)
comb_meta2 = as.data.frame(comb_meta2)
rownames(comb_meta2) = comb_meta2$cell_id
int@meta.data = comb_meta2


FeaturePlot(
  object = int,
  features = "prediction_score",
  min.cutoff = min(int@meta.data$prediction_score),
  max.cutoff = max(int@meta.data$prediction_score),
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle("scC&T integration, H3K4me3 - G4", subtitle = "Prediction scores") +
  scale_color_viridis_c() +
  theme(
    legend.position = 'bottom',
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  "../results/Seurat/callpeaks_GFPsorted/H3K4me3_G4_int_UMAP_preds.pdf",
  plot = last_plot(),
  width = 10,
  height = 7,
  device = "pdf"
)

DimPlot(
  int,
  reduction = "umap",
  pt.size = 1,
  label.size = 7,
  group.by = 'pred_cell_type',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle("scC&T integration, H3K4me3 - G4", subtitle = "Predicted cell types") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  "../results/Seurat/callpeaks_GFPsorted/H3K4me3_G4_int_UMAP_pred_celltypes.pdf",
  plot = last_plot(),
  width = 10,
  height = 7,
  device = "pdf"
)

comb_meta2_filt = comb_meta2 %>% filter(cell_id %in% high_preds_cellids)
int@meta.data = comb_meta2_filt

DimPlot(
  int,
  reduction = "umap",
  pt.size = 1,
  label.size = 7,
  group.by = 'pred_cell_type',
  repel = TRUE,
  na.value = "grey50",
  raster = TRUE
) +
  scale_color_brewer(palette = "Set3") +
  xlim(-15, 15) +
  ylim(-15, 15) +
  ggtitle("scC&T integration, H3K4me3 - G4", subtitle = "pred. score >= 0.75") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes()

ggsave(
  "../results/Seurat/callpeaks_GFPsorted/H3K4me3_G4_int_UMAP_highly_pred_celltypes.pdf",
  plot = last_plot(),
  width = 10,
  height = 7,
  device = "pdf"
)



