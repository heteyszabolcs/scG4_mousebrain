# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("tidyverse")
  library("ggplot2")
  library("GenomicFeatures")
})

set.seed(5)

# path to result folder
result_folder = "../results/Seurat/"

# read CellRanger output files (h5 formats)
mes_mef <- Read10X_h5(filename =
                        "../data/CellRanger/mES-mEF/filtered_peak_bc_matrix.h5")
mes_mef <-
  CreateSeuratObject(counts = mes_mef, project = "mESC-MEF")
mes_mef_cells = names(mes_mef@active.ident)
mes_mef_cells = paste0("mESC-MEF_", mes_mef_cells)

unsorted <- Read10X_h5(filename =
                         "../data/CellRanger/unsorted/filtered_peak_bc_matrix.h5")
unsorted <-
  CreateSeuratObject(counts = unsorted, project = "unsorted")
unsorted_cells = names(unsorted@active.ident)
unsorted_cells = paste0("brain_", unsorted_cells)

sorted <- Read10X_h5(filename =
                       "../data/CellRanger/GFP_sorted/filtered_peak_bc_matrix.h5")
sorted <- CreateSeuratObject(counts = sorted, project = "sorted")
sorted_cells = names(sorted@active.ident)
sorted_cells = paste0("oligo_", sorted_cells)

# combine Seurat objects
combined = merge(
  unsorted,
  y = mes_mef,
  add.cell.ids = c("unsorted brain", "mESC-MEF"),
  project = "combined"
)



# normalization and dim. reduction
all.genes = rownames(combined)
combined = ScaleData(combined, features = all.genes)
combined = FindVariableFeatures(object = combined)
combined = RunPCA(combined, features = VariableFeatures(object = combined))
combined = FindNeighbors(combined, dims = 1:10)
combined = FindClusters(combined, resolution = 0.5)
# run non-linear dimensional reduction
combined = RunUMAP(
  combined,
  dims = 1:10,
  n.neighbors = 250,
  n.components = 2
)

# dimplot
combined_dim =
  DimPlot(
    combined,
    label = FALSE,
    cells.highlight = mes_mef_cells,
    cols.highlight = "#fc9272",
    pt.size = 1.5,
  ) +
  #scale_colour_manual(values = cols, breaks = c("0", "4"), labels = c("GFP+", "unsorted")) +
  scale_colour_manual(
    labels = c("mESC-MEF", "unsorted brain"),
    values = c("#9ecae1", "#fc9272")
  ) +
  ggtitle("merged scCut&Tag datasets") +
  xlim(-25, 25) + ylim(-25, 25) +
  theme(
    text = element_text(size = 30),
    plot.title = element_text(25),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  )
combined_dim

ggsave(
  glue("{result_folder}merged_mESC-MEF_unsorted-UMAP.png"),
  plot = combined_dim,
  width = 7,
  height = 7,
  dpi = 500,
)

ggsave(
  glue("{result_folder}merged_mESC-MEF_unsorted-UMAP.pdf"),
  plot = combined_dim,
  width = 7,
  height = 7
)
