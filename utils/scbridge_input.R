if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("Seurat", "glue", "tidyverse", "data.table")

# result folder
result_folder = "../results/scBridge/"

# scBridge needs:
# A: gene activity scores of scATAC/scCut&Tag
# B: normalized scRNA-Seq counts
# keep those genes that overlap between A and B.

# scRNA-Seq input - Marques et al.
# counts = fread("../data/GSE75330/GSE75330_Marques_et_al_mol_counts2.tab")
# rna =
#   readRDS(
#     "../results/Seurat/final/sorted_brain/res0.8/integration/outputs/scRNA_Seq_Marques_et_el.Rds"
#   )
# 
# mat = counts[, 2:ncol(counts)]
# mat = as.data.frame(mat)
# rownames(mat) = counts$cellid
# 
# meta = rna@meta.data
# as.character(tail(meta$cellid))
# colnames(mat) = meta$cellid

# scRNA-Seq input - Bartosovic et al.
rna_bartosovic = readRDS("../data/GSE163484/brain_object.Rds")
meta = rna_bartosovic@meta.data
rna_bartosovic = rna_bartosovic@assays$RNA@data
rna_bartosovic = as.data.frame(rna_bartosovic)

colnames(rna_bartosovic) = unname(sapply(colnames(rna_bartosovic), function(x) {
  strsplit(x, "_1")[[1]][1]
}))

# G4 input
g4 =
  readRDS("../results/Seurat/final/sorted_brain/res0.8/outputs/Seurat_object.Rds")
ga = g4@assays$GA$counts
ga = as.data.frame(ga)

# intersect
int = intersect(rownames(ga), rownames(rna_bartosovic))
rna_bartosovic = rna_bartosovic[int, ]
ga = ga[int, ]

# export
write.csv(
  rna_bartosovic,
  file = glue("{result_folder}GSE163484_Bartosovic_et_al_counts.csv"),
  quote = FALSE
)
cell_ids = unname(sapply(rownames(meta), function(x) {
  strsplit(x, "_1")[[1]][1]
}))
annot = tibble(cell_id = cell_ids, CellType = meta$cell_type)
write_csv(annot,
          file = glue("{result_folder}GSE163484_Bartosovic_et_al_annot.csv"))
write.csv(ga,
          file = glue("{result_folder}gene_activity_scores.csv"),
          quote = FALSE)
