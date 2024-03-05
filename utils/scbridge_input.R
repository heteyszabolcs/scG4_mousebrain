if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("Seurat",
               "glue",
               "tidyverse",
               "data.table"
)

# result folder
result_folder = "../results/scBridge/"

# scBridge needs:
# A: gene activity scores of scATAC/scCut&Tag
# B: normalized scRNA-Seq counts
# keep those genes that overlap between A and B. 

# scRNA-Seq input
counts = fread("../data/GSE75330/GSE75330_Marques_et_al_mol_counts2.tab")
rna = 
  readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/scRNA_Seq_Marques_et_el.Rds")

mat = counts[,2:ncol(counts)]
mat = as.data.frame(mat)
rownames(mat) = counts$cellid

meta = rna@meta.data
as.character(tail(meta$cellid))
colnames(mat) = meta$cellid

# G4 input
g4 = 
  readRDS("../results/Seurat/final/sorted_brain/res0.8/outputs/Seurat_object.Rds")
ga = g4@assays$GA$counts
ga = as.data.frame(ga)

# intersect
int = intersect(rownames(ga), rownames(mat))
mat = mat[int,]
ga = ga[int,]

# export
write.csv(mat, file = glue("{result_folder}GSE75330_Marques_et_al_mol_counts2.csv"), quote = FALSE)
annot = tibble(cell_id = meta$cellid, CellType = meta$cell_class)
write_csv(annot, file = glue("{result_folder}GSE75330_Marques_et_al_annot.csv"))

write.csv(ga, file = glue("{result_folder}gene_activity_scores.csv"), quote = FALSE)


