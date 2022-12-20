suppressPackageStartupMessages({
  library("SingleCellExperiment")
  library("Matrix")
  library("Seurat")
  library("Signac")
  library("glue")
  library("data.table")
  library("tidyverse")
  library("cccd") 
})


cellranger_output = "../data/CellRanger/GFP_sorted/filtered_peak_bc_matrix/"
output = "../results/Seurat/"

# Marques et al. GSE75330
rna = read.table(
  "../data/GSE75330/GSE75330_Marques_et_al_mol_counts2.tab",
  stringsAsFactors = FALSE,
  header = FALSE
)
annot = readRDS(file = "../data/GSE75330/Marques2016annotation.rds")

colnames(rna) = gsub("-", "_", sub('-', '', gsub("C1-", "", rna[1,])))
rna = rna[-1,]
genes = rna[, 1]
rna = rna[,-1]
rna = as.matrix(rna)
rna = apply(rna, 2, as.numeric)
rownames(rna) = genes

rna = CreateSeuratObject(counts = rna,
                         meta.data = annot,
                         assay = 'RNA')

# PCA analysis
all.genes = rownames(rna)
rna = NormalizeData(rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna = ScaleData(rna, features = all.genes)
rna = RunPCA(rna, features = VariableFeatures(object = rna))


peaks = fread(glue("{cellranger_output}peaks.bed"))
peaks$V4 = sapply(seq(1, length(peaks$V1)), helper)
write_tsv(peaks, "../data/CellRanger/GFP_sorted/filtered_peak_bc_matrix/peaks_labeled.bed", col_names = FALSE)


helper = function(number) {
  y = paste0("peak_", number)
  return(y)
}

# Signac TF-IDF and SVD
process_signac <- function(data_path, ndims = 50) {
  sce <- create_sce(data_path)
  dat <- as.Seurat(sce, counts = "counts", data = "counts")
  dat <- RunTFIDF(dat)
  dat <- FindTopFeatures(dat, min.cutoff = "q0")
  dat <- RunSVD(dat, n = ndims)
  sce <- as.SingleCellExperiment(dat)
  return(sce)
}

create_sce <- function(path) {
  expression_matrix <- ReadMtx(mtx = glue("{cellranger_output}matrix.mtx"),
                               features = glue("{cellranger_output}peaks_labeled.bed"),
                               cells = glue("{cellranger_output}barcodes.tsv"), feature.column = 4)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  return(as.SingleCellExperiment(seurat_object))
}

create_sce(cellranger_output)

# run Signac process and LSI
signac = process_signac(cellranger_output, ndim=10)
write.csv(reducedDim(signac, "LSI"),
          glue("{output}GFP_sorted.csv"))


# read scRNA-Seq data and create Seurat object
annot = readRDS(file = "../data/GSE75330/Marques2016annotation.rds")

colnames(rna) = gsub("-", "_", sub('-', '', gsub("C1-", "", rna[1,])))
rna = rna[-1,]
genes = rna[, 1]
rna = rna[,-1]
rna = as.matrix(rna)
rna = apply(rna, 2, as.numeric)
rownames(rna) = genes

rna = CreateSeuratObject(counts = rna,
                         meta.data = annot,
                         assay = 'RNA')

# Seurat PCA on scRNA-Seq + knn graph
all.genes = rownames(rna)
rna = NormalizeData(rna,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
rna = FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna = ScaleData(rna, features = all.genes)
rna = RunPCA(rna, features = VariableFeatures(object = rna)) # PCA
rna_pca = rna@reductions$pca@feature.loadings

generate_k = function(cellnumber_percentage) {
  cell_number = rna@assays$RNA@counts@Dim[2] # for k value 
  k = cell_number * cellnumber_percentage
  return(k)
}
k = generate_k(0.001)
k
rna_nng = nng(rna_pca, k = k) # knn (package cccd)

# Signac LSI + knn graph
g4_lsi = signac@int_colData@listData$reducedDims@listData$LSI
g4_nng = nng(g4_lsi, k = k)

#
both = sum(as.matrix(rna_nng) & as.matrix(g4_nng))
score = both / (dim(g4_lsi)[1] * k)


dim(as.matrix(rna_nng))
dim(as.matrix(g4_nng))
