library("ggplot2")
library("glue")
library("GenomicRanges")
library("Seurat")

H3K4me3 <-
  readRDS(file = "../data/GSE157637/H3K4me3_seurat_object.Rds")
g4 <- readRDS(file = "../data/merged/Seurat_merged.Rds")

p1 <- DimPlot(H3K4me3, label = TRUE) + NoLegend()
p2 <- DimPlot(g4, label = TRUE) + NoLegend()
p1 + p2

DefaultAssay(H3K4me3) <- "GA"
DefaultAssay(g4) <- "GA"

common.genes <- intersect(rownames(H3K4me3), rownames(g4))

g4[["bins_25000"]] = NULL
g4[["peaks"]] = NULL
g4[["PA"]] = NULL
g4[["bins_5000"]] = NULL

H3K4me3[["peaks"]] = NULL
H3K4me3[["PA"]] = NULL
H3K4me3[["bins_5000"]] = NULL

H3K4me3 = subset(x = H3K4me3, downsample = 500)
g4 = subset(x = g4, downsample = 500)

transfer.anchors <- FindTransferAnchors(
  reference = H3K4me3,
  query = g4,
  reduction = 'cca',
  query.assay = 'GA',
  reference.assay = 'GA',
  k.filter = NA,
  features = common.genes
)

saveRDS(transfer.anchors,
        glue("../results/Seurat/H3K4me3_G4_transfer_anch.rds"))

genes.use <- VariableFeatures(H3K4me3)
refdata <- GetAssayData(H3K4me3, assay = "GA", slot = "data")

imputation <-
  TransferData(
    anchorset = transfer.anchors,
    refdata = refdata,
    weight.reduction = g4[["lsi"]],
    dims = 1:50
  )
g4[['GA']] <- imputation
coembed <- merge(x = H3K4me3, y = g4)


coembed <- ScaleData(coembed, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:50)

saveRDS(transfer.anchors,
        glue("../results/Seurat/H3K4me3_G4_transfer_anch2.rds"))
