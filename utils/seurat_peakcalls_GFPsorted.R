# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("EnsDb.Mmusculus.v79")
  library("GenomicFeatures")
  library("GenomeInfoDb")
  library("Matrix")
})

set.seed(5)

# Cellranger-ATAC folder
cellranger_folder = "../data/CellRanger/GFP_sorted/"
# path to result folder
result_folder = "../results/Seurat/callpeaks_GFPsorted/"

# read and process data
counts = Read10X_h5(filename = glue("{cellranger_folder}filtered_peak_bc_matrix.h5"))
metadata = read.csv(
  file = glue("{cellranger_folder}singlecell.csv"),
  header = TRUE,
  row.names = 1
)

chrom_assay = CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = glue("{cellranger_folder}fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)

g4 = CreateSeuratObject(counts = chrom_assay,
                        assay = "peaks",
                        meta.data = metadata)

# check peaks in GRanges format
granges(g4)

# give EnsDB mm10 annotation
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) = 'UCSC'
Annotation(g4) = annotations

# normalization
g4 = RunTFIDF(g4)
g4 = FindTopFeatures(g4, min.cutoff = 'q0')
g4 = RunSVD(g4)
DepthCor(g4) # 1st component highly correlate with sequencing depth

# Non-linear dimension reduction and clustering
g4 = RunUMAP(object = g4,
             reduction = 'lsi',
             dims = 2:30)
g4 = FindNeighbors(object = g4,
                   reduction = 'lsi',
                   dims = 2:30)
g4 = FindClusters(object = g4,
                  verbose = FALSE,
                  algorithm = 3)
dim = DimPlot(object = g4, label = TRUE, pt.size = 6) + 
  NoLegend() +
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("GFP-sorted scCutnTag")

# tuning
# for(i in seq(5, 50, by = 5)) {
#   tune = RunUMAP(object = g4,
#                reduction = 'lsi',
#                n.neighbors = i,
#                dims = 2:30)
#   print(DimPlot(object = tune, label = TRUE, pt.size = 6) + 
#     NoLegend() +
#     scale_color_brewer(palette = "Set3") +
#     xlim(-10, 10) + 
#     ylim(-10, 10) + 
#     ggtitle(as.character(i)))
#   
# }
  
ggsave(
  glue("{result_folder}Seurat_GFPsorted_UMAP.png"),
  plot = dim,
  width = 10,
  height = 10,
  dpi = 300,
)

# Create a gene activity matrix
# Signac function: https://stuartlab.org/signac/reference/geneactivity
# Compute counts per cell in gene body and promoter region
gene.activities = GeneActivity(g4)

g4[['GA']] = CreateAssayObject(counts = gene.activities)
g4 = NormalizeData(
  object = g4,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(g4$nCount_peaks)
)

# extract cell barcodes per cluster
barcodes = g4@meta.data %>% 
  rownames_to_column("barcodes") %>%
  dplyr::select(barcodes, seurat_clusters)
write_tsv(barcodes, glue("{result_folder}barcodes_per_cluster.tsv"))

for(cluster in unique(barcodes$seurat_clusters)) {
  subset = barcodes %>% dplyr::filter(seurat_clusters == cluster)
  write_tsv(subset, glue("{result_folder}barcodes_cluster_{as.character(cluster)}.tsv"), col_names = FALSE)
}

# export Rds
saveRDS(g4, glue("{result_folder}GFPsorted.Rds"))

# peak calling per Seurat clusters
peaks = CallPeaks(
  object = g4,
  group.by = "seurat_clusters",
  cleanup = FALSE,
  outdir = result_folder,
  effective.genome.size = 2652783500
)

write.table(
  as.data.frame(peaks),
  file = glue("{result_folder}peaks_per_clusters.bed"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = F
)



