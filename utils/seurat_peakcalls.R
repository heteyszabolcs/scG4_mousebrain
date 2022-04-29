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
})

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
dim = DimPlot(object = g4, label = TRUE) + NoLegend()
ggsave(
  glue("{result_folder}Seurat_GFP_sorted_UMAP.png"),
  plot = dim,
  width = 10,
  height = 10,
  dpi = 300,
)

# Create a gene activity matrix
gene.activities = GeneActivity(g4)

g4[['GA']] = CreateAssayObject(counts = gene.activities)
g4 = NormalizeData(
  object = g4,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(g4$nCount_peaks)
)

# export Rds
saveRDS(g4, glue("{result_folder}GFP_sorted.Rds"))

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



