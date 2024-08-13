# packages
suppressPackageStartupMessages({
  library("EnsDb.Mmusculus.v75")
  library("Seurat")
  library("Signac")
  library("ggpubr")
})

# function for creating quality violin plots
qc = function(filtered_peak_mat = "../data/CellRanger/unsorted/filtered_peak_bc_matrix.h5",
              metadata = "../data/CellRanger/unsorted/singlecell.csv",
              fragments = "../data/CellRanger/unsorted/fragments.tsv.gz",
              ensdb = EnsDb.Mmusculus.v75,
              label = "unsorted",
              res = 0.1) {
  counts <- Read10X_h5(filename = filtered_peak_mat)
  metadata <- read.csv(
    file = metadata,
    header = TRUE,
    row.names = 1
  )
  chrom_assay = CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "mm10",
    fragments = fragments,
    min.cells = 10,
    min.features = 200
  )
  
  seurat = CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  # nucleosome signal score
  seurat = NucleosomeSignal(object = seurat)
  
  # compute TSS enrichment score per cell
  annotations = GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevels(annotations) = paste0('chr', seqlevels(annotations))
  genome(annotations) = "mm10"
  
  Annotation(seurat) = annotations
  seurat = TSSEnrichment(object = seurat, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks (FRIP)
  seurat$pct_reads_in_peaks = seurat$peak_region_fragments / seurat$passed_filters * 100
  seurat$blacklist_ratio = seurat$blacklist_region_fragments / seurat$peak_region_fragments
  
  DensityScatter(seurat, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
  ggsave(
    glue("../results/Seurat/{label}_density_scatter.pdf"),
    plot = last_plot(),
    width = 6,
    height = 4,
    dpi = 300,
  )
  
  ## processing
  # normalization
  seurat = RunTFIDF(seurat)
  seurat = FindTopFeatures(seurat, min.cutoff = 'q0')
  seurat = RunSVD(seurat)
  
  # Non-linear dimension reduction and clustering
  seurat = RunUMAP(object = seurat,
                   reduction = 'lsi',
                   dims = 2:30)
  # clustering (community detection)
  seurat = FindNeighbors(object = seurat,
                         reduction = 'lsi',
                         dims = 2:30)
  seurat = FindClusters(object = seurat,
                        verbose = FALSE,
                        resolution = res,
                        algorithm = 3)
  
  
  seurat$high.tss = ifelse(seurat$TSS.enrichment > 3, 'High', 'Low')
  TSSPlot(seurat, group.by = 'high.tss') + NoLegend()
  
  tss_enrich = VlnPlot(
    object = seurat,
    features = 'TSS.enrichment',
    group.by = "seurat_clusters",
    pt.size = 0.1,
    ncol = 1,
    cols = c("#addd8e", "#f0f0f0", "#636363")
  ) +
    xlab("cluster") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoLegend()
  nucl_signal = VlnPlot(
    object = seurat,
    features = 'nucleosome_signal',
    group.by = "seurat_clusters",
    pt.size = 0.1,
    ncol = 1,
    cols = c("#addd8e", "#f0f0f0", "#636363")
  ) +
    xlab("cluster") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoLegend()
  pct_in_reads = VlnPlot(
    object = seurat,
    features = 'pct_reads_in_peaks',
    group.by = "seurat_clusters",
    pt.size = 0.1,
    ncol = 1,
    cols = c("#addd8e", "#f0f0f0", "#636363")
  ) +
    xlab("cluster") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoLegend()
  nF_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "nFeature_peaks",
    pt.size = 0.1,
    ncol = 1,
    cols = c("#addd8e", "#f0f0f0", "#636363")
  ) +
    xlab("cluster") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(
        size = 25,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoLegend()
  nC_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "nCount_peaks",
    pt.size = 0.1,
    ncol = 1,
    cols = c("#addd8e", "#f0f0f0", "#636363")
  ) +
    xlab("cluster") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(
        size = 25,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 25, color = "black")
    )
  TSS_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "TSS_fragments",
    pt.size = 0.1,
    ncol = 1,
    cols = c("#addd8e", "#f0f0f0", "#636363")
  ) +
    xlab("cluster") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(
        size = 25,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 25, color = "black")
    )
  mito_violin = VlnPlot(
    seurat,
    group.by = "seurat_clusters",
    features = "mitochondrial",
    pt.size = 0.1,
    ncol = 1,
    cols = c("#addd8e", "#f0f0f0", "#636363")
  ) +
    xlab("cluster") +
    theme(
      text = element_text(size = 25),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(
        size = 25,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(size = 25, color = "black")
    )

  qc_violins = ggarrange(tss_enrich, nucl_signal, pct_in_reads,
                         nF_violin, nC_violin, TSS_violin, mito_violin)  
  
  ggsave(
    glue("../results/Seurat/{label}_quality_plots.pdf"),
    plot = qc_violins,
    width = 12,
    height = 10,
    dpi = 300,
  )
  
  return(qc_violins)
  
}

# run
qc()













