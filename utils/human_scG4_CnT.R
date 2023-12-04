# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("ggplot2")
  library("EnsDb.Mmusculus.v79")
  library("ensembldb")
  library("GenomicRanges")
  library("dplyr")
  library("glue")
  library("tidyverse")
  library("data.table")
})

# read into Seurat
create_seurat = function(counts,
                         fragments,
                         metadata) {
  chrom_assay = CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "mm10",
    fragments = fragments,
    min.cells = 10,
    min.features = 200
  )
  
  seurat = CreateSeuratObject(counts = chrom_assay,
                              assay = "peaks",
                              meta.data = metadata)
  return(seurat)
}

# FRiP
calculate_frip = function(fragments, seurat) {
  fragments = CountFragments(fragments,
                             cells = NULL,
                             max_lines = NULL,
                             verbose = TRUE)
  cell_ids = fragments$CB
  fragments = data.frame(frequency_count = fragments$frequency_count)
  rownames(fragments) = cell_ids
  fragments = as.matrix(fragments)
  
  
  seurat = AddMetaData(object = seurat,
                       metadata = fragments,
                       col.name = "fragments")
  seurat = FRiP(
    object = seurat,
    assay = 'peaks',
    total.fragments = "fragments",
    col.name = "FRiP"
  )
  
  return(seurat)
  
}

# cellranger output source:
# human G4 scCut&Tag (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181373)
mcf7_counts = Read10X_h5(filename = "../data/GSE181373/MCF7_G4_scCT/cellranger_output/filtered_peak_bc_matrix.h5")
mcf7_meta = read.csv(file = "../data/GSE181373/MCF7_G4_scCT/cellranger_output/singlecell.csv",
                     header = TRUE,
                     row.names = 1)
mcf7_fragments = "../data/GSE181373/MCF7_G4_scCT/cellranger_output/fragments.tsv.gz"

u2os_counts = Read10X_h5(filename = "../data/GSE181373/U2OS_G4_scCT/cellranger_output/filtered_peak_bc_matrix.h5")
u2os_meta = read.csv(file = "../data/GSE181373/U2OS_G4_scCT/cellranger_output/singlecell.csv",
                     header = TRUE,
                     row.names = 1)
u2os_fragments = "../data/GSE181373/U2OS_G4_scCT/cellranger_output/fragments.tsv.gz"


# compute frip
mcf7_seurat = create_seurat(counts = mcf7_counts,
                            metadata = mcf7_meta,
                            fragments = mcf7_fragments)

mcf7_frip = calculate_frip(fragments = mcf7_fragments, 
                             seurat = mcf7_seurat)

u2os_seurat = create_seurat(counts = u2os_counts,
                            metadata = u2os_meta,
                            fragments = u2os_fragments)

u2os_frip = calculate_frip(fragments = u2os_fragments, 
                           seurat = u2os_seurat)

mcf7_tibble = tibble(FRiP = mcf7_frip@meta.data$FRiP, data = "MCF7")
u2os_tibble = tibble(FRiP = u2os_frip@meta.data$FRiP, data = "U2OS")
input = rbind(mcf7_tibble, u2os_tibble)

ggplot(input, aes(x = data, y = FRiP, fill = data)) +
  geom_violin(trim = FALSE, color = "black") + 
  scale_fill_manual(values = c("#bcbddc", "#bcbddc")) +
  ylim(0, 1) +
  labs(
    title = "Fraction of reads in peaks (FRiP)",
    x = "experiment",
    y = " ",
    fill = " "
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.title.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
# stat_compare_means(label.y = 1, label.x = 1.5) +
# stat_compare_means(label = "p.signif", method = "t.test",
#                    ref.group = ".all.")




