# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("EnsDb.Mmusculus.v79")
  library("GenomicFeatures")
  library("GenomeInfoDb")
  library("cowplot")
  library("ggpubr")
})


set.seed(5)

# Cellranger-ATAC folder
cellranger_folder = "../data/CellRanger/unsorted/"
# path to result folder
result_folder = "../results/Seurat/"

# Seurat objects
unsorted_g4 = readRDS("../results/Seurat/callpeaks_unsorted/unsorted.Rds")
mesc_mef_g4 = readRDS("../results/Seurat/callpeaks_mESC-MEF/mESC-MEF.Rds")

calculate_frip = function(fragments, seurat) {
  
  fragments = CountFragments(fragments, cells = NULL, max_lines = NULL, verbose = TRUE)
  cell_ids = fragments$CB 
  fragments = data.frame(frequency_count = fragments$frequency_count)
  rownames(fragments) = cell_ids
  fragments = as.matrix(fragments)
  
  
  seurat = AddMetaData(
    object = seurat,
    metadata = fragments,
    col.name = "fragments"
  )
  seurat = FRiP(object = seurat, assay = 'peaks', total.fragments = "fragments", col.name = "FRiP")
  
  return(seurat)
  
}



unsorted_g4 = calculate_frip(fragments = "../data/CellRanger/unsorted/fragments.tsv.gz", seurat = unsorted_g4)
mesc_mef_g4 = calculate_frip(fragments = "../data/CellRanger/mES-mEF/fragments.tsv.gz", seurat = mesc_mef_g4)

unsorted_tibble = tibble(FRiP = unsorted_g4@meta.data$FRiP, data = "unsorted brain")
mesc_mef_tibble = tibble(FRiP = mesc_mef_g4@meta.data$FRiP, data = "mESC-MEF")
input = rbind(unsorted_tibble, mesc_mef_tibble)

ggplot(input, aes(x = data, y = FRiP, fill = data)) +
  geom_violin(trim = FALSE, color = "black") + 
  scale_fill_manual(values = c("#9ecae1", "#a1d99b")) +
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

ggsave(
  glue("{result_folder}FRiP_violins.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 10,
  height = 10,
  dpi = 300,
)

# unsorted_tibble = tibble(reads_count = CountFragments("../data/CellRanger/unsorted/fragments.tsv.gz", cells = NULL, max_lines = NULL, verbose = TRUE)$mononucleosomal, 
#                          data = "unsorted brain")
# mesc_mef_tibble = tibble(reads_count = CountFragments("../data/CellRanger/mES-mEF/fragments.tsv.gz", cells = NULL, max_lines = NULL, verbose = TRUE)$mononucleosomal, 
#                          data = "mESC-MEF")
# input = rbind(unsorted_tibble, mesc_mef_tibble)
# 
# ggplot(input, aes(x = data, y = log10(reads_count), fill = data)) +
#   geom_violin(trim = FALSE, color = "black") + 
#   scale_fill_manual(values = c("#9ecae1", "#a1d99b")) +
#   ylim(0, 2) +
#   labs(
#     title = "Mononucleosomal",
#     x = "experiment",
#     y = " ",
#     fill = " "
#   ) +
#   theme_classic() +
#   theme(
#     text = element_text(size = 15),
#     plot.title = element_text(size = 25),
#     legend.text = element_text(size = 25),
#     axis.text.x = element_text(size = 25, color = "black"),
#     axis.title.x = element_text(size = 25, color = "black"),
#     axis.text.y = element_text(size = 25, color = "black")
#   )
# stat_compare_means(label.y = 1, label.x = 1.5) +
# stat_compare_means(label = "p.signif", method = "t.test",
#                    ref.group = ".all.")

