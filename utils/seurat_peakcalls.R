# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("data.table")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("GenomicFeatures")
})

# path to result folder
result_folder = "../results/Seurat/"

# merged G4 scCut&Tag object
g4 = readRDS("../data/merged/Seurat_merged.Rds")
DefaultAssay(g4) = "peaks"

# visualize dimplot
DimPlot(g4, reduction = "umap")

# an example showing how to find peaks open in >X% of a given group of cells:
# Idents(g4) shows the samples (cell types)

get_peaks = function(sample) {
  cells.use = WhichCells(g4, ident = sample)
  counts = GetAssayData(g4, assay = 'peaks', slot = 'counts')[, cells.use]
  percent.open = rowSums(BinarizeCounts(object = counts)) / length(cells.use)
  peaks = names(which(percent.open > 0.05))
  peaks = tibble(sample = sample, "peak" = peaks)
  peaks
}

list_of_peaks = lapply(levels(Idents(g4)), get_peaks)
all_peaks = bind_rows(list_of_peaks)

aggr_peaks = all_peaks %>% group_by(peak) %>% summarise(count = n()) %>% arrange(desc(count))
aggr_samples = all_peaks %>% group_by(sample) %>% summarise(count = n()) %>% arrange(desc(count))

# create aggregated bar plot
bar = ggplot(data = aggr_samples, aes(x = reorder(sample,-count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(title = "Aggregated peak counts / G4 scCut&Tag cell types", 
       x = "Seurat Ident", 
       y = "peaks above 5% of cell type") +
  theme_minimal() +
  theme(
    title = element_text(size = 12, face = 'bold'),
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 0,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(size = 10)
  )
bar

ggsave(
  plot = bar,
  filename = glue("{result_folder}aggregated_G4_peak_counts.png"),
  width = 10,
  height = 8,
  dpi = 500,
  device = "png"
)

# check the common peaks
split <- all_peaks %>%
  split(.$sample)

split$`0` %>%
  semi_join(split$`0`, by = "peak") %>%
  semi_join(split$`1`, by = "peak") %>%
  semi_join(split$`2`, by = "peak") %>%
  semi_join(split$`3`, by = "peak") %>%
  semi_join(split$`4`, by = "peak") %>%
  semi_join(split$`5`, by = "peak") %>%
  semi_join(split$`6`, by = "peak") %>%
  semi_join(split$`7`, by = "peak") %>%
  semi_join(split$`8`, by = "peak") %>%
  pull(peak)
