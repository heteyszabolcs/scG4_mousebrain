library("wigglescout")
library("rtracklayer")
library("GenomicRanges")
library("ggplot2")
library("glue")

# path to save images
results = "../results/wigglescout/"

# pairwise comparison of two bigwig files
scatterplot = function(bw1,
                       bw2,
                       bed = "../data/bed/H33_dependent_G4_JL.v2.bed",
                       color = "#a6bddb",
                       bin_size = 10000,
                       export = TRUE) {
  subset = import(bed, format = "BED")
  lab1 = strsplit(strsplit(bw1, "../data/bw/")[[1]][2], ".bw")[[1]][1]
  lab2 = strsplit(strsplit(bw2, "../data/bw/")[[1]][2], ".bw")[[1]][1]
  
  # wigglescout scatterplot
  p = plot_bw_bins_scatter(
    bw1,
    bw2,
    bin_size = bin_size,
    genome = "mm9",
    remove_top = 0.001,
    selection = subset
  )
  
  # modify
  p + geom_point(size = 4,
                 color = color,
                 alpha = 0.5) +
    geom_smooth(
      method = lm,
      linetype = "dashed",
      color = "black",
      fill = "black"
    ) + labs(title = "",
             xlab = lab1,
             ylab = lab2)
  
  if (export) {
    ggsave(
      glue("{results}{lab1}_{lab2}.png"),
      plot = last_plot(),
      width = 10,
      height = 10,
      dpi = 300,
    )
  }
}

# loop through
bigwig_list = list.files("../data/bw/", pattern = "*.bw", full.names = TRUE)

for (i in bigwig_list) {
  lapply(bigwig_list, scatterplot, bw1 = i)
}





