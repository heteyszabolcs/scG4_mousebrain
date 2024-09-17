suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("wigglescout")
  library("ggrastr")
  library("ggpubr")
})

result_folder = "../results/Seurat/callpeaks_unsorted/"

# scatterplot to visualize signal relationship and PQS sites
pqs_scatter = function(bigwig_path_1,
                       bigwig_path_2,
                       pqs_path,
                       peak_path,
                       label_1,
                       label_2,
                       peak_label,
                       output_name,
                       color = "#addd8e") {
  # peaks
  peaks = rtracklayer::import(peak_path)
  
  # PQS sites
  pqs = rtracklayer::import(pqs_path)
  start(pqs) = start(pqs) - 50
  end(pqs) = end(pqs) + 50
  
  # wigglescout
  aggr_signal_1 <- bw_loci(
    bigwig_path_1,
    bg_bwfiles = NULL,
    loci = pqs,
    per_locus_stat = "max",
    labels = "score"
  )
  aggr_signal_2 <- bw_loci(
    bigwig_path_2,
    bg_bwfiles = NULL,
    loci = pqs,
    per_locus_stat = "max",
    labels = "score"
  )
  
  # scatter plot input
  ranges = as.data.frame(aggr_signal_1)
  ranges = ranges %>% mutate(seqnames = paste(seqnames, start, end, sep = "_"))

  input = tibble(cluster_a = aggr_signal_1$score, cluster_b = aggr_signal_2$score)
  input = input %>% mutate(type = "PQS", seqnames = ranges$seqnames)
  
  # the query (PQS) interval must be wholly contained within the G4 peak
  input[unique(queryHits(findOverlaps(
    aggr_signal_2, peaks, minoverlap = 0, type = "within"
  ))), "type"] =
    glue("PQS + {peak_label} peak")
  
  # collect where regions are zero
  excl = input %>% dplyr::filter(cluster_a == 0 &
                                   cluster_b == 0 &
                                   type == glue("PQS + {peak_label} peak"))
  
  input = input %>% dplyr::filter(!seqnames %in% excl$seqnames) %>% 
    mutate(cluster_a = log2(cluster_a + 1),
           cluster_b = log2(cluster_b + 1))
  
  input$type = factor(input$type, levels = c("PQS", glue("PQS + {peak_label} peak")))
  input = input[order(input$type == "PQS", decreasing = TRUE), ]
    
  pqs_and_peak = 
    input %>% dplyr::filter(type == glue("PQS + {peak_label} peak")) %>%
    pull(seqnames) %>% 
    unique %>% 
    length
  
  print(glue("Number of PQS & cluster peak sites: {pqs_and_peak}
              Number of all peaks: {as.character(length(peaks))}"))
  
  # scatter plot
  p = ggplot(input, aes(x = cluster_a, y = cluster_b,
                                  color = type)) + 
    ggrastr::geom_point_rast(size = 0.5) +
    scale_color_manual(values = c("#f0f0f0", color)) +
    labs(
      title = "",
      x = label_1,
      y = label_2,
      color = ""
    ) +
    xlim(0, 15) +
    ylim(0, 15) +
    theme_classic() +
    theme(
      text = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 20, color = "black"),
      axis.title.y = element_text(size = 20, color = "black")
    )
  p = rasterize(p, layers = 'Point', dpi = 150)
  
  ggsave(
    glue("{result_folder}{output_name}.pdf"),
    width = 7,
    height = 5,
    plot = p,
    device = "pdf"
  )
  
  return(print(p))
  
}

## unsorted G4 scCUT&Tag
# res 0.1
pqs_scatter(
  bigwig_path_1 = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/0.bam_RPGC.bigwig",
  bigwig_path_2 = "../data/bw/ATAC-H33WT.mm10.bw",
  label_1 = "cluster 0",
  label_2 = "mESC ATAC-Seq",
  output_name = "unsorted_res0.1_peak_scatterplot-cl0_ATAC",
  pqs_path = "../data/bed/mm10_canonical_G4_PQS-regex.bed",
  peak_path = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak",
  color = "#addd8e",
  peak_label = "cluster 0"
)

pqs_scatter(
  bigwig_path_1 = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_bigwigs/1.bam_RPGC.bigwig",
  bigwig_path_2 = "../data/bw/ATAC-H33WT.mm10.bw",
  label_1 = "cluster 1",
  label_2 = "mESC ATAC-Seq",
  output_name = "unsorted_res0.1_peak_scatterplot-cl1_ATAC",
  pqs_path = "../data/bed/mm10_canonical_G4_PQS-regex.bed",
  peak_path = "../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak",
  color = "#bcbcbc",
  peak_label = "cluster 1"
)


# res 0.8
bigwigs = c("../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/0.bam_RPGC.bigwig",
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/1.bam_RPGC.bigwig",
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/2.bam_RPGC.bigwig",
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/3.bam_RPGC.bigwig",
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/4.bam_RPGC.bigwig")
peaks = c("../results/GenomicRanges/unsorted_outputs/unsorted-unique_cl0_peaks.bed",
            "../results/GenomicRanges/unsorted_outputs/unsorted-unique_cl1_peaks.bed",
            "../results/GenomicRanges/unsorted_outputs/unsorted-unique_cl2_peaks.bed",
            "../results/GenomicRanges/unsorted_outputs/unsorted-unique_cl3_peaks.bed",
            "../results/GenomicRanges/unsorted_outputs/unsorted-unique_cl4_peaks.bed")
peak_labels = c("unique cl. 0", "unique cl. 1", "unique cl. 2", "unique cl. 3", "unique cl. 4")

l = list()
for (i in 1:length(bigwigs)) {
  print(paste0("Working on: ", bigwigs[i]))
  for (j in 1:length(bigwigs)) {
    p = pqs_scatter(
      bigwig_path_1 = bigwigs[i],
      bigwig_path_2 = bigwigs[j],
      label_1 = paste0("cluster ", strsplit(
        strsplit(
          bigwigs[i],
          "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"
        )[[1]][2],
        ".bam_RPGC.bigwig"
      )[[1]][1]),
      label_2 = paste0("cluster ", strsplit(
        strsplit(
          bigwigs[j],
          "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"
        )[[1]][2],
        ".bam_RPGC.bigwig"
      )[[1]][1]),
      output_name = paste0(
        "unsorted_",
        strsplit(
          strsplit(
            bigwigs[i],
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"
          )[[1]][2],
          ".bam_RPGC.bigwig"
        )[[1]][1],
        "_vs_",
        strsplit(
          strsplit(
            bigwigs[j],
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"
          )[[1]][2],
          ".bam_RPGC.bigwig"
        )[[1]][1]
      ),
      pqs_path = "../data/bed/mm10_canonical_G4_PQS-regex.bed",
      peak_path = peaks[i],
      peak_label = peak_labels[i],
      color = "#addd8e"
    )
    l[[paste0(
        "unsorted_",
        strsplit(
          strsplit(
            bigwigs[i],
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"
          )[[1]][2],
          ".bam_RPGC.bigwig"
        )[[1]][1],
        "_vs_",
        strsplit(
          strsplit(
            bigwigs[j],
            "../results/Seurat/final/unsorted_brain/res0.8/cluster_spec_bigwigs/"
          )[[1]][2],
          ".bam_RPGC.bigwig"
        )[[1]][1]
      )]] = p
  }
}

# scatterplot matrix
arranged = ggarrange(plotlist = l[16:25])

ggsave(
  glue("{result_folder}unsorted_signal_at_PQS-sc_matrix_3.pdf"),
  width = 28,
  height = 12,
  plot = arranged,
  device = "pdf"
)


