# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("ggplot2")
  library("data.table")
  library("GenomicFeatures")
  library("GenomicRanges")
  library("cowplot")
})

# result folder
result_folder = "../results/Seurat/callpeaks_unsorted/"
peak_folder = "../results/Seurat/callpeaks_unsorted/peak_sets/"
peak_list = list.files(glue("{peak_folder}"), pattern = "[0-9]_peaks.narrowPeak$")

# customized upset plot function
plot_freq_intersect <-
  function(dat,
           .by = "Group",
           .levels = NA,
           .split = "Category",
           .color = 1:10,
           top_n = 10) {
    # Args:
    #   dat: gene_id (unique)
    #   .by: row feature
    #   .levels: dot plot row names
    #   .split: bar split feature
    
    if (is.na(.levels[1])) {
      .levels <- seq_len(as.character(nchar(dat[1, .by])))
    }
    
    # limit less frequent groups
    top_groups <- names(tail(sort(table(dat[, .by])), top_n))
    dat <- dat[dat[[.by]] %in% top_groups,]
    
    dat$Group <-
      factor(dat[, .by], levels = names(sort(table(dat[, .by]), decreasing = T)))
    dat$Type <- dat[, .split]
    
    # barplot
    dat_g1 <- dplyr::count(dat, Group, Type, sort = TRUE)
    
    frac_tbl <- table(dat[, "Group"], dat[, "Type"])
    frac_tbl <- frac_tbl / rowSums(frac_tbl)
    
    dat_g1_text <- data.frame(Group = unique(dat_g1$Group),
                              #Fraction = paste0(round(frac_tbl[, 1], 2) * 100, "%"),
                              Type = colnames(frac_tbl)[1])
    
    Group_n <- dplyr::count(dat, Group, sort = TRUE)
    dat_g1_text$n  <-
      Group_n[match(dat_g1_text$Group, Group_n$Group), "n"]
    
    g1 <- ggplot(dat_g1, aes(x = Group, y = n, fill = Type)) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_text(
        data = dat_g1_text,
        aes(x = Group, y = n, label = n),
        vjust = -0.25,
        hjust = 0.5,
        size = 5
      ) +
      xlab("") + ylab("Number of regions") +
      theme_classic() +
      scale_fill_manual(name = .split, values = c("#636363", .color)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(c(0, 1, 1, 1))
      )
    
    # dot plot
    dat_tile <- NULL
    for (m in seq_along(levels(dat$Group))) {
      tmp <- cbind(x = m,
                   y = .levels,
                   color = ifelse(strsplit(levels(dat$Group)[m], "")[[1]] == 1, 1, 0))
      dat_tile <- rbind(dat_tile, tmp)
    }
    
    dat_tile <- as.data.frame(dat_tile)
    dat_tile$x <-
      factor(dat_tile$x, levels = order(as.numeric(dat_tile$x)))
    dat_tile$y <- factor(dat_tile$y, levels = rev(.levels))
    
    g2 <-
      ggplot(dat_tile, aes(
        x = x,
        y = y,
        color = factor(color),
        group = x
      )) +
      geom_point(size = 4) +
      scale_color_manual(values = c("grey90", .color)) +
      xlab("") + ylab("") +
      ggpubr::theme_pubclean() +
      theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(c(1, 1, 0, 1))
      )
    
    # add lines
    dat_seg <- NULL
    for (i in unique(dat_tile$x)) {
      tmp <- dat_tile[dat_tile$x == i,]
      tmp <- tmp[tmp$color == 1,]
      if (nrow(tmp) > 1) {
        dat_seg <- rbind(dat_seg,
                         data.frame(
                           x = tmp$x[-nrow(tmp)],
                           y = tmp$y[-nrow(tmp)],
                           xend = tmp$x[-1],
                           yend = tmp$y[-1]
                         ))
      }
    }
    
    if (!is.null(dat_seg)) {
      levels(dat_seg$x) <- levels(dat_seg$xend) <- levels(dat_tile$x)
      levels(dat_seg$y) <- levels(dat_seg$yend) <- levels(dat_tile$y)
      dat_seg$color <- 1
      
      g2 <-
        g2 + geom_segment(data = dat_seg,
                          aes(
                            x = x,
                            y = y,
                            xend = xend,
                            yend = yend
                          ),
                          lty = 2)
    }
    
    cowplot::plot_grid(
      g1,
      g2,
      ncol = 1,
      align = "v",
      rel_heights = c(1, 0.4)
    )
  }

overlap = function(peak_set1, peak_set2) {
  ol = peak_set1[queryHits(findOverlaps(
    peak_set1,
    peak_set2,
    maxgap = -1L,
    minoverlap = 1,
    type = c("any")
  )), ]
  
  ol = as_tibble(ol)
  ol = ol %>% mutate(id = paste0(seqnames, "_", start, "_", end)) %>% pull(id)
  return(ol)
}

# create ID column
add_id = function(peak_set = peak0) {
  peak_set = as_tibble(peak_set)
  ids = peak_set %>% mutate(id = paste0(seqnames, "_", start, "_", end)) %>% pull(id)
  return(ids)
}

# create input data frame for upset plot
upset_input = function(peak_set, bulk_set) {
  peak_set_t = as_tibble(peak_set)
  peak_set_t = peak_set_t %>% mutate(id = paste0(seqnames, "_", start, "_", end))
  
  # create ids for upset plot
  groups = paste0(as.numeric(
    peak_set_t$id %in% overlap(peak_set1 = peak_set, peak_set2 = cluster0)
  ),
  paste0(
    as.numeric(
      peak_set_t$id %in% overlap(peak_set1 = peak_set, peak_set2 = cluster1)
    ),
    paste0(
      as.numeric(
        peak_set_t$id %in% overlap(peak_set1 = peak_set, peak_set2 = bulk_set)
      )
      )
    ))
  
  input = tibble(id = peak_set_t$id,
                 category = "",
                 group = groups)
  input = as.data.frame(input)
  
  return(input)
}

# FindAllMarkers region
# markers = fread("../results/Seurat/final/unsorted_brain/res0.1/outputs/FindAllMarkers_logreg_output.tsv")
# cluster1 = markers %>% dplyr::filter(cluster == 1) %>% dplyr::select(chr, start, end) %>% 
#   dplyr::select(V1 = chr, V2 = start, V3 = end)

# Signac MACS2 peaks
cluster0 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")
cluster0$type = "cluster0"
cluster0 = GRanges(
  seqnames = cluster0$V1,
  ranges = IRanges(
    start = cluster0$V2,
    end = cluster0$V3,
    names = cluster0$type,
  )
)


cluster1 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
cluster1$type = "cluster1"
cluster1 = GRanges(
  seqnames = cluster1$V1,
  ranges = IRanges(
    start = cluster1$V2,
    end = cluster1$V3,
    names = cluster1$type,
  )
)

bulk_npc_peak = fread("../data/bed/bulk_CnT_G4_NPC_mm10.bed")
bulk_npc_peak$type = "NPC_bulk_CnT"
bulk_npc_peak = GRanges(
  seqnames = bulk_npc_peak$V1,
  ranges = IRanges(
    start = bulk_npc_peak$V2,
    end = bulk_npc_peak$V3,
    names = bulk_npc_peak$type,
  )
)

bulk_neuron_peak = fread("../data/bed/bulk_CnT_G4_neuron_mm10.bed")
bulk_neuron_peak$type = "neuron_bulk_CnT"
bulk_neuron_peak = GRanges(
  seqnames = bulk_neuron_peak$V1,
  ranges = IRanges(
    start = bulk_neuron_peak$V2,
    end = bulk_neuron_peak$V3,
    names = bulk_neuron_peak$type,
  )
)

cluster0$type = "cluster 0"
cluster1$type = "cluster 1"
bulk_npc_peak$type = "bulk npc"
bulk_neuron_peak$type = "bulk neuron"

# concatenate peak lists
peak_all =
  GenomicRanges::reduce(c(bulk_npc_peak, cluster1, cluster0))

input = upset_input(peak_all, bulk_set = bulk_npc_peak)
upset_plot = plot_freq_intersect(
  input,
  .by = "group",
  .levels = c("cluster 0", "cluster 1", "bulk NPC"),
  .split = "category",
  .color = "#a1d99b",
  top_n = 10
)
upset_plot

ggsave(
  glue("{result_folder}bulk_NPC-unsorted_brain_cluster0_1_upsetplot.png"),
  plot = upset_plot,
  width = 6,
  height = 6,
  dpi = 500,
)
ggsave(
  glue("{result_folder}bulk_NPC-unsorted_brain_cluster0_1_upsetplot.pdf"),
  plot = upset_plot,
  width = 6,
  height = 6,
  device = "pdf"
)

peak_all =
  GenomicRanges::reduce(c(bulk_neuron_peak, cluster1, cluster0))
input = upset_input(peak_all, bulk_set = bulk_neuron_peak)
upset_plot = plot_freq_intersect(
  input,
  .by = "group",
  .levels = c("cluster 0", "cluster 1", "bulk neuron"),
  .split = "category",
  .color = "#a1d99b",
  top_n = 10
)
upset_plot

# differentiate cluster0, and cluster1
plot_euler(
  list(cluster1, cluster0),
  ignore.strand = TRUE,
  fills = c("white", "#2ca25f"),
  names = c("cluster 1", "cluster 0")
)

plot_euler(
  list(bulk_neuron_peak, cluster1, cluster0),
  ignore.strand = TRUE,
  fills = c("white", "#bdbdbd", "#2ca25f"),
  names = c("bulk neuron", "cluster 1", "cluster 0")
)

ggsave(
  glue("{result_folder}bulk_neuron-unsorted_brain_cluster0_1_upsetplot.png"),
  plot = upset_plot,
  width = 6,
  height = 6,
  dpi = 500,
)
ggsave(
  glue("{result_folder}bulk_neuron-unsorted_brain_cluster0_1_upsetplot.pdf"),
  plot = upset_plot,
  width = 6,
  height = 6,
  device = "pdf"
)

