# packages
suppressPackageStartupMessages({
  library("glue")
  library("ggplot2")
  library("ggrepel")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("cowplot")
})

enhancer_anal = "../results/GenomicRanges/"
collect = fread(glue("{enhancer_anal}enhancer_analysis_output.tsv"))

cm_bar = collect %>%
  mutate(Cruz_Molina_enh = as.numeric(Cruz_Molina_enh)) %>%
  group_by(Seurat_cluster) %>%
  summarize(sum = sum(Cruz_Molina_enh)) %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  ungroup() %>%
  ggplot(data = ., aes(
    x = reorder(Seurat_cluster, -sum),
    y = sum,
    fill = Seurat_cluster
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  ylim(0, 12) +
  scale_y_continuous(breaks = seq(0, 12, 1)) +
  labs(
    title = expression(
      paste(
        "G4 (peaks above 75th perc) overlaps with active enhancers of ",
        italic("Cruz-Molina et al.")
      )
    ),
    x = "Seurat cluster",
    y = "# of G4 - active enhancer overlaps",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
cm_bar

ggsave(
  glue("{enhancer_anal}signG4_overlaps_w_CruzM_active_enh.png"),
  plot = cm_bar,
  width = 10,
  height = 10,
  dpi = 300,
)

gl_bar = collect %>%
  mutate(Glaser_enh = as.numeric(Glaser_enh)) %>%
  group_by(Seurat_cluster) %>%
  summarize(sum = sum(Glaser_enh)) %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  ungroup() %>%
  ggplot(data = ., aes(
    x = reorder(Seurat_cluster, -sum),
    y = sum,
    fill = Seurat_cluster
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  ylim(0, 12) +
  scale_y_continuous(breaks = seq(0, 12, 1)) +
  labs(
    title = expression(
      paste(
        "G4 (peaks above 75th perc) overlaps with active enhancers of ",
        italic("Glaser et al.")
      )
    ),
    x = "Seurat cluster",
    y = "# of G4 - active enhancer overlaps",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
gl_bar

ggsave(
  glue("{enhancer_anal}signG4_overlaps_w_Glaser_active_enh.png"),
  plot = gl_bar,
  width = 10,
  height = 10,
  dpi = 300,
)

tss_rank = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  mutate(Distance_to_TSS = Distance_to_TSS / 1000) %>%
  ggplot(data = ., aes(x = reorder(id,-Distance_to_TSS), y = Distance_to_TSS)) +
  geom_point(
    stat = 'identity',
    aes(col = Seurat_cluster),
    size = 3,
    alpha = 0.4
  ) +
  scale_color_brewer(palette = "YlOrRd") +
  labs(
    title = expression(paste("TSS distances of significant G4 peaks")),
    x = "G4 peaks above 75th percentile",
    y = "Distance to TSS (kb)",
    color = "Seurat cluster"
  ) +
  geom_label_repel(
    aes(label = ifelse(
      Cruz_Molina_enh == 1 |
        Glaser_enh == 1,
      as.character(Gene_name),
      ''
    )),
    box.padding   = 0.9,
    max.overlaps = Inf,
    point.padding = 0.5,
    size = 2,
    segment.color = 'black'
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 20),
    plot.title = element_text(size = 15)
  )
tss_rank

ggsave(
  glue("{enhancer_anal}G4_TSS_distance_rank-labeled.png"),
  plot = tss_rank,
  width = 10,
  height = 5,
  dpi = 300,
)

tss_rank = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  mutate(Distance_to_TSS = Distance_to_TSS / 1000) %>%
  ggplot(data = ., aes(x = reorder(id,-Distance_to_TSS), y = Distance_to_TSS)) +
  geom_point(
    stat = 'identity',
    aes(col = Seurat_cluster),
    size = 3,
    alpha = 0.4
  ) +
  scale_color_brewer(palette = "YlOrRd") +
  labs(
    title = expression(paste("TSS distances of significant G4 peaks")),
    x = "G4 peaks above 75th percentile",
    y = "Distance to TSS (kb)",
    color = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 20),
    plot.title = element_text(size = 15)
  )
tss_rank

ggsave(
  glue("{enhancer_anal}G4_TSS_distance_rank.png"),
  plot = tss_rank,
  width = 10,
  height = 5,
  dpi = 300,
)

peak_rank = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(col = Seurat_cluster),
    size = 3,
    alpha = 0.4
  ) +
  scale_color_brewer(palette = "YlOrRd") +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "Seurat cluster"
  ) +
  geom_label_repel(
    aes(label = ifelse(
      Cruz_Molina_enh == 1 |
        Glaser_enh == 1,
      as.character(Gene_name),
      ''
    )),
    box.padding   = 0.9,
    max.overlaps = Inf,
    point.padding = 0.5,
    size = 2,
    segment.color = 'black'
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 20),
    plot.title = element_text(size = 15)
  )
peak_rank

ggsave(
  glue("{enhancer_anal}G4_peak_score_rank-labeled.png"),
  plot = peak_rank,
  width = 10,
  height = 5,
  dpi = 300,
)

peak_rank = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(col = Seurat_cluster),
    size = 3,
    alpha = 0.4
  ) +
  scale_color_brewer(palette = "YlOrRd") +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "Seurat cluster"
  ) +
  # geom_label_repel(aes(label = ifelse(Cruz_Molina_enh == 1 | Glaser_enh == 1, as.character(Gene_name), '')),
  #                  box.padding   = 0.9,
  #                  max.overlaps = Inf,
  #                  point.padding = 0.5,
  #                  size = 2,
  #                  segment.color = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 20),
    plot.title = element_text(size = 15)
  )
peak_rank

ggsave(
  glue("{enhancer_anal}G4_peak_score_rank.png"),
  plot = peak_rank,
  width = 10,
  height = 5,
  dpi = 300,
)


collect[is.na(collect)] = 0

peak_rank_ac_astro = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(color = H3K27ac_Astrocytes.bw),
    size = 3,
    alpha = 1
  ) +
  scale_color_gradient2(
    midpoint = mean(collect$H3K27ac_Astrocytes.bw),
    low = "blue",
    mid = "white",
    high = "red",
    space = "Lab"
  ) +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "H3K27ac - Astrocytes"
  ) +
  # geom_label_repel(aes(label = ifelse(Cruz_Molina_enh == 1 | Glaser_enh == 1, as.character(Gene_name), '')),
  #                  box.padding   = 0.9,
  #                  max.overlaps = Inf,
  #                  point.padding = 0.5,
  #                  size = 2,
  #                  segment.color = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
peak_rank_ac_astro

peak_rank_ac_astro = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(color = H3K27ac_Astrocytes.bw),
    size = 3,
    alpha = 1
  ) +
  scale_color_gradient2(
    midpoint = mean(collect$H3K27ac_Astrocytes.bw),
    low = "blue",
    mid = "white",
    high = "red",
    space = "Lab",
    limits = c(0, 15)
  ) +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "H3K27ac - Astrocytes"
  ) +
  # geom_label_repel(aes(label = ifelse(Cruz_Molina_enh == 1 | Glaser_enh == 1, as.character(Gene_name), '')),
  #                  box.padding   = 0.9,
  #                  max.overlaps = Inf,
  #                  point.padding = 0.5,
  #                  size = 2,
  #                  segment.color = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
peak_rank_ac_astro

peak_rank_ac_mol = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(color = H3K27ac_mOL.bw),
    size = 3,
    alpha = 1
  ) +
  scale_color_gradient2(
    midpoint = mean(collect$H3K27ac_mOL.bw),
    low = "blue",
    mid = "white",
    high = "red",
    space = "Lab",
    limits = c(0, 15)
  ) +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "H3K27ac - mOL"
  ) +
  # geom_label_repel(aes(label = ifelse(Cruz_Molina_enh == 1 | Glaser_enh == 1, as.character(Gene_name), '')),
  #                  box.padding   = 0.9,
  #                  max.overlaps = Inf,
  #                  point.padding = 0.5,
  #                  size = 2,
  #                  segment.color = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
peak_rank_ac_mol

peak_rank_ac_oec = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(color = H3K27ac_OEC.bw),
    size = 3,
    alpha = 1
  ) +
  scale_color_gradient2(
    midpoint = mean(collect$H3K27ac_OEC.bw),
    low = "blue",
    mid = "white",
    high = "red",
    space = "Lab",
    limits = c(0, 15)
  ) +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "H3K27ac - OEC"
  ) +
  # geom_label_repel(aes(label = ifelse(Cruz_Molina_enh == 1 | Glaser_enh == 1, as.character(Gene_name), '')),
  #                  box.padding   = 0.9,
  #                  max.overlaps = Inf,
  #                  point.padding = 0.5,
  #                  size = 2,
  #                  segment.color = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
peak_rank_ac_oec

peak_rank_ac_opc = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(color = H3K27ac_OPC.bw),
    size = 3,
    alpha = 1
  ) +
  scale_color_gradient2(
    midpoint = mean(collect$H3K27ac_OPC.bw),
    low = "blue",
    mid = "white",
    high = "red",
    space = "Lab",
    limits = c(0, 15)
  ) +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "H3K27ac - OPC"
  ) +
  # geom_label_repel(aes(label = ifelse(Cruz_Molina_enh == 1 | Glaser_enh == 1, as.character(Gene_name), '')),
  #                  box.padding   = 0.9,
  #                  max.overlaps = Inf,
  #                  point.padding = 0.5,
  #                  size = 2,
  #                  segment.color = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
peak_rank_ac_opc

peak_rank_ac_vlmc = collect %>%
  mutate(Seurat_cluster = as.character(Seurat_cluster)) %>%
  mutate(id = as.character(row_number())) %>%
  ggplot(data = ., aes(x = reorder(id,-signalValue), y = signalValue)) +
  geom_point(
    stat = 'identity',
    aes(color = H3K27ac_VLMC.bw),
    size = 3,
    alpha = 1
  ) +
  scale_color_gradient2(
    midpoint = mean(collect$H3K27ac_VLMC.bw),
    low = "blue",
    mid = "white",
    high = "red",
    space = "Lab",
    limits = c(0, 15)
  ) +
  labs(
    title = expression(paste("MACS2 peak scores")),
    x = "G4 peaks above 75th percentile",
    y = "MACS2 signalValue",
    color = "H3K27ac - VLMC"
  ) +
  # geom_label_repel(aes(label = ifelse(Cruz_Molina_enh == 1 | Glaser_enh == 1, as.character(Gene_name), '')),
  #                  box.padding   = 0.9,
  #                  max.overlaps = Inf,
  #                  point.padding = 0.5,
  #                  size = 2,
  #                  segment.color = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
peak_rank_ac_vlmc

peak_rank_grid = plot_grid(
  peak_rank_ac_astro,
  peak_rank_ac_mol,
  peak_rank_ac_oec,
  peak_rank_ac_opc,
  peak_rank_ac_vlmc
)
peak_rank_grid

ggsave(
  glue("{enhancer_anal}G4_peak_score_rank_grid.png"),
  plot = peak_rank_grid,
  width = 10,
  height = 5,
  dpi = 300,
)

scatter_astro = collect %>% mutate(LTR_vicinity = as.character(LTR_vicinity)) %>%
  ggplot(.,
         aes(x = H3K27ac_Astrocytes.bw, y = signalValue)) +
  geom_point(aes(fill = LTR_vicinity, size = LTR_vicinity), 
             colour = "black", pch = 21, size = 5, alpha = 1) +
  scale_fill_manual(values = c('#f0f0f0','#de2d26')) +
  scale_size_manual(values = c(3, 8)) +
  geom_label_repel(
    aes(label = ifelse(
      LTR_vicinity == 1,
      as.character(Gene_name),
      ''
    )),
    box.padding   = 0.9,
    max.overlaps = Inf,
    point.padding = 0.5,
    size = 4,
    segment.color = 'black'
  ) +
  labs(
    title = "",
    x = "H3K27ac - Astrocytes",
    y = "MACS2 signalValue",
    fill = "LTR vicinity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
scatter_astro

scatter_mol = collect %>% mutate(LTR_vicinity = as.character(LTR_vicinity)) %>%
  ggplot(.,
         aes(x = H3K27ac_mOL.bw, y = signalValue)) +
  geom_point(aes(fill = LTR_vicinity), 
             colour = "black", pch = 21, size = 5, alpha = 1) +
  scale_fill_manual(values = c('#f0f0f0','#de2d26')) +
  geom_label_repel(
    aes(label = ifelse(
      LTR_vicinity == 1,
      as.character(Gene_name),
      ''
    )),
    box.padding   = 0.9,
    max.overlaps = Inf,
    point.padding = 0.5,
    size = 4,
    segment.color = 'black'
  ) +
  labs(
    title = "",
    x = "H3K27ac - mOL",
    y = "MACS2 signalValue",
    fill = "LTR vicinity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
scatter_mol

scatter_oec = collect %>% mutate(LTR_vicinity = as.character(LTR_vicinity)) %>%
  ggplot(.,
         aes(x = H3K27ac_OEC.bw, y = signalValue)) +
  geom_point(aes(fill = LTR_vicinity), 
             colour = "black", pch = 21, size = 5, alpha = 1) +
  scale_fill_manual(values = c('#f0f0f0','#de2d26')) +
  geom_label_repel(
    aes(label = ifelse(
      LTR_vicinity == 1,
      as.character(Gene_name),
      ''
    )),
    box.padding   = 0.9,
    max.overlaps = Inf,
    point.padding = 0.5,
    size = 4,
    segment.color = 'black'
  ) +
  labs(
    title = "",
    x = "H3K27ac - OEC",
    y = "MACS2 signalValue",
    fill = "LTR vicinity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
scatter_oec

scatter_opc = collect %>% mutate(LTR_vicinity = as.character(LTR_vicinity)) %>%
  ggplot(.,
         aes(x = H3K27ac_OPC.bw, y = signalValue)) +
  geom_point(aes(fill = LTR_vicinity), 
             colour = "black", pch = 21, size = 5, alpha = 1) +
  scale_fill_manual(values = c('#f0f0f0','#de2d26')) +
  geom_label_repel(
    aes(label = ifelse(
      LTR_vicinity == 1,
      as.character(Gene_name),
      ''
    )),
    box.padding   = 0.9,
    max.overlaps = Inf,
    point.padding = 0.5,
    size = 4,
    segment.color = 'black'
  ) +
  labs(
    title = "",
    x = "H3K27ac - OPC",
    y = "MACS2 signalValue",
    fill = "LTR vicinity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )
scatter_opc

scatter_grid = plot_grid(
  scatter_astro,
  scatter_mol,
  scatter_oec,
  scatter_opc
)
scatter_grid

ggsave(
  glue("{enhancer_anal}G4_H3K27ac_scatter_grid.png"),
  plot = scatter_grid,
  width = 10,
  height = 10,
  dpi = 300,
)
