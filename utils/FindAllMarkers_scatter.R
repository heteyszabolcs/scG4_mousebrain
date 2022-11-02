suppressPackageStartupMessages({
  library("data.table")
  library("tidyverse")
  library("ggplot2")
  library("ggrepel")
  library("glue")
  library("ggpubr")
})

# turn off warnings
options(warn = -1)

# folders
result_folder = "../results/Seurat/callpeaks_unsorted/"

# function for scatterplot
create_scatter =
  function(findallmarker_output, title) {
    prom = fread(findallmarker_output)
    sc = fread("../data/GSE75330/FindAllMarker_fc0.2_wilcox.tsv")
    
    comb = prom %>% left_join(., sc, by = c("Gene Name" = "gene")) %>%
      dplyr::select(
        `Gene Name`,
        G4_fc = "avg_log2FC.x",
        G4_padj = p_val_adj.x,
        "Distance to TSS",
        "scRNA_Seq_fc" = avg_log2FC.y,
        "scRNA_Seq_padj" =  p_val_adj.y,
        "scRNA_Seq_cluster" = cluster
      ) %>% drop_na()
    
    filt = comb %>% dplyr::filter(G4_padj < 0.05 &
                                    scRNA_Seq_padj < 0.05) %>%
      mutate(label = ifelse(abs(scRNA_Seq_fc) > 3, `Gene Name`, NA_character_))
    
    
    plot = ggplot(filt,
                  aes(
                    x = scRNA_Seq_fc,
                    y = G4_fc,
                    color = factor(scRNA_Seq_cluster),
                    label = label
                  )) +
      geom_point(size = 2) +
      scale_color_brewer(palette = "Set3") +
      labs(
        title = title,
        x = "Marques et al. scRNA-Seq log2FC",
        y = "G4 log2FC",
        color = " "
      ) +
      theme_classic() +
      labs(color = "Seurat cluster") +
      ylim(-1, 1) +
      xlim(-4, 4) +
      geom_text_repel(colour = "black", max.overlaps = 100) +
      theme(
        text = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")
      )
    return(print(plot))
    
  }

# loop over the promoter G4s that were marker regions in the Seurat clusters
# creating G4 log2FC vs. scRNA-Seq log2FC scatterplots
proms = list.files("../results/Seurat/callpeaks_unsorted/", pattern = "_promoters.tsv")
titles = c(
  "Promoter G4s - cluster 0",
  "Promoter G4s - cluster 1",
  "Promoter G4s - cluster 2",
  "Promoter G4s - cluster 3",
  "Promoter G4s - cluster 4",
  "Promoter G4s - cluster 5",
  "Promoter G4s - cluster 6"
)

for (i in 1:length(proms)) {
  print(create_scatter(
    findallmarker_output = glue("{result_folder}{proms[i]}"),
    title = titles[i]
  ))
  
  ggsave(
    plot = last_plot(),
    glue("{result_folder}{proms[i]}_scatter.pdf"),
    width = 6,
    height = 6
  )
  ggsave(
    plot = last_plot(),
    glue("{result_folder}{proms[i]}_scatter.png"),
    width = 6,
    height = 6,
    dpi = 300
  )
  
}

# grid 
plots = list()
for (i in 1:length(proms)) {
  plot = create_scatter(
    findallmarker_output = glue("{result_folder}{proms[i]}"),
    title = titles[i]
  )
  plots[[i]] = plot
}

ggarrange(plotlist = plots)

ggsave(
  plot = last_plot(),
  glue("{result_folder}FindAllMarker_proms_scatter.pdf"),
  width = 14,
  height = 12
)
ggsave(
  plot = last_plot(),
  glue("{result_folder}FindAllMarker_proms_scatter.png"),
  width = 12,
  height = 14,
  dpi = 300
)







