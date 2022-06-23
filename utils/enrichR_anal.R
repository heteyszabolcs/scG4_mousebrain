# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("enrichR")
  library("scales")
})

# export folder
result_folder = "../results/enrichment_analysis/"

gene_table = "../results/Seurat/callpeaks_GFPsorted/ArchRSignacLanc_promoter_G4_peaks_genes.tsv"
gene_table = fread(gene_table, header = FALSE)

# enrichR 
dbs = c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enriched = enrichr(gene_table$V1, dbs)
biol = enriched[["GO_Biological_Process_2015"]]
biol = enriched[["GO_Biological_Process_2015"]]
biol = biol %>% filter(0.01 > Adjusted.P.value) 

# horizontal bar
go_bars = ggplot(biol[1:15, ], aes(x = reorder(Term, Adjusted.P.value), y = Adjusted.P.value)) +
  geom_bar(stat = 'identity', fill = "#2ca25f", color = "black") +
  coord_flip() +
  scale_y_continuous(labels = scientific) +
  labs(title = "enrichR output of Promoter G4 peaks",
       y = "adj. p-value",
       x = "GO Biologocial process (2015)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars

# enrichr plots
plotEnrich(enriched[["GO_Biological_Process_2015"]], title = "GO Biological Process (2015)")
plotEnrich(enriched[["GO_Molecular_Function_2015"]], title = "GO Molecular Function (2015)")

# save
ggsave(
  glue("{result_folder}enrichR_GO_biol-promoterG4_peaks.png"),
  plot = go_bars,
  width = 10,
  height = 10,
  dpi = 300,
)
