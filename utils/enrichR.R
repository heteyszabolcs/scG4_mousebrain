# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ComplexHeatmap")
  library("circlize")
  library("enrichR")
  library("ggpubr")
})

# folders
result_folder = "../results/Seurat/callpeaks_unsorted/"

# source annotation script for mm10
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# unsorted G4 scCut&Tag, res 0.1
# FindAllMarker output
marker_regions = fread("../results/Seurat/final/unsorted_brain/res0.1/outputs/FindAllMarkers_logreg_output.tsv")
marker_regions

# FindMarker output - positive enrichments of both clusters
cl1_pos = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_logreg_res0.1-cluster1.tsv")
cl0_pos = fread("../results/Seurat/callpeaks_unsorted/FindAllMarkers_logreg_res0.1-cluster0.tsv")

# MACS2 peaks
peaks_1 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/1_peaks.narrowPeak")
peaks_1 = peaks_1 %>% dplyr::filter(V7 > quantile(V7, .75)) # keep high peaks
peaks_1 = mm10_annotation(regions = peaks_1, seqname_col = "V1", start_col = "V2", end_col = "V3")
promoter_peaks_1 = peaks_1 %>% dplyr::select(gene_symbol = SYMBOL, distanceToTSS) %>% dplyr::filter(abs(distanceToTSS) < 3000)

peaks_0 = fread("../results/Seurat/final/unsorted_brain/res0.1/cluster_spec_peaks/0_peaks.narrowPeak")
peaks_0 = peaks_0 %>% dplyr::filter(V7 > quantile(V7, .75)) # keep high peaks
peaks_0 = mm10_annotation(regions = peaks_0, seqname_col = "V1", start_col = "V2", end_col = "V3")
promoter_peaks_0 = peaks_0 %>% dplyr::select(gene_symbol = SYMBOL, distanceToTSS) %>% dplyr::filter(abs(distanceToTSS) < 3000)

annot = mm10_annotation(regions = marker_regions, seqname_col = "chr", start_col = "start", end_col = "end", feature_1 = "avg_log2FC", feature_2 = "p_val_adj", feature_3 = "cluster")
marker_regions = annot %>% dplyr::select(gene_symbol = SYMBOL, distanceToTSS, avg_log2FC = feature_1, p_val_adj = feature_2, cluster = feature_3)

## enrichR on peaks enriched only in cluster 0
# enrichr: https://maayanlab.cloud/Enrichr/
cluster0_input = marker_regions %>% 
  dplyr::filter(cluster == 0) %>% 
  dplyr::filter(abs(distanceToTSS) < 3000) %>% 
  dplyr::filter(avg_log2FC < -1)

# enrichR for cluster 0
dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
enriched_0 = enrichr(cluster0_input$gene_symbol, dbs)
biol = enriched_0[["GO_Biological_Process_2023"]]
biol = biol %>% filter(0.05 > P.value) 

biol_proc_0 = enriched_0$GO_Biological_Process_2023
biol_proc_0 = biol_proc_0 %>% arrange(P.value)

# horizontal bar for biol. processes
go_bars_cl0 = ggplot(biol_proc_0[1:15, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity', fill = "#9ecae1", color = "black") +
  coord_flip() +
  labs(title = "promoter peaks enriched in cluster 0",
       y = "-log10(p-value)",
       x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_cl0

molfunc = enriched_0[["GO_Molecular_Function_2023"]]
molfunc = molfunc %>% filter(0.05 > P.value) 

molfunc_0 = enriched_0$GO_Molecular_Function_2023
molfunc_0 = molfunc_0 %>% arrange(P.value)

# horizontal bar for mol. functions
go_bars_molfunc_cl0 = ggplot(molfunc_0[1:15, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity', fill = "#9ecae1", color = "black") +
  coord_flip() +
  labs(title = "promoter peaks enriched in cluster 0",
       y = "-log10(p-value)",
       x = "GO Molecular function (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_molfunc_cl0

biol_molfunc_bars = ggarrange(go_bars_cl0, go_bars_molfunc_cl0)
ggsave(
  glue("{result_folder}enrichR_GO-cluster0_enriched_prom_peaks.png"),
  plot = biol_molfunc_bars,
  width = 25,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}enrichR_GO-cluster0_enriched_prom_peaks.pdf"),
  plot = biol_molfunc_bars,
  width = 25,
  height = 10
)

## enrichR on promoter G4s coming from the Seurat's FindMarker analysis
# for cluster 0
dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
enriched_fma_0 = enrichr(cl0_pos$gene_symbol, dbs)
biol = enriched_fma_0[["GO_Biological_Process_2023"]]
biol = biol %>% filter(0.05 > P.value) 

biol_proc_0 = enriched_fma_0$GO_Biological_Process_2023
biol_proc_0 = biol_proc_0 %>% arrange(P.value)

# horizontal bar for biol. processes
go_bars_fma_cl0 = ggplot(biol_proc_0[1:15, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity', fill = "#deebf7", color = "black") +
  coord_flip() +
  labs(title = "promoter peaks enriched in cluster 0",
       y = "-log10(p-value)",
       x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_fma_cl0

molfunc = enriched_fma_0[["GO_Molecular_Function_2023"]]
molfunc = molfunc %>% filter(0.05 > P.value) 

molfunc_0 = enriched_fma_0$GO_Molecular_Function_2023
molfunc_0 = molfunc_0 %>% arrange(P.value)

# horizontal bar for mol. functions
go_bars_molfunc_fma_cl0 = ggplot(molfunc_0[1:15, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity', fill = "#deebf7", color = "black") +
  coord_flip() +
  labs(title = "promoter peaks enriched in cluster 0",
       y = "-log10(p-value)",
       x = "GO Molecular function (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_molfunc_fma_cl0

biol_fma_bars = ggarrange(go_bars_fma_cl0, go_bars_molfunc_fma_cl0)
ggsave(
  glue("{result_folder}enrichR_GO-cluster0_enriched_prom_peaks-findmarker.png"),
  plot = biol_fma_bars,
  width = 25,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}enrichR_GO-cluster0_enriched_prom_peaks-findmarker.pdf"),
  plot = biol_fma_bars,
  width = 25,
  height = 10
)

# for cluster 1
dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
enriched_fma_1 = enrichr(cl1_pos$gene_symbol, dbs)
biol = enriched_fma_1[["GO_Biological_Process_2023"]]
biol = biol %>% filter(0.05 > P.value) 

biol_proc_1 = enriched_fma_1$GO_Biological_Process_2023
biol_proc_1 = biol_proc_1 %>% arrange(P.value)

# horizontal bar for biol. processes
go_bars_fma_cl1 = ggplot(biol_proc_1[1:15, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity', fill = "#9ecae1", color = "black") +
  coord_flip() +
  labs(title = "promoter peaks enriched in cluster 1",
       y = "-log10(p-value)",
       x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_fma_cl1

molfunc = enriched_fma_1[["GO_Molecular_Function_2023"]]
molfunc = molfunc %>% filter(0.05 > P.value) 

molfunc_1 = enriched_fma_1$GO_Molecular_Function_2023
molfunc_1 = molfunc_1 %>% arrange(P.value)

# horizontal bar for mol. functions
go_bars_molfunc_fma_cl1 = ggplot(molfunc_1[1:15, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity', fill = "#9ecae1", color = "black") +
  coord_flip() +
  labs(title = "promoter peaks enriched in cluster 1",
       y = "-log10(p-value)",
       x = "GO Molecular function (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_molfunc_fma_cl1

biol_fma_bars = ggarrange(go_bars_fma_cl1, go_bars_molfunc_fma_cl1)
ggsave(
  glue("{result_folder}enrichR_GO-cluster1_enriched_prom_peaks-findmarker.png"),
  plot = biol_fma_bars,
  width = 25,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}enrichR_GO-cluster1_enriched_prom_peaks-findmarker.pdf"),
  plot = biol_fma_bars,
  width = 25,
  height = 10
)

## enrichR on high promoter proximal peaks - MACS2 based
# enrichR for cluster 0
dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
enriched_promoter_peaks_0 = enrichr(promoter_peaks_0$gene_symbol, dbs)
biol_promoter_peaks_0 = enriched_promoter_peaks_0[["GO_Biological_Process_2023"]]
biol_promoter_peaks_0 = biol_promoter_peaks_0 %>% filter(0.05 > Adjusted.P.value) %>% 
  top_n(-15, wt = Adjusted.P.value)
molfunc_promoter_peaks_0 = enriched_promoter_peaks_0[["GO_Molecular_Function_2023"]]
molfunc_promoter_peaks_0 = molfunc_promoter_peaks_0 %>% filter(0.05 > Adjusted.P.value) %>% 
  top_n(-15, wt = Adjusted.P.value)

# horizontal bar
go_bars_prom_peaks_cl0 = ggplot(biol_promoter_peaks_0, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = "#deebf7", color = "black") +
  coord_flip() +
  labs(title = "cluster 0 - promoter proximal G4 peaks",
       y = "-log10(adj. p-value)",
       x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_prom_peaks_cl0

go_bars_prom_peaks_molfunc_cl0 = ggplot(molfunc_promoter_peaks_0, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = "#deebf7", color = "black") +
  coord_flip() +
  labs(title = "cluster 0 - promoter proximal G4 peaks",
       y = "-log10(adj. p-value)",
       x = "GO Molecular function (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_prom_peaks_molfunc_cl0

# enrichR for cluster 1
enriched_promoter_peaks_1 = enrichr(promoter_peaks_1$gene_symbol, dbs)
biol_promoter_peaks_1 = enriched_promoter_peaks_1[["GO_Biological_Process_2023"]]
biol_promoter_peaks_1 = biol_promoter_peaks_1 %>% filter(0.05 > Adjusted.P.value) %>% 
  top_n(-15, wt = Adjusted.P.value)
molfunc_promoter_peaks_1 = enriched_promoter_peaks_1[["GO_Molecular_Function_2023"]]
molfunc_promoter_peaks_1 = molfunc_promoter_peaks_1 %>% filter(0.05 > Adjusted.P.value) %>% 
  top_n(-15, wt = Adjusted.P.value)

go_bars_prom_peaks_cl1 = ggplot(biol_promoter_peaks_1, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = "#deebf7", color = "black") +
  coord_flip() +
  labs(title = "cluster 1 - promoter proximal G4 peaks",
       y = "-log10(adj. p-value)",
       x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_prom_peaks_cl1

go_bars_prom_peaks_molfunc_cl1 = ggplot(molfunc_promoter_peaks_1, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = 'identity', fill = "#deebf7", color = "black") +
  coord_flip() +
  labs(title = "cluster 1 - promoter proximal G4 peaks",
       y = "-log10(adj. p-value)",
       x = "GO Molecular function (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
go_bars_prom_peaks_molfunc_cl1

biol_bars = ggarrange(go_bars_prom_peaks_cl0, go_bars_prom_peaks_cl1)
ggsave(
  glue("{result_folder}enrichR_GO_biol-promoterG4_peaks.png"),
  plot = biol_bars,
  width = 25,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}enrichR_GO_biol-promoterG4_peaks.pdf"),
  plot = biol_bars,
  width = 25,
  height = 10
)

molfunc_bars = ggarrange(go_bars_prom_peaks_molfunc_cl0, go_bars_prom_peaks_molfunc_cl1)
ggsave(
  glue("{result_folder}enrichR_GO_Mol_func-promoterG4_peaks.png"),
  plot = molfunc_bars,
  width = 25,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}enrichR_GO_Mol_func-promoterG4_peaks.pdf"),
  plot = molfunc_bars,
  width = 25,
  height = 10
)


