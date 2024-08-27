if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("data.table",
               "tidyverse",
               "reshape2",
               "plyr",
               "ComplexHeatmap",
               "ggplot2",
               "viridis",
               "RColorBrewer",
               "circlize"
)


# Seurat object
# scG4 - scRNA integrated
seurat = readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/G4_scRNA_integration.Rds")
meta = seurat@meta.data

# scBridge outputs
rel = fread("../results/scBridge/output/scbridge_reliability.csv")
pred = fread("../results/scBridge/output/scbridge_predictions.csv", header = TRUE)

meta = meta %>% rownames_to_column(., var = "cell_id")
meta = meta %>% inner_join(., pred, by = c("cell_id" = "V1")) %>% 
  dplyr::rename(scBridge_prediction = Prediction)

# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2){
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}

# calculate overlap
# Seurat integration
ident2rna = data.frame(idents = rownames(meta), rna_label = meta$pred_cell_type)
ident2rna = ident2rna[complete.cases(ident2rna), ]

ident2g4 = data.frame(idents = rownames(meta), rna_label = meta$seurat_clusters)
ident2g4 = ident2g4[complete.cases(ident2g4), ]

ovlpScore.df = cal_ovlpScore(ident2rna, ident2g4)

mapSubclass = ovlpScore.df
colnames(mapSubclass) <- c("cell_type", "seurat_cluster", "ovlpScore")
mapSubclass = dcast(mapSubclass, cell_type~seurat_cluster, value.var = "ovlpScore", 
                     fun.aggregate = identity, fill = 0)
rows = mapSubclass$cell_type
mapSubclass = mapSubclass[,-1]
rownames(mapSubclass) = rows


col_fun = colorRamp2(c(0, 0.25, 0.5), c("#452258", "#679b81", "#f0e527"))

pdf(
   file = "../results/Seurat/Seurat_overlap_score_hm.pdf",
   width = 5,
   height = 4
 )
hm = Heatmap(
  mapSubclass,
  name = "overlap score",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "Seurat prediction",
  row_title_side = "right",
  column_title = "Seurat cluster",
  column_title_side = "bottom",
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  column_names_rot = 0
)
hm
dev.off()

# scBridge integration
meta = seurat@meta.data %>% mutate(cell_id = rownames(seurat@meta.data)) %>% 
  left_join(., rel, c("cell_id" = "V1")) %>% dplyr::rename(scBridge_reliability = "Reliability") %>% 
  left_join(., pred, c("cell_id" = "V1")) %>% dplyr::rename(scBridge_prediction = "Prediction") %>% 
  filter(scBridge_reliability > 0.9)

ident2rna = data.frame(idents = rownames(meta), rna_label = meta$scBridge_prediction)
ident2rna = ident2rna[complete.cases(ident2rna), ]
ident2rna = ident2rna[which(ident2rna$rna_label != "Novel (Most Unreliable)"), ]

ident2g4 = data.frame(idents = rownames(meta), rna_label = meta$seurat_clusters)
ident2g4 = ident2g4[complete.cases(ident2g4), ]

ovlpScore.df = cal_ovlpScore(ident2rna, ident2g4)

mapSubclass = ovlpScore.df
colnames(mapSubclass) <- c("cell_type", "seurat_cluster", "ovlpScore")
mapSubclass = dcast(mapSubclass, cell_type~seurat_cluster, value.var = "ovlpScore", 
                    fun.aggregate = identity, fill = 0)
rows = mapSubclass$cell_type
mapSubclass = mapSubclass[,-1]
rownames(mapSubclass) = rows


col_fun = colorRamp2(c(0, 0.25, 0.5), c("#452258", "#679b81", "#f0e527"))

pdf(
  file = "../results/Seurat/scBridge_overlap_score_hm.pdf",
  width = 5,
  height = 4
)
hm2 = Heatmap(
  mapSubclass,
  name = "overlap score",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "scBridge prediction",
  row_title_side = "right",
  column_title = "Seurat cluster",
  column_title_side = "bottom",
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  column_names_rot = 0
)
hm2
dev.off()

hm + hm2

# Seurat - scBridge overlap 
meta = meta %>% dplyr::filter(scBridge_prediction != "Novel (Most Unreliable)") %>% 
  dplyr::filter(pred_max_score > 0.5) %>% 
  mutate(scBridge_prediction = str_replace_all(scBridge_prediction, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(scBridge_prediction = str_replace_all(scBridge_prediction, pattern = "Oligodendrocytes", replacement = "MOL")) %>% 
  mutate(pred_cell_type = str_replace_all(pred_cell_type, pattern = "Astrocytes", replacement = "AST")) %>% 
  mutate(pred_cell_type = str_replace_all(pred_cell_type, pattern = "Oligodendrocytes", replacement = "MOL"))

ident2seurat = data.frame(idents = rownames(meta), rna_label = meta$pred_cell_type)
ident2seurat = ident2seurat[complete.cases(ident2seurat), ]

ident2sc_br = data.frame(idents = rownames(meta), rna_label = meta$scBridge_prediction)
ident2sc_br = ident2sc_br[complete.cases(ident2sc_br), ]

ovlpScore.df = cal_ovlpScore(ident2seurat, ident2sc_br)

mapSubclass = ovlpScore.df
colnames(mapSubclass) <- c("cell_type", "seurat_cluster", "ovlpScore")
mapSubclass = dcast(mapSubclass, cell_type~seurat_cluster, value.var = "ovlpScore", 
                    fun.aggregate = identity, fill = 0)
rows = mapSubclass$cell_type
mapSubclass = mapSubclass[,-1]
rownames(mapSubclass) = rows

pdf(
  file = "../results/Seurat/scBridge_Seurat_overlap_score_hm.pdf",
  width = 5,
  height = 5
)
hm3 = Heatmap(
  mapSubclass,
  name = "overlap score",
  clustering_distance_rows = "pearson",
  col = col_fun,
  row_title = "Seurat prediction (pred. score > 0.5)",
  row_title_side = "right",
  column_title = "scBridge prediction",
  column_title_side = "bottom",
  show_row_names = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  width = unit(7, "cm"),
  height = unit(3, "cm"),
  column_names_rot = 90
)
hm3
dev.off()
