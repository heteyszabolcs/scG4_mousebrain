
coembed.g4 = coembed[,coembed$data_type == "sorted G4 scCut&Tag"]
coembed.scrna = coembed[,coembed$data_type == "scRNA-Seq (Marques et al.)"]
max(coembed.g4@assays$GA["Ptprz1",])
max(coembed.scrna@assays$GA["Ptprz1",])

x = WhichCells(coembed.g4, expression = Ptprz1 > 1.1)
x

coembed.g4@assays$RNA


DimPlot(object = coembed.g4, cells.highlight = x, cols.highlight = "red", cols = "gray", order = TRUE)


predicted_labels_gpr17 = predictions_long %>% dplyr::filter(cell_id %in% x)

x = sorted@meta.data

ggplot(predicted_labels_gpr17,
       aes(x = types, y = pred_score, fill = types)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Set3") +
  ylim(0, 1) +
  labs(title = "",
       x = "",
       y = "",
       fill = "") +
  theme_classic() +
  guides(fill = "none") +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )
