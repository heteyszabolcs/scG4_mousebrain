import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy.external as sce
import pandas as pd
import anndata as ad
from sklearn.metrics import accuracy_score, f1_score, silhouette_score

# folders
scbridge_folder = "C:/Szabolcs/Karolinska/Data/Jing/LTRIS2_BRG1_H33_G4s/results/scBridge/output/"
seurat_folder = "C:/Szabolcs/Karolinska/Data/Jing/LTRIS2_BRG1_H33_G4s/results/Seurat/final/sorted_brain/res0.8/integration/outputs/"
seurat_other = "C:/Szabolcs/Karolinska/Data/Jing/LTRIS2_BRG1_H33_G4s/results/Seurat/callpeaks_GFPsorted/"

# integrated output
g4_output = ad.read_h5ad(scbridge_folder + "scCutTag_gene_activity_scores-integrated.h5ad")
rna_output = ad.read_h5ad(scbridge_folder + "Marques_scRNA-Seq-integrated.h5ad")

seurat_preds = pd.read_csv(seurat_folder + "g4_cell_label_preds.tsv", sep='\t')
seurat_clusters = pd.read_csv(seurat_other + "Seurat_clusters_res0.8.tsv", sep = "\t")

# outcomes
g4_output_preds = g4_output.obs["Prediction"]
g4_output_rel = g4_output.obs["Reliability"]

g4_output_rel.to_csv(scbridge_folder + "scbridge_reliability.csv",index=True)
g4_output_preds.to_csv(scbridge_folder + "scbridge_predictions.csv", index=True)

### UMAPs ###
pastel = sns.color_palette("pastel")
palette = {}
types = ["COP", "MFOL", "MOL", "NFOL", "Novel (Most Unreliable)", "OPC", "PPR"]
for i in range(len(types)):
    if "Novel" in types[i]:
        palette[types[i]] = "#f0f0f0"
    else:
        palette[types[i]] = pastel[i]

# scBridge UMAps
sc.pl.umap(rna_output, color='CellType', palette = palette,  save = "Marques_CellTypes.pdf", show = True)
sc.pl.umap(g4_output, color='Prediction', palette = palette, save = "Predictions.pdf")
sc.pl.umap(g4_output, color='Reliability', save = "Reliabilities.pdf")

# Seurat prediction scores
g4_output.obs["Seurat_predscore"] = pd.Series(seurat_preds["Seurat_predscore"].tolist(),
          index=seurat_preds["id"].tolist())
sc.pl.umap(g4_output, color='Seurat_predscore',
           save = "Seurat_predscore.pdf")

# Seurat predicted labels
g4_output.obs["Seurat_predlabel"] = pd.Series(seurat_preds["Seurat_prediction"].tolist(),
          index=seurat_preds["id"].tolist())
sc.pl.umap(g4_output, color='Seurat_predlabel',
           save = "Seurat_predlabel.pdf", palette = palette)

# Seurat clusters
seurat_cl_strs = [str(i) for i in seurat_clusters["Seurat_cluster"].tolist()]

g4_output.obs["Seurat_cluster"] = pd.Series(seurat_cl_strs,
          index=seurat_clusters["id"].tolist())
domain_pal = {"0" : sns.color_palette("pastel")[0],
              "1" : sns.color_palette("pastel")[1],
              "2" : sns.color_palette("pastel")[2],
              "3" : sns.color_palette("pastel")[3]}
sc.pl.umap(g4_output, color='Seurat_cluster',
           save = "Seurat_cluster.pdf", palette = domain_pal)

comb = ad.concat([g4_output, rna_output])
domain_pal = {"Marques_scRNA-Seq" : sns.color_palette("pastel")[0],
              "scCutTag_gene_activity_scores" : sns.color_palette("pastel")[1]}
sc.pl.umap(comb, color='Domain', palette = domain_pal, save = "Domains.pdf")

comb.write('C:/Szabolcs/Karolinska/Data/Jing/LTRIS2_BRG1_H33_G4s/results/scBridge/output/combined.h5ad',
               compression="gzip")



