import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc

g4_output = ad.read_h5ad('/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/output/scCutTag_gene_activity_scores-integrated.h5ad')
g4_output_preds = g4_output.obs["Prediction"]
g4_output_rel = g4_output.obs["Reliability"]

g4_output_rel.to_csv('/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/output/scbridge_reliability.csv',
                     index=True)
g4_output_preds.to_csv('/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/output/scbridge_predictions.csv',
                       index=True)

# UMAPs
rna_output = ad.read_h5ad('/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/output/Marques_scRNA-Seq-integrated.h5ad')
sc.pl.umap(rna_output, color='CellType', palette = "Set3", save = "Marques_CellTypes.pdf")

sc.pl.umap(g4_output, color='Prediction', palette = 'Set3', save = "Predictions.pdf")
sc.pl.umap(g4_output, color='Reliability', save = "Reliabilities.pdf")
sc.pl.umap(g4_output, color='Domain', save = "Reliabilities.pdf")
