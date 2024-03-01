import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
print(ad.__version__)

# scRNA annotations
annot = pd.read_csv("/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/GSE75330_Marques_et_al_annot.csv",
                    index_col = 0)
annot.loc[annot['CellType'].str.contains('MOL'), 'CellType'] = 'MOL'
annot.loc[annot['CellType'].str.contains('MFOL'), 'CellType'] = 'MFOL'
annot.loc[annot['CellType'].str.contains('NFOL'), 'CellType'] = 'NFOL'

# process scRNA
counts = pd.read_csv("/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/GSE75330_Marques_et_al_mol_counts2.csv", index_col=0)
counts = counts.transpose()
counts_mat = csr_matrix(counts)
adata = ad.AnnData(counts_mat)
adata.obs_names = list(counts.index)
adata.var_names = list(counts.columns.values)
rna_genes = counts.columns.to_list()
adata.obs["CellType"] = annot
adata.obs["Domain"] = "RNA"

# process G4 scCut&Tag GA scores
ga = pd.read_csv("/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/gene_activity_scores.csv", index_col=0)
ga = ga.transpose()
ga_mat = csr_matrix(ga)
ga_adata = ad.AnnData(ga_mat)
ga_adata.obs_names = list(ga.index)
ga_adata.var_names = list(ga.columns.values)
g4_genes = ga.columns.to_list()
ga_adata.obs["Domain"] = "G4"

# export
ga_adata.write('/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/input/scCutTag_gene_activity_scores.h5ad',
               compression="gzip")
adata.write('/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/results/scBridge/input/Marques_scRNA-Seq.h5ad',
            compression="gzip")
