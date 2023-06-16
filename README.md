# Code

## Main steps includes: 

* Seurat processing of G scCut&Tag 
* generating cluster specific bam and bigwig files 
* integration with _Marques et al._ scRNA-Seq data 

## Workflow:

1. _seurat_workflow.R_:

parameters:

-c: path to CellRanger output directory containing the sorted bam files, fragments and annotation table

-w: path to working directory where the script will return the Seurat outputs (including Seurat object) and UMAPs

-r: a numeric value to indicate the resolution for Seurat clustering (FindClusters function)

This script returns the Seurat object to be used for the downstream steps, output of marker analysis, gene activity heatmap and UMAP plots 

2. _sinto.R_:

parameters:

-s: path the Seurat object coming from the _seurat_workflow.R_ script 

-w: path to working directory where the script will place the cluster specific barcodes and bam files

-b: path to sorted bam file computed by CellRanger workflow

sinto makes Seurat cluster specific bam files by its filterbarcode function. It requires a virtual environment with Python 3.8 <. sinto documentation can be found [here](https://timoast.github.io/sinto/). 

3. _bam_process.R_:

parameters:

-s: path the Seurat object coming from the _seurat_workflow.R_ script 

-w: path to working directory where the _sinto.R_ placed the cluster specific barcodes and bam files

Converts cluster specific bam files to RPGC normalized bigwigs and call cluster specific peaks. It depends on [MACS2](https://github.com/macs3-project/MACS) and [deeptools](https://deeptools.readthedocs.io/en/develop/).

4. _seurat_integration.R_:

parameters:

-s: path the Seurat object coming from the _seurat_workflow.R_ script 

-w: path to working directory where the script returns the integration outputs including anchor scores, label predictins, post-integration UMAPs etc.

It makes a multi-omics integration across GFP sorted G4 scCut&Tag and [_Marques et al._ oigodendrocyte scRNA-Seq](https://www.science.org/doi/10.1126/science.aaf6463) data using Seurat's FindTransferAnchors function. It attempts to label the G4 clusters by label imputation using TransferData function. 


# Data

_seurat_integration.R_ file requires scRNA-Seq count table of _Marques et al._ that can be found in the data folder or [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75330). 

All medium or high sized files needed for explorative analysis are available at Uppmax: 
  
* @ Uppmax:
    
external bigwigs: _/proj/snic2020-6-3/GEOPIPE/OUTPUT/bw_
    
our working directory: _/proj/snic2020-6-3/SZABOLCS/LTRIS2_BRG1_H33_and_G4s/

# Task and discussion
   
* @ Trello: 
    
[Marek's analysis](https://trello.com/c/c2a0wan6/12-mareks-analysis)

