suppressMessages({
  library("cicero");
  library("Seurat");
  library("Signac");
  library("ggplot2");
  library("monocle3");
  library("data.table")
})

set.seed(42)

# cicero: https://www.bioconductor.org/packages/release/bioc/vignettes/cicero/inst/doc/website.html
# implementation: https://github.com/Castelo-Branco-lab/scCut-Tag_2020/blob/master/notebooks/H3K27ac/Cicero.Rmd

mm10.chromsize = fread("../data/mm10.chrom.sizes.txt")
g4 = readRDS(file = "../results/Seurat/final/unsorted_brain/res0.8/outputs/Seurat_object.Rds")

input_cds = monocle3::new_cell_data_set(expression_data =  g4[['peaks']]@counts, cell_metadata = g4@meta.data,gene_metadata = data.frame('gene_short_name' = rownames(g4[['peaks']]),row.names=rownames(g4[['peaks']])))
input_cds = detect_genes(input_cds)
input_cds = input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

input_cds = detect_genes(input_cds)
input_cds = estimate_size_factors(input_cds)
input_cds = preprocess_cds(input_cds, method = "LSI")

# dimension reduction using LSI followed by UMAP
input_cds = reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")      


umap_coords = reducedDims(input_cds)$UMAP
cicero_cds = make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# estimate the co-G4 sites in the genome in order to predict cis-regulatory interactions
conns = run_cicero(cicero_cds, mm10.chromsize) 
save(conns, file = "../results/Seurat/callpeaks_unsorted/unsorted_cicero.Rds")

# finding cis-Coaccessibility networks (CCANS), coaccess_cutoff_override > 0.5 strict
CCAN_assigns = generate_ccans(conns, coaccess_cutoff_override = 0.1)
save(CCAN_assigns, file = "../results/Seurat/callpeaks_unsorted/unsorted_cicero_coG4networks.Rds")

# visualization
gene_anno = rtracklayer::readGFF("../data/gencode.vM10.annotation.gff3")

gene_anno$chromosome = paste0("chr", gene_anno$seqid)
gene_anno$gene = gene_anno$gene_id
gene_anno$transcript = gene_anno$transcript_id
gene_anno$symbol = gene_anno$gene_name

plot_connections(conns, "chr1", 106334061, 106404978,
                 gene_model = as.data.frame(gene_anno),
                 coaccess_cutoff = .10,
                 connection_width = .5)

