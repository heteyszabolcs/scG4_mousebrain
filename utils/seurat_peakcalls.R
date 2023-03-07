# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("EnsDb.Mmusculus.v79")
  library("GenomicFeatures")
  library("GenomeInfoDb")
  library("Matrix")
})

## seurat peak caller function: 
# subset sc sorted bam file based on a list of barcodes
# returns a bam file

# path to result folder
seurat_peak_caller = function(seurat_object,
                              initial_table,
                              barcode_column,
                              feature_column,
                              peak_name,
                              result_folder,
                              output) {
  
  # Seurat object
  seurat_object = readRDS(seurat_object)
  # initial table
  initial_table = read_tsv(initial_table)
  # get barcodes from initial table
  barcodes = initial_table %>% pull(barcode_column)
    
  # modify metadata
  meta = seurat_object@meta.data
  rows_as_barcodes = rownames(meta)
  meta = meta %>% mutate(barcode = rows_as_barcodes) %>% 
    mutate(peakcall_group = ifelse(rows_as_barcodes %in% barcodes, peak_name, "rest"))
  seurat_object@meta.data = meta
  
  # write out barcodes (it's necessary for sinto command)
  barcodes = initial_table %>% dplyr::select(barcode_column, feature_column)
  write_tsv(barcodes, glue("{result_folder}barcodes.tsv"), col_names = FALSE)
  
  # peak calling per group of interest (MACS2)
  peaks = CallPeaks(
    object = seurat_object,
    group.by = "peakcall_group",
    cleanup = FALSE,
    outdir = result_folder,
    effective.genome.size = 2652783500
  )
  
  # write out peak files
  write.table(
    as.data.frame(peaks),
    file = glue("{result_folder}{output}"),
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
  
  return(peaks)
  
}

# test
seurat_peak_caller(
  seurat_object = "../results/Seurat/GFP_sorted.Rds",
  initial_table = "../results/Seurat/callpeaks_GFPsorted/high_pred_MOL.tsv",
  barcode_column = "barcode",
  feature_column = "cell_type",
  peak_name = "high_pred_MOL",
  result_folder = "../results/Seurat/callpeaks_GFPsorted/",
  output = "peaks_per_MOL_preds.bed"
)














