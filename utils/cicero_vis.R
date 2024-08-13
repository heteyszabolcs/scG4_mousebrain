suppressMessages({
  library("GenomicRanges")
  library("Seurat")
  library("Signac")
  library("ggplot2")
  library("data.table")
  library("tidyverse")
  library("glue")
  library("ggpubr")
})

set.seed(42)

# result folder
result_folder = "../results/genome_browser/Figure_4/"

# K27ac-ed mouse brain cCREs (Li et al.)
li_enh = readRDS("../results/GenomicRanges/Li_et_al-mousebrain.union.cCRE_with_K27ac.Rds")

# all putative mouse brain cCREs (Li et al.)
# li_enh = "../data/bed/Li_et_al-mousebrain.union.cCRE.bed"
# li_enh = fread(li_enh)
# li_enh$V5 = "Li_enh"
# li_enh = GRanges(
#   seqnames = li_enh$V1,
#   ranges = IRanges(
#     start = li_enh$V2,
#     end = li_enh$V3,
#     names = li_enh$V5
#   )
# )

# G4 scCut&Tag GFP sorted
g4 = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")

# scBridge outputs
pred = fread("../results/scBridge/output/scbridge_predictions.csv",
             header = TRUE)
barcodes_scbridge = pred %>% filter(Prediction == "Astrocytes") %>% pull(V1)
other = setdiff(colnames(g4@assays$GA@counts), barcodes_scbridge)

g4@meta.data = g4@meta.data %>%
  mutate(rownames_to_column(., var = "cell_id")) %>%
  mutate(AST_status = ifelse(cell_id %in% barcodes_scbridge, "AST", "non-AST"))

# cicero outputs
load("../results/Seurat/callpeaks_GFPsorted/sorted_cicero-pred_nonAST.Rds")
load(
  "../results/Seurat/callpeaks_GFPsorted/sorted_cicero_coG4networks-pred_nonAST.Rds"
)
nonAST = conns
nonAST_ccan = CCAN_assigns
rm(conns)
load("../results/Seurat/callpeaks_GFPsorted/sorted_cicero-predAST.Rds")
load("../results/Seurat/callpeaks_GFPsorted/sorted_cicero_coG4networks-predAST.Rds")
AST = conns
AST_ccan = CCAN_assigns
rm(conns)

nonAST_x = nonAST %>%
  dplyr::filter(coaccess > 0.5) %>%
  separate(Peak1,
           sep = "-",
           into = c("chr1", "start1", "end1")) %>%
  separate(Peak2,
           sep = "-",
           into = c("chr2", "start2", "end2")) %>%
  dplyr::select(chr1, start1, end2) %>%
  mutate(start1 = as.numeric(start1), end2 = as.numeric(end2)) %>%
  dplyr::filter(start1 < end2)

nonAST_x$type = "non-AST"
nonAST_peakset = GRanges(
  seqnames = nonAST_x$chr1,
  ranges = IRanges(
    start = nonAST_x$start1,
    end = nonAST_x$end2,
    names = nonAST_x$type,
  )
)
AST_x = AST %>%
  dplyr::filter(coaccess > 0.5) %>%
  separate(Peak1,
           sep = "-",
           into = c("chr1", "start1", "end1")) %>%
  separate(Peak2,
           sep = "-",
           into = c("chr2", "start2", "end2")) %>%
  dplyr::select(chr1, start1, end2) %>%
  mutate(start1 = as.numeric(start1), end2 = as.numeric(end2)) %>%
  dplyr::filter(start1 < end2)

AST_x$type = "AST"
AST_peakset = GRanges(
  seqnames = AST_x$chr1,
  ranges = IRanges(
    start = AST_x$start1,
    end = AST_x$end2,
    names = AST_x$type,
  )
)

# find AST spec G4 regions
ol = findOverlaps(AST_peakset,
                  nonAST_peakset,
                  type = "any",
                  ignore.strand = FALSE)

AST_spec = as_tibble(AST_peakset[-queryHits(ol)])
AST_spec_regions = AST_spec %>% mutate(region = paste(seqnames, start, end, sep = "-")) %>%
  pull(region)

# find AST spec G4 regions over cCREs
ol_w_enh = findOverlaps(AST_peakset[-queryHits(ol)],
                        li_enh,
                        type = "any",
                        ignore.strand = FALSE)
AST_spec_enh_regions = as_tibble(AST_peakset[-queryHits(ol)][queryHits(ol_w_enh)])
AST_spec_enh_regions = distinct_all(AST_spec_enh_regions) %>%
  mutate(region = paste(seqnames, start, end, sep = "-")) %>%
  pull(region)


# make genome browser figures
pred_AST_links = ConnectionsToLinks(conns = AST, ccans = AST_ccan)
Links(g4) = pred_AST_links

make_coverage_plot = function(region,
                              distance = 5000,
                              seurat_object = g4) {
  chr = strsplit(region, "-")[[1]][1]
  start = as.numeric(strsplit(region, "-")[[1]][2]) - distance
  end = as.numeric(strsplit(region, "-")[[1]][3]) + distance
  p = CoveragePlot(
    seurat_object,
    region = glue("{chr}-{start}-{end}"),
    annotation = TRUE,
    show.bulk = TRUE,
    group.by = "AST_status",
    ranges.title = region,
    ranges = li_enh,
    height.tracks = 20
  )
  
  annot = AnnotationPlot(object = g4,
                         region = region)
  ggsave(
    glue("{result_folder}{region}_cicero_covplots.pdf"),
    plot = p,
    width = 12,
    height = 8,
    device = "pdf"
  )
  ggsave(
    glue("{result_folder}{region}_annotplot.pdf"),
    plot = annot,
    width = 12,
    height = 8,
    device = "pdf"
  )
  return(print(p))
}

make_coverage_plot_all = function(region,
                                  distance = 5000) {
  
  g4 = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")
  g4@meta.data = g4@meta.data %>%
    mutate(rownames_to_column(., var = "cell_id")) %>%
    mutate(AST_status = ifelse(cell_id %in% barcodes_scbridge, "AST", "non-AST"))
  
  pred_AST_links = ConnectionsToLinks(conns = AST, ccans = AST_ccan)
  Links(g4) = pred_AST_links
  AST_example = make_coverage_plot(region = region, seurat_object = g4)
  
  g4 = readRDS(file = "../results/Seurat/callpeaks_GFPsorted/GFPsorted.Rds")
  g4@meta.data = g4@meta.data %>%
    mutate(rownames_to_column(., var = "cell_id")) %>%
    mutate(AST_status = ifelse(cell_id %in% barcodes_scbridge, "AST", "non-AST"))
  pred_nonAST_links = ConnectionsToLinks(conns = nonAST, ccans = nonAST_ccan)
  Links(g4) = pred_nonAST_links
  nonAST_example = make_coverage_plot(region = region, seurat_object = g4)
  
  plots = ggarrange(
    plotlist = list(nonAST_example, AST_example),
    ncol = 2,
    nrow = 1
  )
  
  ggsave(
    glue("{result_folder}{region}_cicero_arranged_covplots.pdf"),
    plot = plots,
    width = 12,
    height = 8,
    device = "pdf"
  )
  
}

selected_enhs = sample(AST_spec_enh_regions, 50)
lapply(selected_enhs, make_coverage_plot)
pred_AST_links = ConnectionsToLinks(conns = AST, ccans = AST_ccan)
Links(g4) = pred_AST_links
lapply(selected_enhs, make_coverage_plot)

# make ggarrange CoveragePlots
lapply(selected_enhs, make_coverage_plot_all)



