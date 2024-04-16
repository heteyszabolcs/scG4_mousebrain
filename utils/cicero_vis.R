suppressMessages({
  library("GenomicRanges")
  library("Seurat")
  library("Signac")
  library("ggplot2")
  library("data.table")
  library("tidyverse")
  library("glue")
})

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
barcodes_scbridge = pred %>% filter(Prediction == "MOL") %>% pull(V1)
other = setdiff(colnames(g4@assays$GA@counts), barcodes_scbridge)

g4@meta.data = g4@meta.data %>%
  mutate(rownames_to_column(., var = "cell_id")) %>%
  mutate(MOL_status = ifelse(cell_id %in% barcodes_scbridge, "MOL", "non-MOL"))

# cicero outputs
load("../results/Seurat/callpeaks_GFPsorted/sorted_cicero-pred_nonMOL.Rds")
load(
  "../results/Seurat/callpeaks_GFPsorted/sorted_cicero_coG4networks-pred_nonMOL.Rds"
)
nonMOL = conns
nonMOL_ccan = CCAN_assigns
rm(conns)
load("../results/Seurat/callpeaks_GFPsorted/sorted_cicero-predMOL.Rds")
load("../results/Seurat/callpeaks_GFPsorted/sorted_cicero_coG4networks-predMOL.Rds")
MOL = conns
MOL_ccan = pred_mol
rm(conns)

nonMOL_x = nonMOL %>%
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

nonMOL_x$type = "non-MOL"
nonMOL_peakset = GRanges(
  seqnames = nonMOL_x$chr1,
  ranges = IRanges(
    start = nonMOL_x$start1,
    end = nonMOL_x$end2,
    names = nonMOL_x$type,
  )
)
MOL_x = MOL %>%
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

MOL_x$type = "MOL"
MOL_peakset = GRanges(
  seqnames = MOL_x$chr1,
  ranges = IRanges(
    start = MOL_x$start1,
    end = MOL_x$end2,
    names = MOL_x$type,
  )
)

# find MOL spec G4 regions
ol = findOverlaps(MOL_peakset,
                  nonMOL_peakset,
                  type = "any",
                  ignore.strand = FALSE)

MOL_spec = as_tibble(MOL_peakset[-queryHits(ol)])
MOL_spec_regions = MOL_spec %>% mutate(region = paste(seqnames, start, end, sep = "-")) %>%
  pull(region)

# find MOL spec G4 regions over cCREs
ol_w_enh = findOverlaps(MOL_peakset[-queryHits(ol)],
                        li_enh,
                        type = "any",
                        ignore.strand = FALSE)
MOL_spec_enh_regions = as_tibble(MOL_peakset[-queryHits(ol)][queryHits(ol_w_enh)])
MOL_spec_enh_regions = distinct_all(MOL_spec_enh_regions) %>%
  mutate(region = paste(seqnames, start, end, sep = "-")) %>%
  pull(region)


# make genome browser figures
pred_mol_links = ConnectionsToLinks(conns = MOL, ccans = MOL_ccan)
Links(g4) = pred_mol_links

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
    group.by = "MOL_status",
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

#lapply(MOL_spec_enh_regions, make_coverage_plot)

pred_nonmol_links = ConnectionsToLinks(conns = nonMOL, ccans = nonMOL_ccan)
Links(g4) = pred_nonmol_links
lapply(MOL_spec_enh_regions, make_coverage_plot)

# examples (cherry picking)
# Pth1r
pred_mol_links = ConnectionsToLinks(conns = MOL, ccans = MOL_ccan)
Links(g4) = pred_mol_links

mol_pth1r = make_coverage_plot(region = "chr9-110716521-110759131")

pred_nonmol_links = ConnectionsToLinks(conns = nonMOL, ccans = nonMOL_ccan)
Links(g4) = pred_nonmol_links

nonmol_pth1r = make_coverage_plot(region = "chr9-110716521-110759131")

pth1r_plots = ggarrange(
  plotlist = list(nonmol_pth1r, mol_pth1r),
  ncol = 2,
  nrow = 1
)

ggsave(
  glue("{result_folder}Pth1r_cicero_covplots.pdf"),
  plot = pth1r_plots,
  width = 12,
  height = 8,
  device = "pdf"
)

#
pred_mol_links = ConnectionsToLinks(conns = MOL, ccans = MOL_ccan)
Links(g4) = pred_mol_links
mol_vw1 = make_coverage_plot(region = "chr6-125602999-125605170")

pred_nonmol_links = ConnectionsToLinks(conns = nonMOL, ccans = nonMOL_ccan)
Links(g4) = pred_nonmol_links
nonmol_vw1 = make_coverage_plot(region = "chr6-125602999-125605170")

vw1_plots = ggarrange(
  plotlist = list(nonmol_vw1, mol_vw1),
  ncol = 2,
  nrow = 1
)

ggsave(
  glue("{result_folder}Vw1_cicero_covplots.pdf"),
  plot = vw1_plots,
  width = 12,
  height = 8,
  device = "pdf"
)

#
pred_mol_links = ConnectionsToLinks(conns = MOL, ccans = MOL_ccan)
Links(g4) = pred_mol_links
mol_elavl2 = make_coverage_plot(region = "chr4-91376129-91625283")

pred_nonmol_links = ConnectionsToLinks(conns = nonMOL, ccans = nonMOL_ccan)
Links(g4) = pred_nonmol_links
nonmol_elavl2 = make_coverage_plot(region = "chr4-91376129-91625283")

elavl2_plots = ggarrange(
  plotlist = list(nonmol_elavl2, mol_elavl2),
  ncol = 2,
  nrow = 1
)

ggsave(
  glue("{result_folder}Elavl2_cicero_covplots.pdf"),
  plot = elavl2_plots,
  width = 12,
  height = 8,
  device = "pdf"
)

# chr19
pred_mol_links = ConnectionsToLinks(conns = MOL, ccans = MOL_ccan)
Links(g4) = pred_mol_links
mol_chr19 = make_coverage_plot(region = "chr19-8897652-8993498")

pred_nonmol_links = ConnectionsToLinks(conns = nonMOL, ccans = nonMOL_ccan)
Links(g4) = pred_nonmol_links
nonmol_chr19 = make_coverage_plot(region = "chr19-8897652-8993498")

chr19_plots = ggarrange(
  plotlist = list(nonmol_chr19, mol_chr19),
  ncol = 2,
  nrow = 1
)

#
pred_mol_links = ConnectionsToLinks(conns = MOL, ccans = MOL_ccan)
Links(g4) = pred_mol_links
mol_wdr7 = make_coverage_plot(region = "chr18-63764748-63771732")

pred_nonmol_links = ConnectionsToLinks(conns = nonMOL, ccans = nonMOL_ccan)
Links(g4) = pred_nonmol_links
nonmol_wdr7 = make_coverage_plot(region = "chr18-63764748-63771732")

wdr7_plots = ggarrange(
  plotlist = list(nonmol_wdr7, mol_wdr7),
  ncol = 2,
  nrow = 1
)

ggsave(
  glue("{result_folder}Wdr7_cicero_covplots.pdf"),
  plot = wdr7_plots,
  width = 12,
  height = 8,
  device = "pdf"
)

CoveragePlot(object = g4,
             region = "chr9-110711875-110712480")
CoveragePlot(
  seurat_object,
  region = "chr9-110711875-110712480",
  ranges = li_enh,
  annotation = TRUE,
  show.bulk = TRUE,
  group.by = "MOL_status",
  ranges.title = region,
  height.tracks = 20
)
