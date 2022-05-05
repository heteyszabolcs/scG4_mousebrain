# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("data.table")
  library("GenomicFeatures")
  library("GenomicRanges")
  library("wigglescout")
  library("plyranges")
  library("ArchR")
})

result_folder = "../results/GenomicRanges/"

# data
g4_peaks = "../results/Seurat/callpeaks_GFPsorted/peaks_per_clusters.bed"
cm_enh = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
gl_enh = "../data/bed/GSE171771_FAIRE_STARR_enh_mESC.bed"
ltrs = "../data/bed/RepMasker_lt200bp.LTRIS2.bed"

## data exploration on G4 peak set
# G4 peak set
g4_peaks = fread(g4_peaks)
g4_peaks$V7 = "G4"
# enhancers
cm_enh = fread(cm_enh)
cm_enh$V6 = "Cruz-Molina_active_enh"
cm_enh = GRanges(
  seqnames = cm_enh$V1,
  ranges = IRanges(
    start = cm_enh$V2,
    end = cm_enh$V3,
    names = cm_enh$V6
  )
)

gl_enh = fread(gl_enh)
gl_enh$V6 = "Glaser_active_enh"
gl_enh = GRanges(
  seqnames = gl_enh$V1,
  ranges = IRanges(
    start = gl_enh$V2,
    end = gl_enh$V3,
    names = gl_enh$V6
  )
)

# G4 length distributions
hist = g4_peaks %>% mutate(diff = V3 - V2) %>%
  filter(V6 %in% c("0", "1", "2", "3", "4")) %>%
  ggplot(., aes(x = diff, fill = V6)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "G4 length distributions / Seurat cluster",
       x = "length (bp)",
       y = "Density",
       fill = "Seurat cluster") +
  theme_classic() +
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, color = "black"))
hist

ggsave(
  glue("{result_folder}G4_length_distr_Seurat_clst.png"),
  plot = hist,
  width = 10,
  height = 10,
  dpi = 300,
)

# convert G4 peak ranges to GRanges
g4_peaks =
  lapply(split(g4_peaks, g4_peaks$V6), function(i) {
    GRanges(seqnames = i$V1,
            ranges = IRanges(
              start = i$V2,
              end = i$V3,
              names = i$V7
            ))
  })


## find overlaps between G4 peaks and active enhancers
# query: g4_peaks
# subject: enhancers
g4s_no_ol = g4_peaks[c("0", "1", "2", "3", "4")]
g4_enh_quant = numeric()

for (cluster in names(g4s_no_ol)) {
  ol = suppressWarnings(findOverlaps(g4_peaks[[cluster]], cm_enh, type = "any", minoverlap = 1))
  g4_enh_quant = c(g4_enh_quant, length(ol))
}

bar = tibble(
  overlap = g4_enh_quant,
  cluster = names(g4s_no_ol),
  fill_col = names(g4s_no_ol)
) %>%
  ggplot(data = ., aes(
    x = reorder(cluster, -overlap),
    y = overlap,
    fill = fill_col
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(
    title = expression(paste(
      "G4 overlaps with active enhancers of ",
      italic("Cruz-Molina et al.")
    )),
    x = "Seurat cluster",
    y = "# of G4 - active enhancer overlaps",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, color = "black"))
bar

ggsave(
  glue("{result_folder}G4_overlaps_w_CruzM_active_enh.png"),
  plot = bar,
  width = 10,
  height = 10,
  dpi = 300,
)

## HOMER annotation of G4 peaks
result_folder = "../results/Seurat/callpeaks_GFPsorted/"
bw_folder = "../data/GSE157637/"
reference = "mm10"

annotation = function(narrowpeak, percentile = 0.75) {
  cluster = strsplit(narrowpeak, ".narrowPeak")[[1]][1]
  print(cluster)
  np = fread(glue("{result_folder}{narrowpeak}"))
  
  # create narrowpeak table
  np = np %>% dplyr::select(
    chrom = V1,
    start = V2,
    end = V3,
    name = V4,
    score = V5,
    strand = V6,
    signalValue = V7,
    pValue = V8,
    qValue = V9,
    peak = V10
  ) %>%
    separate(name, sep = "_" , into = c(".", "..", "name")) %>%
    dplyr::select(-., -..) %>%
    mutate(start = as.character(start), end = as.character(end))
  
  # create HOMER input
  # keep peak scores with above 75th percentile
  np = np %>% filter(signalValue >= quantile(signalValue, percentile))
  bed_format = np %>% dplyr::select(chrom, start, end)
  write_tsv(bed_format,
            glue("{result_folder}{cluster}.bed"),
            col_names = FALSE)
  
  # annotate by HOMER
  bed_file = glue("{cluster}.bed")
  system(
    glue(
      "annotatePeaks.pl {result_folder}{bed_file} {reference} > {result_folder}{cluster}_annot.tsv"
    )
  )
  
  # read annotated file
  annot = fread(glue("{result_folder}{cluster}_annot.tsv"))
  id_col = colnames(annot)[1]
  
  # assign MACS2 columns
  annot = annot %>% mutate(Start = as.character(Start), End = as.character(End)) %>%
    dplyr::select("PeakID" = all_of(id_col), everything()) %>%
    mutate(PeakID = as.character(PeakID)) %>%
    inner_join(., np, by = c("PeakID" = "name")) %>% dplyr::select(Chr,
                                                                   Start,
                                                                   End,
                                                                   "Distance to TSS",
                                                                   "Gene Name",
                                                                   signalValue,
                                                                   pValue,
                                                                   qValue,
                                                                   peak) %>%
    
    mutate(Seurat_cluster = cluster) %>%
    separate(Seurat_cluster, sep = "_", into = "Seurat_cluster")
  
  # export
  write_tsv(annot, glue("{result_folder}{cluster}_annot.tsv"))
  
  return(annot)
  
}

# run on all narrowpeaks
narrowpeaks = list.files(glue("{result_folder}"), pattern = "*.narrowPeak")
lapply(narrowpeaks, annotation)

## complete annotated peaks with enhancer signatures
add_enhancer_status = function(annot_peak_file) {
  annot = fread(glue("{result_folder}{annot_peak_file}"))
  cluster = strsplit(annot_peak_file, split = "_peaks_annot.tsv")[[1]][1]
  
  # convert to GRanges object
  annot_gr = GRanges(
    seqnames = annot$Chr,
    ranges = IRanges(
      start = annot$Start,
      end = annot$End,
      names = annot$Seurat_cluster
    )
  )
  
  # intersect with Cruz-Molina enhancers
  cm_int = suppressWarnings(findOverlaps(annot_gr, cm_enh, type = "any", minoverlap = 1))
  annot = annot %>% mutate(Cruz_Molina_enh = "0")
  for (hit in cm_int@from) {
    annot[hit, which(colnames(annot) == "Cruz_Molina_enh")] = "1"
  }
  # intersect with Glaser enhancers
  gl_int = suppressWarnings(findOverlaps(annot_gr, gl_enh, type = "any", minoverlap = 1))
  annot = annot %>% mutate(Glaser_enh = "0")
  for (hit in gl_int@from) {
    annot[hit, which(colnames(annot) == "Glaser_enh")] = "1"
  }
  
  # add Bartosovic et al. scCut&Tag data
  require("wigglescout")
  bws = list.files(bw_folder, pattern = "*.bw", full.names = TRUE)
  bws = bws[grep("H3K27ac", bws)]
  labels = unname(sapply(bws, function(x)
    str_split(x, "/")[[1]][4]))
  
  # convert to GRanges object
  annot_bed = annot[, 1:3]
  annot_bed[["name"]] = cluster
  annot_bed = GRanges(
    seqnames = annot_bed$Chr,
    ranges = IRanges(
      start = annot_bed$Start,
      end = annot_bed$End,
      names = annot_bed$name
    )
  )
  
  # add metadata
  annot_bed$Distance_to_TSS = annot$`Distance to TSS`
  annot_bed$Gene_name = annot$`Gene Name`
  annot_bed$signalValue = annot$signalValue
  annot_bed$pValue = annot$pValue
  annot_bed$qvalue = annot$qValue
  annot_bed$peak = annot$peak
  annot_bed$Seurat_cluster = annot$Seurat_cluster
  annot_bed$Cruz_Molina_enh = annot$Cruz_Molina_enh
  annot_bed$Glaser_enh = annot$Glaser_enh
  
  aggr = suppressWarnings(bw_loci(bws, annot_bed, labels = labels))
  # left join by plyranges
  output = join_overlap_left_directed(annot_bed, aggr)
  output = as_tibble(output)
  
  # export
  write_tsv(output, glue("{result_folder}{cluster}_enh_anal.tsv"))
  
  return(output)
  
}

# run on all annotated peaks
annot_peaks = list.files(glue("{result_folder}"), pattern = "*peaks_annot.tsv")
lapply(annot_peaks, add_enhancer_status)

enh_anals = list.files(glue("{result_folder}"), pattern = "*enh_anal.tsv")
collect = lapply(enh_anals, function(x)
  fread(glue("{result_folder}{x}")))
collect = bind_rows(collect)

## find G4 peaks near to LTR sequences
ltrs = fread(ltrs)
ltrs = GRanges(
  seqnames = ltrs$V1,
  ranges = IRanges(
    start = ltrs$V2,
    end = ltrs$V3,
    names = rep("ltrs", nrow(ltrs))
  )
)

sign_g4s = collect %>% dplyr::select(seqnames, start, end)
sign_g4s = GRanges(
  seqnames = sign_g4s$seqnames,
  ranges = IRanges(
    start = sign_g4s$start,
    end = sign_g4s$end,
    names = rep("G4", nrow(sign_g4s))
  )
)
# expand G4 peaks by 1-1 kbp
sign_g4s = extendGR(sign_g4s, upstream = 500, downstream = 500)

# overlap with LTRs
ltr_ol = suppressWarnings(findOverlaps(sign_g4s, ltrs, type = "any", minoverlap = 1))

collect = collect %>% mutate(LTR_vicinity = "0")
for (hit in ltr_ol@from) {
  collect[hit, which(colnames(collect) == "LTR_vicinity")] = "1"
}

write_tsv(collect,
          glue("{result_folder}enhancer_analysis_output.tsv"))





