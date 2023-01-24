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

# export folder
result_folder = "../results/GenomicRanges/unsorted_outputs/"

# data
g4_peaks = "../results/Seurat/callpeaks_unsorted/peak_sets/peaks_per_clusters.bed"
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
  dplyr::filter(V6 %in% c("0", "1", "2", "3", "4", "5", "6")) %>%
  ggplot(., aes(x = diff, fill = V6)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  geom_density(alpha = 1.0) +
  xlim(200, 750) +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "G4 length distributions / Seurat cluster",
       x = "length (bp)",
       y = "Density",
       fill = "Seurat cluster") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
hist

ggsave(
  glue("{result_folder}G4_length_distr_Seurat_clst.pdf"),
  plot = hist,
  width = 10,
  height = 7,
  device = "pdf"
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
g4s_no_ol = g4_peaks[c("0", "1", "2", "3", "4", "5", "6")]
g4_cm_enh_quant = numeric()
for (cluster in names(g4s_no_ol)) {
  ol = suppressWarnings(findOverlaps(g4_peaks[[cluster]], cm_enh, type = "any", minoverlap = 1))
  g4_cm_enh_quant = c(g4_cm_enh_quant, length(ol))
}




get_signals = function(selected_cluster = "0") {
  require("wigglescout")
  cl_cm_ol = suppressWarnings(subsetByOverlaps(g4_peaks[[selected_cluster]], cm_enh, type = "any", minoverlap = 1))
  cl_cm_ol = as.tibble(cl_cm_ol)
  cl_cm_ol = suppressWarnings(subsetByOverlaps(g4_peaks[[selected_cluster]], cm_enh, type = "any", minoverlap = 1))
  cl_cm_ol = as.tibble(cl_cm_ol)
  
  bigwigs = c(
    "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/0.bw",
    "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/1.bw",
    "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/2.bw",
    "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/3.bw",
    "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/4.bw",
    "../results/Seurat/callpeaks_unsorted/cluster_bigwigs/5.bw",
    "../data/bw/Martire2019_ESC_H33WT_H3K27ac.bw"
  )
  
  bed = cl_cm_ol[, 1:3]
  write_tsv(
    bed,
    glue(
      "../data/bed/Cruz_Molina_enh-unsorted_cluster{selected_cluster}.bed"
    ),
    col_names = FALSE
  )
  bed = glue("../data/bed/Cruz_Molina_enh-unsorted_cluster{selected_cluster}.bed")
  selected_bw = glue("../results/Seurat/callpeaks_unsorted/cluster_bigwigs/{selected_cluster}.bw")
  bigwigs = c(bigwigs[which(bigwigs != selected_bw)], selected_bw)
  
  read_cov = bw_loci(bigwigs, loci = bed)
  read_cov = as.data.frame(read_cov)
  
  columns = colnames(read_cov)[grepl("X", colnames(read_cov)[which(colnames(read_cov) != glue("X{selected_cluster}"))])]
  read_cov_filt = read_cov %>% 
    mutate(average = rowMeans(dplyr::select(., columns))) %>% 
    mutate(enrichment = get(glue("X{selected_cluster}")) / average) %>% 
    dplyr::filter(enrichment >= 10)
  
  return(read_cov_filt)
  
}

cl3_over_cm_enh = get_signals("3")
cl5_over_cm_enh = get_signals("5")

g4_gl_enh_quant = numeric()
for (cluster in names(g4s_no_ol)) {
  ol = suppressWarnings(findOverlaps(g4_peaks[[cluster]], gl_enh, type = "any", minoverlap = 1))
  g4_gl_enh_quant = c(g4_gl_enh_quant, length(ol))
}

clusters = c("0", "1", "2", "3", "4", "5", "6")
cm_bar = tibble(
  overlap = g4_cm_enh_quant,
  peak_count = sapply(clusters, function(x) length(g4_peaks[[which(names(g4_peaks) == x)]])),
  enhancer_ratio = (overlap / peak_count) * 100,
  cluster = names(g4s_no_ol),
  fill_col = names(g4s_no_ol)
) %>%
  ggplot(data = ., aes(
    x = reorder(cluster,-enhancer_ratio),
    y = enhancer_ratio,
    fill = fill_col
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  labs(
    title = expression(paste(
      "G4 overlaps with active enhancers of ",
      italic("Cruz-Molina et al.")
    )),
    x = "Seurat cluster",
    y = "active enhancer %",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
cm_bar

cm_bar_input = tibble(
  overlap = g4_cm_enh_quant,
  peak_count = sapply(clusters, function(x) length(g4_peaks[[which(names(g4_peaks) == x)]])),
  enhancer_ratio = (overlap / peak_count) * 100,
  cluster = names(g4s_no_ol),
  fill_col = names(g4s_no_ol)
)

print(glue("Average enhancer percentage across clusters: {round(mean(cm_bar_input$enhancer_ratio), 2)} % (SD: {round(sd(cm_bar_input$enhancer_ratio), 2)})"))
print(glue("Average peak number over active enhancer: {round(mean(cm_bar_input$overlap), 0)} (SD: {round(sd(cm_bar_input$overlap), 0)})"))

ggsave(
  glue("{result_folder}G4_overlaps_w_CruzM_active_enh.pdf"),
  plot = cm_bar,
  width = 10,
  height = 7,
  device = "pdf"
)

clusters = c("0", "1", "2", "3", "4", "5", "6")
gl_bar = tibble(
  overlap = g4_gl_enh_quant,
  peak_count = sapply(clusters, function(x) length(g4_peaks[[which(names(g4_peaks) == x)]])),
  enhancer_ratio = (overlap / peak_count) * 100,
  cluster = names(g4s_no_ol),
  fill_col = names(g4s_no_ol)
) %>%
  ggplot(data = ., aes(
    x = reorder(cluster,-enhancer_ratio),
    y = enhancer_ratio,
    fill = fill_col
  )) +
  geom_bar(stat = "identity",
           width = 0.5,
           color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title = expression(paste(
      "G4 overlaps with active enhancers of ",
      italic("Glaser et al.")
    )),
    x = "Seurat cluster",
    y = "active enhancer %",
    fill = "Seurat cluster"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  )
gl_bar

ggsave(
  glue("{result_folder}G4_overlaps_w_Glaser_active_enh.pdf"),
  plot = gl_bar,
  width = 10,
  height = 7,
  device = "pdf"
)

## HOMER annotation of G4 peaks
result_folder = "../results/Seurat/callpeaks_unsorted/peak_sets/"
bw_folder = "../data/GSE157637/"
reference = "mm10"

annotation = function(narrowpeak = "0_peaks.narrowPeak",
                      percentile = 0.75) {
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
    dplyr::select(-.,-..) %>%
    mutate(start = as.character(start), end = as.character(end))
  
  # create HOMER input
  # keep peak scores with above 75th percentile
  np = np %>% dplyr::filter(signalValue >= quantile(signalValue, percentile))
  write_tsv(np, glue("{result_folder}{cluster}_robust_peaks.bed"))
  
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

# check package plyranges
if("plyranges" %in% installed.packages()) {
  print(" plyranges is OK")
} else {
  install.packages("plyranges")
  library("plyranges")
}

## complete annotated peaks with enhancer signatures
add_enhancer_status = function(annot_peak_file = "0_peaks_annot.tsv") {
  print(annot_peak_file)
  
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
  bws_ac = list.files(bw_folder, pattern = "*.bw", full.names = TRUE)
  bws_ac = bws_ac[grep("H3K27ac", bws_ac)]
  labels_ac = unname(sapply(bws_ac, function(x)
    str_split(x, "//")[[1]][2]))
  
  bws_k4 = list.files(bw_folder, pattern = "*.bw", full.names = TRUE)
  bws_k4 = bws_k4[grep("H3K4me3", bws_k4)]
  labels_k4 = unname(sapply(bws_k4, function(x)
    str_split(x, "//")[[1]][2]))
  
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
  
  aggr_ac = suppressWarnings(bw_loci(bws_ac, annot_bed, labels = labels_ac))
  # left join by plyranges
  output = join_overlap_left_directed(annot_bed, aggr_ac)
  
  bws_k4 = bws_k4[!str_detect(bws_k4, "H3K4me3_Astrocytes")]
  labels_k4 = labels_k4[!str_detect(labels_k4, "H3K4me3_Astrocytes")]
  aggr_k4 = 
    suppressWarnings(bw_loci(bws_k4, annot_bed, labels = labels_k4))
  output = join_overlap_left_directed(output, aggr_k4)
  output = as_tibble(output)

# export
write_tsv(output, glue("{result_folder}{cluster}_enh_anal.tsv"))

return(output)

}

# run on all annotated peaks

print(glue("####### {result_folder}"))
annot_peaks = list.files(glue("{result_folder}"), pattern = "*peaks_annot.tsv")
print(length(annot_peaks))
lapply(annot_peaks, add_enhancer_status)

enh_anals = list.files(glue("{result_folder}"), pattern = "*enh_anal.tsv")
collect = lapply(enh_anals, function(x)
  fread(glue("{result_folder}{x}")))
collect = bind_rows(collect)

## find G4 peaks near to LTR sequences
ltrs = fread(ltrs)
ltrs = GRanges(seqnames = ltrs$V1,
               ranges = IRanges(
                 start = ltrs$V2,
                 end = ltrs$V3,
                 names = rep("ltrs", nrow(ltrs))
               ))

sign_g4s = collect %>% dplyr::select(seqnames, start, end)
sign_g4s = GRanges(
  seqnames = sign_g4s$seqnames,
  ranges = IRanges(
    start = sign_g4s$start,
    end = sign_g4s$end,
    names = rep("G4", nrow(sign_g4s))
  )
)
# expand G4 peaks
sign_g4s = extendGR(sign_g4s, upstream = 500, downstream = 500)

# overlap with LTRs
ltr_ol = suppressWarnings(findOverlaps(sign_g4s, ltrs, type = "any", minoverlap = 1))

collect = collect %>% mutate(LTR_vicinity = "0")
for (hit in ltr_ol@from) {
  collect[hit, which(colnames(collect) == "LTR_vicinity")] = "1"
}

write_tsv(collect,
          glue("{result_folder}enhancer_analysis_output.tsv"))
