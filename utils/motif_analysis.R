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

# result / peak folder
peak_folder = "../results/Seurat/callpeaks_GFPsorted/peak_sets/"
result_folder = "../results/Seurat/callpeaks_GFPsorted/"

# Signac MACS2 peak set
signac_macs2 = fread(glue("{peak_folder}enhancer_analysis_output.tsv"))

# pull genes close to G4 structure and show elevated H3K4me3 signal (like above median H3K4me3 signal)
tss_proximal = signac_macs2 %>% filter(abs(Distance_to_TSS) < 3000)

k4me3_status = tss_proximal %>% select(Gene_name, starts_with("H3K4me3")) %>% 
  mutate(H3K4me3_mean = rowMeans(select(., starts_with("H3K4me3")), na.rm = TRUE)) %>% 
  arrange(desc(H3K4me3_mean)) %>% 
  filter(H3K4me3_mean >= median(H3K4me3_mean)) %>% select(Gene_name)

write_tsv(k4me3_status, glue("{result_folder}motif_analysis_input.txt"), col_names = FALSE)

# run findMotifs.pl of HOMER
system(
  glue("mkdir -p {result_folder}findMotif_output/"))

system(
  glue(
    "findMotifs.pl {result_folder}motif_analysis_input.txt mouse {result_folder}findMotif_output/ -start -500 -end 100 -len 8,10,12 -p 4"
))
