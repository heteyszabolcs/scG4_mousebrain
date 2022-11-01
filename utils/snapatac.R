library("SnapATAC")
library("viridisLite")
library("ggplot2")

x.sp = createSnap(
  file="../data/CellRanger/unsorted/fragments.snap",
  sample="unsorted",
  num.cores=1
)

barcodes = read.csv(
  "../data/CellRanger/unsorted/singlecell.csv",
  head=TRUE
)

x.sp@barcode = barcodes$barcode
x.sp@bmat

x.sp@mmat

barcodes = barcodes[2:nrow(barcodes),]
promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1)
UMI = log(barcodes$passed_filters+1, 10)
data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio)
barcodes$promoter_ratio = promoter_ratio
p1 = ggplot(
  data, 
  aes(x= UMI, y= promoter_ratio)) + 
  geom_point(size=0.1, col="red") +
  theme_classic() +
  ggtitle("") +
  ylim(0, 1) + xlim(0, 6) +
  labs(x = "log10(UMI)", y="promoter ratio") 
p1 

barcodes.sel = barcodes[which(UMI >= 0.5 & UMI <= 5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6),]
rownames(barcodes.sel) = barcodes.sel$barcode
x.sp = x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),]
x.sp@metaData = barcodes.sel[x.sp@barcode,]
x.sp

x.sp = addBmatToSnap(x.sp, bin.size = 5000, num.cores = 1)


file.name = system.file("extdata", "demo.snap", package = "SnapATAC")
demo.sp = createSnap(file.name, sample="demo", do.par=FALSE)
demo.sp@feature
showBinSizes("../data/CellRanger/unsorted/fragments.snap")

