
# Copy Number Profile Generating

library(tidyverse)
library(readxl)

### PCAWG dataset
setwd("~/HRD/HRDCNA/data/download/cn_pcawg/")
txt <- list.files()
vcfs <- purrr::map(txt, ~data.table::fread(., select = c(1,2,3,4)))
names(vcfs) <- sub(".consensus.20170119.somatic.cna.annotated.txt", "", txt)
vcfs <- data.table::rbindlist(vcfs, use.names = FALSE, idcol = "sample")
names(vcfs) <- c("sample", "chromosome", "start", "end", "segVal")

# saveRDS(vcfs, file = "~/HRD/HRDCNA/data/cndata/cn_pcawg_wgs.rds")

### 560 breast dataset
setwd("~/HRD/HRDCNA/data/download/cn_breast560/")
vcf560 <- list.files()
vcf560s <- purrr::map(vcf560, ~data.table::fread(., select = c(2, 3,4,7)))
names(vcf560s) <- sapply(str_split(string = vcf560, pattern = "[.]"), "[[",1)
vcf560s <- data.table::rbindlist(vcf560s, use.names = FALSE, idcol = "sample")
names(vcf560s) <- c("sample", "chromosome", "start", "end", "segVal")

# saveRDS(vcf560s, file = "~/HRD/HRDCNA/data/cndata/cn_560_snp.rds")

### panel dataset
setwd("~/HRD/HRDCNA/")
cnv_panel501 <- read_excel("./data/download/cn_panel/CNV_12885_2022_9602_MOESM11_ESM.xlsx")
colnames(cnv_panel501) <- cnv_panel501[1, ]
cnv_data <- cnv_panel501[-1, -6]
names(cnv_data) <- c("sample","chromosome","start","end","segVal")
cnv_data$segVal <- as.numeric(cnv_data$segVal)
cnv_panel_85 <- cnv_data[(1:8723), ]
cnv_panel_416 <- cnv_data[-(1:8723), ]

# saveRDS(cnv_data, file = "./data/cndata/cn_panel.rds")
# saveRDS(cnv_panel_85, file = "./data/cndata/cn_panel_85.rds")
# saveRDS(cnv_panel_416, file = "./data/cndata/cn_panel_416.rds")

### 66 breast dataset
setwd("~/HRD/HRDCNA/")
arraycopy <- read.csv("./data/download/cn_breast66/segmentation_array.tsv",sep = "\t")
arraycopysig <- arraycopy[,c(1,2,3,4,6)]
names(arraycopysig) <- c("sample","chromosome","start","end","segVal")

wgscopy <- read.csv("./data/download/cn_breast66/segmentation_original_wgs.tsv",sep = "\t")
wgscopysig <- wgscopy[,c(1,2,3,4,6)]
names(wgscopysig) <- c("sample","chromosome","start","end","segVal")

copy10x<- read.csv("./data/download/cn_breast66/segmentations_10_coverage.tsv",sep = "\t")
copysig10x <- copy10x[,c(1,2,3,4,6)]
names(copysig10x) <- c("sample","chromosome","start","end","segVal")

copy15x<- read.csv("./data/download/cn_breast66/segmentations_15_coverage.tsv",sep = "\t")
copysig15x <- copy15x[,c(1,2,3,4,6)]
names(copysig15x) <- c("sample","chromosome","start","end","segVal")

copy30x <- read.csv("./data/download/cn_breast66/segmentations_30_coverage.tsv",sep = "\t")
copysig30x <- copy30x[,c(1,2,3,4,6)]
names(copysig30x) <- c("sample","chromosome","start","end","segVal")

# saveRDS(arraycopysig, file = "./data/cndata/cn_60_snp.rds")
# saveRDS(wgscopysig, file = "./data/cndata/cn_60_wgs.rds")
# saveRDS(copysig30x, file = "./data/cndata/cn_60_wgs30x.rds")
# saveRDS(copysig15x, file = "./data/cndata/cn_60_wgs15x.rds")
# saveRDS(copysig10x, file = "./data/cndata/cn_60_wgs10x.rds")

### TCGA dataset
setwd("~/HRD/HRDCNA/")
sample <- read.csv("./data/download/cn_tcga/tcga_sample", header = F)
dataall <- data.frame()
for( i in 1:nrow(sample)){
  cancer <- sample[i,]
  rdspath <-paste("./data/download/cn_tcga/cn_tcga/",cancer,".rds", sep = "")
  rds <- readRDS(rdspath)
  dataall <- rbind(dataall, rds)
}

# saveRDS(dataall, file = "./data/cndata/cn_tcga_snp.rds")




