
# Calling CNA Features

library(sigminer)

### data preparing
cn_pcawg_wgs <- readRDS("~/HRD/HRDCNA/data/cndata/cn_pcawg_wgs.rds")
cn_560_snp <- readRDS("~/HRD/HRDCNA/data/cndata/cn_560_snp.rds")
cn_60_wgs <- readRDS("~/HRD/HRDCNA/data/cndata/cn_60_wgs.rds")
cn_60_snp <- readRDS("~/HRD/HRDCNA/data/cndata/cn_60_snp.rds")
cn_60_wgs30x <- readRDS("~/HRD/HRDCNA/data/cndata/cn_60_wgs30x.rds")
cn_60_wgs15x <- readRDS("~/HRD/HRDCNA/data/cndata/cn_60_wgs15x.rds")
cn_60_wgs10x <- readRDS("~/HRD/HRDCNA/data/cndata/cn_60_wgs10x.rds")
cn_tcga_snp <- readRDS("~/HRD/HRDCNA/data/cndata/cn_tcga_snp.rds")

cn_panel <- readRDS("~/HRD/HRDCNA/data/cndata/cn_panel.rds")
cn_panel_85 <- readRDS("~/HRD/HRDCNA/data/cndata/cn_panel_85.rds")
cn_panel_416 <- readRDS("~/HRD/HRDCNA/data/cndata/cn_panel_416.rds")


### PCAWG dataset
callpcawg <- read_copynumber(vcfs,
                             seg_cols = c("chromosome", "start", "end", "segVal"),
                             genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_pcawg <- sig_tally(callpcawg, method = "W")

# saveRDS(tally_W_pcawg, file = "./data/tallydata/tally_W_pcawg.rds")

### 560 breast dataset
call560 <- read_copynumber(vcf560s,
                           seg_cols = c("chromosome", "start", "end", "segVal"),
                           genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_560 <- sig_tally(call560, method = "W")

# saveRDS(tally_W_560, file = "./data/tallydata/snp560/tally_W_560.rds")

### panel dataset
callpanel <- read_copynumber(cnv_data,
                             seg_cols = c("chromosome", "start", "end", "segVal"),
                             genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_panel <- sig_tally(callpanel, method = "W")

callpanel85 <- read_copynumber(cnv_panel_85,
                               seg_cols = c("chromosome", "start", "end", "segVal"),
                               genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_panel85 <- sig_tally(callpanel85, method = "W")

callpanel416 <- read_copynumber(cnv_panel_416,
                                seg_cols = c("chromosome", "start", "end", "segVal"),
                                genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_panel416 <- sig_tally(callpanel416, method = "W")

# saveRDS(tally_W_panel, file = "./data/tallydata/tally_W_panel.rds")
# saveRDS(tally_W_panel85, file = "./data/tallydata/tally_W_panel85.rds")
# saveRDS(tally_W_panel416, file = "./data/tallydata/tally_W_panel416.rds")

### 66 breast dataset
call60wgs <- read_copynumber(copysigwgs,
                             seg_cols = c("chromosome", "start", "end", "segVal"),
                             genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_60wgs <- sig_tally(call60wgs, method = "W")

call60snp <- read_copynumber(copysigsnp,
                             seg_cols = c("chromosome", "start", "end", "segVal"),
                             genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_60array <- sig_tally(call60snp, method = "W")

call60wgs30 <- read_copynumber(copysigwgs30,
                               seg_cols = c("chromosome", "start", "end", "segVal"),
                               genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_60wgs30x <- sig_tally(call60wgs30, method = "W")

call60wgs15 <- read_copynumber(copysigwgs15,
                               seg_cols = c("chromosome", "start", "end", "segVal"),
                               genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_60wgs15x <- sig_tally(call60wgs15, method = "W")

call60wgs10 <- read_copynumber(copysigwgs10,
                               seg_cols = c("chromosome", "start", "end", "segVal"),
                               genome_build = "hg19", complement = FALSE, verbose = TRUE)
tally_W_60wgs10x <- sig_tally(call60wgs10, method = "W")

# saveRDS(tally_W_60wgs, file = "./data/tallydata/tally_W_60wgs.rds")
# saveRDS(tally_W_60array, file = "./data/tallydata/tally_W_60array.rds")
# saveRDS(tally_W_60wgs30x, file = "./data/tallydata/tally_W_60wgs30x.rds")
# saveRDS(tally_W_60wgs15x, file = "./data/tallydata/tally_W_60wgs15x.rds")
# saveRDS(tally_W_60wgs10x, file = "./data/tallydata/tally_W_60wgs10x.rds")

### TCGA dataset
call <- read_copynumber(dataall,
                        seg_cols = c("chromosome", "start", "end", "segVal"),
                        genome_build = "hg38", complement = FALSE, verbose = TRUE)
tally_W_tcga <- sig_tally(cnall, method = "W")

# saveRDS(tally_W_tcga, file = "./data/tallydata/tally_W_tcga.rds")


