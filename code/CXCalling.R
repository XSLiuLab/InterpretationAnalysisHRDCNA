
# CX calling

library(tidyverse)
library(CINSignatureQuantification)

cn_pcawg_wgs <- readRDS("./data/cndata/cn_pcawg_wgs.rds")
cn_560_snp <- readRDS("./data/cndata/cn_560_snp.rds")

### wgs_pcawg
dat = data.frame(cn_pcawg_wgs)
dat[,c(2)][dat[,c(2)]==23] <- "X"
dat[,c(2)][dat[,c(2)]==24] <- "Y"
cn <- dat
cn <- cn[ , c(2,3,4,5,1)]

cn <- na.omit(cn)
mySigs = quantifyCNSignatures(cn)

pcawg_dt <- mySigs@activities[["rawAct0"]]

### snp_560
dat = data.frame(cn_560_snp)
dat[,c(2)][dat[,c(2)]==23] <- "X"
dat[,c(2)][dat[,c(2)]==24] <- "Y"
cn <- dat
cn <- cn[ , c(2,3,4,5,1)]

cn <- na.omit(cn)
mySigs = quantifyCNSignatures(cn)

snp560_dt <- mySigs@activities[["rawAct0"]]

# saveRDS(pcawg_dt, file = "./data/modelselection/data/pcawg_cx.rds")
# saveRDS(snp560_dt, file = "./data/modelselection/data/snp560_cx.rds")


