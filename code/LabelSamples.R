
# Sample labeling

library(tidyverse)
library(readxl)

setwd("~/HRD/HRDCNA")

### PCAWG dataset
pca <- read_xlsx("./data/download/type_pcawg/media-1.xlsx", sheet = 1)
pcawg <- pca[pca$group == "PCAWG",]
pcawgcancertype <- pcawg[,c(2,14)]

pcawghrd <- pcawg %>% filter(response != "none" & used_for_perf_eval == "TRUE")
pcawgbrca1 <- pcawghrd %>% filter(response == "BRCA1" & hr_status == "HR_deficient")
pcawgbrca2 <- pcawghrd %>% filter(response == "BRCA2" & hr_status == "HR_deficient")
pcawgcontro <- pcawg %>% filter(response == "none" & used_for_perf_eval == "TRUE" & hr_status == "HR_proficient")

pcawg_hrd <- rbind(pcawgbrca1, pcawgbrca2)
pcawg_hrr <- pcawgcontro

# saveRDS(pcawg_hrd, file = "./data/typedata/pcawg_hrd.rds")
# saveRDS(pcawg_hrr, file = "./data/typedata/pcawg_hrr.rds")

### 560 breast dataset
supplement1 <- read_excel("./data/download/type_breast560/supplement1.xlsx", sheet = 1, na = "NA")
supplement1 <- supplement1[-1,]
supplement1 <- as.data.frame(supplement1)
colnames(supplement1) <- supplement1[1,]
supplement1 <- supplement1[-1,]
evulate <- supplement1 %>% filter(supplement1$isUsedForEvaluation == "TRUE")

isKnownGermline <- evulate %>% filter(evulate$isKnownGermline == "TRUE")
newgermline <- evulate %>% filter(evulate$isNewGermline == "TRUE")
germlinemu <- rbind(isKnownGermline,newgermline)
brca1ger560 <- germlinemu %>% filter(Gene == "BRCA1")
brca2ger560 <- germlinemu %>% filter(Gene == "BRCA2")

stomameth <- evulate %>% filter(evulate$IsSomaticMeth == "TRUE")
brca1sto560 <- stomameth %>% filter(Gene == "BRCA1")
brca2sto560 <- stomameth %>% filter(Gene == "BRCA2")

isBrcaMonoallelic <- supplement1 %>% filter(supplement1$isBrcaMonoallelic == "TRUE")

control560 <- evulate %>% filter(evulate$isQuiescentGenomeControl == "TRUE")
all77 <- rbind(isKnownGermline,newgermline,stomameth)

# saveRDS(all77, file = "./data/typedata/a560_hrd.rds")
# saveRDS(control560, file = "./data/typedata/a560_hrr.rds")

### panel dataset
panel_85 <- read_excel("./data/download/type_panel/panel_85_12885_2022_9602_MOESM8_ESM.xlsx")
colnames(panel_85) <- panel_85[1, ]
panel_85 <- panel_85[-1, ]

panel_416 <- read_excel("./data/download/type_panel/panel_416_12885_2022_9602_MOESM10_ESM.xlsx")
colnames(panel_416) <- panel_416[1, ]
panel_416 <- panel_416[-1, ]

HRD_panel_85 <- panel_85[ , c(1,9,10,11,12)]
HRD_panel_416 <- panel_416[ , c(1,3,4,5,6)]
HRD_panel <- rbind(HRD_panel_85, HRD_panel_416)
HRD_panel$`HRD score` <- as.numeric(HRD_panel$`HRD score`)

panel_hrd <- HRD_panel %>% filter(HRD_panel$`HRD score` >= 38)
panel_hrr <- HRD_panel %>% filter(HRD_panel$`HRD score` < 38)

# saveRDS(panel_hrd, file = "./data/typedata/panel_all_hrd.rds")
# saveRDS(panel_hrr, file = "./data/typedata/panel_all_hrr.rds")


### 66 breast dataset
weahrd <- read_excel("./data/download/type_breast66/SupplementaryTableS2.xlsx")
colnames(weahrd) <- weahrd[2,]
weahrd <- weahrd[-c(1,2),]
weahrdtype <- weahrd[,c(1,32)]
weahrdtype <- as.data.frame(weahrdtype)

weahrdbrca1ger <- weahrdtype %>% filter(`Germline status` == "BRCA1" | `Germline status`=="BRCA1 + BRCA2")
weahrdbrca2ger <- weahrdtype %>% filter(`Germline status` == "BRCA2")
weahrdstome <- weahrdtype %>% filter(`Germline status` == "Somatic BRCA1 LOH + methylation" | `Germline status` == "Somatic BRCA1 LOH")
weahrdcontrol <- weahrdtype %>% filter(`Germline status` == "non-BRCA1/2" | `Germline status` == "BRCA2 (not biallelic)")
weahrdplab <- weahrdtype %>% filter(`Germline status` == "PALB2")
weauncladd <- weahrdtype %>% filter(`Germline status` == "non-BRCA1/2 (BRCA2 UV)")

wgs67_hrd_PALB2 <- rbind(weauncladd, weahrdstome, weahrdbrca1ger, weahrdbrca2ger)
wgs67_hrr_PALB2 <- rbind(weahrdcontrol)
wgs67_PALB2 <- weahrdplab

# saveRDS(wgs67_hrd_PALB2, file = "./data/typedata/wgs67_hrd.rds")
# saveRDS(wgs67_hrr_PALB2, file = "./data/typedata/wgs67_hrr.rds")



