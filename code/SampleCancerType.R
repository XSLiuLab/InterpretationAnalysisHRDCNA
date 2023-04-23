
library(tidyverse)

### cancer type frequencie

pca <- read_xlsx("~/HRD/HRDCNA0/data/download/type_pcawg/media-1.xlsx", sheet = 1)

pcawg <- pca[pca$group == "PCAWG",]
pcawgcancertype <- pcawg[,c(2,14)]

pcawghrd <- pcawg %>% filter(response != "none" & used_for_perf_eval == "TRUE")
pcawgbrca1 <- pcawghrd %>% filter(response == "BRCA1" & hr_status == "HR_deficient")
pcawgbrca2 <- pcawghrd %>% filter(response == "BRCA2" & hr_status == "HR_deficient")
pcawgcontro <- pcawg %>% filter(response == "none" & used_for_perf_eval == "TRUE" & hr_status == "HR_proficient")

pcawg_hrd <- rbind(pcawgbrca1, pcawgbrca2)
pcawg_hrr <- pcawgcontro

pcawg_all <- rbind(pcawg_hrd, pcawg_hrr)

pcawg_all_cancertype <- pcawg_all[,c(2,14)]
colnames(pcawg_all_cancertype) <- c("Sample", "CancerType")



supplement1 <- read_excel("~/HRD/HRDCNA0/data/download/type_breast560/supplement1.xlsx", sheet = 1, na = "NA")
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

breast560 <- rbind(control560, all77)

breast560 <- breast560[1]
breast560$CancerType <- "Breast"


allsample <- rbind(pcawg_all_cancertype, breast560)


fre <- as.data.frame(table(allsample$CancerType))
colnames(fre) <- c("CancerType", "SampleNumber")
fre$Frequency <- round(fre$SampleNumber/sum(fre$SampleNumber), 3)


library(writexl)

write_xlsx(allsample, "./supp/data/SampleCancerType.xlsx")
write_xlsx(fre, "./supp/data/CancerTypeFrequency.xlsx")




