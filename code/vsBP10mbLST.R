
library(tidyverse)
library(pROC)

### BP10MB[1] & LST

tally_W_80 <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_80.rds")
HRD_80 <- readRDS("~/HRD/yhz_CNHRD/data_new/HRDdata/HRD_80.rds")

a80_hrd <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/a80_hrd.rds")
a80_hrr <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/a80_hrr.rds")


nmf80 <- as.data.frame(tally_W_80$nmf_matrix)
nmf80$type <- ifelse(rownames(nmf80) %in% a80_hrd$Sample, "1",
                     ifelse(rownames(nmf80) %in% a80_hrr$Sample, "0", "null"))
nmf80 <- nmf80 %>% filter(type!="null")

plot.roc(nmf80$type, nmf80$`BP10MB[1]`, col="#bc5148",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="AUC of BP10MB[1]: %.1f%%",
         print.auc.y=50)


HRD_80$type <- ifelse(rownames(HRD_80) %in% a80_hrd$Sample, "1",
                     ifelse(rownames(HRD_80) %in% a80_hrr$Sample, "0", "null"))
HRD_80 <- HRD_80 %>% filter(type!="null")

plot.roc(HRD_80$type, HRD_80$HRD, col="#4a4266",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="AUC of LOH: %.1f%%",
         print.auc.y=45, add=T)

plot.roc(HRD_80$type, HRD_80$TAI, col="#549688",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="AUC of TAI: %.1f%%",
         print.auc.y=40, add=T)

plot.roc(HRD_80$type, HRD_80$LST, col="#d0a727",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="AUC of LST: %.1f%%",
         print.auc.y=35, add=T)




HRD_alldata <- readRDS("~/HRD/yhz_CNHRD/data_new/HRDdata/HRD_alldata.rds")
testall <- readRDS("~/HRD/HRDCNA/data/modeldata/testall.rds")

testall$sample <- rownames(testall)
testall <- inner_join(testall, HRD_alldata, by = "sample")


plot.roc(testall$type, testall$`BP10MB[1]`, col="#bc5148",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="BP10MB[1] : AUC = %.0f%%",
         print.auc.y=50)

plot.roc(testall$type, testall$LOH, col="#4a4266",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="LOH : AUC = %.0f%%",
         print.auc.y=45, add=T)

plot.roc(testall$type, testall$TAI, col="#549688",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="TAI : AUC = %.0f%%",
         print.auc.y=40, add=T)

plot.roc(testall$type, testall$LST, col="#d0a727",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="LST : AUC =  %.0f%%",
         print.auc.y=35, add=T)


