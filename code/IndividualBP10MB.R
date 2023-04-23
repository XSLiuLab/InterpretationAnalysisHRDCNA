

### individual cancers bp10mb


vali_brca_sample <- readRDS("~/HRD/HRDCNA/supp/data/vali_brca_sample.rds")
vali_ov_sample <- readRDS("~/HRD/HRDCNA/supp/data/vali_ov_sample.rds")
vali_paca_sample <- readRDS("~/HRD/HRDCNA/supp/data/vali_paca_sample.rds")
vali_prad_sample <- readRDS("~/HRD/HRDCNA/supp/data/vali_prad_sample.rds")
vali_lymp_sample <- readRDS("~/HRD/HRDCNA/supp/data/vali_lymp_sample.rds")
vali_live_sample <- readRDS("~/HRD/HRDCNA/supp/data/vali_live_sample.rds")


tally_W_pcawg <- readRDS("~/HRD/HRDCNA/data/tallydata/tally_W_pcawg.rds")
nmf_pcawg <- as.data.frame(tally_W_pcawg$nmf_matrix)
hrd_sampleid <- readRDS("~/HRD/HRDCNA/supp/data/hrd_sampleid.rds")
hrp_sampleid <- readRDS("~/HRD/HRDCNA/supp/data/hrp_sampleid.rds")


brca <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_brca_sample$Sample)

brca$type <- ifelse(rownames(brca) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(brca) %in% hrp_sampleid$SampleID, "0", "null"))
brca <- brca %>% filter(type!="null")


ov <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_ov_sample$Sample)

ov$type <- ifelse(rownames(ov) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(ov) %in% hrp_sampleid$SampleID, "0", "null"))
ov <- ov %>% filter(type!="null")


paca <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_paca_sample$Sample)

paca$type <- ifelse(rownames(paca) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(paca) %in% hrp_sampleid$SampleID, "0", "null"))
paca <- paca %>% filter(type!="null")


prad <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_prad_sample$Sample)

prad$type <- ifelse(rownames(prad) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(prad) %in% hrp_sampleid$SampleID, "0", "null"))
prad <- prad %>% filter(type!="null")



lymp <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_lymp_sample$Sample)

lymp$type <- ifelse(rownames(lymp) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(lymp) %in% hrp_sampleid$SampleID, "0", "null"))
lymp <- lymp %>% filter(type!="null")


live <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_live_sample$Sample)

live$type <- ifelse(rownames(live) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(live) %in% hrp_sampleid$SampleID, "0", "null"))
live <- live %>% filter(type!="null")




plot.roc(brca$type, brca$`BP10MB[1]`, col="#bc5148",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Beast Cancer : AUC = %.0f%%",
         print.auc.y=50)

plot.roc(ov$type, ov$`BP10MB[1]`, col="#4a4266",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Ovarian Cancer : AUC = %.0f%%",
         print.auc.y=45, add=T)

plot.roc(paca$type, paca$`BP10MB[1]`, col="#549688",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Pancreatic Cancer : AUC = %.0f%%",
         print.auc.y=40, add=T)

plot.roc(prad$type, prad$`BP10MB[1]`, col="#d0a727",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Prostate Cancer : AUC = %.0f%%",
         print.auc.y=35, add=T)

plot.roc(lymp$type, lymp$`BP10MB[1]`, col="#19419B",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Lymphatic Cancer : AUC = %.0f%%",
         print.auc.y=30, add=T)

plot.roc(live$type, live$`BP10MB[1]`, col="#C15F16",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Liver Cancer : AUC = %.0f%%",
         print.auc.y=25, add=T)





