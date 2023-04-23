


# HRDCNA vs HRD score, HRDetect (new prad cohort)

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")

tally_W_tcga <- readRDS("~/HRD/HRDCNA/data/tallydata/tally_W_tcga.rds")

nmf_tcga <- as.data.frame(tally_W_tcga$nmf_matrix)
churn.pred <- predict(churn.gbmtest4, nmf_tcga, n.trees = bestTree, type = "response")

nmf_tcga$HRDCNAScore <- churn.pred
nmf_tcga$Sample <- substr(rownames(nmf_tcga), 1, 12)
vali_all <- nmf_tcga[ , c(82,81)]
rownames(vali_all) <- NULL


vali_vs <- read_xlsx("./data/download/vs_tcga/hrdetect_brca.xlsx")
vali_all <- inner_join(vali_all,vali_vs, by = "Sample")


LLL <- read.table("./data/download/hrd_tcga/TCGA.HRD_withSampleID.txt")
LLL <- as.data.frame(t(LLL))
colnames(LLL) <- c("Sample", "TAI", "LST", "LOH", "HRD")
LLL <- LLL[-1, ]
rownames(LLL) <- NULL
LLL$Sample <- substr(LLL$Sample, 1, 12)



vali_all <- inner_join(vali_all, LLL, by = "Sample")

vali_all_hrd <- vali_all %>% filter(vali_all$`HR Status` == "HR-deficiency")
vali_all_hrp <- vali_all %>% filter(vali_all$`HR Status` == "HR-proficiency")

vali_all$Type <- ifelse(vali_all$Sample %in% vali_all_hrd$Sample, "1",
                        ifelse(vali_all$Sample %in% vali_all_hrp$Sample, "0", "null"))
vali_all <- vali_all %>% filter(vali_all$Type != "null")
vali_all <- unique(vali_all)



p <- ggplot(vali_all, aes(x=HRDCNAScore,y=HRDetectScore, colour=Type))+
  geom_point(size=3)+
  xlab("HRDCNA Score")+ylab("HRDetect Score")+
  theme(title=element_text(size=15),
        panel.background  = element_rect(fill="white",colour = "black",size = 1),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = c(.75,.25),
        legend.key=element_blank(),
        axis.text.x  =element_text(size=15,color = "black"),
        axis.text.y = element_text(size=15,color = "black"))+
  geom_hline(yintercept = 0.7)+
  geom_vline(xintercept = 0.2)+
  annotate("text",label="HRD",x=0.6,y=0.85,size=5,colour="#D74B4B")+
  annotate("text",label="HRP",x=0.3,y=0.45,size=5,colour="#354B5E")+
  scale_discrete_manual(values = c("#716e77", "#e79686"),
                        aesthetics = "colour",
                        labels=c("HRP","HRD"))
p



churn.roc = pROC::roc(vali_all$Type, vali_all$HRDCNAScore)
churn.roc$auc

churn.roc = pROC::roc(vali_all$Type, vali_all$HRDetectScore)
churn.roc$auc

plot.roc(vali_all$Type, vali_all$HRDCNAScore, col="#bc5148",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="HRDCNA : AUC = %.0f%%",
         print.auc.y=50)

plot.roc(vali_all$Type, vali_all$HRDetectScore, col="#549688",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="HRDetect : AUC = %.0f%%",
         print.auc.y=45, add=T)




### HRDCNA vs HRD score (held-out data)

testall <- readRDS("~/HRD/HRDCNA/data/modeldata/testall.rds")

HRD_560 <- readRDS("~/HRD/yhz_CNHRD/data_new/HRDdata/HRD_560.rds")
HRD_pcawg <- readRDS("~/HRD/yhz_CNHRD/data_new/HRDdata/HRD_pcawg.rds")

HRD_560$HRDALL <- HRD_560$TAI + HRD_560$LST + HRD_560$LOH
HRD_560$sample <- paste(HRD_560$SampleID, "a", sep = "")
HRD_560 <- HRD_560[ , c(8,2,4,3,7)]


HRD_all <- rbind(HRD_pcawg, HRD_560)


churn.pred <- predict(churn.gbmtest4, testall, n.trees = bestTree, type = "response")
testall$churn.pred <- churn.pred
testall$sample <- rownames(testall)
testall <- testall[ , c(83,82,81)]
rownames(testall) <- NULL

data <- inner_join(HRD_all, testall, by="sample")

churn.roc = pROC::roc(data$type, data$churn.pred)
churn.roc$auc
# Area under the curve: 0.9727
churn.roc = pROC::roc(data$type, data$TAI)
churn.roc$auc
# Area under the curve: 0.9413
churn.roc = pROC::roc(data$type, data$LST)
churn.roc$auc
# Area under the curve: 0.9728
churn.roc = pROC::roc(data$type, data$LOH)
churn.roc$auc
# Area under the curve: 0.9380
churn.roc = pROC::roc(data$type, data$HRDALL)
churn.roc$auc
# Area under the curve: 0.9625


plot.roc(data$type, data$churn.pred, col="#bc5148",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="HRDCNA : AUC = %.0f%%",
         print.auc.y=50)

plot.roc(data$type, data$TAI, col="#549688",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="TAI : AUC = %.0f%%",
         print.auc.y=45, add=T)

plot.roc(data$type, data$LST, col="#4a4266",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="LST : AUC = %.0f%%",
         print.auc.y=40, add=T)

plot.roc(data$type, data$LOH, col="#d0a727",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="LOH : AUC = %.0f%%",
         print.auc.y=35, add=T)

plot.roc(data$type, data$HRDALL, col="#5B6049",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="HRD Score : AUC = %.0f%%",
         print.auc.y=30, add=T)



