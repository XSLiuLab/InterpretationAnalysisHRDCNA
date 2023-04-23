
library(tidyverse)
library(ggplot2)


### vali in individual cancer (breast, prostate, pancreas, ovarian – i.e. those with high frequency of HRD)

pca <- read_xlsx("~/HRD/HRDCNA0/data/download/type_pcawg/media-1.xlsx", sheet = 1)

pcawg <- pca%>% filter(pca$group == "PCAWG" & pca$used_for_perf_eval == "TRUE")

pcawghrd <- pcawg %>% filter(hr_status == "HR_deficient")
pcawghrp <- pcawg %>% filter(hr_status == "HR_proficient")

hrd_sampleid <- as.data.frame(pcawghrd$sample)
colnames(hrd_sampleid) <- "SampleID"
hrp_sampleid <- as.data.frame(pcawghrp$sample)
colnames(hrp_sampleid) <- "SampleID"


# saveRDS(hrd_sampleid, file = "./supp/data/hrd_sampleid.rds")
# saveRDS(hrp_sampleid, file = "./supp/data/hrp_sampleid.rds")



### individual cancer

pcawg_all_cancertype <- pcawg[,c(2,14)]
colnames(pcawg_all_cancertype) <- c("Sample", "CancerType")


vali_brca_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Breast")
vali_paca_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Pancreas")
vali_prad_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Prostate")
vali_ov_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Ovary")
vali_esop_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Esophagus")
vali_kidn_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Kidney")
vali_live_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Liver")
vali_lymp_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Lymphoid")
vali_medu_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Medulloblastoma")
vali_skin_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Skin")
vali_stom_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Stomach")


# saveRDS(vali_brca_sample, "./supp/data/vali_brca_sample.rds")
# saveRDS(vali_paca_sample, "./supp/data/vali_paca_sample.rds")
# saveRDS(vali_prad_sample, "./supp/data/vali_prad_sample.rds")
# saveRDS(vali_ov_sample, "./supp/data/vali_ov_sample.rds")
# saveRDS(vali_esop_sample, "./supp/data/vali_esop_sample.rds")
# saveRDS(vali_kidn_sample, "./supp/data/vali_kidn_sample.rds")
# saveRDS(vali_live_sample, "./supp/data/vali_live_sample.rds")
# saveRDS(vali_lymp_sample, "./supp/data/vali_lymp_sample.rds")
# saveRDS(vali_medu_sample, "./supp/data/vali_medu_sample.rds")
# saveRDS(vali_skin_sample, "./supp/data/vali_skin_sample.rds")
# saveRDS(vali_stom_sample, "./supp/data/vali_stom_sample.rds")




### vali

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")

tally_W_pcawg <- readRDS("~/HRD/HRDCNA/data/tallydata/tally_W_pcawg.rds")
nmf_pcawg <- as.data.frame(tally_W_pcawg$nmf_matrix)


#### BRCA

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_brca_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9953

#### ROC PRC
plot_data <- data.frame(sampletype = data$type, probablity = churn.pred,
                    sampleid = rownames(data))

### ROC add confidence interval
roc <- roc(plot_data$sampletype, plot_data$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100)
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255))


pre_obj1 <- mmdata(plot_data$probablity, plot_data$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_single <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "BRCA"

roc <- pre1_df[pre1_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(aes(colour="red"))+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14))+
  scale_color_manual(values='#bc5148')+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p1 <- p + annotate("text", x = .65, y = .35, size=5, colour = "black",
                   label = paste("AUC =", round(auc_single$aucs[1], 4)))

p1



#### PACA

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_paca_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9371

#### ROC PRC
plot_data <- data.frame(sampletype = data$type, probablity = churn.pred,
                        sampleid = rownames(data))

### ROC add confidence interval
roc <- roc(plot_data$sampletype, plot_data$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100)
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255))


pre_obj1 <- mmdata(plot_data$probablity, plot_data$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_single <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "PACA"

roc <- pre1_df[pre1_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(aes(colour="red"))+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14))+
  scale_color_manual(values='#bc5148')+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p1 <- p + annotate("text", x = .65, y = .35, size=5, colour = "black",
                   label = paste("AUC =", round(auc_single$aucs[1], 4)))

p1


#### PRAD

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_prad_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9555

#### ROC PRC
plot_data <- data.frame(sampletype = data$type, probablity = churn.pred,
                        sampleid = rownames(data))

### ROC add confidence interval
roc <- roc(plot_data$sampletype, plot_data$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100)
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255))



pre_obj1 <- mmdata(plot_data$probablity, plot_data$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_single <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "PRAD"

roc <- pre1_df[pre1_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(aes(colour="red"))+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14))+
  scale_color_manual(values='#bc5148')+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p1 <- p + annotate("text", x = .65, y = .35, size=5, colour = "black",
                   label = paste("AUC =", round(auc_single$aucs[1], 4)))

p1



#### OV

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_ov_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 1.0000

#### ROC PRC
plot_data <- data.frame(sampletype = data$type, probablity = churn.pred,
                        sampleid = rownames(data))

pre_obj1 <- mmdata(plot_data$probablity, plot_data$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_single <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "OV"

roc <- pre1_df[pre1_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(aes(colour="red"))+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14))+
  scale_color_manual(values='#bc5148')+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p1 <- p + annotate("text", x = .65, y = .35, size=5, colour = "black",
                   label = paste("AUC =", round(auc_single$aucs[1], 4)))

p1


### Lymphoid

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_lymp_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9769

#### ROC
plot_data <- data.frame(sampletype = data$type, probablity = churn.pred,
                        sampleid = rownames(data))

### ROC add confidence interval
roc <- roc(plot_data$sampletype, plot_data$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100)
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255))



pre_obj1 <- mmdata(plot_data$probablity, plot_data$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_single <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "Lymphoid"

roc <- pre1_df[pre1_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(aes(colour="red"))+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14))+
  scale_color_manual(values='#bc5148')+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p1 <- p + annotate("text", x = .65, y = .35, size=5, colour = "black",
                   label = paste("AUC =", round(auc_single$aucs[1], 4)))

p1


### Liver

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_live_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 1.0000

#### ROC
plot_data <- data.frame(sampletype = data$type, probablity = churn.pred,
                        sampleid = rownames(data))

pre_obj1 <- mmdata(plot_data$probablity, plot_data$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_single <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "Liver"

roc <- pre1_df[pre1_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(aes(colour="red"))+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14))+
  scale_color_manual(values='#bc5148')+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p1 <- p + annotate("text", x = .65, y = .35, size=5, colour = "black",
                   label = paste("AUC =", round(auc_single$aucs[1], 4)))

p1



### show result in table

### Esophagus
data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_esop_sample$Sample)
data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))
data <- data %>% filter(type!="null")
churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
table_data_Esophagus <- data.frame(SampleHRstatus = data$type, ProbablityScore = churn.pred,
                        CancerType = "Esophagus")

### Kidney
data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_kidn_sample$Sample)
data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))
data <- data %>% filter(type!="null")
churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
table_data_Kidney <- data.frame(SampleHRstatus = data$type, ProbablityScore = churn.pred,
                                   CancerType = "Kidney")

### Medulloblastoma
data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_medu_sample$Sample)
data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))
data <- data %>% filter(type!="null")
churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
table_data_Medulloblastoma <- data.frame(SampleHRstatus = data$type, ProbablityScore = churn.pred,
                                   CancerType = "Medulloblastoma")

### Skin
data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_skin_sample$Sample)
data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))
data <- data %>% filter(type!="null")
churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
table_data_Skin <- data.frame(SampleHRstatus = data$type, ProbablityScore = churn.pred,
                                   CancerType = "Skin")

### Stomach
data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_stom_sample$Sample)
data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))
data <- data %>% filter(type!="null")
churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
table_data_Stomach <- data.frame(SampleHRstatus = data$type, ProbablityScore = churn.pred,
                                   CancerType = "Stomach")


### ALL
table_data_allcancer <- rbind(table_data_Esophagus, table_data_Kidney,
                              table_data_Medulloblastoma, table_data_Skin, table_data_Stomach)

library(writexl)
write_xlsx(table_data_allcancer, "./supp/data/table_data_allcancer.xlsx")


