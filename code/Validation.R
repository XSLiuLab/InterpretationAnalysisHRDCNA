rm(list=ls())

# validation

library(tidyverse)
library(ggplot2)
library(caret)
library(pROC)
library(survival)
library(survminer)

setwd("~/HRD/HRDCNA/")

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")

tally_W_panel <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_panel.rds")
tally_W_wgs <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_60wgs.rds")
tally_W_SNP <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_60array.rds")
tally_W_10x <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_60wgs10x.rds")
tally_W_15x <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_60wgs15x.rds")
tally_W_30x <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_60wgs30x.rds")

wgs66_hrd <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/wgs67_hrd.rds")
wgs66_hrr <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/wgs67_hrr.rds")
panel_all_hrr <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/panel_all_hrr.rds")
panel_all_hrd <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/panel_all_hrd.rds")



## different Platforms
#### panel ----
data <- tally_W_panel$nmf_matrix
data <- as.data.frame(data)

data$type <- ifelse(rownames(data) %in% panel_all_hrd$`Sample ID`, "1",
                    ifelse(rownames(data) %in% panel_all_hrr$`Sample ID`, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9580

#### ROC PRC
panel <- data.frame(sampletype = data$type, probablity = churn.pred,
                    sampleid = rownames(data))

pre_obj1 <- mmdata(panel$probablity, panel$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_panel <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "501 panel"

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
                   label = paste("AUC =", round(auc_panel$aucs[1], 2)))

p1


PRC <-  pre1_df[pre1_df$curvetype == "PRC",]

p3 <- ggplot(PRC, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(colour="#bc5148")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=14))+
  scale_color_manual(values='#bc5148')+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text", x = .65, y = .35, size=5, colour = "black",
                    label = paste("PR-AUC =", round(auc_panel$aucs[2],2)))

p4


#### score
data <- tally_W_panel$nmf_matrix
data <- as.data.frame(data)

data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$Sample <- rownames(data)
data.pred$type <- ifelse(data.pred$Sample %in% panel_all_hrd$`Sample ID`, "HRD",
                         ifelse(data.pred$Sample %in% panel_all_hrr$`Sample ID`, "HRP", "null"))
data.pred <- data.pred %>% filter(type != "null")

colnames(data.pred)[1] <- c("probability")

data_reorder <- data.pred

data_reorder$type <- factor(data_reorder$type,
                            levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("501 Panel Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
p5


#### Survival Probability of PFS
data <- tally_W_panel$nmf_matrix
data <- as.data.frame(data)
data$Sample <- rownames(data)
data$type <- ifelse(data$Sample %in% panel_all_hrd$`Sample ID`, "HRD",
                    ifelse(data$Sample %in% panel_all_hrr$`Sample ID`, "HRP", "null"))
data <- data %>% filter(type != "null")
data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data$score <- data.pred
data.roc=pROC::roc(data$type, data.pred)
data.roc$auc
# Area under the curve: 0.9580


panel_85 <- readxl::read_excel("./data/download/type_panel/panel_85_12885_2022_9602_MOESM8_ESM.xlsx")
colnames(panel_85) <- panel_85[1, ]
panel_85 <- panel_85[-1, ]
panel_85 <- as.data.frame(panel_85)

data85 <- data %>% filter(data$Sample %in% panel_85$`Sample ID`)
data85$OS <- panel_85$`PFS status`
data85$OS.time <- panel_85$`PFS (days)`
data85$TYPE <- ifelse(data85$score > 0.16, "HRD", "HRP")

alldatasurvival <- data85[ , c("Sample", "TYPE", "OS", "OS.time")]
alldatasurvival$OS <- as.numeric(alldatasurvival$OS)
alldatasurvival$OS.time <- as.numeric(alldatasurvival$OS.time)

fit <- survfit(Surv(alldatasurvival$OS.time, alldatasurvival$OS)~TYPE, data = alldatasurvival)

ggsurvplot(fit, data = alldatasurvival, surv.median.line = "hv",
           pval = TRUE,
           xlab = "Time(days)", ylab = "Survival Probability of PFS",
           legend.title = "",
           legend.labs = c("HRD ", "HRP"),
           break.x.by = 100,
           color = "strata",
           palette = c("#bc5148", "#3090a1"))


### wgs66 SNP ----
data <- tally_W_SNP$nmf_matrix
data <- as.data.frame(data)

data$type <- ifelse(rownames(data) %in% wgs66_hrd$Sample, "1",
                    ifelse(rownames(data) %in% wgs66_hrr$Sample, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9838


#### ROC PRC
snp66 <- data.frame(sampletype = data$type, probablity = churn.pred,
                    sampleid = rownames(data))

pre_obj1 <- mmdata(snp66$probablity, snp66$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_snp66 <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "66 SNP array"

roc <- pre1_df[pre1_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(aes(colour="#bc5148"))+
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
                   label = paste("AUC =", round(auc_snp66$aucs[1], 2)))

p1


PRC <-  pre1_df[pre1_df$curvetype == "PRC",]

p3 <- ggplot(PRC, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(colour="#bc5148")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=14))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text", x = .65, y = .35, size=5, colour = "black",
                    label = paste("PR-AUC =", round(auc_snp66$aucs[2],2)))

p4

#### score
data <- tally_W_SNP$nmf_matrix
data<- as.data.frame(data)

data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$Sample <- rownames(data)
data.pred$type <- ifelse(data.pred$Sample %in% wgs66_hrd$Sample, "HRD",
                         ifelse(data.pred$Sample %in% wgs66_hrr$Sample, "HRP", "null"))

data.pred <- data.pred %>% filter(type != "null")

colnames(data.pred)[1] <- c("probability")

data_reorder <- data.pred

data_reorder$type <- factor(data_reorder$type,
                            levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("66 SNP array Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
p5


### wgs66 WGS ----
data <- tally_W_wgs$nmf_matrix
data <- as.data.frame(data)

data$type <- ifelse(rownames(data) %in% wgs66_hrd$Sample, "1",
                    ifelse(rownames(data) %in% wgs66_hrr$Sample, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9697


#### ROC PRC
wgs66 <- data.frame(sampletype = data$type, probablity = churn.pred,
                    sampleid = rownames(data))

pre_obj1 <- mmdata(wgs66$probablity, wgs66$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_wgs66 <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "66 WGS"

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
                   label = paste("AUC =", round(auc_wgs66$aucs[1], 2)))

p1


PRC <-  pre1_df[pre1_df$curvetype == "PRC",]

p3 <- ggplot(PRC, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(colour="#bc5148")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=14))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text", x = .65, y = .35, size=5, colour = "black",
                    label = paste("PR-AUC =", round(auc_wgs66$aucs[2],2)))

p4

#### score
data <- tally_W_wgs$nmf_matrix
data <- as.data.frame(data)

data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$Sample <- rownames(data)
data.pred$type <- ifelse(data.pred$Sample %in% wgs66_hrd$Sample, "HRD",
                         ifelse(data.pred$Sample %in% wgs66_hrr$Sample, "HRP", "null"))

data.pred <- data.pred %>% filter(type != "null")

colnames(data.pred)[1] <- c("probability")

data_reorder <- data.pred

data_reorder$type <- factor(data_reorder$type,
                            levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("66 WGS Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
p5


## shallow WGS
### wgs66 30x ----
data <- tally_W_30x$nmf_matrix
data <- as.data.frame(data)

data$type <- ifelse(rownames(data) %in% wgs66_hrd$Sample, "1",
                    ifelse(rownames(data) %in% wgs66_hrr$Sample, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9545

#### ROC PRC
wgs30x66 <- data.frame(sampletype = data$type, probablity = churn.pred,
                       sampleid = rownames(data))

pre_obj1 <- mmdata(wgs30x66$probablity, wgs30x66$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_wgs30x66 <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "66 WGS_30X"

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
                   label = paste("AUC =", round(auc_wgs30x66$aucs[1], 2)))

p1


PRC <-  pre1_df[pre1_df$curvetype == "PRC",]

p3 <- ggplot(PRC, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(colour="#bc5148")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=14))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text", x = .65, y = .35, size=5, colour = "black",
                    label = paste("PR-AUC =", round(auc_wgs30x66$aucs[2],2)))

p4

#### score
data <- tally_W_30x$nmf_matrix
data <- as.data.frame(data)

data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$Sample <- rownames(data)
data.pred$type <- ifelse(data.pred$Sample %in% wgs66_hrd$Sample, "HRD",
                         ifelse(data.pred$Sample %in% wgs66_hrr$Sample, "HRP", "null"))

data.pred <- data.pred %>% filter(type != "null")

colnames(data.pred)[1] <- c("probability")

data_reorder <- data.pred

data_reorder$type <- factor(data_reorder$type,
                            levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("66 WGS 30X Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
p5


### wgs66 15x ----
data <- tally_W_15x$nmf_matrix
data <- as.data.frame(data)

data$type <- ifelse(rownames(data) %in% wgs66_hrd$Sample, "1",
                    ifelse(rownames(data) %in% wgs66_hrr$Sample, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.8969


#### ROC PRC
wgs15x66 <- data.frame(sampletype = data$type, probablity = churn.pred,
                       sampleid = rownames(data))

pre_obj1 <- mmdata(wgs15x66$probablity, wgs15x66$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_wgs15x66 <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "66 WGS_15X"

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
                   label = paste("AUC =", round(auc_wgs15x66$aucs[1], 2)))

p1


PRC <-  pre1_df[pre1_df$curvetype == "PRC",]

p3 <- ggplot(PRC, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(colour="#bc5148")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=14))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text", x = .65, y = .35, size=5, colour = "black",
                    label = paste("PR-AUC =", round(auc_wgs15x66$aucs[2],2)))

p4

#### score
data <- tally_W_15x$nmf_matrix
data <- as.data.frame(data)

data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$Sample <- rownames(data)
data.pred$type <- ifelse(data.pred$Sample %in% wgs66_hrd$Sample, "HRD",
                         ifelse(data.pred$Sample %in% wgs66_hrr$Sample, "HRP", "null"))

data.pred <- data.pred %>% filter(type != "null")

colnames(data.pred)[1] <- c("probability")

data_reorder <- data.pred

data_reorder$type <- factor(data_reorder$type,
                            levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("66 WGS 15X Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
p5


### wgs66 10x ----
data <- tally_W_10x$nmf_matrix
data <- as.data.frame(data)

data$type <- ifelse(rownames(data) %in% wgs66_hrd$Sample, "1",
                    ifelse(rownames(data) %in% wgs66_hrr$Sample, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9489


#### ROC PRC
wgs10x66 <- data.frame(sampletype = data$type, probablity = churn.pred,
                       sampleid = rownames(data))

pre_obj1 <- mmdata(wgs10x66$probablity, wgs10x66$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_wgs10x66 <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

pre1_df$Dataset <- "66 WGS_10X"

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
                   label = paste("AUC =", round(auc_wgs10x66$aucs[1], 2)))

p1


PRC <-  pre1_df[pre1_df$curvetype == "PRC",]

p3 <- ggplot(PRC, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(colour="#bc5148")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=14))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text", x = .65, y = .35, size=5, colour = "black",
                    label = paste("PR-AUC =", round(auc_wgs10x66$aucs[2],2)))

p4

#### score
data <- tally_W_10x$nmf_matrix
data <- as.data.frame(data)

data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$Sample <- rownames(data)
data.pred$type <- ifelse(data.pred$Sample %in% wgs66_hrd$Sample, "HRD",
                         ifelse(data.pred$Sample %in% wgs66_hrr$Sample, "HRP", "null"))

data.pred <- data.pred %>% filter(type != "null")

colnames(data.pred)[1] <- c("probability")

data_reorder <- data.pred

data_reorder$type <- factor(data_reorder$type,
                            levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("66 WGS 10X samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
p5



### alldata ----

nmfpanel <- tally_W_panel$nmf_matrix
nmfpanel <- as.data.frame(nmfpanel)
nmfpanel$type <- ifelse(rownames(nmfpanel) %in% panel_all_hrd$`Sample ID`, "1",
                        ifelse(rownames(nmfpanel) %in% panel_all_hrr$`Sample ID`, "0", "null"))
nmfpanel <- nmfpanel %>% filter(type!="null")

nmfwgs <- tally_W_wgs$nmf_matrix
nmfwgs <- as.data.frame(nmfwgs)
nmfwgs$type <- ifelse(rownames(nmfwgs) %in% wgs66_hrd$Sample, "1",
                      ifelse(rownames(nmfwgs) %in% wgs66_hrr$Sample, "0", "null"))
nmfwgs <- nmfwgs %>% filter(type!="null")

nmfsnp <- tally_W_SNP$nmf_matrix
nmfsnp <- as.data.frame(nmfsnp)
nmfsnp$type <- ifelse(rownames(nmfsnp) %in% wgs66_hrd$Sample, "1",
                      ifelse(rownames(nmfsnp) %in% wgs66_hrr$Sample, "0", "null"))
nmfsnp <- nmfsnp %>% filter(type!="null")

data <- rbind(nmfpanel, nmfwgs, nmfsnp)
data$sample <- rownames(data)

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.roc = pROC::roc(data$type, churn.pred)
churn.roc$auc
# Area under the curve: 0.9613


#### ROC PRC
alldata <- data.frame(sampletype = data$type, probablity = churn.pred,
                      sampleid = rownames(data))

pre_obj1 <- mmdata(alldata$probablity, alldata$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auc_alldata <- auc(pre_obj1)

pre1_df <- fortify(pre_obj1)

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
                   label = paste("AUC =", round(auc_alldata$aucs[1], 2)))

p1


PRC <-  pre1_df[pre1_df$curvetype == "PRC",]

p3 <- ggplot(PRC, aes(x=x, y=y)) +
  theme_bw()+
  geom_line(colour="#bc5148")+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        title=element_text(size=14))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text", x = .65, y = .35, size=5, colour = "black",
                    label = paste("PR-AUC =", round(auc_alldata$aucs[2],2)))

p4

#### score
data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$Sample <- rownames(data)

##### 解决样本名重复问题
wgs66_hrd2 <- wgs66_hrd
wgs66_hrd2$id <- "1"
wgs66_hrd2$Sample <- str_c(wgs66_hrd2$Sample, wgs66_hrd2$id)
wgs66_hrd2 <- wgs66_hrd2[ , -3]

wgs66_hrr2 <- wgs66_hrr
wgs66_hrr2$id <- "1"
wgs66_hrr2$Sample <- str_c(wgs66_hrr2$Sample, wgs66_hrr2$id)
wgs66_hrr2 <- wgs66_hrr2[ , -3]


data.pred$type <- ifelse(data.pred$Sample %in% wgs66_hrd$Sample, "HRD",
                         ifelse(data.pred$Sample %in% wgs66_hrr$Sample, "HRP",
                                ifelse(data.pred$Sample %in% panel_all_hrd$`Sample ID`, "HRD",
                                       ifelse(data.pred$Sample %in% panel_all_hrr$`Sample ID`, "HRP",
                                              ifelse(data.pred$Sample %in% wgs66_hrd2$Sample, "HRD",
                                                     ifelse(data.pred$Sample %in% wgs66_hrr2$Sample, "HRP", "null"))))))

data.pred <- data.pred %>% filter(type != "null")

colnames(data.pred)[1] <- c("probability")

data_reorder <- data.pred

data_reorder$type <- factor(data_reorder$type, levels = c("HRD", "HRP"))

p5 <- ggplot(data = data_reorder, mapping = aes(x = reorder(Sample, probability), y = probability, fill = type)) +
  geom_bar(stat = "identity") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.x  = element_blank()) +
  xlab("633 Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))
p5



