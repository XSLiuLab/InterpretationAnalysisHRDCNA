rm(list=ls())

# validation

library(tidyverse)
library(ggplot2)
library(caret)
library(pROC)
library(survival)
library(survminer)
library(gbm)
library(precrec)

setwd("~/HRD/HRDCNA/")

churn.gbmtest4 <- readRDS("./data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("./data/modeldata/bestTree.rds")

tally_W_panel <- readRDS("./data/tallydata/tally_W_panel.rds")
tally_W_wgs <- readRDS("./data/tallydata/tally_W_60wgs.rds")
tally_W_SNP <- readRDS("./data/tallydata/tally_W_60array.rds")
tally_W_10x <- readRDS("./data/tallydata/tally_W_60wgs10x.rds")
tally_W_15x <- readRDS("./data/tallydata/tally_W_60wgs15x.rds")
tally_W_30x <- readRDS("./data/tallydata/tally_W_60wgs30x.rds")

wgs66_hrd <- readRDS("./data/typedata/wgs67_hrd.rds")
wgs66_hrr <- readRDS("./data/typedata/wgs67_hrr.rds")
panel_all_hrr <- readRDS("./data/typedata/panel_all_hrr.rds")
panel_all_hrd <- readRDS("./data/typedata/panel_all_hrd.rds")



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

### ROC add confidence interval
roc <- roc(panel$sampletype, panel$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100) 
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255)) 


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
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)
p5




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

### ROC add confidence interval
roc <- roc(snp66$sampletype, snp66$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100) 
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255)) 


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
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)
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

### ROC add confidence interval
roc <- roc(wgs66$sampletype, wgs66$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100) 
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255)) 


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
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)
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

### ROC add confidence interval
roc <- roc(wgs30x66$sampletype, wgs30x66$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100) 
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255)) 


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
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)
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

### ROC add confidence interval
roc <- roc(wgs15x66$sampletype, wgs15x66$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100) 
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255)) 


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
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)
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

### ROC add confidence interval
roc <- roc(wgs10x66$sampletype, wgs10x66$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100) 
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255)) 


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
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)
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

### ROC add confidence interval
roc <- roc(alldata$sampletype, alldata$probablity, levels = c(0,1))
plot(roc)
sp.obj <- ci.sp(roc, sensitivities=seq(0, 1, .01), conf.level=0.95, boot.n=100) 
plot(sp.obj, type="shape", col = rgb(192,192,192, max=255)) 


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
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)

p5



