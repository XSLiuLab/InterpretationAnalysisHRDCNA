rm(list=ls())

# GIS vs HRDcna —— validation

library(precrec)
library(dplyr)
library(ggplot2)

setwd("~/HRD/HRDCNA/")

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")

HRD_60array <- readRDS("~/HRD/yhz_CNHRD/data_new/HRDdata/HRD_60array.rds")
HRD_60array$SampleID <- paste(HRD_60array$SampleID, "a", sep = "")
HRD_60wgs <- readRDS("~/HRD/yhz_CNHRD/data_new/HRDdata/HRD_60wgs.rds")
HRD_panel <- readRDS("~/HRD/yhz_CNHRD/data_new/HRDdata/HRD_panel.rds")



tally_W_60array <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_60array.rds")
tally_W_60wgs <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_60wgs.rds")
tally_W_panel <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_panel.rds")

wgs66_hrd <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/wgs67_hrd_PALB2.rds")
wgs66_hrr <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/wgs67_hrr_PALB2.rds")
snp_hrd <- wgs66_hrd
snp_hrd$Sample <- paste(snp_hrd$Sample, "a", sep = "")
snp_hrr <- wgs66_hrr
snp_hrr$Sample <- paste(snp_hrr$Sample, "a", sep = "")
panel_all_hrr <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/panel_all_hrr.rds")
panel_all_hrd <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/panel_all_hrd.rds")


HRD_all <- rbind(HRD_60wgs, HRD_60array, HRD_panel)


nmfwgs <- as.data.frame(tally_W_60wgs$nmf_matrix)
nmfsnp <- as.data.frame(tally_W_60array$nmf_matrix)
rownames(nmfsnp) <- paste(rownames(nmfsnp), "a", sep = "")
nmfpanel <- as.data.frame(tally_W_panel$nmf_matrix)

data <- rbind(nmfwgs, nmfsnp, nmfpanel)
data$type <- ifelse(rownames(data) %in% wgs66_hrd$Sample, "1",
                    ifelse(rownames(data) %in% wgs66_hrr$Sample, "0",
                           ifelse(rownames(data) %in% snp_hrd$Sample, "1",
                                  ifelse(rownames(data) %in% snp_hrr$Sample, "0", "null"))))
data <- data %>% filter(type != "null")


churn.pred <- predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.pred <- as.data.frame(churn.pred)
churn.pred$SampleID <- rownames(data)
churn.pred$type <- data$type

data <- inner_join(HRD_all, churn.pred, by="SampleID")



churn.roc = pROC::roc(data$type, data$churn.pred)
churn.roc$auc
# Area under the curve: 0.9788
churn.roc = pROC::roc(data$type, data$TAI)
churn.roc$auc
# Area under the curve: 0.9609
churn.roc = pROC::roc(data$type, data$LST)
churn.roc$auc
# Area under the curve: 0.9952
churn.roc = pROC::roc(data$type, data$LOH)
churn.roc$auc
# Area under the curve: 0.9243
churn.roc = pROC::roc(data$type, data$ALLHRD)
churn.roc$auc
# Area under the curve: 0.9977

# 十字交叉图
# data$HRDCNA_TYPE <- ifelse(data$churn.pred > 0.16, "HRD", "HRP")
# data$ALLHRD_TYPE <- ifelse(data$ALLHRD > 38, "HRD", "HRP")
#
# # predpcawg$gbm <- ifelse(predpcawg$churnpcawg >= 0.5, "HRD", "HRP")
# # predpcawg$chord <- ifelse(predpcawg$p_hrd >= 0.5, "HRD", "HRP")
# # predpcawg$type <- ifelse(predpcawg$response == "none", "0", "1")
#
# p <- ggplot(data,aes(x=churn.pred,y=ALLHRD,colour=type))+
#   geom_point(size=3)+
#   xlab("HRDCNA")+ylab("GIS")+
#   theme(title=element_text(size=15),
#         panel.background  = element_rect(fill="white",colour = "black",size = 1),
#         panel.grid.major  = element_blank(),
#         panel.grid.minor  = element_blank(),
#         legend.background = element_blank(),
#         legend.text = element_text(size = 15),
#         legend.title = element_blank(),
#         legend.position = c(.75,.25),
#         legend.key=element_blank(),
#         axis.text.x  =element_text(size=15,color = "black"),
#         axis.text.y = element_text(size=15,color = "black"))+
#   geom_hline(yintercept = 38)+
#   geom_vline(xintercept = 0.16)+
#   annotate("text",label="HRD",x=0.8,y=0.25,size=5,colour="#D74B4B")+
#   annotate("text",label="HRP",x=0.1,y=0.25,size=5,colour="#354B5E")+
#   scale_discrete_manual(values = c("#716e77", "#e79686"),
#                         aesthetics = "colour",
#                         labels=c("HRP","HRD"))
# p
#


pre_obj1 <- mmdata(data$churn.pred, data$type)
pre_obj1 <- evalmod(pre_obj1)
auchrdcna <- auc(pre_obj1)

pre_obj2 <- mmdata(data$TAI, data$type)
pre_obj2 <- evalmod(pre_obj2)
auctai <- auc(pre_obj2)

pre_obj3 <- mmdata(data$LST, data$type)
pre_obj3 <- evalmod(pre_obj3)
auclst <- auc(pre_obj3)

pre_obj4 <- mmdata(data$LOH, data$type)
pre_obj4 <- evalmod(pre_obj4)
aucloh <- auc(pre_obj4)

pre_obj5 <- mmdata(data$ALLHRD, data$type)
pre_obj5 <- evalmod(pre_obj5)
aucall <- auc(pre_obj5)


pre1_df <- fortify(pre_obj1)
pre2_df <- fortify(pre_obj2)
pre3_df <- fortify(pre_obj3)
pre4_df <- fortify(pre_obj4)
pre5_df <- fortify(pre_obj5)

pre1_df$Dataset <- "HRDCNA"
pre2_df$Dataset <- "TAI"
pre3_df$Dataset <- "LST"
pre4_df$Dataset <- "LOH"
pre5_df$Dataset <- "LOH+LST+TAI"

performance_df <- Reduce(rbind,list(pre1_df,pre2_df,pre3_df,pre4_df,pre5_df))


roc <- performance_df[performance_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y, group = Dataset)) +
  theme_bw()+
  geom_line(aes(color = Dataset))+
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
  scale_color_manual(values=c('#bc5148','#4a4266',"#549688","#d0a727","#5B6049"))
p1 <- p + annotate("text",x = .65, y = .35, size=5,colour="#bc5148",
                   label = paste("HRDCNA =",round(auchrdcna$aucs[1],4))) +
  annotate("text",x = .65, y = .30,size=5,colour="#4a4266",
           label=paste("TAI =",round(auctai$aucs[1],4)))+
  annotate("text",x = .65, y = .25,size=5,colour="#549688",
           label=paste("LST =", round(auclst$aucs[1],4)))+
  annotate("text",x = .65, y = .20,size=5,colour="#d0a727",
           label=paste("LOH =", round(aucloh$aucs[1],4)))+
  annotate("text",x = .65, y = .15,size=5,colour="#5B6049",
           label=paste("GIS =", round(aucall$aucs[1],4)))
p1


prc <-  performance_df[performance_df$curvetype == "PRC",]

p3 <- ggplot(prc, aes(x=x, y=y, group = Dataset)) +
  theme_bw()+
  geom_line(aes(color = Dataset))+
  xlab("Recall")+
  ylab("Precision")+
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14))+
  scale_color_manual(values=c('#bc5148','#4a4266',"#549688","#d0a727","#5B6049"))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 + annotate("text",x = .65, y = .35, size=5,colour="#bc5148",
                    label = paste("HRDCNA =",round(auchrdcna$aucs[2],4))) +
  annotate("text",x = .65, y = .30,size=5,colour="#4a4266",
           label=paste("TAI =",round(auctai$aucs[2],4)))+
  annotate("text",x = .65, y = .25,size=5,colour="#549688",
           label=paste("LST =", round(auclst$aucs[2],4)))+
  annotate("text",x = .65, y = .20,size=5,colour="#d0a727",
           label=paste("LOH =", round(aucloh$aucs[2],4)))+
  annotate("text",x = .65, y = .15,size=5,colour="#5B6049",
           label=paste("GIS =", round(aucall$aucs[2],4)))
p4



