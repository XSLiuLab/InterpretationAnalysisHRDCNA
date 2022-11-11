rm(list=ls())

# GIS vs HRDcna —— held-out

library(precrec)
library(dplyr)
library(ggplot2)

setwd("~/HRD/HRDCNA/")

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")

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
# Area under the curve: 0.9741
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

pre_obj5 <- mmdata(data$HRDALL, data$type)
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


