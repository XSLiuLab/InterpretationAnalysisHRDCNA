rm(list=ls())

# Performance Verification

library(precrec)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("~/HRD/HRDCNA/")

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")

trainall <- readRDS("~/HRD/HRDCNA/data/modeldata/trainall.rds")
testall <- readRDS("~/HRD/HRDCNA/data/modeldata/testall.rds")


churntrain = predict(churn.gbmtest4, trainall, n.trees = bestTree, type = "response")
churntrain.roc = pROC::roc(trainall$type, churntrain)
churntrain.roc$auc
# Area under the curve: 1.0000

churntest = predict(churn.gbmtest4, testall, n.trees = bestTree,type = "response")
churntest.roc = pROC::roc(testall$type, churntest)
churntest.roc$auc
# Area under the curve: 0.9751

### ROC/PRC

trainscore <- data.frame(sampletype = trainall$type, probablity = churntrain,
                         sampleid = rownames(trainall), datatype = "Training Dataset")

testscore <- data.frame(sampletype = testall$type, probablity = churntest,
                        sampleid = rownames(testall), datatype = "Held-out Dataset")


pre_obj1 <- mmdata(trainscore$probablity, trainscore$sampletype)
pre_obj1 <- evalmod(pre_obj1)
auctrain <- auc(pre_obj1)

pre_obj2 <- mmdata(testscore$probablity, testscore$sampletype)
pre_obj2 <- evalmod(pre_obj2)
auctest <- auc(pre_obj2)

pre1_df <- fortify(pre_obj1)
pre2_df <- fortify(pre_obj2)


pre1_df$Dataset <- "Training Dataset"
pre2_df$Dataset <- "Held-out Dataset"


performance_df <- Reduce(rbind,list(pre1_df,pre2_df))


roc <- performance_df[performance_df$curvetype == "ROC",]

p <- ggplot(roc, aes(x=x, y=y, group = Dataset)) +
  theme_bw() +
  geom_line(aes(color = Dataset)) +
  xlab("1-Specificity") +
  ylab("Sensitivity") +
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14),
  ) +
  scale_color_manual(values=c('#bc5148','#4a4266'))
p1 <- p +
  annotate("text",x = .55, y = .25,size=5, colour="#bc5148",
           label=paste("Held-out Dataset: AUC =",round(auctest$aucs[1],2))) +
  annotate("text",x = .55, y = .15,size=5, colour="#4a4266",
           label=paste("Training Dataset: AUC =", round(auctrain$aucs[1],2)))
p1


prc <- performance_df[performance_df$curvetype == "PRC",]

p3 <- ggplot(prc, aes(x=x, y=y, group = Dataset)) +
  theme_bw() +
  geom_line(aes(color = Dataset))+
  xlab("Recall") +
  ylab("Precision") +
  theme(plot.title = element_text(hjust = 0.5),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  =  element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        axis.line = element_line(colour="black"),
        legend.position = "none",
        title=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14))+
  scale_color_manual(values=c('#bc5148','#4a4266'))+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

p4 <- p3 +
  annotate("text",x = .55, y = .25,size=5,colour="#bc5148",
           label=paste("Held-out Dataset: AUC-PR =",round(auctest$aucs[2],2))) +
  annotate("text",x = .55, y = .15,size=5,colour="#4a4266",
           label=paste("Training Dataset: AUC-PR =", round(auctrain$aucs[2],2)))
p4

# ------------------------------------------------------------------------------

### 1470个样本在模型中的预测得分
alldata <- readRDS("~/HRD/HRDCNA/data/modeldata/alldata.rds")

churn.pred = predict(churn.gbmtest4, alldata, n.trees = bestTree, type = "response")
pred <- as.data.frame(churn.pred)
pred$Sample <- rownames(alldata)

pred$type <- ifelse(alldata$type == 1, "HRD",
                    ifelse(alldata$type == 0, "HRP", "null"))
pred <- pred %>% filter(!type=="null")

colnames(pred)[1] <- c("probability")
data_reorder <- pred
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
  xlab("1470 Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))

p5


### cut-off score
library(parallel)
clus <- makeCluster(18)
clusterExport(clus, c("pred", "churn.pred", "alldata"), envir = environment())
Run <- function(y){
  library(caret)
  predict.class = ifelse(pred$probability >= y, "1", "0")
  predicetd <- as.factor(ifelse(churn.pred >= y, 1, 0))
  x <- as.factor(alldata$type)
  matrix <- confusionMatrix(predicetd, x)
  cutoff <- cbind(matrix$overall[1], matrix$byClass[1], matrix$byClass[2])
  colnames(cutoff) <- c("Accuracy", "Sensitivity", "Specificity")
  write.table(cutoff, file="./data/modeldata/cutoff.txt", append=T, row.names = F)
}
nums <- seq(from=0.1, to=0.7, by=0.01)
parLapply(clus, nums ,fun = Run)


cutoff <- read.table("./data/modeldata/cutoff.txt", header = T)
cutoff <- cutoff %>% filter(cutoff$Accuracy != "Accuracy")
cutoff$score <- nums

# saveRDS(cutoff, file = "./data/modeldata/cutoff.rds")

predicetd <- as.factor(ifelse(churn.pred >= 0.16, 1, 0))
x <- as.factor(alldata$type)
matrix <- confusionMatrix(predicetd, x)
matrix
# Confusion Matrix and Statistics
#
#             Reference
# Prediction    0    1
#          0 1329    6
#          1   11  124
#
# Accuracy : 0.9884
# 95% CI : (0.9815, 0.9932)
# No Information Rate : 0.9116
# P-Value [Acc > NIR] : <2e-16
#
# Kappa : 0.9295
#
# Mcnemar's Test P-Value : 0.332
#
#             Sensitivity : 0.9918
#             Specificity : 0.9538
#          Pos Pred Value : 0.9955
#          Neg Pred Value : 0.9185
#              Prevalence : 0.9116
#          Detection Rate : 0.9041
#    Detection Prevalence : 0.9082
#       Balanced Accuracy : 0.9728
#
#        'Positive' Class : 0
fourfoldplot(matrix$table, color = c("#5B6044", "#CF9B61"),
             conf.level = 0, margin = 1)











