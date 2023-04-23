rm(list=ls())

# SS78

library(pROC)
library(ggpubr)
library(ggplot2)


churn.gbmtest4 <- readRDS("./data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("./data/modeldata/bestTree.rds")

trainall <- readRDS("./data/modeldata/trainall.rds")
testall <- readRDS("./data/modeldata/testall.rds")
alldata <- readRDS("./data/modeldata/alldata.rds")


tally_W_panel <- readRDS("./data/tallydata/tally_W_panel.rds")
tally_W_wgs <- readRDS("./data/tallydata/tally_W_60wgs.rds")
tally_W_snp <- readRDS("./data/tallydata/tally_W_60array.rds")

panel_all_hrd <- readRDS("./data/typedata/panel_all_hrd.rds")
panel_all_hrr <- readRDS("./data/typedata/panel_all_hrr.rds")
wgs67_hrd <- readRDS("./data/typedata/wgs67_hrd.rds")
wgs67_hrr <- readRDS("./data/typedata/wgs67_hrr.rds")


### ROC
testall$sample <- rownames(testall)
plot.roc(testall$type, testall$`SS[>7 & <=8]`, col='#4a4266',
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Held-out Dataset: AUC = %.1f%%",
         print.auc.y=25, print.auc.x=75)

trainall$sample <- rownames(trainall)
plot.roc(trainall$type, trainall$`SS[>7 & <=8]`, col="#bc5148",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Training Dataset: AUC = %.1f%%",
         print.auc.y=20, print.auc.x=75, add=T)



### validation
#### count
nmfpanel <- tally_W_panel$nmf_matrix
nmfpanel <- as.data.frame(nmfpanel)
nmfpanel$type <- ifelse(rownames(nmfpanel) %in% panel_all_hrd$`Sample ID`, "1",
                        ifelse(rownames(nmfpanel) %in% panel_all_hrr$`Sample ID`, "0", "null"))
nmfpanel <- nmfpanel %>% filter(type!="null")

nmfwgs <- tally_W_wgs$nmf_matrix
nmfwgs <- as.data.frame(nmfwgs)
nmfwgs$type <- ifelse(rownames(nmfwgs) %in% wgs67_hrd$Sample, "1",
                      ifelse(rownames(nmfwgs) %in% wgs67_hrr$Sample, "0", "null"))
nmfwgs <- nmfwgs %>% filter(type!="null")

nmfsnp <- tally_W_snp$nmf_matrix
nmfsnp <- as.data.frame(nmfsnp)
nmfsnp$type <- ifelse(rownames(nmfsnp) %in% wgs67_hrd$Sample, "1",
                      ifelse(rownames(nmfsnp) %in% wgs67_hrr$Sample, "0", "null"))
nmfsnp <- nmfsnp %>% filter(type!="null")

data <- rbind(nmfpanel, nmfwgs, nmfsnp)
data$sample <- rownames(data)

data <- data %>% pivot_longer(-c("type", "sample"), names_to = "feature",
                              values_to = "count")

data <- data %>% filter(feature=="SS[>7 & <=8]")
data$type <- ifelse(data$type == "0", "HRP", "HRD")
data_reorder <- data

data_reorder$type <- factor(data_reorder$type, levels = c("HRD","HRP"))


ggplot(data = data_reorder, mapping = aes(x = reorder(sample,count),y=count,fill=type)) +
  geom_bar(stat = "identity")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        panel.background  = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size=15,color = "black"),
        axis.line = element_line(colour="black"),
        axis.ticks.x  = element_blank())+
  xlab("633 Samples")+ylab("SS[>7 & <=8]  Count")+
  theme(title=element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom")+labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))


#### ROC

dtwgs <- nmfwgs
dtwgs$sample <- rownames(dtwgs)
dtsnp <- nmfsnp
dtsnp$sample <- rownames(dtsnp)
dtpanel <- nmfpanel
dtpanel$sample <- rownames(dtpanel)


plot.roc(dtwgs$type, dtwgs$`SS[>7 & <=8]`, col = "#bc5148",
         percent=T,
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="WGS: AUC = %.1f%%",
         print.auc.y=30, print.auc.x = 65)

plot.roc(dtsnp$type, dtsnp$`SS[>7 & <=8]`, col = "#4a4266",
         percent=T,thresholds="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="SNP array: AUC = %.1f%%",
         print.auc.y=25,  print.auc.x = 65, add=T)

plot.roc(dtpanel$type, dtpanel$`SS[>7 & <=8]`, col = "#549688",
         percent=T,thresholds="best",
         lwd=2, print.auc=T, print.auc.cex=1, print.auc.pattern="Panel: AUC = %.1f%%",
         print.auc.y=20, print.auc.x = 65, add=T)

