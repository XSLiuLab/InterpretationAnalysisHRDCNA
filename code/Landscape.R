rm(list=ls())

# HRD Landscape —— TCGA

library(ggplot2)

churn.gbmtest4 <- readRDS("./data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("./data/modeldata/bestTree.rds")

nmftcga_cancertype <- readRDS("./data/download/mut_tcga/nmftcga_cancertype.rds")
tally_W_tcga <- readRDS("./data/tallydata/tally_W_tcga.rds")

nmftcga <- as.data.frame(tally_W_tcga$nmf_matrix)
churn <- predict(churn.gbmtest4, nmftcga, n.trees =bestTree,type ="response")
pred <- data.frame(Predictedscore=churn)
pred$sample <- nmftcga_cancertype$sample
pred$cancertype <- nmftcga_cancertype$cancertype


ggplot(pred, aes(x=cancertype, y=Predictedscore, fill =cancertype )) +
  scale_fill_manual(values = rep("white", length(unique(pred$cancertype))))+
  geom_jitter(position = position_jitterdodge(jitter.width = .95,
                                              jitter.height = 0,
                                              dodge.width = .95),
              aes(fill = cancertype, col = cancertype), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),title = element_text(size=18),
        line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        axis.text.x  = element_text(size = 9,color = "black",angle = 90,vjust = 0.6),
        axis.text.y = element_text(size=10,color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size=15),
        axis.line = element_line(colour="black")) +
  ylab("Probability Score")+xlab("Cancer Type")+
  guides(fill="none", colour="none")+
  ggtitle("")

