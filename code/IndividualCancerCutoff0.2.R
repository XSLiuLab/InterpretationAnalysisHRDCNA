
### individual cancer cut-off score testing

library(tidyverse)
library(ggplot2)
library(readxl)

churn.gbmtest4 <- readRDS("./data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("./data/modeldata/bestTree.rds")

pca <- read_xlsx("./data/download/type_pcawg/media-1.xlsx", sheet = 1)

pcawg <- pca%>% filter(pca$group == "PCAWG" & pca$used_for_perf_eval == "TRUE")

pcawghrd <- pcawg %>% filter(hr_status == "HR_deficient")
pcawghrp <- pcawg %>% filter(hr_status == "HR_proficient")

hrd_sampleid <- as.data.frame(pcawghrd$sample)
colnames(hrd_sampleid) <- "SampleID"
hrp_sampleid <- as.data.frame(pcawghrp$sample)
colnames(hrp_sampleid) <- "SampleID"


### individual cancer

pcawg_all_cancertype <- pcawg[,c(2,14)]
colnames(pcawg_all_cancertype) <- c("Sample", "CancerType")


vali_brca_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Breast")
vali_paca_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Pancreas")
vali_prad_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Prostate")
vali_lymp_sample <- pcawg_all_cancertype %>% filter(pcawg_all_cancertype$CancerType == "Lymphoid")



tally_W_pcawg <- readRDS("~/HRD/HRDCNA/data/tallydata/tally_W_pcawg.rds")
nmf_pcawg <- as.data.frame(tally_W_pcawg$nmf_matrix)
nmf_pcawg <- as.data.frame(tally_W_pcawg$nmf_matrix)


#### BRCA

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_brca_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.pred <- as.data.frame(churn.pred)
churn.pred$Sample <- rownames(data)
churn.pred$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                          ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))

churn.pred <- churn.pred %>% filter(type != "null")

colnames(churn.pred)[1] <- c("probability")

data_reorder <- churn.pred

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
  xlab("Breast Cancer Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)

p5



#### PACA

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_paca_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.pred <- as.data.frame(churn.pred)
churn.pred$Sample <- rownames(data)
churn.pred$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                          ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))

churn.pred <- churn.pred %>% filter(type != "null")

colnames(churn.pred)[1] <- c("probability")

data_reorder <- churn.pred

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
  xlab("Pancreatic Cancer Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)

p5



#### PRAD

data <- nmf_pcawg %>% filter(rownames(nmf_pcawg) %in% vali_prad_sample$Sample)

data$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "1",
                    ifelse(rownames(data) %in% hrp_sampleid$SampleID, "0", "null"))
data <- data %>% filter(type!="null")

churn.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
churn.pred <- as.data.frame(churn.pred)
churn.pred$Sample <- rownames(data)
churn.pred$type <- ifelse(rownames(data) %in% hrd_sampleid$SampleID, "HRD",
                          ifelse(rownames(data) %in% hrp_sampleid$SampleID, "HRP", "null"))

churn.pred <- churn.pred %>% filter(type != "null")

colnames(churn.pred)[1] <- c("probability")

data_reorder <- churn.pred

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
  xlab("Prostate Cancer Samples") + ylab("Probability Score") +
  theme(title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") + labs(fill = "HR Status ")+
  scale_fill_manual(values = c("#e79686", "#716e77"))+
  geom_hline(yintercept = 0.2, linetype = "dashed", size = 0.2)

p5


