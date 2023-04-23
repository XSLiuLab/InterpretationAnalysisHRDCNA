
### The difference in 10 CNA features between HRD and HRP samples in individual cancers


churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")


tally_W_pcawg <- readRDS("./data/tallydata/tally_W_pcawg.rds")
nmfpcawg <- as.data.frame(tally_W_pcawg$nmf_matrix)
pca <- read_xlsx("~/HRD/HRDCNA0/data/download/type_pcawg/media-1.xlsx", sheet = 1)
pcawg <- pca%>% filter(pca$group == "PCAWG")
pcawghrd <- pcawg %>% filter(hr_status == "HR_deficient")
pcawghrp <- pcawg %>% filter(hr_status == "HR_proficient")
nmfpcawg$type <- ifelse(rownames(nmfpcawg) %in% pcawghrd$sample,"1","0")
pcawg <- pcawg[,c(2,14)]
vali_ov_sample <- pcawg %>% filter(pcawg$cancer_type == "Ovary")
vali_brca_sample <- readRDS("./supp/data/vali_brca_sample.rds")
vali_paca_sample <- readRDS("./supp/data/vali_paca_sample.rds")
vali_prad_sample <- readRDS("./supp/data/vali_prad_sample.rds")


feat <- summary(churn.gbmtest4)
feat$var <- gsub("\\`","",feat$var)


fea <- nmfpcawg[ , c(which(colnames(nmfpcawg) %in% feat$var), 81)]
fea$sample <- rownames(fea)
fea$type <- ifelse(fea$type=="0","HRP","HRD")
fea <- fea %>% pivot_longer(-c(sample, type), names_to = "feature", values_to = "count")


fea_brca <- fea %>% filter(fea$sample %in% vali_brca_sample$Sample)
fea_paca <- fea %>% filter(fea$sample %in% vali_paca_sample$Sample)
fea_prad <- fea %>% filter(fea$sample %in% vali_prad_sample$Sample)
fea_ov <- fea %>% filter(fea$sample %in% vali_ov_sample$sample)




ggplot(fea_brca, aes(x = type, y = count, fill = type)) +
  geom_boxplot() + facet_wrap(~feature, ncol=5, scales = "free_y")+
  ylab("Count of Features") + xlab("Breast Cancer Samples")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "HR Status ")+
  geom_signif(comparisons = list(c("HRP","HRD")),
              map_signif_level = T, test = wilcox.test, y_position = c(80,30),
              tip_length = c(0.06,0.06))+
  cowplot::theme_cowplot(font_size = 15,line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))

fea_q <- fea_brca %>% group_split(fea_brca$feature)
for (i in 1:length(fea_q)) {
  dt <- compare_means(count ~ type, data = fea_q[[i]])
  dt$feature <- fea_q[[i]]$feature[i]
  write.table(dt, file="./supp/data/individual_brca_feature_qvalue.txt", append=T, row.names = F)
}

individual_brca_feature_qvalue <- read.table("./supp/data/individual_brca_feature_qvalue.txt", header = T)
individual_brca_feature_qvalue <- individual_brca_feature_qvalue %>% filter(individual_brca_feature_qvalue$.y. != ".y.")
individual_brca_feature_qvalue$cancertype <- "Breast"

saveRDS(individual_brca_feature_qvalue, file = "./supp/data/individual_brca_feature_qvalue.rds")


ggplot(fea_paca, aes(x = type, y = count, fill = type)) +
  geom_boxplot() + facet_wrap(~feature, ncol=5, scales = "free_y")+
  ylab("Count of Features") + xlab("Pancreatic Cancer Samples")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "HR Status ")+
  geom_signif(comparisons = list(c("HRP","HRD")),
              map_signif_level = T, test = wilcox.test, y_position = c(80,30),
              tip_length = c(0.06,0.06))+
  cowplot::theme_cowplot(font_size = 15,line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))


fea_q <- fea_paca %>% group_split(fea_paca$feature)
for (i in 1:length(fea_q)) {
  dt <- compare_means(count ~ type, data = fea_q[[i]])
  dt$feature <- fea_q[[i]]$feature[i]
  write.table(dt, file="./supp/data/individual_paca_feature_qvalue.txt", append=T, row.names = F)
}

individual_paca_feature_qvalue <- read.table("./supp/data/individual_paca_feature_qvalue.txt", header = T)
individual_paca_feature_qvalue <- individual_paca_feature_qvalue %>% filter(individual_paca_feature_qvalue$.y. != ".y.")
individual_paca_feature_qvalue$cancertype <- "Pancreatic"

saveRDS(individual_paca_feature_qvalue, file = "./supp/data/individual_paca_feature_qvalue.rds")



ggplot(fea_prad, aes(x = type, y = count, fill = type)) +
  geom_boxplot() + facet_wrap(~feature, ncol=5, scales = "free_y")+
  ylab("Count of Features") + xlab("Prostate Cancer Samples")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "HR Status ")+
  geom_signif(comparisons = list(c("HRP","HRD")),
              map_signif_level = T, test = wilcox.test, y_position = c(80,30),
              tip_length = c(0.06,0.06))+
  cowplot::theme_cowplot(font_size = 15,line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))

fea_q <- fea_prad %>% group_split(fea_prad$feature)
for (i in 1:length(fea_q)) {
  dt <- compare_means(count ~ type, data = fea_q[[i]])
  dt$feature <- fea_q[[i]]$feature[i]
  write.table(dt, file="./supp/data/individual_prad_feature_qvalue.txt", append=T, row.names = F)
}

individual_prad_feature_qvalue <- read.table("./supp/data/individual_prad_feature_qvalue.txt", header = T)
individual_prad_feature_qvalue <- individual_prad_feature_qvalue %>% filter(individual_prad_feature_qvalue$.y. != ".y.")
individual_prad_feature_qvalue$cancertype <- "Prostate"

saveRDS(individual_prad_feature_qvalue, file = "./supp/data/individual_prad_feature_qvalue.rds")


ggplot(fea_ov, aes(x = type, y = count, fill = type)) +
  geom_boxplot() + facet_wrap(~feature, ncol=5, scales = "free_y")+
  ylab("Count of Features") + xlab("Ovarian Cancer Samples")+
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "HR Status ")+
  geom_signif(comparisons = list(c("HRP","HRD")),
              map_signif_level = T, test = wilcox.test, y_position = c(80,30),
              tip_length = c(0.06,0.06))+
  cowplot::theme_cowplot(font_size = 15,line_size = 1)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E"))


fea_q <- fea_ov %>% group_split(fea_ov$feature)
for (i in 1:length(fea_q)) {
  dt <- compare_means(count ~ type, data = fea_q[[i]])
  dt$feature <- fea_q[[i]]$feature[i]
  write.table(dt, file="./supp/data/individual_ov_feature_qvalue.txt", append=T, row.names = F)
}

individual_ov_feature_qvalue <- read.table("./supp/data/individual_ov_feature_qvalue.txt", header = T)
individual_ov_feature_qvalue <- individual_ov_feature_qvalue %>% filter(individual_ov_feature_qvalue$.y. != ".y.")
individual_ov_feature_qvalue$cancertype <- "Ovarian"

saveRDS(individual_ov_feature_qvalue, file = "./supp/data/individual_ov_feature_qvalue.rds")



individual_cancer_feature_qvalue <- rbind(individual_brca_feature_qvalue, individual_paca_feature_qvalue, individual_prad_feature_qvalue, individual_ov_feature_qvalue)
saveRDS(individual_cancer_feature_qvalue, "./supp/data/individual_cancer_feature_qvalue.rds")










