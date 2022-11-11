
# model vs --- Sig-CNS, Sig-CX, Feature

## ROC

feature <- read.table("./data/modelselection/model_selection_feature.txt", header = T)
sig_cns <- read.table("./data/modelselection/model_selection_cns.txt", header = T)
sig_cx <- read.table("./data/modelselection/model_selection_cx.txt", header = T)

feature$model <- "CNA Feature"
sig_cns$model <- "Sig-CNS"
sig_cx$model <- "Sig-CX"


alldata <- rbind(feature, sig_cns, sig_cx)
alldata_test <- alldata %>% filter(dataset == "testdata")
alldata_test$ROC <- as.numeric(alldata_test$ROC)
alldata_test$PRC <- as.numeric(alldata_test$PRC)

# saveRDS(alldata_test, file = "./data/modelselection/modelvs.rds")

p1 <- alldata_test %>%
  ggplot(aes(x = model, y = ROC, fill = model)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "AUC")+
  scale_x_discrete(name = "") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 2, face = "bold"),
        text = element_text(size = 12),
        axis.text.x=element_text(size = 11))+
  labs(fill = "Model")+
  stat_compare_means(comparisons = list( c("CNA Feature", "Sig-CNS"),
                                         c("Sig-CNS", "Sig-CX"),
                                         c("CNA Feature", "Sig-CX")),
                     method = "wilcox.test", label = "p.signif",
                     size=3,label.y=0.99)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E", "#E8CA78"))+
  coord_cartesian(ylim = c(0.9, 1))
p1




## PRC
p2 <- alldata_test %>%
  ggplot(aes(x = model, y = PRC, fill = model)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "PR-AUC")+
  scale_x_discrete(name = "") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.text.x=element_text(size = 11))+
  labs(fill = "Model")+
  stat_compare_means(comparisons = list( c("Feature", "Sig-CNS"),
                                         c("Sig-CNS", "Sig-CX"),
                                         c("Feature", "Sig-CX")),
                     method = "wilcox.test", label = "p.signif",
                     size=3,label.y=0.95)+
  scale_fill_manual(values = c("#D74B4B", "#354B5E", "#E8CA78"))+
  coord_cartesian(ylim = c(0.6, 1))
p2

