rm(list=ls())

# HR-gene —— Count

library(tidyverse)


churn.gbmtest4 <- readRDS("./data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("./data/modeldata/bestTree.rds")

tally_W_tcga <- readRDS("./data/tallydata/tally_W_tcga.rds")



nmftcga <- as.data.frame(tally_W_tcga$nmf_matrix)
nmftcga$sample <- rownames(nmftcga)

data.pred <- predict(churn.gbmtest4, nmftcga, n.trees = bestTree, type = "response")
data.pred <- as.data.frame(data.pred)
data.pred$sampleid <- nmftcga$sample
data.pred$sample <- substr(data.pred$sampleid, 1, 12)
tcgaalldatapred <- data.pred


allmat <- data.table::fread("./data/download/mut_tcga/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt")##come from http://www.github.com/riazn/biallelic_hr
mat <- t(allmat)
colnames(mat) <- mat[1,]
mat <- mat[-1,]
mat2 <- apply(mat, 2, function(x){as.numeric(x)})
mat2 <- as.data.frame(mat2)
rownames(mat2) <- rownames(mat)
mat2$sample <- rownames(mat2)


TCGAalldata <- tcgaalldatapred %>% filter(tcgaalldatapred$sample %in% rownames(mat2)) # 8035

allbrca <- mat2 %>% filter(BRCA1 != 0 | BRCA2 != 0) %>% select(c(1, 2)) # 369
allbrca$sample <- rownames(allbrca)

nobrca <- mat2 %>% filter(!mat2$sample %in% allbrca$sample) # 7809

nobrca <- inner_join(TCGAalldata, nobrca, by="sample") # 7668

colnames(nobrca)[1] <- "Predictedscore"


# HRD
HRDtcga <- nobrca %>% filter(Predictedscore > 0.2) # 1171

condition <- c(1,5)
# 1- Germline biallelic pathogenic with LOH
# 5- Somatic biallelic pathogenic with LOH
matrix_num <- HRDtcga[4:105]
logical_col <- apply(matrix_num, 2, function(x) all(!(x %in% condition)))
final_matrix <- matrix_num[!logical_col] # 33 genes

HRDtcgaother <- cbind(HRDtcga[1:3], final_matrix)

HRDtcgaother <- HRDtcgaother[ , -2]


genehrd <- HRDtcgaother[ , c(3:35)]

num = as.data.frame(apply(genehrd, 2, function(x){table(x == 1 | x ==5)}))

mun <- as.data.frame(t(num))

mun$Gene <- rownames(mun)
mun <- mun[ , c(3,2)]
rownames(mun) <- NULL
colnames(mun)[2] <- "Count"

mun$Count <- as.numeric(mun$Count)

mun$Relative_Count <- mun$Count/1171
mun$Type <- "HRD"

# saveRDS(mun, file = "./data/brcaness/hrdsample_count.rds")



# HRP
HRPtcga <- nobrca %>% filter(Predictedscore < 0.2) # 6497

condition <- c(1,5)
# 1- Germline biallelic pathogenic with LOH
# 5- Somatic biallelic pathogenic with LOH
matrix_num <- HRPtcga[4:105]
logical_col <- apply(matrix_num, 2, function(x) all(!(x %in% condition)))
final_matrix <- matrix_num[!logical_col] # 60 genes

HRPtcgaother <- cbind(HRPtcga[1:3], final_matrix)

HRPtcgaother <- HRPtcgaother[ , -2]


genehrd <- HRPtcgaother[ , c(3:62)]


num = as.data.frame(apply(genehrd, 2, function(x){table(x == 1 | x ==5)}))

mun_hrr <- as.data.frame(t(num))
mun_hrr$Gene <- rownames(mun_hrr)
mun_hrr <- mun_hrr[ , c(3,2)]
rownames(mun_hrr) <- NULL
colnames(mun_hrr)[2] <- "Count"

mun_hrr$Count <- as.numeric(mun_hrr$Count)
mun_hrr$Relative_Count <- mun_hrr$Count/6497
mun_hrr$Type <- "HRP"

# saveRDS(mun_hrr, file = "./data/brcaness/hrpsample_count.rds")


munall <- rbind(mun, mun_hrr)

# saveRDS(munall, file = "./data/brcaness/allsample_count.rds")



position_dodge(width = NULL, preserve = c("total", "single"))

position_dodge2(
  width = NULL,
  preserve = c("total", "single"),
  padding = 0.1,
  reverse = FALSE
)


p1 <- ggplot(munall, aes(x = reorder(Gene, -Relative_Count), y = Relative_Count, fill = Type)) +
  geom_bar(stat = 'identity', show.legend = T, position = position_dodge(preserve = 'single')) +
  theme_bw() + xlab("")+ylab("Probability of Biallelic Pathogenic with LOH") + theme_classic() +
  scale_fill_manual(name='HR Status', values=c('#a64942','#4d638c')) +
  theme(axis.text.x  = element_text(size = 8,color = "black", angle = 90,vjust = 0.6),
        axis.text.y = element_text(size=10,color = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size=13),
        axis.line = element_line(colour="black"),
        legend.position = "right",legend.text = element_text(size = 11),
        legend.title = element_text(size = 10))
p1



# Sort manually by odds ratio
hrdsample_count <- readRDS("./data/brcaness/hrdsample_count.rds")
hrpsample_count <- readRDS("./data/brcaness/hrpsample_count.rds")

colnames(hrdsample_count) <- c("Gene", "Count_HRD", "Relative_Count_HRD", "Type_HRD")
colnames(hrpsample_count) <- c("Gene", "Count_HRP", "Relative_Count_HRP", "Type_HRP")

all <- full_join(hrdsample_count, hrpsample_count, by = "Gene")
all$OR <- all$Relative_Count_HRD/all$Relative_Count_HRP
all <- all[ , c(1,3,6,8)]

# saveRDS(all, file = "./data/brcaness/allsample_count_OR.rds")


