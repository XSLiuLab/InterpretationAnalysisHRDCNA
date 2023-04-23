rm(list=ls())

# Variants of Uncertain Significance

library(gbm)
library(ggBubbles)


setwd("~/HRD/HRDCNA/")

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")


### 560 breast 244/559 ----

tally_W_560 <- readRDS("./data/tallydata/tally_W_560.rds")

nmf560 <- tally_W_560$nmf_matrix
nmf560 <- as.data.frame(nmf560)
nmf560$sample <- rownames(nmf560)
rownames(nmf560) <- NULL

churn.pred <- predict(churn.gbmtest4, nmf560, n.trees = bestTree, type ="response")
pred <- as.data.frame(churn.pred)
pred$Sample <- nmf560$sample

supplement <- readxl::read_excel("~/HRD/yhz_CNHRD/data/referencepaper/breast560/supplement4.xlsx", sheet = 4, na = 'NA')
colnames(supplement) <- supplement[1, ]

test560 <- inner_join(pred, supplement, by = 'Sample')
colnames(test560)[1] <- c("Predicted_value")
colnames(test560)[9] <- "hgvsc"
test560 <- test560 %>% unite(position, Chromosome, Position, sep = "-")


annoteafile <- data.table::fread("~/HRD/yhz_CNHRD/data/referencepaper/clinvar/variant_summary.txt")
annoteafile37 <- annoteafile %>% filter(GeneSymbol == "BRCA1" | GeneSymbol == "BRCA2") %>%
  filter(Assembly == "GRCh37")
anno <- annoteafile37 %>% unite(position, Chromosome, Start,sep = "-") %>%
  select(c(2,3,5,7,18,19))
anno <- anno %>% separate(col = Name, into = c("Name", "mut"), sep = ":")
anno <- anno %>% separate(col = mut, into = c("hgvsc", "protein"), sep = " ")
anno$protein <- gsub("\\(", "", anno$protein)
anno$protein <- gsub("\\)", "", anno$protein)

testall <- inner_join(test560, anno, by = "hgvsc")

test <- testall %>% filter(position.x == position.y) %>% filter(Loss_of_Alternative_Allele == "Yes") %>%
  filter(ClinicalSignificance == "Uncertain significance")

# saveRDS(test, file = "./data/vus/vus_560_all.rds")


### TCGA ---- 10906

allmat <- data.table::fread("~/HRD/biallelic_hr/biallelic_hr-master/biallelic_hr-master/Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt")##come from http://www.github.com/riazn/biallelic_hr
mat <- t(allmat)
colnames(mat) <- mat[1,]
mat <- mat[-1,]
mat2 <- apply(mat, 2, function(x){as.numeric(x)})
mat2 <- as.data.frame(mat2)
rownames(mat2) <- rownames(mat)

allbrca <- mat2 %>% filter(BRCA1 != 0 | BRCA2 != 0) %>% select(c(1, 2)) # 369 samples
allbrca$sampleid <- rownames(allbrca)

tally_W_tcga <- readRDS("./data/tallydata/tally_W_tcga.rds")
nmftcga <- as.data.frame(tally_W_tcga$nmf_matrix)
nmftcga$sample <- rownames(nmftcga)
rownames(nmftcga) <- NULL


churn.pred <- predict(churn.gbmtest4, nmftcga, n.trees =bestTree, type ="response")
pred <- as.data.frame(churn.pred)
pred$Sample <- nmftcga$sample
pred$cancertype <- nmftcga$cancertype
tcgaalldatapred <- pred
tcgaalldatapred$sampleid <- substr(tcgaalldatapred$Sample, 1, 12)


allbrcatcga <- inner_join(tcgaalldatapred, allbrca, by = "sampleid")
gerbiallbrca <- allbrcatcga %>% filter(BRCA1 %in% c(1,2,3,4) | BRCA2 %in% c(1,2,3,4))
BRCATCGA <- allbrcatcga %>% filter(BRCA1 %in% c(5,6,7,8,9,10) | BRCA2 %in% c(5,6,7,8,9,10))


#### BRCA Somatic Mutations

maf2 <- data.table::fread("./data/download/mut_tcga/mc3.v0.2.8.PUBLIC.maf")
BRCAmaf <- maf2 %>% filter(substr(maf2$Tumor_Sample_Barcode, 1, 16) %in% BRCATCGA$Sample) %>%
  filter(Hugo_Symbol == "BRCA1" | Hugo_Symbol == "BRCA2")
BRCAmaf$Sample <- substr(BRCAmaf$Tumor_Sample_Barcode,1,16)
BRCAmaf<- left_join(BRCAmaf, BRCATCGA, by="Sample")
BRCAmaf <- BRCAmaf %>% unite(position, Chromosome, Start_Position, End_Position, sep = "-") # 354

annoteafile37 <- annoteafile %>% filter(GeneSymbol == "BRCA1" | GeneSymbol == "BRCA2") %>%
  filter(Assembly == "GRCh37")

anno <- annoteafile37 %>% unite(position, Chromosome, Start, Stop, sep = "-") %>%
  select(c(2,3,5,7,18,19))
anno <- anno %>% separate(col = Name, into = c("Name", "mut"), sep = ":")
anno <- anno %>% separate(col = mut, into = c("hgvsc", "protein"), sep = " ")
anno$protein <- gsub("\\(", "", anno$protein)
anno$protein <- gsub("\\)", "", anno$protein)


annBRCAmaf <- right_join(anno, BRCAmaf, by = "position")
annBRCAmaf <- annBRCAmaf %>% filter(hgvsc == HGVSc)

annBRCAuncer <- annBRCAmaf %>%
  filter(ClinicalSignificance == "Uncertain significance" | ClinicalSignificance=="not provided") # 67

path_tcga <- annBRCAmaf %>%
  filter(ClinicalSignificance == "Pathogenic" | ClinicalSignificance=="Pathogenic/Likely pathogenic") # 24

path_uncer <- annBRCAuncer %>% filter(Tumor_Sample_Barcode %in% path_tcga$Tumor_Sample_Barcode) # 5

annBRCAuncer <- annBRCAuncer %>% filter(!Tumor_Sample_Barcode %in% path_uncer$Tumor_Sample_Barcode) # 62

annBRCAuncer$brca1type <- ifelse(annBRCAuncer$BRCA1==1,"Germline biallelic with LOH",
                                 ifelse(annBRCAuncer$BRCA1==2,"Germline biallelic compound heterozygous",
                                        ifelse(annBRCAuncer$BRCA1==3,"Germline biallelic compound heterozygous",
                                               ifelse(annBRCAuncer$BRCA1==4,"Germline monoallelic ",
                                                      ifelse(annBRCAuncer$BRCA1==5,"Somatic biallelic with LOH",
                                                             ifelse(annBRCAuncer$BRCA1==6,"Somatic biallelic compund heterozygous",
                                                                    ifelse(annBRCAuncer$BRCA1==7,"Somatic biallelic compund heterozygous",
                                                                           ifelse(annBRCAuncer$BRCA1==8,"Somatic biallelic with LOH",
                                                                                  ifelse(annBRCAuncer$BRCA1==9,"Somatic monoallelic",
                                                                                         ifelse(annBRCAuncer$BRCA1==10,"Somatic monoallelic","wild type"))))))))))




annBRCAuncer$brca2type <- ifelse(annBRCAuncer$BRCA2==1,"Germline biallelic with LOH",
                                 ifelse(annBRCAuncer$BRCA2==2,"Germline biallelic compound heterozygous",
                                        ifelse(annBRCAuncer$BRCA2==3,"Germline biallelic compound heterozygous",
                                               ifelse(annBRCAuncer$BRCA2==4,"Germline mono allelic ",
                                                      ifelse(annBRCAuncer$BRCA2==5,"Somatic biallelic with LOH",
                                                             ifelse(annBRCAuncer$BRCA2==6,"Somatic biallelic compund heterozygous",
                                                                    ifelse(annBRCAuncer$BRCA2==7,"Somatic biallelic compund heterozygous",
                                                                           ifelse(annBRCAuncer$BRCA2==8,"Somatic biallelic with LOH",
                                                                                  ifelse(annBRCAuncer$BRCA2==9,"Somatic monoallelic",
                                                                                         ifelse(annBRCAuncer$BRCA2==10,"Somatic monoallelic","wild type"))))))))))


annBRCAuncer$brca1type <- paste("BRCA1",annBRCAuncer$brca1type,sep = "-")
annBRCAuncer$brca2type <- paste("BRCA2",annBRCAuncer$brca2type,sep = "-")
annBRCAuncer <- annBRCAuncer %>% unite(type,brca1type,brca2type,sep = "-")
annBRCAuncer$type <- gsub("BRCA1-wild type-","",annBRCAuncer$type)
annBRCAuncer$type <- gsub("-BRCA2-wild type","",annBRCAuncer$type)

colnames(annBRCAuncer)[121] <- "Predictedscore"


alltcgauncer <- annBRCAuncer[ , c(2,3,4,5,6,8,9,14,40,41,91,120,121,122,123,124,125)]
alltcgauncer <- na.omit(alltcgauncer)

alltcgauncer <- alltcgauncer[!duplicated(alltcgauncer$protein)]

alltcgauncerbiabrca1 <- alltcgauncer %>% filter(alltcgauncer$BRCA1 == 5 |
                                                  alltcgauncer$BRCA1 == 6 |
                                                  alltcgauncer$BRCA1 == 7 |
                                                  alltcgauncer$BRCA1 == 8)
alltcgauncerbiabrca2 <- alltcgauncer %>% filter(alltcgauncer$BRCA2 == 5 |
                                                  alltcgauncer$BRCA2 == 6 |
                                                  alltcgauncer$BRCA2 == 7 |
                                                  alltcgauncer$BRCA2 == 8)
alltcgauncerbia <- rbind(alltcgauncerbiabrca1, alltcgauncerbiabrca2)


# saveRDS(alltcgauncerbia, file = "./data/vus/vus_tcga_som_bia.rds") # 14


#### BRCA Germline Mutations

BRCAGERMLINE <- data.table::fread(file="./data/download/mut_tcga/PCA_pathVar_integrated_filtered_adjusted.tsv") # download from:https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Germline-AWG

brcagermline <- BRCAGERMLINE %>% filter(HUGO_Symbol == "BRCA1" | HUGO_Symbol == "BRCA2")
colnames(brcagermline)[1] <- "sampleid"

anngermline <- inner_join(tcgaalldatapred, brcagermline, by = "sampleid")

gerbrca<- allbrca %>% filter(BRCA1 %in% c(1,2,3,4) | BRCA2 %in% c(1,2,3,4))

anngermline2 <- inner_join(gerbrca, anngermline, by = "sampleid")
anngermline2 <- anngermline2[ , c(1,2,3,4,6,7,9:19,21,28,29,30,44,45)]
anngermline2 <- anngermline2 %>% unite(position, Chromosome, Start, Stop, sep="-")


anno <- annoteafile37 %>% unite(position, Chromosome, Start, Stop, sep = "-") %>% select(c(2,3,5,7,15,18,19,29,30))

anngerbrca <- inner_join(anno, anngermline2, by = "position") %>% select(1:7,10:20,8,9,20,21,22)

annuncerger <- anngerbrca %>% filter(ClinicalSignificance == "Uncertain significance" | ClinicalSignificance == "not provided")

path_tcga <- anngerbrca %>%
  filter(ClinicalSignificance == "Pathogenic" | ClinicalSignificance=="Pathogenic/Likely pathogenic")

path_uncer <- annuncerger %>% filter(sampleid %in% path_tcga$sampleid)

annuncerger <- annuncerger %>% filter(!sampleid %in% path_uncer$sampleid)

annuncerger$brca1type <- ifelse(annuncerger$BRCA1==1,"Germline biallelic with LOH",
                                ifelse(annuncerger$BRCA1==2,"Germline biallelic compound heterozygous",
                                       ifelse(annuncerger$BRCA1==3,"Germline biallelic compound heterozygous",
                                              ifelse(annuncerger$BRCA1==4,"Germline monoallelic ",
                                                     ifelse(annuncerger$BRCA1==5,"Somatic biallelic with LOH",
                                                            ifelse(annuncerger$BRCA1==6,"Somatic biallelic compund heterozygous",
                                                                   ifelse(annuncerger$BRCA1==7,"Somatic biallelic compund heterozygous",
                                                                          ifelse(annuncerger$BRCA1==8,"Somatic biallelic with LOH",
                                                                                 ifelse(annuncerger$BRCA1==9,"Somatic monoallelic",
                                                                                        ifelse(annuncerger$BRCA1==10,"Somatic monoallelic","wild type"))))))))))

annuncerger$brca2type <- ifelse(annuncerger$BRCA2==1,"Germline biallelic with LOH",
                                ifelse(annuncerger$BRCA2==2,"Germline biallelic compound heterozygous",
                                       ifelse(annuncerger$BRCA2==3,"Germline biallelic compound heterozygous",
                                              ifelse(annuncerger$BRCA2==4,"Germline mono allelic ",
                                                     ifelse(annuncerger$BRCA2==5,"Somatic biallelic with LOH",
                                                            ifelse(annuncerger$BRCA2==6,"Somatic biallelic compund heterozygous",
                                                                   ifelse(annuncerger$BRCA2==7,"Somatic biallelic compund heterozygous",
                                                                          ifelse(annuncerger$BRCA2==8,"Somatic biallelic with LOH",
                                                                                 ifelse(annuncerger$BRCA2==9,"Somatic monoallelic",
                                                                                        ifelse(annuncerger$BRCA2==10,"Somatic monoallelic","wild type"))))))))))



annuncerger$brca1type <- paste("BRCA1",annuncerger$brca1type,sep = "-")
annuncerger$brca2type <- paste("BRCA2",annuncerger$brca2type,sep = "-")
annuncerger <- annuncerger %>% unite(type,brca1type,brca2type,sep = "-")
annuncerger$type <- gsub("BRCA1-wild type-","",annuncerger$type)
annuncerger$type <- gsub("-BRCA2-wild type","",annuncerger$type)

annuncerger <- annuncerger %>% separate(Name,into = c("name","protein"),sep = " ")
annuncerger$protein <- gsub("\\(","",annuncerger$protein)
annuncerger$protein <- gsub("\\)","",annuncerger$protein)
annuncerger <- annuncerger %>% filter(!is.na(protein))
colnames(annuncerger)[12] <- "Predictedscore"
anngeruncer1 <- annuncerger %>% unite(protein, sampleid, protein, sep = " ")
anngeruncer1 <- anngeruncer1[!duplicated(anngeruncer1$protein)]
annuncerger <- anngeruncer1 %>% separate(protein, into = c("sampleid","protein"),sep = " ")


annuncergerbiabrca1 <- annuncerger %>% filter(annuncerger$BRCA1 == 1 |
                                                annuncerger$BRCA1 == 2 |
                                                annuncerger$BRCA1 == 3)
annuncergerbiabrca2 <- annuncerger %>% filter(annuncerger$BRCA2 == 1 |
                                                annuncerger$BRCA2 == 2 |
                                                annuncerger$BRCA2 == 3)
annuncergerbia <- rbind(annuncergerbiabrca1, annuncergerbiabrca2)

# saveRDS(annuncergerbia, file = "./data/vus/vus_tcga_ger_bia.rds")


som_vus <- alltcgauncerbia[, c(14,3,2,4,5,8,13)]
ger_vus <- annuncergerbia[, c(3,4,23,5,6,22,12)]
ger_vus[1,3] <- "c.5293_5294insA"
colnames(ger_vus)[3] <- "hgvsc"

tcga_allvus <- rbind(som_vus, ger_vus)
tcga_allvus <- tcga_allvus[-13]

# saveRDS(tcga_allvus, file = "./data/vus/vus_tcga_all.rds")



vus_tcga_all <- readRDS("~/HRD/HRDCNA/data/vus/vus_tcga_all.rds")
vus_560_all <- readRDS("~/HRD/HRDCNA/data/vus/vus_560_all.rds")
vus_560_all <- vus_560_all[, c(2,14,8,15,16,10,1)]
colnames(vus_560_all)[6] <- "Variant_Classification"
colnames(vus_560_all)[7] <- "Predictedscore"
colnames(vus_tcga_all)[1] <- "Sample"

allvus <- rbind(vus_560_all, vus_tcga_all)
allvus <- allvus[-19, ]

# saveRDS(vus_560_all, file = "./data/vus/VUS_560.rds")
# saveRDS(vus_tcga_all, file = "./data/vus/VUS_TCGA.rds")
# saveRDS(allvus, file = "./data/vus/VUS_ALL.rds")


all_vus <- allvus %>% unite(protein, Sample, protein, sep = " - ")


ggplot(all_vus, aes(x = Predictedscore, y = reorder(protein, -Predictedscore))) +
  geom_bar(aes(fill = GeneSymbol), stat ='identity') +
  theme_bw() + ylab("") + xlab("Probability Score")+ theme_classic() +
  scale_fill_manual(name='', values=c('#587498','#CC9900'))+
  theme(axis.text.x  = element_text(size = 10,color = "black"),
        axis.text.y = element_text(size=10,color = "black"),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size=10),
        axis.line = element_line(colour="black"),
        legend.position = "right",legend.text = element_text(size = 11),
        legend.title = element_text(size = 10))+
  ggtitle(" ")



