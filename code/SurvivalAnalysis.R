

library(tidyverse)
library(survival)
library(survminer)
library(gbm)


### Kaplan-Meier survival analysis on patients receiving platinum chemotherapy

setwd("~/HRD/HRDCNA/")

churn.gbmtest4 <- readRDS("~/HRD/HRDCNA/data/modeldata/churn.gbmtest4.rds")
bestTree <- readRDS("~/HRD/HRDCNA/data/modeldata/bestTree.rds")

tally_W_panel <- readRDS("~/HRD/yhz_CNHRD/data_new/tallydata/tally_W_panel.rds")
panel_all_hrr <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/panel_all_hrr.rds")
panel_all_hrd <- readRDS("~/HRD/yhz_CNHRD/data_new/typedata/panel_all_hrd.rds")

panel_85 <- readxl::read_excel("./data/download/type_panel/panel_85_12885_2022_9602_MOESM8_ESM.xlsx")
colnames(panel_85) <- panel_85[1, ]
panel_85 <- panel_85[-1, ]
panel_85 <- as.data.frame(panel_85)


data <- tally_W_panel$nmf_matrix
data <- as.data.frame(data)
data$Sample <- rownames(data)
data$type <- ifelse(data$Sample %in% panel_all_hrd$`Sample ID`, "HRD",
                    ifelse(data$Sample %in% panel_all_hrr$`Sample ID`, "HRP", "null"))
data <- data %>% filter(type != "null")
data.pred = predict(churn.gbmtest4, data, n.trees = bestTree, type = "response")
data$score <- data.pred

data85 <- data %>% filter(data$Sample %in% panel_85$`Sample ID`)
data85.pred <- data85$score

data.roc=pROC::roc(data85$type, data85.pred)
data.roc$auc
# Area under the curve: 0.9837

pROC::coords(data.roc, "best")
#   threshold specificity sensitivity
# 1 0.1438813   0.9482759   0.9259259

data85$OS <- panel_85$`PFS status`
data85$OS.time <- panel_85$`PFS (days)`
data85$TYPE <- ifelse(data85$score > 0.14, "HRD", "HRP")

alldatasurvival <- data85[ , c("Sample", "TYPE", "OS", "OS.time")]
alldatasurvival$OS <- as.numeric(alldatasurvival$OS)
alldatasurvival$OS.time <- as.numeric(alldatasurvival$OS.time)

fit <- survfit(Surv(alldatasurvival$OS.time, alldatasurvival$OS)~TYPE, data = alldatasurvival)
ggsurvplot(fit, data = alldatasurvival, surv.median.line = "hv",
           pval = TRUE,
           xlab = "Time(days)", ylab = "Survival Probability",
           legend.title = "",
           legend.labs = c("HRD ", "HRP"),
           break.x.by = 100,
           color = "strata",
           palette = c("#bc5148", "#3090a1"))

