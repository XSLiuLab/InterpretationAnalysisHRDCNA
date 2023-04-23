rm(list=ls())

# Model Building

library(tidyverse)
library(sigminer)
library(pROC)
library(gbm)

### data preparing
snp560_dt <- readRDS("./data/modelselection/data/snp560_cx.rds")
pcawg_dt <- readRDS("./data/modelselection/data/pcawg_cx.rds")
pcawghrd <- readRDS("./data/typedata/pcawg_hrd.rds")
pcawghrr <- readRDS("./data/typedata/pcawg_hrr.rds")
a560hrd <- readRDS("./data/typedata/a560_hrd.rds")
a560hrr <- readRDS("./data/typedata/a560_hrr.rds")

#### all data 1899
alldata <- rbind(pcawg_dt, snp560_dt)
alldata <- as.data.frame(alldata)
alldata$sample <- rownames(alldata)
rownames(alldata) <- NULL

alldata$type <- ifelse(alldata$sample %in% pcawghrd$sample, "1",
                       ifelse(alldata$sample %in% pcawghrr$sample, "0",
                              ifelse(alldata$sample %in% a560hrd$Sample, "1",
                                     ifelse(alldata$sample %in% a560hrr$Sample, "0", "null"))))
alldata <- alldata %>% filter(type != "null")
alldata$sample <- as.factor(alldata$sample)
# 1371

# saveRDS(alldata, file = "alldata.rds")

set.seed(1234) # for reproducibility
ind = sample(2, nrow(alldata), replace = T, prob = c(0.8, 0.2))

trainall = alldata[ind == 1, ] # #the training dataset
testall = alldata[ind == 2, ] # #the test dataset

t <- table(trainall$type)
t[2]/t[1]
#        1
# 0.1539256
t <- table(testall$type)
t[2]/t[1]
#        1
# 0.1759259

saveRDS(trainall, file = "./data/modelselection/modelcxdata/trainall.rds")
saveRDS(testall, file = "./data/modelselection/modelcxdata/testall.rds")



### model building

library(parallel)

clus <- makeCluster(20)


Run <- function(x){

  library(tidyverse)
  library(sigminer)
  library(pROC)
  library(gbm)
  library(precrec)
  library(dplyr)

  trainall <- readRDS("./data/modelselection/modelcxdata/trainall.rds")
  testall <- readRDS("./data/modelselection/modelcxdata/testall.rds")
  alldata <- readRDS("./data/modelselection/modelcxdata/alldata.rds")

  trainall <- trainall[ , -c(7,14,18)]
  testall <- testall[ , -c(7,14,18)]
  alldata <- alldata[ , -c(7,14,18)]

  set.seed(x)
  churn.gbmtest = gbm(formula = type~., distribution = "bernoulli", data = trainall,
                      n.trees = 6000,
                      interaction.depth = 7,
                      shrinkage = 0.01,
                      cv.folds = 5,
                      bag.fraction = 0.8,
                      n.minobsinnode = 5)
  bestTree <- gbm.perf(churn.gbmtest, plot.it = T,
                       oobag.curve = FALSE,
                       overlay = TRUE,
                       method = "cv")

  set.seed(3)
  churn.gbmtest2 = gbm(formula = type~., distribution = "bernoulli",
                       data = trainall,
                       n.trees = bestTree,
                       interaction.depth = 7,
                       shrinkage = 0.01,
                       cv.folds = 5,
                       bag.fraction = 0.8,
                       n.minobsinnode = 3)

  rel_infs <- relative.influence(churn.gbmtest2, n.trees = bestTree)


  featuresall <- gsub("\\`", "", names(rel_infs)[rel_infs/sum(rel_infs) > 0.01])

  choosedata = trainall[ , c(featuresall, "type")]

  set.seed(2)
  churn.gbmtest3 = gbm(formula = type~., distribution = "bernoulli",
                       data = choosedata,
                       n.trees = bestTree,
                       interaction.depth = 7,
                       shrinkage = 0.01,
                       cv.folds = 5,
                       bag.fraction = 0.8,
                       n.minobsinnode = 3)

  set.seed(1)
  churn.gbmtest4 = gbm(formula = type~., distribution = "bernoulli",
                       data = choosedata,
                       n.trees = bestTree,
                       interaction.depth = 7,
                       shrinkage = 0.01,
                       bag.fraction = 0.8,
                       n.minobsinnode = 3)


  churntrain = predict(churn.gbmtest4, trainall, n.trees = bestTree, type = "response")

  churnall = predict(churn.gbmtest4, alldata, n.trees = bestTree, type = "response")

  churntest = predict(churn.gbmtest4, testall, n.trees = bestTree,type = "response")


  trainscore <- data.frame(sampletype = trainall$type, probablity = churntrain,
                           sampleid = rownames(trainall), datatype = "Training dataset")

  testscore <- data.frame(sampletype = testall$type, probablity = churntest,
                          sampleid = rownames(testall), datatype = "Testing dataset")

  alldatascore <- data.frame(sampletype = alldata$type, probablity = churnall,
                             sampleid = rownames(alldata), datatype = "All dataset")



  pre_obj1 <- mmdata(trainscore$probablity, trainscore$sampletype)
  pre_obj1 <- evalmod(pre_obj1)
  auctrain <- auc(pre_obj1)

  pre_obj2 <- mmdata(testscore$probablity, testscore$sampletype)
  pre_obj2 <- evalmod(pre_obj2)
  auctest <- auc(pre_obj2)

  pre_obj3 <- mmdata(alldatascore$probablity, alldatascore$sampletype)
  pre_obj3 <- evalmod(pre_obj3)
  aucall <- auc(pre_obj3)


  model_selection <- data.frame(ROC = c(round(auctrain$aucs[1],4), round(auctest$aucs[1],4), round(aucall$aucs[1],4)),
                                PRC = c(round(auctrain$aucs[2],4), round(auctest$aucs[2],4), round(aucall$aucs[1],4)),
                                dataset = c("traindata","testdata","alldata"))

  write.table(model_selection, file="./data/modelselection/model_selection_cx.txt", append=T, row.names = F)

}

nums <- 1:100

parLapply(clus, nums ,fun = Run)


