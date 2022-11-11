
rm(list=ls())

# Model Building

library(tidyverse)
library(sigminer)
library(pROC)
library(gbm)

### data preparing
all_cn_sigs_8_signature <- readRDS("./data/modelselection/data/all_cn_sigs_8_signature.rds")
pcawghrd <- readRDS("./data/typedata/pcawg_hrd.rds")
pcawghrr <- readRDS("./data/typedata/pcawg_hrr.rds")
a560hrd <- readRDS("./data/typedata/a560_hrd.rds")
a560hrr <- readRDS("./data/typedata/a560_hrr.rds")


#### all data 1470
alldata <- t(all_cn_sigs_8_signature$Exposure)
alldata <- as.data.frame(alldata)

alldata$type <- ifelse(rownames(alldata) %in% pcawghrd$sample, "1",
                       ifelse(rownames(alldata) %in% pcawghrr$sample, "0",
                              ifelse(rownames(alldata) %in% a560hrd$Sample, "1",
                                     ifelse(rownames(alldata) %in% a560hrr$Sample, "0", "null"))))
alldata <- alldata %>% filter(type != "null")


# saveRDS(alldata, file = "./data/modelselection/modelcnsdata/alldata.rds")

set.seed(1234) # for reproducibility
ind = sample(2, nrow(alldata), replace = T, prob = c(0.8, 0.2))

trainall = alldata[ind == 1, ] # #the training dataset 1174
testall = alldata[ind == 2, ] # #the test dataset 296

t <- table(trainall$type)
t[2]/t[1]
# 1
# 0.0961718
t <- table(testall$type)
t[2]/t[1]
# 1
# 0.1003717


saveRDS(trainall, file = "./data/modelselection/modelcnsdata/trainall.rds")
saveRDS(testall, file = "./data/modelselection/modelcnsdata/testall.rds")



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

  trainall <- readRDS("./data/modelselection/modelcnsdata/trainall.rds")
  testall <- readRDS("./data/modelselection/modelcnsdata/testall.rds")
  alldata <- readRDS("./data/modelselection/modelcnsdata/alldata.rds")

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

  write.table(model_selection, file="./data/modelselection/model_selection_cns.txt", append=T, row.names = F)

}

nums <- 1:100

parLapply(clus, nums ,fun = Run)


