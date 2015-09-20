library(foreign)
library(glmertree)

# DATA PREPARATION
metadata <- read.dta("Database IPDMA CBT PHA Version 11.dta")
metadata[metadata == 999] <- NA
metadata[metadata == 888] <- NA
vars <- c("studyid", "Tx_group", "Age", "Gender", "education", "ComorbidAnxietyDisorder", "HRSDt0", "HRSDt1")
factors <- c("studyid", "Tx_group", "Gender", "education", "ComorbidAnxietyDisorder")
metadata$education <- factor(metadata$education, ordered = T)
for (i in 1:length(factors)) {metadata[,factors[i]] <- factor(metadata[,factors[i]])}
metadata <- metadata[vars] # select only relevant variables
metadata <- metadata[complete.cases(metadata[,vars]),] # select only complete data
metadata <- metadata[!metadata$Tx_group == "placebo",] # remove placebo observations
metadata$Tx_group <- factor(metadata$Tx_group)
summary(metadata)

## calculate trees
hrsd <- lm(HRSDt1 ~ HRSDt0, data = metadata) 
metadata$HRSDfit <- fitted(hrsd) 
lm_f <- lmtree(HRSDt1 ~ Tx_group + offset(HRSDfit) | Age + Gender +
                 education + ComorbidAnxietyDisorder + HRSDt0, data = metadata)
lm_o <- lmtree(HRSDt1 ~ Tx_group | Age + Gender +
                 education + ComorbidAnxietyDisorder + HRSDt0, offset = HRSDfit, data = metadata)
## when joint=TRUE, prediction for random effects becomes cumbersome
#lmer_f <- lmertree(HRSDt1 ~ Tx_group | (1 | studyid) + offset(HRSDfit) |
#                     Age + Gender + education + ComorbidAnxietyDisorder + HRSDt0,
#                   data = metadata, ranefstart = metadata$HRSDfit)
lmer_f <- lmertree(HRSDt1 ~ Tx_group | (1 | studyid) + offset(HRSDfit) |
                     Age + Gender + education + ComorbidAnxietyDisorder + HRSDt0,
                   data = metadata, ranefstart = metadata$HRSDfit, joint = F)
plot(lm_f); print(lm_f)
plot(lm_o); print(lm_o)
plot(lmer_f); print(lmer_f)

## calculate predictions
metadata$lmtreepred <- predict(lm_f, newdata = metadata)
metadata$lmertreepred <- predict(lmer_f$tree, newdata = metadata) + predict(lmer_f$lmer, newdata = metadata)
mean(metadata$HRSDt1); hist(metadata$HRSDt1)
mean(metadata$lmtreepred); hist(metadata$lmtreepred)
mean(metadata$lmertreepred); hist(metadata$lmertreepred)
cor(metadata$lmertreepred, metadata$HRSDt1)
cor(metadata$lmtreepred, metadata$HRSDt1)
cor(metadata$lmertreepred, metadata$lmtreepred)

## calculate amount of variance explained
lmer_f$varcorr[[1]][[1]] # variance component for random effects
var(metadata$HRSDt1) # total variance
lmer_f$varcor[[1]][[1]] / var(metadata$HRSDt1) # variance explained by random effects
var(metadata$lmertreepred) / var(metadata$HRSDt1) # Total variance explained by lmertree
var(metadata$lmtreepred) / var(metadata$HRSDt1) # variance explained by the lmtree
var(metadata$hrs) / var(metadata$HRSDt1) # variance explained by HRSDt0
cor(metadata$HRSDt1, metadata$lmertreepred)^2 # lmertree explains 11% of variance in outcome variable
cor(metadata$HRSDt1, metadata$lmtreepred)^2 # lmtree explains 9% of variance outcome variable



## 50 fold CV
metadata <- metadata[,vars]
library(peperr)
set.seed(32)
samplething <- resample.indices(n = 694, sample.n = 50, method = "cv")
testdata <- list()
traindata <- list()
lmertrees <- list()
lmtrees <- list()
for (i in 1:50) {
  print(i)
  trainids <- samplething$sample.index[[i]]
  testids <- samplething$not.in.sample[[i]]
  traindata[[i]] <- metadata[trainids,]
  testdata[[i]] <- metadata[testids,]
  # fit lm(HRSDt1~HRSDt0) and predictions, based on training data 
  traindata[[i]]$HRSDfit <- fitted(lm(HRSDt1 ~ HRSDt0, data = traindata[[i]]))
  testdata[[i]]$HRSDfit <- predict(lm(HRSDt1 ~ HRSDt0, data = traindata[[i]]), newdata = testdata[[i]])

  # grow lmtree
  lmtrees[[i]] <- lmtree(HRSDt1 ~ Tx_group + offset(HRSDfit) | Age + Gender + 
                           education + ComorbidAnxietyDisorder + HRSDt0, data = traindata[[i]])  
  # calculate lmtree predictions
  testdata[[i]]$lmtreepreds <- predict(lmtrees[[i]], newdata = testdata[[i]])

  ## grow lmertree
  lmertrees[[i]] <- lmertree(HRSDt1 ~ Tx_group | (1 | studyid) + offset(HRSDfit) | 
                              Age + Gender + education + ComorbidAnxietyDisorder + HRSDt0,
                            data = traindata[[i]], ranefstart = traindata[[i]]$HRSDfit, joint = F)
 
  # calculate lmertree predictions
  testdata[[i]]$lmertreepreds <- predict(lmertrees[[i]]$tree, newdata = testdata[[i]]) + 
    predict(lmertrees[[i]]$lmer, newdata = testdata[[i]])   
}  


## Assess 50 fold CV results
lmtreecors <- list()
lmertreecors <- list()
for (i in 1:50){
  lmertreecors[[i]] <- cor(testdata[[i]]$lmertreepreds, testdata[[i]]$HRSDt1)
  lmtreecors[[i]] <- cor(testdata[[i]]$lmtreepreds, testdata[[i]]$HRSDt1)  
}
mean(unlist(lmtreecors))
mean(unlist(lmertreecors))
var(unlist(lmtreecors))
var(unlist(lmertreecors))