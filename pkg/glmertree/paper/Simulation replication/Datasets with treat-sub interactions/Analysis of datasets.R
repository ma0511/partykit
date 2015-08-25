library(lme4)
library(partykit)
source("lmertree.R")

# Fit trees on training data
for (c in 1:50) {
  load(paste("datasets", c, sep="")) # load training data
  load(paste("descriptions", c, sep="")) # load trainig data descriptions
  REEMobtrees <- list()
  Mobtrees <- list()
  for (i in 1:length(datasets)) {
    print(i)
    if (descriptions[[i]][[3]] == "number of covariates = 5") {
      formula <- Y ~ T | X1 + X2 + X3 + X4 + X5
    }
    if (descriptions[[i]][[3]] == "number of covariates = 15") {
      formula <- Y ~ T | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + 
                  X11 + X12 + X13 + X14 + X15
    }
    REEMobtrees[[i]] <- lmertree(lmtreeformula = formula, data=datasets[[i]], 
                                  randomformula = Y~(1|bi), verbose = F, plotting=F, 
                                  ErrorTolerance = 0.001, maxdepth=4)
    Mobtrees[[i]] <- lmtree(formula = formula, data=datasets[[i]], maxdepth=4)
  }
  save(REEMobtrees, file=paste("REEMobtrees", c, sep=""))
  save(Mobtrees, file=paste("Mobtrees", c, sep=""))
}

# get tree characteristics
for(c in 1:50){
  print(c)
  load(paste("Mobtrees",c,sep=""))
  load(paste("REEMobtrees",c,sep=""))

  # Summarize tree characteristics of REEMob trees
  treesize.REEMobs <- list()
  for (i in 1:length(datasets)) {
    treesize.REEMobs[[i]] <- length(REEMobtrees[[i]]$Tree)
  }
  save(treesize.REEMobs, file=paste("treesize.REEMobs", c, sep=""))

  breaks.REEMobs <- list()
  varids.REEMobs <- list()
  for (i in 1:length(datasets)) {
    breaks.REEMobs[[i]] <- list()
    varids.REEMobs[[i]] <- list()
    for (j in 1:length(REEMobtrees[[i]][[1]])) {
      varids.REEMobs[[i]][[j]] <- REEMobtrees[[i]][[1]][[j]]$node$split$varid
      breaks.REEMobs[[i]][[j]] <- REEMobtrees[[i]][[1]][[j]]$node$split$breaks
    }    
  }    
  splits.REEMobs <- list()
  for (i in 1:length(datasets)) {
    splits.REEMobs[[i]] <- list()
    for (j in 1:length(varids.REEMobs[[i]])) {
      splits.REEMobs[[i]][[j]] <- c(varids.REEMobs[[i]][[j]], breaks.REEMobs[[i]][[j]])
    }
  }
  splits.REEMobs <- data.frame(varids=unlist(splits.REEMobs)[1:length(unlist(splits.REEMobs))%%2==1],
                               breaks=unlist(splits.REEMobs)[1:length(unlist(splits.REEMobs))%%2==0])
  save(splits.REEMobs, file=paste("splits.REEMobs", c, sep=""))

  # Summarize tree characteristics of Mob trees
  treesize.Mobs <- list()
  for (i in 1:length(datasets)) {
    treesize.Mobs[[i]] <- length(Mobtrees[[i]])
  }
  save(treesize.Mobs, file=paste("treesize.Mobs", c, sep=""))

  breaks.Mobs <- list()
  varids.Mobs <- list()
  for (i in 1:length(datasets)) {
    breaks.Mobs[[i]] <- list()
    varids.Mobs[[i]] <- list()
    for (j in 1:length(Mobtrees[[i]][[1]])) {
      varids.Mobs[[i]][[j]] <- Mobtrees[[i]][[1]][[j]]$node$split$varid
      breaks.Mobs[[i]][[j]] <- Mobtrees[[i]][[1]][[j]]$node$split$breaks
    }    
  }    
  splits.Mobs <- list()
  for (i in 1:length(datasets)) {
    if(length(varids.Mobs[[i]])>0) {
      splits.Mobs[[i]] <- list()
      for (j in 1:length(varids.Mobs[[i]])) {  
        splits.Mobs[[i]][[j]] <- c(varids.Mobs[[i]][[j]], breaks.Mobs[[i]][[j]])
      }
    }
  }
  splits.Mobs <- data.frame(varids=unlist(splits.Mobs)[1:length(unlist(splits.Mobs))%%2==1],
                            breaks=unlist(splits.Mobs)[1:length(unlist(splits.Mobs))%%2==0])
  save(splits.Mobs, file=paste("splits.Mobs", c, sep=""))
}


# Evaluate performance of trees with test data
for(c in 1:50){
  print(c)
  load(paste("Mobtrees",c,sep=""))
  load(paste("REEMobtrees",c,sep=""))
  load(paste("testdatasets", c, sep=""))
  load(paste("testdescriptions", c, sep=""))

  # compare treatment difference estimates of Mob and REEMob with true difference
  treatdiffs <- list()
  for (i in 1:length(testdata)) {
    diff <- as.numeric(substr(testdescriptions[[i]][[4]], 31, 33))
    tmp <- testdata[[i]][[1]]
    tmp$true_d_hat <- NA
    tmp$true_d_hat[tmp$X2<=30 & tmp$X1<=17] <- diff
    tmp$true_d_hat[tmp$X2<=30 & tmp$X1>17] <- 0
    tmp$true_d_hat[tmp$X2>30 & tmp$X1<=63] <- 0 
    tmp$true_d_hat[tmp$X2>30 & tmp$X5>63] <- -diff
    treatdiffs[[i]] <- data.frame(true_d_hat=tmp$true_d_hat)
    treatdiffs[[i]]$REEMob_d_hat <- NA
    treatdiffs[[i]]$REEMob_d_hat <- predict(REEMobtrees[[i]]$Tree, newdata=testdata[[i]][[1]], 
              type="response") - predict(REEMobtrees[[i]]$Tree, 
              newdata=testdata[[i]][[2]], type="response")
    treatdiffs[[i]]$Mob_d_hat <- predict(Mobtrees[[i]], newdata=testdata[[i]][[1]], 
              type="response") - predict(Mobtrees[[i]], 
              newdata=testdata[[i]][[2]], type="response")
  }
  save(treatdiffs, file=paste("treatdiffs", c, sep=""))

  # Compute true and predicted Ys
  true_pred_Y <- list()
  for (i in 1:length(testdata)) {
    true_pred_Y[[i]] <- data.frame(
      Ytrue = testdata[[1]][[3]]$Y,
      REEMobYpred = predict(REEMobtrees[[1]]$Tree, newdata=testdata[[1]][[3]]), 
      MobYpred = predict(Mobtrees[[1]], newdata=testdata[[1]][[3]]))
    }
  save(true_pred_Y, file=paste("true_pred_Y", c, sep=""))
}