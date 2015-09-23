library(lme4)
library(partykit)
source("lmertree.R")

for (c in 9:50) {
  load(paste("datasets", c, sep="")) # loads training data
  load(paste("descriptions", c, sep="")) # loads trainig data descriptions

  # Fit btrees on training data
  REEMobtrees <- list()
  Mobtrees <- list()
  for (i in 1:length(datasets)) {
  print(i)
  if (descriptions[[i]][[3]] == "number of covariates = 5") {
    formula <- Y ~ T | X1 + X2 + X3 + X4 + X5
  }
  if (descriptions[[i]][[3]] == "number of covariates = 15") {
    formula <- Y ~ T | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15
  }
  REEMobtrees[[i]] <- lmertree(lmtreeformula = formula, data=datasets[[i]], 
                              randomformula = Y~(1|bi), verbose = F, plotting=F, 
                              ErrorTolerance = 0.001, maxdepth=4)
  Mobtrees[[i]] <- lmtree(formula = formula, data=datasets[[i]],
                          maxdepth=4)
  }
  save(REEMobtrees, file=paste("REEMobtrees", c, sep=""))
  save(Mobtrees, file=paste("Mobtrees", c, sep=""))
}

for(c in 1:50){
  load(paste("Mobtrees",c,sep=""))
  load(paste("REEMobtrees",c,sep=""))
  # Summarize tree characteristics of REEMob trees
  treesize.REEMobs <- list()
  for (i in 1:length(datasets)) {
    treesize.REEMobs[[i]] <- length(REEMobtrees[[i]]$Tree)
  }
  save(treesize.REEMobs, file=paste("treesize.REEMobs", c, sep=""))
  # Summarize tree characteristics of Mob trees
  treesize.Mobs <- list()
  for (i in 1:length(datasets)) {
    treesize.Mobs[[i]] <- length(Mobtrees[[i]])
  }
  save(treesize.Mobs, file=paste("treesize.Mobs", c, sep=""))
}

# Evaluate tree sizes
REEMobtreesize <- vector()
Mobtreesize <- vector()
for (c in 1:50){
  load(paste("treesize.REEMobs", c, sep=""))
  REEMobtreesize <- c(REEMobtreesize, unlist(treesize.REEMobs))
  load(paste("treesize.Mobs", c, sep=""))
  Mobtreesize <- c(Mobtreesize, unlist(treesize.Mobs))  
}
treesizes <- data.frame(REEMobtreesize, Mobtreesize)

load("descriptions1")
for (i in 1:length(descriptions)) {
  treesizes[i,"N"] <- substr(descriptions[[i]][[1]], 5, 8) # N
  treesizes[i,"rho"] <- substr(descriptions[[i]][[2]], 7, 9) # is rho (correlation between covariates)
  treesizes[i,"np"] <- substr(descriptions[[i]][[3]], 24, 25) # is number of covariates
  treesizes[i,"treatdiff"] <- substr(descriptions[[i]][[4]], 31, 33) # is treatment effect difference
  treesizes[i,"corUbi"] <- substr(descriptions[[i]][[5]], 32, 64) # is correlation between U and bi
  treesizes[i,"numbclus"] <- substr(descriptions[[i]][[6]], 34, 35) # is number of random intercept values
  treesizes[i,"sigmabi"] <- substr(descriptions[[i]][[7]], 12, 14) # max and -min random intercept value    
}
save(treesizes,file="treesizes")
load("treesizes")

treesizes$N <- factor(treesizes$N)
treesizes$rho <-  factor(treesizes$rho)            
treesizes$np <- factor(treesizes$np)              
treesizes$treatdiff <- factor(treesizes$treatdiff)          
treesizes$corUbi <- factor(treesizes$corUbi)           
treesizes$numbclus <- factor(treesizes$numbclus)        
treesizes$sigmabi <- factor(treesizes$sigmabi)

treesizes$N <- rep(treesizes$N[1:length(descriptions)], times=1)
treesizes$rho <- rep(treesizes$rho[1:length(descriptions)], times=1)
treesizes$np <- rep(treesizes$np[1:length(descriptions)], times=1)
treesizes$treatdiff <- rep(treesizes$treatdiff[1:length(descriptions)], times=1)
treesizes$corUbi <- rep(treesizes$corUbi[1:length(descriptions)], times=1)
treesizes$numbclus <- rep(treesizes$numbclus[1:length(descriptions)], times=1)
treesizes$sigmabi <- rep(treesizes$sigmabi[1:length(descriptions)], times=1)

table(treesizes[,1]) # clustermob
prop.table(table(treesizes[,1])) # clustermob
table(treesizes[,2]) # mob
prop.table(table(treesizes[,2])) # mob

treesizes.long <- data.frame(treesize=stack(treesizes[,1:2]), 
                             N=rep(treesizes$N, 2), rho=rep(treesizes$rho, 2), np=rep(treesizes$np, 2), 
                             treatdiff=rep(treesizes$treatdiff, 2),  corUbm=rep(treesizes$corUbi, 2), 
                             numbclus=rep(treesizes$numbclus, 2), sigmabm=rep(treesizes$sigmabi, 2),
                             datasetID=factor(rep(1:nrow(treesizes), 2)))
tmp <- as.character(treesizes.long$treesize.ind)
tmp[treesizes.long$treesize.ind=="REEMobtreesize"] <- "GLMM tree"
tmp[treesizes.long$treesize.ind=="Mobtreesize"] <- "GLM tree"
treesizes.long$treesize.ind <- factor(tmp)

mean(treesizes[,1]);sd(treesizes[,1]) # REEMob
mean(treesizes[,2]);sd(treesizes[,2]) # mob

# create lattice xyplots
library(lattice)
treesize.anova <- aov(treesize.values ~ treesize.ind + N + rho + np + treatdiff + numbclus + sigmabm + corUbm +
                        treesize.ind*(N + rho + np + treatdiff + numbclus + sigmabm + corUbm), 
                      data=treesizes.long)
summary(treesize.anova)[[1]]["Sum Sq"] / sum(summary(treesize.anova)[[1]]["Sum Sq"] )
# N, sigmabm & corUbi have main and interaction effects with eta squared > .01
treesizes.long$sigmabm <- factor(as.numeric(as.character(treesizes.long$sigmabm)), ordered=T)
treesizes.long$N <- factor(as.numeric(as.character(treesizes.long$N)), ordered=T)
aggdata <- aggregate(formula=treesize.values ~ treesize.ind + N + sigmabm + corUbm, FUN=mean, data=treesizes.long)

# indirect greyscale plot
pdf("xy_treesizes_maineff.pdf")
xyplot(treesize.values ~ sigmabm | N + corUbm, data = aggdata, groups=treesize.ind, type="b",
       ylab="tree size", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), abline=c(1,0), 
       auto.key=list(space="top", columns=2, title="Algorithm type", cex.title=1,lines=T, points=T))
dev.off()

