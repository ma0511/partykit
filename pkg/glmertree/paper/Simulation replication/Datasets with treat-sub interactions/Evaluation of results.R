library(lme4)
library(partykit)

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

load("testdescriptions1")
for (i in 1:length(testdescriptions)) {
  treesizes[i,"N"] <- substr(testdescriptions[[i]][[1]], 5, 8) # N
  treesizes[i,"rho"] <- substr(testdescriptions[[i]][[2]], 7, 9) # is rho (correlation between covariates)
  treesizes[i,"np"] <- substr(testdescriptions[[i]][[3]], 24, 25) # is number of covariates
  treesizes[i,"treatdiff"] <- substr(testdescriptions[[i]][[4]], 31, 33) # is treatment effect difference
  treesizes[i,"corUbi"] <- substr(testdescriptions[[i]][[5]], 32, 64) # is correlation between U and bi
  treesizes[i,"numbclus"] <- substr(testdescriptions[[i]][[6]], 34, 35) # is number of random intercept values
  treesizes[i,"sigmabi"] <- substr(testdescriptions[[i]][[7]], 12, 14) # max and -min random intercept value    
}

treesizes$N <- factor(treesizes$N)
treesizes$rho <-  factor(treesizes$rho)            
treesizes$np <- factor(treesizes$np)              
treesizes$treatdiff <- factor(treesizes$treatdiff)          
treesizes$corUbi <- factor(treesizes$corUbi)           
treesizes$numbclus <- factor(treesizes$numbclus)        
treesizes$sigmabi <- factor(treesizes$sigmabi)

treesizes$N <- rep(treesizes$N[1:length(testdescriptions)], times=1)
treesizes$rho <- rep(treesizes$rho[1:length(testdescriptions)], times=1)
treesizes$np <- rep(treesizes$np[1:length(testdescriptions)], times=1)
treesizes$treatdiff <- rep(treesizes$treatdiff[1:length(testdescriptions)], times=1)
treesizes$corUbi <- rep(treesizes$corUbi[1:length(testdescriptions)], times=1)
treesizes$numbclus <- rep(treesizes$numbclus[1:length(testdescriptions)], times=1)
treesizes$sigmabi <- rep(treesizes$sigmabi[1:length(testdescriptions)], times=1)

save(treesizes,file="treesizes")
load("treesizes")

table(treesizes[,1])
table(treesizes[,2])
prop.table(table(treesizes[,1]))
prop.table(table(treesizes[,2]))

mean(treesizes[,1]);sd(treesizes[,1]) # REEMob
mean(treesizes[,2]);sd(treesizes[,2]) # mob

treesizes.long <- data.frame(treesize=stack(treesizes[,1:2]), 
  N=rep(treesizes$N, 2), rho=rep(treesizes$rho, 2), np=rep(treesizes$np, 2), 
  treatdiff=rep(treesizes$treatdiff, 2),  corUb=rep(treesizes$corUbi, 2), 
  numbclus=rep(treesizes$numbclus, 2), sigmab=rep(treesizes$sigmabi, 2),
  datasetID=factor(rep(1:nrow(treesizes), 2)))
tmp <- as.character(treesizes.long$treesize.ind)
tmp[treesizes.long$treesize.ind=="REEMobtreesize"] <- "GLMM tree"
tmp[treesizes.long$treesize.ind=="Mobtreesize"] <- "GLM tree"
treesizes.long$treesize.ind <- factor(tmp)

# calculate correlations of true and predicted treatment diffs (Rsquareds) 
treespecs <- treesizes
treespecs$cor_REEMob_true <- NA
treespecs$cor_mob_true <- NA
treespecs$cor_mob_REEMob <- NA
treespecs$mean_true_diffs <- NA
treespecs$mean_REEMob_diffs <- NA
treespecs$mean_mob_diffs <- NA
treespecs$sd_true_diffs <- NA
treespecs$sd_REEMob_diffs <- NA
treespecs$sd_mob_diffs <- NA

for(c in 1:50) {
  print(c)
  load(paste("treatdiffs", c, sep=""))
  treespecs[(c-1)*length(testdescriptions)+(1:length(testdescriptions)),10:12] <- matrix(unlist(lapply(treatdiffs, cor)), ncol=9, byrow=T)[,c(2,3,6)]
  treespecs[(c-1)*length(testdescriptions)+(1:length(testdescriptions)),13:15] <- t(sapply(treatdiffs, apply, 2, mean))
  treespecs[(c-1)*length(testdescriptions)+(1:length(testdescriptions)),16:18] <- t(sapply(treatdiffs, apply, 2, sd))  
}

apply(treespecs[,10:18], 2, mean)
apply(treespecs[,10:18], 2, sd)

treatdiffcors.long <- data.frame(correlation=stack(treespecs[,c(10,11)]), 
                             N=rep(treesizes$N, 2), rho=rep(treesizes$rho, 2), np=rep(treesizes$np, 2), 
                             treatdiff=rep(treesizes$treatdiff, 2),  corUb=rep(treesizes$corUbi, 2), 
                             numbclus=rep(treesizes$numbclus, 2), sigmab=rep(treesizes$sigmabi, 2),
                             datasetID=factor(rep(1:nrow(treesizes), 2)))
tmp <- as.character(treatdiffcors.long$correlation.ind)
tmp[treatdiffcors.long$correlation.ind=="cor_REEMob_true"] <- "GLMM tree"
tmp[treatdiffcors.long$correlation.ind=="cor_mob_true"] <- "GLM tree"
treatdiffcors.long$correlation.ind <- factor(tmp)




# check splitting variables and split points

tmp <- data.frame()
for (c in 1:50){
  load(paste("splits.Mobs",c,sep=""))
  tmp <- data.frame(rbind(tmp, splits.Mobs))
}
Mobsplits <- data.frame(mobsplitvar=tmp[,1], mobsplitval=tmp[,2])

tmp <- data.frame()
for (c in 1:50){
  load(paste("splits.REEMobs",c,sep=""))
  tmp <- data.frame(rbind(tmp, splits.REEMobs))
}
REEMobsplits <- data.frame(REEMobsplitvar=tmp[,1], REEMobsplitval=tmp[,2])

tmp <- vector() 
tmpnodeno <- vector()
for (j in 1:nrow(treesizes)) {
  print(j)
  tmp <- c(tmp, rep(treesizes$Mobtreesize[j], times=(treesizes$Mobtreesize[j]-1)/2))
  tmpnodeno <- c(tmpnodeno, seq(1,(treesizes$Mobtreesize[j]-1)/2))
} 
Mobsplits$treesize <- tmp 
Mobsplits$nodenumber <- tmpnodeno

tmp <- vector() 
tmpnodeno <- vector()
for (j in 1:nrow(treesizes)) {
  print(j)
  tmp <- c(tmp, rep(treesizes$REEMobtreesize[j], times=(treesizes$REEMobtreesize[j]-1)/2))
  tmpnodeno <- c(tmpnodeno, seq(1,(treesizes$REEMobtreesize[j]-1)/2))
} 
REEMobsplits$treesize <- tmp
REEMobsplits$nodenumber <- tmpnodeno

REEMobsplits$REEMobsplitvar <- factor(REEMobsplits$REEMobsplitvar)
Mobsplits$mobsplitvar <- factor(Mobsplits$mobsplitvar)
save(REEMobsplits, file="REEMobsplits");save(Mobsplits, file="Mobsplits")
load("REEMobsplits"); load("Mobsplits")

# Assess accuracy of first split
table(REEMobsplits[REEMobsplits$nodenumber==1, "REEMobsplitvar"])
mean(REEMobsplits[REEMobsplits$nodenumber==1, "REEMobsplitval"])
sd(REEMobsplits[REEMobsplits$nodenumber==1, "REEMobsplitval"])
table(Mobsplits[Mobsplits$nodenumber==1, "mobsplitvar"])
mean(Mobsplits[Mobsplits$nodenumber==1, "mobsplitval"])
sd(Mobsplits[Mobsplits$nodenumber==1, "mobsplitval"])

# identify which split belongs to which tree
mobtreeno <- rep(1:32400, times=(treesizes$Mobtreesize-1)/2)
reemobtreeno <- rep(1:32400, times=(treesizes$REEMobtreesize-1)/2)
REEMobsplits$treeno <- reemobtreeno
Mobsplits$treeno <- mobtreeno

# Find the true trees (that is, true value plus or minus 5 = .5SD)
REEMobtruefirstsplit <- REEMobsplits[REEMobsplits$treesize==7 & REEMobsplits$REEMobsplitvar==4 & 
                                       REEMobsplits$REEMobsplitval>25 & REEMobsplits$REEMobsplitval<35,"treeno"]
REEMobtruesecsplit <- REEMobsplits[REEMobsplits$treesize==7 & REEMobsplits$REEMobsplitvar==3 & 
                                     REEMobsplits$REEMobsplitval>12 & REEMobsplits$REEMobsplitval<23,"treeno"]
REEMobtruethirdsplit <- REEMobsplits[REEMobsplits$treesize==7 & REEMobsplits$REEMobsplitvar==7 & 
                                       REEMobsplits$REEMobsplitval>58 & REEMobsplits$REEMobsplitval<68,"treeno"]
REEMobtruetreenos <- REEMobtruefirstsplit[REEMobtruefirstsplit %in% REEMobtruesecsplit][
  REEMobtruefirstsplit[REEMobtruefirstsplit %in% REEMobtruesecsplit] %in% REEMobtruethirdsplit]
REEMobsplits$truetree <- REEMobsplits$treeno %in% REEMobtruetreenos
table(REEMobsplits$truetree)
# get statistics for right trees 
names(REEMobsplits)
tapply(REEMobsplits[REEMobsplits$truetree,"REEMobsplitval"],
       REEMobsplits[REEMobsplits$truetree,"REEMobsplitvar"], length)
tapply(REEMobsplits[REEMobsplits$truetree,"REEMobsplitval"],
       REEMobsplits[REEMobsplits$truetree,"REEMobsplitvar"], mean)
tapply(REEMobsplits[REEMobsplits$truetree,"REEMobsplitval"],
       REEMobsplits[REEMobsplits$truetree,"REEMobsplitvar"], sd)
# get statistics for wrong trees
tapply(REEMobsplits[!REEMobsplits$truetree,"REEMobsplitval"],
       REEMobsplits[!REEMobsplits$truetree,"REEMobsplitvar"], length)

Mobtruefirstsplit <- Mobsplits[Mobsplits$treesize==7 & Mobsplits$mobsplitvar==4 & 
                                 Mobsplits$mobsplitval>25 & Mobsplits$mobsplitval<35,"treeno"]
Mobtruesecsplit <- Mobsplits[Mobsplits$treesize==7 & Mobsplits$mobsplitvar==3 & 
                               Mobsplits$mobsplitval>12 & Mobsplits$mobsplitval<22,"treeno"]
Mobtruethirdsplit <- Mobsplits[Mobsplits$treesize==7 & Mobsplits$mobsplitvar==7 & 
                                 Mobsplits$mobsplitval>58 & Mobsplits$mobsplitval<68,"treeno"]
Mobtruetreenos <- Mobtruefirstsplit[Mobtruefirstsplit %in% Mobtruesecsplit][
  Mobtruefirstsplit[Mobtruefirstsplit %in% Mobtruesecsplit] %in% Mobtruethirdsplit]
Mobsplits$truetree <- Mobsplits$treeno %in% Mobtruetreenos
table(Mobsplits$truetree)
# get statistics for right trees
tapply(Mobsplits[Mobsplits$truetree,"mobsplitval"],
       Mobsplits[Mobsplits$truetree,"mobsplitvar"], length)
tapply(Mobsplits[Mobsplits$truetree,"mobsplitval"],
       Mobsplits[Mobsplits$truetree,"mobsplitvar"], mean)
tapply(Mobsplits[Mobsplits$truetree,"mobsplitval"],
       Mobsplits[Mobsplits$truetree,"mobsplitvar"], sd)
# get statistics for wrong trees
tapply(Mobsplits[!Mobsplits$truetree,"mobsplitval"],
       Mobsplits[!Mobsplits$truetree,"mobsplitvar"], length)

treesizes$treeno <- 1:32400
treesizes$truemobtree <- treesizes$treeno %in% unique(Mobsplits$treeno[Mobsplits$truetree])
treesizes$truereemobtree <- treesizes$treeno %in% unique(REEMobsplits$treeno[REEMobsplits$truetree])

treeacc.long <- data.frame(treesize=stack(treesizes[,1:2]), 
                             N=rep(treesizes$N, 2), rho=rep(treesizes$rho, 2), np=rep(treesizes$np, 2), 
                             treatdiff=rep(treesizes$treatdiff, 2),  corUb=rep(treesizes$corUbi, 2), 
                             numbclus=rep(treesizes$numbclus, 2), sigmab=rep(treesizes$sigmabi, 2),
                             datasetID=factor(rep(1:nrow(treesizes), 2)), truetree=stack(treesizes[,11:12]))
tmp <- as.character(treeacc.long$treesize.ind)
tmp[treeacc.long$truetree.ind=="truereemobtree"] <- "GLMM tree"
tmp[treeacc.long$truetree.ind=="truemobtree"] <- "GLM tree"
treeacc.long$truetree.ind <- factor(tmp)
treeacc.long$truetree.values <- factor(treeacc.long$truetree.values)
prop.table(table(treesizes$truemobtree))
prop.table(table(treesizes$truereemobtree))





# create lattice xyplots
library(lattice)

treesize.anova <- aov(treesize.values ~ treesize.ind + N + rho + np + treatdiff + numbclus + sigmab + 
                        corUb + treesize.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                      data=treesizes.long)
summary(treesize.anova)[[1]]["Sum Sq"] / sum(summary(treesize.anova)[[1]]["Sum Sq"] )
# N, sigmab & corUbi have main and interaction effects with eta squared > .01
treesizes.long$N <- factor(as.numeric(as.character(treesizes.long$N)), ordered=T)
treesizes.long$sigmab <- factor(as.numeric(as.character(treesizes.long$sigmab)), ordered=T)
aggdata.size <- aggregate(formula=treesize.values ~ treesize.ind + N + sigmab + corUb, FUN=mean, 
                     data=treesizes.long)
levels(aggdata.size$corUb)[levels(aggdata.size$corUb)=="bi and splitting U correlated"] <- "b correlated with splitting U"
levels(aggdata.size$corUb)[levels(aggdata.size$corUb)=="bi and non-splitting U correlated"] <- "b correlated with non-splitting U"
levels(aggdata.size$corUb)[levels(aggdata.size$corUb)=="uncorrelated"] <- "b and U uncorrelated"
pdf("xy_treesizes_treatsubs.pdf", width=9)
xyplot(treesize.values ~ sigmab | N + corUb, data = aggdata.size, groups=treesize.ind, type="b", 
       ylab="tree size", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), abline=c(7,0),
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))
dev.off()

correlation.anova <- aov(correlation.values ~ correlation.ind + N + rho + np + treatdiff + numbclus + 
                          sigmab + corUb + correlation.ind*(N + rho + np + treatdiff + numbclus + 
                            sigmab + corUb), data=treatdiffcors.long)
summary(correlation.anova)[[1]]["Sum Sq"] / sum(summary(correlation.anova)[[1]]["Sum Sq"] )
# N, treatdiffs & sigmab have main and/or interaction effects with eta squared > .01
treatdiffcors.long$N <- factor(as.numeric(as.character(treatdiffcors.long$N)), ordered=T)
treatdiffcors.long$sigmab <- factor(as.numeric(as.character(treatdiffcors.long$sigmab)), ordered=T)
treatdiffcors.long$treatdiff <- factor(as.numeric(as.character(treatdiffcors.long$treatdiff)), ordered=T)
aggdata.cor <- aggregate(formula=correlation.values ~ correlation.ind + N + sigmab + treatdiff, FUN=mean, 
                     data=treatdiffcors.long)
pdf("xy_correlations.pdf")
xyplot(correlation.values ~ sigmab | N + treatdiff, data = aggdata.cor, groups=correlation.ind, type="b",
       ylab="correlation", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))
dev.off()

plot(glmtree(truetree.values ~ truetree.ind | N + rho + np + treatdiff + corUb + numbclus + sigmab, 
             data=treeacc.long, maxdepth=4, family="binomial"), type="simple")
treeacc.glm <- glm(truetree.values ~ truetree.ind + N + rho + np + treatdiff + numbclus + sigmab + corUb +
                     truetree.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                   data=treeacc.long, family="binomial")
summary(treeacc.glm)
# N, corUb & sigmab have main and/or interaction effects
treeacc.long$truetree.values <- as.numeric(treeacc.long$truetree.values)-1
treeacc.long$N <- factor(as.numeric(as.character(treeacc.long$N)), ordered=T)
treeacc.long$sigmab <- factor(as.numeric(as.character(treeacc.long$sigmab)), ordered=T)
aggdata.acc <- aggregate(formula=truetree.values ~ truetree.ind + N + sigmab + corUb, FUN=mean, 
                     data=treeacc.long)
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="bi and splitting U correlated"] <- "b correlated with splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="bi and non-splitting U correlated"] <- "b correlated with non-splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="uncorrelated"] <- "b and U uncorrelated"
pdf("xy_accuracy.pdf", width=9)
xyplot(truetree.values ~ sigmab | N + corUb, data = aggdata.acc, groups=truetree.ind, type="b",
       ylab="tree accuracy", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))
dev.off()
