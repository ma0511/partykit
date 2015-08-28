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
treespecs <- data.frame(REEMobtreesize, Mobtreesize)

load("descriptions1")
for (i in 1:length(descriptions)) {
  treespecs[i,"N"] <- substr(descriptions[[i]][[1]], 5, 8)
  treespecs[i,"rho"] <- substr(descriptions[[i]][[2]], 7, 9)
  treespecs[i,"np"] <- substr(descriptions[[i]][[3]], 24, 25)
  treespecs[i,"treatdiff"] <- substr(descriptions[[i]][[4]], 31, 33)
  treespecs[i,"corUbi"] <- substr(descriptions[[i]][[5]], 32, 64)
  treespecs[i,"numbclus"] <- substr(descriptions[[i]][[6]], 34, 35)
  treespecs[i,"sigmabi"] <- substr(descriptions[[i]][[7]], 12, 14)  
}

treespecs$N <- rep(factor(treespecs$N[1:length(descriptions)]), times=50)
treespecs$rho <- rep(factor(treespecs$rho[1:length(descriptions)]), times=50)
treespecs$np <- rep(factor(treespecs$np[1:length(descriptions)]), times=50)
treespecs$treatdiff <- rep(factor(treespecs$treatdiff[1:length(descriptions)]), times=50)
treespecs$corUbi <- rep(factor(treespecs$corUbi[1:length(descriptions)]), times=50)
treespecs$numbclus <- rep(factor(treespecs$numbclus[1:length(descriptions)]), times=50)
treespecs$sigmabi <- rep(factor(treespecs$sigmabi[1:length(descriptions)]), times=50)

table(treespecs$REEMobtreesize)
prop.table(table(treespecs$REEMobtreesize))
mean(treespecs$REEMobtreesize)
sd(treespecs$REEMobtreesize)

table(treespecs$Mobtreesize)
prop.table(table(treespecs$Mobtreesize))
mean(treespecs$Mobtreesize)
sd(treespecs$Mobtreesize)

treespecs.long <- data.frame(treesize=stack(treespecs[c("REEMobtreesize","Mobtreesize")]), 
  N=rep(treespecs$N, 2), rho=rep(treespecs$rho, 2), np=rep(treespecs$np, 2), 
  treatdiff=rep(treespecs$treatdiff, 2),  corUb=rep(treespecs$corUbi, 2), 
  numbclus=rep(treespecs$numbclus, 2), sigmab=rep(treespecs$sigmabi, 2),
  datasetID=factor(rep(1:nrow(treespecs), 2)))
tmp <- as.character(treespecs.long$treesize.ind)
tmp[treespecs.long$treesize.ind=="REEMobtreesize"] <- "GLMM tree"
tmp[treespecs.long$treesize.ind=="Mobtreesize"] <- "GLM tree"
treespecs.long$treesize.ind <- factor(tmp)

# calculate correlations of true and predicted treatment diffs (Rsquareds) 
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
  treespecs[(c-1)*length(descriptions)+(1:length(descriptions)),10:12] <- matrix(unlist(lapply(treatdiffs, cor)), ncol=9, byrow=T)[,c(2,3,6)]
  treespecs[(c-1)*length(descriptions)+(1:length(descriptions)),13:15] <- t(sapply(treatdiffs, apply, 2, mean))
  treespecs[(c-1)*length(descriptions)+(1:length(descriptions)),16:18] <- t(sapply(treatdiffs, apply, 2, sd))  
}

apply(treespecs[,10:18], 2, mean)
apply(treespecs[,10:18], 2, sd)

treespecs.long <- data.frame(correlation=stack(treespecs[c("cor_REEMob_true","cor_mob_true")]), treespecs.long)
tmp <- as.character(treespecs.long$correlation.ind)
tmp[treespecs.long$correlation.ind=="cor_REEMob_true"] <- "GLMM tree"
tmp[treespecs.long$correlation.ind=="cor_mob_true"] <- "GLM tree"
treespecs.long$correlation.ind <- factor(tmp)




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
for (j in 1:nrow(treespecs)) {
  print(j)
  tmp <- c(tmp, rep(treespecs$Mobtreesize[j], times=(treespecs$Mobtreesize[j]-1)/2))
  tmpnodeno <- c(tmpnodeno, seq(1,(treespecs$Mobtreesize[j]-1)/2))
} 
Mobsplits$treesize <- tmp 
Mobsplits$nodenumber <- tmpnodeno

tmp <- vector() 
tmpnodeno <- vector()
for (j in 1:nrow(treespecs)) {
  print(j)
  tmp <- c(tmp, rep(treespecs$REEMobtreesize[j], times=(treespecs$REEMobtreesize[j]-1)/2))
  tmpnodeno <- c(tmpnodeno, seq(1,(treespecs$REEMobtreesize[j]-1)/2))
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
mean(Mobsplits[Mobsplits$nodenumber==1 & Mobsplits$mobsplitvar==4, "mobsplitval"])
sd(Mobsplits[Mobsplits$nodenumber==1 & Mobsplits$mobsplitvar==4, "mobsplitval"])

# identify which split belongs to which tree
mobtreeno <- rep(1:32400, times=(treespecs$Mobtreesize-1)/2)
reemobtreeno <- rep(1:32400, times=(treespecs$REEMobtreesize-1)/2)
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
# get statistics for right trees 
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

treespecs$treeno <- 1:32400
treespecs$truereemobtree <- treespecs$treeno %in% unique(REEMobsplits$treeno[REEMobsplits$truetree])
treespecs$truemobtree <- treespecs$treeno %in% unique(Mobsplits$treeno[Mobsplits$truetree])

treespecs.long <- data.frame(truetree=stack(treespecs[c("truereemobtree","truemobtree")]), treespecs.long)
tmp <- as.character(treespecs.long$truetree.ind)
tmp[treespecs.long$truetree.ind=="truereemobtree"] <- "GLMM tree"
tmp[treespecs.long$truetree.ind=="truemobtree"] <- "GLM tree"
treespecs.long$truetree.ind <- factor(tmp)
treespecs.long$truetree.values <- factor(treespecs.long$truetree.values)
prop.table(table(treespecs$truemobtree))
prop.table(table(treespecs$truereemobtree))
save("treespecs.long", file="treespecs_long.dat")
save("treespecs", file="treespecs.dat")




# create lattice xyplots
library(lattice)

treesize.anova <- aov(treesize.values ~ treesize.ind + N + rho + np + treatdiff + numbclus + sigmab + 
                        corUb + treesize.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                      data=treespecs.long)
summary(treesize.anova)[[1]]["Sum Sq"] / sum(summary(treesize.anova)[[1]]["Sum Sq"] )
# N, sigmab & corUbi have main and interaction effects with eta squared > .01
treespecs.long$N <- factor(as.numeric(as.character(treespecs.long$N)), ordered=T)
treespecs.long$sigmab <- factor(as.numeric(as.character(treespecs.long$sigmab)), ordered=T)
aggdata.size <- aggregate(formula=treesize.values ~ treesize.ind + N + sigmab + corUb, FUN=mean, 
                     data=treespecs.long)
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
                            sigmab + corUb), data=treespecs.long)
summary(correlation.anova)[[1]]["Sum Sq"] / sum(summary(correlation.anova)[[1]]["Sum Sq"] )
# N, treatdiffs & sigmab have main and/or interaction effects with eta squared > .01
treespecs.long$treatdiff <- factor(as.numeric(as.character(treespecs.long$treatdiff)), ordered=T)
aggdata.cor <- aggregate(formula=correlation.values ~ correlation.ind + N + sigmab + treatdiff, FUN=mean, 
                     data=treespecs.long)
pdf("xy_correlations.pdf")
xyplot(correlation.values ~ sigmab | N + treatdiff, data = aggdata.cor, groups=correlation.ind, type="b",
       ylab="correlation", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))
dev.off()

plot(glmtree(truetree.values ~ truetree.ind | N + rho + np + treatdiff + corUb + numbclus + sigmab, 
             data=treespecs.long, maxdepth=4, family="binomial"), type="simple")
treeacc.glm <- glm(truetree.values ~ truetree.ind + N + rho + np + treatdiff + numbclus + sigmab + corUb +
                     truetree.ind*(N + rho + np + treatdiff + numbclus + sigmab + corUb), 
                   data=treespecs.long, family="binomial")
summary(treeacc.glm)
# N, corUb & sigmab have main and/or interaction effects
treespecs.long$truetree.values <- as.numeric(treespecs.long$truetree.values)-1
aggdata.acc <- aggregate(formula=truetree.values ~ truetree.ind + N + sigmab + corUb, FUN=mean, 
                     data=treespecs.long)
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="bi and splitting U correlated"] <- "b correlated with splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="bi and non-splitting U correlated"] <- "b correlated with non-splitting U"
levels(aggdata.acc$corUb)[levels(aggdata.acc$corUb)=="uncorrelated"] <- "b and U uncorrelated"
pdf("xy_accuracy.pdf", width=9)
xyplot(truetree.values ~ sigmab | N + corUb, data = aggdata.acc, groups=truetree.ind, type="b",
       ylab="tree accuracy", xlab="sigma_b", par.settings=standard.theme("pdf",color=F), 
       auto.key=list(space="top", columns=2, title=" ", cex.title=1,lines=T, points=T))
dev.off()
