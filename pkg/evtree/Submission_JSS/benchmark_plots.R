###################################################
### Benchmark plots for Section: 5. Performance comparisons 
###################################################

## load results
rm(list = ls())
for(i in Sys.glob("results/*.RData")) load(i)
for(i in Sys.glob("results_j48/*.RData")) load(i)

## preprocess for reference evtree
preprocess <- function(d, dname = "datasetname", isclassification = TRUE){
    if(isclassification){
        colnames(d) <- c("evtree", "rpart", "ctree", "J48","evtree", "rpart", "ctree", "J48")
		d[, 1:4] <- 1 - d[ ,1:4]
    }else{
    	colnames(d) <- c("evtree", "rpart", "ctree","evtree", "rpart", "ctree")    	
   	}
    d <- as.data.frame(d)
	nAlgorithms = dim(d)[2]/2
    for(i in nAlgorithms:1) d[, i] <- d[, i] / d[, 1] * 100
    if(isclassification)  # for J48 the total number of nodes is used
	    d[, nAlgorithms*2] <- d[, nAlgorithms*2] / (d[, nAlgorithms+1]*2+1) * 100
	else
		d[, nAlgorithms*2] <- d[, nAlgorithms*2] / d[, nAlgorithms+1] * 100
    for(i in (nAlgorithms*2-1):(nAlgorithms+1)) d[, i] <- d[, i] / d[, nAlgorithms+1] * 100
    
    x <- d[, 1:nAlgorithms]
    y <- d[, (nAlgorithms+1):(nAlgorithms*2)]
    rval <- reshape(x, idvar="samp", times=names(x), timevar = "alg",varying= list(names(x)), direction="long")
    names(rval)[2] <- "accuracy"
    rval$complexity <- reshape(y, idvar="samp", times=names(y), timevar = "alg",varying= list(names(y)), direction="long")[,2]
    if(isclassification) 
    	rval$alg <- factor(rval$alg, levels = c("evtree", "ctree", "rpart", "J48"))
    else 
    	rval$alg <- factor(rval$alg, levels = c("evtree", "ctree", "rpart"))
    rval$ds <- dname
    rval
}

## collect results for all datasets
r <- rbind(
preprocess(d = cbind(rglass[,1:3], rglass2[,3], rglass[,4:6], rglass2[,4]), dname = "Glass identification", isclassification = TRUE),
preprocess(d = cbind(rheart[,1:3], rheart2[,3], rheart[,4:6], rheart2[,4]), dname = "Statlog heart", isclassification = TRUE),
preprocess(d = cbind(rionosphere[,1:3], rionosphere2[,3], rionosphere[,4:6], rionosphere2[,4]), dname = "Ionosphere", isclassification = TRUE),
preprocess(d = cbind(rmusk[,1:3], rmusk2[,3], rmusk[,4:6], rmusk2[,4]), dname = "Musk", isclassification = TRUE),
preprocess(d = cbind(rbreastcancer[,1:3], rbreastcancer2[,3], rbreastcancer[,4:6], rbreastcancer2[,4]), dname = "Breast cancer database", isclassification = TRUE),
preprocess(d = cbind(rpima[,1:3], rpima2[,3], rpima[,4:6], rpima2[,4]), dname = "Pima Indians diabetes", isclassification = TRUE),
preprocess(d = cbind(rvowel[,1:3], rvowel2[,3], rvowel[,4:6], rvowel2[,4]), dname = "Vowel", isclassification = TRUE),
preprocess(d = cbind(rcredit[,1:3], rcredit2[,3], rcredit[,4:6], rcredit2[,4]), dname = "Statlog German credit", isclassification = TRUE),
preprocess(d = cbind(rcontraceptive[,1:3], rcontraceptive2[,3], rcontraceptive[,4:6], rcontraceptive2[,4]), dname = "Contraceptive method", isclassification = TRUE),
preprocess(d = cbind(rdna[,1:3], rdna2[,3], rdna[,4:6], rdna2[,4]), dname = "DNA", isclassification = TRUE),
preprocess(d = cbind(rspam[,1:3], rspam2[,3], rspam[,4:6], rspam2[,4]), dname = "Spam", isclassification = TRUE),
preprocess(d = cbind(rmagicgamma[,1:3], rmagicgamma2[,3], rmagicgamma[,4:6], rmagicgamma2[,4]), dname = "Magic gamma telescope", isclassification = TRUE),
preprocess(d = rservo, dname = "Servo", isclassification = FALSE),
preprocess(d = rbostonhousing, dname = "Boston housing", isclassification = FALSE),
preprocess(d = rmel0101, dname = "MEL0101", isclassification = FALSE),
preprocess(d = rhdg0202, dname = "HDG0202", isclassification = FALSE),
preprocess(d = rhdg0502, dname = "HDG0502", isclassification = FALSE)
)

r$ds <- factor(r$ds)
r$samp <- factor(r$samp)
r$dssamp <- r$ds:r$samp

## compute multiple comparisons
library("multcomp")
cstats <- function(alg = "evtree", value = "accuracy", data = r) {
  dlab <- rev(unique(data$ds))
  if(alg == "J48"){ 
  	dlab <- dlab[-c(1:5)] ## J48: skip regression datasets
  }
  k <- length(dlab)  
  mean  <- numeric(k)
  lower <- numeric(k)
  upper <- numeric(k)
  names(data)[names(data) == value] <- "value"
  firstDS <- 1
  for(i in 1:k) {
  	dsub <- subset(data, ds == dlab[i])
	dsub$alg <- factor(dsub$alg)
    mod1 <- lm(value ~ alg, data = dsub)
    pt <- glht(mod1, linfct = mcp(alg = "Dunnett"))
    w <- confint(pt)$confint
    d <- which(levels(dsub$alg) == alg) - 1
    mean[i]  <-  w[d]
    lower[i] <-  w[d + length(levels(dsub$alg))-1]
    upper[i] <-  w[d + (length(levels(dsub$alg))-1)*2]
  }
  rval <- data.frame(mean, lower, upper)
  rownames(rval) <- dlab
  return(rval)
}

acc_rpart <- cstats("rpart", "accuracy")
com_rpart <- cstats("rpart", "complexity")
acc_ctree <- cstats("ctree", "accuracy")
com_ctree <- cstats("ctree", "complexity")
acc_J48 <- cstats("J48", "accuracy") 
com_J48 <- cstats("J48", "complexity") 

## function for visualization
ciplot <- function(x, xlim = NULL, main = "", xlab = "", ylab = TRUE) {
  nam <- rownames(x)
  k <- length(nam)
  plot(x$mean, 1:k, xlim = xlim, axes = FALSE, xlab = "", ylab = "", pch = 19)
  arrows(x$lower, 1:k, x$upper, 1:k, angle = 90, length = 0.05, code = 3)
  if(xlab == "") axis(1, labels = FALSE) else axis(1)
  if(ylab) ylab <- nam
  axis(2, at = 1:k, labels = ylab, las = 1, cex = 0.8)  
  axis(2, at = k + 1.5, labels = main, tick = FALSE, las = 1, outer = TRUE, cex.axis = 1.5, xpd = TRUE)
  mtext(xlab, side = 1, line = 3, xpd = TRUE)
  if (dim(x) >= 17) abline(h = 5.5)
  abline(v = 0, lty = 2)  
  box()
}

## plot the results if evtree vs. rpart and evtree vs. ctree
par(mfrow = c(2, 2), oma = c(5, 10, 2, 0), mar = c(1, 1, 2, 1))

xlim1 <- range(cbind(acc_rpart, acc_ctree))
xlim2 <- range(cbind(com_rpart, com_ctree))

ciplot(acc_rpart, xlim = xlim1, main = "rpart", ylab = TRUE, xlab = "")
ciplot(com_rpart, xlim = xlim2, main = "",      ylab = FALSE, xlab = "")
ciplot(acc_ctree, xlim = xlim1, main = "ctree", ylab = TRUE,
  xlab = "relative difference in predictive accuracy (%)")
ciplot(com_ctree, xlim = xlim2, main = "",      ylab = FALSE,
  xlab = "relative difference in complexity (%)")

## plot the results of evtree vs. J48
par(mfrow = c(1, 2), oma = c(5, 10, 2, 0), mar = c(1, 1, 2, 1))

xlim1 <- range(acc_J48)
xlim2 <- range(com_J48)

ciplot(acc_J48, xlim = xlim1, main = "J48", ylab = TRUE,
  xlab = "relative difference in predictive accuracy (%)")
ciplot(com_J48, xlim = xlim2, main = "",      ylab = FALSE,
  xlab = "relative difference in complexity (%)")



###################################################
### Benchmark plots for Section: 6. Parameter optimization 
###################################################
library("plotrix")
library("lme4")
library("multcomp")

## load results
for(i in Sys.glob("results_parameter/*.RData")) load(i)

panel.mean <- function(x,y,...){
	x <- as.numeric(x)
	x.unique <- unique(x)
	for(X in x.unique) {
		Y <- y[x == X]
		if (!length(Y)) next
		mean.value <- list(y = mean(Y), x = X)
		do.call("lpoints", c(mean.value, pch = 20, col= "red"))
	}
}

# functions for the visualisation of different operatorprobabilities 

preprocess_op <- function(d, dname = "datasetname"){
	d <- as.data.frame(d)
    colnames(d) <- c("c0m50sp50", "c20m40sp40","c40m30sp30","c20m20sp60","c20m60sp20")
    x <- d*100
    x <- reshape(x, idvar="samp", times=names(x), timevar = "operatorprob",varying= list(names(x)), direction="long")
	names(x)[[2]] <- "value"
	x$ds <- dname
    return(x)
}

preprocess_op_comp <- function(ntrees, colum){	
	rt <- preprocess_op(d = rheart[,colum], dname = "Statlog heart")
	rt <- rbind(rt, preprocess_op(d = rcredit[,colum], dname = "Statlog German credit"))
	rt <- rbind(rt, preprocess_op(d = rspam[,colum], dname = "Spam"))
	rt <- rbind(rt, preprocess_op(d = rchessboard44_5[,colum], dname = "Chessboard 4x4 with 5% noise"))
	rt <- cbind(rt, rep( ntrees, dim(rt)[1]))
	colnames(rt)[5] <- "nIter"
	rt
}

r2 <- preprocess_op_comp("200 iterations", c(11:15))
r2 <- rbind(r2, preprocess_op_comp("500 iterations", c(16:20)))
r2 <- rbind(r2, preprocess_op_comp("10000 iterations", c(21:25)))


sort_op <- function(x){
	x$ds <- relevel(x$ds, "Statlog heart")
	x$ds <- relevel(x$ds, "Statlog German credit")
	x$ds <- relevel(x$ds, "Spam")
	x$ds <- relevel(x$ds, "Chessboard 4x4 with 5% noise")
	x
}

r2$operatorprob <- factor(r2$operatorprob)
r2$ds <- factor(r2$ds)
r2$nIter <- factor(r2$nIter)
r2 <- sort_op(r2)

b1 <- bwplot (value ~ factor(operatorprob)| nIter+ds , data= as.data.frame(r2), 
horizontal = FALSE,  
ylab= list("Accuracy (%)", cex=1.1),
pch= '|',
layout = c(3,4),
ylim=as.data.frame(matrix(c(
rep(c(55,100),3),
rep(c(85,94),3),
rep(c(60,80),3),
rep(c(58,90),3)
), nrow=2)),
scales= list(x= list(rot=60), y=list(relation="free"), alternating = F), 
panel=function(x,y,...) {
	panel.bwplot(x,y,...)
	panel.mean(x,y,...)
}
)
b1


# functions for the visualisation of different population sizes 
preprocess_ntrees <- function(d, dname = "datasetname"){
	d <- as.data.frame(d)
    colnames(d) <- c("25 trees", "50 trees","100 trees","250 trees","500 trees")
    x <- d*100
    x <- reshape(x, idvar="samp", times=names(x), timevar = "ntrees",varying= list(names(x)), direction="long")
    names(x)[[2]] <- "value"
	x$ds <- dname
    return(x)
}

r <- preprocess_ntrees(d = rheart[,1:5], dname = "Statlog heart")
r <- rbind(r, preprocess_ntrees(d = rcredit[,1:5], dname = "Statlog German credit"))
r <- rbind(r, preprocess_ntrees(d = rspam[,1:5], dname = "Spam"))
r <- rbind(r, preprocess_ntrees(d = rchessboard44_5[,1:5], dname = "Chessboard 4x4 with 5% noise"))

sort_ntrees <- function(x){
	x$ds <- relevel(x$ds, "Chessboard 4x4 with 5% noise")
	x$ds <- relevel(x$ds, "Spam")
	x$ds <- relevel(x$ds, "Statlog German credit")
    x$ds <- relevel(x$ds, "Statlog heart")
	x$ntrees <- relevel(x$ntrees, "250 trees")
	x$ntrees <- relevel(x$ntrees, "100 trees")
	x$ntrees <- relevel(x$ntrees, "50 trees")
	x$ntrees <- relevel(x$ntrees, "25 trees")
	x
}

r$ntrees <- factor(r$ntrees) 
r$ds <- factor(r$ds)
r <- sort_ntrees(r) 
par.settings = list(cex=1.2)
b1 <- bwplot (value ~ ntrees | ds, data= as.data.frame(r), 
horizontal = FALSE,  
ylab= list("Accuracy (%)", cex=1.1),
pch= '|',
ylim=as.data.frame(matrix(c(
c(60,90),
c(65,80),
c(89,94),
c(60,95)
), nrow=2)),
scales= list(x= list(rot=60), y=list(relation="free"), alternating = F), 
layout = c(4,1),
panel=function(x,y,...) {
	panel.bwplot(x,y,...)
	panel.mean(x,y,...)
}
)
b1














