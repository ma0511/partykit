###################################################
### This file reproduces the benchmarks experiments and plots of the evtree paper.
### Note that the full analysis takes rather long (our analysis on a mainframe with 
### 6 parallel processes took 1-2 weeks). 
###################################################

## required libraries
library("rpart")
library("party")
library("partykit")
library("kernlab")
library("mlbench")
library("evtree")
library("RWeka")

###################################################
### Calculates results of J48 for the benchmark experiments.
### Input: A benchmark dataset (the dependent variable has to be in the last row).
### Output: Results of the benchmark.
###################################################
bootstrap_j48 <- function(kdata, nboots = 250, seed = 1000){
    set.seed(seed)
    accuracy_j48 <- array(-999999, nboots)
    accuracy_j48_up <- array(-999999, nboots)
    NN_j48 <- array(-999999, nboots)
    NN_j48_up <- array(-999999, nboots)

    print("bootstrap no.:")
    for(f in 1:nboots){
        print(f)
        rand <- sample(nrow(kdata), replace = T)
        rand <- sort(rand)
        nbootTraining <- kdata[rand, ]
        s <- 1
        k <- 1
        flag <- FALSE
        while(k <= dim(kdata)[1]){
            if(rand[k] == s ){
                s <- s + 1
                k <- k + 1
            }else if(rand[k] < s ){
                k <- k + 1
            }else if(rand[k] > s ){
                if(flag == FALSE){
                  nbootTest <- kdata[s,]
                  flag <- TRUE
                }else{
                  nbootTest <- rbind(nbootTest, kdata[s,])
                }
                s <- s + 1
            }
        }

        xtrain <- nbootTraining[ ,1:(dim(nbootTraining)[2]-1)]
        ytrain <- nbootTraining[ ,dim(nbootTraining)[2]]
        xpredict <- nbootTest[ ,1:(dim(nbootTest)[2]-1)]
        ypredict <- nbootTest[ ,dim(nbootTest)[2]]
		
		out_j48 <- J48(
                ytrain ~ ., data = cbind(xtrain, ytrain),
                control = Weka_control(
                	R = TRUE,
					M = 7
				)
		)
		out_j48_up <- J48(
                ytrain ~ ., data = cbind(xtrain, ytrain),
                control = Weka_control(
					M = 7
				)
		)

		fit_j48 <- predict(out_j48, newdata= xpredict)
		fit_j48_up <- predict(out_j48_up, newdata= xpredict)
        accuracy_j48[f] <- sum(fit_j48 == ypredict) / length(ypredict)
        accuracy_j48_up[f] <- sum(fit_j48_up == ypredict) / length(ypredict)
	    NN_j48[f] <- length(as.party(out_j48))
	    NN_j48_up[f] <- length(as.party(out_j48_up))
		
    } # end for
    accuracy <- cbind(accuracy_j48, NN_j48, accuracy_j48_up, NN_j48)
    accuracy
}


###################################################
### Calculate the results for the chessboard problems.
### Input: number of points, number of noise variable, noise level, number of simulations and seed.
### Output: results of the simulation.
###################################################
chessboard44 <- function(n = 4000, noisevariables = 6, noise = 0) {
  chess44 <- array(0,c(n,noisevariables+3))
  for(i in 1:(noisevariables+2))
      chess44[,i] <- as.numeric(runif(dim(chess44)[1]))*4

   x <- chess44[,1]
   y <- chess44[,2]
   chess44[, ncol(chess44)] <- 0
   for(k in 1:4)  
      chess44[(x <= k & x > k-1 & y <= k & y > k-1), ncol(chess44)] <- 1
   for(k in 1:2)  
      chess44[(x <= k & x > k-1 & y <= k+2 & y > k+1), ncol(chess44)] <- 1
   for(k in 1:2)  
      chess44[(y <= k & y > k-1 & x <= k+2 & x > k+1), ncol(chess44)] <- 1

   if(noise > 0) {
      flipclasslist <- sample(n, n * (noise / 100), replace = FALSE)

      for(i in 1:length(flipclasslist)){
	  if(chess44[flipclasslist[i], ncol(chess44)] == 1)
	      chess44[flipclasslist[i], ncol(chess44)] = 0
	  else if(chess44[flipclasslist[i], ncol(chess44)] == 0)
	      chess44[flipclasslist[i], ncol(chess44)] = 1
      }
  }

  chess44 <- as.data.frame(chess44)
  chess44[,ncol(chess44)] <- as.factor(chess44[,ncol(chess44)])
  names(chess44) <- c(paste("X", 1:8, sep = ""), "Y")
  chess44
}

benchchessboard44_j48 <- function(n = 4000, noisevariables = 6, noise = 0, nrealizations = 250, seed = 1000){
    set.seed(seed)
    accuracy_j48 <- array(-999999, nboots)
    accuracy_j48_up <- array(-999999, nboots)
    NN_j48 <- array(-999999, nboots)
    NN_j48_up <- array(-999999, nboots)


    print("bootstrap no.:")
    for(f in 1:nrealizations){
        print(f)

        ch <- chessboard44(n = n, noisevariables = noisevariables, noise = noise )
        ch_train <- ch[1:(n/2), ]
        ch_test <-  ch[(n/2+1):n, ]
        xtrain <-   ch_train[, 1:(ncol(ch_train)-1)]
        ytrain <-   ch_train[, ncol(ch_train)]
        xpredict <- ch_test[, 1:(ncol(ch_test)-1)]
        ypredict <- ch_test[, ncol(ch_test)]

		out_j48 <- J48(
                ytrain ~ ., data = cbind(xtrain, ytrain),
                control = Weka_control(
                	R = TRUE,
					M = 7
				)
		)
		out_j48_up <- J48(
                ytrain ~ ., data = cbind(xtrain, ytrain),
                control = Weka_control(
					M = 7
				)
		)

		fit_j48 <- predict(out_j48, newdata= xpredict)
		fit_j48_up <- predict(out_j48_up, newdata= xpredict)
        accuracy_j48[f] <- sum(fit_j48 == ypredict) / length(ypredict)
        accuracy_j48_up[f] <- sum(fit_j48_up == ypredict) / length(ypredict)
	    NN_j48[f] <- length(as.party(out_j48))
	    NN_j48_up[f] <- length(as.party(out_j48_up))
		
    } # end for
    accuracy <- cbind(accuracy_j48, NN_j48, accuracy_j48_up, NN_j48)
    accuracy
}

###################################################
### Processing of datasets
###################################################

### 1. Glass ###
data(Glass)
rglass <- bootstrap_j48(Glass, nboots=250)
apply(apply(rglass, 2, as.numeric),2,mean)
save(rglass, file="results2/glass_j48.RData")

### 2. statlogHeart ###
data(StatlogHeart)
rheart <- bootstrap_j48(StatlogHeart, nboots=250)
apply(apply(rheart, 2, as.numeric),2,mean)
save(rheart, file="results2/heart_j48.RData")

### 3. Ionosphere ###
data(Ionosphere)
rionosphere <- bootstrap_j48(Ionosphere, nboots=250)
apply(apply(rionosphere, 2, as.numeric),2,mean)
save(rionosphere, file="results2/ionosphere_j48.RData")

### 4. musk ###
data(musk)
rmusk <- bootstrap_j48(musk, nboots=250)
apply(apply(rmusk, 2, as.numeric),2,mean)
save(rmusk, file="results2/musk.RData_j48")

### 5. BreastCancer ###
data(BreastCancer)
BreastCancer <- BreastCancer[complete.cases(BreastCancer),2:11]
for(i in 1:5)
	BreastCancer[,i] <- as.numeric(BreastCancer[,i]) 
for(i in 6:10)
	BreastCancer[,i] <- as.factor(BreastCancer[,i]) 
rbreastcancer <- bootstrap_j48(BreastCancer, nboots=250)
apply(apply(rbreastcancer, 2, as.numeric),2,mean)
save(rbreastcancer, file="results2/breastcancer_j48.RData")

### 6. PimaIndiansDiabetes ###
data(PimaIndiansDiabetes)
rpima <- bootstrap_j48(PimaIndiansDiabetes, nboots=250)
apply(apply(rpima, 2, as.numeric),2,mean)
save(rpima, file="results2/pima_j48.RData")

### 7. Vowel ###
data(Vowel)
rvowel <- bootstrap_j48(Vowel, nboots=250)
apply(apply(rvowel, 2, as.numeric),2,mean)
save(rvowel, file="results2/vowel_j48.RData")

### 8. statlogGermanCredit ###
data(GermanCredit)
rcredit <- bootstrap_j48(GermanCredit, nboots=250)
apply(apply(rcredit, 2, as.numeric),2,mean)
save(rcredit, file="results2/credit_j48.RData")

### 9. contraceptiveMethodChoice ###
data(ContraceptiveChoice)
rcontraceptive <- bootstrap_j48(ContraceptiveChoice, nboots = 250)
apply(apply(rcontraceptive, 2, as.numeric),2,mean)
save(rcontraceptive, file="results2/contraceptive_j48.RData")

### 10. DNA ###
data(DNA)
rdna <- bootstrap_j48(DNA, nboots = 250)
apply(apply(rdna, 2, as.numeric),2,mean)
save(rdna, file="results2/dna_j48.RData")

### 11. spam ###
data(spam)
rspam <- bootstrap_j48(spam, nboots = 250)
apply(apply(rspam, 2, as.numeric),2,mean)
save(rspam, file="results2/spam_j48.RData")

### 12. MAGICGammaTelescope ###
data(MAGICGammaTelescope)
mg1 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 1000)
apply(apply(mg1, 2, as.numeric),2,mean)
save(mg1, file="results2/mg1_")

data(MAGICGammaTelescope)
mg2 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 2000)
apply(apply(mg2, 2, as.numeric),2,mean)
save(mg2, file="results2/mg2_")

data(MAGICGammaTelescope)
mg3 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 3000)
apply(apply(mg3, 2, as.numeric),2,mean)
save(mg3, file="results2/mg3_")

data(MAGICGammaTelescope)
mg4 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 4000)
apply(apply(mg4, 2, as.numeric),2,mean)
save(mg4, file="results2/mg4_")

data(MAGICGammaTelescope)
mg5 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 5000)
apply(apply(mg5, 2, as.numeric),2,mean)
save(mg5, file="results/mg5_")

load(file="results2/mg1_")
load(file="results2/mg2_")
load(file="results2/mg3_")
load(file="results2/mg4_")
load(file="results2/mg5_")
magicgamma <- rbind(mg1_, mg2_, mg3_, mg4_, mg5_)
save(magicgamma, file="results2/magicgamma_j48.RData")

### simulatied chessboard examples ###
rchessboard44_0 <- benchchessboard44_j48(noise = 0, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_0, 2, as.numeric),2,mean)
save(rchessboard44_0, file="results/chessboard44_0_j48.RData")

rchessboard44_5 <- benchchessboard44_j48(noise = 5, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_5, 2, as.numeric),2,mean)
save(rchessboard44_5, file="results/chessboard44_5_j48.RData")

rchessboard44_10 <- benchchessboard44_j48(noise = 10, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_10, 2, as.numeric),2,mean)
save(rchessboard44_10, file="results/chessboard44_10_j48.RData")


###################################################
### Benchmark plots
###################################################

### under construction 

## load results
for(i in Sys.glob("results/*.RData")) load(i)

## preprocess for reference evtree
preprocess <- function(d, dname = "datasetname", isclassification = TRUE){
    if(isclassification) d[, 1:3] <- 1 - d[ ,1:3]
    d <- as.data.frame(d)
    colnames(d) <- c("evtree", "J48", "J48_up","evtree", "J48", "J48_up")
    for(i in 3:1) d[, i] <- d[, i] / d[, 1] * 100
    for(i in 6:4) d[, i] <- d[, i] / d[, 4] * 100
    x <- d[, 1:3]
    y <- d[, 4:6]
    rval <- reshape(x, idvar="samp", times=names(x), timevar = "alg",varying= list(names(x)), direction="long")
    names(rval)[2] <- "accuracy"
    rval$complexity <- reshape(y, idvar="samp", times=names(y), timevar = "alg",varying= list(names(y)), direction="long")[,2]
    rval$alg <- factor(rval$alg, levels = c("evtree", "ctree", "rpart"))
    rval$ds <- dname
    rval
}

## collect results for all datasets
r <- rbind(
  preprocess(d = rglass, dname = "Glass identification", isclassification = TRUE),
  preprocess(d = rheart, dname = "Statlog heart", isclassification = TRUE),
  preprocess(d = rionosphere, dname = "Ionosphere", isclassification = TRUE),
  preprocess(d = rmusk, dname = "Musk", isclassification = TRUE),
  preprocess(d = rbreastcancer, dname = "Breast cancer database", isclassification = TRUE),
  preprocess(d = rpima, dname = "Pima Indians diabetes", isclassification = TRUE),
  preprocess(d = rvowel, dname = "Vowel", isclassification = TRUE),
  preprocess(d = rcredit, dname = "Statlog German credit", isclassification = TRUE),
  preprocess(d = rcontraceptive, dname = "Contraceptive method", isclassification = TRUE),
  preprocess(d = rdna, dname = "DNA", isclassification = TRUE),
  preprocess(d = rspam, dname = "Spam", isclassification = TRUE),
  preprocess(d = rmagicgamma, dname = "Magic gamma telescope", isclassification = TRUE),
)

r$ds <- factor(r$ds)
r$samp <- factor(r$samp)
r$dssamp <- r$ds:r$samp

## compute multiple comparisons
library("multcomp")
cstats <- function(alg = "rpart", value = "accuracy", data = r) {
  dlab <- rev(unique(data$ds))
  k <- length(dlab)
  mean  <- numeric(k)
  lower <- numeric(k)
  upper <- numeric(k)
  names(data)[names(data) == value] <- "value"
  for(i in 1:k) {
    mod1 <- lm(value ~ alg, data = subset(data, ds == dlab[i]))
    pt <- glht(mod1, linfct = mcp(alg = "Dunnett"))
    w <- confint(pt)$confint
    d <- which(levels(r$alg) == alg) - 2
    mean[i]  <-  w[1+d]
    lower[i] <-  w[3+d]
    upper[i] <-  w[5+d]
  }
  rval <- data.frame(mean, lower, upper)
  rownames(rval) <- dlab
  return(rval)
}

acc_rpart <- cstats("J48", "accuracy")
com_rpart <- cstats("J48_up", "complexity")
acc_ctree <- cstats("J48", "accuracy")
com_ctree <- cstats("J48_up", "complexity")

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
  abline(h = 5.5)
  abline(v = 0, lty = 2)  
  box()
}

## plot the results
par(mfrow = c(2, 2), oma = c(5, 10, 2, 0), mar = c(1, 1, 2, 1))

xlim1 <- range(cbind(acc_rpart, acc_ctree))
xlim2 <- range(cbind(com_rpart, com_ctree))

ciplot(acc_rpart, xlim = xlim1, main = "J48", ylab = TRUE, xlab = "")
ciplot(com_rpart, xlim = xlim2, main = "",      ylab = FALSE, xlab = "")
ciplot(acc_ctree, xlim = xlim1, main = "J48_up", ylab = TRUE,
  xlab = "relative difference in predictive accuracy (%)")
ciplot(com_ctree, xlim = xlim2, main = "",      ylab = FALSE,
  xlab = "relative difference in complexity (%)")

