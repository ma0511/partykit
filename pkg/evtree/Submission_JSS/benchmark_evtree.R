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

###################################################
### Calculates results of benchmark experiments.
### Input: A benchmark dataset (the dependent variable has to be in the last row).
### Output: Results of the benchmark.
###################################################
bootstrap_evtree <- function(kdata, nboots = 250, seed = 1000){
    set.seed(seed)
    accuracy_evtree <- array(-999999, nboots)
    accuracy_crpart <- array(-999999, nboots)
    accuracy_ctree <- array(-999999, nboots)

    NN_evtree <- array(-999999, nboots)
    NN_rpart   <- array(-999999, nboots)
    NN_ctree   <- array(-999999, nboots)

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
        out_ctree <- party::ctree(
                ytrain ~ ., data = cbind(xtrain, ytrain),
                control = party::ctree_control(
                       minbucket = 7,
                       minsplit = 20,
                       maxdepth = 9
                )
        )

        out_rpart <- rpart(
                ytrain ~ ., data = xtrain,
                        minbucket = 7,
                        minsplit = 20,
                        maxdepth = 9,
                        maxsurrogate = 0,
                        maxcompete = 0,
                        )

        out_evtree <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(
                        minbucket = 7,
                        minsplit = 20,
                        maxdepth = 9,
                        seed = seed
                        )
        )

        fit_ctree <- predict(out_ctree, newdata= xpredict)
        fit_evtree <- predict(out_evtree, newdata= xpredict)
        if(is.factor(ypredict)){
            fit_rpart <- predict(out_rpart, newdata= xpredict, type = c("class"))
            accuracy_evtree[f] <- sum(fit_evtree == ypredict) / length(ypredict)
            accuracy_crpart[f] <- sum(fit_rpart == ypredict) / length(ypredict)
            accuracy_ctree[f] <- sum(fit_ctree == ypredict) / length(ypredict)
        }else{
            fit_rpart <- predict(out_rpart, newdata= xpredict)
            accuracy_evtree[f] <- sum((fit_evtree - ypredict)^2) / length(ypredict)
            accuracy_crpart[f] <- sum((fit_rpart - ypredict)^2) / length(ypredict)
            accuracy_ctree[f] <- sum((fit_ctree - ypredict)^2) / length(ypredict)
        }

        NN_evtree[f] <- width(node_party(out_evtree))
        NN_rpart[f] <- length(out_rpart$splits[,4])+1
        NN_ctree[f] <- (max(where(out_ctree))+1)/2
               
    } # end for
    accuracy <- cbind(accuracy_evtree, accuracy_crpart, accuracy_ctree, NN_evtree, NN_rpart, NN_ctree)
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

benchchessboard44 <- function(n = 4000, noisevariables = 6, noise = 0, nrealizations = 250, seed = 1000){
    set.seed(seed)
    accuracy_evtree <- array(-999999, nrealizations)
    accuracy_rpart <- array(-999999, nrealizations)
    accuracy_ctree <- array(-999999, nrealizations)

    NN_evtree <- array(-999999, nrealizations)
    NN_rpart   <- array(-999999, nrealizations)
    NN_ctree   <- array(-999999, nrealizations)

    print("realization no.:")
    for(f in 1:nrealizations){
        print(f)

        ch <- chessboard44(n = n, noisevariables = noisevariables, noise = noise )
        ch_train <- ch[1:(n/2), ]
        ch_test <-  ch[(n/2+1):n, ]
        xtrain <-   ch_train[, 1:(ncol(ch_train)-1)]
        ytrain <-   ch_train[, ncol(ch_train)]
        xpredict <- ch_test[, 1:(ncol(ch_test)-1)]
        ypredict <- ch_test[, ncol(ch_test)]

        out_ctree <- party::ctree(
                ytrain ~ ., data = cbind(xtrain, ytrain),
                control = party::ctree_control(
                       minbucket = 7,
                       minsplit = 20,
                       maxdepth = 9
                )
        )

        out_rpart <- rpart(
                ytrain ~ ., data = xtrain,
                        minbucket = 7,
                        minsplit = 20,
                        maxdepth = 9,
                        maxsurrogate = 0,
                        maxcompete = 0
                        )

        out_evtree <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(
                        minbucket = 7,
                        minsplit = 20,
                        maxdepth = 9,
                        alpha=1,
                        seed = seed
                        )
        )

        fit_ctree <- predict(out_ctree, newdata = xpredict)
        fit_evtree <- predict(out_evtree, newdata = xpredict)
        fit_rpart <- predict(out_rpart, newdata = xpredict, type = c("class"))

        accuracy_evtree[f] <- sum(fit_evtree == ypredict) / length(ypredict)
        accuracy_rpart[f] <- sum(fit_rpart == ypredict) / length(ypredict)
        accuracy_ctree[f] <- sum(fit_ctree == ypredict) / length(ypredict)

        NN_evtree[f] <- width(node_party(out_evtree))
        NN_rpart[f] <- length(out_rpart$splits[,4])+1
        NN_ctree[f] <- (max(where(out_ctree))+1)/2
    } # end for
    accuracy <- cbind(accuracy_evtree, accuracy_rpart, accuracy_ctree, NN_evtree, NN_rpart, NN_ctree)
    accuracy
}

###################################################
### Processing of datasets
###################################################

### 1. Glass ###
data(Glass)
rglass <- bootstrap_evtree(Glass, nboots=250)
apply(apply(rglass, 2, as.numeric),2,mean)
save(rglass, file="results/glass.RData")

### 2. statlogHeart ###
data(StatlogHeart)
rheart <- bootstrap_evtree(StatlogHeart, nboots=250)
apply(apply(rheart, 2, as.numeric),2,mean)
save(rheart, file="results/heart.RData")

### 3. Ionosphere ###
data(Ionosphere)
rionosphere <- bootstrap_evtree(Ionosphere, nboots=250)
apply(apply(rionosphere, 2, as.numeric),2,mean)
save(rionosphere, file="results/ionosphere.RData")

### 4. musk ###
data(musk)
rmusk <- bootstrap_evtree(musk, nboots=250)
apply(apply(rmusk, 2, as.numeric),2,mean)
save(rmusk, file="results/musk.RData")

### 5. BreastCancer ###
data(BreastCancer)
BreastCancer <- BreastCancer[complete.cases(BreastCancer),2:11]
for(i in 1:5)
	BreastCancer[,i] <- as.numeric(BreastCancer[,i]) 
for(i in 6:10)
	BreastCancer[,i] <- as.factor(BreastCancer[,i]) 
rbreastcancer <- bootstrap_evtree(BreastCancer, nboots=250)
apply(apply(rbreastcancer, 2, as.numeric),2,mean)
save(rbreastcancer, file="results/breastcancer.RData")

### 6. PimaIndiansDiabetes ###
data(PimaIndiansDiabetes)
rpima <- bootstrap_evtree(PimaIndiansDiabetes, nboots=250)
apply(apply(rpima, 2, as.numeric),2,mean)
save(rpima, file="results/pima.RData")

### 7. Vowel ###
data(Vowel)
rvowel <- bootstrap_evtree(Vowel, nboots=250)
apply(apply(rvowel, 2, as.numeric),2,mean)
save(rvowel, file="results/vowel.RData")

### 8. statlogGermanCredit ###
data(GermanCredit)
rcredit <- bootstrap_evtree(GermanCredit, nboots=250)
apply(apply(rcredit, 2, as.numeric),2,mean)
save(rcredit, file="results/credit.RData")

### 9. contraceptiveMethodChoice ###
data(ContraceptiveChoice)
rcontraceptive <- bootstrap_evtree(ContraceptiveChoice, nboots = 250)
apply(apply(rcontraceptive, 2, as.numeric),2,mean)
save(rcontraceptive, file="results/contraceptive.RData")

### 10. DNA ###
data(DNA)
rdna <- bootstrap_evtree(DNA, nboots = 250)
apply(apply(rdna, 2, as.numeric),2,mean)
save(rdna, file="results/dna.RData")

### 11. spam ###
data(spam)
rspam <- bootstrap_evtree(spam, nboots = 250)
apply(apply(rspam, 2, as.numeric),2,mean)
save(rspam, file="results/spam.RData")

### 12. MAGICGammaTelescope ###
# The analysis of this dataset takes extremly long. We therefore splitted the work into 5 seperate runs.
data(MAGICGammaTelescope)
mg1 <- bootstrap_evtree(MAGICGammaTelescope, nboots = 50, seed = 1000)
apply(apply(mg1, 2, as.numeric),2,mean)
save(mg1, file="results/mg1")

data(MAGICGammaTelescope)
mg2 <- bootstrap_evtree(MAGICGammaTelescope, nboots = 50, seed = 2000)
apply(apply(mg2, 2, as.numeric),2,mean)
save(mg2, file="results/mg2")

data(MAGICGammaTelescope)
mg3 <- bootstrap_evtree(MAGICGammaTelescope, nboots = 50, seed = 3000)
apply(apply(mg3, 2, as.numeric),2,mean)
save(mg3, file="results/mg3")

data(MAGICGammaTelescope)
mg4 <- bootstrap_evtree(MAGICGammaTelescope, nboots = 50, seed = 4000)
apply(apply(mg4, 2, as.numeric),2,mean)
save(mg4, file="results/mg4")

data(MAGICGammaTelescope)
mg5 <- bootstrap_evtree(MAGICGammaTelescope, nboots = 50, seed = 5000)
apply(apply(mg5, 2, as.numeric),2,mean)
save(mg5, file="results/mg5")

load(file="results/mg1")
load(file="results/mg2")
load(file="results/mg3")
load(file="results/mg4")
load(file="results/mg5")
magicgamma <- rbind(mg1, mg2, mg3, mg4, mg5)
save(magicgamma, file="results/magicgamma.RData")

### 13. Servo ###
data("Servo")
rservo <- bootstrap_evtree(Servo, nboots = 250)
apply(apply(rservo, 2, as.numeric),2,mean)
save(rservo, file="results/servo.RData")

### 14. BostonHousing ###
data("BostonHousing")
rbostonhousing <- bootstrap_evtree(BostonHousing, nboots = 250)
apply(apply(rbostonhousing, 2, as.numeric),2,mean)
save(rbostonhousing, file="results/bostonhousing.RData")

### Datasets from the Austrian DRG system. We are not allowed to make them available.

#load("/home/grubinger/Desktop/mel0101.rda" )
#rmel0101 <- bootstrap_evtree(adata[,-c(4,6)], nboots = 250, seed = 1000)
#save(rmel0101, file="results/mel0101.RData")
#
#load("/home/grubinger/Desktop/hdg0502.rda" )
#rhdg0502 <- bootstrap_evtree(adata, nboots= 250, seed=1000)
#save(rhdg0502, file="results/hdg0502.RData")
#
#load("/home/grubinger/Desktop/hdg0202.rda" )
#rhdg0202 <- bootstrap_evtree(adata[,-c(10,46,47,50)], nboots = 250, seed = 1000)
#save(rhdg0202, file = "results/hdg0202.RData")

### simulatied chessboard examples ###
rchessboard44_0 <- benchchessboard44(noise = 0, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_0, 2, as.numeric),2,mean)
save(rchessboard44_0, file="results/chessboard44_0.RData")

rchessboard44_5 <- benchchessboard44(noise = 5, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_5, 2, as.numeric),2,mean)
save(rchessboard44_5, file="results/chessboard44_5.RData")

rchessboard44_10 <- benchchessboard44(noise = 10, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_10, 2, as.numeric),2,mean)
save(rchessboard44_10, file="results/chessboard44_10.RData")
