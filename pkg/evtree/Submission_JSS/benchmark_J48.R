###################################################
### This file reproduces the benchmarks experiments and plots of the evtree paper.
### Note that the full analysis takes rather long (our analysis on a mainframe with 
### 6 parallel processes took 1-2 weeks). 
###################################################

## required libraries
library("J48")
library("party")
library("partykit")
library("kernlab")
library("mlbench")
library("evtree")
library("RWeka")

length_j48 <- function(x) {
	x <- capture.output(print(x))
	x <- x[grep("Size of the tree", x, fixed = TRUE)]
	x <- tail(strsplit(x, "\t")[[1]], 1)
	as.numeric(x)
}

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
	    NN_j48[f] <- length_j48(out_j48)
	    NN_j48_up[f] <- length_j48(out_j48_up)
		
    } # end for
    accuracy <- cbind(accuracy_j48, NN_j48, accuracy_j48_up, NN_j48_up)
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
    accuracy_j48 <- array(-999999, nrealizations)
    accuracy_j48_up <- array(-999999, nrealizations)
    NN_j48 <- array(-999999, nrealizations)
    NN_j48_up <- array(-999999, nrealizations)


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
	    NN_j48[f] <- length_j48(out_j48)
	    NN_j48_up[f] <- length_j48(out_j48_up)
		
    } # end for
    accuracy <- cbind(accuracy_j48, NN_j48, accuracy_j48_up, NN_j48_up)
    accuracy
}

###################################################
### Processing of datasets
###################################################

### 1. Glass ###
data(Glass)
rglass2 <- bootstrap_j48(Glass, nboots=250)
apply(apply(rglass2, 2, as.numeric),2,mean)
save(rglass2, file="results_j48/glass_j48.RData")

### 2. statlogHeart ###
data(StatlogHeart)
rheart2 <- bootstrap_j48(StatlogHeart, nboots=250)
apply(apply(rheart2, 2, as.numeric),2,mean)
save(rheart2, file="results_j48/heart_j48.RData")

### 3. Ionosphere ###
data(Ionosphere)
rionosphere2 <- bootstrap_j48(Ionosphere, nboots=250)
apply(apply(rionosphere2, 2, as.numeric),2,mean)
save(rionosphere2, file="results_j48/ionosphere_j48.RData")

### 4. musk ###
data(musk)
rmusk2 <- bootstrap_j48(musk, nboots=250)
apply(apply(rmusk2, 2, as.numeric),2,mean)
save(rmusk2, file="results_j48/musk.RData_j48")

### 5. BreastCancer ###
data(BreastCancer)
BreastCancer <- BreastCancer[complete.cases(BreastCancer),2:11]
for(i in 1:5)
	BreastCancer[,i] <- as.numeric(BreastCancer[,i]) 
for(i in 6:10)
	BreastCancer[,i] <- as.factor(BreastCancer[,i]) 
rbreastcancer2 <- bootstrap_j48(BreastCancer, nboots=250)
apply(apply(rbreastcancer2, 2, as.numeric),2,mean)
save(rbreastcancer2, file="results_j48/breastcancer_j48.RData")

### 6. PimaIndiansDiabetes ###
data(PimaIndiansDiabetes)
rpima2 <- bootstrap_j48(PimaIndiansDiabetes, nboots=250)
apply(apply(rpima2, 2, as.numeric),2,mean)
save(rpima2, file="results_j48/pima_j48.RData")

### 7. Vowel ###
data(Vowel)
rvowel2 <- bootstrap_j48(Vowel, nboots=250)
apply(apply(rvowel2, 2, as.numeric),2,mean)
save(rvowel2, file="results_j48/vowel_j48.RData")

### 8. statlogGermanCredit ###
data(GermanCredit)
rcredit2 <- bootstrap_j48(GermanCredit, nboots=250)
apply(apply(rcredit2, 2, as.numeric),2,mean)
save(rcredit2, file="results_j48/credit_j48.RData")

### 9. contraceptiveMethodChoice ###
data(ContraceptiveChoice)
rcontraceptive2 <- bootstrap_j48(ContraceptiveChoice, nboots = 250)
apply(apply(rcontraceptive2, 2, as.numeric),2,mean)
save(rcontraceptive2, file="results_j48/contraceptive_j48.RData")

### 10. DNA ###
data(DNA)
rdna2 <- bootstrap_j48(DNA, nboots = 250)
apply(apply(rdna2, 2, as.numeric),2,mean)
save(rdna2, file="results_j48/dna_j48.RData")

### 11. spam ###
data(spam)
rspam2 <- bootstrap_j48(spam, nboots = 250)
apply(apply(rspam2, 2, as.numeric),2,mean)
save(rspam2, file="results_j48/spam_j48.RData")

### 12. MAGICGammaTelescope ###
data(MAGICGammaTelescope)
mg1 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 1000)
apply(apply(mg1, 2, as.numeric),2,mean)
save(mg1, file="results_j48/mg1_")

data(MAGICGammaTelescope)
mg2 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 2000)
apply(apply(mg2, 2, as.numeric),2,mean)
save(mg2, file="results_j48/mg2_")

data(MAGICGammaTelescope)
mg3 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 3000)
apply(apply(mg3, 2, as.numeric),2,mean)
save(mg3, file="results_j48/mg3_")

data(MAGICGammaTelescope)
mg4 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 4000)
apply(apply(mg4, 2, as.numeric),2,mean)
save(mg4, file="results_j48/mg4_")

data(MAGICGammaTelescope)
mg5 <- bootstrap_j48(MAGICGammaTelescope, nboots = 50, seed = 5000)
apply(apply(mg5, 2, as.numeric),2,mean)
save(mg5, file="results_j48/mg5_")

load(file="results_j48/mg1_")
load(file="results_j48/mg2_")
load(file="results_j48/mg3_")
load(file="results_j48/mg4_")
load(file="results_j48/mg5_")
rmagicgamma2 <- rbind(mg1, mg2, mg3, mg4, mg5)
save(rmagicgamma2, file="results_j48/magicgamma_j48.RData")

### simulatied chessboard examples ###
rchessboard44_02 <- benchchessboard44_j48(noise = 0, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_02, 2, as.numeric),2,mean)
save(rchessboard44_02, file="results_j48/chessboard44_0_j48.RData")

rchessboard44_52 <- benchchessboard44_j48(noise = 5, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_52, 2, as.numeric),2,mean)
save(rchessboard44_52, file="results_j48/chessboard44_5_j48.RData")

rchessboard44_102 <- benchchessboard44_j48(noise = 10, seed = 1000, nrealizations = 250)
apply(apply(rchessboard44_102, 2, as.numeric),2,mean)
save(rchessboard44_102, file="results_j48/chessboard44_10_j48.RData")
