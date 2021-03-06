###################################################
### 
###################################################

## required libraries
library("rpart")
library("party")
library("partykit")
library("kernlab")
library("mlbench")
library("evtree")

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

###################################################
### Calculates results of benchmark experiments 
### Input: A benchmark dataset (the dependent variable has to be in the last row).
### Output: Results of the benchmark.
###################################################

# experiments number of iterations and population size
# kdata=1 for simulation of 4x4 chesboard examples
bootstrap_evtree_parameter <- function(kdata=1, nboots = 50, seed = 1000){
    set.seed(seed)
   	# results evtree with different population size and number of iterations: _#trees_#iterations
	accuracy_evtree_25  <- array(-999999, nboots)	    
	accuracy_evtree_50  <- array(-999999, nboots)	
    accuracy_evtree_100 <- array(-999999, nboots)
    accuracy_evtree_250 <- array(-999999, nboots)
    accuracy_evtree_500 <- array(-999999, nboots)
	
	# results evtree with different operator probabilities: _%crossover_%mutation_%split/prune_#iterations
    accuracy_evtree_0_50_50_200 <- array(-999999, nboots)
    accuracy_evtree_20_40_40_200 <- array(-999999, nboots)
    accuracy_evtree_40_30_30_200 <- array(-999999, nboots)
    accuracy_evtree_20_20_60_200 <- array(-999999, nboots)
    accuracy_evtree_20_60_20_200 <- array(-999999, nboots)
    
    accuracy_evtree_0_50_50_500 <- array(-999999, nboots)
    accuracy_evtree_20_40_40_500 <- array(-999999, nboots)
    accuracy_evtree_40_30_30_500 <- array(-999999, nboots)
    accuracy_evtree_20_20_60_500 <- array(-999999, nboots)
    accuracy_evtree_20_60_20_500 <- array(-999999, nboots)
    
    accuracy_evtree_0_50_50_10000 <- array(-999999, nboots)
    accuracy_evtree_20_40_40_10000 <- array(-999999, nboots)
    accuracy_evtree_40_30_30_10000 <- array(-999999, nboots)
    accuracy_evtree_20_20_60_10000 <- array(-999999, nboots)
    accuracy_evtree_20_60_20_10000 <- array(-999999, nboots)
	
	NN_evtree_25 <- array(-999999, nboots)	
    NN_evtree_50 <- array(-999999, nboots)
    NN_evtree_100 <- array(-999999, nboots)
    NN_evtree_250 <- array(-999999, nboots)
    NN_evtree_500 <- array(-999999, nboots)
    
    NN_evtree_0_50_50_200 <- array(-999999, nboots)
    NN_evtree_20_40_40_200 <- array(-999999, nboots)
    NN_evtree_40_30_30_200 <- array(-999999, nboots)
    NN_evtree_20_20_60_200 <- array(-999999, nboots)
    NN_evtree_20_60_20_200 <- array(-999999, nboots)
    
    NN_evtree_0_50_50_500 <- array(-999999, nboots)
    NN_evtree_20_40_40_500 <- array(-999999, nboots)
    NN_evtree_40_30_30_500 <- array(-999999, nboots)
    NN_evtree_20_20_60_500 <- array(-999999, nboots)
    NN_evtree_20_60_20_500 <- array(-999999, nboots)
    
    NN_evtree_0_50_50_10000 <- array(-999999, nboots)
    NN_evtree_20_40_40_10000 <- array(-999999, nboots)
    NN_evtree_40_30_30_10000 <- array(-999999, nboots)
    NN_evtree_20_20_60_10000 <- array(-999999, nboots)
    NN_evtree_20_60_20_10000 <- array(-999999, nboots)
	
    print("bootstrap no.:")
    for(f in 1:nboots){
        print(f)
        if(length(kdata == 1) && kdata == 1){ #dataset missing; simulate 4x4 chessboard example
          print("chessboard")
            n = 4000
            ch <- chessboard44(n = 4000, noise = 5 )
     	    ch_train <- ch[1:(n/2), ]
       	    ch_test <-  ch[(n/2+1):n, ]
     	    xtrain <-   ch_train[, 1:(ncol(ch_train)-1)]
        	ytrain <-   ch_train[, ncol(ch_train)]
        	xpredict <- ch_test[, 1:(ncol(ch_test)-1)]
        	ypredict <- ch_test[, ncol(ch_test)]
         }else{ #real dataset
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
		}

        out_evtree_25 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 25L
						)
        )
		out_evtree_50 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 50L
						)
        )
		out_evtree_100 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 100L
						)
        )
		out_evtree_250 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 250L
						)
        )
		out_evtree_500 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 500L
						)
        )		
		        
      	out_evtree_0_50_50_200 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 200L, ntrees = 100L,
						operatorprob = list(pcrossover = 0., pmutatemajor = 0.25, pmutateminor = 0.25,  psplit = 0.25, pprune = 0.25),
						)
        )
       	out_evtree_20_40_40_200 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 200L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.2, pmutateminor = 0.2,  psplit = 0.2, pprune = 0.2),
						)
        )
   		out_evtree_40_30_30_200 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 200L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.4, pmutatemajor = 0.15, pmutateminor = 0.15,  psplit = 0.15, pprune = 0.15),
						)
        )
        out_evtree_20_20_60_200 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 200L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.1, pmutateminor = 0.1,  psplit = 0.3, pprune = 0.3),
						)
        )  
        out_evtree_20_60_20_200 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 200L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.3, pmutateminor = 0.3,  psplit = 0.1, pprune = 0.1),
						)
        )                   
		
		
	    out_evtree_0_50_50_500 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 500L, ntrees = 100L,
						operatorprob = list(pcrossover = 0., pmutatemajor = 0.25, pmutateminor = 0.25,  psplit = 0.25, pprune = 0.25),
						)
        )
       	out_evtree_20_40_40_500 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 500L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.2, pmutateminor = 0.2,  psplit = 0.2, pprune = 0.2),
						)
        )
   		out_evtree_40_30_30_500 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 500L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.4, pmutatemajor = 0.15, pmutateminor = 0.15,  psplit = 0.15, pprune = 0.15),
						)
        )
        out_evtree_20_20_60_500 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 500L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.1, pmutateminor = 0.1,  psplit = 0.3, pprune = 0.3),
						)
        )  
        out_evtree_20_60_20_500 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						niterations = 500L, ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.3, pmutateminor = 0.3,  psplit = 0.1, pprune = 0.1),
						)
        )        
		
		
	    out_evtree_0_50_50_10000 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 100L,
						operatorprob = list(pcrossover = 0., pmutatemajor = 0.25, pmutateminor = 0.25,  psplit = 0.25, pprune = 0.25),
						)
        )
       	out_evtree_20_40_40_10000 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.2, pmutateminor = 0.2,  psplit = 0.2, pprune = 0.2),
						)
        )
   		out_evtree_40_30_30_10000 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 100L,
						operatorprob = list(pcrossover = 0.4, pmutatemajor = 0.15, pmutateminor = 0.15,  psplit = 0.15, pprune = 0.15),
						)
        )
        out_evtree_20_20_60_10000 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.1, pmutateminor = 0.1,  psplit = 0.3, pprune = 0.3),
						)
        )  
        out_evtree_20_60_20_10000 <- evtree(ytrain ~ .,xtrain,
                        control = evtree.control(minbucket = 7, minsplit = 20, maxdepth = 9, seed = seed,
						ntrees = 100L,
						operatorprob = list(pcrossover = 0.2, pmutatemajor = 0.3, pmutateminor = 0.3,  psplit = 0.1, pprune = 0.1),
						)
        )        		
                
		
        accuracy_evtree_25[f] <- sum( predict(out_evtree_25, newdata= xpredict) == ypredict) / length(ypredict)
        accuracy_evtree_50[f] <- sum( predict(out_evtree_50, newdata= xpredict) == ypredict) / length(ypredict)
        accuracy_evtree_100[f] <- sum( predict(out_evtree_100, newdata= xpredict) == ypredict) / length(ypredict)
        accuracy_evtree_250[f] <- sum( predict(out_evtree_250, newdata= xpredict) == ypredict) / length(ypredict)
        accuracy_evtree_500[f] <- sum( predict(out_evtree_500, newdata= xpredict) == ypredict) / length(ypredict)

        accuracy_evtree_0_50_50_200[f] <- sum( predict(out_evtree_0_50_50_200, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_40_40_200[f] <- sum( predict(out_evtree_20_40_40_200, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_40_30_30_200[f] <- sum( predict(out_evtree_40_30_30_200, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_20_60_200[f] <- sum( predict(out_evtree_20_20_60_200, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_60_20_200[f] <- sum( predict(out_evtree_20_60_20_200, newdata= xpredict) == ypredict)/ length(ypredict)

        accuracy_evtree_0_50_50_500[f] <- sum( predict(out_evtree_0_50_50_500, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_40_40_500[f] <- sum( predict(out_evtree_20_40_40_500, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_40_30_30_500[f] <- sum( predict(out_evtree_40_30_30_500, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_20_60_500[f] <- sum( predict(out_evtree_20_20_60_500, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_60_20_500[f] <- sum( predict(out_evtree_20_60_20_500, newdata= xpredict) == ypredict)/ length(ypredict)

        accuracy_evtree_0_50_50_10000[f] <- sum( predict(out_evtree_0_50_50_10000, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_40_40_10000[f] <- sum( predict(out_evtree_20_40_40_10000, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_40_30_30_10000[f] <- sum( predict(out_evtree_40_30_30_10000, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_20_60_10000[f] <- sum( predict(out_evtree_20_20_60_10000, newdata= xpredict) == ypredict)/ length(ypredict)
        accuracy_evtree_20_60_20_10000[f] <- sum( predict(out_evtree_20_60_20_10000, newdata= xpredict) == ypredict)/ length(ypredict)
             
        NN_evtree_25[f] <- width(node_party(out_evtree_25))
	    NN_evtree_50[f] <- width(node_party(out_evtree_50))
        NN_evtree_100[f] <- width(node_party(out_evtree_100))
        NN_evtree_250[f] <- width(node_party(out_evtree_250))
        NN_evtree_500[f] <- width(node_party(out_evtree_500))
		
		NN_evtree_0_50_50_200[f] <- width(node_party(out_evtree_0_50_50_200))
	    NN_evtree_20_40_40_200[f] <- width(node_party(out_evtree_20_40_40_200))
        NN_evtree_40_30_30_200[f] <- width(node_party(out_evtree_40_30_30_200))
        NN_evtree_20_20_60_200[f] <- width(node_party(out_evtree_20_20_60_200))
        NN_evtree_20_60_20_200[f] <- width(node_party(out_evtree_20_60_20_200))
        
      	NN_evtree_0_50_50_500[f] <- width(node_party(out_evtree_0_50_50_500))
	    NN_evtree_20_40_40_500[f] <- width(node_party(out_evtree_20_40_40_500))
        NN_evtree_40_30_30_500[f] <- width(node_party(out_evtree_40_30_30_500))
        NN_evtree_20_20_60_500[f] <- width(node_party(out_evtree_20_20_60_500))
        NN_evtree_20_60_20_500[f] <- width(node_party(out_evtree_20_60_20_500))
        
		NN_evtree_0_50_50_10000[f] <- width(node_party(out_evtree_0_50_50_10000))
	    NN_evtree_20_40_40_10000[f] <- width(node_party(out_evtree_20_40_40_10000))
        NN_evtree_40_30_30_10000[f] <- width(node_party(out_evtree_40_30_30_10000))
        NN_evtree_20_20_60_10000[f] <- width(node_party(out_evtree_20_20_60_10000))
        NN_evtree_20_60_20_10000[f] <- width(node_party(out_evtree_20_60_20_10000))
               
    } # end for
	
    accuracy <- cbind(
	accuracy_evtree_25, accuracy_evtree_50, accuracy_evtree_100,
    accuracy_evtree_250, accuracy_evtree_500,
	
	NN_evtree_25, NN_evtree_25, NN_evtree_50,
    NN_evtree_100, NN_evtree_250, NN_evtree_500,
	
	accuracy_evtree_0_50_50_200, accuracy_evtree_20_40_40_200, accuracy_evtree_40_30_30_200,
    accuracy_evtree_20_20_60_200, accuracy_evtree_20_60_20_200, 
    
	accuracy_evtree_0_50_50_500, accuracy_evtree_20_40_40_500, accuracy_evtree_40_30_30_500,
    accuracy_evtree_20_20_60_500, accuracy_evtree_20_60_20_500, 
    
    accuracy_evtree_0_50_50_10000, accuracy_evtree_20_40_40_10000, accuracy_evtree_40_30_30_10000,
    accuracy_evtree_20_20_60_10000, accuracy_evtree_20_60_20_10000, 
    
	NN_evtree_0_50_50_200, NN_evtree_20_40_40_200, NN_evtree_40_30_30_200,
    NN_evtree_20_20_60_200, NN_evtree_20_60_20_200, 
    
	NN_evtree_0_50_50_500, NN_evtree_20_40_40_500, NN_evtree_40_30_30_500,
    NN_evtree_20_20_60_500, NN_evtree_20_60_20_500, 
    
    NN_evtree_0_50_50_10000, NN_evtree_20_40_40_10000, NN_evtree_40_30_30_10000,
    NN_evtree_20_20_60_10000, NN_evtree_20_60_20_10000 
	)
    accuracy
}

###################################################
### Processing of datasets
###################################################

### 2. statlogHeart ###
data(StatlogHeart)
rheart <- bootstrap_evtree_parameter(StatlogHeart, nboots = 250)
apply(apply(rheart, 2, as.numeric),2,mean)
save(rheart, file="results_parameter/heart_parameter.RData")

### 8. statlogGermanCredit ###
data(GermanCredit)
rcredit <- bootstrap_evtree_parameter(GermanCredit, nboots = 250)
apply(apply(rcredit, 2, as.numeric),2,mean)
save(rcredit, file="results_parameter/credit_parameter.RData")

### 11. spam ###
data(spam)
rspam <- bootstrap_evtree_parameter(spam, nboots = 250)
apply(apply(rspam, 2, as.numeric),2,mean)
save(rspam, file="results_parameter/spam_parameter.RData")

### 3. 4x4 chessboard with 5% noise ###
rchessboard44_5 <- bootstrap_evtree_parameter(kdata = 1, nboots = 250)
apply(apply(rchessboard44_5, 2, as.numeric),2,mean)
save(rchessboard44_5, file="results_parameter/chessboard44_5_parameter.RData")
