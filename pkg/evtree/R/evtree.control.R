evtree.control <- function(minbucket = 7L, minsplit = 20L, maxdepth = 9L, niterations = 10000L, ntrees = 100L, alpha = 1,
  operatorprob= list(pmutatemajor = 20, pmutateminor = 20, pcrossover = 20, psplit = 20, pprune = 20), seed = NULL){
    list(minbucket = minbucket,
         minsplit = minsplit,
         maxdepth = maxdepth,
         niterations = niterations,
         ntrees = ntrees,
         alpha = alpha,
         operatorprob = operatorprob,
	 seed = seed 
    )
}
