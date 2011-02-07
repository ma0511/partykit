evtree.control <- function(minbucket = 20L, minsplit = 7L, maxdepth = 9L, niterations = 10000L, ntrees = 100L, alpha = 1,
  operatorprob= list(pmutatemajor = 20, pmutateminor = 20, pcrossover = 20, psplit = 20, pprune = 20)){
    list(minbucket = minbucket,
         minsplit = minsplit,
         maxdepth = maxdepth,
         niterations = niterations,
         ntrees = ntrees,
         operatorprob = operatorprob,
         alpha = alpha
    )
}
