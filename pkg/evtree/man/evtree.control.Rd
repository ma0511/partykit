\name{evtree.control}
\alias{ evtree.control}
\title{ Control for evtree}
\description{
  Various parameters that control aspects of the `evtree' fit.
}
\usage{
evtree.control(minbucket = 7L, minsplit = 20L, maxdepth = 9L, niterations = 10000L, ntrees = 100L, alpha = 1, operatorprob = list(pmutatemajor = 0.2, pmutateminor = 0.2, pcrossover = 0.2, psplit = 0.2, pprune = 0.2), seed = NULL, ...)
}

\arguments{
  \item{minbucket}{ the minimum sum of weights in a terminal node. }
  \item{minsplit}{ the minimum sum of weights in a node in order to be considered for splitting. }
  \item{maxdepth}{ maximum depth of the tree. Note, that the memory requirements increase by the square of the maximum tree depth. The default \code{maxdepth = 9}. }
  \item{niterations}{ in case the run does not converge, it terminates after a specified number of iterations defined by \code{niterations}. The default \code{niterations = 10000}. }
  \item{ntrees}{ the number of trees in the population. The default \code{ntrees = 100}. }
  \item{alpha}{ regulates the complexity part of the cost function. Increasing values of alpha encourage decreasing tree sizes. 
The default \code{alpha = 1}.}
  \item{operatorprob}{ probabilities for the selection of variation operators. The default 
\code{list(pmutatemajor = 0.2, pmutateminor = 0.2, pcrossover = 0.2, psplit = 0.2, pprune = 0.2)}. 
The probabilities must sum up to 1.} 
  \item{seed}{ a seed to initialize the random number generator. By default the random number generator is initialized by the system time. Results can be reproduced by re-using the same seed}
}

\value{
  A list.
}

\keyword{misc}
