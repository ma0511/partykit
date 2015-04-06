\name{lmertree}
\alias{lmertree}

\title{Linear Mixed Model Trees}

\description{
  Model-based recursive partitioning based on linear mixed model.
}

\usage{
lmertree(lmtreeformula, randomformula, data,
  initialRandomEffects = NULL, 
  ErrorTolerance = 0.001, MaxIterations = 1000, 
  verbose = TRUE, plotting = FALSE,
  lmer.control = lmerControl(), \dots)
}



\arguments{
  \item{lmtreeformula}{Model formula for building the \code{lmtree}. See examples below, or \code{lmtree} documentation.} 
  \item{randomformula}{Formula for random-effects model. See examples below, or \code{lmer} and \code{glmer} documentation.}
  \item{data}{Dataset to be used for estimating the lmertree.}
  \item{initialRandomEffects}{A vector of length \code{nrow(data)}, to be used as an offset in estimation of the first lmtree. \code{NULL} by default, which results in no offset being used.}
  \item{ErrorTolerance}{The convergence criterion used for estimation of the lmertree. When the difference in log likelihoods of the random-effects model of two consecutive iterations is smaller than \code{ErrorTolerance}, estimation of the lmertree has converged.} 
  \item{MaxIterations}{The maximum number of iterations to be performed in executing the lmertree function.}
  \item{verbose}{Should the log likelihood value of the estimated random-effects model be printed for every iteration during execution of the lmertree function?} 
  \item{plotting}{Should the lmertree be plotted at every iteration during execution of the lmertree function? Note that plotting of a tree at every iteration slows down execution of the function, and using \code{plotting=F} is recommended.}
  \item{lmer.control}{An optional list (of correct class, resulting from \code{lmerControl()}), containing control parameters to be passed to \code{lmer()}. See \code{lmerControl} documentation for details.} 
  \item{\dots}{Additional arguments to be passed to \code{lmtree()}. See \code{mob_control} documentation for details.}
}
  

\details{
  This code is the first implementation and will probably change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{Tree}{The final linear model tree.}
  \item{RandomEffectModel}{The final random-effects model.}
  \item{RandomEffects}{ranef(lme)} 
  \item{BetweenMatrix}{VarCorr(lme)}
  \item{ErrorVariance}{attr(VarCorr(lme),"sc")^2} 
  \item{data}{The dataset specified with the \code{dataset} argument.}
  \item{logLik}{The log likelihood value of the last iteration.}
  \item{IterationsUsed}{The number of iterations used to estimate the lmertree.} 
  \item{MaxIterations}{The maximum number of iterations specified with the \code{MaxIterations} argument.}
  \item{initialRandomEffects}{The random effects used as an offset, as specified with the \code{initialRandomEffects} argument.} 
  \item{TreeFormula}{The formula as specified with the \code{lmtreeformula} argument.}
  \item{RandomFormula}{The formula as specified with the \code{randomformula} argument.}
  \item{ErrorTolerance}{The prespecified value for the change in log likelihood to evaluate convergence, as specified with the \code{ErrorTolerance} argument.} 
  \item{tree.control}{A list containing control parameters passed to \code{lmtree()}, as specified with \dots.}
  \item{lmer.control}{A list containing control parameters passed to \code{lmer()}, as specified in the \code{lmer.control} argument.}
}

\seealso{\code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\link[partykit]{lmtree}}, \code{\link[partykit]{glmtree}}}

\examples{
## lmertree example
data(dataset.lmer)
lmer_ex <- lmertree(Y ~ treatment | U1 + U2 + U3 + U4 + U5, Y ~ 1|cluster, 
  data=dataset.lmer, plotting=FALSE)
plot(lmer_ex$Tree)
coef(lmer_ex$Tree)
lmer_ex$RandomEffectModel
lmer_ex$RandomEffects
}

\keyword{tree}
