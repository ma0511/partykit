\name{lmertree}
\alias{lmertree}
\alias{print.lmertree}
\alias{plot.lmertree}
\alias{coef.lmertree}
\alias{ranef.lmertree}

\title{Linear Mixed Model Trees}

\description{
  Model-based recursive partitioning based on linear mixed model.
}

\usage{
lmertree(formula, randomformula, data,
  ranefstart = NULL, abstol = 0.001, maxit = 1000, 
  verbose = FALSE, plot = FALSE, lmer.control = lmerControl(), \dots)
}

\arguments{
  \item{formula}{Model formula for building the \code{lmtree}.
    See examples below, or \code{\link[partykit]{lmtree}} documentation.} 
  \item{randomformula}{Formula for random-effects model.
    See examples below, or \code{\link[lme4]{lmer}} and
    \code{\link[lme4]{glmer}} documentation.}
  \item{data}{Dataset to be used for estimating the lmertree.}
  \item{ranefstart}{A vector of length \code{nrow(data)}, to be used as
    an offset in estimation of the first lmtree. \code{NULL} by default,
    which results in a zero offset initialization.}
  \item{abstol}{The convergence criterion used for estimation of the \code{lmertree}.
    When the difference in log likelihoods of the random-effects model of two
    consecutive iterations is smaller than \code{abstol}, estimation of the
    \code{lmertree} has converged.} 
  \item{maxit}{The maximum number of iterations to be performed in executing
    the \code{lmertree} function.}
  \item{verbose}{Should the log likelihood value of the estimated
    random-effects model be printed for every iteration during execution of
    the \code{lmertree} function?} 
  \item{plot}{Should the lmertree be plotted at every iteration during
    execution of the lmertree function? Note that plotting of a tree at
    every iteration slows down execution of the function.}
  \item{lmer.control}{An optional list (of correct class, resulting from
    \code{lmerControl()}), containing control parameters to be passed to
    \code{lmer()}. See \code{\link[lme4]{lmerControl}} documentation for
    details.} 
  \item{\dots}{Additional arguments to be passed to \code{lmtree()}.
    See \code{\link[partykit]{mob_control}} documentation for details.}
}
  
\details{
  This code is the first implementation and might change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{tree}{The final \code{lmtree}.}
  \item{lmer}{The final \code{lmer} random-effects model.}
  \item{ranef}{The corresponding random effects of \code{lmer}.} 
  \item{varcorr}{The corresponding \code{VarCorr(lmer)}.}
  \item{variance}{The corresponding \code{attr(VarCorr(lmer), "sc")^2}.}
  \item{data}{The dataset specified with the \code{data} argument
    including added auxiliary variables \code{.ranef} and \code{.treeresponse}
    from the last iteration.}
  \item{loglik}{The log-likelihood value of the last iteration.}
  \item{iterations}{The number of iterations used to estimate the \code{lmertree}.} 
  \item{maxit}{The maximum number of iterations specified with the \code{maxit} argument.}
  \item{ranefstart}{The random effects used as an offset, as specified with
    the \code{ranefstart} argument.}
  \item{formula}{The formula as specified with the \code{formula} argument.}
  \item{randomformula}{The formula as specified with the \code{randomformula} argument.}
  \item{abstol}{The prespecified value for the change in log-likelihood to evaluate
    convergence, as specified with the \code{abstol} argument.} 
  \item{mob.control}{A list containing control parameters passed to
    \code{lmtree()}, as specified with \dots.}
  \item{lmer.control}{A list containing control parameters passed to
    \code{lmer()}, as specified in the \code{lmer.control} argument.}
}

\seealso{\code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\link[partykit]{lmtree}}, \code{\link[partykit]{glmtree}}}

\examples{
## lmertree example
data("TreatmentSubgroups", package = "glmertree")

## fit normal linear regression LMM tree
lt <- lmertree(ynum ~ treatment | u1 + u2 + u3 + u4 + u5,
  ynum ~ 1 | cluster, data = TreatmentSubgroups)

## results
lt
plot(lt)
coef(lt)
ranef(lt)
}

\keyword{tree}
