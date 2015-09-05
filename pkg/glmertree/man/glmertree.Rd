\name{glmertree}
\alias{glmertree}
\alias{print.glmertree}
\alias{plot.glmertree}
\alias{coef.glmertree}
\alias{ranef.glmertree}

\title{Generalized Linear Mixed Model Trees}

\description{
  Model-based recursive partitioning based on generalized linear mixed models.
}

\usage{
glmertree(formula, randomformula, data,
  family = "binomial", ranefstart = NULL, 
  abstol = 0.001, maxit = 1000, verbose = FALSE, plot = FALSE,
  glmer.control = glmerControl(), \dots)
}

\arguments{
  \item{formula}{Model formula for building the \code{glmtree}. See
    examples below, or \code{\link[partykit]{glmtree}} documentaton.}
  \item{randomformula}{Formula for random-effects model. See examples
    below, or \code{lmer} and \code{\link[lme4]{glmer}} documentation.}
  \item{data}{Dataset to be used for estimating the \code{glmertree}.}
  \item{family}{Specification of a family for \code{glmtree} and \code{glmer}.
    See \code{\link[stats]{glm}} documentation for families.}
  \item{ranefstart}{A vector of length \code{nrow(data)}, to be used as an
    offset in estimation of the first \code{glmtree}. \code{NULL} by default,
    which results in a zero offset being used.}
  \item{abstol}{The convergence criterion used for estimation of the
    \code{glmertree}. When the difference in log-likelihoods of the
    random-effects model of two consecutive iterations is smaller than
    \code{abstol}, estimation of the \code{glmertree} has converged.} 
  \item{maxit}{The maximum number of iterations to be performed in
    executing the \code{glmertree} function.}
  \item{verbose}{Should the log-likelihood value of the estimated
    random-effects model be printed for every iteration during execution
    of the \code{glmertree} function?}
  \item{plot}{Should the \code{glmertree} be plotted at every iteration during
    execution of the glmertree function? Note that plotting of a tree at
    every iteration slows down execution of the function.}
  \item{glmer.control}{An optional list (of correct class, resulting from
    \code{glmerControl()}), containing control parameters to be
    passed to \code{glmer()}. See \code{\link[lme4]{glmerControl}} documentation
    for details.}
  \item{\dots}{Additional arguments to be passed to \code{glmtree()}.
    See \code{\link[partykit]{mob_control}} documentation for details.}
}
  

\details{
  This code is the first implementation and might change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{tree}{The final \code{glmtree}.}
  \item{glmer}{The final \code{glmer} random-effects model.}
  \item{ranef}{The corresponding random effects of \code{glmer}.} 
  \item{varcorr}{The corresponding \code{VarCorr(glmer)}.}
  \item{variance}{The corresponding \code{attr(VarCorr(glmer), "sc")^2}.}
  \item{data}{The dataset specified with the \code{data} argument
    including added auxiliary variables \code{.ranef} and \code{.treeresponse}
    from the last iteration.}
  \item{loglik}{The log-likelihood value of the last iteration.}
  \item{iterations}{The number of iterations used to estimate the \code{glmertree}.} 
  \item{maxit}{The maximum number of iterations specified with the \code{maxit} argument.}
  \item{ranefstart}{The random effects used as an offset, as specified with
    the \code{ranefstart} argument.}
  \item{formula}{The formula as specified with the \code{formula} argument.}
  \item{randomformula}{The formula as specified with the \code{randomformula} argument.}
  \item{abstol}{The prespecified value for the change in log-likelihood to evaluate
    convergence, as specified with the \code{abstol} argument.} 
  \item{mob.control}{A list containing control parameters passed to
    \code{glmtree()}, as specified with \dots.}
  \item{glmer.control}{A list containing control parameters passed to
    \code{glmer()}, as specified in the \code{glmer.control} argument.}
}

\seealso{\code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\link[partykit]{lmtree}}, \code{\link[partykit]{glmtree}}}

\examples{
## glmertree example (binary outcome)
data("TreatmentSubgroups", package = "glmertree")

## fit logistic regression GLMM tree
gt <- glmertree(ybin ~ treatment | u1 + u2 + u3 + u4 + u5,
  ybin ~ 1 | cluster, data = TreatmentSubgroups)

## results
gt
plot(gt)
coef(gt)
ranef(gt)
}

\keyword{tree}
