\name{glmertree}
\alias{lmertree}
\alias{glmertree}

\title{Generalized Linear Mixed Model Trees}

\description{
  Model-based recursive partitioning based on (generalized) linear mixed models.
}

\usage{
lmertree(lmtreeformula, randomformula, data,
  initialRandomEffects = NULL, 
  ErrorTolerance = 0.001, MaxIterations = 1000, 
  verbose = TRUE, plotting = TRUE,
  lmer.control = lmerControl(), \dots)

glmertree(glmtreeformula, randomformula, data,
  family = "binomial",
  initialRandomEffects = NULL, 
  ErrorTolerance = 0.001, MaxIterations = 1000, 
  verbose = TRUE, plotting = TRUE,
  glmer.control = glmerControl(), \dots)
}

\arguments{
  \item{lmtreeformula, glmtreeformula, randomformula, data, initialRandomEffects,
    ErrorTolerance, MaxIterations, verbose, plotting, lmer.control, glmer.control,
    family, \dots}{various arguments that might change in future versions.}
}
  

\details{
This code is the first implementation and will probably change in future versions.
}

\value{
  A list with various components.
}

\seealso{\code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\link[partykit]{lmtree}}, \code{\link[partykit]{glmtree}}}

\examples{
## currently no example yet
1 + 1
}

\keyword{tree}
