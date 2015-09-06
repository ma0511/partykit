\name{TreatmentSubgroups}
\alias{TreatmentSubgroups}

\title{Artificial Data for Illustrating (G)LMM Trees}

\description{
Artificial dataset for illustrating (G)LMM trees for detecting
treatment-subgroup interactions.
}

\usage{
data("TreatmentSubgroups")
}

\format{
A data frame containing 1000 observations on 8 variables.
\describe{
  \item{u1,u2,u3,u4,u5}{numeric. Continuous partitioning variables.}
  \item{treatment}{factor. Binary treatment variable.}
  \item{cluster}{factor. Indicator for cluster with 10 levels.}
  \item{ynum}{numeric. Continuous treatment outcome.}
  \item{ybin}{factor. Binary treatment outcome.}
}
}

\details{
The data contains four treatment effect subgroups with respect to the
continuous partitioning variables (\code{u1}--\code{u5}). The \code{cluster}
variable has an additional random intercept that should be accounted for.
The outcome is assessed either by a continuous variable (\code{ynum}) or
a binary variable (\code{ybin}).

There are two large subgroups where Treatment 1 is better than Treatment 2 and
vice versa, respectively. Additionally, there are two smaller subgroups where
both treatments lead to comparable outcomes. For the corresponding (G)LMM trees
see \code{\link{lmertree}}.
}

\seealso{\code{\link{lmertree}}}

\examples{
data("TreatmentSubgroups", package = "glmertree")
summary(TreatmentSubgroups)
}

\keyword{datasets}
