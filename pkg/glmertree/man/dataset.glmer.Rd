\name{dataset.glmer}
\alias{dataset.glmer}
\docType{data}

\title{Example dataset for glmertree function}

\description{The dataset contains values on 8 variables for 1000 simulated observations.}

\usage{data("dataset.glmer")}

\format{
  The format is:
  chr "dataset.glmer"
}

\details{
  The dataset contains values on 8 variables for 1000 simulated observations. Variables are:
  U1 through U5: continuous partitioning covariates 
  treatment: a binary treatment indicator
  cluster: an indicator for cluster with 10 levels
  Y: a binary treatment outcome
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\references{
%%  ~~ possibly secondary sources and usages ~~
}

\examples{
  data(dataset.glmer)
  str(dataset.glmer)
}

\keyword{datasets}
