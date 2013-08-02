\name{mob_control}
\alias{mob_control}

\title{Control Parameters for Model-Based Partitioning}

\description{
  Various parameters that control aspects the fitting algorithm
  for recursively partitioned \code{\link{mob}} models.
}

\usage{
mob_control(alpha = 0.05, bonferroni = TRUE, minsplit = NULL, trim = 0.1,
  breakties = FALSE, parm = NULL, verbose = FALSE, caseweights = TRUE,
  ytype = "vector", xtype = "matrix", terminal = "object", inner = terminal,
  model = TRUE)
}

\arguments{
  \item{alpha}{numeric significance level. A node is splitted when
    the (possibly Bonferroni-corrected) \eqn{p} value for any parameter
    stability test in that node falls below \code{alpha}.}
  \item{bonferroni}{logical. Should \eqn{p} values be Bonferroni corrected?}
  \item{minsplit}{integer. The minimum number of observations in a node.
    If \code{NULL}, the default is to use 10 times the number of parameters
    to be estimated (divided by the number of responses per observation
    if that is greater than 1).}
  \item{trim}{numeric. This specifies the trimming in the parameter instability
    test for the numerical variables. If smaller than 1, it is interpreted
    as the fraction relative to the current node size.}
  \item{breakties}{logical. Should ties in numeric variables be broken
    randomly for computing the associated parameter instability test?}
  \item{parm}{numeric or character. Number or name of model parameters
    included in the parameter instability tests (by default all parameters
    are included).}
  \item{verbose}{logical. Should information about the fitting process
    of \code{\link{mob}} (such as test statistics, \eqn{p} values, selected
    splitting variables and split points) be printed to the screen?}
  \item{caseweights}{logical. Should weights be interpreted as case weights?
    If \code{TRUE}, the number of observations is \code{sum(weights)},
    otherwise it is \code{sum(weights > 0)}.}
  \item{ytype, xtype}{character. Specification of how \code{mob} should
    preprocess \code{y} and \code{x} variables. Possible choice are:
    \code{"vector"} (for \code{y} only), i.e., only one variable;
    \code{"matrix"}, i.e., the model matrix of all variables;
    \code{"data.frame"}, i.e., a data frame of all variables.}
  \item{terminal, inner}{character. Specification of which additional
    information (\code{"estfun"}, \code{"object"}, or both) should be
    stored in each node. If \code{NULL}, no additional information is
    stored.}
  \item{model}{logical. Should the full model frame be stored in the
    resulting object?}
}

\details{
  See \code{\link{mob}} for more details and references.
}

\seealso{\code{\link{mob}}}

\value{
  A list of class \code{mob_control} containing the control parameters.
}

\keyword{misc}
