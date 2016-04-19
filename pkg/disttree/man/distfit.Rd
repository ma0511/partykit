\name{distfit}
\alias{distfit}
\alias{bread.distfit}
\alias{coef.distfit}
\alias{confint.distfit}
\alias{estfun.distfit}
\alias{logLik.distfit}
\alias{nobs.distfit}
\alias{vcov.distfit}

\title{Maximum-Likelihood Fitting of Parametric Distributions}

\description{
  The function \code{distfit} carries out maximum-likelihood
  estimation of parameters for distributions from the GAMLSS
  family (for generalized additive models for location, scale,
  and shape). The parameters can be transformed through
  link functions but do not depend on further covariates
  (i.e., are constant across observations).
}

\usage{
  distfit(y, family, weights = NULL, start = NULL,
    vcov. = TRUE, estfun = TRUE, \dots)
}

\arguments{
  \item{y}{numeric vector.}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object or
    a function generating such a family.}
  \item{weights}{optional numeric vector of case weights.}
  \item{start}{FIXME}
  \item{vcov.}{logical or character. FIXME}
  \item{estfun}{logical. Should the matrix of observation-wise
    score contributions (or empirical estimating functions) be
    returned?}
  \item{\dots}{further arguments passed to \code{control} for
    \code{\link[stats]{optim}}.}
}

\details{
  The function \code{distfit} fits distributions,
  similar to \code{fitdistr} from MASS (Venables and Ripley 2002)
  but based on GAMLSS families (Stasinopoulos and Rigby 2007).
  
  Provides analytical gradients and hessian, can be plugged into
  \code{\link[partykit]{mob}}.

  The resulting object of class \code{"distfit"} comes with a set of 
  standard methods to generic functions including \code{coef}, FIXME.
}

\value{
  \code{distfit} returns an object of class \code{"distfit"} which is a list with
  the following components:
  \item{y}{y,}
  \item{weights}{weights,}
  \item{family}{family,}
  \item{startpar}{startpar,}
  \item{starteta}{starteta,}
  \item{opt}{opt,}
  \item{par}{par,}
  \item{eta}{eta,}
  \item{hess}{hessian,}
  \item{call}{cl,}
  \item{ny}{ny,}	
  \item{nobs}{nobs,}	
  \item{vcov}{vc,}
  \item{estfun}{ef,}
  \item{ddist}{ddist,}
  \item{pdist}{pdist,}
  \item{qdist}{qdist,}
  \item{rdist}{rdist.}
}

\references{
Stasinopoulos DM, Rigby RA (2007). 
  Generalized Additive Models for Location Scale and Shape (GAMLSS) in R,
  \emph{Journal of Statistical Software}, \bold{23}(7), 1-46.
  \doi{10.18637/jss.v023.i07}

Venables WN, Ripley BD (2002).
  \emph{Modern Applied Statistics with S}.
  4th Edition. Springer-Verlag, New York.
}

\seealso{\code{\link[gamlss.dist]{gamlss.family}}, \code{\link[stats]{optim}}}

\examples{
## simulate artifical negative binomial data
set.seed(0)
y <- rnbinom(1000, size = 1, mu = 2)

## simple distfit
df <- distfit(y, family = NBI)
coef(df)
confint(df)
logLik(df)

## using tabulated data
ytab <- table(y)
df2 <- distfit(as.numeric(names(ytab)), family = NBI, weights = ytab)
coef(df2)
confint(df2)
logLik(df2)

## coefficients tests
if(require("lmtest")) {
  coeftest(df)
  coeftest(df2)
}
}

\keyword{distribution}
