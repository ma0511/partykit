\name{lmtree}
\alias{lmtree}

\alias{plot.lmtree}
\alias{predict.lmtree}
\alias{print.lmtree}

\title{Linear Model Trees}

\description{
  Model-based recursive partitioning based on least squares
  regression.
}

\usage{
lmtree(formula, data, subset, na.action, weights, offset, \dots)
}

\arguments{
  \item{formula}{symbolic description of the model (of type
    \code{y ~ z1 + \dots + zl} or \code{y ~ x1 + \dots + xk | z1 + \dots + zl};
    for details see below).}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights. By default these are
    treated as case weights but the default can be changed in
    \code{\link{mob_control}}.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the model \code{y ~ x1 + \dots + xk} (i.e., only when
    \code{x} variables are specified).}
  \item{\dots}{optional control parameters passed to
    \code{\link{mob_control}}.}
}

\details{
Convenience interface for fitting MOBs (model-based recursive partitions) via
the \code{\link{mob}} function. \code{lmtree} internally sets up a model
\code{fit} function for \code{mob}, using either \code{\link[stats]{lm.fit}}
or \code{\link[stats]{lm.wfit}} (depending on whether weights are used or not).
Then \code{mob} is called using the residual sum of squares as the objective
function.

The implementation makes tries to be avoid making unnecessary computations
while growing the tree. Also, it provides a more elaborate plotting function.
}

\value{
  An object of class \code{lmtree} inheriting from \code{\link{modelparty}}.
  The \code{info} element of the overall \code{party} and the individual
  \code{node}s contain various informations about the models.
}

\references{ 
  Achim Zeileis, Torsten Hothorn, and Kurt Hornik (2008). Model-Based
  Recursive Partitioning. \emph{Journal of Computational and Graphical Statistics}, 
  \bold{17}(2), 492--514.
}

\seealso{\code{\link{mob}}, \code{\link{mob_control}}, \code{\link{glmtree}}}

\examples{
if(require("mlbench")) {

## Boston housing data
data("BostonHousing", package = "mlbench")
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)

## linear model tree
mb_bh <- lmtree(medv ~ log(lstat) + I(rm^2) | zn +
  indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  data = BostonHousing, minsplit = 40)

## printing whole tree or individual nodes
print(mb_bh)
print(mb_bh, node = 7)

## plotting
plot(mb_bh)
plot(mb_bh, tp_args = list(which = "log(lstat)"))
plot(mb_bh, terminal_panel = NULL)

## estimated parameters
coef(mb_bh)
coef(mb_bh, node = 9)
summary(mb_bh, node = 9)

## various ways for computing the mean squared error (on the training data)
mean((BostonHousing$medv - fitted(mb_bh))^2)
mean(residuals(mb_bh)^2)
deviance(mb_bh)/sum(weights(mb_bh))
deviance(mb_bh)/nobs(mb_bh)

## log-likelihood and information criteria
logLik(mb_bh)
AIC(mb_bh)
BIC(mb_bh)
## (Note that this penalizes estimation of error variances, which
## were treated as nuisance parameters in the fitting process.)

## different types of predictions
bh <- BostonHousing[c(1, 10, 50), ]
predict(mb_bh, newdata = bh, type = "node")
predict(mb_bh, newdata = bh, type = "response")
predict(mb_bh, newdata = bh, type = function(object) summary(object)$r.squared)

}
}
\keyword{tree}
