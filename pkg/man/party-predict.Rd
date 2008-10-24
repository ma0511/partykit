\name{party-predict}
\alias{party-predict}
\alias{predict.party}
\alias{predict_party}
\alias{predict_party.default}
\alias{predict_party.constparty}
\alias{predict_party.simpleparty}
\title{ Tree Predictions }
\description{

  Compute predictions from \code{party} objects.

}
\usage{
\method{predict}{party}(object, newdata = NULL, ...)
predict_party(party, id, newdata = NULL, ...)
\method{predict_party}{default}(party, id, newdata = NULL, ...) 
\method{predict_party}{constparty}(party, id, newdata = NULL,
    type = c("response", "prob", "node"), FUN = NULL, simplify = TRUE, ...)
\method{predict_party}{simpleparty}(party, id, newdata = NULL,
    type = c("response", "prob", "node"), ...)
}
\arguments{
  \item{object}{ objects of class \code{\link{party}}. }
  \item{newdata}{ an optional data frame in which to look for variables with
          which to predict,if omitted, the fitted values are used.}
  \item{party}{ objects of class \code{\link{party}}. }
  \item{id}{ a vector of terminal node identifiers. }
  \item{type}{}
  \item{FUN}{ a function to compute summary statistics, i.e., constant
              predictions for each node with argument list \code{(y, w)}
              where \code{y} is the response and \code{w} are case weights.}
  \item{simplify}{}
  \item{\dots}{ additional arguments. }
}
\details{

  The \code{\link{predict}} method for \code{\link{party}} objects
  computes the identifiers of the predicted terminal nodes, either
  for new data in \code{newdata} or for the learning samples
  (only possible for objects of class \code{constparty}).
  These identifiers are delegated to the corresponding
  \code{predict_party} method which computes (via 
  \code{FUN} for class \code{constparty})
  or extracts (class \code{simpleparty}) the actual predictions.

}
\value{
  A list of predictions, possibly simplified to a numeric vector or
  character.
}
\examples{

  ## fit tree using rpart
  library("rpart")
  rp <- rpart(skips ~ Opening + Solder + Mask + PadType + Panel,
              data = solder, method = 'anova')

  ## coerse to `constparty'
  pr <- as.party(rp)

  ## mean predictions
  predict(pr, newdata = solder[c(3, 541, 640),])

  ## terminal node identifiers
  predict(pr, newdata = solder[c(3, 541, 640),], type = "node")

  ## median predictions
  predict(pr, newdata = solder[c(3, 541, 640),], 
          FUN = function(y, w = 1) median(y))

}
\keyword{tree}
