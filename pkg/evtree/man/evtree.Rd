\name{evtree}
\alias{evtree}
\alias{evtree-package}
\title{
evtree
}
\description{
 	evtree: Evolutionary Learning of Globally Optimal Classification and Regression Trees in R
}
\usage{evtree(formula, data=list(), weights=NULL, subset = NULL, control = evtree.control(...), ...)}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. Only non-negative integer valued weights are
                  allowed.}
  \item{control}{ a list with control parameters, see
                 \code{\link{evtree.control}}.}
}

\details{
Commonly used classification and regression tree methods like the CART algorithm are recursive methods that build the model in a forward stepwise search. Although this approach is known to be an efficient heuristic, the results of recursive tree methods are only locally optimal, as splits are chosen to maximize homogeneity at the next step only. One way to build trees that are more globally optimal is to use stochastic optimization methods, like evolutionary algorithms. This method implements and evolutionary algorithm to build classification and regression trees. 

First, a set of trees is initialized with random split rules in the root nodes. Secondly, the trees are iteratively modified by mutation and crossover operators. After each modification a survivor selection mechanism selects the best candidate models for the next iteration. The algorithm terminates when the quality of the best trees does not improve further, but not later than a maximum number of iterations specified by \code{niterations} in \code{\link{evtree.control}}. 
}


\value{
  An object of class \code{\link{party}}.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Thomas Grubinger, Achim Zeileis, Karl-Peter Pfeiffer
}

\examples{
    ### regression
    airq <- subset(airquality, !is.na(Ozone) & complete.cases(airquality))
    airct <- evtree(Ozone ~ ., data = airq)
    airct
    plot(airct)
    mean((airq$Ozone - predict(airct))^2)

    ### classification
    irisct <- evtree(Species ~ .,data = iris)
    irisct
    plot(irisct)
    table(predict(irisct), iris$Species)
}

\keyword{tree}

