\name{evtree}
\alias{evtree}
\alias{evtree-package}
\title{
Evolutionary Learning of Globally Optimal Trees
}
\description{
	evtree: learning of globally optimal classiffication and regression trees by using evolutionary algorithms	
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
  \item{...}{ arguments to \code{evtree.control} may also be specified in the
          call to \code{evtree}.}
  }

\details{
Commonly used classification and regression tree methods like the CART algorithm are recursive partitioning methods that build the model in a forward stepwise search. Although this approach is known to be an efficient heuristic, the results of recursive tree methods are only locally optimal, as splits are chosen to maximize homogeneity at the next step only. An alternative way to search over the parameter space of trees is to use global optimization methods like evolutionary algorithms. The evtree package implements an evolutionary algorithm for learning globally optimal classification and regression trees in R. CPU and memory-intensive tasks are fully computed in C++ while the partykit package is leveraged to represent the resulting trees in R, providing unified infrastructure for summaries, visualizations, and predictions.
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

