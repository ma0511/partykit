\name{evtree}
\alias{evtree}
\alias{evtree-package}
\title{
Evolutionary Learning of Globally Optimal Trees
}
\description{
	evtree: learning of globally optimal classification and regression trees by using evolutionary algorithms.	
}
\usage{evtree(formula, data=list(), weights=NULL, subset = NULL, 
control = evtree.control(...), ...)}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit.}
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
Globally optimal classification and regression trees are learned by using evolutionary algorithm. Roughly, the algorithm works as follows. First, a set of trees is initialized with random split rules in the root nodes. Second, mutation and crossover operators are applied to modify the trees' structure and the tests that are applied in the internal nodes. After each modification step a survivor selection mechanism selects the best candidate models for the next iteration. In this evolutionary process the mean quality of the population increases over time. The algorithm terminates when the quality of the best trees does not improve further, but not later than a maximum number of iterations specified by \code{niterations} in \code{\link{evtree.control}}.


More details on the algorithm are provided by \code{vignette("evtree")}. 


The resulting trees can be summarized and visualized by \code{print()}, \code{summary()}, and \code{plot()} functions provided by the \code{partykit} package. A \code{predict()} function is utilized to predict the class labels of new data. 
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

