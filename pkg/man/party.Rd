\name{party}
\alias{party}
\alias{names.party}
\alias{names<-.party}
\alias{node_party}
\title{ Recursive Partytioning }
\description{

  A class for representing decision trees
  and accessor functions.

}
\usage{
party(node, data, fitted = NULL, terms = NULL, names = NULL, info = NULL)
\method{names}{party}(x)
\method{names}{party}(x) <- value
node_party(party)
}
\arguments{
  \item{node}{ an object of class \code{\link{partynode}}.}
  \item{data}{ a (potentially empty) \code{\link{data.frame}}.}
  \item{fitted}{ an optional \code{\link{data.frame}} with \code{nrow(data)} 
                 rows (only if \code{nrow(data) != 0} and 
                 containing at least the fitted terminal node ids
                 as element \code{(fitted)}.}
  \item{terms}{ an optional \code{\link{terms}} object.  }
  \item{names}{ an optional vector of names to be assigned to each node of \code{node}. }
  \item{info}{ additional information. }
  \item{x}{ an object of class \code{party}.}
  \item{party}{ an object of class \code{party}.}
  \item{value}{a character vector of up to the same length as \code{x}, or
    \code{NULL}.}
}
\details{

}
\value{
  The constrcutor returns an object of class \code{party}:
  \item{node}{ an object of class \code{\link{partynode}}.}
  \item{data}{ a (potentially empty) \code{\link{data.frame}}.}
  \item{fitted}{ an optional \code{\link{data.frame}} with \code{nrow(data)}
                 rows (only if \code{nrow(data) != 0} and
                 containing at least the fitted terminal node ids
                 as element \code{(fitted)}.}
  \item{terms}{ an optional \code{\link{terms}} object.  }
  \item{names}{ an optional vector of names to be assigned to each node of \code{node}. }
  \item{info}{ additional information. }
  \code{names} can be used to set and retrieve names of nodes and 
  \code{node_party} returns an object of class \code{\link{partynode}}.
}
\examples{

  data("iris")

  ## a stump defined by a binary split in Sepal.Length
  stump <- partynode(id = 1L,
      split = partysplit(which(names(iris) == "Sepal.Length"),
                         breaks = 5),
      kids = lapply(2:3, partynode))

  party(stump, iris, 
        fitted = data.frame("(fitted)" = fitted_node(stump, data = iris), 
                            check.names = FALSE),
        names = c("root", "left", "right"))

}
\keyword{tree}
