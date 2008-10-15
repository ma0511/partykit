\name{partynode}
\alias{partynode}
\alias{kidids_node}
\alias{fitted_node}
\alias{id_node}
\alias{split_node}
\alias{surrogates_node}
\alias{kids_node}
\alias{info_node}
\title{ Nodes in Regression Trees }
\description{

  A class for representing recursive nodes in decision trees
  and for computing with nodes.

}
\usage{
  partynode(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL)
  kidids_node(node, data, vmatch = 1:ncol(data), obs = NULL)
  fitted_node(node, data, vmatch = 1:ncol(data), obs = 1:nrow(data))
  id_node(node)
  split_node(node)
  surrogates_node(node)
  kids_node(node)
  info_node(node)
}
\arguments{
  \item{id}{ a unique integer identifier for a node. }
  \item{split}{ an object of class \code{\link{partysplit}}. }
  \item{kids}{ a list of \code{node} objects. }
  \item{surrogates}{ a list of \code{partysplit} objects.}
  \item{info}{ additional information. }
  \item{node}{ an object of class \code{partynode}.}
  \item{data}{ a \code{\link{list}} or \code{\link{data.frame}}.}
  \item{vmatch}{ a permutation of the variable numbers in \code{data}.}
  \item{obs}{ a logical or integer vector indicating a subset of the
              observations in \code{data}.}
}
\details{

  A node represents a recursive regression tree and
  consists of an identifier. A node consisting only
  of an identifier (and possibly additional information
  in \code{info}) is a terminal node.

  Inner nodes consist of a primary split 
  (an object of class \code{\link{partysplit}})
  and at least two kids. Kid nodes are objects of class \code{node}
  itself, so the structure is defined recursively.
  In addition, a list of \code{partysplit} objects 
  offering surrogate splits can be supplied. Like
  \code{\help{partysplit}} objects, \code{partynode}
  objects aren't connected to the actual data.

  Function \code{kidids_node()} determines how
  the observation in \code{data[obs,]} are assigned 
  to the kid nodes (and returns the kid node number).
  This is done by evaluating \code{split} (and possibly
  all surrogate splits) on \code{data}.

  Function \code{fitted_node()} performs all 
  splits recursively and returns the \code{id}
  of the terminal node each observation in
  \code{data[obs,]} belongs to.

  All other functions are accessor function for
  extracting information from objects of class \code{partynode}.

}
\value{
  An object of class \code{partynode}:
  \item{id}{ a unique integer identifier for a node. }
  \item{split}{ an object of class \code{\link{partysplit}}. }
  \item{kids}{ a list of \code{node} objects. }
  \item{surrogates}{ a list of \code{partysplit} objects.}
  \item{info}{ additional information. }
}
\examples{

  data("iris")

  ## a stump defined by a binary split in Sepal.Length
  stump <- partynode(id = 1L, 
      split = partysplit(which(names(iris) == "Sepal.Length"),
                         breaks = 5),
      kids = lapply(2:3, partynode))

  ## textual representation
  print(stump, data = iris)  

  ## number and node id of the two terminal nodes
  table(kidids_node(stump, iris), 
        fitted_node(stump, data = iris))

}
\keyword{tree}
