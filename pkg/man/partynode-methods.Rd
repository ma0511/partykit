\name{partynode-methods}
\alias{partynode-methods}
\alias{is.partynode}
\alias{as.partynode}
\alias{as.partynode.partynode}
\alias{as.partynode.list}
\alias{as.list.partynode}
\alias{length.partynode}
\alias{[.partynode}
\alias{[[.partynode}
\alias{is.terminal}
\alias{is.terminal.partynode}
\alias{depth}
\alias{depth.partynode}
\alias{width}
\alias{width.partynode}
\title{ Methods for Node Objects}
\description{

  Methods for computing on \code{partynode} objects.

}
\usage{
  is.partynode(x)
  as.partynode(x, ...)
  \method{as.partynode}{partynode}(x, from = NULL, ...)
  \method{as.partynode}{list}(x, ...)
  \method{as.list}{partynode}(x, ...)
  \method{length}{partynode}(x)
  \method{[}{partynode}(x, i, ...)
  \method{[[}{partynode}(x, i, ...)
  is.terminal(x, ...)
  \method{is.terminal}{partynode}(x, ...)
  depth(x, ...)
  \method{depth}{partynode}(x, ...)
  width(x, ...)
  \method{width}{partynode}(x, ...)
}
\arguments{
  \item{x}{ an object of class \code{partynode} or \code{list}.}
  \item{from}{ an integer giving the identifier of the root node.}
  \item{i}{ an integer specifying the kid to extract.}
  \item{\dots}{ additional arguments.}
}
\details{

  \code{is.partynode} checks if the argument is a valid \code{partynode}
  objects. \code{is.terminal} is \code{TRUE} for terminal nodes
  and \code{FALSE} for inner nodes. The subset methods
  return the \code{partynode} object corresponding to the \code{i}th
  kid.

  The \code{as.partynode} and \code{as.list} methods can be used
  to convert flat list structures into recursive \code{partynode}
  objects and vice versa.

  \code{length} gives the number of kid nodes, \code{depth}
  the depth of the tree and \code{width} the number of terminal
  nodes.

}
\examples{

  ## a tree as flat list structure
  nodelist <- list(
    list(id = 1L, split = partysplit(varid = 4L, breaks = 1.9), 
         kids = 2:3),
    list(id = 2L),
    list(id = 3L, split = partysplit(varid = 5L, breaks = 1.7), 
         kids = c(4L, 7L)),
    list(id = 4L, split = partysplit(varid = 4L, breaks = 4.8), 
         kids = 5:6),
    list(id = 5L),
    list(id = 6L),
    list(id = 7L)
  )

  ## convert to a recursive structure
  node <- as.partynode(nodelist)

  ## print tree
  data("iris")
  print(node, data = iris)

  ## print subtree
  print(as.partynode(node[2]), data = iris)

  ## number of kids
  length(node)

  ## depth of tree
  depth(node)

  ## number of terminal nodes
  width(node)

  ## convert back to flat structure
  as.list(node)

}
\keyword{tree}
