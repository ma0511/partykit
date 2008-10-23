\name{party-methods}
\alias{party-methods}
\alias{length.party}
\alias{[.party}
\alias{[[.party}
\alias{depth.party}
\alias{width.party}
\title{ Methods for Party Objects}
\description{

  Methods for computing on \code{party} objects.

}
\usage{
  \method{length}{party}(x)
  \method{[}{party}(x, i, ...)
  \method{[[}{party}(x, i, ...)
  \method{depth}{party}(x, ...)
  \method{width}{party}(x, ...)
}
\arguments{
  \item{x}{ an object of class \code{\link{party}}.}
  \item{i}{ an integer specifying the root of the subtree to extract.}
  \item{\dots}{ additional arguments.}
}
\details{

  \code{length} gives the number of nodes in the tree (in contrast
  to the \code{length} method for \code{\link{partynode}} objects),
  \code{depth} the depth of the tree and \code{width} 
  the number of terminal nodes.

}
\examples{

  ## a tree as flat list structure
  nodelist <- list(
    # root node
    list(id = 1L, split = partysplit(varid = 4L, breaks = 1.9), 
         kids = 2:3),
    # V4 <= 1.9, terminal node
    list(id = 2L),
    # V4 > 1.9
    list(id = 3L, split = partysplit(varid = 5L, breaks = 1.7), 
         kids = c(4L, 7L)),
    # V5 <= 1.7
    list(id = 4L, split = partysplit(varid = 4L, breaks = 4.8), 
         kids = 5:6),
    # V4 <= 4.8, terminal node
    list(id = 5L),
    # V4 > 4.8, terminal node
    list(id = 6L),
    # V5 > 1.7, terminal node
    list(id = 7L)
  )

  ## convert to a recursive structure
  node <- as.partynode(nodelist)

  ## set up party object
  data("iris")
  tree <- party(node, data = iris, 
                fitted = data.frame("(fitted)" = fitted_node(node, data = iris),
                                    check.names = FALSE),
)
  ## number of kids in root node
  length(tree)

  ## depth of tree
  depth(tree)

  ## number of terminal nodes
  width(tree)

  tree[2]

}
\keyword{tree}
