\name{nodeapply}
\alias{nodeapply}
\alias{nodeapply.party}
\alias{nodeapply.partynode}
\title{ Apply Functions Over Nodes }
\description{

  Returns a list of values obtained by applying a function to
  \code{partynode} objects.

}
\usage{
nodeapply(obj, ids = 1, FUN = NULL, ...)
\method{nodeapply}{partynode}(obj, ids = 1, FUN = NULL, ...)
\method{nodeapply}{party}(obj, ids = 1, FUN = NULL, by_node = TRUE, ...)
}
\arguments{
  \item{obj}{ an object of class \code{\link{partynode}} 
             or \code{\link{party}}.}
  \item{ids}{ integer vector of node identifiers to apply over.}
  \item{FUN}{ a function to be applied to nodes. By default, the 
              node itself is returned.}
  \item{by_node}{ ??? }
  \item{\dots}{ additional arguments.}
}
\details{

  The method for \code{party} essentially calls the \code{nodeapply}
  method on it's \code{node} slot.

}
\value{
  A list.
}
\examples{

  ## a tree as flat list structure
  nodelist <- list(
    # root node
    list(id = 1L, split = partysplit(varid = 4L, breaks = 1.9),
         kids = 2:3),
    # V4 <= 1.9, terminal node
    list(id = 2L, info = "terminal A"),
    # V4 > 1.9
    list(id = 3L, split = partysplit(varid = 5L, breaks = 1.7),
         kids = c(4L, 7L)),
    # V5 <= 1.7
    list(id = 4L, split = partysplit(varid = 4L, breaks = 4.8),
         kids = 5:6),
    # V4 <= 4.8, terminal node
    list(id = 5L, info = "terminal B"),
    # V4 > 4.8, terminal node
    list(id = 6L, info = "terminal C"),
    # V5 > 1.7, terminal node
    list(id = 7L, info = "terminal D")
  )

  ## convert to a recursive structure
  node <- as.partynode(nodelist)

  ## return root node 
  nodeapply(node)

  ## return info slots of terminal nodes
  nodeapply(node, ids = nodeids(node, terminal = TRUE), 
            FUN = function(x) info_node(x))

}
\keyword{tree}
