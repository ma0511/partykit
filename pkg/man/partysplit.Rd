\name{partysplit}
\alias{partysplit}
\alias{kidids_split}
\alias{character_split}
\alias{varid_split}
\alias{breaks_split}
\alias{index_split}
\alias{right_split}
\alias{prob_split}
\alias{info_split}
\title{ Multiway Splits }
\description{

  A class for representing multiway splits in nominal or 
  ordered variables and functions for computing with splits.

}
\usage{
  partysplit(varid, breaks = NULL, index = NULL, right = TRUE, 
             prob = NULL, info = NULL)
  kidids_split(split, data, vmatch = 1:ncol(data), obs = NULL)
  character_split(split, data = NULL, digits = getOption("digits") - 2)
  varid_split(split)
  breaks_split(split)
  index_split(split)
  right_split(split)
  prob_split(split)
  info_split(split)
}
\arguments{
  \item{varid}{ an integer specifying the variable to split in, i.e., 
                a column in \code{data}. }
  \item{breaks}{ a numeric vector of split points. }
  \item{index}{ an integer vector containing a contiguous sequence 
                of 1 to number of kid nodes.}
  \item{right}{ a logical, indicating if the intervals defined by 
                \code{breaks} should be closed on the
                right (and open on the left) or vice versa.}
  \item{prob}{ a numeric vector representing a probability distribution 
               over kid nodes. }
  \item{info}{ additional information. }
  \item{split}{ an object of class \code{partysplit}.}
  \item{data}{ a \code{\link{list}} or \code{\link{data.frame}}.}
  \item{vmatch}{ a permutation of the variable numbers in \code{data}.}
  \item{obs}{ a logical or integer vector indicating a subset of the 
              observations in \code{data}.}
  \item{digits}{ minimal number of significant digits.}
}}
\details{

  A split is basically a function that maps data, 
  more specifically a partitioning variable,
  to a set of integers indicating the kid nodes. Objects of
  class \code{partysplit} describe such a function and can
  be set-up via the \code{partysplit()} constructor.
  The variables are available in a \code{list} or \code{data.frame} 
  (here called \code{data}) and \code{varid} specifies the 
  partitioning variable, i.e., the variable to split in. 
  The constructor \code{partysplit()} doesn't have access 
  to the actual data, i.e., doesn't \emph{estimate} splits.

  \code{kidids_split()} actually performs the split
  represented by its \code{split} argument on observations 
  \code{data} and assigns an integer (giving the
  kid node number) to each observation in \code{data[obs,]}.

  \code{character_split()} returns a character representation
  of its \code{split} argument. The remaining functions 
  defined here are accessor functions for \code{partysplit} 
  objects.

  The numeric vector \code{breaks} defines how the range of 
  the partitioning variable
  is divided into intervals (like in \code{\link{cut}}) and may be 
  \code{NULL}. The intervals are represented by integers 
  \code{1:length(breaks) + 1}.

  \code{index} assigns these \code{length(breaks) + 1} 
  intervals to one of at least two kid nodes.
  Kid nodes are represented by an integer identifiers (starting with 1).
  Either \code{breaks} or \code{index} must be given.

  \code{prob} defines a probability distribution over
  all kid nodes which is used for random splitting
  when a deterministic split isn't possible (due to missing 
  values, for example).

  \code{info} takes arbitrary information.

}
\value{
  An object of class \code{partysplit}:
  \item{varid}{ an integer specifying the variable to split in, i.e.,
                a column in \code{data}. }
  \item{breaks}{ a numeric vector of split points.}
  \item{index}{ an integer vector containing a contiguous sequence
                of 1 to number of kid nodes.}
  \item{right}{ a logical, indicating if the intervals defined by
                \code{breaks} should be closed on the
                right (and open on the left) or vice versa.}
  \item{prob}{ a numeric vector representing a distribution
               over kid nodes.}
  \item{info}{ arbitrary information.}
}
\seealso{\code{\link{cut}}}
\examples{

  data("iris")

  ## binary split in Sepal.Length
  sl5 <- partysplit(which(names(iris) == "Sepal.Length"),
                    breaks = 5)
  character_split(sl5, data = iris)
  table(kidids_split(sl5, data = iris), iris$Sepal.Length <= 5)

  ## multiway split in Sepal.Width, higher values
  ## go to the first kid, smallest values to the last kid
  sw23 <- partysplit(which(names(iris) == "Sepal.Width"),    
                     breaks = c(2, 3), index = 3:1)    
  character_split(sw23, data = iris)    
  table(kidids_split(sw23, data = iris), 
        cut(iris$Sepal.Width, breaks = c(-Inf, 2, 3, Inf)))   

  ## binary split in Species
  sp <- partysplit(which(names(iris) == "Species"),
                   index = c(1L, 1L, 2L))
  character_split(sp, data = iris)
  table(kidids_split(sp, data = iris), iris$Species)

}
\keyword{tree}
