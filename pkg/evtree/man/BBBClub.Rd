\name{BBBClub}
\alias{BBBClub}

\title{Bookbinder's Book Club}

\description{
  Marketing case study about a (fictitious) American book club to whose
  customers a book about \dQuote{The Art History of Florence} was advertised.
}

\usage{data("BBBClub")}

\format{
A data frame containing 1,300 observations on 11 variables.
\describe{
  \item{choice}{factor. Did the customer buy the advertised book?}
  \item{gender}{factor indicating gender.}  
  \item{amount}{total amount of money spent at the BBB Club.}  
  \item{freq}{number of books purchased at the BBB Club.}  
  \item{last}{number of months since the last purchase.}  
  \item{first}{number of months since the first purchase.}  
  \item{child}{number of children's books purchased.}  
  \item{youth}{number of youth books purchased.}  
  \item{cook}{number of cookbooks purchased.}  
  \item{diy}{number of do-it-yourself books purchased.}  
  \item{art}{number of art books purchased.}  
}
}

\details{
  The data is a marketing case study about a (fictitious) American book club. 
  taken from the \emph{Marketing Engineering} textbook of Lilien and Rangaswamy (2004).  
  In this case study, a brochure of the book \dQuote{The Art History of Florence} was
  sent to 20,000 customers and 1,806 of which bought the book. A subsample of 1,300
  customers is provided in \code{BBBClub} for building a predictive model for \code{choice}.

  The use of a cost matrix is suggested for this dataset.	
  Classifying a customer that purchased the book as a non-buyer is worse (cost = 5),
  than it is to classify a custumer that did not purchase the book as a buyer (cost = 1). 
}

\source{
  Complements to Lilien and Rangaswamy (2004).
}

\references{
  Lilien GL, Rangaswamy A (2004).
    \emph{Marketing Engineering: Computer-Assisted Marketing Analysis and Planning},
    2nd edition. Victoria, BC: Trafford Publishing. 
}

\examples{
data("BBBClub")
summary(BBBClub)
ev <- evtree(choice ~ ., BBBClub, control = evtree.control(maxdepth = 3, 
minbucket = 65, seed = 1000))
ev
plot(ev)
table(predict(ev), BBBClub$choice)
# the use of a different seeds can result in different solutions.  
ev2 <- evtree(choice ~ ., BBBClub, control = evtree.control(maxdepth = 3, 
minbucket = 65, seed = 2000))
ev2
plot(ev2)
table(predict(ev2), BBBClub$choice)
#with case weigths
gcw <- array(1, nrow(BBBClub))
gcw[BBBClub$choice == "yes"] <- 5
ev <- evtree(choice ~ ., BBBClub, control = evtree.control(maxdepth = 3, 
minbucket = ceiling(sum(gcw) * 0.05), seed = 100), weights = gcw)
ev
plot(ev)
table(predict(ev), BBBClub$choice)
# the use of a different seeds can result in different solutions.  
ev2 <- evtree(choice ~ ., BBBClub, control = evtree.control(maxdepth = 3, 
minbucket = ceiling(sum(gcw) * 0.05), seed = 200), weights = gcw)
ev2
plot(ev2)
table(predict(ev2), BBBClub$choice)
}

\keyword{datasets}
