\name{BBBClub}
\docType{data}
\alias{BBBClub}
\title{Bookbinder's Book Club}

\description{Data from a case study of an (fictitious) American book club. In this case study, a brochure of the book "The Art History of Florence" was sent to 20,000 customers. 1,806 of the brochure's addresses bought the book. 
The Bookbinder's book club collected several variables of its customers with the purpose of building a prognosis model for the customer's purchase decision. This data-set is a subset of the described case study, containing 1,300 observations on 11 Variables.
}

\usage{data("BBBClub")}

\format{
A data frame containing 1,300 observations on 11 variables.
\describe{
  \item{choice}{binary variable indicating the customers purchase decision, with levels \code{yes} and \code{no}.}
  \item{gender}{binary variable indicating gender.}  
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

\source{
Lilien GL, Rangaswamy A (2004). Marketing Engineering: Computer-Assisted Marketing Analysis and Planning, 2nd edition. Victoria, BC: Trafford Publishing. 
}

\details{
The use of a cost matrix is suggested for this data-set.	Classifying a customer that purchased the book as a non-buyer is worse (cost = 5), than it is to classify a costumer that did not purchase the book as a buyer (cost = 1). 
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





