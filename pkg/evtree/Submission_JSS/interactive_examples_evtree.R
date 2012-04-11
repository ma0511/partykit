###################################################
### setup
###################################################
library("rpart")
library("evtree")
data("BBBClub", package = "evtree")


###################################################
### 2x2 chessboard example 
###################################################

## generate artificial data
X1 <- rep(seq(0.25, 1.75, 0.5), each = 4)
X2 <- rep(seq(0.25, 1.75, 0.5), 4)
Y <- rep(1, 16)
Y[(X1 < 1 & X2 < 1) | (X1 > 1 & X2 > 1)] <- 2
Y <- factor(Y, labels = c("O", "X"))
chess22 <- data.frame(Y, X1, X2)

## (Figure 1)
plot(X2 ~ X1, data = chess22, xlim = c(0, 2), ylim = c(0, 2), pch = c(1, 4)[Y], col = c("black", "slategray")[Y])

## evtree
set.seed(1090)
print(evtree(Y ~ ., data = chess22, minbucket = 1, minsplit = 2))


###################################################
### BBBClub example
###################################################

## data
data("BBBClub", package = "evtree")

## rpart and ctree
library("rpart")
rp  <- as.party(rpart(choice ~ ., data = BBBClub, minbucket = 10))
rp2 <- as.party(rpart(choice ~ ., data = BBBClub, minbucket = 10, maxdepth = 2))
ct  <- ctree(choice ~ ., data = BBBClub, minbucket = 10, mincrit = 0.99)
ct2 <- ctree(choice ~ ., data = BBBClub, minbucket = 10, mincrit = 0.99, maxdepth = 2)
## Figure 2 (upper and lower panel)
plot(rp)
plot(ct)

## evtree
set.seed(1090)
ev <- evtree(choice ~ ., data = BBBClub, minbucket = 10, maxdepth = 2)
## Figure 3
plot(ev)
ev

## performance comparison
mc <- function(obj) 1 - mean(predict(obj) == BBBClub$choice)
evalfun <- function(obj) 2 * nrow(BBBClub) * mc(obj) +
  width(obj) * log(nrow(BBBClub))
trees <- list("evtree" = ev, "rpart" = rp, "ctree" = ct, "rpart2" = rp2,
  "ctree2" = ct2)
round(sapply(trees, function(obj) c("misclassification" = mc(obj),
  "evaluation function" = evalfun(obj))), digits = 3)
  
## structure comparison
ftable(tab <- table(evtree = predict(ev), rpart  = predict(rp),
  ctree  = predict(ct), observed = BBBClub$choice))
sapply(c("evtree", "rpart", "ctree"), function(nam) {
  mt <- margin.table(tab, c(match(nam, names(dimnames(tab))), 4))
  c(abs = as.vector(rowSums(mt))[2],
    rel = round(100 * prop.table(mt, 1)[2, 2], digits = 3))
})


###################################################
### 4x4 chessboard example 
###################################################

## 4x4 chessboard data generating process
chessboard44 <- function(n = 4000, noisevariables = 6, noise = 0) {
  chess44 <- array(0,c(n,noisevariables+3))
  for(i in 1:(noisevariables+2))
      chess44[,i] <- as.numeric(runif(dim(chess44)[1]))*4

   x <- chess44[,1]
   y <- chess44[,2]
   chess44[, ncol(chess44)] <- 0
   for(k in 1:4)  
      chess44[(x <= k & x > k-1 & y <= k & y > k-1), ncol(chess44)] <- 1
   for(k in 1:2)  
      chess44[(x <= k & x > k-1 & y <= k+2 & y > k+1), ncol(chess44)] <- 1
   for(k in 1:2)  
      chess44[(y <= k & y > k-1 & x <= k+2 & x > k+1), ncol(chess44)] <- 1

   if(noise > 0) {
      flipclasslist <- sample(n, n * (noise / 100), replace = FALSE)

      for(i in 1:length(flipclasslist)){
	  if(chess44[flipclasslist[i], ncol(chess44)] == 1)
	      chess44[flipclasslist[i], ncol(chess44)] = 0
	  else if(chess44[flipclasslist[i], ncol(chess44)] == 0)
	      chess44[flipclasslist[i], ncol(chess44)] = 1
      }
  }

  chess44 <- as.data.frame(chess44)
  chess44[,ncol(chess44)] <- as.factor(chess44[,ncol(chess44)])
  names(chess44) <- c(paste("X", 1:8, sep = ""), "Y")
  chess44
}

## draw data
chess44 <- chessboard44(2000)

## Figure 5
plot(X2 ~ X1, data = chess44, xlim = c(0, 4), ylim = c(0, 4), pch = c(1, 4)[Y], col = c("black", "slategray")[Y])
