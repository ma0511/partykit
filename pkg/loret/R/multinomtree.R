#This here is done with the formula interface
## wrapper function to specify fitter and return class
multinomtree <- function(formula, data, subset, na.action, weights, offset, epsilon = 1e-8, maxit = 25, ...)
{
  ## use dots for setting up mob_control
  control <- mob_control(...)

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## extend formula if necessary
  fo <- Formula::Formula(formula)
  if(length(fo)[2L] == 1L) {
    attr(fo, "rhs") <- c(list(1), attr(fo, "rhs"))
    formula[[3L]] <- formula(fo)[[3L]]
  } else { 
    fo <- NULL
  }
  
  ## call mob
  m <- match.call(expand.dots = FALSE)
  if(!is.null(fo)) m$formula <- formula
  m$fit <- multinomfit
  m$control <- control
  m$epsilon <- epsilon
  m$maxit <- maxit
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("multinomtree", class(rval))
  return(rval)
}


#' A fitting function for multinomial logit models. Based on nnet:multinom
#'
#' @param y 
#' @param x 
#' @param weights case weights
#' @param start  an optional list of starting values of the form c(alpha, beta, zeta) for the thresholds and nominal effects (alpha), regression parameters (beta) and scaleparameters (zeta).
#' @param subset
#' @param na.action
#' @param contrasts
#' @param Hess
#' @param summ
#' @param censored
#' @param model 
#' @param ... Additional control parameters passed to the fitting function, a list of control parameteres or a call to \link{clm.control}  
##' @param estfun estimating functions
##' @param object should the object be saved
## '
## ' @return An object of class, a list with components
## ' \itemize{
## '         \item{}
## '         \item{}
## '         \item{}
## '         \item{}
## '         \item{}
## '         \item{}
## '         \item{}
## ' }
## ' 
## ' @export

## library(MASS)
## bwt <- with(birthwt, {
## race <- factor(race, labels = c("white", "black", "other"))
## ptd <- factor(ptl > 0)
## ftv <- factor(ftv)
## levels(ftv)[-(1:2)] <- "2+"

## data.frame(low = factor(low), age, lwt, race, smoke = (smoke > 0), ptd, ht = (ht > 0), ui = (ui > 0), ftv)
## })

## f1 <- multinom(ftv~age+lwt+race+smoke+ptd+ht+ui+low,data=bwt)

## f2 <- multinom(ftv~age+lwt+race+smoke+ptd+ht+ui+low,data=bwt,softmax=TRUE)

## x <- model.matrix(f)
## y <- model.response(f)

## xbu <- x

## df <- cbind(y,x)

## colnames(x) <- sub(":","",colnames(x))
## colnames(x) <- sub("(","",colnames(x))
## colnames(x) <- sub(")","",colnames(x))


## colnames(x) <- c("boatintercept" ,   "charterintercept", "pierintercept"   , "price"          ,    "catch")

## class(x) <- "mlogit.data"

## forms <- formula(paste(colnames(y),"~0+",paste(colnames(x),collapse="+")))

## test1 <- mlogit(mode ~ price + catch,data=df)

## test1 <- mlogit(mode ~ price + catch,data=Fish)

## mlogit(mode ~ 0 + I(boat(intercept)) + I(charter(intercept)) + I(pier(intercept)) + 
##     price + catch,data=df)


## multinom(mode ~ boatintercept + charterintercept + pierintercept + price + catch,data=df)

## mnlogit(mode ~ boatintercept + charterintercept + pierintercept + price + catch,data=df,choiceVar="mode")

## df$mode <- df$mode*1


#multinom.fits <- function(y, x, weights, start = NULL, offset=NULL,...)
#    {
#       multinom(y ~ 0 + x, start = start, Hess=TRUE, ...)
#    }

## formos <- low~age+lwt+race+smoke+ptd+ht+ui+ftv

## formos <- ftv ~ low | age+lwt+race+smoke+ptd+ht+ui

## test_tree <- partykit::mob(formos, data = bwt, fit = multinom.fits)

## test_tree


## test_tree <- partykit::mob(formos, data = bwt, fit = multinomfit)

## y <- bwt$ftv
## x <- model.matrix(f1)

## teso <- multinomfit(y,x)


 #to get coefficients in a vector other than the multinom standard 
coefficients.multinom <- function(object,...)
    {                                 
    return(coef(object))
    }

#coefficients(f1)
    
 #to get coefficients in a vector other than the multinom standard 
coef.multinom <- function(object,...)
    {
    origs <- nnet:::coef.multinom(object)
    coefs <- as.numeric(origs)
    names(coefs) <- Reduce(function(x,y) paste(x,y,sep=":"),expand.grid(attr(origs,"dimnames")[[1]],attr(origs,"dimnames")[[2]]))
    if(length(object$lev)==2L) names(coefs) <- paste(rep(object$lev[2],length(coefs)),names(origs),sep=":")
    return(coefs)
   }


#deleted many arguments - check if we need them
multinomfit <- function(y, x, weights, start = NULL, offset=NULL, Hess=TRUE, summ=0, model=FALSE, trace=FALSE,..., estfun = FALSE, object = FALSE)
  {  

 if(missing(weights)) weights <- rep(1,nrow(x))
  
 ## catch control arguments
 addargs <- list(...)

 #set up design matrix 
  xo <- x

  x <- x[,-which(colnames(x)%in%"(Intercept)"),drop=FALSE]
  
  y <- as.data.frame(y,ncol=1)
 
  mf <- data.frame(y,x)
   
  forms <- formula(paste(colnames(y),"~",paste(colnames(x),collapse="+")))
  
  ## call multinom fitting function
  argos <- c(list(formula=forms, data=mf, weights = weights,  start = start, Hess=Hess, summ=summ, model=model,trace=trace),addargs) #addargs has issues here
  z <- do.call("multinom", argos)

  coefs <- coefficients.multinom(z)

  ## list structure
  rval <- list(
    coefficients = coefs,
    objfun=  z$deviance/2, #(z$AIC/2-length(coefficients(z)))
    estfun = NULL,
    object = NULL
  )

  #estimating functions
  if(estfun) {
    #browser()  
    xmat <- model.matrix(z)
    m <- length(z$lev) #hier ist das Problem - der droppt die Kategorien dann stimmt das alles nimma; er kann dann nicht mt nur 2 Katgeorien unten weitermachen bei 
    n <- nrow(xmat)
    k <- ncol(xmat)
    w <- model.weights(model.frame(z))
    if(is.null(w)) w <- rep(1, n)
    res <- residuals(z) #get resid
    if(m==2L) res <- cbind(-res,res)  #if there are only 2 categories, residuals() returns a nx1 matrix instead a nx2 matrix; in case of more than 2 categories it is nxk  so we augment it here accordingly that we can do it the same way later on - funtion returns the nx k-1 estimating functions though, the baselines are kicked out 
    wres <- res*w #weighted residuals
   #calculate the estimating functions as wres_k*xmat for each k    
    wresxmat <- lapply(1:ncol(wres), function(x) wres[,x] * xmat) #wres*xmat for each category individually
    wresxmat <- matrix(unlist(wresxmat), ncol = m*ncol(xmat)) #turn it into an nxk matrix
    colnam <- paste(rep(z$lev,each=ncol(xmat)),rep(colnames(xmat),m),sep=":") #naming
   # browser()
    #problem!
    attr(wresxmat,"dimnames") <- list(rownames(xmat), colnam) #naming
#    coefs <- coefficients.multinom(z) #take the overloaded coefficients function here
    rval$estfun <- wresxmat[,-grep(paste(z$lev[1],":",sep=""),dimnames(wresxmat)[[2]])] #kick out those columns that are baseline category
   # rval$estfun <- estfun.multinom(z)
}

  ## add model (if desired)
  if(object) {
    class(z) <- c("multinom","nnet")
    #z$offset <- offset #made this to be the offset right from the start
    z$contrasts <- attr(xo, "contrasts")
    z$xlevels <- attr(xo, "xlevels")
    z$x <- x
    cl <- as.call(expression(multinom))
    cl$formula <- attr(xo, "formula")	
    z$call <- cl
    z$terms <- attr(xo, "terms")
    rval$object <- z
  }

  return(rval)
}


## multinom.fit <- function(y, x, data, weights, start = NULL, reflevel=NULL, nests=NULL, un.nest.el=FALSE, unscaled =FALSE, hetersc=FALSE, rpar=NULL, probit=FALSE, R=40, correlation=FALSE, halton=NULL, random.nb=NULL, panel=FALSE, estiamte=TRUE,seed=10, ..., estfun = FALSE, object = FALSE)
##   {  

##  if(missing(weights)) weights <- rep(1,nrow(x))
  
##  ## catch control arguments
##   args <- list(...)
##   ctrl <- list()
##   for(n in c("fTol", "iterlim")) {
##     if(n %in% names(args)) {
##       ctrl[[n]] <- args[[n]]
##       args[[n]] <- NULL
##     }
##   }
##   args$control <- do.call("mlogit.optim", ctrl)



##  x <- x[,-which(colnames(x)%in%"(Intercept)")]


##  ##   x <- x[,-which(colnames(x)%in%"(Intercept)")]
  
## ##   formula <- formula(mf)$formula
## ##   mf <- data.frame(y,x)
## ##   colnames(mf)[1] <- as.character(formula)[2]
  
##   mf <- cbind(y,x)
   
##   forms <- formula(paste(colnames(y),"~",paste(paste("'",colnames(x),"'",sep=""),collapse="+")))
  
##   ## call clm fitting function
##   args <- c(list(formula=formula, data=mf, weights = weights, model=model, start = start, na.action=na.action, contrasts=contrasts, link=link,threshold=threshold), args)
##   z <- do.call("mlogit", args)

##   ## z2 <- rapply(z, function(x) {
##   ##                   if(is.character(x)) x <- sub("^x","", x) else x
##   ##               }
##   ##             , how="replace")

##   ## z2 <- rapply(z2, function(x) {
##   ##                        if(!is.null(attr(x,"names"))) attr(x,"names") <- sub("^x","", attr(x,"names")) else x
##   ##               }
##   ##             , how="replace")

##   ## z2 <- rapply(z2, function(x) {
##   ##                  if(!is.null(attr(x,"dimnames"))) lapply(attr(x,"dimnames"), function(l) l <- sub("^x","", l)) else x }
##   ##             , how="replace")

##   ## z2
##   ## attr(z$coefficients,"names") <- sub("^x", "", attr(z$coefficients,"names"))

      
##   ## degrees of freedom
##   df <- z$edf 

##   ## list structure
##   rval <- list(
##     coefficients = z$coefficients,
##     objfun= -z$logLik,  
##     estfun = NULL,
##     object = NULL
##   )

##   #add estimating functions (if desired)
##   #based on estfun.clm in sandwich
##   #only works for flexible threshold and no scale regression
##   if(estfun) {
##    rval$estfun <- sandwich::estfun(z)
## }

##   ## add model (if desired)
##   #TODO: add more info 
##   if(object) {
##    # class(z) <- c("clm")
##    # z$offset <- offset #made this to be the offset right from the start
##   #  z$contrasts <- attr(x, "contrasts")
##   #  z$xlevels <- attr(x, "xlevels")
##     #z$x <- x
##   #  cl <- as.call(expression(clm))
##    # cl$formula <- attr(x, "formula")	
##    # z$call <- cl
##    # z$terms <- attr(x, "terms")
##     #should we collect link information?
##     rval$object <- z
##   }

##   return(rval)
## }
