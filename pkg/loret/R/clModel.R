## simple wrapper function to specify fitter and return class
clmtree <- function(formula, data, subset, na.action, weights, offset, link=c("logit","probit","cloglog","loglog","cauchit"),threshold=c("flexible","symmetric","symmetric2","equidistant"),epsilon = 1e-8, maxit = 25, ...)
{
  if(missing(threshold)) threshold <- "flexible"
  if(missing(link)) link <- "logit"
  if (threshold != "flexible") stop("only flexible thresholds implemented at the moment")
  ## use dots for setting up mob_control
  control <- mob_control(...)

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## extend formula if necessary
  fo <- Formula::Formula(formula)
  if(length(fo)[2L] == 1L) {
    attr(fo, "rhs") <- c(list(1), attr(fo, "rhs"))
    formula[[3L]] <- formula(fo)[[3L]]
  } else {    #not sure this is what we want?
    fo <- NULL
  }

  #Check what to do with link/family as in glmtree
  
  ## call mob
  m <- match.call(expand.dots = FALSE)
  if(!is.null(fo)) m$formula <- formula
  m$fit <- clmfit
  m$control <- control
  m$epsilon <- epsilon
  m$maxit <- maxit
  if("..." %in% names(m)) m[["..."]] <- NULL
  m$link <- link
  m$threshold <- threshold
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  rval$info$link <- link
  rval$info$threshold <- threshold
  class(rval) <- c("clmtree", class(rval))
  return(rval)
}

#' A fitting function for cumulative link logit models. Based on ordinal:clm. In loret use it with y=factor(y,ordered=TRUE)
#'
#' @param y the response variable; an ordered factor.
#' @param x design matrix for regression  parameters
#' @param s design matrix for scale  parameters  
#' @param n design matrix for nominal parameters 
#' @param start  an optional list of starting values of the form c(alpha, beta, zeta) for the thresholds and nominal effects (alpha), regression parameters (beta) and scaleparameters (zeta).
#' @param weights case weights
#' @param offset an optional offset
#' @param S.offset an optional offset for the scale
#' @param link the link function
#' @param threshold the threshold structure, see \link{clm}
#' @param ... Additional control parameters passed to the fitting function, a list of control parameteres or a call to \link{clm.control}  
#' @param estfun estimating functions
#' @param object should the object be saved
#'
#' @return An object of class, a list with components
#' \itemize{
#'         \item{}
#'         \item{}
#'         \item{}
#'         \item{}
#'         \item{}
#'         \item{}
#'         \item{}
#' }
#' 
#' @export
#'
#'       if (missing(X)) 
#'        X <- cbind(`(Intercept)` = rep(1, length(y)))
#'  stopifnot(is.factor(y), is.matrix(X))
#'  if (!missing(S) && !is.null(S)) {
#'        stopifnot(is.matrix(S), length(y) == nrow(S))
#'    }
#'    if (!missing(N) && !is.null(N)) {
#'        stopifnot(is.matrix(N), length(y) == nrow(N))
#'    }

clmfit <- function(y, x, s=NULL, n=NULL, start = NULL, weights = rep(1,nrow(x)), offset = rep(0,nrow(x)),S.offset=rep(0,nrow(x)), link=c("logit","probit","cloglog","loglog","cauchit"),threshold=c("flexible","symmetric","symmetric2","equidistant"), ..., estfun = TRUE, object = TRUE)
  {

  if (!is.null(s)) stop("Scale regression not implemented yet")
  if (!is.null(n)) stop("Nominal regression not implemented as cumulative link model. Use mlogittree().")
  
 ## catch control arguments
  args <- list(...)
  ctrl <- list()
  for(n in c("relTol", "maxIter")) {
    if(n %in% names(args)) {
      ctrl[[n]] <- args[[n]]
      args[[n]] <- NULL
    }
  }
  args$control <- do.call("clm.control", ctrl)
  
  ## call clm fitting function
  args <- c(list(y = y, X = x, weights = weights, offset = offset, S.offset=S.offset,start = start,link=link,threshold=threshold), args)
  z <- do.call("clm.fit", args)

  ## degrees of freedom
  df <- z$edf 

  ## list structure
  rval <- list(
    coefficients = z$coefficients,
    objfun= -z$logLik,  
    estfun = NULL,
    object = NULL
  )

  #add estimating functions (if desired)
  #based on estfun.clm in sandwich
  #only works for flexible threshold and no scale regression
  if(estfun) {
   if(length(link)>1L) link <- link[1]
  # cat("link:", link,"\n")   
   mueta <- make.link(link)$mu.eta
   xmat <- x[, -1L, drop = FALSE]
   n <- nrow(xmat)
   k <- ncol(xmat)
   m <- length(z$alpha)
   y <- as.numeric(y)
   prob <- z$fitted.values
   xb <- if (k >= 1L) 
        as.vector(xmat %*% z$beta)
   else rep(0, n)
   zeta <- z$alpha
   lp <- cbind(0, mueta(matrix(zeta, nrow = n, ncol = m, byrow = TRUE) - 
        xb), 0)
   rrval <- matrix(0, nrow = n, ncol = k + m + 2L)
   if (k >= 1L) 
        rrval[, 1L:k] <- (-xmat * as.vector(lp[cbind(1:n, y + 
            1L)] - lp[cbind(1:n, y)]))
   rrval[cbind(1:n, k + y)] <- -as.vector(lp[cbind(1:n, y)])
   rrval[cbind(1:n, k + y + 1L)] <- as.vector(lp[cbind(1:n, y + 
        1L)])
   rrval <- rrval[, -c(k + 1L, k + m + 2L), drop = FALSE]
   rrval <- weights/prob * rrval
   dimnames(rrval) <- list(rownames(xmat), c(colnames(xmat), names(z$alpha)))
   ix <- if (k >= 1L) 
        c((k + 1L):(k + m), 1L:k)
   else 1L:m
   rval$estfun <- rrval[, ix, drop = FALSE]
}

  ## add model (if desired)
  #TODO: add more info 
  if(object) {
    class(z) <- c("clm.fit", "clm")
    z$offset <- offset #made this to be the offset right from the start
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- attr(x, "xlevels")    

    cl <- as.call(expression(clm))
    cl$formula <- attr(x, "formula")	
    z$call <- cl
    z$terms <- attr(x, "terms")
    #should we collect link information?
    rval$object <- z
  }

  return(rval)
}

## methods
print.clmtree <- function(x,
  title = NULL, objfun = "negative log-likelihood", ...)
{
  if(is.null(title)) title <- paste("Cumulative link model tree for ordinal data (link: ", x$info$link,", threshold: ", x$info$threshold,")",sep="")
  print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.clmtree <- function(object, newdata = NULL, type = "response", ...)
{
  ## FIXME: possible to get default?
  if(is.null(newdata) & !identical(type, "node")) stop("newdata has to be provided")
  predict.modelparty(object, newdata = newdata, type = type, ...)
}



plot.clmtree <- function(x, terminal_panel = node_bivplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  ###FIXME: Should be rewritten; like in psychotree?   
  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(nreg < 1L & missing(terminal_panel)) {
    plot.constparty(as.constparty(x),
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  } else {
    if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
    if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
    plot.modelparty(x, terminal_panel = terminal_panel,
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  }
}
