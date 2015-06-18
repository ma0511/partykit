#There are two ways to do this: use a fitter function or a formula interface function.
#With CLM: clm.fit interface works, but there is no predict.clm.fit -> needs to be written first; problem here is that the clm.fit object does not return a model frame with the response or x matrix; need to find a god workaround that does not need to much memory (so appending x toe very fit is not acceptable)
#Second way: use clm and write a wrapper that allows to use it as clm(y,x); also how its done with mlogit 

## wrapper function to specify fitter and return class
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


#' A fitting function for cumulative link logit models. Based on ordinal:clm. 
#'
#' @param formula the model formula
#' @param scale 
#' @param nominal 
#' @param weights case weights
#' @param start  an optional list of starting values of the form c(alpha, beta, zeta) for the thresholds and nominal effects (alpha), regression parameters (beta) and scaleparameters (zeta).
#' @param subset
#' @param na.action
#' @param contrasts
#' @param model
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

## clmfit <- function(y, x, scale=NULL, nominal=NULL, weights, start = NULL, na.action="na.exclude", contrasts=NULL, model=TRUE, link=c("logit","probit","cloglog","loglog","cauchit"), threshold=c("flexible","symmetric","symmetric2","equidistant"), ..., estfun = FALSE, object = FALSE)
##   {
##    #got rid of argument subset, data, scale and nominal 
##   if (!is.null(scale)) stop("Scale regression not implemented yet.")
##   if (!is.null(nominal)) stop("Nominal regression not implemented as cumulative link model. Use mlogittree().")
  
##   if(missing(weights)) weights <- rep(1,nrow(x))
  
##  ## catch control arguments
##   args <- list(...)
##   ctrl <- list()
##   for(n in c("relTol", "maxIter")) {
##     if(n %in% names(args)) {
##       ctrl[[n]] <- args[[n]]
##       args[[n]] <- NULL
##     }
##   }
##   args$control <- do.call("clm.control", ctrl)

##   x <- x[,-which(colnames(x)%in%"(Intercept)")]
  
##   formula <- formula(mf)$formula
##   mf <- data.frame(y,x)
##   colnames(mf)[1] <- as.character(formula)[2]
  
##   ## call clm fitting function
##   args <- c(list(formula=formula, data=mf, weights = weights, model=model, start = start, na.action=na.action, contrasts=contrasts, link=link,threshold=threshold), args)
##   z <- do.call("clm", args)

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


#' A fitting function for cumulative link logit models. Based on ordinal:clm.fit. In loret use it with y=factor(y,ordered=TRUE)
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

clmfit <- function(y, x, s=NULL, n=NULL, start = NULL, weights = rep(1,nrow(x)), offset = rep(0,nrow(x)),S.offset=rep(0,nrow(x)), link=c("logit","probit","cloglog","loglog","cauchit"),threshold=c("flexible","symmetric","symmetric2","equidistant"), ..., estfun = FALSE, object = FALSE)
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
    class(z) <- c("clm.fit")
    z$offset <- offset #made this to be the offset right from the start
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- attr(x, "xlevels")
    z$x <- x #added for predict
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


predict.clm.fit <- function(object,newdata=NULL,type = c("prob", "class"),se.fit=FALSE,interval=FALSE,level=0.95,na.action=na.pass,...)
    {
    #TODO look at predict.clm.fit from predict.clm
    type <- match.arg(type)
    se.fit <- as.logical(se.fit)[1]
    interval <- as.logical(interval)[1]
    stopifnot(length(level) == 1 && is.numeric(level) && level < 1 && level > 0)
    if (type == "class" && (se.fit || interval)) {
        warning("se.fit and interval set to FALSE for type = 'class'")
        se.fit <- interval <- FALSE
    }
    cov <- if (se.fit || interval) 
        object$object$vcov
    else NULL
    has.response <- TRUE
    if (type == "class" && missing(newdata)) 
        newdata <- object$object$x
    has.newdata <- !(missing(newdata) || is.null(newdata))
    if (has.newdata || type == "class") {
        if (has.newdata && sum(unlist(object$aliased)) > 0) 
            warning("predictions from column rank-deficient fit may be misleading")
        newdata <- as.data.frame(newdata)
        resp <- object$object$y 
        if (type == "class") 
            newdata <- newdata[!names(newdata) %in% resp]
        has.response <- resp %in% names(newdata)
        if (!has.response) {
            ylev <- object$y.levels
            nlev <- length(ylev)
            nnd <- nrow(newdata)
            newdata <- cbind(newdata[rep(1:nnd, each = nlev), 
                , drop = FALSE], factor(rep(ylev, nnd), levels = ylev, 
                ordered = TRUE))
            names(newdata)[ncol(newdata)] <- resp
        }
        if (is.null(attr(object$terms, "predvars"))) 
            warning(paste0("terms object does not have a predvars attribute: ", 
                "predictions may be misleading"))
        mf <- model.frame(object$terms, newdata, na.action = na.action, 
            xlev = object$xlevels)
        if (nrow(mf) != nrow(newdata)) 
            stop("length of variable(s) found do not match nrow(newdata)")
        if (!is.null(cl <- attr(object$terms, "dataClasses"))) 
            .checkMFClasses(cl, mf)
        X <- model.matrix(object$terms, mf, contrasts = object$contrasts)
        Xint <- match("(Intercept)", colnames(X), nomatch = 0L)
        n <- nrow(X)
        if (Xint <= 0) 
            X <- cbind(`(Intercept)` = rep(1, n), X)
        if (sum(object$aliased$beta) > 0) 
            X <- X[, !c(FALSE, object$aliased$beta), drop = FALSE]
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(object$terms, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(object$terms, 
                "variables")[[i + 1]], newdata)
        y <- model.response(mf)
        if (any(!levels(y) %in% object$y.levels)) 
            stop(gettextf("response factor '%s' has new levels", 
                response.name(object$terms)))
        if (is.nom <- !is.null(object$nom.terms)) {
            nom.mf <- model.frame(object$nom.terms, newdata, 
                na.action = na.action, xlev = object$nom.xlevels)
            if (nrow(nom.mf) != nrow(newdata)) 
                stop("length of variable(s) found do not match nrow(newdata)")
            if (!is.null(cl <- attr(object$nom.terms, "dataClasses"))) 
                .checkMFClasses(cl, nom.mf)
            NOM <- model.matrix(object$nom.terms, nom.mf, contrasts = object$nom.contrasts)
            NOMint <- match("(Intercept)", colnames(NOM), nomatch = 0L)
            if (NOMint <= 0) 
                NOM <- cbind(`(Intercept)` = rep(1, n), NOM)
            alias <- t(matrix(object$aliased$alpha, nrow = length(object$y.levels) - 
                1))[, 1]
            if (sum(alias) > 0) 
                NOM <- NOM[, !c(FALSE, alias), drop = FALSE]
        }
        if (is.scale <- !is.null(object$S.terms)) {
            S.mf <- model.frame(object$S.terms, newdata, na.action = na.action, 
                xlev = object$S.xlevels)
            if (nrow(S.mf) != nrow(newdata)) 
                stop("length of variable(s) found do not match nrow(newdata)")
            if (!is.null(cl <- attr(object$S.terms, "dataClasses"))) 
                .checkMFClasses(cl, S.mf)
            S <- model.matrix(object$S.terms, S.mf, contrasts = object$S.contrasts)
            Sint <- match("(Intercept)", colnames(S), nomatch = 0L)
            if (Sint <= 0) 
                S <- cbind(`(Intercept)` = rep(1, n), S)
            if (sum(object$aliased$zeta) > 0) 
                S <- S[, !c(FALSE, object$aliased$zeta), drop = FALSE]
            Soff <- rep(0, nrow(S))
            if (!is.null(off.num <- attr(object$S.terms, "offset"))) 
                for (i in off.num) Soff <- Soff + eval(attr(object$S.terms, 
                  "variables")[[i + 1]], newdata)
        }
        tJac <- object$tJac
        dimnames(tJac) <- NULL
        env <- clm.newRho(parent.frame(), y = y, X = X, NOM = if (is.nom) 
            NOM
        else NULL, S = if (is.scale) 
            S
        else NULL, weights = rep(1, n), offset = offset, S.offset = if (is.scale) 
            Soff
        else rep(0, n), tJac = tJac)
        setLinks(env, link = object$link)
    }
    else {
        env <- get_clmRho.clm(object)
    }
    env$par <- as.vector(coef(object))
    env$par <- env$par[!is.na(env$par)]
    pred <- switch(type, prob = prob.predict.clm(env = env, cov = cov, 
        se.fit = se.fit, interval = interval, level = level), 
        class = prob.predict.clm(env = env, cov = cov, se.fit = se.fit, 
            interval = interval, level = level), cum.prob = cum.prob.predict.clm(env = env, 
            cov = cov, se.fit = se.fit, interval = interval, 
            level = level), linear.predictor = lin.pred.predict.clm(env = env, 
            cov = cov, se.fit = se.fit, interval = interval, 
            level = level))
    if (!has.response || type == "class") {
        pred <- lapply(pred, function(x) {
            x <- matrix(unlist(x), ncol = nlev, byrow = TRUE)
            dimnames(x) <- list(1:nrow(x), ylev)
            x
        })
        if (type == "class") 
            pred <- lapply(pred, function(x) {
                factor(max.col(x), levels = seq_along(ylev), 
                  labels = ylev)
            })
    }
    if (missing(newdata) && !is.null(object$na.action)) 
        pred <- lapply(pred, function(x) napredict(object$na.action, 
            x))
    return(pred)
}

predict.clmtree <- function(object, newdata = NULL, type = "node", ...)
{
  ## FIXME: This does not work currently
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
