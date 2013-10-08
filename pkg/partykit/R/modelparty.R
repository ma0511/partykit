mob <- function(formula, data, subset, na.action, weights, offset,
  fit, control = mob_control(), ...)
{
  ## required packages
  stopifnot(require("Formula"))
  
  ## control parameters (used repeatedly)
  minsplit <- control$minsplit
  if(!is.null(minsplit) && !is.integer(minsplit)) minsplit <- as.integer(minsplit)
  verbose <- control$verbose
  ytype <- control$ytype
  xtype <- control$xtype
  rnam <- c("estfun", "object")
  terminal <- lapply(rnam, function(x) x %in% control$terminal)
  inner    <- lapply(rnam, function(x) x %in% control$inner)
  names(terminal) <- names(inner) <- rnam

  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula FIXME: y ~ . or y ~ x | .
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula(formula(as.Formula(formula(formula), ~ 0), rhs = 2L:1L))
    xreg <- FALSE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1L:2L))
      warning("formula must not have more than two RHS parts")
    }
    xreg <- TRUE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, response, regressor matrix (if any), partitioning variables
  mt <- terms(formula, data = data)
  mtY <- terms(formula, data = data, rhs = if(xreg) 1L else 0L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- switch(ytype,
    "vector" = model.part(formula, mf, lhs = 1L)[[1L]],
    "matrix" = model.matrix(~ 0 + ., model.part(formula, mf, lhs = 1L)),
    "data.frame" = model.part(formula, mf, lhs = 1L)
  )
  X <- if(!xreg) NULL else switch(xtype,
    "matrix" = model.matrix(mtY, mf),
    "data.frame" = model.part(formula, mf, rhs = 1L)
  )
  if(!is.null(X) && ncol(X) < 1L) {
    X <- NULL
    xreg <- FALSE
  }
  if(xreg) {
    attr(X, "formula") <- formula(formula, rhs = 1L)
    attr(X, "terms") <- mtY
  }
  Z <- model.part(formula, mf, rhs = 2L)
  n <- nrow(Z)
  nyx <- length(mf) - length(Z) - as.numeric("(weights)" %in% names(mf)) - as.numeric("(offset)" %in% names(mf))

  ## weights and offset
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1L
  if(length(weights) == 1L) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  offset <- if(xreg) model.offset(mf) else NULL

  ## check fitting function
  fitargs <- names(formals(fit))
  if(!all(c("y", "x", "start", "weights", "offset") %in% fitargs)) {
    stop("no fitting function")
  }

  ## augment fitting function (if necessary)
  if(!all(c("estfun", "object") %in% fitargs)) {
    stopifnot(require("sandwich"))
    afit <- function(y,
      x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
      estfun = FALSE, object = FALSE)
    {
      obj <- fit(y = y, x = x, start = start, weights = weights, offset = offset, ...)
      list(
        coefficients = coef(obj),
        objfun = -as.numeric(logLik(obj)),
	## FIXME: vcov/info/?
        estfun = if(estfun) estfun(obj) else NULL,
        object = if(object) obj else NULL
      )
    }
  } else {
    afit <- fit
  }

  ## convenience functions
  w2n <- function(w) if(control$caseweights) sum(w) else sum(w > 0)
  suby <- function(y, index) {
    if(ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  subx <- if(xreg) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels")   <- attr(x, "xlevels")
      attr(sx, "formula")   <- attr(x, "formula")
      attr(sx, "terms")     <- attr(x, "terms")
      sx
    }
  } else {
    function(x, index) NULL
  }
  subz <- function(z, index) z[index, , drop = FALSE]
  ## from strucchange
  root.matrix <- function(X) {
    if((ncol(X) == 1L)&&(nrow(X) == 1L)) return(sqrt(X)) else {
      X.eigen <- eigen(X, symmetric = TRUE)
      if(any(X.eigen$values < 0)) stop("matrix is not positive semidefinite")
      sqomega <- sqrt(diag(X.eigen$values))
      V <- X.eigen$vectors
      return(V %*% sqomega %*% t(V))
    }
  }

  ## core mob_grow_* functions

  ## variable selection: given model scores, conduct
  ## all M-fluctuation tests for orderins in z
  mob_grow_fluctests <- function(estfun, z, weights)
  {  
    ## set up return values
    m <- NCOL(z)
    pval <- rep.int(0, m)
    stat <- rep.int(0, m)
    ifac <- rep.int(FALSE, m)

    ## estimating functions (dropping zero weight observations)
    process <- as.matrix(estfun)
    ww0 <- (weights > 0)
    process <- process[ww0, , drop = FALSE]
    z <- z[ww0, , drop = FALSE]
    k <- NCOL(process)
    n <- NROW(process)

    ## scale process
    process <- process/sqrt(n)
    J12 <- root.matrix(crossprod(process))
    process <- t(chol2inv(chol(J12)) %*% t(process))  

    ## select parameters to test
    if(!is.null(control$parm)) process <- process[, control$parm, drop = FALSE]
    k <- NCOL(process)

    ## get critical values for supLM statistic
    from <- if(control$trim > 1) control$trim else ceiling(n * control$trim)
    from <- max(from, minsplit)
    to <- n - from
    lambda <- ((n - from) * to)/(from * (n - to))

    beta <- mob_beta_suplm
    logp.supLM <- function(x, k, lambda)
    {
      if(k > 40L) {
        ## use Estrella (2003) asymptotic approximation
        logp_estrella2003 <- function(x, k, lambda)
  	  -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * (1 - k/x) + 2/x))
        ## FIXME: Estrella only works well for large enough x
        ## hence require x > 1.5 * k for Estrella approximation and
        ## use an ad hoc interpolation for larger p-values
        p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
      } else {
        ## use Hansen (1997) approximation
        m <- ncol(beta) - 1L
        if(lambda<1) tau <- lambda
        else tau <- 1/(1 + sqrt(lambda))
        beta <- beta[(((k - 1) * 25 + 1):(k * 25)),]
        dummy <- beta[,(1L:m)]%*%x^(0:(m-1))
        dummy <- dummy*(dummy>0)
        pp <- pchisq(dummy, beta[,(m+1)], lower.tail = FALSE, log.p = TRUE)
        if(tau == 0.5)
          p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
        else if(tau <= 0.01)
          p <- pp[25L]
        else if(tau >= 0.49)
          p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
        else
        {
  	  taua <- (0.51 - tau) * 50
    	  tau1 <- floor(taua)
  	  p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1 + 1L]))
        }
      }
      return(as.vector(p))
    }

    ## compute statistic and p-value for each ordering
    for(i in 1L:m) {
      zi <- z[,i]
      if(is.factor(zi)) {
        proci <- process[order(zi), , drop = FALSE]
        ifac[i] <- TRUE

        # re-apply factor() added to drop unused levels
        zi <- factor(zi[order(zi)])
        # compute segment weights
        segweights <- as.vector(table(zi))/n

        # compute statistic only if at least two levels are left
        if(length(segweights) < 2L) {
          stat[i] <- 0
	  pval[i] <- NA
        } else {      
          stat[i] <- sum(sapply(1L:k, function(j) (tapply(proci[,j], zi, sum)^2)/segweights))
          pval[i] <- pchisq(stat[i], k*(length(levels(zi))-1), log.p = TRUE, lower.tail = FALSE)
        }
      } else {
        oi <- if(control$breakties) {
          mm <- sort(unique(zi))
	  mm <- ifelse(length(mm) > 1L, min(diff(mm))/10, 1)
  	  order(zi + runif(length(zi), min = -mm, max = +mm))
        } else {
          order(zi)
        }    
        proci <- process[oi, , drop = FALSE]
        proci <- apply(proci, 2L, cumsum)
        stat[i] <- if(from < to) {
	  xx <- rowSums(proci^2)
	  xx <- xx[from:to]
	  tt <- (from:to)/n
	  max(xx/(tt * (1 - tt)))	  
	} else {
	  0
	}
        pval[i] <- if(from < to) logp.supLM(stat[i], k, lambda) else NA
      }
    }

    ## select variable with minimal p-value
    best <- which.min(pval)
    if(length(best) < 1L) best <- NA
    rval <- list(pval = exp(pval), stat = stat, best = best)
    names(rval$pval) <- names(z)
    names(rval$stat) <- names(z)
    if(!all(is.na(rval$best)))
      names(rval$best) <- names(z)[rval$best]
    return(rval)
  }

  ### split in variable zselect, either ordered (numeric or ordinal) or nominal
  mob_grow_findsplit <- function(y, x, zselect, weights, offset, ...)
  {
    ## process minsplit (to minimal number of observations)
    if(minsplit > 0.5 & minsplit < 1) minsplit <- 1 - minsplit
    if(minsplit < 0.5) minsplit <- ceiling(w2n(weights) * minsplit)
  
    if(is.numeric(zselect)) {
    ## for numerical variables
      uz <- sort(unique(zselect))
      if (length(uz) == 0L) stop("cannot find admissible split point in z")
      dev <- vector(mode = "numeric", length = length(uz))

      for(i in 1L:length(uz)) {
        zs <- zselect <= uz[i]
        start_left <- NULL
        start_right <- NULL
        if(w2n(weights[zs]) < minsplit || w2n(weights[!zs]) < minsplit) {
          dev[i] <- Inf
        } else {
          fit_left <- afit(y = suby(y, zs), x = subx(x, zs), start = start_left,
	    weights = weights[zs], offset = offset[zs], ...)
          fit_right <- afit(y = suby(y, !zs), x = subx(x, !zs), start = start_right,
	    weights = weights[!zs], offset = offset[!zs], ...)
  	  start_left <- fit_left$coef
	  start_right <- fit_right$coef
	  dev[i] <- fit_left$objfun + fit_right$objfun
        }
      }

      ## maybe none of the possible splits is admissible
      if(all(!is.finite(dev))) {
        split <- list(
          breaks = NULL,
	  index = NULL
        )
      } else {
        split <- list(
          breaks = as.double(uz[which.min(dev)]),
          index = NULL
        )
      }

    } else {

      ## for categorical variables
      al <- mob_grow_getlevels(zselect)
      dev <- apply(al, 1L, function(w) {
        zs <- zselect %in% levels(zselect)[w]
        if(w2n(weights[zs]) < minsplit || w2n(weights[!zs]) < minsplit) {
          return(Inf)
        } else {
	  if(nrow(al) == 1L) 1 else {
            fit_left <- afit(y = suby(y, zs), x = subx(x, zs), start = NULL,
	      weights = weights[zs], offset = offset[zs], ...)
            fit_right <- afit(y = suby(y, !zs), x = subx(x, !zs), start = NULL,
	      weights = weights[!zs], offset = offset[!zs], ...)
    	    fit_left$objfun + fit_right$objfun
	  }
        }
      })

      if(all(!is.finite(dev))) {
        split <- list(
          breaks = NULL,
	  index = NULL
        )
      } else {
        if(is.ordered(zselect)) {
          split <- list(
            breaks = which.min(dev),
	    index = NULL
          )
        } else {
          split <- list(
            breaks = NULL,
	    index = as.integer(!al[which.min(dev),]) + 1L
          )
        }
      }
    }
  
    return(split)
  }

  ## grow tree by combining fluctuation tests for variable selection
  ## and split selection recursively
  mob_grow <- function(id = 1L, y, x, z, weights, offset, ...)
  {
    if(verbose) {
      if(id == 1L) cat("\n")
      cat(sprintf("-- Node %i %s\n", id, paste(rep("-", 32 - floor(log10(id)) + 1L), collapse = "")))
      cat(sprintf("Number of observations: %i\n", w2n(weights)))
    }

    ## fit model
    mod <- afit(y, x, weights = weights, offset = offset, ...,
      estfun = TRUE, object = terminal$object)
    mod$test <- NULL
    mod$nobs <- w2n(weights)
    mod$p.value <- NULL

    ## set default for minsplit if not specified
    if(is.null(minsplit)) minsplit <<- as.integer(ceiling(10L * length(mod$coefficients)/NCOL(y)))

    ## if too few observations: no split = return terminal node
    if(w2n(weights) < 2 * minsplit) {
      if(verbose) {
        cat(sprintf("Too few observations, stop splitting (minsplit = %i)\n\n", minsplit))
      }
      if(!terminal$estfun) mod$estfun <- NULL
      return(partynode(id = id, info = mod))
    }

    ## conduct all parameter instability tests
    test <- try(mob_grow_fluctests(mod$estfun, z, weights))

    if(!inherits(test, "try-error")) {
      if(control$bonferroni) {
        pval1 <- pmin(1, sum(!is.na(test$pval)) * test$pval)
        pval2 <- 1 - (1 - test$pval)^sum(!is.na(test$pval))
        test$pval <- ifelse(!is.na(test$pval) & (test$pval > 0.01), pval2, pval1)
      }

      best <- test$best
      TERMINAL <- is.na(best) || test$pval[best] > control$alpha
      mod$p.value <- as.numeric(test$pval[best])

      if (verbose) {
        cat("\nParameter instability tests:\n")
        print(rbind(statistic = test$stat, p.value = test$pval))
        cat(sprintf("\nBest splitting variable: %s", names(test$stat)[best]))
        cat(sprintf("\nPerform split? %s", ifelse(TERMINAL, "no\n\n", "yes\n")))
      }
    } else {
      if(verbose) cat("Parameter instability tests failed\n\n")
      TERMINAL <- TRUE
      test <- list(stat = NA, pval = NA)
    }
    
    ## update model information
    if(!terminal$estfun) mod$estfun <- NULL
    mod$test <- rbind("statistic" = test$stat, "p.value" = test$pval)

    if(TERMINAL) {
      return(partynode(id = id, info = mod))
    } else {
      zselect <- z[[best]]
      sp <- mob_grow_findsplit(y, x, zselect, weights, offset, ...)
    
      ## split successful?
      if(is.null(sp$breaks) & is.null(sp$index)) {
        if(verbose)
          cat(sprintf("No admissable split found in %s\n\n", sQuote(names(test$stat)[best])))
        return(partynode(id = id, info = mod))
      } else {
        sp <- partysplit(as.integer(best), breaks = sp$breaks, index = sp$index)
        if(verbose) cat(sprintf("Selected split: %s\n\n",
	  paste(character_split(sp, data = z)$levels, collapse = " | ")))
      }
    }
  
    ## actually split the data
    kidids <- kidids_split(sp, data = z)
    
    ## set-up all daugther nodes
    kids <- vector(mode = "list", length = max(kidids))
    for(kidid in 1L:max(kidids)) {
      ## select obs for current next node
      nxt <- kidids == kidid

      ## get next node id
      if(kidid > 1L) {
        myid <- max(nodeids(kids[[kidid - 1L]]))
      } else {
        myid <- id
      }

      ## start recursion on this daugther node
      kids[[kidid]] <- mob_grow(id = myid + 1L,
        suby(y, nxt), subx(x, nxt), subz(z, nxt), weights[nxt], offset[nxt], ...)
    }

    ## shift split varid from z to mf
    sp$varid <- as.integer(sp$varid + nyx)
    
    ## return nodes
    return(partynode(id = id, split = sp, kids = kids, info = mod))
  }

  ## grow tree
  nodes <- mob_grow(id = 1L, Y, X, Z, weights, offset, ...)

  ## compute terminal node number for each observation
  fitted <- fitted_node(nodes, data = mf)
  fitted <- data.frame(
      "(fitted)" = fitted,
      ## "(response)" = Y, ## probably not really needed
      check.names = FALSE,
      row.names = rownames(mf))
  if(!identical(weights, rep.int(1L, n))) fitted[["(weights)"]] <- weights
  if(!is.null(offset)) fitted[["(offset)"]] <- offset

  ## update minsplit (in case it was set internally)
  control$minsplit <- minsplit

  ## return party object
  rval <- party(nodes, 
    data = if(control$model) mf else mf[0,],
    fitted = fitted,
    terms = mt,
    info = list(
      call = cl,
      formula = oformula,
      Formula = formula,
      terms = list(response = mtY, partitioning = mtZ),
      fit = afit,
      control = control,
      dots = list(...),
      nreg = as.integer(xreg) * (nyx - NCOL(Y)))
  )
  class(rval) <- c("modelparty", class(rval))
  return(rval)
}

## determine all possible splits for a factor, both nominal and ordinal
mob_grow_getlevels <- function(z) {
  nl <- nlevels(z)
  if(inherits(z, "ordered")) {
    indx <- diag(nl)
    indx[lower.tri(indx)] <- 1
    indx <- indx[-nl,]
    rownames(indx) <- levels(z)[-nl]
  } else {
    mi <- 2^(nl - 1L) - 1L
    indx <- matrix(0, nrow = mi, ncol = nl)
    for (i in 1:mi) {
      ii <- i
      for (l in 1:nl) {
        indx[i, l] <- ii %% 2L
	ii <- ii %/% 2L   
      }
    }
    rownames(indx) <- apply(indx, 1L, function(x) paste(levels(z)[x > 0], collapse = "+"))
  }
  colnames(indx) <- as.character(levels(z))
  storage.mode(indx) <- "logical"
  indx
}

## control splitting parameters
mob_control <- function(alpha = 0.05, bonferroni = TRUE, minsplit = NULL, trim = 0.1,
  breakties = FALSE, parm = NULL, verbose = FALSE, caseweights = TRUE,
  ytype = "vector", xtype = "matrix", terminal = "object", inner = terminal,
  model = TRUE)
{
  ytype <- match.arg(ytype, c("vector", "data.frame", "matrix"))
  xtype <- match.arg(xtype, c("data.frame", "matrix"))

  if(!is.null(terminal)) terminal <- as.vector(sapply(terminal, match.arg, c("estfun", "object")))
  if(!is.null(inner))    inner    <- as.vector(sapply(inner,    match.arg, c("estfun", "object")))

  rval <- list(alpha = alpha, bonferroni = bonferroni, minsplit = minsplit,
    trim = ifelse(is.null(trim), minsplit, trim),
    breakties = breakties, parm = parm, verbose = verbose,
    caseweights = caseweights, ytype = ytype, xtype = xtype,
    terminal = terminal, inner = inner, model = model)
  return(rval)
}

## methods concerning call/formula/terms/etc.
## (default methods work for terms and update)

formula.modelparty <- function(x, extended = FALSE, ...)
  if(extended) x$info$Formula else x$info$formula

getCall.modelparty <- function(x, ...) x$info$call

model.frame.modelparty <- function(formula, ...)
{
  mf <- formula$data
  if(nrow(mf) > 0L) return(mf)

  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- formula$info$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf[names(nargs)] <- nargs
  if(is.null(env <- environment(formula$info$terms))) env <- parent.frame()
  mf$formula <- Formula(as.formula(mf$formula))
  eval(mf, env)
}

weights.modelparty <- function(object, ...) {
  fit <- object$fitted
  ww <- if(!is.null(w <- fit[["(weights)"]])) w else rep.int(1L, NROW(fit))
  structure(ww, .Names = rownames(fit))
}

## methods concerning model/parameters/loglik/etc.

coef.modelparty <- function(object, node = NULL, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = TRUE)
  drop(do.call("rbind",
    nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)))
}

refit.modelparty <- function(object, node = NULL, drop = TRUE, ...)
{
  ## by default use all ids
  if(is.null(node)) node <- nodeids(object)
  
  ## estimated coefficients
  cf <- nodeapply(object, ids = node, FUN = function(n) info_node(n)$coefficients)
  
  ## model.frame
  mf <- model.frame(object)
  weights <- weights(object)
  offset <- model.offset(mf)
  
  ## fitted ids
  fitted <- object$fitted[["(fitted)"]]

  ## response variables
  Y <- switch(object$info$control$ytype,
    "vector" = model.part(object$info$Formula, mf, lhs = 1L)[[1L]],
    "matrix" = model.matrix(~ 0 + ., model.part(object$info$Formula, mf, lhs = 1L)),
    "data.frame" = model.part(object$info$Formula, mf, lhs = 1L)
  )
  X <- if(object$info$nreg < 1L) NULL else switch(object$info$control$xtype,
    "matrix" = model.matrix(object$info$terms$response, mf),
    "data.frame" = model.part(object$info$Formula, mf, rhs = 1L)
  )
  if(!is.null(X)) {
    attr(X, "formula") <- formula(object$info$Formula, rhs = 1L)
    attr(X, "terms") <- object$info$terms$response
  }

  suby <- function(y, index) {
    if(object$info$control$ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  subx <- if(object$info$nreg > 0L) {
    function(x, index) {
      sx <- x[index, , drop = FALSE]
      attr(sx, "contrasts") <- attr(x, "contrasts")
      attr(sx, "xlevels")   <- attr(x, "xlevels")
      attr(sx, "formula")   <- attr(x, "formula")
      attr(sx, "terms")     <- attr(x, "terms")
      sx
    }
  } else {
    function(x, index) NULL
  }
  
  ## fitting function
  afit <- object$info$fit
  
  ## refit
  rval <- lapply(seq_along(node), function(i) {
    ix <- fitted %in% nodeids(object, from = node[i], terminal = TRUE)
    args <- list(y = suby(Y, ix), x = subx(X, ix), start = cf[[i]],
      weights = weights[ix], offset = offset[ix], object = TRUE)
    args <- c(args, object$info$dots)
    do.call("afit", args)$object
  })
  names(rval) <- node

  ## drop?
  if(drop & length(rval) == 1L) rval <- rval[[1L]]
  
  ## return
  return(rval)
}

apply_to_models <- function(object, node = NULL, FUN = NULL, drop = FALSE, ...) {
  if(is.null(node)) node <- nodeids(object, terminal = FALSE)
  if(is.null(FUN)) FUN <- function(object, ...) object  
  rval <- if("object" %in% object$info$control$terminal) {
    nodeapply(object, node, function(n) FUN(info_node(n)$object))
  } else {
    lapply(refit.modelparty(object, node, drop = FALSE), FUN)
  }
  names(rval) <- node
  if(drop & length(node) == 1L) rval <- rval[[1L]]
  return(rval)
}

logLik.modelparty <- function(object, splits = TRUE, ...)
{
  ids <- nodeids(object, terminal = TRUE)
  ll <- apply_to_models(object, node = ids, FUN = logLik)
  df <- if(splits) length(ll) - 1L else 0L  
  structure(
    sum(as.numeric(ll)),
    df = sum(sapply(ll, function(x) attr(x, "df"))) + df,
    nobs = nobs(object),
    class = "logLik"
  )
}

nobs.modelparty <- function(object, ...) {
  sum(unlist(nodeapply(object,
    nodeids(object, terminal = TRUE),
    function(n) info_node(n)$nobs
  )))
}

deviance.modelparty <- function(object, ...)
{
  ids <- nodeids(object, terminal = TRUE)
  dev <- apply_to_models(object, node = ids, FUN = deviance)
  sum(unlist(dev))
}

summary.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = TRUE) else node
  rval <- apply_to_models(object, node = ids, FUN = summary)
  if(length(ids) == 1L) rval[[1L]] else rval
}

sctest.modelparty <- function(object, node = NULL, ...)
{
  ids <- if(is.null(node)) nodeids(object, terminal = FALSE) else node
  rval <- nodeapply(object, ids, function(n) info_node(n)$test)
  names(rval) <- ids
  if(length(ids) == 1L) rval[[1L]] else rval  
}

print.modelparty <- function(x, node = NULL,
  FUN = NULL, digits = getOption("digits") - 4L,
  header = TRUE, footer = TRUE,
  title = "Model-based recursive partitioning", objfun = "", ...)
{
  digits <- max(c(0, digits))
  if(objfun != "") objfun <- paste(" (", objfun, ")", sep = "")

  if(is.null(node)) {
    header_panel <- if(header) function(party) {      
      c(title, "", "Model formula:", deparse(party$info$formula), "", "Fitted party:", "")
    } else function(party) ""
  
    footer_panel <- if(footer) function(party) {
      n <- width(party)
      n <- format(c(length(party) - n, n))
      info <- nodeapply(x, ids = nodeids(x, terminal = TRUE),
        FUN = function(n) c(length(info_node(n)$coefficients), info_node(n)$objfun))
      k <- mean(sapply(info, "[", 1L))
      of <- format(sum(sapply(info, "[", 2L)), digits = getOption("digits"))

      c("", paste("Number of inner nodes:   ", n[1L]),
        paste("Number of terminal nodes:", n[2L]),
        paste("Number of parameters per node:", format(k, digits = getOption("digits"))),
        paste("Objective function", objfun, ": ", of, sep = ""), "")
    } else function (party) ""

    if(is.null(FUN)) {
      FUN <- function(x) c(sprintf(": n = %i", x$nobs), capture.output(print(x$coefficients)))
    }
    terminal_panel <- function(node) formatinfo_node(node,
      default = "*", prefix = NULL, FUN = FUN)

    print.party(x, terminal_panel = terminal_panel,
      header_panel = header_panel, footer_panel = footer_panel, ...)
  } else {
    node <- as.integer(node)
    info <- nodeapply(x, ids = node,
      FUN = function(n) info_node(n)[c("coefficients", "objfun", "test")])    
    for(i in seq_along(node)) {
      if(i == 1L) {
        cat(paste(title, "\n", collapse = ""))
      } else {
        cat("\n")
      }
      cat(sprintf("-- Node %i --\n", node[i]))
      cat("\nEstimated parameters:\n")
      print(info[[i]]$coefficients)
      cat(sprintf("\nObjective function:\n%s\n", format(info[[i]]$objfun)))
      cat("\nParameter instability tests:\n")
      print(info[[i]]$test)
    }
  }
  invisible(x)
}

predict.modelparty <- function(object, newdata = NULL, type = "node", ...)
{
  ## predicted node ids
  node <- predict.party(object, newdata = newdata)
  if(identical(type, "node")) return(node)

  ## obtain fitted model objects
  ids <- sort(unique(as.integer(node)))
  mod <- apply_to_models(object, node = ids)

  ## obtain predictions
  pred <- if(is.character(type)) {
    function(object, newdata = NULL, ...) predict(object, newdata = newdata, type = type, ...)
  } else {
    type
  }
  if("newdata" %in% names(formals(pred))) {    
    ix <- lapply(seq_along(ids), function(i) which(node == ids[i]))
    preds <- lapply(seq_along(ids), function(i)
      pred(mod[[i]], newdata = newdata[ix[[i]], , drop = FALSE], ...))
    nc <- NCOL(preds[[1L]])
    rval <- if(nc > 1L) {
      matrix(0, nrow = length(node), ncol = nc, dimnames = list(names(node), colnames(preds[[1L]])))
    } else {
      rep(preds[[1L]], length.out = length(node))
    }      
    for(i in seq_along(ids)) {
      if(nc > 1L) {
        rval[ix[[i]], ] <- preds[[i]]
	rownames(rval) <- names(node)
      } else {
        rval[ix[[i]]] <- preds[[i]]
	names(rval) <- names(node)
      }
    }
  } else {
    rval <- lapply(mod, pred, ...)
    if(NCOL(rval[[1L]]) > 1L) {
      rval <- do.call("rbind", rval)
      rownames(rval) <- ids
      rval <- rval[as.character(node), , drop = FALSE]
      rownames(rval) <- names(node)
      rval <- drop(rval)
    } else {
      ## provide a c() method for factors locally
      c.factor <- function(...) {
        args <- list(...)
	lev <- levels(args[[1L]])
	args[[1L]] <- unclass(args[[1L]])
	rval <- do.call("c", args)
	factor(rval, levels = 1L:length(lev), labels = lev)
      }
      rval <- do.call("c", rval)
      names(rval) <- ids
      rval <- rval[as.character(node)]
      names(rval) <- names(node)
    }
  }
  return(rval)
}

fitted.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  fit <- apply_to_models(object, node = ids, FUN = fitted)

  nc <- NCOL(fit[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(fit[[1L]])))
  } else {
    rep(fit[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- fit[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- fit[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

residuals.modelparty <- function(object, ...)
{
  ## fitted nodes
  node <- predict.party(object, type = "node")

  ## obtain fitted model objects
  ids <- nodeids(object, terminal = TRUE)
  res <- apply_to_models(object, node = ids, FUN = residuals)

  nc <- NCOL(res[[1L]])
  rval <- if(nc > 1L) {
    matrix(0, nrow = length(ids), ncol = nc, dimnames = list(names(node), colnames(res[[1L]])))
  } else {
    rep(res[[1L]], length.out = length(node))
  }	 
  for(i in seq_along(ids)) {
    if(nc > 1L) {
      rval[node == ids[i], ] <- res[[i]]
      rownames(rval) <- names(node)
    } else {
      rval[node == ids[i]] <- res[[i]]
      names(rval) <- names(node)
    }
  }

  return(rval)
}

plot.modelparty <- function(x, terminal_panel = NULL, FUN = NULL, ...) {
  if(is.null(terminal_panel)) {
    if(is.null(FUN)) {
      FUN <- function(x) {
        cf <- x$coefficients
	cf <- matrix(cf, ncol = 1, dimnames = list(names(cf), ""))
        c(sprintf("n = %i", x$nobs), "Estimated parameters:",
          strwrap(capture.output(print(cf, digits = 4L))[-1L]))
      }
    }
    terminal_panel <- node_terminal(x, FUN = FUN)
  }
  plot.party(x, terminal_panel = terminal_panel, ...)
}
