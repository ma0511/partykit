utils::globalVariables(c(".tree", ".ranef"))

lmertree <- function(formula, data,
  ranefstart = NULL, abstol = 0.001, maxit = 100, 
  joint = TRUE, dfsplit = TRUE, verbose = FALSE, plot = FALSE,
  lmer.control = lmerControl(), ...)
{
  ## remember call
  cl <- match.call()
  
  ## formula processing (full, tree, random)
  ff <- Formula::as.Formula(formula)
  tf <- formula(ff, lhs = 1L, rhs = c(1L, 3L))
  if(length(attr(ff, "rhs")[[2L]]) == 1L) {
    rf <- (. ~ (1 | id))[[3L]]
    rf[[2L]][[3L]] <- attr(ff, "rhs")[[2L]]
    attr(ff, "rhs")[[2L]] <- rf
  }
  if(joint) {
    rf <- formula(ff, lhs = 1L, rhs = 1L)
    rf <- update(rf, . ~ .tree / .)
    rf <- formula(Formula::as.Formula(rf, formula(ff, lhs = 0L, rhs = 2L)),
      lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
  } else {
    rf <- formula(ff, lhs = 1L, rhs = 2L)
  }

  ## initialization
  iteration <- 0L
  data$.ranef <- if (is.null(ranefstart)) {
    rep(0, times = dim(data)[1L])
  } else {
    ranefstart  
  }
  continue <- TRUE
  oldloglik <- -Inf

  ## iterate between lmer and lmtree estimation
  while (continue) {
    iteration <- iteration + 1L

    ## lmtree
    tree <- lmtree(tf, data = data, offset = .ranef, dfsplit = FALSE, ...)
    if(plot) plot(tree)
    data$.tree <- if(joint) {
      factor(predict(tree, newdata = data, type = "node"))
    } else {
      predict(tree, newdata = data, type = "response")
    }

    ## lmer
    if(joint) {
      ## estimate full lmer model but force all coefficients from the
      ## .tree (and the overall intercept) to zero for the prediction
      lme <- lmer(rf, data = data)
      b <- structure(lme@beta, .Names = names(coef(lme)[[1L]]))
      b[substr(names(b), 1L, 5L) %in% c("(Inte", ".tree")] <- 0
      data$.ranef <- suppressWarnings(suppressMessages(predict(lme, newdata = data, newparams = list(beta = b))))
    } else {
      ## estimate only a partial lmer model using the .tree fitted
      ## values as an offset
      lme <- lmer(rf, data = data, offset = .tree)
      data$.ranef <- predict(lme, newdata = data)    
    }

    ## iteration information
    newloglik <- logLik(lme)    
    continue <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if(verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    formula = formula,
    call = cl,
    tree = tree,
    lmer = lme,
    ranef = ranef(lme), 
    varcorr = VarCorr(lme),
    variance = attr(VarCorr(lme),"sc")^2, 
    data = data,
    nobs = nrow(data),
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    dfsplit = dfsplit,
    iterations = iteration, 
    maxit = maxit,
    ranefstart = ranefstart, 
    abstol = abstol,
    mob.control = list(...),
    lmer.control = lmer.control
  )
  class(result) <- "lmertree"
  return(result)
}

glmertree <- function(formula, data, family = "binomial",
  ranefstart = NULL, abstol = 0.001, maxit = 100, 
  joint = TRUE, dfsplit = TRUE, verbose = FALSE, plot = FALSE,
  glmer.control = glmerControl(), ...)
{
  ## remember call
  cl <- match.call()
  
  ## formula processing (full, tree, random)
  ff <- Formula::as.Formula(formula)
  tf <- formula(ff, lhs = 1L, rhs = c(1L, 3L))
  if(length(attr(ff, "rhs")[[2L]]) == 1L) {
    rf <- (. ~ (1 | id))[[3L]]
    rf[[2L]][[3L]] <- attr(ff, "rhs")[[2L]]
    attr(ff, "rhs")[[2L]] <- rf
  }
  if(joint) {
    rf <- formula(ff, lhs = 1L, rhs = 1L)
    rf <- update(rf, . ~ .tree / .)
    rf <- formula(Formula::as.Formula(rf, formula(ff, lhs = 0L, rhs = 2L)),
      lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
  } else {
    rf <- formula(ff, lhs = 1L, rhs = 2L)
  }

  ## initialization
  iteration <- 0L
  data$.ranef <- if (is.null(ranefstart)) {
    rep(0, times = dim(data)[1L])
  } else {
    ranefstart  
  }
  continue <- TRUE
  oldloglik <- -Inf

  ## iterate between glmer and glmtree estimation
  while (continue) {
    iteration <- iteration + 1L

    ## glmtree
    tree <- glmtree(tf, data = data, family = family, offset = .ranef, dfsplit = FALSE, ...)
    if(plot) plot(tree)
    data$.tree <- if(joint) {
      factor(predict(tree, newdata = data, type = "node"))
    } else {
      predict(tree, newdata = data, type = "link")
    }

    ## glmer
    if(joint) {
      ## estimate full glmer model but force all coefficients from the
      ## .tree (and the overall intercept) to zero for the prediction
      glme <- glmer(rf, data = data, family = family)
      b <- structure(glme@beta, .Names = names(coef(glme)[[1L]]))
      b[substr(names(b), 1L, 5L) %in% c("(Inte", ".tree")] <- 0
      data$.ranef <- suppressWarnings(suppressMessages(predict(glme, newdata = data, type = "link", newparams = list(beta = b))))
    } else {
      ## estimate only a partial glmer model using the .tree fitted
      ## values as an offset
      glme <- glmer(rf, data = data, family = family, offset = .tree)
      data$.ranef <- predict(glme, newdata = data, type = "link")
    }

    ## iteration information
    newloglik <- logLik(glme)    
    continue <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if(verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    formula = formula,
    call = cl,
    tree = tree,
    glmer = glme,
    ranef = ranef(glme), 
    varcorr = VarCorr(glme),
    variance = attr(VarCorr(glme),"sc")^2, 
    data = data,
    nobs = nrow(data),
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    dfsplit = dfsplit,
    iterations = iteration, 
    maxit = maxit,
    ranefstart = ranefstart, 
    abstol = abstol,
    mob.control = list(...),
    glmer.control = glmer.control
  )
  class(result) <- "glmertree"
  return(result)
}

coef.lmertree <- coef.glmertree <- function(object, ...) {
  coef(object$tree, ...)
}

plot.lmertree <- plot.glmertree <- function(x, ...) {
  plot(x$tree, ...)
}

ranef.lmertree <- ranef.glmertree <- function(object, ...) {
  object$ranef
}

logLik.lmertree <- logLik.glmertree <- function(object, dfsplit = NULL, ...)
{
  if(is.null(dfsplit)) dfsplit <- object$dfsplit
  dfsplit <- as.integer(dfsplit) * (length(object$tree) - length(nodeids(object$tree, terminal = TRUE)))
  structure(object$loglik, df = object$df + dfsplit, nobs = object$nobs, class = "logLik")
}

print.lmertree <- function(x, title = "Linear mixed model tree", ...) {
  print(x$tree, title = title, ...)
  cat("\nRandom effects:\n")
  print(x$ranef)
  invisible(x)
}

print.glmertree <- function(x, title = "Generalized linear mixed model tree", ...) {
  print(x$tree, title = title, ...)
  cat("\nRandom effects:\n")
  print(x$ranef)
  invisible(x)
}
