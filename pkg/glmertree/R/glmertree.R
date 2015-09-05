lmertree <- function(formula, randomformula, data,
  ranefstart = NULL, abstol = 0.001, maxit = 1000, 
  verbose = FALSE, plot = FALSE, lmer.control = lmerControl(), ...)
{
  ## initialization
  iteration <- 0L
  data$.ranef <- .ranefstart <- if (is.null(ranefstart)) {
    rep(0, times = dim(data)[1L])
  } else {
    ranefstart  
  }
  continuecondition <- TRUE
  oldloglik <- -Inf

  ## iterate between lmer and lmtree estimation  
  while (continuecondition) {
    iteration <- iteration + 1L

    ## lmtree
    tree <- lmtree(formula, data, offset = .ranef, ...)
    if(plot) plot(tree)
    data$.treeresponse <- .treeresponse <- predict(tree, newdata = data)

    ## lmer
    lme <- lmer(randomformula, data = data, offset = .treeresponse)
    data$.ranef <- .ranef <- predict(lme, newdata = data)

    ## iteration information
    newloglik <- logLik(lme)    
    continuecondition <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if(verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    tree = tree,
    lmer = lme,
    ranef = ranef(lme), 
    varcorr = VarCorr(lme),
    variance = attr(VarCorr(lme),"sc")^2, 
    data = data,
    logLik = newloglik,
    iterations = iteration, 
    maxit = maxit,
    ranefstart = ranefstart, 
    formula = formula,
    randomformula = randomformula,
    abstol = abstol,
    mob.control = list(...),
    lmer.control = lmer.control
  )
  class(result) <- "lmertree"
  return(result)
}

glmertree <- function(formula, randomformula, data,
  family = "binomial", ranefstart = NULL, abstol = 0.001, maxit = 1000, 
  verbose = FALSE, plot = FALSE, glmer.control = glmerControl(), ...)
{
  ## initialization
  iteration <- 0L
  data$.ranef <- .ranef <- if (is.null(ranefstart)) {
    rep(0, times = dim(data)[1L])
  } else {
    ranefstart
  }
  continuecondition <- TRUE
  oldloglik <- -Inf
  
  ## iterate between glmer and glmtree estimation
  while (continuecondition) {
    iteration <- iteration + 1L

    ## glmtree
    tree <- glmtree(formula, data, family = family, offset = .ranef, ...)
    if(plot) plot(tree)
    data$.treeresponse <- .treeresponse <- predict(tree, newdata = data, type = "link")

    ## glmer
    glme <- glmer(randomformula, family = family, data = data, offset = .treeresponse)
    data$.ranef <- .ranef <- predict(glme, newdata = data, type = "link")

    ## iteration information
    newloglik <- logLik(glme)    
    continuecondition <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if(verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    tree = tree,
    glmer = glme,
    ranef = ranef(glme), 
    varcorr = VarCorr(glme),
    variance = attr(VarCorr(glme), "sc")^2, 
    data = data,
    logLik = newloglik,
    iterations = iteration, 
    maxit = maxit,
    ranefstart = ranefstart, 
    formula = formula,
    randomformula = randomformula,
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

print.lmertree <- function(x, title = "Linear mixed model tree", ...) {
  print(x$tree, title = title, ...)
  cat("\nRandom effects:\n")
  print(x$ranef)
  invisible(x)
}

print.glmertree <- function(x, title = "General linear mixed model tree", ...) {
  print(x$tree, title = title, ...)
  cat("\nRandom effects:\n")
  print(x$ranef)
  invisible(x)
}
