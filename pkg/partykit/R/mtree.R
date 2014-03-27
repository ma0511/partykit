## auxiliary function for MLR which checks for separation
mconstfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
  ..., estfun = FALSE, object = FALSE)
{
  if(is.factor(y)) {

  ## tables and probabilities
  tab <- tapply(weights, y, sum)
  tab[is.na(tab)] <- 0L
  pr <- tab/sum(tab)
  alias <- tab == 0L
  ix1 <- which(!alias)[1L]
  if(estfun) ef <- matrix(0, nrow = length(y), ncol = length(tab),
    dimnames = list(names(y), names(tab)))
  
  if(sum(!alias) < 2L) {
    return(list(
      coefficients = log(pr[-ix1]) - log(pr[ix1]),
      objfun = 0,
      estfun = NULL,
      object = NULL
    ))
  }
  
  ## information required for mob()
  rval <- list(
    coefficients = log(pr[-ix1]) - log(pr[ix1]),
    objfun = -sum(tab[tab > 0L] * log(pr[tab > 0L])),
    estfun = NULL,
    object = NULL
  )
  if(estfun) {
    cf <- log(pr) - log(pr[ix1])
    ef[] <- rep(-pr, each = nrow(ef))
    ef[cbind(1:nrow(ef), as.numeric(y))] <- (1 - pr[y])
    ef <- ef[, !alias, drop = FALSE]
    ef <- ef[, -1L, drop = FALSE]
    rval$estfun <- ef * weights
  }  
  
  } else {
  
  if(!is.matrix(y)) y <- matrix(y, ncol = 1L)
  cf <- colMeans(y * weights)/mean(weights)
  res <- y - rep(cf, each = NROW(y))
  rval <- list(
    coefficients = cf,
    objfun = sum(res^2 * weights),
    estfun = if(estfun) res * weights else NULL,
    object = NULL
  )
  
  }

  return(rval)
}

mtree <- function(formula, data, na.action, weights, subset, minsplit = 7L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(minsplit = minsplit, ...)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  stopifnot(require("Formula"))
  m$formula <- formula(as.Formula(formula, ~ 1), rhs = 2:1)
  m$fit <- mconstfit
  m$control <- control
  m$minsplit <- NULL
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  rval <- as.constparty(rval)
  class(rval) <- c("constparty", "party", "modelparty")
  return(rval)
}


if(FALSE) {

## multinomial logistic regression MOB
mt <- mtree(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = iris)
plot(mt)

## least squares MOB
airq <- subset(airquality, !is.na(Ozone) & !is.na(Solar.R))
mt2 <- mtree(Ozone ~ Solar.R + Wind + Temp + Month + Day, data = airq)

## gaussian MOB
}
