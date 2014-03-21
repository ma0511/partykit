meanvarfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
  ..., estfun = FALSE, object = FALSE)
{
  m <- mean(y * weights)/mean(weights)
  res <- y - m
  s2 <- mean(res^2 * weights)/mean(weights)

  rval <- list(
    coefficients = structure(c(m, log(s2)), .Names = c("mean", "log(variance)")),
    objfun = -sum(weights * dnorm(y, mean = m, sd = sqrt(s2), log = TRUE)),
    estfun = if(estfun) cbind(res, (res^2 - s2)/2) * (1/s2) * weights else NULL,
    object = NULL
  )
}

meanvartree <- function(formula, data, na.action, weights, subset, minsplit = 10L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(minsplit = minsplit, ...)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  stopifnot(require("Formula"))
  m$formula <- formula(as.Formula(formula, ~ 1), rhs = 2:1)
  m$fit <- meanvarfit
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

mixforest <- function(formula, data, na.action, weights, subset,
  ntree = 500L, mtry = 3L, minsplit = 10L, replace = TRUE, fraction = 0.632, ...)
{
  ## call and formula
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  ## call model.frame()
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  n <- nrow(mf)
  
  ## extract terms, model matrices, response, weights
  mt <- terms(formula, data = data)
  y <- mf[[1L]]
  z <- mf[, -1L, drop = FALSE]
  w <- model.weights(mf)
  if(is.null(w)) w <- rep.int(1L, n)
  z[["(weights)"]] <- NULL

  ## resampling weights (hmm, there must be a more elegant solution...)
  rw <- matrix(0L, nrow = n, ncol = ntree)
  for(i in 1:ntree) {
    tab <- table(sample(1:n, n, replace = TRUE))
    rw[as.numeric(names(tab)), i] <- as.integer(tab)
  }

  ## fit trees
  nodes <- lapply(1:ntree, function(i) partykit:::mob_partynode(Y = y, Z = z, weights = w * rw[, i],
    fit = meanvarfit, control = mob_control(mtry = mtry, minsplit = minsplit, ...), nyx = 1L))

  list(
    nodes = nodes,
    data = mf,
    weights = rw
  )
}

if(FALSE) {
library("partykit")

## single tree
airq <- subset(airquality, !is.na(Ozone) & !is.na(Solar.R))
mv <- meanvartree(Ozone ~ Solar.R + Wind + Temp + Month + Day, data = airq)
plot(mv)
coef(mv)

coef(mv)[as.character(predict(mv, type = "node")), ]


## forest
data("BostonHousing", package = "mlbench")
BostonHousing <- transform(BostonHousing,
  chas = factor(chas, levels = 0:1, labels = c("no", "yes")),
  rad = factor(rad, ordered = TRUE))

set.seed(1)
bh <- mixforest(medv ~ lstat + rm + zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  data = BostonHousing)

plot(party(bh$nodes[[1]], bh$data),
  tp_args = list(FUN = function(i) paste("n =", i$nobs)))
plot(party(bh$nodes[[2]], bh$data),
  tp_args = list(FUN = function(i) paste("n =", i$nobs)))
}
