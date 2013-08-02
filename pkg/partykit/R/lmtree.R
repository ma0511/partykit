## simple wrapper function to specify fitter and return class
lmtree <- function(formula, data, subset, na.action, weights, offset, ...)
{
  ## use dots for setting up mob_control
  control <- mob_control(...)

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- lmfit
  m$control <- control
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("lmtree", class(rval))
  return(rval)
}

## actual fitting function for mob()
lmfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  ## call lm fitting function
  if(is.null(weights) || identical(as.numeric(weights), rep.int(1, length(weights)))) {
    z <- lm.fit(x, y, offset = offset, ...)
  } else {
    weights <- 1
    z <- lm.wfit(x, y, weights = weights, offset = offset, ...)
  }

  ## list structure
  rval <- list(
    coefficients = z$coefficients,
    objfun = sum(z$residuals^2),
    estfun = NULL,
    object = NULL
  )

  ## add estimating functions (if desired)
  if(estfun) {
    rval$estfun <- as.vector(z$residuals) * weights * x
  }

  ## add model (if desired)
  if(object) {
    class(z) <- c(if(is.matrix(z$fitted)) "mlm", "lm")
    z$offset <- if(is.null(offset)) 0 else offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- attr(x, "xlevels")    

    cl <- as.call(expression(lm))
    cl$formula <- attr(x, "formula")	
    z$call <- cl
    z$terms <- attr(x, "terms")

    rval$object <- z
  }

  return(rval)
}

## methods
print.lmtree <- function(x,
  title = "Linear model tree", objfun = "residual sum of squares", ...)
{
  print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.lmtree <- function(object, newdata = NULL, type = "response", ...)
{
  ## FIXME: possible to get default?
  if(is.null(newdata) & !identical(type, "node")) stop("newdata has to be provided")
  predict.modelparty(object, newdata = newdata, type = type, ...)
}

plot.lmtree <- function(x, terminal_panel = node_bivplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}



if(FALSE) {

data("Journals", package = "AER")
jour <- Journals
jour$age <- 2000 - Journals$foundingyear
jour$citeprice <- with(Journals, price/citations)
jour$lsubs <- log(jour$subs)
jour$lciteprice <- log(jour$citeprice)

mb_jour <- lmtree(log(subs) ~ log(citeprice) | price + citations + age + charpp + society,
  data = jour, minsplit = 20)

myfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) lm(y ~ 0 + x)
mb <- mob(log(subs) ~ log(citeprice) | price + citations + age + charpp + society,
  data = jour, fit = myfit)

data("BostonHousing", package = "mlbench")
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)

mb_bh <- lmtree(medv ~ log(lstat) + I(rm^2) | zn +
  indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  data = BostonHousing, minsplit = 40)

mb_bh
plot(mb_bh)
coef(mb_bh)
summary(mb_bh, node = 7)
sctest(mb_bh, node = 7)
mean((BostonHousing$medv - fitted(mb_bh))^2)
mean(residuals(mb_bh)^2)
deviance(mb_bh)/sum(weights(mb_bh))
deviance(mb_bh)/nobs(mb_bh)
logLik(mb_bh)
AIC(mb_bh)

}
