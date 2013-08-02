## simple wrapper function to specify fitter and return class
lmtree <- function(formula, data, subset, na.action, weights, offset, ...)
{
  ## call mob and return
  m <- match.call(expand.dots = TRUE)
  m$fit <- lmfit
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  class(rval) <- c("lmtree", class(rval))
  return(rval)
}

## fitting function for mob()
lmfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...,
  estfun = FALSE, object = FALSE)
{
  ## call lm fitting function
  if(is.null(weights) || identical(as.numeric(weights), rep.int(1, length(weights)))) {
    z <- lm.fit(x, y, offset = offset)
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
print.lmtree <- function(x, ...) {
  print.modelparty(x, name = "Linear model tree", ...)
}

predict.lmtree <- function(object, newdata = NULL, type = "response", ...)
{
  ## FIXME: possible to get default?
  if(is.null(newdata) & !identical(type, "node")) stop("newdata has to be provided")
  predict.modelparty(object, newdata = newdata, type = type, ...)
}

plot.lmtree <- function(x, terminal_panel = node_bivplot, tnex = NULL, ...)
{
  if(is.null(tnex)) tnex <- as.integer(!is.null(terminal_panel)) + 1L
  plot.modelparty(x, terminal_panel = terminal_panel, tnex = tnex, ...)
}



if(FALSE) {

data("Journals", package = "AER")
jour <- Journals
jour$age <- 2000 - Journals$foundingyear
jour$citeprice <- with(Journals, price/citations)
jour$lsubs <- log(jour$subs)
jour$lciteprice <- log(jour$citeprice)

mb <- lmtree(log(subs) ~ log(citeprice) | price + citations + age + charpp + society,
  data = jour, control = mob_control(minsplit = 10))
}
