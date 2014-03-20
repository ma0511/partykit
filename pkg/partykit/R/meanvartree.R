meanvarfit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
  ..., estfun = FALSE, object = FALSE)
{
  m <- mean(y * weights)/mean(weights)
  res <- y - m
  s2 <- mean(res^2 * weights)/mean(weights)

  rval <- list(
    coefficients = c(m, log(s2)),
    objfun = -sum(weights * dnorm(y, mean = m, sd = sqrt(s2), log = TRUE)),
    estfun = if(estfun) cbind(res, res^2 - s2) * weights else NULL,
    object = NULL
  )
}

meanvartree <- function(formula, data, na.action, weights, subset, minsplit = 7L, ...)
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

if(FALSE) {
airq <- subset(airquality, !is.na(Ozone) & !is.na(Solar.R))
mv <- meanvartree(Ozone ~ Solar.R + Wind + Temp + Month + Day, data = airq)
plot(mv)
}
