mtree_ols <- list(
  fit = function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
    ..., estfun = FALSE, object = FALSE)
  {
    if(!is.null(offset)) y <- y - offset
    m <- mean(y * weights)/mean(weights)
    res <- y - m
    rss <- sum(res^2 * weights)
    
    list(
      coefficients = c("mean" = m),
      objfun = rss,
      estfun = if(estfun) res * weights else NULL,
      object = if(object) list(nuisance = c("log(variance)" = log(rss/sum(weights)))) else NULL
    )
  },
  
  ddist = function(x, par, log = FALSE) {
    par <- matrix(par, ncol = 2L)
    dnorm(x, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), log = log)
  },

  pdist = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    par <- matrix(par, ncol = 2L)
    pnorm(q, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), lower.tail = lower.tail, log = log)
  },

  qdist = function(p, par, lower.tail = TRUE, log.p = FALSE) {
    par <- matrix(par, ncol = 2L)
    qnorm(p, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), lower.tail = lower.tail, log = log)
  }
)

mtree_norm <- list(
  fit = function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
    ..., estfun = FALSE, object = FALSE)
  {
    if(!is.null(offset)) y <- y - offset
    m <- mean(y * weights)/mean(weights)
    res <- y - m
    s2 <- mean(res^2 * weights)/mean(weights)

    list(
      coefficients = structure(c(m, log(s2)), .Names = c("mean", "log(variance)")),
      objfun = -sum(weights * dnorm(y, mean = m, sd = sqrt(s2), log = TRUE)),
      estfun = if(estfun) cbind(res, (res^2 - s2)/2) * (1/s2) * weights else NULL,
      object = NULL
    )
  },
  
  ddist = function(x, par, log = FALSE) {
    par <- matrix(par, ncol = 2L)
    dnorm(x, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), log = log)
  },

  pdist = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    par <- matrix(par, ncol = 2L)
    pnorm(q, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), lower.tail = lower.tail, log = log)
  },

  qdist = function(p, par, lower.tail = TRUE, log.p = FALSE) {
    par <- matrix(par, ncol = 2L)
    qnorm(p, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), lower.tail = lower.tail, log = log)
  }
)

mtree_multinom <- list(
  fit = function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
    ..., estfun = FALSE, object = FALSE)
  {
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
    
    return(rval)
  },
  
  ddist = function(x, par, log = FALSE) {
    par <- matrix(par, ncol = NCOL(par))
    par <- cbind(par, 1 - rowSums(par))
    dmultinom(x, prob = par, log = log)
  },

  pdist = NULL,

  qdist = NULL
)

mtree_binom <- list(
  fit = function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
    ..., estfun = FALSE, object = FALSE)
  {
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
    
    return(rval)
  },
  
  ddist = function(x, par, log = FALSE) {
    par <- matrix(par, ncol = 2L)
    dbinom(x, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), log = log)
  },

  pdist = function(q, par, lower.tail = TRUE, log.p = FALSE) {
    par <- matrix(par, ncol = 2L)
    pbinom(q, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), lower.tail = lower.tail, log = log)
  },

  qdist = function(p, par, lower.tail = TRUE, log.p = FALSE) {
    par <- matrix(par, ncol = 2L)
    qbinom(p, mean = par[, 1L], sd = sqrt(exp(par[, 2L])), lower.tail = lower.tail, log = log)
  }
)
  

mtree <- function(formula, data, na.action, weights, subset, type = NULL, minsplit = 7L, ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## use dots for setting up mob_control
  control <- mob_control(minsplit = minsplit, ...)

  ## model type
  stopifnot(!is.null(type))
  if(is.character(type)) getAnywhere(paste("mtree", type, sep = "_"))
  if(!is.list(type) || !all(c("fit", "ddist", "pdist", "qdist") %in% names(type))) {
    stop("invalid specification of 'type'")
  }

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$formula <- formula(Formula::as.Formula(formula, ~ 1), rhs = 2:1)
  m$fit <- type$fit
  m$control <- control
  m$minsplit <- NULL
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  rval$info$ddist <- type$ddist
  rval$info$pdist <- type$pdist
  rval$info$qdist <- type$qdist
  rval <- as.constparty(rval)
  class(rval) <- c("mtree", "constparty", "party", "modelparty")
  return(rval)
}

mforest <- function(formula, data, na.action, weights, subset, type = NULL,
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
  mf0 <- mf[0L, ]
  
  ## extract terms, model matrices, response, weights
  mt <- terms(formula, data = data)
  y <- mf[[1L]]
  z <- mf[, -1L, drop = FALSE]
  w <- model.weights(mf)
  if(is.null(w)) w <- rep.int(1L, n)
  z[["(weights)"]] <- NULL

  ## default model type
  if(is.null(type)) {
    if(is.factor(y)) {
      if(length(levels(y)) > 2L) {
        type <- "multinom"
      } else {
        type <- "binom"
      }
    } else {
      type <- "norm"
    }
  }
  if(is.character(type)) getAnywhere(paste("mtree", type), sep = "_")
  if(!is.list(type) || !all(c("fit", "ddist", "pdist", "qdist") %in% names(type))) {
    stop("invalid specification of 'type'")
  }

  ## resampling weights (hmm, there must be a more elegant solution...)
  rw <- matrix(0L, nrow = n, ncol = ntree)
  for(i in 1:ntree) {
    tab <- table(sample(1:n, n, replace = TRUE))
    rw[as.numeric(names(tab)), i] <- as.integer(tab)
  }

  ## fit trees
  nodes <- lapply(1:ntree, function(i) party(
    partykit:::mob_partynode(Y = y, Z = z, weights = w * rw[, i],
      fit = type$fit, control = mob_control(mtry = mtry, minsplit = minsplit, ...), nyx = 1L),
    mf0))

  rval <- list(
    nodes = nodes,
    data = mf,
    weights = rw,
    dist = list(
      ddist = type$ddist,
      pdist = type$pdist,
      qdist = type$qdist
    )
  )
  class(rval) <- "mforest"
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
