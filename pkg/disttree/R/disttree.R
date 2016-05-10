## high-level convenience interface to mob()
disttree <- function(formula, data, na.action, cluster, family = NO,
  control = mob_control(...), ocontrol = list(), ...)
{
  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## process family if necessary
  if(is.function(family)) family <- family()

  ## glue code for calling distfit() with given family in mob()
  dist_family_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
    cluster = NULL, ..., estfun = FALSE, object = FALSE)
  {
    if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
    if(!is.null(offset)) warning("offset not used")
    rval <- distfit(y, family = family, weights = weights, start = start,
      vcov. = object, estfun = estfun, ...)
    rval <- list(
      coefficients = rval$par,
      objfun = -rval$opt$value,
      estfun = if(estfun) rval$estfun else NULL,
      object = if(object) rval else NULL
    )
    return(rval)
  }

  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$fit <- dist_family_fit
  m$family <- m$ocontrol <- NULL
  for(n in names(ocontrol)) m[[n]] <- ocontrol[[n]]
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.name("mob")
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call/family/control
  rval$info$call <- cl
  rval$info$family <- family
  rval$info$ocontrol <- ocontrol
  class(rval) <- c("disttree", class(rval))
  return(rval)
}


## methods
print.disttree <- function(x,
  title = NULL, objfun = "negative log-likelihood", ...)
{
  if(is.null(title)) title <- sprintf("Distributional regression tree (%s)", x$info$family$family[2L])
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

## predict.disttree <- function(object, newdata = NULL,
##   type = c("worth", "rank", "best", "node"), ...)
## {
##   ## type of prediction
##   type <- match.arg(type)
##   
##   ## nodes can be handled directly
##   if(type == "node") return(partykit::predict.modelparty(object, newdata = newdata, type = "node", ...))
##   
##   ## get default newdata otherwise
##   if(is.null(newdata)) newdata <- model.frame(object)
##   
##   pred <- switch(type,
##     "worth" = worth,
##     "rank" = function(obj, ...) rank(-worth(obj)),
##     "best" = function(obj, ...) {
##       wrth <- worth(obj)
##       factor(names(wrth)[which.max(wrth)], levels = names(wrth))
##     }
##   )
##   partykit::predict.modelparty(object, newdata = newdata, type = pred, ...)
## }

## plot.disttree <- function(x, terminal_panel = node_histogram,
##   tp_args = list(...), tnex = NULL, drop_terminal = NULL, ...)
## {
##   if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L
##   if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
##   partykit::plot.modelparty(x, terminal_panel = terminal_panel,
##     tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
## }
