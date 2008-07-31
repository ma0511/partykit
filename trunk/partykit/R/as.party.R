as.party <- function(obj, ...)
    UseMethod("as.party")

as.party.rpart <- function(obj, ...) {

    ff <- obj$frame
    n  <- nrow(ff)
    if (n==1) return(node(as.integer(1)))  # special case of no splits

    is.leaf <- (ff$var == "<leaf>")
    vnames <- ff$var[!is.leaf]  #the variable names for the primary splits

    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
    splitindex <- list()
    splitindex$primary <- numeric(n)
    splitindex$primary[!is.leaf] <- index[c(!is.leaf, FALSE)]
    splitindex$surrogate <- lapply(1:n, function(i) {
        prim <- splitindex$primary[i]
        if (prim < 1 || ff[i, "nsurrogate"] == 0) return(NULL)
        else return(prim + ff[i, "ncompete"] + 1:ff[i, "nsurrogate"])
    })
    
    rpart_fitted <- function() {
        weights <- NULL
        dc <- attr(obj$terms, "dataClasses")
        y <- obj$y
        if (obj$method == "exp") y <- NULL
        if (is.null(y) | "(weights)" %in% names(dc)) {
            mf <- model.frame(obj)
            y <- model.response(mf)
            weights <- model.weights(mf)
        } else {
            y <- obj$y
            if (!is.null(attr(obj, "ylevels")))
                y <- factor(y, labels = attr(obj, "ylevels"))
        }
        ret <- data.frame("(fitted)" = obj$where, "(response)" = y, check.names = FALSE)
        if (!is.null(weights)) ret[["(weights)"]] <- weights
        ret
    }
    fitted <- rpart_fitted()

    rpart_modelframe0 <- function() {
        ### extract metadata from terms
        varnames <- names(attr(obj$terms, "dataClasses"))
        classes <- as.vector(attr(obj$terms, "dataClasses"))
        ylev <- list(attr(obj, "ylevels"))
        names(ylev) <- metadata$varnames[1]
        lev <- c(attr(obj, "xlevels"), ylev)
        levels <- lapply(varnames, function(var) lev[[var]])

        ### create data frame without obs from description
	## currently we handle
	stopifnot(all(classes %in% c("numeric", "integer", "factor", "ordered", "Surv")))
        
	mf <- rep(list(numeric(0)), length(varnames))
	names(mf) <- varnames
	wi <- which(classes %in% "integer")
	mf[wi] <- rep(list(integer(0)), length(wi))
	wi <- which(classes %in% "factor")
	mf[wi] <- lapply(wi, function(i) factor(integer(0),
	    levels = seq_along(levels[[i]]), labels = levels[[i]]))
	wi <- which(classes %in% "ordered")
	mf[wi] <- lapply(wi, function(i) factor(integer(0),
	    levels = seq_along(levels[[i]]), labels = levels[[i]], ordered = TRUE))
	wi <- which(classes %in% "Surv")
        ## blank Surv object
	surv <- structure(numeric(0), .Dim = c(0L, 2L), .Dimnames = list(NULL, c("time", "status")),
	                  class = "Surv", type = "right")
	mf[wi] <- rep(list(surv), length(wi))
	mf <- structure(mf, row.names = integer(0), class = "data.frame")
	mf
    }
    ### won't work because of dataClass problem in model.frames
    ### mf <- rpart_modelframe0()
    
    mf <- model.frame(obj)[0,]

    rpart_kids <- function(i) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1, 
            which((cumsum(!is.leaf[-(1:i)]) + 1) == cumsum(is.leaf[-(1:i)]))[1] + 1 + i))
    }

    rpart_onesplit <- function(j) {
        if (j < 1) return(NULL)
        ### numeric
        if (abs(obj$split[j, "ncat"]) == 1) {
            ret <- split(varid = which(rownames(obj$split)[j] == names(mf)),
                      breaks = as.double(obj$split[j, "index"]),
                      right = FALSE,
                      index = if(obj$split[j, "ncat"] > 0) 2:1)
        } else {
            index <- obj$csplit[obj$split[j, "index"],]
            index[index == 2] <- NA
            index[index == 3] <- 2
            ret <- split(varid = which(rownames(obj$split)[j] == names(mf)),
                      index = as.integer(index))
        }
        ret
    }
                      
    rpart_split <- function(i)
        rpart_onesplit(splitindex$primary[i])
    
    rpart_surrogates <- function(i)
        lapply(splitindex$surrogate[[i]], rpart_onesplit)

    rpart_node <- function(i) {
        if (is.null(rpart_kids(i))) return(node(as.integer(i)))
        node(as.integer(i), split = rpart_split(i),
	     kids = lapply(rpart_kids(i), rpart_node),
	     surrogates = rpart_surrogates(i))
    }

    node <- rpart_node(1)

    rval <- party(node = node, data = mf, fitted = fitted, terms = obj$terms)
    class(rval) <- c("cparty", class(rval))
    return(rval)
}

model.frame.rpart <- function(formula, ...) {
  mf <- formula$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L))]
  if (is.null(mf$na.action)) mf$na.action <- na.rpart
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  env <- if (is.null(environment(formula$terms))) environment(formula$terms) 
             else parent.frame()
  mf <- eval(mf, env)
  return(mf)
}

## FIXME: put into RWeka
model.frame.Weka_classifier <- model.frame.rpart

as.party.J48 <- function(obj, ...) {

  ## construct metadata
  mf <- model.frame(obj)
  mf_class <- sapply(mf, class)
  mf_class <- lapply(mf, levels)

  x <- .jcall(obj$classifier, "S", "graph")
  x <- RWeka:::parse_Weka_digraph(x, plainleaf = TRUE)
  nodes <- x$nodes
  edges <- x$edges
  is.leaf <- x$nodes[, "splitvar"] == ""

  j48_kids <- function(i) {
    if (is.leaf[i]) return(NULL)
      else return(which(nodes[,"name"] %in% edges[nodes[i,"name"] == edges[,"from"], "to"]))
  }

  j48_split <- function(i) {
    if(is.leaf[i]) return(NULL)
    
    var_id <- which(nodes[i, "splitvar"] == names(mf))
    split <- strsplit(edges[nodes[i,"name"] == edges[,"from"], "label"], " ")

    if(mf_class[var_id] %in% c("ordered", "factor")) {
      stopifnot(all(sapply(split, head, 1) == "="))
      stopifnot(all(sapply(split, tail, 1) %in% mf_levels[[var_id]]))
      
      split <- split(varid = as.integer(var_id),
        index = match(mf_levels[[var_id]], sapply(split, tail, 1)))
    } else {
      breaks <- unique(as.numeric(sapply(split, tail, 1)))
      breaks <- if(mf_class[var_id] == "integer") as.integer(breaks) else as.double(breaks) ## FIXME: check
      
      stopifnot(length(breaks) == 1 && !is.na(breaks))
      stopifnot(all(sapply(split, head, 1) %in% c("<=", ">")))
      
      split <- split(varid = as.integer(var_id),
        breaks = breaks, right = TRUE,
	index = if(split[[1]][1] == ">") 2:1)
    }
    return(split)
  }

  j48_node <- function(i) {
    if(is.null(j48_kids(i))) return(node(as.integer(i)))
    node(as.integer(i), split = j48_split(i), kids = lapply(j48_kids(i), j48_node))
  }

  node <- j48_node(1)

  j48 <- party(node = node,
               data = mf[0,],
               fitted = data.frame("(fitted)" = fitted_node(node, mf),
	                           "(response)" = model.response(mf),
				   check.names = FALSE),
               terms = obj$terms)

  class(j48) <- c("cparty", class(j48))
  return(j48)
}

## FIXME: small convenience function (just temporary)
summary.cparty <- function(object) {
  j48summary <- function(x)
    paste(levels(x)[which.max(table(x))], " (", length(x), "/", length(x) - max(table(x)), ")", sep = "")

  info <- info_node(object)
  tapply(info$responses, info$fitted, j48summary)
}

