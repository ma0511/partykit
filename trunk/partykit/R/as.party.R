as.party <- function(obj, ...)
    UseMethod("as.party")

as.party.rpart <- function(obj, ...) {

    ff <- obj$frame
    n  <- nrow(ff)
    if (n==1) return(new_node(as.integer(1)))  # special case of no splits

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
    
    rpart_metadata <- function() {
        ### extract metadata from terms
        varnames <- names(attr(obj$terms, "dataClasses"))
        class <- as.vector(attr(obj$terms, "dataClasses"))
        ylev <- list(attr(obj, "ylevels"))
        names(ylev) <- metadata$varnames[1]
        lev <- c(attr(obj, "xlevels"), ylev)
        levels <- lapply(varnames, function(var) lev[[var]])
        new_metadata(varnames = varnames, class = class, levels = levels)
    }
    objmeta <- rpart_metadata()

    rpart_info <- function() {
        if (is.null(obj$y)) {
            y <- model.response(model.frame(obj))
        } else {
            y <- obj$y
            if (!is.null(attr(obj, "ylevels")))
                y <- factor(y, levels = attr(obj, "ylevels"))
        }
        list(responses = y, fitted = obj$where, terms = obj$terms)
    }
    objinfo <- rpart_info()

    rpart_kids <- function(i) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1, 
            which((cumsum(!is.leaf[-(1:i)]) + 1) == cumsum(is.leaf[-(1:i)]))[1] + 1 + i))
    }

    rpart_onesplit <- function(j) {
        if (j < 1) return(NULL)
        ### numeric
        if (abs(obj$split[j, "ncat"]) == 1) {
            ret <- new_split(fun = which(rownames(obj$split)[j] == objmeta$varnames),
                      breaks = as.double(obj$split[j, "index"]),
                      right = FALSE,
                      index = if(obj$split[j, "ncat"] > 0) 2:1)
        } else {
            index <- obj$csplit[obj$split[j, "index"],]
            index[index == 2] <- NA
            index[index == 3] <- 2
            ret <- new_split(fun = which(rownames(obj$split)[j] == objmeta$varnames),
                      index = as.integer(index))
        }
        ret
    }
                      
    rpart_splits <- function(i)
        lapply(c(splitindex$primary[i], splitindex$surrogate[[i]]), rpart_onesplit)
    
    rpart_node <- function(i) {
        if (is.null(rpart_kids(i))) return(new_node(as.integer(i)))
        new_node(as.integer(i), split = rpart_splits(i), kids = 
            lapply(rpart_kids(i), rpart_node))
    }

    node <- rpart_node(1)

    new_party(node = node, metadata = objmeta, info = objinfo)
}

## FIXME: put into RWeka
model.frame.Weka_classifier <- function(formula, ...) {
  mf <- formula$call
  mf <- mf[c(1, match(c("formula", "data", "subset", "na.action"), names(mf), 0))]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, environment(formula))
  return(mf)
}

model.frame.rpart <- function(formula, ...) {
  mf <- formula$call
  mf <- mf[c(1, match(c("formula", "data", "subset", "na.action"), names(mf), 0))]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, environment(formula))
  return(mf)
}

as.party.J48 <- function(obj, ...) {

  ## construct metadata
  mf <- model.frame(obj)
  meta <- metadata(mf)

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
    
    var_id <- which(nodes[i, "splitvar"] == meta$varnames)
    split <- strsplit(edges[nodes[i,"name"] == edges[,"from"], "label"], " ")

    if(meta$class[var_id] %in% c("ordered", "factor")) {
      stopifnot(all(sapply(split, head, 1) == "="))
      stopifnot(all(sapply(split, tail, 1) %in% meta$levels[[var_id]]))
      
      split <- new_split(fun = as.integer(var_id),
        index = match(meta$levels[[var_id]], sapply(split, tail, 1)))
    } else {
      breaks <- unique(as.numeric(sapply(split, tail, 1)))
      breaks <- if(meta$class[var_id] == "integer") as.integer(breaks) else as.double(breaks) ## FIXME: check
      
      stopifnot(length(breaks) == 1 && !is.na(breaks))
      stopifnot(all(sapply(split, head, 1) %in% c("<=", ">")))
      
      split <- new_split(fun = as.integer(var_id),
        breaks = breaks, right = TRUE,
	index = if(split[[1]][1] == ">") 2:1)
    }
    return(split)
  }

  j48_node <- function(i) {
    if(is.null(j48_kids(i))) return(new_node(as.integer(i)))
    new_node(as.integer(i), split = list(j48_split(i)), kids = lapply(j48_kids(i), j48_node))
  }

  node <- j48_node(1)

  j48 <- new_party(node = node, metadata = meta,
      info = list(responses = model.response(mf), fitted = do_nodeid(node, mf), terms = terms(obj)))

  class(j48) <- c("R48", class(j48))
  return(j48)
}

## FIXME: small convenience function (just temporary)
summary.R48 <- function(object) {
  j48summary <- function(x)
    paste(levels(x)[which.max(table(x))], " (", length(x), "/", length(x) - max(table(x)), ")", sep = "")

  info <- get_info(object)
  tapply(info$responses, info$fitted, j48summary)
}

predict.R48 <- function(object, newdata = NULL,
    type = c("response", "prob", "node"), ...)
{
  ## match type
  type <- match.arg(type)
  
  ## extract info slot (response and fitted node ids)
  info <- get_info(object)
  
  ## special case: fitted ids
  if(is.null(newdata) & type == "node")
    return(structure(info$fitted, .Names = names(info$responses)))
  
  ## empirical distribution in each leaf
  tab <- tapply(info$responses, info$fitted, function(x) prop.table(table(x)))
  tab <- t(structure(as.vector(unlist(tab)),
    .Dim = c(length(tab[[1]]), length(tab)),
    .Dimnames = list(names(tab[[1]]), names(tab))))

  ## get predicted leaf
  if(!is.null(newdata)) {

    ### get all relevant variables
    terminal <- nodeids(object, terminal = TRUE)
    inner <- 1:max(terminal)
    inner <- inner[-terminal]
    fun <- function(node)
        sapply(get_split(node), is.functional)

    ### we can't handle functional splits this way
    if (!any(unlist(nodeapply(object, ids = inner, FUN = fun)))) {

        used_vars <- nodeapply(object, ids = inner, FUN = function(node) {
            sapply(get_split(node), get_fun)
        })
        vnames <- object$metadata$varnames[unique(unlist(used_vars))]
    } else {
        vnames <- object$metadata$varnames[-1]
    }
        
    ## FIXME: Does this handle functional splits correctly?
    stopifnot(all(vnames %in% names(newdata)))

    ### determine correct matching (passed to do_nodeid)
    vmatch <- match(object$metadata$varnames, names(newdata))

    ## FIXME: Is there a better way for this?
    # newdata <- model.frame(delete.response(info$terms), newdata)
    # newdata[[object$metadata$varnames[1]]] <- FALSE
    # newdata <- newdata[, object$metadata$varnames, drop = FALSE]
    pred <- do_nodeid(object$node, newdata, vmatch)
    nam <- rownames(newdata)
  } else {
    pred <- info$fitted
    nam <- names(info$responses)
  }

  ## handle different types
  pred <- switch(type,
    "node" = pred,
    "prob" = tab[as.character(pred),],
    "response" = factor(apply(tab, 1, which.max)[as.character(pred)],
      levels = 1:ncol(tab), labels = colnames(tab))
  )  
  
  ## observation names
  if(NCOL(pred) > 1)
    rownames(pred) <- nam
  else
    names(pred) <- nam
  
  return(pred)
}
