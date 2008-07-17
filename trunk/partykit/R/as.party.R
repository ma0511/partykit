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
        metadata <- list(varnames = names(attr(obj$terms, "dataClasses")),
                         class = as.vector(attr(obj$terms, "dataClasses")))
        ylev <- list(attr(obj, "ylevels"))
        names(ylev) <- metadata$varnames[1]
        lev <- c(attr(obj, "xlevels"), ylev)
        metadata$levels <- lapply(metadata$varnames, function(var) lev[[var]])
        class(metadata) <- "metadata"
        metadata
    }
    objmeta <- rpart_metadata()

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

    new_party(node = node, metadata = objmeta)
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

as.party.J48 <- function(obj, ...) {

  ## construct metadata
  meta <- metadata(model.frame(obj))

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

    if(meta$class[var_id] %in% c("ordered", "factor")) { ## FIXME: ordered splits should be handled differently
      stopifnot(all(sapply(split, head, 1) == "="))
      stopifnot(all(sapply(split, tail, 1) %in% meta$levels[[var_id]]))
      
      split <- new_split(fun = as.integer(var_id), index = match(meta$levels[[var_id]], sapply(split, tail, 1)))
    } else {
      breaks <- unique(as.numeric(sapply(split, tail, 1)))
      breaks <- if(meta$class[var_id] == "integer") as.integer(breaks) else as.double(breaks) ## FIXME: check
      
      stopifnot(length(breaks) == 1 && !is.na(breaks))
      stopifnot(all(sapply(split, head, 1) %in% c("<=", ">")))
      
      split <- new_split(fun = as.integer(var_id), breaks = breaks, right = TRUE, index = if(split[[1]][1] == ">") 2:1)
    }
    return(split)
  }

  j48_node <- function(i) {
    if(is.null(j48_kids(i))) return(new_node(as.integer(i)))
    new_node(as.integer(i), split = list(j48_split(i)), kids = lapply(j48_kids(i), j48_node))
  }

  node <- j48_node(1)

  ## FIXME: node IDs
  j48 <- new_party(node = node, metadata = meta)
  return(j48)
}

## FIXME: small convenience function (just temporary)
j48fit <- function(fit) {
  mf <- model.frame(fit)
  fit <- as.party(fit)

  j48summary <- function(x)
    paste(levels(x)[which.max(table(x))], " (", length(x), "/", length(x) - max(table(x)), ")", sep = "")

  tapply(mf[,1], get_node_id(fit$node, mf), j48summary)
}


