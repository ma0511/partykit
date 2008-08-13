## FIXME: data in party
##   - currently assumed to be a data.frame
##   - potentially empty
##   - the following are all assumed to work:
##     dim(data), names(data)
##     sapply(data, class), lapply(data, levels)
##   - potentially these need to be modified if data/terms
##     should be able to deal with data bases

party <- function(node, data, fitted, terms = NULL, names = NULL) {

    stopifnot(inherits(node, "node"))
    stopifnot(inherits(data, "data.frame"))
    
    if(!is.null(fitted)) {
        stopifnot(inherits(fitted, "data.frame"))
        stopifnot("(fitted)" == names(fitted)[1])
        stopifnot(nrow(data) == 0 | nrow(data) == nrow(fitted))

        nt <- nodeids(node, terminal = TRUE)
        stopifnot(all(fitted[["(fitted)"]] %in% nt))

        node <- as.node(node, from = 1L)
        nt2 <- nodeids(node, terminal = TRUE)
        fitted[["(fitted)"]] <- nt2[match(fitted[["(fitted)"]], nt)]
    } else {
        node <- as.node(node, from = 1L)
    }
    
    party <- list(node = node, data = data, fitted = fitted, 
                  terms = NULL, names = NULL)
    class(party) <- "party"

    if(!is.null(terms)) {
        stopifnot(inherits(terms, "terms"))
	party$terms <- terms
    }

    if (!is.null(names)) {
        n <- length(nodeids(party, terminal = FALSE))
        if (length(names) != n)
            stop("invalid", " ", sQuote("names"), " ", "argument")
        party$names <- names
    }

    party
}

length.party <- function(x)
    length(nodeids(x))

names.party <- function(x)
    x$names

"names<-.party" <- function(x, value) {
     n <- length(nodeids(x, terminal = FALSE))
     if (!is.null(value) && length(value) != n)
         stop("invalid", " ", sQuote("names"), " ", "argument")
     x$names <- value
     x
}

names_party <- function(party) {
    names <- party$names
    if (is.null(names))
        names <- as.character(nodeids(party, terminal = FALSE))
    names
}

node_party <- function(party) {
    stopifnot(inherits(party, "party"))
    party$node
}

"[.party" <- "[[.party" <- function(x, i, ...) {
    if (is.character(i) && !is.null(names(x)))
        i <- which(names(x) %in% i)
    stopifnot(length(i) == 1 & is.numeric(i))
    stopifnot(i <= length(x) & i >= 1)
    i <- as.integer(i)
    dat <- data_party(x, i)
    if (!is.null(x$fitted)) {
        findx <- which("(fitted)" == names(dat))[1]
        fit <- dat[,findx:ncol(dat), drop = FALSE]
        dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
        if (ncol(dat) == 0)
            dat <- x$data
    } else {
        fit <- NULL
        dat <- x$data
    }
    nam <- names(x)[nodeids(x, from = i, terminal = FALSE)]

    recFun <- function(node) {
        if (id_node(node) == i) return(node)
        kid <- sapply(kids_node(node), id_node)
        return(recFun(node[[max(which(kid <= i))]]))
    }
    node <- recFun(node_party(x))

    ret <- party(node = node, data = dat, fitted = fit, terms = x$terms, names = nam)
    class(ret) <- class(x)
    ret
}

nodeids <- function(obj, ...)
    UseMethod("nodeids")

nodeids.node <- function(obj, from = NULL, terminal = FALSE, ...) {

    if(is.null(from)) from <- id_node(obj)

    id <- function(node, record = TRUE, terminal = FALSE) {
      if(!record) return(NULL)
      if(!terminal)
          return(id_node(node))
      else
          if(is.terminal(node)) return(id_node(node)) else return(NULL)
    }

    rid <- function(node, record = TRUE, terminal = FALSE) {  
        myid <- id(node, record = record, terminal = terminal)
        if(is.terminal(node)) return(myid)
        kids <- kids_node(node)    
        kids_record <- if(record)  
            rep(TRUE, length(kids))
        else
            sapply(kids, id_node) == from
        return(c(myid,
            unlist(lapply(1:length(kids), function(i)
                rid(kids[[i]], record = kids_record[i], terminal = terminal)))
        ))
    }

    return(rid(obj, from == id_node(obj), terminal))
}

nodeids.party <- function(obj, from = NULL, terminal = FALSE, ...)
    nodeids(node_party(obj), from = from, terminal = terminal, ...)


nodeapply <- function(party, ids = 1, FUN = NULL, by_node = TRUE, ...) {

    stopifnot(inherits(party, "party"))
    stopifnot(isTRUE(all.equal(ids, round(ids))))
    ids <- as.integer(ids)

    if(is.null(FUN)) FUN <- function(x, ...) x

    if (length(ids) == 0)
        return(NULL)

    if (!by_node) {
        rval <- lapply(ids, function(i) FUN(party[[i]], ...))
    } else {
        rval <- vector(mode = "list", length = length(ids))
	rval_id <- rep(0, length(ids))
	i <- 1
	
	recFUN <- function(node, ...) {
	    if(id_node(node) %in% ids) {
	        rval_id[i] <<- id_node(node)
	        rval[[i]] <<- FUN(node, ...)
	        i <<- i + 1
	    }
	    kids <- kids_node(node)
	    if(length(kids) > 0) {
	        for(j in 1:length(kids)) recFUN(kids[[j]])
	    }
	    	  
	    invisible(TRUE)
	}
	foo <- recFUN(node_party(party))
	
        rval <- rval[match(rval_id, ids)]
    }

    names(rval) <- names_party(party)[ids]
    return(rval)
}


predict.party <- function(object, newdata = NULL, ...)
{
    fitted <- if(is.null(newdata)) object$fitted[["(fitted)"]] else {

        terminal <- nodeids(object, terminal = TRUE)
        inner <- 1:max(terminal)
        inner <- inner[-terminal]

        primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
            varid_split(split_node(node))
        })
        surrogate_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
            surr <- surrogates_node(node)
            if(is.null(surr)) return(NULL) else return(sapply(surr, varid_split))
        })
        vnames <- names(object$data)

        ## ## FIXME: the is.na() call takes loooong on large data sets
        ## unames <- if(any(sapply(newdata, is.na))) 
        ##     vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
        ## else 
        ##     vnames[unique(unlist(primary_vars))]
	unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
	
        vclass <- structure(lapply(object$data, class), .Names = vnames)
        ndnames <- names(newdata)
        ndclass <- structure(lapply(newdata, class), .Names = ndnames)
        if(all(unames %in% ndnames) &&
           all(unlist(lapply(unames, function(x) vclass[[x]] == ndclass[[x]])))) {
            vmatch <- match(vnames, ndnames)
            fitted_node(node_party(object), newdata, vmatch)
        } else {
            if (!is.null(object$terms)) {
                mf <- model.frame(delete.response(object$terms), newdata)
                fitted_node(node_party(object), mf, match(vnames, names(mf)))
            } else
                stop("") ## FIXME: write error message
        }
    }
    predict_party(object, fitted, newdata, ...)
}

predict_party <- function(party, id, newdata = NULL, ...)
    UseMethod("predict_party")

predict_party.default <- function(party, id, newdata = NULL, ...) {
  ## FIXME: improve warning message
  if(length(as.list(...)) > 0) warning("too many arguments") 
  return(structure(id, .Names = if(is.null(newdata)) 
      rownames(object$fitted) else rownames(newdata)))
}

### response: at scale of the response
### prob: discrete densities, KM for survival
predict_party.cparty <- function(object, id, newdata = NULL,
    type = c("response", "prob", "node"), FUN = NULL, simplify = TRUE, ...)
{
    ## extract fitted information
    response <- object$fitted[["(response)"]]
    weights <- object$fitted[["(weights)"]]
    fitted <- object$fitted[["(fitted)"]]
    if (is.null(weights)) weights <- rep(1, NROW(response))

    ## get observation names: either node names or
    ## observation names from newdata
    nam <- if(is.null(newdata)) names_party(object)[id] else rownames(newdata)
    if(length(nam) != length(id)) nam <- NULL

    ## match type
    type <- match.arg(type)

    ## special case: fitted ids
    if(type == "node")
      return(structure(id, .Names = nam))

    if (is.null(FUN)) {
        ### FIXME: multivariate response
#        if (is.list(response)) return(lapply(response, function(r) {
#            object$fitted[["(response)"]] <- r
#            predict_party(object, id = id, newdata = newdata, type = type, 
#                          simplify = simplify, ...)
#        }))

        rtype <- class(response)[1]
        if (rtype == "ordered") rtype <- "factor"    
        if (rtype == "integer") rtype <- "numeric"

        FUN <- switch(rtype,
                    "Surv" = function(y, w) {
                        sf <- survival:::survfit(y, weights = w, subset = w > 0)
                        if (type == "response")
                            sf <- party:::mst(sf) ### FIXME: Surv(mst)???, copy
                        sf
                    },
                    "factor" = function(y, w) {
                        sumw <- tapply(w, y, sum)
                        sumw[is.na(sumw)] <- 0
                        prob <- sumw / sum(w)
                        names(prob) <- levels(response)
                        if (type == "response")
                            return(factor(which.max(prob), levels = 1:nlevels(response),
                                          labels = levels(response), 
                                          ordered = is.ordered(response)))
                        return(prob)
                    },
                    "numeric" = {
                        if (type == "prob")
                            stop(sQuote("type = \"prob\""), " ", "is not available")
                        function(y, w)
                            weighted.mean(y, w)
                    })
    }
      
    ## empirical distribution in each leaf
    if (all(id %in% fitted)) {
        tab <- tapply(1:NROW(response), fitted, 
                      function(i) FUN(response[i], weights[i]), simplify = FALSE)
    } else {
        ### id may also refer to inner nodes
        tab <- as.array(lapply(sort(unique(id)), function(i) {
            index <- fitted %in% nodeids(object, i, terminal = TRUE)
            FUN(response[index], weights[index])
        }))
        names(tab) <- as.character(sort(unique(id)))
    }
    tn <- names(tab)
    dim(tab) <- NULL
    names(tab) <- tn

    if (simplify) {
        if (all(sapply(tab, length) == 1) & all(sapply(tab, is.atomic))) {
            ret <- do.call("c", tab)
            names(ret) <- names(tab)
            ret <- if (is.factor(tab[[1]]))
                factor(ret[as.character(id)], levels = 1:length(levels(tab[[1]])),
		       labels = levels(tab[[1]]), ordered = is.ordered(tab[[1]]))
            else 
                ret[as.character(id)]
            names(ret) <- nam
        } else if (length(unique(sapply(tab, length))) == 1 & 
                   all(sapply(tab, is.numeric))) {
            ret <- matrix(unlist(tab), nrow = length(tab), byrow = TRUE)
            colnames(ret) <- names(tab[[1]])
            rownames(ret) <- names(tab)
            ret <- ret[as.character(id),, drop = FALSE]
            rownames(ret) <- nam
        } else {
            ret <- tab[as.character(id)]
            names(ret) <- nam
        }
    } else {
        ret <- tab[as.character(id)]
        names(ret) <- nam
    }
    
    return(ret)
}

print_party <- function(party, id, ...)
    UseMethod("print_party")

data_party <- function(party, id)
    UseMethod("data_party")

print_party.default <- function(party, id, newdata = NULL) {

}

data_party.default <- function(party, id) {
    
    extract <- function(id) {
        if(is.null(party$fitted))
            if(nrow(party$data) == 0) return(NULL)
        else
            stop("cannot subset data without fitted ids")

        ### which terminal nodes follow node number id?
        nt <- nodeids(party, id, terminal = TRUE)
        wi <- party$fitted[["(fitted)"]] %in% nt

        ret <- if (nrow(party$data) == 0)
            subset(party$fitted, wi)
        else
            subset(cbind(party$data, party$fitted), wi)
        ret
    }
    if (length(id) > 1)
        return(lapply(id, extract))
    else 
        return(extract(id))
}
