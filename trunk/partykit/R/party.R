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
    stopifnot(inherits(fitted, "data.frame"))
    stopifnot("(fitted)" %in% names(fitted))

    party <- list(node = as.node(node, from = 1L), data = data, fitted = fitted, 
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
    stopifnot(length(i) == 1 & is.numeric(i))
    ## FIXME: extract subtree, subset fitted, relabel fitted$"(fitted)"
    x
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

predict_party <- function(party, id, newdata = NULL, ...)
    UseMethod("predict_party")

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
        unames <- if(any(sapply(newdata, is.na))) vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
            else vnames[unique(unlist(primary_vars))]
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
                stop("")
        }
    }
    # return(fitted)
    predict_party(object, fitted, newdata, ...)
}

predict_party.default <- function(party, id, newdata = NULL) {

}

predict_party.cparty <- function(object, id, newdata = NULL,
    type = c("response", "prob", "node"), FUN = NULL, ...)
{
    ## match type
    type <- match.arg(type)
    if (type == "node") {
         names(id) <- rownames(newdata)
         return(id)
    }
  
    response <- object$fitted[["(response)"]]
    rname <- names(object$fitted["(response)"])
    weights <- object$fitted[["(weights)"]]
    fitted <- object$fitted[["(fitted)"]]
    if (is.null(weights)) weights <- rep(1, NROW(response))

    ## special case: fitted ids
    if(type == "node")
      return(structure(id, .Names = rname))

    myFUN <- FUN    
    if (is.null(myFUN)) {
        myFUN <- switch(class(response)[1],
                      "Surv" = function(y, w)
                          survival:::survfit(y, weights = w, subset = w > 0),
                      "factor" = function(y, w) {
                          sumw <- tapply(w, y, sum)
                          sumw[is.na(sumw)] <- 0
                          sumw / sum(w)
                      },
                      ### FIXME: really?
                      "ordered" = function(y, w) {
                          sumw <- tapply(w, y, sum)
                          sumw[is.na(sumw)] <- 0
                          sumw / sum(w)
                      },
                      "numeric" = function(y, w)
                          weighted.mean(y, w),
                      "integer" = function(y, w)
                          weighted.mean(y, w))
    }
      
    ## empirical distribution in each leaf
    tab <- tapply(1:NROW(response), fitted, function(i) myFUN(response[i], weights[i]))
    if (!is.null(FUN)) return(tab[as.character(id)])

    if (inherits(response, "Surv")) {
        ret <- tab[as.character(id)]
        names(ret) <- rownames(newdata)
        return(ret)
    }

    if (is.numeric(response)) {
        ret <- as.vector(tab[as.character(id)])
        names(ret) <- rownames(newdata)
        return(ret)
    }

    ## handle different types
    switch(type,
        "prob" = {
             ret <- matrix(unlist(tab[as.character(id)]), nrow = length(id), byrow = TRUE)
             colnames(ret) <- levels(response)
             rownames(ret) <- rownames(newdata)
        },
        "response" = {
             ret <- sapply(tab, which.max)
             names(ret) <- names(tab)
             ret <- factor(ret[as.character(id)],
                           labels = levels(response))
             names(ret) <- rownames(newdata)
        }
    )  

    return(ret)
}


print_party <- function(party, id, ...)
    UseMethod("print_party")

data_party <- function(party, id)
    UseMethod("data_party")

print_party.default <- function(party, id, newdata = NULL) {

}

data_party.default <- function(party, id) {
    ## FIXME: merge "data" and "fitted" (if possible)
    ## extract data/fitted for node "id"
}

