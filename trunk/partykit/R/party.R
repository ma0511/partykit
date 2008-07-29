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
    stopifnot(length(i) == 1 & is.numeric(i))
    ## FIXME: extract subtree, subset fitted, relabel fitted$"(fitted)"
    x
}


nodeids <- function(party, from = NULL, terminal = TRUE) {

    if(is.null(from)) from <- id_node(node_party(party))

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

    return(rid(node_party(party), from == id_node(node_party(party)), terminal))
}

## FIXME: remove when [[.party is available
nodeapply <- function(party, ids = 1, FUN = NULL, by_node = TRUE, ...) {

    stopifnot(inherits(party, "party"))
    stopifnot(isTRUE(all.equal(ids, round(ids))))
    ids <- as.integer(ids)

    if(is.null(FUN)) FUN <- function(x, ...) x

    if (length(ids) == 0)
        return(NULL)

    if (!by_node) {
        rval <- lapply(ids, function(i) FUN(party, i, ...))
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

predict_party <- function(party, id, newdata = NULL)
    UseMethod("predict_party")

print_party <- function(party, id, ...)
    UseMethod("print_party")

data_party <- function(party, id)
    UseMethod("data_party")

predict_party.default <- function(party, id, newdata = NULL) {

}

print_party.default <- function(party, id, newdata = NULL) {

}

data_party.default <- function(party, id) {
    ## FIXME: merge "data" and "fitted" (if possible)
    ## extract data/fitted for node "id"
}

