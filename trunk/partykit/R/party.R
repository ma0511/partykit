new_party <- function(node, metadata, names = NULL, info = NULL) {

    stopifnot(inherits(node, "node"))
    stopifnot(inherits(metadata, "metadata"))

    if (is.flat(node))
        stop(sQuote("node"), " ", "is not a recursive node")

    party <- list(node = node, metadata = metadata, 
                  names = NULL, info = info)
    class(party) <- "party"

    if (!is.null(names)) {
        n <- length(nodeids(party, terminal = FALSE))
        if (length(names) != n)
            stop("invalid", " ", sQuote("names"), " ", "argument")
        party$names <- names
    }

    party
}

names.party <- function(x)
    x$names

"names<-.party" <- function(x, value) {
     n <- length(nodeids(x, terminal = FALSE))
     if (!is.null(value) && length(value) != n)
         stop("invalid", " ", sQuote("names"), " ", "argument")
     x$names <- value
     x
}

get_names <- function(party) {
    names <- party$names
    if (is.null(names))
        names <- as.character(nodeids(party, terminal = FALSE))
    names
}

nodeids <- function(party, from = 1, terminal = TRUE) {

    id <- function(node, record = TRUE, terminal = FALSE) {
      if(!record) return(NULL)
      if(!terminal)
          return(get_id(node))
      else
          if(is.terminal(node)) return(get_id(node)) else return(NULL)
    }

    rid <- function(node, record = TRUE, terminal = FALSE) {
        myid <- id(node, record = record, terminal = terminal)
        if(is.terminal(node)) return(myid)
        kids <- get_kids(node)
        kids_record <- if(record)
            rep(TRUE, length(kids))
        else
            sapply(kids, get_id) == from
        return(c(myid,
            unlist(lapply(1:length(kids), function(i)
	        rid(kids[[i]], record = kids_record[i], terminal = terminal)))
        ))
    }

    return(rid(party$node, from == 1, terminal))
}

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
	    if(get_id(node) %in% ids) {
	        rval_id[i] <<- get_id(node)
	        rval[[i]] <<- FUN(node, ...)
	        i <<- i + 1
	    }
	    kids <- get_kids(node)
	    if(length(kids) > 0) {
	        for(j in 1:length(kids)) recFUN(kids[[j]])
	    }
	    	  
	    invisible(TRUE)
	}
	foo <- recFUN(party$node)
	
        rval <- rval[match(rval_id, ids)]
    }

    names(rval) <- get_names(party)[ids]
    return(rval)
}
