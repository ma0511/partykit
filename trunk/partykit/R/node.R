
new_node <- function(id, split = NULL, kids = NULL, info = NULL) {

    if (!is.integer(id) || length(id) != 1)
        stop(sQuote("id"), " ", "is not an integer")

    if (is.null(split) != is.null(kids)) {
        stop(sQuote("split"), " ", "and", " ", sQuote("kids"), " ", 
             "must either both be specified or unspecified")
    }

    if (!is.null(split)) {
        if (!is.list(split) || any(!sapply(split, inherits, "split")))
            stop(sQuote("split"), " ", "is not a list of", " ", sQuote("split"), 
                        " ", "objects")
    }

    if (!is.null(kids)) {
        if (!(is.integer(kids) | 
              (is.list(kids) && all(sapply(kids, inherits, "node")))) 
            || length(kids) < 2)
            stop(sQuote("kids"), " ", "must be an integer vector or a list of", 
                 " ", sQuote("node"), " ", "objects")
    }

    node <- list(id = id, split = split, kids = kids, info = info)
    class(node) <- "node"
    return(node)
}

flat2rec <- function(obj) {

    if (!all(sapply(obj, inherits, "node"))
        || !all(sapply(obj, function(node) is.null(node$kids) | is.integer(node$kids))))
        stop(sQuote("obj"), " ", "is not a list of flat", " ", 
             sQuote("node"), " ", "objects")
    
    obj <- obj[order(sapply(obj, function(node) node$id))]
    if (length(obj) == 1) return(obj)

    new_recnode <- function(id) {
        if (is.null(obj[[id]]$kids))
            new_node(id = id, info = obj[[id]]$info)
        else
            new_node(id = id, split = obj[[id]]$split, 
                     kids = lapply(obj[[id]]$kids, new_recnode), 
                     info = obj[[id]]$info)
    }
        
    node <- new_node(id = as.integer(1), split = obj[[1]]$split, 
                     kids = lapply(obj[[1]]$kids, new_recnode), info = obj[[1]]$info)
    return(node)
}

get_node_id <- function(node, data) {

    if (is.null(node$split))
        return(rep(node$id, nrow(data)))
    retid <- nextid <- do_splitlist(node$split, data)
    for (i in unique(nextid)) {
        retid[nextid == i] <- get_node_id(node$kids[[i]], data[nextid == i, , drop = FALSE])
    }
    return(retid)
}

rec2flat <- function(node) {

    if (is.null(node$split) && is.null(node$kids))
        return(list(node))

    obj <- list()
    
    new_flatnode <- function(node) {
        if (is.null(node$kids))
            obj[[node$id]] <<- new_node(id = node$id, info = node$info)
        else {
            obj[[node$id]] <<- new_node(id = node$id, split = node$split,
                 kids = sapply(node$kids, function(k) k$id), info = node$info)
            lapply(node$kids, new_flatnode)
        }
    }
    new_flatnode(node)
    return(obj)
}

length.node <- function(x)
    length(x$kids)


"[.node" <- "[[.node" <- function(x, i, ...) {
    stopifnot(length(i) == 1 & is.numeric(i))
    rval <- x$kids[[i]]
    if (!inherits(rval, "node"))
        warning(sQuote("x"), " ", "is not a recursive node")
    return(rval)
}

