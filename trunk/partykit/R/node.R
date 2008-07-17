
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

get_id <- function(node) {
    stopifnot(inherits(node, "node"))
    node$id
}

get_split <- function(node) {
    stopifnot(inherits(node, "node"))
    node$split
}

get_kids <- function(node) {
    stopifnot(inherits(node, "node"))
    node$kids
}

is.terminal <- function(node) {
    kids <- is.null(get_kids(node))
    split <- is.null(get_split(node))
    stopifnot(kids == split)
    kids
}

is.flat <- function(node)
    is.terminal(node) | is.integer(get_kids(node))


flat2rec <- function(obj) {

    if (!all(sapply(obj, is.flat)))
        stop(sQuote("obj"), " ", "is not a list of flat", " ", 
             sQuote("node"), " ", "objects")
    
    obj <- obj[order(sapply(obj, function(node) get_id(node)))]
    if (length(obj) == 1) return(obj)

    new_recnode <- function(id) {
        if (is.null(get_kids(obj[[id]])))
            new_node(id = id, info = obj[[id]]$info)
        else
            new_node(id = id, split = get_split(obj[[id]]),
                     kids = lapply(get_kids(obj[[id]]), new_recnode), 
                     info = obj[[id]]$info)
    }
        
    node <- new_node(id = as.integer(1), split = get_split(obj[[1]]),
                     kids = lapply(get_kids(obj[[1]]), new_recnode), 
                     info = obj[[1]]$info)
    return(node)
}

do_nodeid <- function(node, data) {

    if (is.terminal(node))
        return(rep(get_id(node), nrow(data)))
    retid <- nextid <- do_splitlist(get_split(node), data)
    for (i in unique(nextid))
        retid[nextid == i] <- do_nodeid(get_kids(node)[[i]], 
                                          data[nextid == i, , drop = FALSE])
    return(retid)
}

rec2flat <- function(node) {

    if (is.terminal(node))
        return(list(node))

    if (is.flat(node))
        stop(sQuote("node"), " ", "is not a recursive node")

    obj <- list()
    
    new_flatnode <- function(node) {
        if (is.terminal(node))
            obj[[node$id]] <<- new_node(id = get_id(node), info = node$info)
        else {
            obj[[node$id]] <<- new_node(id = get_id(node), split = get_split(node),
                 kids = sapply(get_kids(node), function(k) get_id(k)), info = node$info)
            lapply(get_kids(node), new_flatnode)
        }
    }
    new_flatnode(node)
    return(obj)
}

length.node <- function(x)
    length(get_kids(x))


"[.node" <- "[[.node" <- function(x, i, ...) {
    if (is.flat(x))
        warning(sQuote("x"), " ", "is not a recursive node")
    stopifnot(length(i) == 1 & is.numeric(i))
    get_kids(x)[[i]]
}
