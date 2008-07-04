
new_node <- function(id, split = NULL, kids = NULL, ...) {

    if (!is.integer(id) || length(id) != 1)
        stop(sQuote("id"), " ", "is not an integer")

    if (!is.null(split)) {
        if (!is.list(split) || any(!sapply(split, inherits, "split")))
            stop(sQuote("split", " ", "is not a list of", " ", sQuote("split"), 
                        " ", "objects")
    }

    if (!(is.integer(kids) | 
          (is.list(kids) && all(sapply(kids, inherits, "node")))) 
        || length(kids) < 2)
        stop(sQuote("kids"), " ", "must be an integer vector or a list of", 
             " ", sQuote("node"), " ", "objects")

    node <- list(id = id, split = split, kids = kids, info = list(...))
    class(node) <- "node"
    return(node)
}

