node <- function(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL) {

    if (!is.integer(id) || length(id) != 1)
        stop(sQuote("id"), " ", "must be a single integer")

    if (is.null(split) != is.null(kids)) {
        stop(sQuote("split"), " ", "and", " ", sQuote("kids"), " ", 
             "must either both be specified or unspecified")
    }

    if (!is.null(kids)) {
        if (!(is.integer(kids) | 
              (is.list(kids) && all(sapply(kids, inherits, "node")))) 
            || length(kids) < 2)
            stop(sQuote("kids"), " ", "must be an integer vector or a list of", 
                 " ", sQuote("node"), " ", "objects")
    }

    if (!is.null(surrogates)) {
        if (!is.list(surrogates) || any(!sapply(surrogates, inherits, "split")))
            stop(sQuote("split"), " ", "is not a list of", " ", sQuote("split"), 
                 " ", "objects")
    }

    node <- list(id = id, split = split, kids = kids, surrogates = surrogates, 	info = info)
    class(node) <- "node"
    return(node)
}

id_node <- function(node) {
    stopifnot(inherits(node, "node"))
    node$id
}

kids_node <- function(node) {
    stopifnot(inherits(node, "node"))
    node$kids
}

info_node <- function(obj) {
    stopifnot(inherits(obj, "node"))
    obj$info
}

kidids_node <- function(node, data, vmatch = 1:ncol(data), obs = NULL) {

    primary <- split_node(node)
    surrogates <- surrogates_node(node)

    ### perform primary split
    x <- kidids_split(primary, data, vmatch, obs)

    ### surrogate / random splits if needed
    if (any(is.na(x))) {
        ### surrogate splits
        if (length(surrogates) >= 1) {
            for (surr in surrogates) {
                nax <- is.na(x)
                if (!any(nax)) break;
                x[nax] <- kidids_split(surr, data, vmatch, obs = obs[nax])
            }
        }
        nax <- is.na(x)
        ### random splits
        if (any(nax)) {
            prob <- prob_split(primary)
            x[nax] <- sample(1:length(prob), sum(nax), prob = prob)
        }
    }
    return(x)
}

fitted_node <- function(node, data, vmatch = 1:ncol(data), obs = 1:nrow(data)) {
    if (is.terminal(node))
        return(rep(id_node(node), length(obs)))
    retid <- nextid <- kidids_node(node, data, vmatch, obs)
    for (i in unique(nextid))
        retid[nextid == i] <- fitted_node(kids_node(node)[[i]], data,
                                          vmatch, obs[nextid == i])
    return(retid)
}


length.node <- function(x)
    length(kids_node(x))

"[.node" <- "[[.node" <- function(x, i, ...) {
    stopifnot(length(i) == 1 & is.numeric(i))
    kids_node(x)[[i]]
}

split_node <- function(node) {
    stopifnot(inherits(node, "node"))
    node$split
}

surrogates_node <- function(node) {
    stopifnot(inherits(node, "node"))
    node$surrogates
}

is.terminal <- function(x, ...)
    UseMethod("is.terminal")

is.terminal.node <- function(x, ...) {
    kids <- is.null(kids_node(x))
    split <- is.null(split_node(x))
    stopifnot(kids == split)
    kids
}

depth <- function(x, ...)
    UseMethod("depth")

depth.node <- function(x) {
    if (is.terminal(x)) return(1)
    max(sapply(kids_node(x), depth)) + 1
}

width <- function(x, ...)
    UseMethod("width")

width.node <- function(x) {
    if (is.terminal(x)) return(1)
    sum(sapply(kids_node(x), nterminal))
}

