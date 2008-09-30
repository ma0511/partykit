
### cd ../src
### R CMD SHLIB *.c -o partykit.so
### dyn.load("../src/partykit.so")

partynode <- function(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL) {

    if (!is.integer(id) || length(id) != 1)
        stop(sQuote("id"), " ", "must be a single integer")

    if (is.null(split) != is.null(kids)) {
        stop(sQuote("split"), " ", "and", " ", sQuote("kids"), " ", 
             "must either both be specified or unspecified")
    }

    if (!is.null(kids)) {
        if (!(is.integer(kids) | 
              (is.list(kids) && all(sapply(kids, inherits, "partynode")))) 
            || length(kids) < 2)
            stop(sQuote("kids"), " ", "must be an integer vector or a list of", 
                 " ", sQuote("partynode"), " ", "objects")
    }

    if (!is.null(surrogates)) {
        if (!is.list(surrogates) || any(!sapply(surrogates, inherits, "partysplit")))
            stop(sQuote("split"), " ", "is not a list of", " ", sQuote("partysplit"), 
                 " ", "objects")
    }

    node <- list(id = id, split = split, kids = kids, surrogates = surrogates, info = info)
    class(node) <- "partynode"
    return(node)
}

is.partynode <- function(x) {
    if (!inherits(x, "partynode")) return(FALSE)
    rval <- diff(nodeids(x, terminal = FALSE))
    isTRUE(all.equal(unique(rval), 1))
}

as.partynode <- function(x, ...)
    UseMethod("as.partynode")

as.partynode.partynode <- function(x, from = NULL, ...) {
    if(is.null(from)) from <- id_node(x)
    from <- as.integer(from)
    if(is.partynode(x) & id_node(x) == from) return(x)
    id <- from - 1L
     
    new_node <- function(x) {
        id <<- id + 1L
        if(is.terminal(x)) return(partynode(id, info = info_node(x)))
        partynode(id,
             split = split_node(x),
             kids = lapply(kids_node(x), new_node),
             surrogates = surrogates_node(x),
             info = info_node(x))
    }
    
    return(new_node(x))    
}

id_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$id
}

kids_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$kids
}

info_node <- function(obj) {
    stopifnot(inherits(obj, "partynode"))
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
            x[nax] <- sample(1:length(prob), sum(nax), prob = prob, 
                             replace = TRUE)
        }
    }
    return(x)
}

fitted_node <- function(node, data, vmatch = 1:ncol(data), obs = 1:nrow(data)) {

    ### should be equivalent to:
    #return(.Call("R_fitted_node", node, data, vmatch, as.integer(obs)))

    if (is.terminal(node))
        return(rep(id_node(node), length(obs)))
    retid <- nextid <- kidids_node(node, data, vmatch, obs)
    for (i in unique(nextid)) {
        indx <- nextid == i
        retid[indx] <- fitted_node(kids_node(node)[[i]], data,
                                   vmatch, obs[indx])
    }
    return(retid)
}


length.partynode <- function(x)
    length(kids_node(x))

"[.partynode" <- "[[.partynode" <- function(x, i, ...) {
    stopifnot(length(i) == 1 & is.numeric(i))
    kids_node(x)[[i]]
}

split_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$split
}

surrogates_node <- function(node) {
    stopifnot(inherits(node, "partynode"))
    node$surrogates
}

is.terminal <- function(x, ...)
    UseMethod("is.terminal")

is.terminal.partynode <- function(x, ...) {
    kids <- is.null(kids_node(x))
    split <- is.null(split_node(x))
    stopifnot(kids == split)
    kids
}

depth <- function(x, ...)
    UseMethod("depth")

depth.partynode <- function(x) {
    if (is.terminal(x)) return(1)
    max(sapply(kids_node(x), depth)) + 1
}

width <- function(x, ...)
    UseMethod("width")

width.partynode <- function(x) {
    if (is.terminal(x)) return(1)
    sum(sapply(kids_node(x), width.partynode))
}
