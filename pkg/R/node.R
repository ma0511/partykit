
partynode <- function(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL) {

    if (!is.integer(id) || length(id) != 1)
        stop(sQuote("id"), " ", "must be a single integer")

    if (is.null(split) != is.null(kids)) {
        stop(sQuote("split"), " ", "and", " ", sQuote("kids"), " ", 
             "must either both be specified or unspecified")
    }

    if (!is.null(kids)) {
        if (!(is.list(kids) && all(sapply(kids, inherits, "partynode"))) 
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

as.partynode.list <- function(obj, ...) {

    if (!all(sapply(obj, inherits, what = "list")))
        stop("'obj' has to be a list of lists")

    if (!all(sapply(obj, function(x) "id" %in% names(x))))
        stop("each list in 'obj' has to define a node 'id'")

    ok <- sapply(obj, function(x) 
              all(names(x) %in% c("id", "split", "kids", "surrogates", "info")))
    if (any(!ok))
        sapply(which(!ok), function(i) 
            warning(paste("list element", i, "defines additional elements:", 
                          paste(names(x[[i]])[!(names(x[[i]]) %in% 
                                c("id", "split", "kids", "surrogates", "info"))], 
                                collapse = ", "))))
    
    ids <- sapply(obj, function(node) node$id)
    if (!all(ids %in% 1:length(obj)))
        stop("ids must match 1:length(obj)")

    obj <- obj[order(ids)]
    if (length(obj) == 1) return(do.call("partynode", obj[[1]]))

    new_recnode <- function(id) {
        if (is.null(obj[[id]]$kids))
            partynode(id = id, info = obj[[id]]$info)
        else
            partynode(id = id, split = obj[[id]]$split,
                 kids = lapply(obj[[id]]$kids, new_recnode),
		 surrogates = obj[[id]]$surrogates,
                 info = obj[[id]]$info)
    }
        
    node <- partynode(id = as.integer(1), split = obj[[1]]$split,
                 kids = lapply(obj[[1]]$kids, new_recnode), 
		 surrogates = obj[[1]]$surrogates,
                 info = obj[[1]]$info)
    return(node)
}

as.list.partynode <- function(x, ...) {

    obj <- list()
    
    nodelist <- function(node) {
        if (is.terminal(node))
            obj[[node$id]] <<- list(id = id_node(node), info = info_node(node))
        else {
            obj[[node$id]] <<- list(id = id_node(node), split = split_node(node),
                 kids = sapply(kids_node(node), function(k) id_node(k)))
             if (!is.null(surrogates_node(node)))
		 obj[[node$id]]$surrogates <- surrogates_node(node)
             if (!is.null(info_node(node)))
		 obj[[node$id]]$info <- info_node(node)
            lapply(kids_node(node), nodelist)
        }
    }
    nodelist(x)
    return(obj)
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
