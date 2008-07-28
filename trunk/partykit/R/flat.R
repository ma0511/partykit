flat2rec <- function(obj) {

    if (!all(sapply(obj, is.flat)))
        stop(sQuote("obj"), " ", "is not a list of flat", " ", 
             sQuote("node"), " ", "objects")
    
    obj <- obj[order(sapply(obj, function(node) get_id(node)))]
    if (length(obj) == 1) return(obj)

    new_recnode <- function(id) {
        if (is.null(kids_node(obj[[id]])))
            node(id = id, info = obj[[id]]$info)
        else
            node(id = id, split = split_node(obj[[id]]),
                 kids = lapply(kids_node(obj[[id]]), new_recnode),
		 surrogates = surrogates_node(obj[[id]]),
                 info = obj[[id]]$info)
    }
        
    node <- node(id = as.integer(1), split = split_node(obj[[1]]),
                 kids = lapply(kids_node(obj[[1]]), new_recnode), 
		 surrogates = surrogates_node(obj[[1]]),
                 info = obj[[1]]$info)
    return(node)
}

rec2flat <- function(node) {

    if (is.terminal(node))
        return(list(node))

    if (is.flat(node))
        stop(sQuote("node"), " ", "is not a recursive node")

    obj <- list()
    
    new_flatnode <- function(node) {
        if (is.terminal(node))
            obj[[node$id]] <<- node(id = get_id(node), info = node$info)
        else {
            obj[[node$id]] <<- node(id = get_id(node), split = split_node(node),
                 kids = sapply(kids_node(node), function(k) get_id(k)), 
		 surrogates = surrogates_node(node),
		 info = node$info)
            lapply(kids_node(node), new_flatnode)
        }
    }
    new_flatnode(node)
    return(obj)
}

is.flat <- function(node)
    is.terminal(node) | is.integer(kids_node(node))
