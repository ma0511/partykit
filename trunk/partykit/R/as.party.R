
as.party <- function(obj, ...)
    UseMethod("as.party")

as.party.rpart <- function(obj, ...) {

    ff <- obj$frame
    n  <- nrow(ff)
    if (n==1) return(new_node(as.integer(1)))  # special case of no splits

    is.leaf <- (ff$var == "<leaf>")
    vnames <- ff$var[!is.leaf]  #the variable names for the primary splits

    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
    splitindex <- list()
    splitindex$primary <- numeric(n)
    splitindex$primary[!is.leaf] <- index[c(!is.leaf, FALSE)]
    splitindex$surrogate <- lapply(1:n, function(i) {
        prim <- splitindex$primary[i]
        if (prim < 1 || ff[i, "nsurrogate"] == 0) return(NULL)
        else return(prim + ff[i, "ncompete"] + 1:ff[i, "nsurrogate"])
    })
    
    rpart_metadata <- function() {
        ### extract metadata from terms
        metadata <- list(varnames = names(attr(obj$terms, "dataClasses")),
                         class = as.vector(attr(obj$terms, "dataClasses")))
        ylev <- list(attr(obj, "ylevels"))
        names(ylev) <- metadata$varnames[1]
        lev <- c(attr(obj, "xlevels"), ylev)
        metadata$levels <- lapply(metadata$varnames, function(var) lev[[var]])
        class(metadata) <- "metadata"
        metadata
    }
    objmeta <- rpart_metadata()

    rpart_kids <- function(i) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1, 
            which((cumsum(!is.leaf[-(1:i)]) + 1) == cumsum(is.leaf[-(1:i)]))[1] + 1 + i))
    }

    rpart_onesplit <- function(j) {
        if (j < 1) return(NULL)
        ### numeric
        if (abs(obj$split[j, "ncat"]) == 1) {
            ret <- new_split(fun = which(rownames(obj$split)[j] == objmeta$varnames),
                      breaks = as.double(obj$split[j, "index"]),
                      right = FALSE,
                      index = if(obj$split[j, "ncat"] > 0) 2:1)
        } else {
            index <- obj$csplit[obj$split[j, "index"],]
            index[index == 2] <- NA
            index[index == 3] <- 2
            ret <- new_split(fun = which(rownames(obj$split)[j] == objmeta$varnames),
                      index = as.integer(index))
        }
        ret
    }
                      
    rpart_splits <- function(i)
        lapply(c(splitindex$primary[i], splitindex$surrogate[[i]]), rpart_onesplit)
    
    rpart_node <- function(i) {
        if (is.null(rpart_kids(i))) return(new_node(as.integer(i)))
        new_node(as.integer(i), split = rpart_splits(i), kids = 
            lapply(rpart_kids(i), rpart_node))
    }

    node <- rpart_node(1)

    new_party(node = node, metadata = objmeta)
}
