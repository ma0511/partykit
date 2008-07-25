
new_metadata <- function(varnames, class, levels) {

    rval <- list(varnames = varnames,
                 class = class,
                 levels = levels)
    class(rval) <- "metadata"
    rval
}

metadata <- function(x, ...)
    UseMethod("metadata")

metadata.data.frame <- function(x, ...) {

    new_metadata(varnames = colnames(x),
                 class = sapply(x, function(z) class(z)[1]),
                 levels = lapply(x, levels))
}
