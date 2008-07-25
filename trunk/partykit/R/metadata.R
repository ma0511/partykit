
new_metadata <- function(varnames, class, levels, responses = 1,
                         inputs = 2:length(varnames)) {

    stopifnot(is.character(varnames))
    stopifnot(is.character(class))
    stopifnot(length(varnames) == length(class))
    stopifnot(length(varnames) == length(levels))
    indx <- 1:length(varnames)

    if (is.character(responses))
        responses <- indx[varnames %in% responses]
    if (is.character(inputs))
        inputs <- indx[varnames %in% inputs]
    stopifnot((all(responses %in% indx) & all(inputs %in% indx)))
    stopifnot(!any(responses %in% inputs))

    rval <- list(varnames = varnames,
                 class = class,
                 levels = levels,
                 responses = as.integer(responses),
                 inputs = as.integer(inputs))
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
