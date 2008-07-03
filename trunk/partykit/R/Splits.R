
## Basic idea:
## A split is a function that maps data (specifying the partitioning variables)
## to a set of integers (specifying the daughter nodes).
## 
## Extensions:
##   - The function does not need to go all the way itself,
##     we can help in doing the splitting of some numeric score of the data.
##   - The function might be very simple: the identity of a single variable
##     from the data.
##
## Arguments:
##  - fun: Either an R function computing some splitting score from the data,
##    or an integer specifying the column in the data.
##    The resulting splitting score must be
##  - breaks: A numeric vector of breaks for splitting the score via cut().
##    Needs to be double for numeric scores and integer for factor scores.
##    

newsplit <- function(fun, breaks = NULL, ids = NULL, right = TRUE) {

    ### informal class for splits
    split <- vector(mode = "list", length = 4)
    names(split) <- c("fun", "breaks", "ids", "right")

    ### split is either a function or an id referring to a variable
    if (is.function(fun) || is.integer(fun)) {
        split$fun <- fun
    } else {
        stop(sQuote("fun"), " ", "should be a function or an integer")
    }

    ### vec
    if (!is.null(breaks)) {
        if (is.double(breaks) && (length(breaks) >= 1)) {
            split$breaks <- breaks
        } else {
            stop(sQuote("break"), " ",
                 "should be a double with at least one element")
        }
    }

    if (!is.null(ids)) {
        if (is.integer(ids)) {
            if (!(length(ids) >= 2)) 
                stop(sQuote("ids"), " ", "has less than two elements")
            if (!(min(ids) == 1))
                stop("minimum of", " ", sQuote("ids"), " ", "is not equal to 1")
            if (!all.equal(diff(sort(unique(ids))), rep(1, max(ids)-1)))
                stop(sQuote("ids"), " ", "is not a contiguous sequence")
            split$ids <- ids
        } else {
            stop(sQuote("ids"), " ", "is not a class", " ", sQuote("integer"))
        }
        if (!is.null(breaks)) {
            if (length(breaks) != (length(ids) - 1))
                stop("length of", " ", sQuote("breaks"), " ", 
                     "does not match length of", " ", sQuote("ids"))
        }
    }

    if (!isTRUE(right))
        split$right <- FALSE

    return(split)
}

split <- function(data, split) {

    if (is.function(split$fun))
        stop("")

    x <- data[[split$fun]]

    if (is.null(split$breaks)) {
        stopifnot(storage.mode(x) == "integer")
    } else {
        stopifnot(storage.mode(x) == "double")
        x <- as.integer(cut(x, breaks = c(-Inf, split$breaks, Inf), 
                            right = is.null(split$right)))
    }
    if (!is.null(split$ids))
        x <- split$ids[x]
    return(x)
}

metadata <- function(data) {

    list(varnames = colnames(data),
         class = sapply(data, function(x) class(x)[1]),
         levels = lapply(data, levels))
}

nodelabels <- function(split, meta) {

    ## determine type
    type <- ifelse(is.function(split$fun), 
                   "function", meta$class[split$fun])
    type[!(type %in% c("function", "factor", "ordered"))] <- "numeric"

    ## process defaults for breaks and ids
    breaks <- split$breaks
    ids <- split$ids
    if (is.null(breaks) && is.null(ids)) {
        if (type %in% c("factor", "ordered")) {
            ids <- seq_along(meta$levels[[split$fun]])
        } else {
            stop("")
        }
    }
    if (is.null(breaks)) breaks <- 1:(length(ids) - 1)
    if (is.null(ids)) ids <- 1:(length(breaks) + 1)

    right <- is.null(split$right)

    ## check whether ordered are really ordered
    if (type == "ordered") {
        if (length(breaks) > 1)
            type <- "factor"
    }
    ### <FIXME> format ordered multiway splits? </FIXME>

    switch(type, 
        "factor" = {
            lev <- meta$levels[[split$fun]]
            nids <- as.integer(as.character(cut(seq_along(lev), c(-Inf, breaks, Inf),
                labels = ids, right = right)))
            dlab <- as.vector(tapply(lev, nids, paste, collapse = ", "))
            mlab <- meta$varnames[split$fun]
        },
        "ordered" = {
            lev <- meta$levels[[split$fun]]
            if (length(breaks) == 1) {
                if (right)
                    dlab <- paste(c("<=", ">"), lev[breaks], sep = " ")
                else
                    dlab <- paste(c("<", ">="), lev[breaks], sep = " ")
            } else {
                stop("") ### see above
            }
            dlab <- dlab[ids]
            mlab <- meta$varnames[split$fun]
        },
        "numeric" = {
            if (length(breaks) == 1) {
                if (right)
                    dlab <- paste(c("<=", ">"), breaks, sep = " ")
                else
                    dlab <- paste(c("<", ">="), breaks, sep = " ")
            } else {
                dlab <- levels(cut(0, breaks = c(-Inf, breaks, Inf), right = right))
            }
            dlab <- as.vector(tapply(dlab, ids, paste, collapse = " | "))
            mlab <- meta$varnames[split$fun]
        }
    )

    return(list(dlab = dlab, mlab = mlab))
} 


dat <- data.frame(x = gl(3, 30, labels = LETTERS[1:3]), y = rnorm(90), 
                  z = gl(9, 10, labels = LETTERS[1:9], ordered = TRUE))
csp <- newsplit(as.integer(1), ids = as.integer(c(1, 2, 1)))
split(dat, csp)

nsp <- newsplit(as.integer(2), breaks = c(-1, 0, 1), ids = as.integer(c(1, 2, 1, 3)))
split(dat, nsp)

osp <- newsplit(as.integer(3), breaks = c(3, 6), ids = as.integer(c(2, 1, 2)))

nodelabels(csp, metadata(dat))
nodelabels(nsp, metadata(dat))
nodelabels(osp, metadata(dat))

