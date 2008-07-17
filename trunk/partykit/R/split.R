
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

new_split <- function(fun, breaks = NULL, index = NULL, right = TRUE, 
                      prob = NULL, ...) {

    ### informal class for splits
    split <- vector(mode = "list", length = 6)
    names(split) <- c("fun", "breaks", "index", "right", "prob", "info")

    ### split is either a function or an id referring to a variable
    if (is.function(fun) || is.integer(fun)) {
        split$fun <- fun
    } else {
        stop(sQuote("fun"), " ", "should be a function or an integer")
    }

    if (is.null(breaks) && is.null(index))
        stop("either", " ", sQuote("breaks"), " ", "or", " ",
             sQuote("index"), " ", "must be given")

    ### vec
    if (!is.null(breaks)) {
        if (is.numeric(breaks) && (length(breaks) >= 1)) {
            split$breaks <- breaks
        } else {
            stop(sQuote("break"), " ",
                 "should be a numeric vector containing at least one element")
        }
    }

    ### FIXME: index for factors may contain NA
    if (!is.null(index)) {
        if (is.integer(index)) {
            if (!(length(index) >= 2)) 
                stop(sQuote("index"), " ", "has less than two elements")
            if (!(min(index, na.rm = TRUE) == 1))
                stop("minimum of", " ", sQuote("index"), " ", "is not equal to 1")
            if (!all.equal(diff(sort(unique(index))), rep(1, max(index, na.rm = TRUE) - 1)))
                stop(sQuote("index"), " ", "is not a contiguous sequence")
            split$index <- index
        } else {
            stop(sQuote("index"), " ", "is not a class", " ", sQuote("integer"))
        }
        if (!is.null(breaks)) {
            if (length(breaks) != (length(index) - 1))
                stop("length of", " ", sQuote("breaks"), " ", 
                     "does not match length of", " ", sQuote("index"))
        }
    }

    if (is.logical(right) & !is.na(right))
        split$right <- right
    else
        stop(sQuote("right"), " ", "is not a logical")

    if (!is.null(prob)) {
        if (!is.double(prob) || 
            (any(prob < 0) | any(prob > 1) | !isTRUE(all.equal(sum(prob), 1))))
            stop(sQuote("prob"), " ", "is not a vector of probabilities")
        if (!is.null(index))
            stopifnot(max(index, na.rm = TRUE) == length(prob))
        if (!is.null(breaks) && is.null(index))
            stopifnot(length(breaks) == (length(prob) - 1))
        split$prob <- prob
    }

    split$info <- list(...)

    class(split) <- "split"

    return(split)
}

get_fun <- function(split) {
    stopifnot(inherits(split, "split"))
    split$fun
}

get_varid <- function(split) {
    if (is.functional(split))
        stop(sQuote("split"), " ", "is a functional split")
    get_fun(split)
}

get_breaks <- function(split) {
    stopifnot(inherits(split, "split"))
    split$breaks
}

get_index <- function(split) {
    stopifnot(inherits(split, "split"))
    split$index
}

get_right <- function(split) {
    stopifnot(inherits(split, "split"))
    split$right
}

is.functional <- function(split) {
    stopifnot(inherits(split, "split"))
    is.function(get_fun(split))
}

get_primary <- function(splitlist) {
    stopifnot(all(sapply(splitlist, inherits, "split")))
    splitlist[[1]]
}

get_surrogates <- function(splitlist) {
    stopifnot(all(sapply(splitlist, inherits, "split")))
    splitlist[-1]
}

do_split <- function(split, data) {

    id <- get_varid(split)
    x <- data[[id]]

    if (is.null(get_breaks(split))) {
        stopifnot(storage.mode(x) == "integer")
    } else {
        ### FIXME: rpart may have double splits for integer variables
        # stopifnot(storage.mode(x) == storage.mode(split$breaks))
        x <- as.integer(cut(as.numeric(x), breaks = c(-Inf, get_breaks(split), Inf), 
                            right = get_right(split)))
    }
    index <- get_index(split)
    if (!is.null(index))
        x <- index[x]
    return(x)
}

do_splitlist <- function(splitlist, data) {

    p <- ncol(data)
    ### replace functional splits by splitting score id
    for (i in 1:length(splitlist)) {
        if (is.functional(splitlist[[i]])) {
            p <- p + 1
            data[p] <- splitlist[[i]]$fun(data)
            splitlist[[i]]$fun <- as.integer(p)
        }
    }
    primary <- get_primary(splitlist)
    surrogates <- get_surrogates(splitlist)

    ### perform primary split
    x <- do_split(primary, data)

    ### surrogate / random splits if needed
    if (any(is.na(x))) {
        ### surrogate splits
        if (length(surrogates) >= 1) {
            for (surr in surrogates) {
                nax <- is.na(x)
                if (!any(nax)) break;
                x[nax] <- do_split(surr, data[nax,, drop = FALSE])
            }
        }
        nax <- is.na(x)
        ### random splits
        if (any(nax)) {
            index <- get_index(primary)
            if (is.null(index)) {
                nd <- length(levels(data[[get_varid(primary)]]))
            } else {
                nd <- max(index, na.rm = TRUE)
            }
            prob <- get_prob(primary)
            if (is.null(prob))
                prob <- rep(1, nd) / nd
            x[nax] <- sample(1:nd, sum(nax), prob = prob)
        }
    }
    return(x)
}

metadata <- function(x, ...)
    UseMethod("metadata")

metadata.data.frame <- function(x, ...) {

    rval <- list(varnames = colnames(x),
         class = sapply(x, function(z) class(z)[1]),
         levels = lapply(x, levels))
    class(rval) <- "metadata"
    rval
}

nodelabels <- function(split, meta, digits = getOption("digits") - 2) {

    ## determine type
    type <- ifelse(is.function(split$fun), 
                   "function", meta$class[split$fun])
    type[!(type %in% c("function", "factor", "ordered"))] <- "numeric"

    ## process defaults for breaks and index
    breaks <- get_breaks(split)
    index <- get_index(split)
    right <- get_right(split)
    if (is.null(breaks)) breaks <- 1:(length(index) - 1)
    if (is.null(index)) index <- 1:(length(breaks) + 1)

    ## check whether ordered are really ordered
    if (type == "ordered") {
        if (length(breaks) > 1)
            type <- "factor"
    }
    ### <FIXME> format ordered multiway splits? </FIXME>

    switch(type, 
        "factor" = {
            lev <- meta$levels[[get_varid(split)]]
            nindex <- as.integer(as.character(cut(seq_along(lev), c(-Inf, breaks, Inf),
                labels = index, right = right)))
            dlab <- as.vector(tapply(lev, nindex, paste, collapse = ", "))
            mlab <- meta$varnames[split$fun]
        },
        "ordered" = {
            lev <- meta$levels[[get_varid(split)]]
            if (length(breaks) == 1) {
                if (right)
                    dlab <- paste(c("<=", ">"), lev[breaks], sep = " ")
                else
                    dlab <- paste(c("<", ">="), lev[breaks], sep = " ")
            } else {
                stop("") ### see above
            }
            dlab <- dlab[index]
            mlab <- meta$varnames[get_varid(split)]
        },
        "numeric" = {
            breaks <- round(breaks, digits)
            if (length(breaks) == 1) {
                if (right)
                    dlab <- paste(c("<=", ">"), breaks, sep = " ")
                else
                    dlab <- paste(c("<", ">="), breaks, sep = " ")
            } else {
                dlab <- levels(cut(0, breaks = c(-Inf, breaks, Inf), 
                                   right = right))
            }
            dlab <- as.vector(tapply(dlab, index, paste, collapse = " | "))
            mlab <- meta$varnames[get_varid(split)]
        },
        "function" = {
            ### FIXME: use possible names attribute for split$fun
            mlab <- "functional split"
            if (!is.integer(breaks)) {
                breaks <- round(breaks, digits)
                if (length(breaks) == 1) {
                     if (right)  
                         dlab <- paste(c("<=", ">"), breaks, sep = " ")
                     else
                         dlab <- paste(c("<", ">="), breaks, sep = " ")
                } else {
                    dlab <- levels(cut(0, breaks = c(-Inf, breaks, Inf), right = right))
                }
                dlab <- as.vector(tapply(dlab, index, paste, collapse = " | "))
            } else if (length(breaks) == 1 && length(index) > 2) {
                if (right)  
                    dlab <- paste(c("<=", ">"), breaks, sep = " ")
                else
                    dlab <- paste(c("<", ">="), breaks, sep = " ")
            } else {
                lev <- 1:max(index, na.rm = TRUE)
                nindex <- as.integer(as.character(cut(seq_along(lev), c(-Inf, breaks, Inf),
                    labels = index, right = right)))
                dlab <- as.vector(tapply(lev, nindex, paste, collapse = ", "))
            }
        }
    )

    return(list(splitname = mlab, splitlevels = dlab))
} 
