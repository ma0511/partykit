
mst <- function(x) {
    minmin <- function(y, xx) {

        ### don't complain about Inf
        ww <- getOption("warn")
        on.exit(options(warn = ww))
        options(warn = -1)

        if (any(!is.na(y) & y==.5)) {
            if (any(!is.na(y) & y <.5))
                .5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
            else
                .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
        } else   min(xx[!is.na(y) & y<=.5])
    }
    med <- minmin(x$surv, x$time)
    return(med)
}

dostep <- function(x, y) {

    ### create a step function based on x, y coordinates
    ### modified from `survival:print.survfit'
    if (is.na(x[1] + y[1])) {
        x <- x[-1]
        y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {
        # replace verbose horizonal sequences like
        # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
        # with (1, .2), (3, .1).  They are slow, and can smear the looks
        # of the line type.
        dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
        n2 <- sum(dupy)

        #create a step function
        xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
        yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
        RET <- list(x = xrep, y = yrep)
    } else {
        if (n == 1) {
            RET <- list(x = x, y = y)
        } else {
            RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
        }
    }
    return(RET)
}
