
median_survival_time <- function(x) {
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
