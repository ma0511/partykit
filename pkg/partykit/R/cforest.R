
cforest <- function(formula, data, weights, subset, na.action = na.pass, 
                    control = ctree_control(...), ytrafo = NULL, 
                    ntree = 500L, replace = TRUE, fraction = 0.632,
                    scores = NULL, ...) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    
    ### only necessary for extended model formulae 
    ### e.g. multivariate responses
    if (require("Formula")) {
        formula <- Formula(formula)
    } else {
        if (length(formula[[2]]) > 1)
            stop("Package ", sQuote("Formula"),
                 " not available for handling extended model formula ",
                 sQuote("formula"))
    }
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    response <- names(mf)[1]
    if (inherits(formula, "Formula"))
        response <- names(model.part(formula, mf, lhs = 1))
    weights <- model.weights(mf)
    dat <- mf[, colnames(mf) != "(weights)"]
    if (!is.null(scores)) {
        for (n in names(scores)) {
            sc <- scores[[n]]
            if (is.ordered(dat[[n]]) && 
                nlevels(dat[[n]]) == length(sc)) {
                attr(dat[[n]], "scores") <- as.numeric(sc)
            } else {
                warning("scores for variable ", sQuote(n), " ignored")
            }
        }
    }

    if (is.null(weights))
        weights <- rep(1, nrow(mf))
    storage.mode(weights) <- "integer"

    nvar <- sum(!(colnames(dat) %in% response))

    control$cfun <- function(...) {
        if (control$teststat == "quad")
            p <- partykit:::.pX2(..., pval = (control$testtype != "Teststatistic"))
        if (control$teststat == "max")
            p <- partykit:::.pmaxT(..., pval = (control$testtype != "Teststatistic"))
        names(p) <- c("teststat", "pval")

        if (control$testtype == "Bonferroni")
            p["pval"] <- p["pval"] * min(nvar, control$mtry)
        crit <-  p["teststat"]
        if (control$testtype != "Teststatistic")
        crit <- p["pval"]
        c(crit, p)
    }

    if (replace) {
        rw <- rmultinom(ntree, size = nrow(dat), prob = weights / sum(weights))
    } else {
        idx <- 1:nrow(dat)
        rw <- matrix(0L, nrow = nrow(dat), ncol = ntree)
        for (b in 1:ntree) {
            i <- sample(idx, floor(fraction * nrow(dat)), prob = weights / sum(weights))
            rw[i, b] <- 1L
        }
    }

    forest <- lapply(1:ntree, function(b) {
        tree <- partykit:::.ctree_fit(dat, response, weights = rw[,b], ctrl = control, 
                           ytrafo = ytrafo)
        fitted <- data.frame("(fitted)" = fitted_node(tree, dat),
                             "(weights)" = rw[,b],
                             check.names = FALSE)
        ret <- party(tree, data = dat[-(1:nrow(dat)),,drop = TRUE], fitted = fitted)   
        class(ret) <- c("constparty", class(ret))
        ret
    })

    ret <- list(
        nodes = forest,
        "(response)" = dat[,response, drop = length(response) == 1]
    )
    class(ret) <- "cforest"

    ret$terms <- terms(mf)
    return(ret)
}


predict.cforest <- function(object, newdata, OOB = FALSE, fun = weighted.mean) {

    w <- matrix(0L, nrow = NROW(object[["(response)"]]), ncol = NROW(newdata))
    for (b in 1:length(object$nodes)) {

        ids <- nodeids(object$nodes[[b]], terminal = TRUE)
        fitted <- object$nodes[[b]]$fitted
        f <- fitted[["(fitted)"]]
        rw <- fitted[["(weights)"]]
        if (OOB) rw <- as.integer(rw == 0)
        pw <- sapply(ids, function(i) rw * (f == i))

        nd <- predict(object$nodes[[b]], newdata = newdata, type = "node")
        w <- w + pw[, match(nd, ids)]
    }

    ret <- vector(mode = "list", length = ncol(w))
    for (j in 1:ncol(w))
        ret[[j]] <- fun(object[["(response)"]], w[,j])
    ret
}

if (FALSE) {
library("partykit")
cf <- cforest(dist ~ speed, data = cars)

p <- predict(cf, newdata = cars)

plot(dist ~ speed, data = cars)
lines(cars$speed, unlist(p))
}
