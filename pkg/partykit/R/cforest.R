
perturbe <- function(replace = TRUE, fraction = .632) {
    ret <- function(prob) {
        if (replace) {
            rw <- rmultinom(1, size = length(prob), prob = prob)
        } else {
            rw <- integer(length(prob))
            i <- sample(1:length(prob), floor(fraction * length(prob)), prob = prob)
            rw[i] <- 1L
        }
        rw
    }
    ret
}

boot <- function() perturbe(replace = TRUE)
sampsplit <- function(fraction = 0.632) perturbe(replace = FALSE, fraction = fraction)

cforest <- function(formula, data, weights, subset, na.action = na.pass, 
                    control = ctree_control(...), ytrafo = NULL, 
                    ntree = 500L, perturbe = sampsplit(fraction = 0.632),
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

    probw <- weights / sum(weights)
    emptydat <- dat[-(1:nrow(dat)),,drop = TRUE]
    
    forest <- lapply(1:ntree, function(b) {

        rw <- perturbe(probw)

        tree <- partykit:::.ctree_fit(dat, response, weights = rw, ctrl = control, 
                                      ytrafo = ytrafo)
        fitted <- data.frame("(fitted)" = fitted_node(tree, dat),
                             "(weights)" = rw,
                             check.names = FALSE)
        ret <- party(tree, data = emptydat, fitted = fitted)   
        class(ret) <- c("constparty", class(ret))
        ret
    })

    ret <- list(
        forest = forest,
        "(response)" = dat[,response, drop = length(response) == 1]
    )
    class(ret) <- "cforest"

    ret$terms <- terms(mf)
    return(ret)
}


predict.cforest <- function(object, newdata = NULL, type = c("weights", "response", "prob"), 
                            OOB = FALSE, FUN = NULL, simplify = TRUE) {

    responses <- object[["(response)"]]
    if (is.null(newdata)) {
        nam <- rownames(object$forest[[1]]$fitted)
    } else {
        nam <- rownames(newdata)
    }

    w <- matrix(0L, nrow = NROW(object[["(response)"]]), ncol = NROW(newdata))
    forest <- object$forest
    for (b in 1:length(forest)) {

        ids <- nodeids(forest[[b]], terminal = TRUE)
        fitted <- forest[[b]]$fitted
        f <- fitted[["(fitted)"]]
        rw <- fitted[["(weights)"]]
        if (OOB) rw <- as.integer(rw == 0)
        pw <- sapply(ids, function(i) rw * (f == i))

        if (is.null(newdata)) {
            nd <- f
        } else {
            nd <- predict(forest[[b]], newdata = newdata, type = "node")
        }
        w <- w + pw[, match(nd, ids)]
    }

    type <- match.arg(type)
    if (type == "weights") {
        ret <- w
        rownames(ret) <- rownames(newdata)
        return(ret)
    }
    
    pfun <- function(response) {

        if (is.null(FUN)) {

            rtype <- class(response)[1]
            if (rtype == "ordered") rtype <- "factor"
            if (rtype == "integer") rtype <- "numeric"

            FUN <- switch(rtype,
                "Surv" = if (type == "response") partykit:::.pred_Surv_response else partykit:::.pred_Surv,
                "factor" = if (type == "response") partykit:::.pred_factor_response else partykit:::.pred_factor,
                "numeric" = if (type == "response") partykit:::.pred_numeric_response else partykit:::.pred_ecdf)
        }

        ret <- vector(mode = "list", length = ncol(w))
        for (j in 1:ncol(w))
            ret[[j]] <- FUN(response, w[,j])
        ret <- as.array(ret)
        dim(ret) <- NULL
        names(ret) <- rownames(newdata)
         
        if (simplify)
            ret <- partykit:::.simplify_pred(ret, names(ret), names(ret))
        ret
    }
    if (NCOL(responses) == 1) {
        ret <- pfun(responses)
    } else {
        ret <- lapply(responses, pfun)
        if (all(sapply(ret, is.atomic)))
            ret <- as.data.frame(ret)
        names(ret) <- colnames(response)
    }
    ret
}

.pred_quantile <- function(probs = c(0.1, 0.5, 0.9), ...)
    function(y, w) quantile(rep(y, w), probs = probs, ...)

.pred_density <- function(y, w)
    density(y, weights = w / sum(w))

if (FALSE) {

library("partykit")
cf <- cforest(dist ~ speed, data = cars, ytrafo = list(dist = function(x) cbind(x, x^2)))

p <- predict(cf, newdata = cars, type = "response")

plot(dist ~ speed, data = cars)
lines(cars$speed, unlist(p))

p <- predict(cf, newdata = cars, type = "response", FUN = .pred_quantile())

plot(dist ~ speed, data = cars)
lines(cars$speed, p[,1])
lines(cars$speed, p[,2])
lines(cars$speed, p[,3])

p <- predict(cf, newdata = cars, type = "response", FUN = .pred_density)

for (i in 1:length(p)) 
    plot(p[[i]], xlim = c(-10, 150), ylim = c(0, .025))

}
