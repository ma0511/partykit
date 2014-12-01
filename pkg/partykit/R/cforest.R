
### constructor for forest objects
constparties <- function(nodes, data, weights, fitted = NULL, terms = NULL, info = NULL) {

    stopifnot(all(sapply(nodes, function(x) inherits(x, "partynode"))))
    stopifnot(inherits(data, "data.frame"))
    stopifnot(inherits(weights, "list"))

    if(!is.null(fitted)) {
        stopifnot(inherits(fitted, "data.frame"))
        stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
        if (nrow(data) == 0L)
            stopifnot("(response)" %in% names(fitted))
    } else {
        stopifnot(nrow(data) > 0L)
        stopifnot(!is.null(terms))
        fitted <- data.frame("(response)" = model.response(model.frame(terms, data = data)),
                             check.names = FALSE)
    }

    ret <- list(nodes = nodes, data = data, weights = weights, fitted = fitted)
    class(ret) <- c("constparties", "parties")

    if(!is.null(terms)) {
        stopifnot(inherits(terms, "terms"))
        ret$terms <- terms
    }

    if (!is.null(info))
        ret$info <- info

    ret
}

.perturbe <- function(replace = TRUE, fraction = .632) {
    ret <- function(prob) {
        if (replace) {
            rw <- rmultinom(1, size = length(prob), prob = prob)
        } else {
            rw <- integer(length(prob))
            i <- sample(1:length(prob), floor(fraction * length(prob)), prob = prob)
            rw[i] <- 1L
        }
        as.integer(rw)
    }
    ret
}

boot <- function() .perturbe(replace = TRUE)
sampsplit <- function(fraction = 0.632) 
    .perturbe(replace = FALSE, fraction = fraction)

cforest <- function(formula, data, weights, subset, na.action = na.pass, 
                    control = ctree_control(...), ytrafo = NULL, scores = NULL,
                    ntree = 500L, perturbe = sampsplit(fraction = 0.632),
                    mtry = ceiling(sqrt(nvar)), myapply = lapply, ...) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    
    ### only necessary for extended model formulae 
    ### e.g. multivariate responses
    formula <- Formula::Formula(formula)
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    response <- names(Formula::model.part(formula, mf, lhs = 1))
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

    nvar <- sum(!(colnames(dat) %in% response)) 
    control$mtry <- mtry

    control$cfun <- function(...) {
        if (control$teststat == "quad")
            p <- .pX2(..., pval = (control$testtype != "Teststatistic"))
        if (control$teststat == "max") 
            p <- .pmaxT(..., pval = (control$testtype != "Teststatistic"))
        names(p) <- c("statistic", "p.value")

        if (control$testtype == "Bonferroni")
            p["p.value"] <- p["p.value"] * min(nvar, control$mtry)
        crit <-  p["statistic"]
        if (control$testtype != "Teststatistic")
        crit <- p["p.value"]
        c(crit, p)
    }

    if (!is.matrix(weights)) {
        probw <- weights / sum(weights)
        rw <- replicate(ntree, perturbe(probw), simplify = FALSE)
    } else {
        rw <- as.data.frame(weights)
    }

    
    forest <- myapply(1:ntree, function(b) {
        .ctree_fit(dat, response, weights = rw[[b]], 
                   ctrl = control, ytrafo = ytrafo)
    })

    ret <- constparties(nodes = forest, data = dat, weights = rw,
        fitted = data.frame("(response)" = dat[,response, drop = length(response) == 1], 
                            check.names = FALSE),
        terms = terms(mf))
    class(ret) <- c("cforest", class(ret))

    return(ret)
}


predict.cforest <- function(object, newdata = NULL, type = c("weights", "response", "prob"), 
                            OOB = FALSE, FUN = NULL, simplify = TRUE, ...) {

    responses <- object$fitted[["(response)"]]
    nd <- object$data
    if (!is.null(newdata)) {
        nd <- model.frame(object$terms, newdata)
        OOB <- FALSE
    }
    nam <- rownames(nd)

    ### regenerate weights
    rw <- object$weights

    w <- matrix(0L, nrow = NROW(responses), ncol = length(nam))
    forest <- object$nodes

    for (b in 1:length(forest)) {
        ids <- nodeids(forest[[b]], terminal = TRUE)
        f <- fitted_node(forest[[b]], nd)
        tw <- rw[[b]]
        if (OOB) tw <- as.integer(tw == 0)
        pw <- sapply(ids, function(i) tw * (f == i))
        w <- w + pw[, match(f, ids)]
    }

    type <- match.arg(type)
    if (type == "weights") {
        ret <- w
        colnames(ret) <- nam
        rownames(ret) <- rownames(responses)
        return(ret)
    }
    
    pfun <- function(response) {

        if (is.null(FUN)) {

            rtype <- class(response)[1]
            if (rtype == "ordered") rtype <- "factor"
            if (rtype == "integer") rtype <- "numeric"

            FUN <- switch(rtype,
                "Surv" = if (type == "response") .pred_Surv_response else .pred_Surv,
                "factor" = if (type == "response") .pred_factor_response else .pred_factor,
                "numeric" = if (type == "response") .pred_numeric_response else .pred_ecdf)
        }

        ret <- vector(mode = "list", length = ncol(w))
        for (j in 1:ncol(w))
            ret[[j]] <- FUN(response, w[,j])
        ret <- as.array(ret)
        dim(ret) <- NULL
        names(ret) <- nam
         
        if (simplify)
            ret <- .simplify_pred(ret, names(ret), names(ret))
        ret
    }
    if (!is.data.frame(responses)) {
        ret <- pfun(responses)
    } else {
        ret <- lapply(responses, pfun)
        if (all(sapply(ret, is.atomic)))
            ret <- as.data.frame(ret)
        names(ret) <- colnames(responses)
    }
    ret
}
