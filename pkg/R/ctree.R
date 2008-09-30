
### calculate p-value
pmaxT <- function(lin, exp, cov) {
    v <- diag(V <- matrix(cov, ncol = sqrt(length(cov))))
    lin <- as.vector(lin)[v > 0]
    exp <- as.vector(exp)[v > 0]
    V <- V[v > 0, v > 0, drop = FALSE]
    v <- v[v > 0]
    if (length(v) == 0) return(c(0, 0))
    maxT <- as.vector(max(abs(lin - exp) / sqrt(v)))
    if (is.na(maxT)) return(c(0, 0))
    c(maxT, pmvnorm(lower = rep(-maxT, length(v)), 
                    upper = rep(maxT, length(v)), 
                    sigma = cov2cor(V)))
}

### set up new node for conditional inference tree
cnode <- function(id = 1, data, response, inputs, weights, ctrl) {

    weights <- as.integer(weights)
    if (sum(weights) < ctrl$minsplit) return(partynode(as.integer(id)))

    p <- sapply(inputs, function(i) {
        x <- data[[i]]
        lin <- .Call("R_LinstatExpCov", x, response, weights)
        do.call("pmaxT", lin[-1])
    })
    colnames(p) <- colnames(data)[inputs]
    rownames(p) <- c("maxT", "pval")
    p["pval",] <- p["pval",]^length(inputs)
    if (any(p["pval",] > 0.95)) {
        isel <- which.max(p["pval",])
        if (p["pval",isel] == 1) isel <- which.max(p["maxT",])
        isel <- inputs[isel]
    } else {
        return(partynode(as.integer(id), info = p))
    }

    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    storage.mode(mb) <- "integer"

    swp <- ceiling(sum(weights[!is.na(data[[isel]])]) * mp)
    if (mb < swp) mb <- as.integer(swp)

    sp <- .Call("R_split", data[[isel]], response, weights, mb)
    if (any(is.na(sp))) return(partynode(as.integer(id), info = p))
    if (length(sp) == 1) {
        thissplit <- partysplit(as.integer(isel), breaks = sp)
    } else {
        thissplit <- partysplit(as.integer(isel), index = sp)
    }
    kidids <- kidids_split(thissplit, data)

    left <- weights
    left[kidids == 2] <- 0
    leftnode <- cnode(id + 1, data, response, inputs, left, ctrl)

    right <- weights
    right[kidids == 1] <- 0
    rightnode <- cnode(max(nodeids(leftnode)) + 1, data, response, 
                       inputs, right, ctrl)

    return(partynode(as.integer(id), split = thissplit, 
                kids = list(leftnode, rightnode), info = p))
}

ctree <- function(formula, data, weights, subset, na.action, 
                  ctrl = list(minsplit = 20L, minbucket = 7L, minprob = 0.01)) {

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), 
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- FALSE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
 
    response <- names(mf)[1] ### model.response(mf)
    weights <- model.weights(mf)
    dat <- mf[, colnames(mf) != "(weights)"]
    ret <- ctree_fit(dat, response, weights = weights, ctrl = ctrl)
    ret$terms <- terms(formula, data = mf)
    return(ret)
}

ctree_fit <- function(data, response, weights = NULL, 
                      ctrl = list(minsplit = 20L, minbucket = 7L, minprob = 0.01)) {

    inputs <- which(!(colnames(data) %in% response))

    infl <- y2infl(data, response)

    if (is.null(weights))
        weights <- rep(1, nrow(data))
    storage.mode(weights) <- "integer"

    tree <-  cnode(1L, data, infl, inputs, weights, ctrl)
    fitted <- data.frame("(fitted)" = fitted_node(tree, data), 
                         "(response)" = data[ , response, drop = TRUE], 
                         check.names = FALSE)
    ret <- party(tree, data = data, fitted = fitted)
    class(ret) <- c("constparty", class(ret))
    return(ret)
}

logrank_trafo <- function(x, ties.method = c("logrank", "HL")) {
    ties.method <- match.arg(ties.method)
    time <- x[,1]
    event <- x[,2]
    n <- length(time)
    ot <- order(time, event)
    rt <- rank(time, ties.method = "max")
    mt <- rank(time, ties.method = "min") - 1
    fact <- switch(ties.method, "logrank" = event / (n - mt),
                                "HL" = event/(n - rt + 1)
                  )   
    event - cumsum(fact[ot])[rt]
}

### convert response y to influence function h(y)
y2infl <- function(data, response) {

    if (length(response) == 1) {
        response <- data[[response]]
        rtype <- class(response)[1]
        if (rtype == "integer") rtype <- "numeric"

        infl <- switch(rtype,
            "factor" = model.matrix(~ response - 1),
            "ordered" = (1:nlevels(response))[as.integer(response)],
            "numeric" = response,
            "Surv" = logrank_trafo(response)
        )
    } else {
        ### multivariate response
        infl <- lapply(response, y2infl, data = data)
        infl <- do.call("cbind", infl)
    }
    storage.mode(infl) <- "double"
    return(infl)
}
