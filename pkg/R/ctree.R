
### calculate p-value
.pmaxT <- function(lin, exp, cov) {
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

.MPinv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    Xplus <- if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2:1])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
    list(Xplus = Xplus, rank = sum(Positive))
}

.pX2 <- function(lin, exp, cov) {
    tmp <- matrix(lin - exp, ncol = 1)
    Xplus <- .MPinv(matrix(cov, ncol = sqrt(length(cov))))
    X2 <- crossprod(tmp, Xplus$Xplus) %*% tmp
    c(X2, pchisq(X2, df = Xplus$rank, lower.tail = TRUE))
}

### set up new node for conditional inference tree
.cnode <- function(id = 1, data, response, inputs, weights, ctrl) {

    weights <- as.integer(weights)
    if (sum(weights) < ctrl$minsplit) return(partynode(as.integer(id)))
    if (id > 1 && ctrl$stump) return(partynode(as.integer(id)))

    p <- sapply(inputs, function(i) {
        x <- data[[i]]
        lin <- .Call("R_LinstatExpCov", x, response, weights)
        do.call(ctrl$teststatfun, lin[-1])
    })
    colnames(p) <- colnames(data)[inputs]
    rownames(p) <- c("teststat", "pval")
    
    if (ctrl$testtype == "Bonferroni")
        p["pval",] <- p["pval",]^length(inputs)
    crit <-  p["teststat",]
    if (ctrl$testtype != "Teststatistic" && max(p["pval",]) < 1)
        crit <- p["pval",]

    if (any(crit > ctrl$mincriterion)) {
        isel <- which.max(crit)
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
    leftnode <- .cnode(id + 1, data, response, inputs, left, ctrl)

    right <- weights
    right[kidids == 1] <- 0
    rightnode <- .cnode(max(nodeids(leftnode)) + 1, data, response, 
                       inputs, right, ctrl)

    return(partynode(as.integer(id), split = thissplit, 
                kids = list(leftnode, rightnode), info = p))
}

ctree_control <- function(teststat = c("quad", "max"), 
    testtype = c("Bonferroni", "Univariate", "Teststatistic"), 
    mincriterion = 0.95, minsplit = 20L, minbucket = 7L, minprob = 0.01,
    stump = FALSE) {

    teststat <- match.arg(teststat)
    teststatfun <- ifelse(teststat == "quad", ".pX2", ".pmaxT")
    testtype <- match.arg(testtype)
    list(teststat = teststat, 
         teststatfun = teststatfun, 
         testtype = testtype, mincriterion = mincriterion, 
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump)
}

ctree <- function(formula, data, weights, subset, na.action, 
                  control = ctree_control()) {

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
    ret <- .ctree_fit(dat, response, weights = weights, ctrl = control)
    ret$terms <- terms(formula, data = mf)
    return(ret)
}

.ctree_fit <- function(data, response, weights = NULL, 
                      ctrl = ctree_control()) {

    inputs <- which(!(colnames(data) %in% response))

    infl <- .y2infl(data, response)

    if (is.null(weights))
        weights <- rep(1, nrow(data))
    storage.mode(weights) <- "integer"

    tree <-  .cnode(1L, data, infl, inputs, weights, ctrl)
    fitted <- data.frame("(fitted)" = fitted_node(tree, data), 
                         "(response)" = data[ , response, drop = TRUE], 
                         check.names = FALSE)
    ret <- party(tree, data = data, fitted = fitted)
    class(ret) <- c("constparty", class(ret))
    return(ret)
}

.logrank_trafo <- function(x, ties.method = c("logrank", "HL")) {
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
.y2infl <- function(data, response) {

    if (length(response) == 1) {
        response <- data[[response]]
        rtype <- class(response)[1]
        if (rtype == "integer") rtype <- "numeric"

        infl <- switch(rtype,
            "factor" = model.matrix(~ response - 1),
            "ordered" = (1:nlevels(response))[as.integer(response)],
            "numeric" = response,
            "Surv" = .logrank_trafo(response)
        )
    } else {
        ### multivariate response
        infl <- lapply(response, .y2infl, data = data)
        infl <- do.call("cbind", infl)
    }
    storage.mode(infl) <- "double"
    return(infl)
}
