
### compare partykit::ctree and party::ctree wrt predictions
### and computing time. But not during R CMD check

library("partykit")
library("party")

set.seed(29)

n <- 200
dgp <- function(n, response = c("numeric", "ordered", "binary", "factor"), na = FALSE) {

    response <- match.arg(response)
    x <- runif(n * 10)
    x <- matrix(x, ncol = 10)
    x <- cbind(x, x[,1:2] + cbind(rnorm(n, sd = 0.1), rnorm(n, sd = 0.1)))
    x1 <- sample(gl(4, n / 4))
    x2 <- sample(ordered(gl(4, n / 4)))
    y <- as.integer(interaction(x[,1] > 0.5, x[,2] > 0.25))
    y <- rnorm(n, mean = y)
    if (na) x[sample(1:(n * 10), n)] <- NA
    if (response == "numeric")
        return(data.frame(x, x1, x2, y = y))
    if (response == "binary") {
        y <- factor(y < median(y))
        return(data.frame(x, x1, x2, y = y))
    } else {
        y <- cut(y, 
            breaks = c(-Inf, quantile(y, c(0.25, 0.5, 0.75)), Inf))
    }
    if (response == "ordered")
        return(data.frame(x, x1, x2, y = ordered(y)))
    if (response == "factor")
        return(data.frame(x, x1, x2, y = y))
}

args <- list(response = c("numeric", "ordered", "binary", "factor"),
             na = c(FALSE, TRUE),
             teststat = c("quad", "max"),
             testtype = c("Bonferroni", "Univariate", "Teststatistic"),
             mincriterion = c(0.8, 0.9, 0.95),
             maxsurrogate = 0:2, maxdepth = 3)

gr <- do.call("expand.grid", args)
tgr <- gr
gr[["response"]] <- as.character(gr[["response"]])
gr[["teststat"]] <- as.character(gr[["teststat"]])
gr[["testtype"]] <- as.character(gr[["testtype"]])

### random splits are done differently in party::ctree
gr <- subset(gr, !na | (na & maxsurrogate > 0))


fun <- function(args) {

    learn <- dgp(n, response = args$response, na = args$na)
    test <- dgp(n, response = args$response, na = args$na)

    ct <- party:::ctree_control
    ctrl <- do.call("ct", args[-(1:2)])
    ot <- system.time(oldmod <- party:::ctree(y ~ ., data = learn, control = ctrl))[1]
    old <- treeresponse(oldmod, newdata = test)
    ct <- partykit:::ctree_control
    ctrl <- do.call("ct", args[-(1:2)])
    nt <- system.time(newmod <- try(partykit:::ctree(y ~ ., data = learn, control = ctrl)))[1]
    if (inherits(newmod, "try-error")) {
        new <- as.list(rep(1, length(old)))
    } else {
        new <- predict(newmod, newdata = test, simp = FALSE, 
            type = ifelse(args$response == "numeric", "response", "prob"))
    }
    list(error = max(abs(unlist(old) - unlist(new))), time = c(ot, nt))
}

ret <- lapply(1:nrow(gr), function(i) { print(i); fun(gr[i, ,drop = FALSE]); })

tm <- t(sapply(ret, function(x) x$time))
err <- sapply(ret, function(x) x$error)

save(ret, tm, tgr, fun, gr, dgp, err, file = "results-regtest.Rda")

### scores
y <- gl(3, 10, ordered = TRUE)
x <- rnorm(length(y))
x <- ordered(cut(x, 3))
d <- data.frame(y = y, x = x)

### partykit with scores
ct11 <- partykit::ctree(y ~ x, data = d)
ct12 <- partykit::ctree(y ~ x, data = d, 
                        scores = list(y = c(1, 4, 5)))
ct13 <- partykit::ctree(y ~ x, data = d, 
                        scores = list(y = c(1, 4, 5), x = c(1, 5, 6)))

### party with scores
ct21 <- party::ctree(y ~ x, data = d)
ct22 <- party::ctree(y ~ x, data = d, 
                     scores = list(y = c(1, 4, 5)))
ct23 <- party::ctree(y ~ x, data = d, 
                     scores = list(y = c(1, 4, 5), x = c(1, 5, 6)))

stopifnot(all.equal(ct11$node$info$p.value, 
          1 - ct21@tree$criterion$criterion, check.attr = FALSE))
stopifnot(all.equal(ct12$node$info$p.value, 
          1 - ct22@tree$criterion$criterion, check.attr = FALSE))
stopifnot(all.equal(ct13$node$info$p.value, 
          1 - ct23@tree$criterion$criterion, check.attr = FALSE))

### ytrafo
y <- runif(100, max = 3)
x <- rnorm(length(y))
d <- data.frame(y = y, x = x)

### partykit with scores
ct11 <- partykit::ctree(y ~ x, data = d)
ct12 <- partykit::ctree(y ~ x, data = d,
                        ytrafo = list(y = sqrt))

### party with scores
ct21 <- party::ctree(y ~ x, data = d)
f <- function(data) ptrafo(data, numeric_trafo = sqrt)
ct22 <- party::ctree(y ~ x, data = d,
                     ytrafo = f)

stopifnot(all.equal(ct11$node$info$p.value,
          1 - ct21@tree$criterion$criterion, check.attr = FALSE))
stopifnot(all.equal(ct12$node$info$p.value,
          1 - ct22@tree$criterion$criterion, check.attr = FALSE))


### spotted by Peter Philip Stephensen (DREAM) <PSP@dreammodel.dk>
### splits x >= max(x) where possible in partykit::ctree
library("partykit")
nAge = 30
d = data.frame(Age=rep(1:nAge,2),y=c(rep(1,nAge),rep(0,nAge)),
n=rep(0,2*nAge))
ntot = 100
alpha = .5
d[d$y==1,]$n = floor(ntot * alpha * d[d$y==1,]$Age / nAge)
d[d$y==0,]$n = ntot - d[d$y==1,]$n
ctrl = partykit::ctree_control(maxdepth=3, minbucket = min(d$n) + 1)
tree = partykit::ctree(y ~ Age, weights=n, data=d, control=ctrl)
tree

(w1 <- predict(tree, type = "node"))

detach(package:partykit)
library("party")
ctrl = party::ctree_control(maxdepth=3, minbucket = min(d$n) + 1)
tree2 = party::ctree(y ~ Age, weights=d$n, data=d, control=ctrl)
tree2

(w2 <- where(tree2))

stopifnot(max(abs(w1 - w2)) == 0)

gr <- cbind(gr, err = err)
subset(gr, !na & err > 0.01)
