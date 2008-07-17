
sapply(dir(path = "../R", pattern = "R$", full = TRUE), source)

library("rpart")
fit <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis)

pfit <- as.party(fit)

all(do_nodeid(pfit$node, kyphosis[, pfit$metadata$varnames]) == fit$where)

