
sapply(dir(path = "../R", pattern = "R$", full = TRUE), source)

library("rpart")
fit <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis)

pfit <- as.party(fit)

all(do_nodeid(pfit$node, kyphosis[, pfit$metadata$varnames]) == fit$where)

library("RWeka")

itree <- J48(Species ~ ., data = iris)
pitree <- as.party(itree)
all(predict(pitree) == predict(pitree, newdata = iris[, 3:4]))