
sapply(dir(path = "../R", pattern = "R$", full = TRUE), source)

library("rpart")
fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)

pfit <- as.party(fit)

all(do_nodeid(pfit$node, kyphosis[, pfit$metadata$varnames]) == fit$where)

library("RWeka")

itree <- J48(Species ~ ., data = iris)
pitree <- as.party(itree)
all(predict(pitree) == predict(pitree, newdata = iris[, 3:4]))


data("GlaucomaM", package = "ipred")

fit <- rpart(Class ~ ., data = GlaucomaM)
pfit <- as.party(fit)
all(do_nodeid(pfit$node, GlaucomaM[, pfit$metadata$varnames]) == fit$where)

itree <- J48(Class ~ ., data = GlaucomaM)
pitree <- as.party(itree)
all(predict(pitree) == predict(pitree, newdata = GlaucomaM))

