
sapply(dir(path = "../R", pattern = "R$", full = TRUE), source)

library("rpart")
fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)

pfit <- as.party(fit)

all(predict(pfit, newdata = kyphosis, type = "node") == fit$where)

library("RWeka")

itree <- J48(Species ~ ., data = iris)
pitree <- as.party(itree)
all(predict(pitree) == predict(pitree, newdata = iris[, 3:4]))

all.equal(predict(itree, type = "prob", newdata = iris), predict(pitree, type = "prob", newdata = iris))
all.equal(predict(itree,  newdata = iris), predict(pitree, newdata = iris))

data("GlaucomaM", package = "ipred")

w <- runif(nrow(GlaucomaM))
fit <- rpart(Class ~ ., data = GlaucomaM, weights = w)
pfit <- as.party(fit)
all(predict(pfit, type = "node") == fit$where)
tmp <- GlaucomaM[sample(1:nrow(GlaucomaM), 100),]
all.equal(predict(fit, type = "prob", newdata = tmp), predict(pfit, type = "prob", newdata = tmp))
all.equal(predict(fit, type = "class", newdata = tmp), predict(pfit, newdata = tmp))                             

itree <- J48(Class ~ ., data = GlaucomaM)
pitree <- as.party(itree)
all.equal(predict(itree, newdata = tmp, type = "prob"), predict(pitree, newdata = tmp, type = "prob"))


data("airquality")
aq <- subset(airquality, !is.na(Ozone))

w <- runif(nrow(aq), max = 3)
aqr <- rpart(Ozone ~ ., data = aq, weights = w)
aqp <- as.party(aqr)

tmp <- subset(airquality, is.na(Ozone))
all.equal(predict(aqr, newdata = tmp), predict(aqp, newdata = tmp))