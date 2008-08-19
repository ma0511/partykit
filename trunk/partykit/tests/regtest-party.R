
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

data("GBSG2", package = "ipred")
library("survival")
fit <- rpart(Surv(time, cens) ~ ., data = GBSG2)
pfit <- as.party(fit)
pfit$fitted
predict(pfit, newdata = GBSG2[1:100,], type = "prob")
predict(pfit, newdata = GBSG2[1:100,], type = "response")

dgp <- function(n)
    data.frame(y = gl(4, n), x1 = rnorm(4 * n), x2 = rnorm(4 * n))

learn <- dgp(100)
fit <- as.party(rpart(y ~ ., data = learn))
test <- dgp(100000)
system.time(id <- fitted_node(node_party(fit), test))
system.time(yhat <- predict_party(fit, id = id, newdata = test))

### multiple responses
f <- fitted(pfit)
f[["(response)"]] <- data.frame(srv = f[["(response)"]], hansi = runif(nrow(f)))
mp <- party(node_party(pfit), fitted = f, data = pfit$data)
class(mp) <- c("cparty", "party")

predict(mp, newdata = GBSG2[1:10,])

### predictions in info slots
tmp <- data.frame(x = rnorm(100))
pfit <- party(node = node(1L, split = split(1L, breaks = 0), 
              kids = list(node(2L, info = -0.5), node(3L, info = 0.5))), data = tmp)
pfit
p <- predict(pfit, newdata = tmp)
p
table(p, sign(tmp$x))
