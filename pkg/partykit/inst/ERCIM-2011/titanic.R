## packages
library("partykit")
library("rpart")
library("RWeka")
library("CHAID")
library("evtree")

## data
data("Titanic", package = "datasets")
ttnc <- as.data.frame(Titanic)
ttnc <- ttnc[rep(1:nrow(ttnc), ttnc$Freq), 1:4]
names(ttnc)[2] <- "Gender"
ttnc <- transform(ttnc, First = factor(Gender == "Female" | Age == "Child",
  levels = c(FALSE, TRUE), labels = c("Male&Adult", "Female|Child")))

## rpart
rp <- rpart(Survived ~ Gender + Age + Class, data = ttnc)
plot(as.party(rp))

## ctree
ct <- ctree(Survived ~ Gender + Age + Class, data = ttnc,
  control = ctree_control(mincrit = 0.99, teststat = "max"))
plot(ct)

## J4.8
j48 <- J48(Survived ~ Gender + Age + Class, data = ttnc)
plot(as.party(j48))

## QUEST
qu <- pmmlTreeModel(file.path(system.file("pmml", package = "partykit"), "ttnc.pmml"))
plot(qu)

## CHAID
chd <- chaid(Survived ~ Gender + Age + Class, data = ttnc)
plot(chd)

## evtree
set.seed(1071)
ev <- evtree(Survived ~ Gender + Age + Class, data = ttnc, maxdepth = 3)
plot(ev)

## performance
mc <- function(obj) 1 - mean(predict(obj) == ttnc$Survived)
evalfun <- function(obj) 2 * nrow(ttnc) * mc(obj) + width(obj) * log(nrow(ttnc))
trees <- list(evtree = ev, rpart = as.party(rp), ctree = ct, J48 = as.party(j48))
round(sapply(trees, function(obj) c(
  "complexity" = width(obj),
  "misclassification" = mc(obj),
  "evaluation function" = evalfun(obj))), digits = 3)

## compare predictions
nd <- ttnc[rep(1:nrow(ttnc), 100), ]
system.time(p1_rp <- predict(rp, newdata = nd, type = "class"))
system.time(p2_rp <- predict(as.party(rp), newdata = nd))
table(p1_rp, p2_rp)

system.time(p1_j48 <- predict(j48, newdata = nd))
system.time(p2_j48 <- predict(as.party(j48), newdata = nd))
table(p1_j48, p2_j48)

## mob
mb <- glmtree(Survived ~ First | Class + Gender + Age, data = ttnc,
  family = binomial, alpha = 0.01)
plot(mb)

mb2 <- glmtree(Survived ~ First | Class + Gender + Age, data = ttnc,
  family = binomial, alpha = 0.01, catsplit = "multiway")
plot(mb2)
