
library("rpart")
library("RWeka")
library("partykit")

data("LetterRecognition", package = "mlbench")

### fit rpart and J48 trees
rp <- rpart(lettr ~ ., data = LetterRecognition)
j48 <- J48(lettr ~ ., data = LetterRecognition)

### convert to party
system.time(party_rp <- as.party(rp))
system.time(party_j48 <- as.party(j48))

### check depths
depth(party_rp)
depth(party_j48)

### compare object sizes
object.size(rp)
object.size(party_rp)

object.size(j48) ### tree stored in external pointer
object.size(party_j48)

### set-up large sample
i <- sample(1:nrow(LetterRecognition), 1000000, replace = TRUE)
x <- LetterRecognition[i,]

### compare predictions (speed and accuracy)
system.time(p_rp <- predict(rp, newdata = x, type = "prob"))
system.time(p_party_rp <- predict(party_rp, newdata = x, type = "prob"))
all.equal(p_rp, p_party_rp)

system.time(p_j48 <- predict(j48, newdata = x))
system.time(p_party_j48 <- predict(party_j48, newdata = x))
all.equal(p_j48, p_party_j48, check.attributes = FALSE)
