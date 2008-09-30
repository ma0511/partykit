
sapply(dir(path = "../R", pattern = "R$", full.names = TRUE), source)

foo <- function(x) partysplit(1L, breaks = as.double(x))
x <- vector(mode = "list", length = 3)
x[[1]] <- partynode(id = 1L, split = foo(1 / 3), kids = 2:3, info = "one")
x[[2]] <- partynode(id = 2L, info = "two")
x[[3]] <- partynode(id = 3L, split = foo(2 / 3), kids = 4:5, info = "three")
x[[4]] <- partynode(id = 4L, info = "four")
x[[5]] <- partynode(id = 5L, info = "five")

rx <- flat2rec(x)
stopifnot(identical(rec2flat(rx), x))

dat <- data.frame(x = runif(100))
kidids_node(rx, dat)
