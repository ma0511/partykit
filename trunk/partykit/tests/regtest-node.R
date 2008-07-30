
sapply(dir(path = "../R", pattern = "R$", full.names = TRUE), source)

foo <- function(x) list(split(1L, breaks = as.double(x)))
x <- vector(mode = "list", length = 3)
x[[1]] <- node(id = 1L, split = foo(1 / 3), kids = 2:3, info = "one")
x[[2]] <- node(id = 2L, info = "two")
x[[3]] <- node(id = 3L, split = foo(2 / 3), kids = 4:5, info = "three")
x[[4]] <- node(id = 4L, info = "four")
x[[5]] <- node(id = 5L, info = "five")

rx <- flat2rec(x)
stopifnot(identical(rec2flat(rx), x))

dat <- data.frame(x = runif(100))
kidids_node(rx, dat)
