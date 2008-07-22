
sapply(dir(path = "../R", pattern = "R$", full.names = TRUE), source)

foo <- function(x) list(new_split(as.integer(1), breaks = as.double(x)))
x <- vector(mode = "list", length = 3)
x[[1]] <- new_node(id = 1:1, split = foo(1 / 3), kids = 2:3, info = "one")
x[[2]] <- new_node(id = 2:2, info = "two")
x[[3]] <- new_node(id = 3:3, split = foo(2 / 3), kids = 4:5, info = "three")
x[[4]] <- new_node(id = 4:4, info = "four")
x[[5]] <- new_node(id = 5:5, info = "five")

rx <- flat2rec(x)
stopifnot(identical(rec2flat(rx), x))

dat <- data.frame(x = runif(100))
do_nodeid(rx, dat)
