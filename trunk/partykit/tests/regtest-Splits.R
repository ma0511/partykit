
source("../R/Splits.R")

dat <- data.frame(v1 = as.double(1:100))

sv1 <- new_split(as.integer(1), breaks = as.double(50))
nodelabels(sv1, metadata(dat))
stopifnot(all(do_split(dat, sv1) == ((dat$v1 > 50) + 1)))

sv1 <- new_split(as.integer(1), breaks = as.double(50), 
                index = as.integer(c(2, 1)))
nodelabels(sv1, metadata(dat))
stopifnot(all(do_split(dat, sv1) == ((dat$v1 <= 50) + 1)))

sv1 <- new_split(as.integer(1), breaks = as.double(50), right = FALSE)
nodelabels(sv1, metadata(dat))
stopifnot(all(do_split(dat, sv1) == ((dat$v1 >= 50) + 1)))

sv1 <- new_split(as.integer(1), breaks = as.double(50), 
                index = as.integer(c(2, 1)), right = FALSE)
nodelabels(sv1, metadata(dat))
stopifnot(all(do_split(dat, sv1) == ((dat$v1 < 50) + 1)))

sv1 <- new_split(as.integer(1), breaks = as.double(c(25, 75)))
nodelabels(sv1, metadata(dat))
stopifnot(all(do_split(dat, sv1) == 
              as.integer(cut(dat$v1, c(-Inf, 25, 75, Inf)))))

sv1 <- new_split(as.integer(1), breaks = as.double(c(25, 75)), right = FALSE)
nodelabels(sv1, metadata(dat))
stopifnot(all(do_split(dat, sv1) == 
              as.integer(cut(dat$v1, c(-Inf, c(25, 75), Inf), right = FALSE))))

sv1 <- new_split(as.integer(1), breaks = as.double(c(25, 75)), 
                index = as.integer(3:1), right = FALSE)
nodelabels(sv1, metadata(dat))
stopifnot(all(do_split(dat, sv1) == 
              (3:1)[as.integer(cut(dat$v1, c(-Inf, c(25, 75), Inf), right = FALSE))]))


dat$v2 <- gl(4, 25)

sv2 <- new_split(as.integer(2), index = as.integer(c(1, 2, 1, 2)))
nodelabels(sv2, metadata(dat))
do_split(dat, sv2)

sv2 <- new_split(as.integer(2), breaks = as.integer(c(1, 3)))
nodelabels(sv2, metadata(dat))
do_split(dat, sv2)













dat <- data.frame(x = gl(3, 30, labels = LETTERS[1:3]), y = rnorm(90), 
                  z = gl(9, 10, labels = LETTERS[1:9], ordered = TRUE))
csp <- new_split(as.integer(1), index = as.integer(c(1, 2, 1)))
do_split(dat, csp)
do_splitlist(dat, list(csp))

nsp <- new_split(as.integer(2), breaks = c(-1, 0, 1), index = as.integer(c(1, 2, 1, 3)))
do_split(dat, nsp)

osp <- new_split(as.integer(3), breaks = as.integer(c(3, 6)), index = as.integer(c(2, 1, 2)))
do_split(dat, osp)

nadat <- dat
nadat$x[1:10] <- NA
nadat$y[11:20] <- NA
do_splitlist(nadat, list(csp, nsp, osp))

nodelabels(csp, metadata(dat))
nodelabels(nsp, metadata(dat))
nodelabels(osp, metadata(dat))

