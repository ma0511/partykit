
setClassUnion("integer_OR_NULL", c("integer", "NULL"))

setClass("CatSplit", 
    representation = representation(
        node_ids = "integer_OR_NULL"),
    prototype = list(node_ids = NULL)
)

setClass("SimCatSplit", 
    contains = "CatSplit",
    representation = representation(
        var_id = "integer")
)

setClass("FunCatSplit", 
    contains = "CatSplit",
    representation = representation(
        fun = "function")
)

setClass("SimNumSplit", 
    contains = "SimCatSplit",
    representation = representation(
        partition = "numeric", 
        right     = "logical"),
        prototype(right = TRUE),
)

setClass("FunNumSplit", 
    contains = "FunCatSplit",
    representation = representation(
        partition = "numeric",
        right     = "logical"),
    prototype(right = TRUE),
    validity = function(object) {
        if (!is.null(object@node_ids))
            length(object@node_ids) == length(object@partition) + 1
    }
)


split <- function(data, splits, NAprob = NULL) {

    if (!is.list(splits))
        splits <- list(splits)

    Fun2Sim <- function(split) {

        if (extends(class(split), "FunCatSplit") ||
            extends(class(split), "FunNumSplit")) {

            p <- ncol(data)
            data[[p + 1]] <<- split@fun(data)
            if (is.factor(data[[p + 1]]))
                split <- new("SimCatSplit", var_id = as.integer(p + 1),
                             node_ids = split@node_ids)
            if (is.numeric(data[[p + 1]]))
                split <- new("SimNumSplit", var_id =  as.integer(p + 1),
                             partition = split@partition, node_ids = split@node_ids)
            split@var_id <- as.integer(p + 1)
        }
        return(split)
   }

   splits <- lapply(splits, Fun2Sim)

   x <- simsplit(data, splits[[1]])
   i <- 2

   for (i in 2:length(splits)) {
       if (!any(is.na(x))) break;
       x[is.na(x)] <- simsplit(data[is.na(x),], splits[[i]])
   }

   if (any(is.na(x))) {
       if (is.null(NAprob)) {
           x[is.na(x)] <- 1
       } else {
           x[is.na(x)] <- sample(1:nmax, sum(is.na(x)), FALSE, prob = NAprob)
       }
   }
   return(x)
}

simsplit <- function(data, split) {

    x <- data[[split@var_id]]
    if (extends(class(split), "SimNumSplit"))
        x <- cut(x, breaks = c(-Inf, split@partition, Inf), right = split@right)

    x <- as.integer(x)
    node_ids <- split@node_ids
    if (!is.null(node_ids))
        x <- node_ids[x]
    return(x)
}

NumSplit2char <- function(split) {

    partition <- split@partition
    node_ids <- split@node_ids

    if (is.null(node_ids)) 
        node_ids <- 1:(length(partition) + 1)

    r <- split@right
    if (length(partition) == 1) {
        if (node_ids[1] < node_ids[2]) {
            nodetxt <- c(paste(ifelse(r, "<=", "<"), partition), 
                         paste(ifelse(r, ">", ">="), partition))
        } else {
            nodetxt <- c(paste(ifelse(r, ">", ">="), partition), 
                         paste(ifelse("r", "<=", "<"), partition))
        }
    } else {
        txt <- cbind(c(-Inf, partition), c(partition, Inf))
        txt <- apply(txt, 1, function(x) 
            paste("(", paste(as.vector(x), collapse = ", "), "]", sep = ""))
        nodetxt <- paste("%in%", txt[node_ids])
    }
    nodetxt
}   

CatSplit2char <- function(split, levels = NULL) {

    node_ids <- split@node_ids

    if (is.null(levels)) 
        levels <- paste("c", 1:length(node_ids), sep = "")

    sapply(sort(unique(node_ids)), function(n) 
        paste("%in%", paste("{", 
              paste(levels[node_ids == n], collapse = ", "), "}", sep = ""))
    )
}



dat <- data.frame(x = rnorm(100), y = gl(4, 25), z = ordered(gl(4, 25)))

### categorical split
cs <- new("SimCatSplit", var_id = as.integer(2), 
          node_ids = as.integer(c(1, 2, 2, 1)))
split(dat, cs)

### binary numeric split
ns <- new("SimNumSplit", var_id = as.integer(1),
          partition = as.numeric(1.2))

split(dat, ns)

### binary split with direction
ns <- new("SimNumSplit", var_id = as.integer(1),
          partition = as.numeric(1.2), 
          node_ids = as.integer(c(2, 1)))

split(dat, ns)

### multiway split
ns <- new("SimNumSplit", var_id = as.integer(1),
          partition = as.numeric(c(0, 1.2)),
          node_ids = as.integer(c(1, 3, 2)))
split(dat, ns)

### functional split
fs <- new("FunNumSplit", fun = function(data) data[[1]]^2,
          partition = as.numeric(2))
split(dat, fs)

