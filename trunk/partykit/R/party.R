
new_party <- function(node, metadata) {

    stopifnot(inherits(node, "node"))
    stopifnot(inherits(metadata, "metadata"))

    if (is.flat(node))
        stop(sQuote("node"), " ", "is not a recursive node")

    party <- list(node = node, metadata = metadata)
    class(party) <- "party"
    party
}
