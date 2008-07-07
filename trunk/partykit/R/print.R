

print.node <- function(x, metadata, prefix = "", first = TRUE, ...) {

    ### FIXME: process info slot
    if (first)
        cat(paste(prefix, "[", x$id, "] root\n", sep = ""))

    if (length(x) > 0) {

        labs <- nodelabels(x$split[[1]], metadata, ...)
        labs <- ifelse(substr(labs$splitlevels, 1, 1) %in% c("<", ">"),
                       paste(labs$splitname, labs$splitlevels), 
                       paste(labs$splitname, "in", labs$splitlevels))

        ### FIXME: use id labels instead of raw ids
        terminal <- sapply(1:length(x), function(z) length(x[[z]]) < 1)
        labs <- paste("  ", prefix, "[", sapply(1:length(x), 
           function(z) x[[z]]$id), "] ", labs, 
           ifelse(terminal, " *", ""), "\n", sep = "")
                  
        for (i in 1:length(x)) {
            cat(labs[i])
            print(x[i], metadata, prefix = paste(prefix, "  ", sep = ""), 
                  first = FALSE, ...)
        }
    }
}

print.party <- function(x, ...)
    print(x$node, x$metadata, ...)
