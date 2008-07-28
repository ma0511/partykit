## FIXME: think about whether this should be visible or not,
## or maybe an S3 method or not.
print_node <- function(x, metadata, names, prefix = "", leaf = " *", first = TRUE, ...) {

    ### FIXME: process info slot
    if (first)
        cat(paste(prefix, "[", names[get_id(x)], "] root\n", sep = ""))

    if (length(x) > 0) {
        labs <- nodelabels(get_primary(get_split(x)), metadata, ...)
        labs <- ifelse(substr(labs$splitlevels, 1, 1) %in% c("<", ">"),
                       paste(labs$splitname, labs$splitlevels), 
                       paste(labs$splitname, "in", labs$splitlevels))

        ### FIXME: use id labels instead of raw ids
        terminal <- sapply(1:length(x), function(z) length(x[[z]]) < 1)
        labs <- paste("|   ", prefix, "[", sapply(1:length(x), 
           function(z) names[get_id(x[[z]])]), "] ", labs, 
           ifelse(terminal, leaf, ""), "\n", sep = "")
                  
        for (i in 1:length(x)) {
            cat(labs[i])
            print_node(x[i], metadata, names = names, prefix = paste(prefix, "|   ", sep = ""), 
                  leaf = leaf, first = FALSE, ...)
        }
    }
}

print.party <- function(x, ...)
    print_node(x$node, x$metadata, names = get_names(x), ...)
