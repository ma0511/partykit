## FIXME: think about whether this should be visible or not,
## or maybe an S3 method or not.
character_node <- function(x, data, names, prefix = "", leaf = " *", first = TRUE, ...) {

    ### FIXME: process info slot
    if (first)
        cat(paste(prefix, "[", names[id_node(x)], "] root\n", sep = ""))

    if (length(x) > 0) {
        labs <- character_split(split_node(x), data, ...)
        labs <- ifelse(substr(labs$levels, 1, 1) %in% c("<", ">"),
                       paste(labs$name, labs$levels), 
                       paste(labs$name, "in", labs$levels))

        ### FIXME: use id labels instead of raw ids
        terminal <- sapply(1:length(x), function(z) length(x[[z]]) < 1)
        labs <- paste("|   ", prefix, "[", sapply(1:length(x), 
           function(z) names[id_node(x[[z]])]), "] ", labs, 
           ifelse(terminal, leaf, ""), "\n", sep = "")
                  
        for (i in 1:length(x)) {
            cat(labs[i])
            character_node(x[i], data, names = names, prefix = paste(prefix, "|   ", sep = ""), 
                  leaf = leaf, first = FALSE, ...)
        }
    }
}

print.party <- function(x, ...)
    character_node(node_party(x), x$data, names = names_party(x), ...)
