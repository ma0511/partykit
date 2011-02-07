.initializeNode <- function(mtree){
        # initializes the tree structure stored in node
	temp <- NULL
	if( max(ls(all.names = TRUE) == ".evtree_gid") == 1 ) # just in case .evtree_gid is already defined
		temp <- .evtree_gid 
        .evtree_gid <<- NULL
        node <- .initNode(id = 1L, mtree)
        evtree_gid <- .evtree_gid
        rm(.evtree_gid, inherits = TRUE)
	if(is.null(temp) == FALSE)	
	       .evtree_gid <<- temp
        return(list(node, evtree_gid))
}

.initNode <- function(id, mtree){
        ## is called by initializeNode() 
        ## should not be used directly
        if(is.null(.evtree_gid)){
	     .evtree_gid <<- 1
        }else{
             .evtree_gid <<- c(.evtree_gid, id)
	}
        if(mtree$varType[mtree$splitV[id]] > 0){
              split = partysplit(as.integer(mtree$splitV[id]), breaks = mtree$splitP[id], right = FALSE)
        }else{
              index <- as.integer(mtree$csplit[( (id-1)*(mtree$maxCat)+1 ):( (id-1)*(mtree$maxCat)+1  + abs(mtree$varType[mtree$splitV[id]])-1 )])
              index[index == 2] <- NA
              index[index == 3] <- as.integer(2)
              breaks = NULL
              split = partysplit(as.integer(mtree$splitV[id]), index = index)
        }
        node <- partynode(id = length(.evtree_gid), split = split,
        if( id*2 < 2^(mtree$maxdepth) ){
              kids = list(
                        if( mtree$splitN[id*2] == id*2 ){
                            .initNode(id = id*2, mtree)
                        }else{
                            .evtree_gid <<- c(.evtree_gid, id*2)
                            partynode(as.integer(length(.evtree_gid)))
                        },
                        if( mtree$splitN[id*2+1] == id*2+1 ) {
                            .initNode(id = id*2+1, mtree)
                        }else{
                            .evtree_gid <<- c(.evtree_gid, id*2+1)
                            partynode(as.integer(length(.evtree_gid)))
                        }
                        )
        }else{
              .evtree_gid <<- c(.evtree_gid, id*2)
              .evtree_gid <<- c(.evtree_gid, id*2+1)
              kids = list(
                        partynode(as.integer(length(.evtree_gid)-1)),
                        partynode(as.integer(length(.evtree_gid)  ))
                        )
        }
    	)
}

