evtree.control <- function(minbucket = 7L, minsplit = 20L, maxdepth = 9L, niterations = 10000L, ntrees = 100L, alpha = 1,
  operatorprob= list(pmutatemajor = 0.2, pmutateminor = 0.2, pcrossover = 0.2, psplit = 0.2, pprune = 0.2), seed = NULL, ...){
    args <- list(...)
    if( length(args) > 0 ){   
	for(a in 1:length(args))
		if(names(args)[a] == "minbucket")
			minbucket <- args[[a]]
		else if(names(args)[a] == "minsplit")
			minsplit <- args[[a]]
		else if(names(args)[a] == "maxdepth")
			maxdepth <- args[[a]]
		else if(names(args)[a] == "niterations")
			niterations <- args[[a]]
		else if(names(args)[a] == "ntrees")
			ntrees <- args[[a]]
		else if(names(args)[a] == "alpha")
			alpha <- args[[a]]
		else if(names(args)[a] == "seed")
			seed <- args[[a]]
		else if(names(args)[a] == "operatorprob"){
		  	if(is.list(args)[a]){
		  	          op <- args[[a]]
		  	      for(i in 1:length(op))
			      if(names(op)[a] == "pmutatemajor")
		    	          operatorprob$pmutatemajor <- op[i]
			      else if(names(op)[a] == "pmutateminor")
			          operatorprob$pmutateminor <- op[i]
			      else if(names(op)[a] == "pcrossover")
				  operatorprob$pcrossover <- op[i]
			      else if(names(op)[a] == "psplit")
				  operatorprob$psplit <- op[i]
			      else if(names(op)[a] == "pprune")
			          operatorprob$pprune <- op[i]	  	      
		  	      else
			          warning(paste("extra argument ", names(op)[i], " in \"operatorprob\" is just disregarded"))
		  	}else if(length(args)[a] == 5){
  		  	      op <- list()
    	  		      op$pmutatemajor <- operatorprob[[1]]
    	    		      op$pmutateminor <- operatorprob[[2]]
    	     	  	      op$psplit <- operatorprob[[3]]
    	         	      op$pprune <- operatorprob[[4]]
    	          	      op$pcrossover <- operatorprob[[5]]
    	          	      operatorprob <- op
			}else{
   	   		      stop("argument operatorprob must be of form: operatorprob= list(pmutatemajor = 0.2, pmutateminor = 0.2, pcrossover = 0.2, psplit = 0.2, pprune = 0.2), or a vector with 5 elements")
   	   	        }
   	        }
		else if(names(args)[a] == "pmutatemajor")
			operatorprob$pmutatemajor <- args[[a]]
		else if(names(args)[a] == "pmutateminor")
			operatorprob$pmutateminor <- args[[a]]
		else if(names(args)[a] == "pcrossover")
			operatorprob$pcrossover <- args[[a]]
		else if(names(args)[a] == "psplit")
			operatorprob$psplit <- args[[a]]
		else if(names(args)[a] == "pprune")
			operatorprob$pprune <- args[[a]]
		else
			warning(paste("extra argument ", names(args)[a], " is just disregarded"))
     }
  
      if( is.null(operatorprob$pmutatemajor) )
  	   operatorprob$pmutatemajor <- 0.2
      if( is.null( operatorprob$pmutateminor) )
            operatorprob$pmutateminor<- 0.2
      if( is.null( operatorprob$pprune) )
	    operatorprob$pprune<- 0.2
      if( is.null(operatorprob$pcrossover) )
    	    operatorprob$pcrossover<- 0.2
      if( is.null( operatorprob$psplit) )
  	    operatorprob$psplit<- 0.2
 

      sumop <- operatorprob$pmutatemajor + operatorprob$pmutateminor + operatorprob$pprune + operatorprob$pcrossover + operatorprob$psplit
      if(sumop != 1 && sumop != 100)      
             warning('user defined operator probabilities do not sum up to 1')
      operatorprob$pmutatemajor <- floor(100* operatorprob$pmutatemajor / sumop)
      operatorprob$pmutateminor <- floor(100* operatorprob$pmutateminor / sumop)
      operatorprob$pprune <-  floor(100*operatorprob$pprune / sumop)
      operatorprob$pcrossover <-  floor(100*operatorprob$pcrossover  / sumop)
      operatorprob$psplit <- 100 - (operatorprob$pmutatemajor + operatorprob$pmutateminor + operatorprob$pprune + operatorprob$pcrossover)
  	    	     
      if(minsplit < 2 | minbucket < 1 )
             stop("parameters \"minsplit\" and \"minbucket\" must be defined with at least 1 and 2 observations")

      if(minsplit < minbucket)
 	     minsplit <- minbucket+1

      if(is.null(seed)){
   	     seed = -1
     }else if(seed < 0){
  	     stop("parameter \"seed\" in evtree.control must be a non-negative integer")
     }

     if(niterations < 100)
             stop("parameter \"niterations\" must at least be 100")

     if(ntrees < 10)
             stop("parameter \"ntrees\" must at least be 10")

     if(alpha < 0)
             stop("parameter \"alpha\" must be equal or larger than 0")

     if(maxdepth > 12)
            stop("\"maxdepth > 12\" is not supported")
  
     list(minbucket = minbucket,
   	    minsplit = minsplit,
   	    maxdepth = maxdepth,
   	    niterations = niterations,
            ntrees = ntrees,
            alpha = alpha,
            operatorprob = operatorprob,
	    seed = seed
        )
}

