
evtree <- function(formula, data = list(), weights = NULL, subset = NULL, control = evtree.control()){
    start <- proc.time()
    call <- match.call()
    call[[1L]] <- as.name("model.frame")
    call$control <- NULL
    call <- eval(call, parent.frame())
    terms <- (attr(call,"terms"))
    if(any(attr(terms, "order") > 1L))
	stop("Trees cannot handle interaction terms")

    mf <- model.frame(formula = formula, data = data, drop.unused.levels = FALSE, na.action = NULL)
    mf <- mf[,c(2:ncol(mf),1)]

    if(length(unique(mf[,ncol(mf)])) == 1)
        stop('dependend variable has no variation')

    mf <- mf[,sapply(sapply(mf, unique), length) > 1] ## drop variables with only one level
    if(is.null(dim(mf)))
        stop('independend variables are all constant')

    if(is.null(weights)){
        weights <- array(1,nrow(mf))
    }else{
       if(sum(weights %% 1) != 0)
         stop('weights must be integers')
       weights <- as.integer(weights)
       if(is.numeric(weights) == FALSE)
         stop('weights must be of type numeric')
       if(min(weights < 0) == TRUE)
         stop('weights must be defined as a positive number')
       if((min(weights) == 0) == TRUE){
         mf <- subset(mf, weights > 0)
         weights <- subset(weights, weights > 0)
       }
       if(nrow(mf) != length(weights))
         stop('length of "weights" does not match training data')
    }

    if(is.null(subset) == FALSE){
       mf <- mf[subset,]
       weights <- weights[subset]
    }

    if( sum(is.na(mf)) > 0 ){
        warning(sum(FALSE==(complete.cases(mf))), " rows are removed due to missing values" )
        # remove weights of incomplete cases
        weights <- weights[complete.cases(mf)]
        # remove missing values
        mf <- na.omit(mf)
    }

    nVariables <- ncol(mf)
    nInstances <- nrow(mf)
    prediction <- array(as.integer(0),nInstances)

    if(is.null(control$method)){
         if(is.factor(mf[,nVariables])){
            control$method <- 1
         }else{
            control$method <- 6
         }
    }

     if(control$method < 6){
        if(is.factor(mf[,nVariables]) == FALSE)
            stop('dependend variable is not a factor')
        if(length(levels(mf[,nVariables]))<2)
            stop("dependend variable has only", length(levels(mf[,nVariables])), " level(s)")
    }else{
        if(is.factor(mf[,nVariables]) == TRUE)
            stop('dependend variable is a factor')
        if(var(mf[,nVariables]) <= 0)
            stop("variance of the denpendend variable is 0")
    }

    if(is.null(control$operatorprob$pmutatemajor) == FALSE && is.null(control$operatorprob$pmutateminor) == FALSE &&
        is.null(control$operatorprob$psplit) == FALSE && is.null(control$operatorprob$pprune) == FALSE &&
        is.null(control$operatorprob$pcrossover) == FALSE
        ){
            op <- control$operatorprob
    }else if(is.null(control$operatorprob$pmutatemajor) && is.null(control$operatorprob$pmutateminor)&&
        is.null(control$operatorprob$psplit) && is.null(control$operatorprob$pprune) && is.null(control$operatorprob$pcrossover)
        ){
            op <- list()
            op$pmutatemajor <- control$operatorprob[[1]]
            op$pmutateminor <- control$operatorprob[[2]]
            op$psplit <- control$operatorprob[[3]]
            op$pprune <- control$operatorprob[[4]]
            op$pcrossover <- control$operatorprob[[5]]
    }else{
       stop('argument operatorprob must be of form: operatorprob=
list(pmutatemajor = 20, pmutateminor = 20, pcrossover = 20, psplit = 20, pprune = 20)')
    }

    if(op$pmutatemajor + op$pmutateminor + op$psplit + op$pprune+op$pcrossover != 100)
        stop("elements of operatorprob must sum up to 100%")

    if(control$minsplit < 2 | control$minbucket < 1 )
        stop("parameters \"minsplit\" and \"minbucket\" must be defined with at least 1 and 2 observations")

    if(control$minsplit < control$minbucket)
        control$minsplit <- control$minbucket+1

    if(control$minsplit > sum(weights))
        stop(paste("no split could be found \n \"minsplit\" is larger than the weighted number of observations in the training data"))

    if(2*control$minbucket > sum(weights)-1)
        stop(paste("no split could be found \n \"minbucket\" is larger than half the weighted number of observations in the training data"))

    if(is.null(control$seed)){
   	control$seed = -1
    }else if(control$seed < 0){
        stop("parameter \"seed\" in evtree.control must be a non-negative integer")
    }

    if(control$niterations < 100)
        stop("parameter \"niterations\" must at least be 100")

    if(control$ntrees < 10)
        stop("parameter \"ntrees\" must at least be 100")

    if(control$alpha < 0)
        stop("parameter \"alpha\" must be equal or larger than 0")

    if(control$maxdepth > 12)
        stop("\"maxdepth > 12\" is not supported")

    for (i in 1:(nVariables-1)){
        if(is.character(mf[,i])){
           mf[,i] <- as.factor( mf[,i] )
           warning(paste("character variable", names(mf[i]),"was converted to a factor"))
        }
        if(is.logical(mf[,i]))
           mf[,i] <- as.factor( mf[,i] )
        if(class(mf[,i])[1] == "Date")
           stop(paste("variable type \"Date\" of variable \"", names(mf[i]), "\" is not supported" ), sep = '""')
    }

    # calculate the maximum number of internal nodes a tree with size control$maxdepth can have
    maxNode <- 1
    if (control$maxdepth > 1)
    for(i in 1:(control$maxdepth-1) ){
        maxNode <- maxNode*2+1
    }

    # splitnumbers, splitvariables and splitpoints
    splitN <- array(-999999, maxNode)
    splitV <- array(-999999, maxNode)
    splitP <- array(-999999, maxNode)

    # mf is transformed into a one dimensional vector such it can be passed to the C++ routines
    ndata <- array(-1, (nInstances*nVariables))
    # maxCat is the maximum number of categories of the variables
    maxCat <- 0
    # varType is negative for nominal variables and positive for ordered variables
    # the absolute value of varType[i] is the number of destinct values of variable i
    varType <-  array(-1, nVariables)
    j <- 1
    k <- 1

    # calculation of varType and maxCat
    for (i in 1:nVariables ){
        ndata[  (1+   nInstances*(i-1) ) : ( nInstances[1]+    nInstances*(i-1) )  ] <- mf[,i]
        if( is.factor(mf[,i]) ){
             if(is.ordered(mf[,i])){
                varType[i] <- length(levels(mf[,i]))
                j <- j + j
            }else{
                varType[i] <- -length(levels(mf[,i]))
                k <- k + 1
            }
            if(abs(varType[i]) > maxCat & i < nVariables)
                maxCat <- abs(varType[i])
        }else{
            varType[i] <- length(unique(mf[,i]))
            j <- j +1
        }
    }
    # for categorical splits.
    csplit <- array(2, maxCat*maxNode)

    #specification of function tree
    tree <- function(nInstances,nVariables,varType,ndata,weights,prediction,splitN,splitV,splitP,
                csplit,maxNode, op, control){
            out <-  .C( "tree", PACKAGE="evtree",
                as.integer(nInstances), 
                as.integer(nVariables), 
                as.integer(varType), 
                as.double(ndata), 
                as.integer(weights),
                as.integer(prediction), 
                as.integer(splitN), 
                as.integer(splitV), 
                as.double(splitP), 
                as.integer(csplit), 
                as.integer(maxNode), 
                as.integer(control$minbucket), 
                as.integer(control$minsplit), 
                as.integer(control$niterations), 
                as.integer(control$ntrees),
                as.integer(op$pmutatemajor),
                as.integer(op$pmutateminor),
                as.integer(op$pcrossover),
                as.integer(op$psplit),
                as.integer(op$pprune),
                as.integer(control$method),
                as.double(control$alpha),
		as.integer(control$seed)
                )
            return(out)
    }
    #Call of the tree function
    out <- tree(nInstances, nVariables, varType, ndata, weights, prediction, splitN, splitV, splitP, csplit,
            maxNode, op, control)
        mtree = list()
        mtree$varType <-varType
        mtree$splitN <- out[[7]]
        mtree$splitV <- out[[8]]
        mtree$splitP <- out[[9]]
        mtree$csplit <- out[[10]]
        mtree$maxdepth <- control$maxdepth
        mtree$weights <- weights
        mtree$prediction <- out[[6]]+1
        mtree$maxCat <- maxCat
        mtree$seed <- out[[23]]

        init <- .initializeNode(mtree)
        node <- init[[1]]
        gid  <- init[[2]]

        prediction <- array(-999999, length(mtree$prediction))
        for(i in 1:length(gid))
            prediction[ mtree$prediction == gid[i] ] <- i

        fitted <- data.frame(prediction, mtree$weights, mf[nVariables])
        names(fitted) <- c("(fitted)", "(weights)" ,"(response)")
        partyObject <- party(node, mf, fitted = fitted, terms = terms, info = list(method = "evtree", nIterations = out[[14]], seed = mtree$seed))
        class(partyObject) <- c("constparty", "party")
        return(partyObject)
}




