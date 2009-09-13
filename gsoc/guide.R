
library(partykit)

ChiSquare <- function(x){
	# chi-square test on a table x 
	cs <- colSums(x) > 0
	rs <- rowSums(x) > 0
    	if (sum(cs) < 2 || sum(rs) < 2) return(1)
    	if (min(x) < 10 && sum(x) < 100) {
        	ctest <- chisq.test(x[rs, cs], correct = FALSE, 
                            simulate.p.value = TRUE, B = 9999)
        	X2 <- ctest$statistic
        	ret <- ctest$p.value
    	} else {
        	suppressWarnings(ctest <- chisq.test(x[rs, cs], correct = FALSE))
        	X2 <- ctest$statistic
        	df <- ctest$parameter
        	ret <- pchisq(X2, df = df, lower = FALSE)
    	}
	attr(ret, "Chisq") <- X2
	ret
}


splitvar <- function(response, xvars, weights, complexity, vartype) { 
	# implementation of p. 366 ff in Loh (2002)
	
	switch( complexity, 
                "const" = { # fitting a piecewise constant model:
	                const.fit <- lm( response[ weights > 0] ~ 1)
	                res <- resid(const.fit)	# Here: res_i = Y_i - mean(Y)
                },
                "mult" = { # fitting a multiple linear model:
	                df.nf <- data.frame(response, xvars[,vartype %in% c("n","f")])
	                mult.fit <- lm(df.nf[ weights > 0,])
	                res <- resid(mult.fit)
                },
                "simple" = { # fitting a best simple linear model:
                    ly <- sum(weights)  #assuming weights are either 0 (unused data) or 1 (used data)
                    firstcol <- rep.int(1,ly)
                    xvars.to.fit <- xvars[,vartype %in% c("n","f")]
                    resxvars <- apply( xvars.to.fit[weights > 0,], 2,
                                function(x){ lm.fit( cbind( firstcol,x),response[weights > 0])$residuals})
                    res <- resxvars[,which.min( colSums(resxvars^2))]
                }
    )

    to.split <- vartype %in% c("n","s","c")
    split.num <- sum(to.split)
    index.counter <- 1:length(xvars)
    index.counter <- index.counter[to.split]

    if( sum(to.split) == 1) return( list(id=as.integer(index.counter), residuals=res))

	p.curv <- numeric(split.num)
	X2.curv <- rep(NA, split.num)
	p.interact <- matrix( rep(NA, split.num^2), nrow=split.num)
	X2.interact <- matrix( rep(NA, split.num^2), nrow=split.num)
	pairwise <- lapply( 1:(split.num-1), function(i) (i + 1):split.num)
    for( i in 1:split.num ){
        tmp <- curv.test( xvars[,to.split][[i]], res, weights)
        p.curv[i] <- tmp
        if( !is.null(attr(tmp, "Chisq")) ){ X2.curv[i] <- attr(tmp, "Chisq")}
		if( i < split.num ){
			for( j in pairwise[[i]] ){
                tmp.interact <- interact.test( xvars[,to.split][[i]], xvars[,to.split][[j]], res, weights)
				p.interact[i,j] <- tmp.interact
				if( !is.null(attr(tmp.interact, "Chisq")) ){ X2.interact[i,j] <- attr(tmp.interact, "Chisq")}
			}
		}
    }

    names(p.curv) <- names(xvars[,to.split])
   	attr(p.curv, "Chisq") <- X2.curv	#either NA or value of chisq-statistic
	dimnames(p.interact) <- list( names(xvars[,to.split]), names(xvars[,to.split]))	
    attr(p.interact, "Chisq") <- X2.interact	
    # X2.interact and p.interact contain NA's in and below diagonal
	
	if( min(p.curv) <= min(p.interact,na.rm=TRUE) ){
        id <- index.counter[which.min(p.curv)]
		return( list( id = as.integer(id), residuals=res ) )
	} 
	else{ # smallest p value from interaction test
		wm <- which.min(p.interact)
		# extracting variables whose interaction is most significant
		cn <- ceiling( wm/split.num )   # column number of interacting pair in p.interact
		rn <- wm %% split.num           # row number of interacting pair in p.interact
		# wm %% split.num = 0 cannot occur because of NA's in last row of p.interact

        idcn <- index.counter[cn]    # index of first variable in xvars
        idrn <- index.counter[rn]    # index of second variable in xvars

		x1 <- xvars[[idcn]]
		x2 <- xvars[[idrn]]

		if( all(vartype[c(idcn,idrn)] %in% "n") ){
		# node is split in turn along sample mean of each variable;
		# for each split, the SSE for a fitted model is obtained for each subnode; 
		# the variable yielding the split with the smaller total SSE is selected
			m1 <- mean( x1[weights > 0], na.rm=TRUE )
			m2 <- mean( x2[weights > 0], na.rm=TRUE )
			# Split 1
			x1m1 <- partysplit(varid = as.integer(idcn), breaks = m1)
			kidids1 <- kidids_split( x1m1, data = xvars)
            # Split 2
			x2m2 <- partysplit(varid = as.integer(idrn), breaks = m2)
			kidids2 <- kidids_split( x2m2, data=xvars)

            log.ids11 <- (weights > 0) & (kidids1 == 1)
            log.ids12 <- (weights > 0) & (kidids1 == 2)
            log.ids21 <- (weights > 0) & (kidids2 == 1)
            log.ids22 <- (weights > 0) & (kidids2 == 2)

            if( sum(log.ids11) < 1 || sum(log.ids12) < 1 || 
                sum(log.ids21) < 1 || sum(log.ids22) < 1) {
                # no model to fit in respective childnode because of lack of data,
                # select variable with smaller curvature p-value 
                if( p.curv[cn] <= p.curv[rn]){ 
                    return( list( id = as.integer(idcn), residuals=res) )}
                else{ return( list( id = as.integer(idrn), residuals=res) )}
            }                   

            switch( complexity, 
                "const" = {
                    res.x1m1.1 <- lm( response[log.ids11] ~ 1)$residuals
                    res.x1m1.2 <- lm( response[log.ids12] ~ 1)$residuals
                    res.x2m2.1 <- lm( response[log.ids21] ~ 1)$residuals
                    res.x2m2.2 <- lm( response[log.ids22] ~ 1)$residuals
                },
                "mult" = {
                    res.x1m1.1 <- lm( df.nf[log.ids11, ] )$residuals
                    res.x1m1.2 <- lm( df.nf[log.ids12, ] )$residuals
                    res.x2m2.1 <- lm( df.nf[log.ids21, ] )$residuals
                    res.x2m2.2 <- lm( df.nf[log.ids22, ] )$residuals
                },
                "simple" = {                    
                    res.xvars11 <- apply( xvars.to.fit[log.ids11,], 2,
                        function(x){ lm.fit( cbind(rep.int(1,sum(log.ids11)),x),
                        response[log.ids11])$residuals} )           
                    if(!is.null(dim(res.xvars11)) && length(dim(res.xvars11))>1){
                        res.x1m1.1 <- res.xvars11[,which.min(colSums(res.xvars11^2))]
                    } else { 
                        res.x1m1.1 <- res.xvars11[which.min(res.xvars11^2)]
                    }

                    res.xvars12 <- apply( xvars.to.fit[log.ids12,], 2,
                        function(x){ lm.fit( cbind(rep.int(1,sum(log.ids12)),x),
                        response[log.ids12])$residuals} )
                    if(!is.null(dim(res.xvars12)) && length(dim(res.xvars12))>1){
                        res.x1m1.2 <- res.xvars12[,which.min(colSums(res.xvars12^2))]
                    } else { 
                        res.x1m1.2 <- res.xvars12[which.min(res.xvars12^2)]
                    }

                    res.xvars21 <- apply( xvars.to.fit[log.ids21,], 2,
                        function(x){ lm.fit( cbind(rep.int(1,sum(log.ids21)),x),
                        response[log.ids21])$residuals} )
                    if(!is.null(dim(res.xvars21)) && length(dim(res.xvars21))>1){
                        res.x2m2.1 <- res.xvars21[,which.min(colSums(res.xvars21^2))]
                    } else { 
                        res.x2m2.1 <- res.xvars21[which.min(res.xvars21^2)]
                    }

                    res.xvars22 <- apply( xvars.to.fit[log.ids22,], 2,
                        function(x){ lm.fit( cbind(rep.int(1,sum(log.ids22)),x),
                        response[log.ids22])$residuals} )
                    if(!is.null(dim(res.xvars22)) && length(dim(res.xvars22))>1){
                        res.x2m2.2 <- res.xvars22[,which.min(colSums(res.xvars22^2))]
                    } else { 
                        res.x2m2.2 <- res.xvars22[which.min(res.xvars22^2)]
                    }
                }
			)
            s1 <- sum(res.x1m1.1^2) + sum(res.x1m1.2^2) 
			s2 <- sum(res.x2m2.1^2) + sum(res.x2m2.2^2)
			if(s1 <= s2){ return( list( id= as.integer(idcn), residuals=res))}
			else{ return( list( id= as.integer(idrn), residuals=res))}
		} 
        else if(complexity=="const"){
            # select variable with smaller curvature p-value 
            if( p.curv[cn] <= p.curv[rn]){ return( list( id = as.integer(idcn), residuals=res) )}
            else{ return( list( id = as.integer(idrn), residuals=res) )}
        } 
        else {#i.e. complexity is "mult" or "simple", not both variables of type "n"   
            if( all( !(vartype[c(idcn,idrn)] %in% "n") )){
                # select variable with smaller curvature p-value 
                if( p.curv[cn] <= p.curv[rn]){ return( list( id = as.integer(idcn), residuals=res) )}
                else{ return( list( id = as.integer(idrn), residuals=res) )}
            } else { # return variable that is not of vartype "n"
                return( list(id = as.integer( c(idcn,idrn)[!(vartype[c(idcn,idrn)] %in% "n")]) ,
                        residuals=res) )
            }
        }
	}
}


curv.test <- function(x, resid, weights){
	if( is.numeric(x) ){ # Omitting values of x where weights=0 for grouping according to quantiles
        q <- quantile( x[weights > 0], names=FALSE, na.rm=TRUE)
        i <- 4
        while( any(diff(q)==0) && (i > 2) ){
            i <- i-1
            q <- quantile( x[weights > 0], probs = seq(0, 1, 1/i), names=FALSE, na.rm=TRUE)
        }
        if( i==2 && any(diff(q)==0) ){ #all data belong to the same class, here labeled 1
            x <- factor(rep.int(1,length(x)))
        }
        else {
		    x <- cut( x, breaks= q, include.lowest=TRUE) 
        }
    } 
	else{ # x is factor 
		x <- x[ , drop=TRUE]
    }
	xytabs <- xtabs( ~ (resid >= 0) + x[weights>0])
    # chi-square test on xytabs 
	ChiSquare(xytabs)
}


interact.test <- function(x1, x2, resid, weights){
	# if x1, x2 are both numeric:
	if( is.numeric(x1) && is.numeric(x2)){
		m1 <- median( x1[weights > 0],na.rm=TRUE)
		m2 <- median( x2[weights > 0],na.rm=TRUE)	
		x <- factor( 1*(x1 > m1)+ 2*(x2 > m2))
		xytabs <- xtabs( ~ (resid >= 0) + x[weights>0])
		ChiSquare(xytabs)
	}
	else if( is.factor(x1) && is.factor(x2)){
		# if x1, x2 are both factors:
		x <- factor( paste( x1, x2, sep=""))
		xytabs <- xtabs( ~ (resid >= 0) + x[weights>0])
		ChiSquare(xytabs)
	} 
	else if( is.numeric(x1)){  # if exactly one of x1, x2 is numeric (and the other factor)
		m1 <- median( x1[weights > 0],na.rm=TRUE)
		x <- factor( paste( (x1 > m1), x2, sep=""))
		xytabs <- xtabs( ~ (resid >= 0) + x[weights>0])
		ChiSquare(xytabs)
	}
	else{
		m2 <- median( x2[weights > 0],na.rm=TRUE)
		x <- factor( paste( (x2 > m2), x1, sep=""))
		xytabs <- xtabs( ~ (resid >= 0) + x[weights>0])
		ChiSquare(xytabs)
	}
}


splitpoint <- function(response, xvars, resid, weights, splitid, splitmethod="G", complexity, vartype){
	x <- xvars[[splitid]][weights > 0]
	if( is.numeric(x)){
		if( splitmethod=="G") {
			# G-method: split points by greedy search
			uni <- sort( unique( x)) # similar values would lead to similar splits 
			sp <- uni[-length(uni)] + diff(uni)/2
			a <- numeric( length(sp))
			for( i in 1:length(sp)){	
				xsp <- partysplit(varid = as.integer(splitid), breaks = sp[i])
				kidids <- kidids_split( xsp, data = xvars)
                log.ids1 <- (weights > 0) & (kidids == 1)                        
                log.ids2 <- (weights > 0) & (kidids == 2)                        

                switch( complexity,
                    "const" = {
                        res.xsp1 <- lm( response[log.ids1] ~ 1)$residuals
                        res.xsp2 <- lm( response[log.ids2] ~ 1)$residuals
				    },
                    "mult" = {
                        df.nf <- data.frame(response, xvars[,vartype %in% c("n","f")])
                        res.xsp1 <- lm( df.nf[log.ids1, ] )$residuals
                        res.xsp2 <- lm( df.nf[log.ids2, ] )$residuals
                    },
                    "simple" = {
                        xvars.to.fit <- xvars[, vartype %in% c("n","f")]
                        res.xvars1 <- sapply( xvars.to.fit[log.ids1,],
                            function(x){ lm.fit( cbind(rep.int(1,sum(log.ids1)),x),
                            response[log.ids1])$residuals} )
                        if(!is.null(dim(res.xvars1)) && length(dim(res.xvars1))>1){
                            res.xsp1 <- res.xvars1[,which.min(colSums(res.xvars1^2))]
                        } else { 
                            res.xsp1 <- res.xvars1[which.min(res.xvars1^2)]
                        }
                        res.xvars2 <- sapply( xvars.to.fit[log.ids2,],
                            function(x){ lm.fit( cbind(rep.int(1,sum(log.ids2)),x),
                            response[log.ids2])$residuals} )
                        if(!is.null(dim(res.xvars2)) && length(dim(res.xvars2))>1){
                            res.xsp2 <- res.xvars2[,which.min(colSums(res.xvars2^2))]
                        } else {
                            res.xsp2 <- res.xvars2[which.min(res.xvars2^2)]
                        }                    
                    }
                )
				a[i] <- sum(res.xsp1^2) + sum(res.xsp2^2)
			}
			return( list( breaks=sp[ which.min(a)], index=NULL) )
		} 
		if(splitmethod=="M") {		
			# M-method: split point is sample median
			return( list( breaks=median( x,na.rm=TRUE), index=NULL) )
		} else {
			stop("unknown splitmethod")
		}
	} 
	else{# i.e. x is factor with more than one level!
        class12 <-  resid >= 0 
	    A <- levels( x[,drop=TRUE])[sort.list( by( class12, x[, drop=TRUE], mean) )]	
	    nm <- length(A)-1
	    binom.var <- numeric(nm)
	    for(i in 1:nm){
		    nodeL <- class12[ x %in% A[1:i]]
		    varnodeL <- sum( nodeL) *(length(nodeL)-sum(nodeL))/length(nodeL)^2
		    nodeR <- class12[ !(x %in% A[1:i])]
		    varnodeR <- sum( nodeR) *(length(nodeR)-sum(nodeR))/length(nodeR)^2
		    binom.var[i] <- varnodeL + varnodeR
        }
	    index <- (!(levels(x) %in% A[1:which.min( binom.var)])) + 1L 
	    return( list( breaks= NULL, index=index) )
    }
}


guide_intern <- function(id = 1L, response, x, weights = NULL, complexity, vartype, ctrl = guide_control()) {

    if (is.null(weights)) weights <- rep.int(1, length(response))
    if (sum(weights) < ctrl$minsplit) return(partynode(id = id)) 

	splitvar <- splitvar(response, x, weights, complexity, vartype)

    if( is.factor(x[[splitvar$id]]) && nlevels( x[[splitvar$id]][ weights > 0, drop=TRUE]) < 2) 
        return(partynode(id = id))
    if( is.numeric(x[[splitvar$id]]) && length(unique(x[[splitvar$id]][weights > 0])) < 2)
        return(partynode(id = id))

	splitpoint <- splitpoint(response, x, splitvar$residuals, weights, splitvar$id, 
					splitmethod = ctrl$splitmethod, complexity, vartype)

	sp <- partysplit(varid = splitvar$id, breaks = splitpoint$breaks, index = splitpoint$index)

    kidids <- kidids_split(sp, data = x)
	ln <- max(unique(kidids[weights > 0]))
   	kids <- vector(mode = "list", length = ln  )
    	for (kidid in 1:ln) {
          w <- weights
          w[kidids != kidid] <- 0
          if (kidid > 1) {
              myid <- max(nodeids(kids[[kidid - 1]]))
          } else {
              myid <- id
          }
          kids[[kidid]] <- guide_intern(id = as.integer(myid + 1), response, x, 
                                 weights = w, complexity, vartype, ctrl = guide_control())
    }
    return(partynode(id = as.integer(id), split = sp, kids = kids))
}


guide <- function(formula, data, vartype = NULL, complexity=c("const","mult","simple"), subset, 
            weights, na.action = na.omit, control = guide_control())
{	
    complexity <- match.arg(complexity)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- FALSE
    mf[[1L]] <- as.name("model.frame") 
    m <- eval.parent(mf)
    y <- model.response(m)
    if(any(attr(attr(m, "terms"), "order") > 1L))
	    stop("No interaction terms in formula allowed")

	if(!is.numeric(y))
		stop("Response must be numeric for regression tree")
	x <- m[, c(-1, -which(names(m) == "(weights)")), drop = FALSE]

	categorical <- sapply(x, function(x) class(x)[1]) %in% c("factor", "ordered")
        
	if (!is.null(vartype)) {
		stopifnot( length(vartype) == length(x))
        if (length(vartype[categorical])>0)
		  stopifnot( all(vartype[categorical] == "c"))
        if (length(vartype[!categorical])>0){
            if(complexity == "const"){ stopifnot( all(vartype[!categorical] == "n" ) )}
            else{ stopifnot( all(vartype[!categorical] %in% c("n","f","s") ))
                stopifnot( any(vartype %in% c("n","f"))) 
            }
        }          
	} else {
		vartype <- character(length(x))
		vartype[categorical] <- "c"
		vartype[!categorical] <- "n"
	}

    w <- model.weights(m)
	guidetree <- guide_intern(1L, y, x, weights = w, complexity, vartype, ctrl = control)
	tree <- party(guidetree, data = x, 
                    fitted = data.frame("(fitted)" = fitted_node(guidetree, data = x),
                                        "(response)" = y, check.names = FALSE),
                    terms = terms(formula, data = data))
    if (!missing(weights))
        tree$fitted[["(weights)"]] <- w
    class(tree) <- c("constparty", "party")
    tree
}


guide_control <- function( splitmethod = "G", minsplit = 20) {
    ret <- list( splitmethod = splitmethod, minsplit = minsplit)
    class(ret) <- "guide_control"
    return(ret)
}






