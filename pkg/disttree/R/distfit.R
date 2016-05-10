distfit <- function(y, family, weights = NULL, start = NULL, vcov. = TRUE, estfun = TRUE, ...)
{
  ## match call
  cl <- match.call()
  
  ## number of distribution parameters (mu, sigma, nu, tau)
  if(is.function(family)) family <- family()
  np <- sum(family$parameter == TRUE)
  
  ## number of observations
  ny <- length(y)
  
  ## weights
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1L) weights <- rep.int(weights, ny)
  weights <- as.vector(weights)
  
  ## store y and select observations with weight > 0 
  #FIXME# y.store <- y          
  #FIXME# y <- y[weights > 0]
  
  ## number of observations = sum of weights (i.e., case weights)
  ## FIXME ## also need proportionality weights, i.e., weights = sum(weights > 0) ?
  nobs <- sum(weights)
  
  ## notation:
  # par ... distribution parameters (mu, sigma, nu, tau)
  # eta ... coefficients of the linear predictor, here: intercept (g(mu)=eta[1], g(sigma)=eta[2], g(nu)=eta[3], g(tau)=eta[4])
  
  # if(np > 0L) m <- family$mu.linkinv(eta[1L])          # m ... mu           eta[1] ... g(mu)        g ... link function
  # if(np > 1L) s <- family$sigma.linkinv(eta[2L])       # s ... sigma        eta[2] ... g(sigma)     g ... link function
  # if(np > 2L) v <- family$nu.linkinv(eta[3L])          # v ... nu           eta[3] ... g(nu)        g ... link function
  # if(np > 3L) t <- family$tau.linkinv(eta[4L])         # t ... tau          eta[4] ... g(tau)       g ... link function
  
  
  
  ## Define all necessary functions depending on the number of parameters
  
  # get parameters of a function f, return vector with the indices of the necessary input parameters 
  getpar <- function(f){
    arguments <- names(formals(f))
    
    ## 2 cases: with or without y as input
    
    ## with y
    # if y is one of the input arguments, the derivative function returns a vector of the length of y
    
    if(arguments[1]=="y"){
      # 4 Parameter
      if(length(arguments)==5L) par.id <- c(0,1,2,3,4)        # f(y, mu=m, sigma=s, nu=v, tau=t)
      
      # 3 Parameter
      if(length(arguments)==4L){
        if(arguments[2]=="mu"){
          if(arguments[3]=="sigma"){
            if(arguments[4]=="nu"){
              par.id <- c(0,1,2,3)                            # f(y, mu=m, sigma=s, nu=v)
            } else {
              par.id <- c(0,1,2,4)                            # f(y, mu=m, sigma=s, tau=t)
            }
          } else {
            par.id <- c(0,1,3,4)                              # f(y, mu=m, nu=v, tau=t)
          }
        } else {
          par.id <- c(0,2,3,4)                                # f(y, sigma=s, nu=v, tau=t)
        }
      }
      
      # 2 Parameter
      if(length(arguments)==3L){
        if(arguments[2]=="mu"){
          if(arguments[3]=="sigma") par.id <- c(0,1,2)        # f(y, mu=m, sigma=s)
          if(arguments[3]=="nu")    par.id <- c(0,1,3)        # f(y, mu=m, nu=v)
          if(arguments[3]=="tau")   par.id <- c(0,1,4)        # f(y, mu=m, tau=t)
        }
        if(arguments[2]=="sigma"){
          if(arguments[3]=="nu")    par.id <- c(0,2,3)        # f(y, sigma=s, nu=v)
          if(arguments[3]=="tau")   par.id <- c(0,2,4)        # f(y, sigma=s, tau=t)
        }
        if(arguments[2]=="nu")      par.id <- c(0,3,4)        # f(y, nu=v, tau=t)
      }
      
      # 1 Parameter
      if(length(arguments)==2L){
        if(arguments[2]=="mu")      par.id <- c(0,1)          # f(y, mu=m)
        if(arguments[2]=="sigma")   par.id <- c(0,2)          # f(y, sigma=s)
        if(arguments[2]=="nu")      par.id <- c(0,3)          # f(y, nu=v)
        if(arguments[2]=="tau")     par.id <- c(0,4)          # f(y, tau=t)
      }
      
      # 0 Parameter
      if(length(arguments)==1L)     par.id <- c(0)            # f(y)
    } else {
      
      ## without y
      # in this case the functions return only a single value -> create vector by replicating this value ny times (necessary for the matrix (using sum and *1/ny)) 
      
      # 4 Parameter
      if(length(arguments)==4L) par.id <- c(1,2,3,4)          # f(mu=m, sigma=s, nu=v, tau=t)
      
      # 3 Parameter
      if(length(arguments)==3L){
        if(arguments[1]=="mu"){
          if(arguments[2]=="sigma"){
            if(arguments[3]=="nu"){
              par.id <- c(1,2,3)                              # f(mu=m, sigma=s, nu=v)
            } else {
              par.id <- c(1,2,4)                              # f(mu=m, sigma=s, tau=t)
            }
          } else {
            par.id <- c(1,3,4)                                # f(mu=m, nu=v, tau=t)
          }
        } else {
          par.id <- c(2,3,4)                                  # f(sigma=s, nu=v, tau=t)
        }
      }
      
      # 2 Parameter
      if(length(arguments)==2L){
        if(arguments[1]=="mu"){
          if(arguments[2]=="sigma") par.id <- c(1,2)          # f(mu=m, sigma=s)
          if(arguments[2]=="nu")    par.id <- c(1,3)          # f(mu=m, nu=v)
          if(arguments[2]=="tau")   par.id <- c(1,4)          # f(mu=m, tau=t)
        }
        if(arguments[1]=="sigma"){ 
          if(arguments[2]=="nu")    par.id <- c(2,3)          # f(sigma=s, nu=v)
          if(arguments[2]=="tau")   par.id <- c(2,4)          # f(sigma=s, tau=t)
        }
        if(arguments[1]=="nu")      par.id <- c(3,4)          # f(nu=v, tau=t)
      }
      
      # 1 Parameter
      if(length(arguments)==1L){
        if(arguments[1]=="mu")    par.id <- c(1)              # f(mu=m)
        if(arguments[1]=="sigma") par.id <- c(2)              # f(sigma=s)
        if(arguments[1]=="nu")    par.id <- c(3)              # f(nu=v)
        if(arguments[1]=="tau")   par.id <- c(4)              # f(tau=t)
      }
      
      # 0 Parameter
      if(length(arguments)==0L)  par.id <- NULL       ## fix: possible case? of which class is f when the derivative is a constant?
    }
    
    return(par.id)
  }
  
  
  
  ## define inner and outer derivative functions
  
  if(np > 0L){
    
    # inner derivative functions (dmdeta, d2mdeta2)
    dmdeta <- function(eta) return(family$mu.dr(eta[1]))
    if(family$mu.link=="identity") d2mdeta2 <- function(eta) return(0)
    if(family$mu.link=="log")      d2mdeta2 <- function(eta) return(exp(eta[1]))
    if(family$mu.link=="logit")    d2mdeta2 <- function(eta) return(exp(eta[1]) * (exp(eta[1])-1) / ((1+exp(eta[1]))^3)) 
    
    # outer derivative functions (dldm, d2ldm2)
    par.id.dldm <- getpar(family$dldm)
    if(par.id.dldm[1] == 0L) {
      if(length(par.id.dldm) == 1L){
        dldm <- function(par) return(family$dldm(y))
      } else {
        dldm <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.dldm)) input[[i]] <- rep.int(par[par.id.dldm[i]], length(y))
          return(do.call(family$dldm, input))
        }
      }
    } else {
      dldm <- function(par) {
        input <- list()
        for (i in 1:length(par.id.dldm)) input[[i]] <- par[par.id.dldm[i]]
        return(rep.int(do.call(family$dldm, input), length(y)))
      }
    }
    
    par.id.d2ldm2 <- getpar(family$d2ldm2)
    if(par.id.d2ldm2[1] == 0L){
      if(length(par.id.d2ldm2) == 1L){
        d2ldm2 <- function(par) return(family$d2ldm2(y))
      } else {
        d2ldm2 <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldm2)) input[[i]] <- rep.int(par[par.id.d2ldm2[i]], length(y))
          return(do.call(family$d2ldm2, input))
        }
      }
    } else {
      d2ldm2 <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldm2)) input[[i]] <- par[par.id.d2ldm2[i]]
        return(rep.int(do.call(family$d2ldm2, input), length(y)))
      }
    }
  }
  
  
  if(np > 1L){
    
    # inner derivative functions (dddeta, d2ddeta2)     
    dddeta <- function(eta) return(family$sigma.dr(eta[2]))
    if(family$sigma.link=="identity") d2ddeta2 <- function(eta) return(0)
    if(family$sigma.link=="log")      d2ddeta2 <- function(eta) return(exp(eta[2]))
    if(family$sigma.link=="logit")    d2ddeta2 <- function(eta) return(exp(eta[2]) * (exp(eta[2])-1) / ((1+exp(eta[2]))^3)) 
    
    # outer derivative functions (dldd, d2ldd2, d2ldmdd)
    par.id.dldd <- getpar(family$dldd)
    if(par.id.dldd[1] == 0L){
      if(length(par.id.dldd) == 1L){
        dldd <- function(par) return(family$dldd(y))
      } else {
        dldd <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.dldd)) input[[i]] <- rep.int(par[par.id.dldd[i]], length(y))
          return(do.call(family$dldd, input))
        }
      }
    } else {
      dldd <- function(par) {
        input <- list()
        for (i in 1:length(par.id.dldd)) input[[i]] <- par[par.id.dldd[i]]
        return(rep.int(do.call(family$dldd, input), length(y)))
      }
    }
    
    par.id.d2ldd2 <- getpar(family$d2ldd2)
    if(par.id.d2ldd2[1] == 0L){
      if(length(par.id.d2ldd2) == 1L){
        d2ldd2 <- function(par) return(family$d2ldd2(y))
      } else {
        d2ldd2 <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldd2)) input[[i]] <- rep.int(par[par.id.d2ldd2[i]], length(y))
          return(do.call(family$d2ldd2, input))
        }
      }
    } else {
      d2ldd2 <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldd2)) input[[i]] <- par[par.id.d2ldd2[i]]
        return(rep.int(do.call(family$d2ldd2, input), length(y)))
      }
    }
    
    par.id.d2ldmdd <- getpar(family$d2ldmdd)
    if(par.id.d2ldmdd[1] == 0L){
      if(length(par.id.d2ldmdd) == 1L){
        d2ldmdd <- function(par) return(family$d2ldmdd(y))
      } else {
        d2ldmdd <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldmdd)) input[[i]] <- rep.int(par[par.id.d2ldmdd[i]], length(y))
          return(do.call(family$d2ldmdd, input))
        }
      }
    } else {
      d2ldmdd <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldmdd)) input[[i]] <- par[par.id.d2ldmdd[i]]
        return(rep.int(do.call(family$d2ldmdd, input), length(y)))
      }
    }
  }
  
  
  if(np > 2L){
    
    # inner derivative functions (dvdeta, d2vdeta2)
    dvdeta <- function(eta) return(family$nu.dr(eta[3]))
    if(family$nu.link=="identity") d2vdeta2 <- function(eta) return(0)
    if(family$nu.link=="log")      d2vdeta2 <- function(eta) return(exp(eta[3]))
    if(family$nu.link=="logit")    d2vdeta2 <- function(eta) return(exp(eta[3]) * (exp(eta[3])-1) / ((1+exp(eta[3]))^3)) 
    
    # outer derivatives (dldv, d2ldv2, d2ldmdv, d2ldddv)
    par.id.dldv <- getpar(family$dldv)
    if(par.id.dldv[1] == 0L){
      if(length(par.id.dldv) == 1L){
        dldv <- function(par) return(family$dldv(y))
      } else {
        dldv <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.dldv)) input[[i]] <- rep.int(par[par.id.dldv[i]], length(y))
          return(do.call(family$dldv, input))
        }
      }
    } else {
      dldv <- function(par) {
        input <- list()
        for (i in 1:length(par.id.dldv)) input[[i]] <- par[par.id.dldv[i]]
        return(rep.int(do.call(family$dldv, input), length(y)))
      }
    }
    
    par.id.d2ldv2 <- getpar(family$d2ldv2)
    if(par.id.d2ldv2[1] == 0L){
      if(length(par.id.d2ldv2) == 1L){
        d2ldv2 <- function(par) return(family$d2ldv2(y))
      } else {
        d2ldv2 <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldv2)) input[[i]] <- rep.int(par[par.id.d2ldv2[i]], length(y))
          return(do.call(family$d2ldv2, input))
        }
      }
    } else {
      d2ldv2 <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldv2)) input[[i]] <- par[par.id.d2ldv2[i]]
        return(rep.int(do.call(family$d2ldv2, input), length(y)))
      }
    }
    
    par.id.d2ldmdv <- getpar(family$d2ldmdv)
    if(par.id.d2ldmdv[1] == 0L) { 
      if(length(par.id.d2ldmdv) == 1L){
        d2ldmdv <- function(par) return(family$d2ldmdv(y))
      } else {
        d2ldmdv <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldmdv)) input[[i]] <- rep.int(par[par.id.d2ldmdv[i]], length(y))
          return(do.call(family$d2ldmdv, input))
        }
      }
    } else {
      d2ldmdv <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldmdv)) input[[i]] <- par[par.id.d2ldmdv[i]]
        return(rep.int(do.call(family$d2ldmdv, input), length(y)))
      }
    }
    
    par.id.d2ldddv <- getpar(family$d2ldddv)
    if(par.id.d2ldddv[1] == 0L) { 
      if(length(par.id.d2ldddv) == 1L){
        d2ldddv <- function(par) return(family$d2ldddv(y))
      } else {
        d2ldddv <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldddv)) input[[i]] <- rep.int(par[par.id.d2ldddv[i]], length(y))
          return(do.call(family$d2ldddv, input))
        }
      }
    } else {
      d2ldddv <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldddv)) input[[i]] <- par[par.id.d2ldddv[i]]
        return(rep.int(do.call(family$d2ldddv, input), length(y)))
      }
    }
  }
  
  
  if(np > 3L){
    
    # inner derivatives (dtdeta, d2tdeta2)    
    dtdeta <- function(eta) return(family$tau.dr(eta[4]))
    if(family$tau.link=="identity")  d2tdeta2 <- function(eta) return(0)
    if(family$tau.link=="log")       d2tdeta2 <- function(eta) return(exp(eta[4]))
    if(family$tau.link=="logit")     d2tdeta2 <- function(eta) return(exp(eta[4]) * (exp(eta[4])-1) / ((1+exp(eta[4]))^3)) 
    
    # outer derivatives (dldt, d2ldt2, d2ldmdt, d2ldddt, d2ldvdt)
    par.id.dldt <- getpar(family$dldt)
    if(par.id.dldt[1] == 0L){
      if(length(par.id.dldt) == 1L){
        dldt <- function(par) return(family$dldt(y))
      } else {
        dldt <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.dldt)) input[[i]] <- rep.int(par[par.id.dldt[i]], length(y))
          return(do.call(family$dldt, input))
        }
      }
    } else {
      dldt <- function(par) {
        input <- list()
        for (i in 1:length(par.id.dldt)) input[[i]] <- par[par.id.dldt[i]]
        return(rep.int(do.call(family$dldt, input), length(y)))
      }
    }
    
    par.id.d2ldt2 <- getpar(family$d2ldt2)
    if(par.id.d2ldt2[1] == 0L){
      if(length(par.id.d2ldt2) == 1L){
        d2ldt2 <- function(par) return(family$d2ldt2(y))
      } else {
        d2ldt2 <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldt2)) input[[i]] <- rep.int(par[par.id.d2ldt2[i]], length(y))
          return(do.call(family$d2ldt2, input))
        }
      }
    } else {
      d2ldt2 <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldt2)) input[[i]] <- par[par.id.d2ldt2[i]]
        return(rep.int(do.call(family$d2ldt2, input), length(y)))
      }
    }
    
    par.id.d2ldmdt <- getpar(family$d2ldmdt)
    if(par.id.d2ldmdt[1] == 0L) { 
      if(length(par.id.d2ldmdt) == 1L){
        d2ldmdt <- function(par) return(family$d2ldmdt(y))
      } else {
        d2ldmdt <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldmdt)) input[[i]] <- rep.int(par[par.id.d2ldmdt[i]], length(y))
          return(do.call(family$d2ldmdt, input))
        }
      }
    } else {
      d2ldmdt <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldmdt)) input[[i]] <- par[par.id.d2ldmdt[i]]
        return(rep.int(do.call(family$d2ldmdt, input), length(y)))
      }
    } 
    
    par.id.d2ldddt <- getpar(family$d2ldddt)
    if(par.id.d2ldddt[1] == 0L) { 
      if(length(par.id.d2ldddt) == 1L){
        d2ldddt <- function(par) return(family$d2ldddt(y))
      } else {
        d2ldddt <- function(par) {
          d2ldddt <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldddt)) input[[i]] <- rep.int(par[par.id.d2ldddt[i]], length(y))
          return(do.call(family$d2ldddt, input))
        }
      }
    } else {
      d2ldddt <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldddt)) input[[i]] <- par[par.id.d2ldddt[i]]
        return(rep.int(do.call(family$d2ldddt, input), length(y)))
      }
    }
    
    par.id.d2ldvdt <- getpar(family$d2ldvdt)
    if(par.id.d2ldvdt[1] == 0L) { 
      if(length(par.id.d2ldvdt) == 1L){
        d2ldvdt <- function(par) return(family$d2ldvdt(y))
      } else {
        d2ldvdt <- function(par) {
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldvdt)) input[[i]] <- rep.int(par[par.id.d2ldvdt[i]], length(y))
          return(do.call(family$d2ldvdt, input))
        }
      }
    } else {
      d2ldvdt <- function(par) {
        input <- list()
        for (i in 1:length(par.id.d2ldvdt)) input[[i]] <- par[par.id.d2ldvdt[i]]
        return(rep.int(do.call(family$d2ldvdt, input), length(y)))
      }
    }
  }
  
  
  
  ## define complete derivative functions dpardeta, d2pardeta2, dldpar, d2ldpar2 according to the number of parameters
  
  if(np == 1L){
    
    # define function for the calculation of initial values
    ## FIXME ## use weights?
    initialize <- function(y) {
      mu <- NULL
      eval(family$mu.initial)
      family$mu.linkfun(mean(mu))
    }

    
    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]))
      names(par) <- c("mu")
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / lsit of matrices:
    dldpar <- function(par){
      dmatrix <- cbind(dldm(par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(par){
      
      d2matrix <- rbind(cbind(d2ldm2(par)))
      
      # d2matrix is of size (1*ny x 1) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (1x1) is stored in d2list
      
      d2list <- list()
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i),]
      }
      
      return(d2list)
    }
  }
  
  
  if(np == 2L){
    
    # define function for the calculation of initial values
    initialize <- function(y) {
      mu <- sigma <- NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)))
    }
    
    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]))
      names(par) <- c("mu","sigma")
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list of matrices:
    dldpar <- function(par){
      dmatrix <- cbind(dldm(par), dldd(par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(par){
      
      d2matrix <- rbind(cbind(d2ldm2(par), d2ldmdd(par)),
                        cbind(d2ldmdd(par), d2ldd2(par)))
      
      # d2matrix is of size (2*ny x 2) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (2x2) is stored in d2list
      
      d2list <- list()
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i),]
        }

      return(d2list)
    }
  }
  
  
  if(np == 3L){
    
    # define function for the calculation of initial values
    initialize <- function(y) {
      mu <- sigma <- nu <-  NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      eval(family$nu.initial)
      c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)))
    }
    
    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]))
      names(par) <- c("mu","sigma","nu")
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta), dvdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta), d2vdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list:
    dldpar <- function(par){
      dmatrix <- cbind(dldm(par), dldd(par), dldv(par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(par){
      
      d2matrix <- rbind(cbind(d2ldm2(par), d2ldmdd(par), d2ldmdv(par)),
                        cbind(d2ldmdd(par), d2ldd2(par), d2ldddv(par)),
                        cbind(d2ldmdv(par), d2ldddv(par), d2ldv2(par)))
      
      # d2matrix is of size (3*ny x 3) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (3x3) is stored in d2list
      
      d2list <- list()
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i),]
      }

      return(d2list)
    }
  }
  
  
  if(np == 4L){
    
    # define function for the calculation of initial values
    initialize <- function(y) {
      mu <- sigma <- nu <- tau <- NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      eval(family$nu.initial)
      eval(family$tau.initial)
      c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)), family$tau.linkfun(mean(tau)))
    }

    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]), family$tau.linkinv(eta[4]))
      names(par) <- c("mu","sigma","nu","tau")
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta), dvdeta(eta), dtdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta), d2vdeta2(eta), d2tdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list :
    dldpar <- function(par){
      dmatrix <- cbind(dldm(par), dldd(par), dldv(par), dldt(par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(par){
      
      d2matrix <- rbind(cbind(d2ldm2(par), d2ldmdd(par), d2ldmdv(par), d2ldmdt(par)),
                        cbind(d2ldmdd(par), d2ldd2(par), d2ldddv(par), d2ldddt(par)),
                        cbind(d2ldmdv(par), d2ldddv(par), d2ldv2(par), d2ldvdt(par)),
                        cbind(d2ldmdt(par), d2ldddt(par), d2ldvdt(par), d2ldt2(par)))
      
      # d2matrix is of size (4*ny x 4) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (4x4) is stored in d2list
      
      d2list <- list()
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i, 3*ny+i),]
      }
      
      return(d2list)
    }
  }
    
  

  # starting values for the optimization: parameters of the linear predictor (here: no covariables -> linear predictor = const = intercept) 
  # -> link functions evaluated at the starting values of the distribution parameters
  
  # if start as input: initial value for distribution parameter -> transform adequately to starting values for the parameters of the linear predictor (here: intercepts)
  transstart <- function(startpar){
    starteta <- NA
    if(np > 0) starteta[1] <- family$mu.linkfun(startpar[1])
    if(np > 1) starteta[2] <- family$sigma.linkfun(startpar[2])
    if(np > 2) starteta[3] <- family$nu.linkfun(startpar[3])
    if(np > 3) starteta[4] <- family$tau.linkfun(startpar[4])
    return(starteta)
  }

  
  
  ## set up negative log-likelihood
  nll <- function(eta) {
    # define input list
    par <- distpar(eta)
    input <- list()
    input[[1]] <- y
    for(i in 2:(length(par)+1)) input[[i]] <- par[i-1]         # <- rep.int(par[i-1], length(y)   (FIX?)
    
    # G.dev.incr ... global deviance function = -2*logLik
    nloglik <- do.call(family$G.dev.incr, input)
    nloglik <- sum(weights * nloglik/2)
    return(nloglik)
  }
    
  
  ## set up gradient/scores 
  grad <- function(eta, sum = TRUE) {
    par <- distpar(eta)                                  # get distribution parameters
    gr.out <- dldpar(par)                                # outer derivatives
    gr <- -weights * t(t(gr.out) * dpardeta(eta))        # multiplied with the inner derivatives (componentwise)
    # gr <- as.matrix(gr)                                # for 1 parameter
    if(sum) gr <- colSums(gr)     # *1/nobs ? scale, doesn't influence optimization
    return(gr)
  }
  
  

  ## set up analytical hessian 

  hess <- function(eta){
    
    ## get distribution parameter
    par <- distpar(eta)
    
    ## calculate derivative vectors / matrices / lists
    d2ldpar2.list <- d2ldpar2(par)
    dldpar.mat <- dldpar(par)
    
    dpardeta.vec <- dpardeta(eta)
    d2pardeta2.vec <- d2pardeta2(eta)
    
    ## calculation is split up in 2 parts: 
    # 2nd outer derivatives times first inner derivatives and a diagonal matrix with the first outer and the second inner derivatives
    
    hess <- list()
    length(hess) <- length(d2ldpar2.list)
    for(i in 1:ny){
      hess[[i]] <- weights[i] * (t(d2ldpar2.list[[i]] * dpardeta.vec) * dpardeta.vec + diag(np) * as.vector(dldpar.mat[i,]) * d2pardeta2.vec)
    }
    

    ## calculate the sum over all matrices in the list (each for one observation)  
    
    sumhess <- Reduce('+', hess)
    
    return(-sumhess)   # negative because the optimized function is the neg. logLik (nll) but the derivatives are of the (positive) logLik => *(-1)
    # the entries of sumhes are the sums of the entries of the hessian matrix evaluated at the observations stored in y and the input parameters

  }
  


  ## calculate initial values if necessary or otherwise transform initial values for the distribution parameters to initial values for the intercepts
  if(is.null(start)){
    ## FIXME ## (1) initialization currently doesn't work
    ## FIXME ## (2) initialization with weights? E.g. rep.int(y, round(weights)) ?
    starteta <- initialize(y = rep.int(y, round(weights)))
    startpar <- distpar(starteta)
  } else {
    startpar <- start
    starteta <- transstart(startpar)
  }
  
  
  ## optimize log-likelihood
  opt <- optim(par = starteta, fn = nll, gr = grad, method = "BFGS", hessian = vcov., control = list(...))

  ## extract parameters
  eta <- opt$par
  par <- distpar(eta)


  ## variance-covariance matrix estimate
  if(vcov.) {
    vc <- solve(hess(eta))
    colnames(vc) <- rownames(vc) <- paste0(names(par), ".par")
  } else {
    vc <- NULL
  }
  
  ## estfun
  # each column represents one distribution parameter (1.col -> dldm * dmdpar = "dldmu.par", 2.col -> dldd * dddpar = "dldsigma.par", ...)
  if(estfun) {
    ef <- -grad(eta, sum = FALSE)
    ef <- as.matrix(ef)                                   # FIX: in case ef is a vector (for np=1)
    colnames(ef) <- paste("dld", names(par),".par", sep = "")
  } else {
    ef <- NULL                    
  }

  ## density function
  #ddist <- get(paste("d",family$family[1], sep = ""))
  ddist <- function(x, log = FALSE){
    if(np == 1L) fy <- get(paste("d",family$family[1], sep = ""))(x, mu = par[1], log = FALSE)
    if(np == 2L) fy <- get(paste("d",family$family[1], sep = ""))(x, mu = par[1], sigma = par[2], log = FALSE)
    if(np == 3L) fy <- get(paste("d",family$family[1], sep = ""))(x, mu = par[1], sigma = par[2], nu = par[3], log = FALSE)
    if(np == 4L) fy <- get(paste("d",family$family[1], sep = ""))(x, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], log = FALSE)
    fy
  }
  
  ## cumulative distribution function
  #pdist <- get(paste("p",family$family[1], sep = ""))
  pdist <- function(q, log.p = FALSE){
    if(np == 1L) cdf <- get(paste("p",family$family[1], sep = ""))(q, mu = par[1], log.p = FALSE)
    if(np == 2L) cdf <- get(paste("p",family$family[1], sep = ""))(q, mu = par[1], sigma = par[2], log.p = FALSE)
    if(np == 3L) cdf <- get(paste("p",family$family[1], sep = ""))(q, mu = par[1], sigma = par[2], nu = par[3], log.p = FALSE)
    if(np == 4L) cdf <- get(paste("p",family$family[1], sep = ""))(q, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], log.p = FALSE)
    cdf
  }
  
  ## quantile function
  #qdist <- get(paste("q",family$family[1], sep = ""))
  qdist <- function(p, log.p = FALSE){
    if(np == 1L) q <- get(paste("q",family$family[1], sep = ""))(p, mu = par[1], log.p = FALSE)
    if(np == 2L) q <- get(paste("q",family$family[1], sep = ""))(p, mu = par[1], sigma = par[2], log.p = FALSE)
    if(np == 3L) q <- get(paste("q",family$family[1], sep = ""))(p, mu = par[1], sigma = par[2], nu = par[3], log.p = FALSE)
    if(np == 4L) q <- get(paste("q",family$family[1], sep = ""))(p, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], log.p = FALSE)
    q
  }
  
  ## random function
  #rdist <- get(paste("r",family$family[1], sep = ""))
  rdist <- function(n){
    if(np == 1L) r <- get(paste("r",family$family[1], sep = ""))(n, mu = par[1])
    if(np == 2L) r <- get(paste("r",family$family[1], sep = ""))(n, mu = par[1], sigma = par[2])
    if(np == 3L) r <- get(paste("r",family$family[1], sep = ""))(n, mu = par[1], sigma = par[2], nu = par[3])
    if(np == 4L) r <- get(paste("r",family$family[1], sep = ""))(n, mu = par[1], sigma = par[2], nu = par[3], tau = par[4])
    r
  }
  
  
  ## return value 
  rval <- list(
    y = y,
    weights = weights,
    family = family,
    startpar = startpar,
    starteta = starteta,
    opt = opt,
    par = par,
    eta = eta,
    hess = hess(eta),
    call = cl,
    ny = ny,        
    nobs = nobs,    
    vcov = vc,
    estfun = ef,
    ddist = ddist,
    pdist = pdist,
    qdist = qdist,
    rdist = rdist
  )
  class(rval) <- "distfit"
  return(rval)
}









## print, summary?, predict?
nobs.distfit <- function(object, ...) {
  object$nobs
}

coef.distfit <- function(object, ...) {
  object$opt$par
  # object$eta
}

vcov.distfit <- function(object, ...) {
  object$vcov
}

estfun.distfit <- function(object, ...) {                         
  object$estfun
}

logLik.distfit <- function(object, ...) {
  structure(-object$opt$value, df = length(object$opt$par), class = "logLik")
}

bread.distfit <- function(object, ...) {
  object$vcov * object$nobs
}


confint.distfit <- function(object, ...) {
  vcov <- object$vcov
  eta <- object$eta
  par <- object$par
  np <- length(par)
  if(np > 0){
    mupar.confint <- c(eta[1] + qnorm(0.025) * sqrt(vcov[1,1]), eta[1] + qnorm(0.975) * sqrt(vcov[1,1]))
    confint <- rbind(mupar.confint)
  }
  if(np > 1){
    sigmapar.confint <- c(eta[2] + qnorm(0.025) * sqrt(vcov[2,2]), eta[2] + qnorm(0.975) * sqrt(vcov[2,2]))
    confint <- rbind(confint, sigmapar.confint)
  }
  if(np > 2){
    nupar.confint <- c(eta[3] + qnorm(0.025) * sqrt(vcov[3,3]), eta[3] + qnorm(0.975) * sqrt(vcov[3,3]))
    confint <- rbind(confint, nupar.confint)
  }
  if(np > 3){ 
    taupar.confint <- c(eta[4] + qnorm(0.025) * sqrt(vcov[4,4]), eta[4] + qnorm(0.975) * sqrt(vcov[4,4]))
    confint <- rbind(confint, taupar.confint)
  }
  colnames(confint) <- c("2.5 %", "97.5 %")
  rownames(confint) <- paste(names(par),".par",sep="")
  confint
}


