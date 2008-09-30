pmmlTreeModel <- function(file, ...) {
  stopifnot(require("XML"))
  as.party(xmlRoot(xmlTreeParse(file)))
}

as.party.XMLNode <- function(x) {
  stopifnot(require("XML"))
  ## check whether XML specifies a TreeModel
  stopifnot(c("DataDictionary", "TreeModel") %in% names(x))
  if(any(c("MiningBuildTask", "TransformationDictionary", "Extension") %in% names(x)))
    warning("not yet implemented")
  
  ## needed? x[["Header"]]
  
  ## parse data dictionary
  extract_empty_model_frame <- function(x) {

    ## extract DataDictionary
    dd <- x[["DataDictionary"]]

    ## currently we can only look at DataField  
    if(!all(names(dd) == "DataField")) warning("not yet implemented")
  
    ## check columns
    nc <- as.numeric(xmlAttrs(dd)["numberOfFields"])
    if(!is.na(nc)) stopifnot(nc == length(dd))

    ## set up data frame (only numeric variables)
    mf <- as.data.frame(rep(list(1), nc))[0,]
    names(mf) <- xmlSApply(dd, function(x) xmlAttrs(x)["name"])

    ## modify class if necessary
    for(i in 1:nc) {
      optype <- xmlAttrs(dd[[i]])["optype"]
      switch(optype,
        "categorical" = {
           mf[[i]] <- factor(integer(0),
	     levels = xmlSApply(dd[[i]], function(x) xmlAttrs(x)["value"]))
        },
        "ordinal" = {
           mf[[i]] <- factor(integer(0), ordered = TRUE,
	     levels = xmlSApply(dd[[i]], function(x) xmlAttrs(x)["value"]))
        },
        "continuous" = {
          dataType <- xmlAttrs(dd[[i]])["dataType"]
          if(dataType == "integer") mf[[i]] <- integer(0)
        }
      )
    }
    
    return(mf)
  }
  mf <- extract_empty_model_frame(x)
  mf_names <- names(mf)
  mf_levels <- lapply(mf, levels)

  ## parse MiningSchema
  extract_terms <- function(x) {

    ## extract MiningSchema
    ms <- x[["TreeModel"]]
    stopifnot("MiningSchema" %in% names(ms))
    ms <- ms[["MiningSchema"]]
    
    ## currently we can only look at MiningField  
    if(!all(names(ms) == "MiningField")) warning("not yet implemented")
    
    ## extract variable info
    vars <- t(xmlSApply(ms, xmlAttrs))
    if(sum(vars[,2] == "predicted") > 1) stop("multivariate responses not yet implemented")
    if(!all(vars[,2] %in% c("predicted", "active"))) warning("not yet implemented")
    
    ## set up formula
    ff <- as.formula(paste(vars[vars[,2] == "predicted",1], "~",
      paste(vars[vars[,2] != "predicted",1], collapse = " + ")))

    return(terms(ff))
  }
  trms <- extract_terms(x)
  
  ## parse TreeModel
  tm <- x[["TreeModel"]]
  tm_info <- xmlAttrs(tm)

  ## check response
  stopifnot(tm_info["functionName"] %in% c("classification", "regression"))
  mf_response <- model.response(model.frame(trms, data = mf))
  if(tm_info["functionName"] == "classification") stopifnot(inherits(mf_response, "factor"))
  if(tm_info["functionName"] == "regression") stopifnot(is.numeric(mf_response))
  
  ## convenience functions for parsing nodes
  is_terminal <- function(xnode) !("Node" %in% names(xnode))
  is_root <- function(xnode) "True" %in% names(xnode)
  n_kids <- function(xnode) sum("Node" == names(xnode))
  n_obs <- function(xnode) as.numeric(xmlAttrs(xnode)["recordCount"])
  n_splits <- function(xnode) {
    wi <- which("Node" == names(xnode))
    rval <- unique(sapply(wi, function(i) {
      if(any(c("SimplePredicate", "SimpleSetPredicate") %in% names(xnode[[i]]))) return(1)
      sum(names(xnode[[i]][["CompoundPredicate"]]) %in% c("SimplePredicate", "SimpleSetPredicate"))
    }))
    if(length(rval) > 1) stop("malformatted XML")
    return(rval)
  }
  kid_ids <- function(xnode) {
    wi <- which("Node" == names(xnode))
    rval <- sapply(wi, function(j) {
      as.vector(xmlAttrs(xnode[[j]])["id"])
    })  
  }
  get_pred <- function(xnode) {
    pred <- as.vector(xmlAttrs(xnode)["score"])
    if(is.na(pred)) return(NULL)
    if(is.numeric(mf_response)) as.numeric(pred)
      else factor(pred, levels = levels(mf_response))
  }
  get_dist <- function(xnode) {
    wi <- which("ScoreDistribution" == names(xnode))
    if(length(wi) < 1) return(NULL)
    rval <- sapply(wi, function(i) as.numeric(xmlAttrs(xnode[[i]])["recordCount"]))
    names(rval) <- sapply(wi, function(i) xmlAttrs(xnode[[i]])["value"])
    if(inherits(mf_response, "factor")) rval <- rval[levels(mf_response)]
    return(rval)
  }
  get_error <- function(xnode) {
    if(tm_info["functionName"] != "classification") return(NULL)
    tab <- get_dist(xnode)
    if(is.null(tab)) return(NULL)
    c("%" = sum(100 * prop.table(tab)[names(tab) != get_pred(xnode)]))
  }
  node_info <- function(xnode) list(prediction = get_pred(xnode), n = n_obs(xnode),
    error = get_error(xnode), distribution = get_dist(xnode))
  get_split_prob <- function(xnode) {
    rval <- rep(0, n_kids(xnode))
    wi <- xmlAttrs(xnode)["defaultChild"]
    if(is.na(wi)) rval <- NULL
      else rval[which(kid_ids(xnode) == wi)] <- 1
    return(rval)
  }
  get_split <- function(xnode, i, surrogates) {
    wi <- which("Node" == names(xnode))
    rval <- sapply(wi, function(j) {
      nj <- if(surrogates) xnode[[j]][["CompoundPredicate"]] else xnode[[j]]
      wii <- which(names(nj) %in% c("SimplePredicate", "SimpleSetPredicate"))[i]      
      c("predicateType" = names(nj)[wii], xmlAttrs(nj[[wii]]))
    })
    stopifnot(length(unique(rval["predicateType",])) == 1)
    stopifnot(length(unique(rval["field",])) == 1)
    if(rval["predicateType", 1] == "SimplePredicate") {
      stopifnot(length(unique(rval["value",])) == 1)
      if(ncol(rval) != 2) stop("not yet implemented")
      if(!(identical(sort(rval["operator",]), c("greaterThan", "lessOrEqual")) |
           identical(sort(rval["operator",]), c("greaterOrEqual", "lessThan")))
      ) stop("not yet implemented")
      partysplit(
        varid = which(rval["field", 1] == mf_names),
	breaks = as.numeric(rval["value", 1]),
	index = if(substr(rval["operator", 1], 1, 1) != "l") 2:1 else NULL,
	right = "lessOrEqual" %in% rval["operator",],
	prob = if(i == 1) get_split_prob(xnode) else NULL
      )      
    } else {
      varid <- which(rval["field", 1] == mf_names)
      lev <- mf_levels[[varid]]
      stopifnot(length(lev) > 1)
      idx <- rep(0, length(lev))
      lab <- lapply(wi, function(j) {
        nj <- if(surrogates) xnode[[j]][["CompoundPredicate"]] else xnode[[j]]
        wii <- which(names(nj) %in% c("SimplePredicate", "SimpleSetPredicate"))[i]      
        ar <- nj[[wii]][["Array"]]
	stopifnot(xmlAttrs(ar)["type"] == "string")
	rv <- strsplit(xmlValue(ar), " ")[[1]]
	stopifnot(length(rv) == as.numeric(xmlAttrs(ar)["n"]))
	return(rv)
      })
      for(j in 1:ncol(rval)) {
        if(rval["booleanOperator",j] == "isIn") idx[which(lev %in% lab[[j]])] <- j
	  else idx[which(!(lev %in% lab[[j]]))] <- j
      }
      stopifnot(all(idx > 0))
      partysplit(
        varid = varid,
	breaks = NULL,
	index = as.integer(idx),
	prob = if(i == 1) get_split_prob(xnode) else NULL
      )
    }
  }
  
  ## function for setting up nodes
  ## (using global index ii)
  pmml_node <- function(xnode) {
    ii <<- ii + 1
    if(is_terminal(xnode)) return(partynode(as.integer(ii),
      info = node_info(xnode)
    ))
    wi <- which("Node" == names(xnode))
    ns <- n_splits(xnode)    
    nd <- partynode(as.integer(ii),
      split = get_split(xnode, 1, ns > 1),
      kids = lapply(wi, function(j) pmml_node(xnode[[j]])),
      surrogates = if(ns < 2) NULL else lapply(2:ns, function(j) get_split(xnode, j, TRUE)),
      info = node_info(xnode)
    )
    nd
  }
  
  ## set up node
  ii <- 0
  if(is_root(tm[["Node"]])) nd <- pmml_node(tm[["Node"]]) else stop("mal-formed XML")

  ## set up party
  pt <- party(node = nd, data = mf, fitted = NULL, terms = trms, names = NULL) ## FIXME: info slot
  class(pt) <- c("simple_party", class(pt))

  return(pt)
}
