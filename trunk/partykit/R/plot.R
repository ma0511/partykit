nobs_party <- function(party, id = 1) {
  dat <- data_party(party, id = id)
  if("(weights)" %in% names(dat)) sum(dat[["(weights)"]]) else NROW(dat)
}

node_inner <- function(obj, id = TRUE, abbreviate = FALSE, fill = "white", gp = gpar())
{
  meta <- obj$data
  nam <- names_party(obj)

  extract_label <- function(node) {
    if(is.terminal(node)) return(rep.int("", 2))

    varlab <- character_split(split_node(node), meta)$name
    if(abbreviate > 0) varlab <- abbreviate(varlab, as.numeric(abbreviate))

    ## FIXME: info processing
    plab <- ""
    return(c(varlab, plab))
  }

  maxstr <- function(node) {
      lab <- extract_label(node)
      klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
      lab <- c(lab, klab)
      lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
      return(lab[which.max(nchar(lab))])
  }

  nstr <- maxstr(node_party(obj))

  ### panel function for the inner nodes
  rval <- function(node) {  
    node_vp <- viewport(
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "strwidth", nstr) * 1.3, 
      height = unit(3, "lines"),
      name = paste("node_inner", id_node(node), sep = ""),
      gp = gp
    )
    pushViewport(node_vp)

    xell <- c(seq(0, 0.2, by = 0.01),
  	      seq(0.2, 0.8, by = 0.05),
  	      seq(0.8, 1, by = 0.01))
    yell <- sqrt(xell * (1-xell))

    lab <- extract_label(node)
    fill <- rep(fill, length.out = 2)

    grid.polygon(x = unit(c(xell, rev(xell)), "npc"),
        	 y = unit(c(yell, -yell)+0.5, "npc"),
        	 gp = gpar(fill = fill[1]))

    ## FIXME: something instead of pval ?
    grid.text(lab[1], y = unit(1.5 + 0.5 * FALSE, "lines"))
    if(FALSE) grid.text(lab[2], y = unit(1, "lines"))

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
        width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
        height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2]))
      grid.text(nam[id_node(node)])
      popViewport()
    }
    upViewport()
  }
  
  return(rval)
}
class(node_inner) <- "grapcon_generator"

node_terminal <- function(obj,
                          digits = 3,
		          abbreviate = FALSE,
		          fill = c("lightgray", "white"),
		          id = TRUE,
			  gp = gpar())
{
  nam <- names_party(obj)

  extract_label <- function(node) {
    return(c("terminal", "node"))
    ## FIXME: re-use print method
  }

  maxstr <- function(node) {
      lab <- extract_label(node)
      klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
      lab <- c(lab, klab)
      lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
      return(lab[which.max(nchar(lab))]) ## FIXME: nchar?
  }

  nstr <- maxstr(node_party(obj))

  ### panel function for simple n, Y terminal node labeling
  rval <- function(node) {
    fill <- rep(fill, length.out = 2)	    

    node_vp <- viewport(x = unit(0.5, "npc"),	
      y = unit(0.5, "npc"),   
      width = unit(1, "strwidth", nstr) * 1.1,
      height = unit(3, "lines"),
      name = paste("node_terminal", id_node(node), sep = ""),
      gp = gp
    )
    pushViewport(node_vp)

    lab <- extract_label(node)

    grid.rect(gp = gpar(fill = fill[1]))
    grid.text(y = unit(2, "lines"), lab[1])
    grid.text(y = unit(1, "lines"), lab[2])

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
        width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
        height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2], lty = "solid"))
      grid.text(nam[id_node(node)])
      popViewport()
    }
    upViewport()
  }
  return(rval)
}
class(node_terminal) <- "grapcon_generator"

edge_simple <- function(obj, digits = 3, abbreviate = FALSE)
{
  meta <- obj$data

  ### panel function for simple edge labelling
  function(node, i) {
    split <- character_split(split_node(node), meta)$levels[i]
    ## FIXME: can this be improved?
    if(any(grep(">=", split) > 0) | any(grep("<=", split) > 0))
      split <- parse(text = paste("phantom(0)", split))
    grid.rect(gp = gpar(fill = "white", col = 0), width = unit(1, "strwidth", split)) 
    grid.text(split, just = "center")
  }
}
class(edge_simple) <- "grapcon_generator"

plot_node <- function(node, xlim, ylim, nx, ny, 
               terminal_panel, inner_panel, edge_panel,
	       tnex = 2, drop_terminal = TRUE, debug = FALSE) {

    ### the workhorse for plotting trees
 
    ### set up viewport for terminal node
    if (is.terminal(node)) {
        x <- xlim[1] + diff(xlim)/2
        y <- ylim[1] + 0.5
       
        tn_vp <- viewport(x = unit(x, "native"),
                          y = unit(y, "native") - unit(0.5, "lines"),
                          width = unit(1, "native"), 
                          height = unit(tnex, "native") - unit(1, "lines"),
			  just = c("center", "top"),
                          name = paste("Node", id_node(node), sep = ""))
        pushViewport(tn_vp)
        if (debug)
            grid.rect(gp = gpar(lty = "dotted", col = 4))
        terminal_panel(node) 
        upViewport()
        return(NULL)
    }    

    ## convenience function for computing relative position of splitting node
    pos_frac <- function(node) {
      if(is.terminal(node)) 0.5 else {
        width_kids <- sapply(kids_node(node), width)
        nk <- length(width_kids)
        rval <- if(nk %% 2 == 0) sum(width_kids[1:(nk/2)]) else
	  mean(cumsum(width_kids)[nk/2 + c(-0.5, 0.5)])
	rval/sum(width_kids)
      }
    }

    ## extract information
    split <- split_node(node)
    kids <- kids_node(node)
    width_kids <- sapply(kids, width)
    nk <- length(width_kids)

    ### position of inner node
    x0 <- xlim[1] + pos_frac(node) * diff(xlim)
    y0 <- max(ylim)

    ### relative positions of kids
    xfrac <- sapply(kids, pos_frac)
    x1lim <- xlim[1] + cumsum(c(0, width_kids))/sum(width_kids) * diff(xlim)
    x1 <- x1lim[1:nk] + xfrac * diff(x1lim)
    if (!drop_terminal) {
        y1 <- rep(y0 - 1, nk)
    } else {
        y1 <- ifelse(sapply(kids, is.terminal), tnex - 0.5, y0 - 1)
    }

    ### draw edges
    for(i in 1:nk) grid.lines(x = unit(c(x0, x1[i]), "native"), y = unit(c(y0, y1[i]), "native"))

    ### create viewport for inner node
    in_vp <- viewport(x = unit(x0, "native"),
                      y = unit(y0, "native"),
                      width = unit(1, "native"),
                      height = unit(1, "native") - unit(1, "lines"), 
                      name = paste("Node", id_node(node), sep = ""))
    pushViewport(in_vp)
    if(debug) grid.rect(gp = gpar(lty = "dotted"))
    inner_panel(node)
    upViewport()

    ### position of labels
    y1max <- max(y1)
    ypos <- y0 - (y0 - y1max) * 0.5
    xpos <- x0 - (x0 - x1) * 0.5 * (y0 - y1max)/(y0 - y1)

    ### setup labels
    for(i in 1:nk) {
      sp_vp <- viewport(x = unit(xpos[i], "native"),
                        y = unit(ypos, "native"),
                        width = unit(diff(x1lim)[i], "native"),
                        height = unit(1, "lines"), 
                        name =  paste("edge", id_node(node), "-", i, sep = ""))
      pushViewport(sp_vp)
      if(debug) grid.rect(gp = gpar(lty = "dotted", col = 2))
      edge_panel(node, i)
      upViewport()
    }

    ## call workhorse for kids
    for(i in 1:nk) plot_node(kids[[i]],
      c(x1lim[i], x1lim[i+1]), c(y1[i], 1), nx, ny, 
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}


plot.party <- function(x, main = NULL, type = "simple", ## FIXME: remove, was: c("extended", "simple"),
                       terminal_panel = NULL, tp_args = list(),
		       inner_panel = node_inner, ip_args = list(),
                       edge_panel = edge_simple, ep_args = list(),
		       drop_terminal = (type[1] == "extended"),
		       tnex = (type[1] == "extended") + 1, 
		       newpage = TRUE,
		       pop = TRUE,
		       gp = gpar(),
		       ...) {

    ### extract tree
    node <- node_party(x)
    ### total number of terminal nodes
    nx <- width(node)
    ### maximal depth of the tree
    ny <- depth(node)

    ### compute default settings
    type <- match.arg(type)
    if (type == "simple") {
        if (is.null(terminal_panel)) 
            terminal_panel <- node_terminal
        if (is.null(tnex)) tnex <- 1
    } else {
        if (is.null(terminal_panel))
            terminal_panel <- switch(class(x$fitted[["(response)"]]), ## FIXME: re-structure into plot.cparty
	                             "Surv" = node_surv,
                                     "factor" = node_barplot,
                                     "ordered" = node_barplot,
                                     node_boxplot)
        if (is.null(tnex)) tnex <- 2
    }

    ## setup newpage
    if (newpage) grid.newpage()

    ## setup root viewport
    root_vp <- viewport(layout = grid.layout(3, 3, 
    			height = unit(c(ifelse(is.null(main), 0, 3), 1, 1), 
                                      c("lines", "null", "lines")),
    			width = unit(c(1, 1, 1), 
                                     c("lines", "null", "lines"))), 
    			name = "root",
			gp = gp)       
    pushViewport(root_vp)
  
    ## viewport for main title (if any)
    if (!is.null(main)) {
        main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1, 
                            name = "main")
        pushViewport(main_vp)
        grid.text(y = unit(1, "lines"), main, just = "center")
        upViewport()
    }

    ## setup viewport for tree
    tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
    			xscale = c(0, nx), yscale = c(0, ny + (tnex - 1)), 
                        name = "tree")
    pushViewport(tree_vp)

    ### setup panel functions (if necessary)
    ### the heuristic is as follows: If the first argument
    ### is `obj' than we assume a panel generating function, 
    ### otherwise the function is treated as a panel function
    if(inherits(terminal_panel, "grapcon_generator"))
      terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
    if(inherits(inner_panel, "grapcon_generator"))
      inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
    if(inherits(edge_panel, "grapcon_generator"))
      edge_panel <- do.call("edge_panel", c(list(x), as.list(ep_args)))


    if((nx <= 1 & ny <= 1)) {
      pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", id_node(node), sep = "")))
      terminal_panel(node)
    } else {
      ## call the workhorse
      plot_node(node,
        xlim = c(0, nx), ylim = c(0, ny - 0.5 + (tnex - 1)),
        nx = nx, ny = ny, 
        terminal_panel = terminal_panel,
        inner_panel = inner_panel,
        edge_panel = edge_panel,
        tnex = tnex,
        drop_terminal = drop_terminal,
        debug = FALSE)
    }
    upViewport()
    if (pop) popViewport() else upViewport()
}
