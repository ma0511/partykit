node_inner <- function(obj, id = TRUE, abbreviate = FALSE, fill = "white")
{
  meta <- metadata(obj)

  extract_label <- function(node) {
    if(is.terminal(node)) return(rep.int("", 2))

    varlab <- nodelabels(split_node(node), meta)$splitname
    if(abbreviate > 0) varlab <- abbreviate(varlab, as.numeric(abbreviate))

    if(FALSE) { ## FIXME: p-value processing
        pvalue <- 1 - node$criterion$maxcriterion
        plab <- ifelse(pvalue < 10^(-digits),
        	       paste("p <", 10^(-digits)),
        	       paste("p =", round(pvalue, digits = digits)))
    } else {
  	plab <- ""
    }
    return(c(varlab, plab))
  }

  maxstr <- function(node) {
      lab <- extract_label(node)
      klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
      lab <- c(lab, klab)
      return(lab[which.max(nchar(lab))]) ## FIXME: nchar?
  }

  nstr <- maxstr(get_node(obj))

  ### panel function for the inner nodes
  rval <- function(node) {  
    node_vp <- viewport(
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "strwidth", nstr) * 1.3, 
      height = unit(3, "lines"),
      name = paste("node_inner", id_node(node), sep = "")
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
        width = max(unit(1, "lines"), unit(1.3, "strwidth", as.character(id_node(node)))),
        height = max(unit(1, "lines"), unit(1.3, "strheight", as.character(id_node(node)))))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2]))
      grid.text(id_node(node))
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
		          id = TRUE)
{
  extract_label <- function(node) {
    return(c("", ""))
    
    if(!is.terminal(node)) return(rep.int("", 2))
    nlab <- paste("n =", sum(node$weights))
    ylab <- if(length(node$prediction) > 1)
        	paste("y =", paste("(", paste(round(node$prediction, digits),
  		      collapse = ", "), ")", sep =""))
            else
  		paste("y =", round(node$prediction, digits))
    return(c(nlab, ylab))
  }

  maxstr <- function(node) {
      lab <- extract_label(node)
      klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
      lab <- c(lab, klab)
      return(lab[which.max(nchar(lab))]) ## FIXME: nchar?
  }

  nstr <- maxstr(get_node(obj))

  ### panel function for simple n, Y terminal node labeling
  rval <- function(node) {
    fill <- rep(fill, length.out = 2)	    

    node_vp <- viewport(x = unit(0.5, "npc"),	
      y = unit(0.5, "npc"),   
      width = unit(1, "strwidth", nstr) * 1.1,
      height = unit(3, "lines"),
      name = paste("node_terminal", id_node(node), sep = "")
    )
    pushViewport(node_vp)

    lab <- extract_label(node)

    grid.rect(gp = gpar(fill = fill[1]))
    grid.text(y = unit(2, "lines"), lab[1])
    grid.text(y = unit(1, "lines"), lab[2])

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
        width = max(unit(1, "lines"), unit(1.3, "strwidth", as.character(id_node(node)))),
        height = max(unit(1, "lines"), unit(1.3, "strheight", as.character(id_node(node)))))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2], lty = "solid"))
      grid.text(id_node(node))
      popViewport()
    }
    upViewport()
  }
  return(rval)
}
class(node_terminal) <- "grapcon_generator"

edge_simple <- function(obj, digits = 3, abbreviate = FALSE)
{
  meta <- metadata(obj)

  ### panel function for simple edge labelling
  function(node, i) {
    split <- nodelabels(split_node(node), meta)$splitlevels[i]

	### <FIXME> phantom and . functions cannot be found by
	###	    codetools
	### </FIXME>
	## if (left) split <- as.expression(bquote(phantom(0) <= .(split)))
	##    else split <- as.expression(bquote(phantom(0) > .(split)))

    grid.rect(gp = gpar(fill = "white", col = 0), width = unit(1, "strwidth", split)) 
    grid.text(split, just = "center")
  }
}
class(edge_simple) <- "grapcon_generator"

plot_party <- function(node, xlim, ylim, nx, ny, 
               terminal_panel, inner_panel, edge_panel,
	       tnex = 2, drop_terminal = TRUE, debug = FALSE) {

    ### the workhorse for plotting trees

    ### set up viewport for terminal node
    if (node$terminal) {
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

    ### number of left leafs
    nl <- width(node$left)

    ### number of right leafs
    nr <- width(node$right)

    ### position of inner node
    x0 <- xlim[1] + (nl / (nl + nr)) * diff(xlim)
    y0 <- max(ylim)

    ### proportion of left terminal nodes in left node
    if (node$left$terminal) {
        lf <- 1/2
    } else {
        lf <- width(node$left$left) / (width(node$left$left) + 
                                           width(node$left$right))
    }

    ### proportion of left terminal nodes in right node
    if (node$right$terminal) {
        rf <- 1/2
    } else {
        rf <- width(node$right$left) / (width(node$right$left) + 
                                            width(node$right$right))
    }

    ### position of left and right daugher node
    x1l <- xlim[1] + (x0 - xlim[1]) * lf
    x1r <- x0 + (xlim[2] - x0) * rf
    
    if (!drop_terminal) {
        y1l <- y1r <- y0 - 1
    } else {
        y1l <- if (node$left$terminal) tnex - 0.5 else y0 - 1
        y1r <- if (node$right$terminal) tnex - 0.5 else y0 - 1
    }

    ### draw edges
    grid.lines(x = unit(c(x0, x1l), "native"), 
               y = unit(c(y0, y1l), "native"))
    grid.lines(x = unit(c(x0, x1r), "native"), 
               y = unit(c(y0, y1r), "native"))

    ### create viewport for inner node
    in_vp <- viewport(x = unit(x0, "native"),
                      y = unit(y0, "native"),
                      width = unit(1, "native"),
                      height = unit(1, "native") - unit(1, "lines"), 
                      name = paste("Node", id_node(node), sep = ""))
    pushViewport(in_vp)
    if (debug)
        grid.rect(gp = gpar(lty = "dotted"))
    inner_panel(node)
    upViewport()

    ps <- node$psplit
    if (ps$ordered) {
        if (!is.null(attr(ps$splitpoint, "levels"))) {
            split <- attr(ps$splitpoint, "levels")[ps$splitpoint]
        } else {
            split <- ps$splitpoint
        }
    } else {
        ### <FIXME>: always to the left? </FIXME>
        split <- attr(ps$splitpoint, "levels")[as.logical(ps$splitpoint) & (ps$table > 0)]
    }


    ### position of labels
    y1lr <- max(y1l, y1r)
    ypos <- y0 - (y0 - y1lr) * 0.5
    xlpos <- x0 - (x0 - x1l) * 0.5 * (y0 - y1lr)/(y0 - y1l)
    xrpos <- x0 - (x0 - x1r) * 0.5 * (y0 - y1lr)/(y0 - y1r)

    ### setup left label
    lsp_vp <- viewport(x = unit(xlpos, "native"),
                       y = unit(ypos, "native"),
                       width = unit(xlpos - xrpos, "native"),
                       height = unit(1, "lines"), 
                       name =  paste("lEdge", id_node(node), sep = ""))
    pushViewport(lsp_vp)
    if (debug)
        grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_panel(split, ordered = ps$ordered, left = TRUE)
    upViewport()

    ### setup right label
    if (ps$ordered) {
        if (!is.null(attr(ps$splitpoint, "levels"))) {
            split <- attr(ps$splitpoint, "levels")[ps$splitpoint]
        } else {
            split <- ps$splitpoint
        }
    } else {
        split <- attr(ps$splitpoint, "levels")[!as.logical(ps$splitpoint) & (ps$table > 0)]
    }

    rsp_vp <- viewport(x = unit(xrpos, "native"),
                       y = unit(ypos, "native"),
                       width = unit(xlpos - xrpos, "native"),
                       height = unit(1, "lines"),
                       name =  paste("rEdge", id_node(node), sep = ""))
    pushViewport(rsp_vp) 
    if (debug)
        grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_panel(split, ordered = ps$ordered, left = FALSE)
    upViewport()

    plotTree(node$left, c(xlim[1], x0), c(y1l, 1), nx, ny, 
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
    plotTree(node$right, c(x0, xlim[2]), c(y1r, 1), nx, ny,
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}


plot.party <- function(x, main = NULL, type = c("extended", "simple"),
                            terminal_panel = NULL, tp_args = list(),
			    inner_panel = node_inner, ip_args = list(),
                            edge_panel = edge_simple, ep_args = list(),
			    drop_terminal = (type[1] == "extended"),
			    tnex = (type[1] == "extended") + 1, 
			    newpage = TRUE,
			    pop = TRUE,
			    ...) {

    ### plot BinaryTree objects

    ### extract tree
    ptr <- x@tree
    ### total number of terminal nodes
    nx <- width(ptr)
    ### maximal depth of the tree
    ny <- maxdepth(ptr)

    ### compute default settings
    type <- match.arg(type)
    if (type == "simple") {
        if (is.null(terminal_panel)) 
            terminal_panel <- node_terminal
        if (is.null(tnex)) tnex <- 1
    } else {
        if (is.null(terminal_panel))
            terminal_panel <- switch(class(response(x)[[1]])[1],
	                             "Surv" = node_surv,
                                     "factor" = node_barplot,
                                     "was_ordered" = node_barplot,
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
    			name = "root")       
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
      pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", ptr$nodeID, sep = "")))
      terminal_panel(ptr)    
    } else {
      ## call the workhorse
      plotTree(ptr,
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
