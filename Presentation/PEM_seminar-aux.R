
## Auxiliary functions

## Shows the vertices of a graph, optionally highlighting an arbitrary subset
## among them.
showVertices <- function(x, highlight = NULL, label = FALSE,
                         lab.trim = c(0,0.08)) {
  
  bg <- rep(0, attr(x,"ev")[2L])
  
  if(!is.null(highlight))
    bg[highlight] <- 1
  
  points(x=x$vertex$x, y=x$vertex$y, pch=21L, bg=bg)
  
  if(label)
    text(
      x = x$vertex$x + lab.trim[1L],
      y = x$vertex$y + lab.trim[2L],
      labels = attr(x,"vlabel"),
      xpd = TRUE
    )
  
  invisible(NULL)
  
}

## Shows the edges of a graph as dotted lines between the vertices, optionally
## highlighting an arbitrary subset among them.
showEdges <- function(x, highlight = NULL, label = FALSE, lab.trim = c(0,0)) {
  
  lty <- rep(3L, attr(x,"ev")[1L])
  
  if(!is.null(highlight))
    lty[highlight] <- 1L
  
  segments(
    x0 = x$vertex$x[x$edge[[1L]]],
    x1 = x$vertex$x[x$edge[[2L]]],
    y0 = x$vertex$y[x$edge[[1L]]],
    y1 = x$vertex$y[x$edge[[2L]]],
    lty = lty
  )
  
  if(label) {
    text(
      x = 0.5*(x$vertex$x[x$edge[[1L]]] + x$vertex$x[x$edge[[2L]]]) + lab.trim[1L],
      y = 0.5*(x$vertex$y[x$edge[[1L]]] + x$vertex$y[x$edge[[2L]]]) + lab.trim[2L],
      labels = attr(x,"elabel"),
      xpd = TRUE
    )
  }
  
  invisible(NULL)
  
}

## Shows the influence matrix of a graph, optionally highlighting an arbitrary
## subset of its elements.
showMatrix <- function(x, highlight = NULL) {
  
  font <- matrix(1L, nrow=nrow(x), ncol=ncol(x))
  
  if(!is.null(highlight))
    font[highlight] <- 2L
  
  if(!is.null(rownames(x)))
    text(
      x = rep(0,nrow(x)),
      y = rev((1L:nrow(x))/nrow(x)),
      label = rownames(x)
    )
  
  if(!is.null(colnames(x)))
    text(
      x = (1L:ncol(x))/ncol(x),
      y = rep(1 + 1/nrow(x), ncol(x)),
      label = colnames(x),
      xpd = TRUE
    )
  
  for(i in 1L:nrow(x))
    for(j in 1L:ncol(x))
      text(
        x = j/ncol(x),
        y = (nrow(x) - i + 1L)/nrow(x),
        label = x[i,j],
        font = font[i,j]
      )
  
  invisible(NULL)
  
}

## Draws a graph and its influence matrix as an adjacent inset.
plot1_graph <- function(x, placeHolder) {
  
  par(mfrow = c(1L,2L), mar=c(1,1,2,1))
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  showEdges(x$graph, label = TRUE)
  showVertices(x$graph, label = TRUE)
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  tmp <- x$Infl
  storage.mode(tmp) <- "character"
  if(!missing(placeHolder))
    tmp[] <- placeHolder
  showMatrix(tmp)
  
  invisible(NULL)
  
}

## Draws a drop line going from a the specified vertex towards the corresponding
## row of an adjacently displayed influence matrix.
droplineVertex <- function(x, v, ...) {
  
  segments(
    x0 = x$vertex$x[v],
    y0 = x$vertex$y[v],
    x1 = x$vertex$x[v],
    y1 = (attr(x,"ev")[2L] - v + 1L)/attr(x,"ev")[2L],
    xpd = TRUE
  )
  
  arrows(
    x0 = x$vertex$x[v],
    y0 = (attr(x,"ev")[2L] - v + 1L)/attr(x,"ev")[2L],
    x1 = 1.1,
    y1 = (attr(x,"ev")[2L] - v + 1L)/attr(x,"ev")[2L],
    xpd = TRUE,
    ...
  )
  
  invisible(NULL)
  
}

## Draws a drop line going from the specified edge towards the lower end of
## the plotting area.
droplineEdgeGraph <- function(x, e, ...) {
  
  segments(
    x0 = 0.5*(x$vertex$x[x$edge[[1L]][e]] + x$vertex$x[x$edge[[2L]][e]]),
    x1 = 0.5*(x$vertex$x[x$edge[[1L]][e]] + x$vertex$x[x$edge[[2L]][e]]),
    y0 = 0.5*(x$vertex$y[x$edge[[1L]][e]] + x$vertex$y[x$edge[[2L]][e]]),
    y1 = -0.1,
    xpd = TRUE,
    ...
  )
  
  segments(
    x0 = 0.5*(x$vertex$x[x$edge[[1L]][e]] + x$vertex$x[x$edge[[2L]][e]]),
    x1 = 1.1,
    y0 = -0.1,
    y1 = -0.1,
    xpd = TRUE,
    ...
  )
  
  invisible(NULL)
  
}

## Draws a drop line going from the lower end of the plotting area towards the
## the specified column of an adjacently displayed influence matrix.
droplineEdgeMatrix <- function(x, e, ...) {
  
  segments(
    x0 = -0.1,
    x1 = e/attr(x,"ev")[1L],
    y0 = -0.1,
    y1 = -0.1,
    xpd = TRUE,
    ...
  )
  
  arrows(
    x0 = e/attr(x,"ev")[1L],
    x1 = e/attr(x,"ev")[1L],
    y0 = -0.1,
    y1 = 0.05,
    xpd = TRUE,
    length = 0.1
  )
  
  invisible(NULL)
  
}

## Plot a graph and where a given vertex is located.
plot1_vertex <- function(x, v, placeHolder = "\"", highlight = "x") {
  
  par(mfrow = c(1L,2L), mar=c(1,1,2,1))
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  showEdges(x$graph)
  droplineVertex(x$graph, v, length=0.1)
  showVertices(x$graph, highlight = v)
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  tmp <- x$Infl
  storage.mode(tmp) <- "character"
  tmp[-v,] <- placeHolder
  tmp[v,] <- highlight
  showMatrix(tmp, highlight = cbind(v,1L:attr(x$graph,"ev")[1L]))
  
  invisible(NULL)
  
}

## Plot a graph and where a given edge is located.
plot1_edge <- function(x, e, placeHolder = "\"", highlight = "x") {
  
  par(mfrow = c(1L,2L), mar=c(1,1,2,1))
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  showEdges(x$graph, e)
  showVertices(x$graph)
  droplineEdgeGraph(x$graph, e, lty=2L)
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  tmp <- x$Infl
  storage.mode(tmp) <- "character"
  tmp[,-e] <- placeHolder
  tmp[,e] <- highlight
  showMatrix(tmp, highlight = cbind(1L:attr(x$graph,"ev")[2L],e))
  droplineEdgeMatrix(x$graph, e, lty=2L)
  
  invisible(NULL)
  
}

plot1_influence <- function(x, v, e, placeHolder = "\"") {
  
  par(mfrow = c(1L,2L), mar=c(1,1,2,1))
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  if(missing(e)) {
    showEdges(x$graph, which(!!x$Infl[v,]), label=TRUE)
  } else
    showEdges(x$graph, e)
  droplineVertex(x$graph, v, length=0.1)
  showVertices(x$graph, highlight = v)
  if(!missing(e))
    droplineEdgeGraph(x$graph, e, lty=2L)
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  tmp <- x$Infl
  storage.mode(tmp) <- "character"
  tmp[-(1L:v),] <- placeHolder
  if(!missing(e)) {
    tmp[v,-(1L:e)] <- placeHolder
    showMatrix(tmp, highlight = cbind(v,e))
    droplineEdgeMatrix(x$graph, e, lty=2L)
  } else
    showMatrix(tmp, highlight = cbind(v,1L:attr(x$graph,"ev")[1L]))
  
  invisible(NULL)
  
}
