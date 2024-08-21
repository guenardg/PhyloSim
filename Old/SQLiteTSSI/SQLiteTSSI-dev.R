
## SQLite-based reticulated evolution API

## rm(list=ls())
## library(yaml)
## library(Rcpp)
library(PhyloSim)  ## detach("package:PhyloSim", unload=TRUE)
## library(RSQLite)
## library(rgl)

## recalculateExamples <- FALSE
## recalculateExamples <- TRUE
## source("TSSI-examples.R")
## rm(recalculateExamples)

## source("aux.R")

openNetwork(
  filename = system.file(
    package = "PhyloSim",
    "extdata",
    "reticulated-ex2.sim"
  ),
  load = TRUE
) -> net

## net <- openNetwork("reticulated-ex1.sim", load=TRUE)

## closeNetwork(net, FALSE)

##
### How to display a network?
### Assuming time in the abscissa from left to right (from Western bias).
##

### We have the edge information (origin, destination, and lengths) from the
### graph object:
gr <- getGraph(net)

### We have the average patristic distances:
dst <- getPatristic(net)

### This is the distances from the graph root. This is the positions of the
### vertices in abscissa:
dd <- c(0,dst[dst_idx(attr(gr,"ev")[2L],1L,)])

### We have the vertices that are leaves:
leaf <- getLeaves(net)


tol <- .Machine$double.eps^0.5
apply(
  MPSEM::PEMInfluence(gr),
  1L,
  function(x, y)
    x*y,
  y = sqrt(gr$edge$distance)
) %>%
  t %>%
  scale(scale = FALSE) -> ii

ss <- svd(ii, nu=1L, nv=0L)

a <- as.numeric(ss$u)
a[order(a)] <- seq(0, 1, length.out=length(a))

data.frame(
  x = dd,
  y = a,
  leaf = FALSE
) -> coords
coords$leaf[leaf] <- TRUE

par(mar=c(1,1,1,1))
plot(NA, xlim=range(coords$x), ylim=range(coords$y), axes=FALSE, xlab="",
     ylab="")
points(coords, pch=21L, bg="black", cex=ifelse(coords$leaf, 0.75, 0.25))
for(i in 1L:attr(gr,"ev")[1L])
  arrows(
    x0 = coords[gr$edge[[1L]][i],"x"],
    y0 = coords[gr$edge[[1L]][i],"y"],
    x1 = coords[gr$edge[[2L]][i],"x"],
    y1 = coords[gr$edge[[2L]][i],"y"],
    length = 0.05
  )




ss <- svd(ii, nu=2L, nv=0L)

a <- atan2(ss$u[,2L],ss$u[,1L])
a[order(a)] <- head(seq(-pi,pi,length=101),100L)

data.frame(
  x = dd*cos(a),
  y = dd*sin(a),
  leaf = FALSE
) -> coords2
coords2$leaf[leaf] <- TRUE

par(mar=c(1,1,1,1))
plot(NA, xlim=range(coords2$x), ylim=range(coords2$y), axes=FALSE, xlab="",
     ylab="")
points(coords2, pch=21L, bg="black", cex=ifelse(coords2$leaf, 1, 0.5))
for(i in 1L:attr(gr,"ev")[1L])
  arrows(
    x0 = coords2[gr$edge[[1L]][i],"x"],
    y0 = coords2[gr$edge[[1L]][i],"y"],
    x1 = coords2[gr$edge[[2L]][i],"x"],
    y1 = coords2[gr$edge[[2L]][i],"y"],
    length = 0.08
  )











### We have the lineages for each of the vertices:
lin <- getLineage(net, v = 1L:attr(gr,"ev")[2L])


## ll[leaf]
## lin[leaf]


## dd[leaf]









G <- -0.5*as.matrix(dst)^2
delta <- t(t(G - rowMeans(G)) - colMeans(G)) + mean(G)
eig <- eigen(delta)

data.frame(
  x = dd,
  y = apply(eig$vectors[,1L:10] %*% diag(eig$values[1L:10L]^-1), 1L, sum),
  leaf = FALSE
) -> coords
coords$leaf[leaf] <- TRUE

### Getting the best bet for the positions in ordinates is a harder call.

par(mar=c(4,1,1,1))
plot(NA, xlim=range(coords$x), ylim=range(coords$y), axes=FALSE,
     xlab="Time", ylab="")
points(coords, pch=21L, bg="black", cex=ifelse(coords$leaf, 1, 0.5))
for(i in 1L:attr(gr,"ev")[1L])
  arrows(
    x0 = coords[gr$edge[[1L]][i],"x"],
    y0 = coords[gr$edge[[1L]][i],"y"],
    x1 = coords[gr$edge[[2L]][i],"x"],
    y1 = coords[gr$edge[[2L]][i],"y"],
    length = 0.1
  )








## Not overly impressive but works somehow

tol <- .Machine$double.eps^0.5
Infl <- MPSEM::PEMInfluence(gr)
ii <- scale(Infl, scale = FALSE)
ss <- svd(ii, nu=2L, nv=0L)

data.frame(
  x = dd,
  y = NA,
  leaf = FALSE
) -> coords
coords$leaf[leaf] <- TRUE
coords$y <- ss$u[,1L]

### Getting the best bet for the positions in ordinates is a harder call.

par(mar=c(4,1,1,1))
plot(NA, xlim=range(coords$x), ylim=range(coords$y), axes=FALSE,
     xlab="Time", ylab="")
points(coords, pch=21L, bg="black", cex=ifelse(coords$leaf, 1, 0.5))
for(i in 1L:attr(gr,"ev")[1L])
  arrows(
    x0 = coords[gr$edge[[1L]][i],"x"],
    y0 = coords[gr$edge[[1L]][i],"y"],
    x1 = coords[gr$edge[[2L]][i],"x"],
    y1 = coords[gr$edge[[2L]][i],"y"],
    length = 0.1
  )




a <- atan2(ss$u[,2L],ss$u[,1L])
a[order(a)] <- head(seq(-pi,pi,length=101),100L)

data.frame(
  x = dd*cos(a),
  y = dd*sin(a),
  leaf = FALSE
) -> coords2
coords2$leaf[leaf] <- TRUE

par(mar=c(1,1,1,1))
plot(NA, xlim=range(coords2$x), ylim=range(coords2$y), axes=FALSE, xlab="",
     ylab="")
points(coords2, pch=21L, bg="black", cex=ifelse(coords2$leaf, 1, 0.5))
for(i in 1L:attr(gr,"ev")[1L])
  arrows(
    x0 = coords2[gr$edge[[1L]][i],"x"],
    y0 = coords2[gr$edge[[1L]][i],"y"],
    x1 = coords2[gr$edge[[2L]][i],"x"],
    y1 = coords2[gr$edge[[2L]][i],"y"],
    length = 0.08
  )
