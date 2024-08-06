##
### Pre-created networks
##
## rm(list=ls())

library(magrittr)

NS <- 1000L

net <- list(edge = list(), weight = list())
lnk <- integer(NS)
asc <- integer(NS)
dsc <- integer(NS)

for(i in 1L:NS) {
  net$edge[[i]] <- integer()
  net$weight[[i]] <- integer()
}

lnk[1L] <- 1L + rpois(1L, 3)

for(i in 2L:NS) {
  wh <- which(lnk > 0L)
  np <- 1L + rpois(1L,2.5)
  p <- if(length(wh) <= np) wh else sample(wh, np)
  for(j in p) {
    net$edge[[j]] %<>% c(i)
    net$weight[[j]] %<>% c(rbeta(1L,1,3))
    dsc[j] <- dsc[j] + 1L
    lnk[j] <- lnk[j] - 1L
  }
  asc[i] <- length(p)
  lnk[i] <- 1L + rpois(1L, 3)
}

## all(unlist(lapply(net$edge,length)) == dsc)
## all(c(0,unname(tapply(unlist(net$edge),unlist(net$edge),length))) == asc)

con <- logical(NS)
for(i in 1L:NS) {
  if(length(net$edge[[i]])) {
    con[i] <- TRUE
  } else {
    for(j in (1L:NS)[-i])
      if(any(net$edge[[j]] == i)) {
        con[i] <- TRUE
        break
      }
  }
}
all(con)

## unlist(lapply(net$edge,length))

w <- matrix(0, NS, NS)
for(i in 1L:(NS - 1L))
  for(j in (i + 1L):NS) {
    if(any(net$edge[[i]] == j)) {
      w[i,j] <- w[j,i] <- net$weight[[i]][net$edge[[i]] == j]
    } else if(any(net$edge[[j]] == i)) {
      w[i,j] <- w[j,i] <- net$weight[[j]][net$edge[[j]] == i]
    }
  }

d <- sqrt(1 - w)
diag(d) <- 0

g <- -0.5*d^2

delta <- t(t(g - rowMeans(g)) - colMeans(g)) + mean(g)
## rowMeans(delta)
## colMeans(delta)
## mean(delta)

eig <- eigen(delta)

coords <- eig$vector[,1L:2L] %*% diag(sqrt(eig$values[1L:2L]))

plot(coords)

for(i in 1L:length(net$edge))
  if(length(net$edge[[i]])) {
    arrows(
      x0 = coords[i,1L],
      x1 = coords[net$edge[[i]],1L],
      y0 = coords[i,2L],
      y1 = coords[net$edge[[i]],2L],
      length = 0.05,
      lty = 3L,
      lwd = 0.25
    )
    ## if(is.null(locator(1L))) break
  }

## which(unlist(lapply(net$edge,length)) == 0)[2L] -> n

n <- 200L
wh <- which(unlist(lapply(net$edge,function(x, y) any(x == y), y=n)))

if(length(wh)) {
  segments(x0=coords[wh,1L], x1=coords[n,1L], y0=coords[wh,2L],
           y1=coords[n,2L], col = "green")
  points(x=coords[wh,1L], y=coords[wh,2L], pch=21L, bg="green")
}
if(length(net$edge[[n]])) {
  segments(x0=coords[n,1L], x1=coords[net$edge[[n]],1L], y0=coords[n,2L],
           y1=coords[net$edge[[n]],2L], col = "red")
  points(x=coords[net$edge[[n]],1L], y=coords[net$edge[[n]],2L], pch=21L,
         bg="red")
}
points(x=coords[n,1L], y=coords[n,2L], pch=21L, bg="blue")

NN <- 1000L
NSEG <- 5L




diff(c(1L, 1L + sort(sample(NN - 2L, NSEG - 1L)), NN))


sample(2, 1L)
