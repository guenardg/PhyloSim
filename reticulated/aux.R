
## Plots a network with nodes arranged around a circle:
network_round_plot <- function(x, ...) {
  
  xy <- seq(0, 2*pi, length.out = length(x$cld) + 1L)
  xy <- head(xy, length(x$cld))
  data.frame(
    x=cos(xy),
    y=sin(xy),
    lab=1L:length(x$cld),
    cld=x$cld
  ) -> xy
  
  par(mar=c(1,1,1,1))
  plot(NA, xlim=1.3*c(-1,1), ylim=1.3*c(-1,1), asp=1, axes=FALSE, xlab="",
       ylab="", ...)
  points(xy, pch=21, bg="grey")
  
  apply(
    xy,
    1L,
    function(x) {
      text(x=1.1*x[1L], y=1.1*x[2L], labels=x[3L],
           srt=180*atan(x[2L]/x[1L])/pi, cex=0.8)
      invisible(NULL)
    }
  )
  
  apply(
    xy,
    1L,
    function(x) {
      text(x=1.2*x[1L], y=1.2*x[2L], labels=x[4L],
           srt=180*atan(x[2L]/x[1L])/pi, cex=0.5)
      invisible(NULL)
    }
  )
  
  for(i in 1:length(x$edge))
    if(length(x$edge[[i]]))
      arrows(
        x0 = xy[i,1L],
        x1 = xy[x$edge[[i]],1L],
        y0 = xy[i,2L],
        y1 = xy[x$edge[[i]],2L],
        length = 0.05
      )
  
  invisible(NULL)
  
}

## Function for showing sequence:
show.sequences <- function(
    x, xlim, ylim, xlab="Position", text = FALSE,
    col=c(A="red",`T`="blue",C="green",G="yellow", `-`="grey"), ...) {
  
  par(mar=c(4,7,1,1))
  
  if(missing(xlim))
    xlim <- c(0,10*ceiling(max(sapply(x,nchar))/10))
  
  if(missing(ylim))
    ylim <- c(length(x), 0)
  
  dev.hold()
  
  plot(NA, xlim=xlim, ylim=ylim, axes=FALSE, xlab=xlab, ylab="", ...)
  
  axis(1L)
  
  axis(2L, at=1:length(x), labels = names(x), las=1)
  
  for(i in 1L:length(x)) {
    cc <- unlist(strsplit(x[i],""))
    for(j in 1L:length(cc)) {
      rect(j - 1, i - 0.5, j, i + 0.5, col = col[cc[j]], border = col[cc[j]])
      if(text) text(j - 0.5, i, cc[j], ...)
    }
  }
  
  dev.flush()
  
  invisible(NULL) 
  
}

## Function to obtain the most direct evolutionary trail:
getDescTrail <- function(net, desc) {
  trail <- desc
  while(TRUE) {
    net$edge %>%
      lapply(function(x, y) y %in% x, y = desc) %>%
      unlist %>%
      which -> asc
    if(!length(asc)) break
    if(length(asc) > 1L)
      sapply(
        asc,
        function(x) net$cont[[x]][net$edge[[x]] == desc]
      ) %>%
      which.min %>%
      asc[.] -> asc
    trail <- c(asc, trail)
    desc <- asc
  }
  trail
}
