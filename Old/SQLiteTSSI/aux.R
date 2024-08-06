
PCplot <- function(net, ...) {
  
  dbGetQuery(
    net$con,
    "SELECT * FROM diss"
  ) -> dd
  
  NS <- max(dd$j)
  
  d <- matrix(0, NS, NS)
  
  for(k in 1L:nrow(dd))
    d[dd$i[k],dd$j[k]] <- d[dd$j[k],dd$i[k]] <- dd$d[k]
  
  G <- -0.5*d^2
  delta <- t(t(G - rowMeans(G)) - colMeans(G)) + mean(G)
  ## rowMeans(delta)
  ## colMeans(delta)
  ## mean(delta)
  
  eig <- eigen(delta)
  coords <- eig$vectors[,1L:2L] %*% diag(sqrt(eig$values[1L:2L]))
  
  plot(coords, ...)
  
  dbGetQuery(
    net$con,
    "SELECT i, j FROM edge"
  ) -> edge
  
  for(k in 1L:nrow(edge))
    segments(
      x0 = coords[edge$i[k],1L],
      x1 = coords[edge$j[k],1L],
      y0 = coords[edge$i[k],2L],
      y1 = coords[edge$j[k],2L],
      col="grey",
      lwd = 0.5,
      ...
    )
  
  invisible(NULL)
  
}

PCplot3d <- function(net, ...) {
  
  dbGetQuery(
    net$con,
    "SELECT * FROM diss"
  ) -> dd
  
  NS <- max(dd$j)
  
  d <- matrix(0, NS, NS)
  
  for(k in 1L:nrow(dd))
    d[dd$i[k],dd$j[k]] <- d[dd$j[k],dd$i[k]] <- dd$d[k]
  
  G <- -0.5*d^2
  delta <- t(t(G - rowMeans(G)) - colMeans(G)) + mean(G)
  
  eig <- eigen(delta)
  coords <- eig$vectors[,1L:3L] %*% diag(sqrt(eig$values[1L:3L]))
  
  dbGetQuery(
    net$con,
    "SELECT i, j FROM edge"
  ) -> edge
  
  dbGetQuery(
    net$con,
    "SELECT cld FROM node"
  ) -> node
  
  wire3d(
    mesh3d(
      x = coords,
      segments = rbind(edge$i, edge$j)
    ),
    col="grey"#,
    #...
  )
  
  points3d(coords[1L,], col="orange")#, ...)
  points3d(coords[node$cld == 0L,], col="red")#, ...)
  points3d(coords[node$cld > 0L,], col="green")#, ...)
  
  invisible(NULL)
  
}
