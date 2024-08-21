
## Next: functions openNetwork and closeNetwork have to be canned.

## Opening procedure:
## 
## The function take a an optional file name if no file name is provided, a
## connection is opened in the memory, the database is initialized, and the
## function returns the data base connection.
## 
## If a file name is provided, it is used to open a connection (with the
## appended .sim suffix). Then, the function verifies if the database is empty.
## If the database is empty, it initializes the database with the yaml data file
## provided with the package and returns the data base connection. If the data
## base is not empty, the function verifies whether the data table are all
## present. If they are all present, it returns the connection. If any is
## missing, it issues an error message and stops without returning anything.
## 

##
## closing procedure:
## The function verifies whether the storage is in the memory. If it is not in
## the memory, it calls dbDisconnect right away and returns NULL invisibly. If
## it is in the memory, the function verify if a file with the specified (or
## default) file name is already present. If it is present, it verify if it is
## allowed to overwrite an existing file. If not, it aborts issuing an error
## message. If it may overwrite, it deletes the existing file. Then, it creates
## an empty SQLite connection with the file name, copies the in-memory database
## into it, and close both the in-memory and the in-file data base connections.
## 

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
