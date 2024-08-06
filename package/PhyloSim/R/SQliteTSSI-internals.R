## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Trait Evolution Simulation - SQLite TSI Internal functions **
##
##    This file is part of PhyloSim
##
##    PhyloSim is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    PhyloSim is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with PhyloSim. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##

#' @importFrom DBI dbGetQuery dbExecute
#' @importFrom utils tail

addNode <- function(net, cld) {
  
  dbExecute(
    net$con,
    sprintf("INSERT INTO node (cld) VALUES (%d)",cld[1L])
  )
  
  dbGetQuery(
    net$con,
    "SELECT MAX(rowid) AS node FROM node"
  ) %>% unlist
  
}

getAsc <- function(net) {
  
  dbGetQuery(
    net$con,
    "SELECT rowid AS asc FROM node WHERE cld > 0"
  ) %>%
    unlist
  
}

addEdge <- function(net, i, j, contrib) {
  
  n <- max(length(i),length(j))
  i <- rep(i, length.out=n)
  j <- rep(j, length.out=n)
  contrib <- rep(contrib, length.out=n)
  
  dbExecute(
    net$con,
    sprintf(
      "INSERT INTO edge(i,j,contrib) VALUES %s",
      paste(sprintf("(%d,%d,%f)",i,j,contrib),collapse=",")
    )
  )
  
  invisible(NULL)
  
}

decrementChildren <- function(net, node) {
  
  dbGetQuery(
    net$con,
    sprintf(
      "SELECT rowid, cld FROM node WHERE rowid IN (%s)",
      paste(node, collapse=",")
    )
  ) -> q
  
  for(k in 1L:nrow(q))
    dbExecute(
      net$con,
      sprintf(
        "UPDATE node SET cld = %d WHERE rowid = %d",
        if(q$cld[k]) q$cld[k] - 1L else 0, q$rowid[k]
      )
    )
  
  invisible(NULL)
  
}

addDiss <- function(net, i, j, value) {
  
  n <- max(length(i),length(j))
  i <- rep(i,length.out=n)
  j <- rep(j,length.out=n)
  value <- rep(value, length.out=n)
  
  for(k in 1L:n)
    if(i[k] > j[k]) {
      ii <- i[k]
      i[k] <- j[k]
      j[k] <- ii
    }
  
  ok <- i != j
  
  dbExecute(
    net$con,
    sprintf(
      "INSERT INTO diss (i,j,d) VALUES %s",
      paste(sprintf("(%d,%d,%f)",i[ok],j[ok],value[ok]),collapse=",")
    )
  )
  
  if(sum(!ok))
    warning(
      "Attempt at changing ", sum(!ok), " diagonal value(s): \n",
      paste(sprintf("d[%d,%d] = %f",i[!ok],j[!ok],value[!ok]),collapse="\n")
    )
  
  invisible(NULL)
  
}

getDiss <- function(net, i, j) {
  
  n <- max(length(i),length(j))
  i <- rep(i, length.out=n)
  j <- rep(j, length.out=n)
  out <- numeric(n)
  
  for(k in 1L:n)
    if(i[k] > j[k]) {
      ii <- i[k]
      i[k] <- j[k]
      j[k] <- ii
    }
  
  ok <- i != j
  
  dbGetQuery(
    net$con,
    sprintf(
      "SELECT i, j, d FROM diss WHERE i IN (%s) AND j IN (%s)",
      paste(unique(i),collapse=","), paste(unique(j),collapse=",")
    )
  ) -> q
  
  for(k in which(ok)) {
    tmp <- q$d[q$i == i[k] & q$j == j[k]]
    if(!length(tmp)) {
      out[k] <- NA
      warning("Missing dissimilarity values for ",i[k]," and ",j[k])
    } else if(length(tmp) > 1L) {
      out[k] <- tail(tmp,1L)
      warning("Multiple dissimilarity values found for ",i[k]," and ",j[k])
    } else
      out[k] <- tmp
  }
  
  out
  
}
