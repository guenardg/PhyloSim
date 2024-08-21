## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Trait Evolution Simulation - SQLite TSSI Utilities **
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
#' 
#' SQLite Trait Simulation Utilities
#' 
#' Utility functions for the trait and sequence simulation interface, a database
#' interface for simulating traits along evolutionary time.
#' 
#' @name TSSI-utils
#' 
#' @aliases getLeaves
#' 
#' @param net A two elements list: the name of the network and an SQLite
#' database connection from function \code{\link[DBI]{dbConnect}}.
#' @param v An integer vector: indices of the vertices for which to
#' determine the lineage (or most direct lineage for a reticulated graph).
#' @param leaf A logical: whether or not to declare the leaves of the graph
#' (i.e., the vertices that have no descendants) as species in the graph.
#' 
#' @details The lineages returned by function \code{getLineage} are unique in
#' the case of a phylogenetic tree (including linear trees). For a reticulated
#' network, however, there is oftentimes several possible lineages since any
#' vertex may have multiple ancestor vertices. Whenever necessary, function
#' \code{getLineage} works around that matter by selecting, among the multiple
#' ancestors of a vertex, that which is the closest to the root of the graph.
#' 
#' @return
#' \describe{
#'   \item{getLeaves}{The indices of the vertices having no descendents.}
#'   \item{getLineage}{A list of integer vectors, each of which contains the
#'   indices of the vertices of its (most direct) lineage.}
#'   \item{getPatristic}{An object of class "dist" as obtained from function
#'   \code{\link[stats]{dist}} containing patristic distance values. The
#'   patristic distance is expressed in term of the simulation time step.}
#'   \item{getGraph}{A \code{\link[MPSEM]{graph-class}} object.}
#' }
#' 
#' @author \packageAuthor{PhyloSim}
#' 
#' @examples ## Loading the example graph (reticulated; N = 100) as follows:
#' openNetwork(
#' filename = system.file(
#'   package = "PhyloSim",
#'   "extdata",
#'   "reticulated-ex2.sim"
#' ),
#' load = TRUE
#' ) -> net2
#' 
#' ## Which vertex is a leaf is determined as follows:
#' leaf <- getLeaves(net2)
#' leaf
#' 
#' ## The pairwise average patristic distances are obtained as follows:
#' dst <- getPatristic(net2)
#' 
#' ## The structure of the "dist" object can be viewed as follows:
#' str(dst)
#' 
#' ## The first elements of the contained vector:
#' head(dst)
#' 
#' ## There is no convenient slicing operator for "dist" objects, one can
#' ## display a slice limited to the first six objects as follows:
#' as.dist(as.matrix(dst)[1:6,1:6])
#' 
#' ## The phylogenetic graph is extracted as follows:
#' gr <- getGraph(net2)
#' 
#' ## The structure of the graph-class object can be viewed as follows:
#' str(gr)
#' 
#' ## Note: function from package MPSEM is needed to further process graph-class
#' ## objects. For instance, it can be used to calculate phylogenetic
#' ## eigenvector maps, which can, in turn, be used to predict trait values.
#' 
NULL
#' 
#' @describeIn TSSI-utils
#' 
#' Get Network Leaves
#' 
#' Extract the leaves (i.e., the vertices having no descendents) from a network.
#' 
#' @export
getLeaves <- function(net) {
  
  dbGetQuery(
    net$con,
    "SELECT i,j FROM edge ORDER BY rowid"
  ) -> edge
  
  unique(edge$j[!(edge$j %in% unique(edge$i))])
}
#' 
#' @describeIn TSSI-utils
#' 
#' Get Vertex Lineage
#' 
#' Get the direct lineage of a vertex (node of leaf) in a phylogenetic tree, or
#' the most direct lineage of a vertex in a phylogenetic network.
#' 
#' @export
getLineage <- function(net, v) {
  
  dbGetQuery(
    net$con,
    "SELECT i,j FROM edge ORDER BY rowid"
  ) -> e
  
  lin <- list()
  
  for(i in 1L:length(v)) {
    
    x <- v[i]
    
    while(length(nd <- e$i[which(e$j == tail(x,1))])) {
      
      if(length(nd) > 1L) {
        
        dbGetQuery(
          net$con,
          sprintf(
            "SELECT j, d FROM diss WHERE i = 1 AND j IN(%s)",
            paste(nd, collapse=",")
          )
        ) -> d
        
        nd <- d$j[which.min(d$d)]
      }
      
      x <- c(x, nd)
    }
    
    lin[[i]] <- unname(rev(x))
    
  }
  
  lin
}
#' 
#' @describeIn TSSI-utils
#' 
#' Get Dissimilarity
#' 
#' Obtain the pairwise dissimilarities (linear, tree) for average pairwise
#' dissimilarities (reticulated) between the vertices of a phylogenetic graph.
#' 
#' @export
getPatristic <- function(net) {
  
  unname(
    unlist(
      dbGetQuery(
        net$con,
        "SELECT COUNT(*) FROM node"
      )
    )
  ) -> N
  
  structure(
    numeric(N*(N - 1L)/2),
    Size = N,
    Diag = FALSE,
    Upper = FALSE,
    method = "patristic",
    class = "dist",
    call = match.call()
  ) -> dst
  
  dbGetQuery(
    net$con,
    "SELECT * FROM diss"
  ) -> d
  
  dst[dst_idx(N, d$i, d$j)] <- d$d
  
  dst
  
}
#' 
#' @describeIn TSSI-utils
#' 
#' Get Graph
#' 
#' Obtain a graph object from a simulated phylogenetic graph.
#' 
#' @importFrom MPSEM pop.graph add.edge
#' 
#' @export
getGraph <- function(net, leaf = TRUE) {
  
  unname(
    unlist(
      dbGetQuery(
        net$con,
        "SELECT COUNT(*) FROM node"
      )
    )
  ) -> N
  
  as.matrix(
    dbGetQuery(
      net$con,
      "SELECT i, j FROM edge"
    )
  ) -> e
  
  as.matrix(
    dbGetQuery(
      net$con,
      sprintf(
        "SELECT * FROM diss
     WHERE %s",
        paste(
          sprintf("(i = %d AND j = %d)", e[,1L], e[,2L]),
          collapse = " OR "
        )
      )
    )
  ) -> d
  
  apply(
    e,
    1L,
    function(x)
      d[(x[1L] == d[,1L]) & (x[2L] == d[,2L]), 3L]
  ) -> dd
  
  if(leaf) {
    
    sp <- logical(N)
    sp[unique(e[!(e[,2L] %in% unique(e[,1L])),2L])] <- TRUE
    
    pop.graph(
      n = N,
      vertex = list(species = sp),
      label = sprintf("V%d",1L:N)
    ) -> gr
    
  } else
    pop.graph(
      n = N,
      label = sprintf("V%d",1L:N)
    ) -> gr
  
  add.edge(
    gr,
    from = e[,1L],
    to = e[,2L],
    edge = list(distance = dd),
    label = sprintf("E%d",1L:NROW(e))
  ) -> gr
  
  gr
}
#' 
