## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Molecular Evolution Model (Markov Process) **
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
#' @name molEvolMod
#' @aliases read.mutationMat
#' 
#' @title Molecular Evolution Simulator
#' 
#' @description Functions to simulate the evolution of DNA sequences following a
#' Markov process.
#' 
#' @param par A list containing the parameters needed to construct the model's
#' transition intensity matrix.
#' @param I A 5 x 5 transition intensity matrix such as the one obtained from
#' function \code{read.mutationMat} (see details).
#' @param step The simulation time step (arbitrary units of time).
#' @param rho An evolution rate factor to be used on top of the transition
#' intensity matrix.
#' 
#' @returns
#' \describe{
#'   \item{ \code{read.mutationMat} }{...}
#'   \item{ \code{molEvolMod} }{...}
#' }
#' 
#' @details The molecular evolution model is based on a set of locations that
#' can take one of five states, namely gap ('-'), adenine ('A'), thymine ('T')
#' cytosine ('C'), or guanine ('G'). A seed sequence is evolved one location at
#' a time. Changes from one nucleotides to another appear as nucleotide
#' transitions or transversions, whereas changes from a gap to one of the
#' nucleotides appear as an insertion and changes from one the the nucleotides
#' to a gap as a deletion.
#' 
#' The changes are simulated as a simple Markov process, using a transition
#' probability matrix, which is calculated using a transition intensity matrix.
#' The off-diagonal elements of this matrix can be provided as a list from a
#' YAML configuration file to function \code{read.mutationMat} (see examples
#' below). The diagonal elements are defined by construct (i.e., the negative of
#' the sum of each off-diagonal element's row; user-provided transition
#' intensity matrices must have this property). Multiple transition intensity
#' matrices can be employed for various simulated sequences, it necessary.
#' 
#' Gap-only locations are discarded at the outset of the
#' process, yielding sets of sequences that are more or less shorter than the
#' prescribed number of nucleotides depending on the gap frequency used in the
#' initial sequence. The initial (root) sequence can be drawn from a uniform
#' distribution with user-defined frequencies (i.e., using function
#' \code{\link[PhyloSim]{drawDNASequence}}), whereas each location's evolution
#' rate can be drawn from a gamma distribution (i.e., using function
#' \code{\link[PhyloSim]{drawEvolRate}}).
#' 
#' The molecular evolution model can be instantiated multiple times, by calling
#' \code{molEvolMod} multiple times and storing the results into a list. For
#' instance, if every single location is assigned its own evolution rate, a
#' model is implemented for every single location. Then, member function
#' \code{$evolve(N)} is called to evolve a location 'N', with the returned value
#' being the new location value. When the time step of evolution rate change,
#' member function \code{$recalculate(step, rho)} is called to update the
#' transition probability matrix with time step 'step' and evolution rate 'rho'.
#' The transition probability matrix is obtained using member function
#' \code{$getMt()}.
#' 
#' @author \packageAuthor{PhyloSim}
#' 
#' @importFrom yaml yaml.load_file
#' 
NULL
#' 
#' @rdname molEvolMod
#' 
#' @examples ## Load the example of a configuration file provided with the
#' ## package:
#' yaml.load_file(
#'   system.file(
#'     package = "PhyloSim",
#'     "extdata",
#'     "evolMod.yml"
#'   )
#' ) -> dnaParam
#' 
#' ## The configuration list is located in member `$DNA`
#' read.mutationMat(
#'   dnaParam$DNA 
#' ) -> I
#' 
#' ## The transition intensity matrix:
#' I
#' 
#' @export
read.mutationMat <- function(par) {
  
  matrix(
    NA, 5L, 5L,
    dimnames=list(
      c("-","A","C","G","T"),
      c("-","A","C","G","T")
    )
  ) -> I
  
  I["-",names(par$insertion)] <- unlist(par$insertion)
  I[names(par$deletion),"-"] <- unlist(par$deletion)
  
  unlist(par$transition)[c("AG","GA","CT","TC")] ->
    I[rbind(c("A","G"),c("G","A"),c("C","T"),c("T","C"))]
  
  for(i in names(par$transversion))
    I[i,names(par$transversion[[i]])] <- unlist(par$transversion[[i]])
  
  diag(I) <- -rowSums(I, na.rm = TRUE)
  
  I
}
#' 
#' @rdname molEvolMod
#' 
#' @examples ## Implement the molecular evolution model for a single nucleotide:
#' molEvolMod(I, 1, 1) -> em1
#' 
#' ## Get the transition probability matrix as follows:
#' em1$getMt()
#' 
#' ## A vector of raw as examples of initial traits:
#' tr <- charToRaw("-ACGT")
#' 
#' ## Simulate molecular evolution from:
#' rawToChar(em1$evolve(tr[1L]))    ## a gap.
#' rawToChar(em1$evolve(tr[2L]))    ## an adenine base.
#' rawToChar(em1$evolve(tr[3L]))    ## a cytosine base.
#' rawToChar(em1$evolve(tr[4L]))    ## a guanine base.
#' rawToChar(em1$evolve(tr[5L]))    ## a thymine base.
#' 
#' ## Recalculate the probabilities for a lower mean evolution rate (one tenth
#' ## the previous one):
#' em1$recalculate(1, 0.1)
#' 
#' em1$getMt()        ## The recalculated transition probability matrix.
#' 
#' ## Simulate molecular evolution from:
#' rawToChar(em1$evolve(tr[1L]))    ## a gap.
#' rawToChar(em1$evolve(tr[2L]))    ## an adenine base.
#' rawToChar(em1$evolve(tr[3L]))    ## a cytosine base.
#' rawToChar(em1$evolve(tr[4L]))    ## a guanine base.
#' rawToChar(em1$evolve(tr[5L]))    ## a thymine base.
#' 
#' ## Base changes are now less probable.
#' 
#' ## Clean up:
#' rm(em1, tr)
#' 
#' @export
molEvolMod <- function(I, step, rho) {
  expMat <- function(x) {
    eig <- eigen(x)
    y <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
    dimnames(y) <- dimnames(x)
    y
  }
  Mt <- expMat(step*rho*I)
  nbytes <- charToRaw(paste(rownames(I),collapse=""))
  list(
    recalculate = function(step, rho)
      Mt <<- expMat(step*rho*I),
    evolve = function(N)
      sample(x=nbytes, size=1L, prob=Mt[which(N == nbytes),]),
    getMt = function() Mt
  )
}
#' 
