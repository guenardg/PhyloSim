## **************************************************************************
##
##    (c) 2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Distance vector indexing function **
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
#' Distance Vector Indexing
#' 
#' Function to calculate the indices of a distance vector that contains the
#' pairwise distance between two objects.
#' 
#' @param n An integer giving the number of objects from which the distances
#' were calculated.
#' @param i An integer giving the index of the first object.
#' @param j An integer giving the index of the second object (can be omitted,
#' see details).
#' 
#' @returns The index of the distances in the vector of distance.
#' 
#' @details A distance vector contains the lower triangle of a pairwise distance
#' matrix concatenated in a column by column manner. While storing distance in
#' that way saves substantial space (i.e., by avoiding the repeated upper
#' triangle values and the all-zero diagonal values), accessing the values is
#' somewhat more tedious than for a square matrix, for which the R language
#' provide extraction and replacement operators.
#' 
#' When argument \code{j} is omitted, all the indices between object \code{i}
#' and the other observations are returned. In that case, argument \code{i}
#' accepts only a single value.
#' 
#' @author \packageAuthor{PhyloSim}
#' 
#' @examples
#' ### Here...
#' 
#' 
#' 
#' 
#' @export
dst_idx <- function(n, i, j) {
  if(missing(j)) {
    if(length(i) != 1L)
      stop("Only a single object can be selected when argument 'j' is omitted.")
    j <- (1L:n)[-i]
  }
  .Call("PhyloSim_dstIdx", PACKAGE="PhyloSim", n, i, j)
}
#'
