## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Simulation of DNA sequences evolution - general network case **
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
#' @name simulate_network
#' 
#' @title Simulation of DNA Sequences Evolution Along the Edge of a Phylogenetic
#' Network
#' 
#' @description A function to simulate a set of DNA sequences along the edge of
#' a random phylogenetic network.
#' 
#' @param I A 5 x 5 transition intensity matrix such as the one obtained from
#' function \code{read.mutationMat}.
#' @param NS The number of sequences to generate.
#' @param NN The number of nucleotides for each of the sequence.
#' @param SID DNA segments identifiers; a vector (of any type) whose elements
#' represent the different segments whose nucleotides segregate together during
#' hybridization events.
#' @param traitMod Trait evolution parameters used to generate trait values
#' alongside the sequences (optional, no trait is generated when missing).
#' @param initState If argument \code{traitMod} is given, the initial state of
#' traits with multiple optima is specified by this argument (mandatory when
#' argument \code{traitMod} is given, see details).
#' @param initValue If argument \code{traitMod} is given, the initial state of
#' traits with multiple optima is specified by this argument (see details).
#' @param NC A function (with no argument) returning the number of descendents
#' (children) for each new node created (default: a function returning the value
#' 3).
#' @param NP A function (with no argument) returning the number of ascendants
#' (parents) for each new node created (default: a function returning the value
#' 2, see details).
#' @param gamma.shape Shape parameter of the gamma distribution from which the
#' evolution rates are drawn (see details).
#' @param gamma.scale Scale parameter of the gamma distribution from which the
#' evolution rates are drawn (see details).
#' @param timestep Time step of the simulation (default: 1 arbitrary unit of
#' time per generation).
#' @param maxDst A function (with no argument) returning the maximum distance
#' between the first parent drawn and any other parents drawn when more than one
#' parents is involved with the creation of a descendent node (default: a
#' function returning 4 times \code{timestep}).
#' @param prob A vector of length 5 giving the probabilities for drawing a gap,
#' or any one of the four bases (gap, adenine, cytosine, guanine, thymine, in
#' that order) when generating the first sequence (the ancestor of all the
#' others; default: \code{c(0.3, 0.175, 0.175, 0.175, 0.175)}).
#' @param removeGapOnly A boolean specifying whether gap-only locations should
#' by removed (default: \code{FALSE}).
#' @param verbose A boolean specifying whether to display information about the
#' links as the sequences are being created (default: \code{TRUE}).
#' 
#' @returns A \code{\link{list}} containing the following members:
#' \describe{
#' \item{trait}{ A list containing the simulation data for the trait value (see
#' details). }
#' \item{rate}{ A numeric vector of the evolution rates of each location (bases
#' or gaps). }
#' \item{em}{ A list of molecular evolution models, one for each location. }
#' \item{sqn}{ A raw matrix containing the sequences, whereby sequences are the
#' rows and locations are the columns). }
#' \item{cld}{ An integer vector giving the number of remaining descendents
#' per node at the end of the simulation process. }
#' \item{dst}{ A matrix of the distances among the sequences, in terms of the
#' numbers of time steps between them. }
#' \item{edge}{ A list of \code{NS} integer vectors giving the index of the
#' descendents for each of the nodes. }
#' \item{cont}{ A list of \code{NS} numeric vectors giving the parent's
#' relative contribution for each of their descendents. }
#' }
#' 
#' @details The function generate a random phylogenetic network with simulated
#' DNA data by evolving a set of nucleotides using Markov chain based molecular
#' evolution models as implemented using function \code{\link{molEvolMod}}.
#' Optionally, the function can simulate the evolution of a set of quantitative
#' traits using the approach implemented by function \code{\link{traitEvolMod}}.
#' 
#' The process begins by drawing a set of evolution rates from a gamma
#' distribution with shape and scale parameters provided from arguments
#' \code{gamma.shape} and \code{gamma.scale}, respectively, for each of a set of
#' locations whose number of given by argument \code{NS}. When applicable,
#' Molecular evolution models generated by function \code{\link{molEvolMod}} are
#' then implemented individually for each location from the randomly drawn rates
#' and the transition intensity matrix obtained from a YAML configuration file
#' using function \code{\link{read.mutationMat}}. Then, an initial DNA sequence
#' is generated using the probabilities (for gaps or any of the ACGT bases)
#' provided by argument \code{prob}. Also, a vector of numbers of descendents
#' per node is is initialized with zeros. Non-zero elements of this vectors
#' indicate nodes that are allowed to have descendents. With each reproduction
#' event, the value is decremented by one until the value goes back to zero.
#' Initial values of this descendant vector are generated by the function passed
#' through argument \code{NC}.
#' 
#' Once these staging steps are done, each sequence is generated by evolving
#' each location of an ancestral sequence using the molecular evolution models
#' until the prescribed number of sequences, provided by argument \code{NS}, is
#' reached. The ancestral sequence is either one of the previously generated
#' sequences, or an hybrid of two or more of these sequences. The number of
#' ascendants is generated by the function passed through argument \code{NP}.
#' When more than one ascendant node is needed, a first node is drawn randomly
#' among the ones allowed to bear descendents, then one or more node is drawn
#' randomly among the one located within a distance specified by argument
#' \code{maxDst} from the first node. When the number of ascendant nodes
#' satisfying these aforementioned conditions is smaller than the demanded
#' number, that smaller number is used instead.
#' 
#' The edges connecting each parental sequence to its immediate descendant is
#' written at each step into the edge list. Random
#' parent selection is done by randomly selecting one of the non-zero element of
#' the vector of numbers of descendents. As its name suggests this function
#' generate evolutionary sequences arranged as a tree rather than sequences
#' evolving linearly or along networks.
#' 
#' @author \packageAuthor{PhyloSim}
#' 
#' @import magrittr
#' @importFrom stats rgamma
#' 
#' @examples ## Reticulated sequence simulation
#' 
#' ## Load the example of a configuration file provided with the
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
#' ## Setting seed:
#' set.seed(12345)
#' 
#' ## Generate the sequences:
#' simulate_network(
#'   I = I,
#'   NS = 300,
#'   NN = 150,
#'   SID = rep(c("S1","S2","S3"), each=50),
#'   NC = function() 1 + rpois(1, 1.5),
#'   NP = function() 1 + rpois(1, 0.5),
#'   removeGapOnly = TRUE,
#'   verbose = FALSE
#' ) -> res
#' 
#' ## Get the concatenated sequences
#' original <- concatenate(res$sqn)
#' nogaps <-  concatenate(res$sqn, discard = "-")
#' 
#' ## To save the sequences:
#' ## write.fasta(original, "Simulated sequences - aligned.fst",
#' ##             linebreak = 70)
#' ## write.fasta(nogaps, "Simulated sequences - raw.fst", linebreak = 70)
#' 
#' ## Function for showing sequence:
#' show.sequences <- function(
#'     x, xlim, ylim, xlab = "Position", text = FALSE,
#'     col=c(A="red",`T`="blue",C="green",G="yellow", `-`="grey"), ...) {
#'   
#'   par(mar=c(4,7,1,1))
#'   
#'   if(missing(xlim))
#'     xlim <- c(0,10*ceiling(max(sapply(x,nchar))/10))
#'   
#'   if(missing(ylim))
#'     ylim <- c(length(x), 0)
#'   
#'   dev.hold()
#'   
#'   plot(NA, xlim=xlim, ylim=ylim, axes=FALSE, xlab=xlab, ylab="", ...)
#'   
#'   axis(1)
#'   
#'   axis(2, at=1:length(x), labels = names(x), las=1)
#'   
#'   for(i in 1:length(x)) {
#'     cc <- unlist(strsplit(x[i],""))
#'     for(j in 1:length(cc)) {
#'       rect(j - 1, i - 0.5, j, i + 0.5, col = col[cc[j]], border = col[cc[j]])
#'       if(text) text(j - 0.5, i, cc[j], ...)
#'     }
#'   }
#'   
#'   dev.flush()
#'   
#'   invisible(NULL) 
#'   
#' }
#' 
#' ## Showing the resulting sequences:
#' show.sequences(original, cex=0.75, font=2)   ## Natively aligned
#' show.sequences(nogaps, cex=0.75, font=2)     ## Raw
#' 
#' ## Sequence simulation with traits
#' 
#' ## Resetting RNG seed:
#' set.seed(12345)
#' 
#' ## The trait evolution model parameters:
#' dnaParam$trait %>%
#'   lapply(read.TraitEvolMod) -> traitMod
#' 
#' ## Generate the sequences with the traits:
#' simulate_network(
#'   I = I,
#'   NS = 300,
#'   NN = 150,
#'   SID = rep(c("S1","S2","S3"), each=50),
#'   traitMod = traitMod,
#'   initState = c(2,NA,1,1),
#'   initValue = c(50,0,15,-25),
#'   NC = function() 1 + rpois(1, 1.5),
#'   NP = function() 1 + rpois(1, 0.5),
#'   removeGapOnly = TRUE,
#'   verbose = FALSE
#' ) -> res2
#' 
#' ## These are the trait evolution models that were used:
#' res2$trait$tem
#' 
#' ## Obtain the nodes that are tips:
#' res2$edge %>%
#'   lapply(function(x) !length(x)) %>%
#'   unlist %>%
#'   which -> tips
#' 
#' getDescTrail <- function(net, desc) {
#'   trail <- desc
#'   while(TRUE) {
#'     net$edge %>%
#'       lapply(function(x, y) y %in% x, y = desc) %>%
#'       unlist %>%
#'       which -> asc
#'     if(!length(asc)) break
#'     if(length(asc) > 1L)
#'       sapply(
#'         asc,
#'         function(x) net$cont[[x]][net$edge[[x]] == desc]
#'       ) %>%
#'       which.min %>%
#'       asc[.] -> asc
#'     trail <- c(asc, trail)
#'     desc <- asc
#'   }
#'   trail
#' }
#' 
#' ## Get the concatenated sequences
#' original2 <- concatenate(res2$sqn)
#' 
#' ## Get the trail from the first tip:
#' show.sequences(
#'   original2[getDescTrail(res2, tips[1])],
#'   cex=0.75,
#'   font=2L
#' )
#' 
#' ## Get the trail from the tenth tip:
#' show.sequences(
#'   original2[getDescTrail(res2, tips[10])],
#'   cex=0.75,
#'   font=2L
#' )
#' 
#' ## Plot the trait simulation results for the first trail:
#' i <- getDescTrail(res2, tips[1])
#' 
#' par(mar=c(4,4,1,1))
#' plot(NA, xlim=c(0,(length(i) - 1)), ylim=range(res2$trait$traitLog[i,]),
#'      xlab="Time", ylab="Trait value", las=1)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,1], col="black")
#' if(!is.na(res2$trait$tem[[1]]$getOptima())[1])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,1], col="black",
#'         lty=3)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,2], col="red")
#' if(!is.na(res2$trait$tem[[2]]$getOptima())[1L])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,2], col="red", lty=3)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,3], col="blue")
#' if(!is.na(res2$trait$tem[[3]]$getOptima())[1])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,3], col="blue", lty=3)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,4], col="green")
#' if(!is.na(res2$trait$tem[[4]]$getOptima())[1])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,4], col="green",
#'         lty=3)
#' 
#' ## Plot the trait simulation results for the tenth trail:
#' i <- getDescTrail(res2, tips[10])
#' 
#' par(mar=c(4,4,1,1))
#' plot(NA, xlim=c(0,(length(i) - 1)), ylim=range(res2$trait$traitLog[i,]),
#'      xlab="Time", ylab="Trait value", las=1)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,1], col="black")
#' if(!is.na(res2$trait$tem[[1]]$getOptima())[1])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,1], col="black",
#'         lty=3)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,2], col="red")
#' if(!is.na(res2$trait$tem[[2]]$getOptima())[1])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,2], col="red", lty=3)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,3], col="blue")
#' if(!is.na(res2$trait$tem[[3]]$getOptima())[1])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,3], col="blue", lty=3)
#' lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,4], col="green")
#' if(!is.na(res2$trait$tem[[4]]$getOptima())[1])
#'   lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,4], col="green",
#'         lty=3)
#' 
#' ## Clean up:
#' rm(dnaParam,I,traitMod,res,res2,nogaps,original,original2,show.sequences,
#'    getDescTrail,i,tips)
#' 
#' @export
simulate_network <- function(I, NS, NN, SID, traitMod, initState, initValue,
                             NC = function() 3L, NP = function() 2L,
                             gamma.shape = 5, gamma.scale = 5e-4, timestep = 1,
                             maxDst = function() 4*timestep,
                             prob = c(0.3,rep(0.175,4L)), removeGapOnly = FALSE,
                             verbose = TRUE) {
  
  ## Check the number of segment identifiers:
  if(!missing(SID)) {
    if(length(SID) != NN) {
      warning("Number of segment identifiers (argument 'SID': ", length(SID),
              ") does not match the value of argument argument 'NN' (", NN,
              "); the segment identifiers had to be adjusted.")
      SID <- rep(SID, length.out=NN)
    }
  } else
    SID <- rep("A",NN)
  
  ## Nucleotide segments:
  nseg <- unique(SID)
  
  ## Result list:
  out <- list()
  
  ## Should we also simulate trait evolution alongside the sequences?
  simTrait <- !missing(traitMod)
  
  ## Create the trait evolution models and the log files:
  if(simTrait) {
    
    out$trait <- list(tem = lapply(traitMod, traitEvolMod))
    trNms <- lapply(out$trait$tem, function(x) x$getName())
    NT <- length(trNms)
    out$trait$stateLog <- matrix(NA, NS, NT, dimnames = list(NULL, trNms))
    out$trait$optimLog <- matrix(NA, NS, NT, dimnames = list(NULL, trNms))
    out$trait$traitLog <- matrix(NA, NS, NT, dimnames = list(NULL, trNms))
    
    if(!missing(initState)) {
      out$trait$stateLog[1L,] <- initState
    } else
      stop("Missing initial trait state(s)!")
    
    out$trait$tem %>%
      lapply(
        function(x)
          x$getOptima()[out$trait$stateLog[1L,x$getName()]]
      ) %>%
      unlist -> out$trait$optimLog[1L,]
    
    if(!missing(initValue)) {
      out$trait$traitLog[1L,] <- initValue
    } else {
      warning("Initial trait values missing! Assumed to be 0")
      out$trait$traitLog[1L,] <- 0
    }
    
  }
  
  ## Draw an evolution rate for each locus:
  rgamma(
    NN,
    shape = gamma.shape,
    scale = gamma.scale
  ) -> out$rate
  
  ## Create the nucleotide evolution models from the mutation intensity matrix
  ## I a the specified time step, and the evolution rates drawn previously:
  out$em <- list()
  for(i in 1L:NN)
    out$em[[i]] <- molEvolMod(I, timestep, out$rate[i])
  
  ## A matrix to store NS sequences of NN loci each:
  out$sqn <- matrix(raw(), NS, NN)
  
  ## Vector of the numbers of remaining children per node:
  out$cld <- integer(NS)
  
  ## Generate an initial sequence. The probabilities of drawing a gap or any of
  ## the four nucleotides is given by argument 'prob'. The role of the gaps is
  ## to allow the simulation of nucleotide insertion or deletion. A gap turning
  ## into a nucleotide simulate an insertion and a nucleotide turning into a gap
  ## simulates a deletion.
  sample(
    x = charToRaw("-ACGT"),
    size = NN,
    prob = prob,
    replace = TRUE
  ) -> out$sqn[1L,]
  
  ## Number of children for the initial node:
  out$cld[1L] <- NC()
  
  ## Distance matrix:
  out$dst <- rep(NA,NS*(NS - 1L)/2)
  
  ## Edge and parents contributions lists:
  out$edge <- list()
  out$cont <- list()
  
  ## Initialization of the lists:
  for(i in 1L:NS) {
    out$edge[[i]] <- integer()
    out$cont[[i]] <- numeric()
  }
  
  ## Sequence generation loops:
  for(i in 2L:NS) {
    ## 'i' is the new child node to be created
    
    ## If we are at the first step:
    if(i == 2L) {
      ## For the first step, only the initial node is available as a parent:
      p <- 1L
    } else {
      ## For the remaining steps, 
      np <- NP()
      
      ## 'wh' is the indices of the nodes that are allowed to be parent.
      wh <- which(out$cld > 0L)
      
      ## Selects a maximum of NP parent(s) among the nodes that are allowed to
      ## reproduce:
      p <- sample(wh, 1L)
      wh <- wh[!(wh %in% p)]
      
      ## Here, the potential parents must be checked for sufficient similarity
      wh <- wh[out$dst[dst_idx(NS, p, wh)] <= maxDst()]
      
      while(length(wh) && (length(p) < np)) {
        
        if(length(wh) == 1L) {
          p <- c(p, wh)
          wh <- integer()
        } else {
          ppp <- sample(wh, 1L)
          p <- c(p, ppp)
          wh <- wh[!(wh %in% ppp)]
        }
        
      }
      
    }
    
    ## 1. Generate the new sequence from the parental sequences:
    if(length(p) > 1L) {
      
      ssel <- sample(p, length(nseg), TRUE)
      
      sel <- integer(NN)
      
      for(j in 1L:length(nseg))
        sel[SID == nseg[j]] <- ssel[j]
      
      for(j in 1L:NN)
        out$sqn[i,j] <- out$em[[j]]$evolve(out$sqn[sel[j],j])
      
      pp <- sapply(p, function(x,sel) mean(sel == x), sel=sel)
      
    } else {
      
      for(j in 1L:NN)
        out$sqn[i,j] <- out$em[[j]]$evolve(out$sqn[p,j])
      
      pp <- 1
    }
    ## paste(sapply(out$sqn[i,], rawToChar), collapse="")
    
    ## 
    ## 2. Insert the time in the locations:
    out$dst[dst_idx(NS, i, p)] <- timestep
    
    ## 3. Initialize the parental status of the newly created node:
    out$cld[i] <- NC()
    
    ## 4. Decrease the parental status of the parent node:
    out$cld[p] <- out$cld[p] - 1L
    
    ## 5. Add the new edge to the edge list and updating the distances:
    for(j in 1L:length(p)) {
      out$edge[[p[j]]] %<>% c(i)
      out$cont[[p[j]]] %<>% c(timestep*pp[j])
    }
    
    ## 6. Update the distances:
    others <- (1L:i)[-c(i,p)]
    
    if(length(others)) {
      dst <- pp[1L]*out$dst[dst_idx(NS, p[1L], others)]
      if(length(p) > 1L)
        for(j in 2L:length(p))
          dst <- dst + pp[j]*out$dst[dst_idx(NS, p[j], others)]
      out$dst[dst_idx(NS, i, others)] <- dst + timestep
    }
    
    ## 7. Evolve each traits, if required
    if(simTrait) {
      
      out$trait$tem %>%
        lapply(
          function(x, p) {
            if(length(p) > 1L) {
              x$updateState(
                sample(out$trait$stateLog[p,x$getName()], 1L, prob = pp)
              )
            } else
              x$updateState(out$trait$stateLog[p,x$getName()])
          },
          p = p
        ) %>%
        unlist -> out$trait$stateLog[i,]
      
      out$trait$tem %>%
        lapply(
          function(x)
            x$getOptima()[out$trait$stateLog[i,x$getName()]]
        ) %>%
        unlist -> out$trait$optimLog[i,]
      
      out$trait$tem %>%
        lapply(
          function(x, p)
            if(length(p) > 1L) {
              x$updateValue(
                state = out$trait$stateLog[i,x$getName()],
                value = sum(pp*out$trait$traitLog[p,x$getName()])
              )
            } else
              x$updateValue(
                state = out$trait$stateLog[i,x$getName()],
                value = out$trait$traitLog[p,x$getName()]
              ),
          p = p
        ) %>%
        unlist -> out$trait$traitLog[i,]
      
    }
    
    if(verbose) {
      cat(sprintf("Link: %04d -> %04d (%f)\n",p,i,pp),sep="")
      if(simTrait) {
        cat("Trait(s): ",
            paste(sprintf("%0.3f",out$trait$traitLog[i,]),collapse=", "),"\n")
      } else
        cat("\n")
    }
    
  }
  
  ## Transform the dissimilarity vector into a dist-class object:
  list(
    Size = NS,
    Diag = FALSE,
    Upper = FALSE,
    method = "patristic",
    call = match.call(),
    class = "dist"
  ) -> attributes(out$dst)
  
  ## Purge gap-only locations:
  if(removeGapOnly)
    out$sqn %<>% .[,apply(., 2L, function(x) !all(x == charToRaw("-")))]
  
  ## Give a name to the sequences
  rownames(out$sqn) <- sprintf("sequence%04d", 1L:NS)
  
  out
  
}
#'
