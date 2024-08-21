## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Trait Evolution Model **
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
#' @name traitEvolMod
#' @aliases read.TraitEvolMod TEM-class
#' 
#' @title Trait Evolution Simulator
#' 
#' @description Functions to simulate the evolution of traits as an
#' Ornstein-Uhlenbeck process, optionally with trait optimal value evolving
#' following a Markov process.
#' 
#' @param cfg A list containing the configuration parameters needed to construct
#' an evolution model a single quantitative trait.
#' @param par A list containing the parameters, such as the ones provided by
#' \code{read.TraitEvolMod}, that are needed to build a quantitative trait
#' evolution model (see details).
#' @param x A TEM-class object.
#' @param ... further arguments to be passed to other functions and methods
#' 
#' @return Function \code{read.TraitEvolMod} returns a list with one or more
#' members, each of which being the parameters of a trait evolution model. Each
#' trait evolution parameter set is made of the following member:
#' \describe{
#'   \item{\code{$name}}{the name of the trait,}
#'   \item{\code{$sigma}}{the trait variance parameter,}
#'   \item{\code{$step}}{the trait evolution initial step size,}
#'   \item{\code{$alpha}}{the selection rate,}
#'   \item{\code{optima}}{(may be omitted when \code{alpha = 0}) one or more the
#'   trait optima to be used when \code{alpha > 0},}
#'   \item{\code{$transition}}{(may be omitted when only one optimum is
#'   provided) a n-by-n transition intensity matrix, where n is the number of
#'   optima.}
#' }
#' Function \code{traitEvolMod} returns a TEM-class object, which consists of 8
#' member functions:
#' \describe{
#'   \item{\code{$getName}}{has no argument and returns the name of the trait,}
#'   \item{\code{$getStep}}{has no argument and returns the size of the
#'   simulation time step}
#'   \item{\code{$setStep}}{is given the time step as an argument and sets,
#'   the simulation time step, including the recalculation of the transition
#'   intensity matrix in the case of a trait with multiple optima,}
#'   \item{\code{$getOptima}}{has no argument and returns the trait optimum
#'   values,}
#'   \item{\code{getTransition}}{has no argument and returns the transition
#'   density matrix,}
#'   \item{\code{getProb}}{has no argument and returns the transition
#'   probability matrix,}
#'   \item{\code{$updateState}}{is given the index of the effective trait
#'   optimum and returns the updated index (allows shifts in trait optimum),}
#'   \item{\code{updateValue}}{is given the value of the trait and, optionally,
#'   the index of the optimum state, at returns the updated trait value,}
#'   \item{\code{dumpConfig}}{has no argument and returns the trait evolution
#'   model parameters as a list, and}
#'   \item{\code{print}}{has no argument and prints the parameters of the trait
#'   model on the screen (returns \code{NULL} invisibly).}
#' }
#' 
#' @details The trait evolution model is based on a random walk process of one
#' of categories: Weiner process and the Ornstein-Uhlenbeck process. The Weiner
#' process, also known as 'Brownian motion' is a random neutral process with
#' unbounded traits values. It is characterized using a single variance
#' parameter called \code{sigma} that dictates to what extent the trait value
#' fluctuates randomly in time.
#' 
#' The Ornstein-Uhlenbeck process is a random process with attractors (optima),
#' with values varying about these optima (or single optimum). It is notable
#' that whenever multiple optima are possible, only a single optimum is
#' effective at any given time. The Ornstein-Uhlenbeck process involve two more
#' values: the trait optimum effective at a particular time and the selection
#' rate (\code{alpha}). The optima represents the value toward which the
#' trait is attracted, whereas the value of \code{alpha} dictates to what
#' extent such an attraction occurs. When \code{alpha = 0}, the optimum stops
#' being an attractor and the Ornstein-Uhlenbeck process becomes purely a Weiner
#' process, whereas when \code{alpha} is very large, the trait value detracts
#' only slightly from its optimal value.
#' 
#' The trait simulation interface enabled to handle traits evolving neutrally
#' (pure Weiner process), and according to an Ornstein-Uhlenbeck process with
#' one or many optima. When many optima are used, transitions among the values
#' are simulated as a Markov process. For that purpose, the simulator needs to
#' be provided with a transition density matrix. 
#' 
#' 
#' 
#' @author \packageAuthor{PhyloSim}
#' 
#' @importFrom yaml yaml.load_file
#' @importFrom stats rnorm
#' 
NULL
#' 
#' @rdname traitEvolMod
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
#' ## This file contains four trait configurations.
#' 
#' ## The first trait file appears as follows:
#' 
#' dnaParam$trait[[1]]
#' 
#' ## The second, third, fourth appears as follows:
#' 
#' dnaParam$trait[[2]]
#' 
#' dnaParam$trait[[3]]
#' 
#' ## and
#' 
#' dnaParam$trait[[4]]
#' 
#' ## The configuration list is obtained from the raw file as follows:
#' 
#' cfg <- read.TraitEvolMod(dnaParam$trait[[1]])
#' cfg
#' 
#' ## All four lists are obtained as follows:
#' 
#' cfg <- lapply(dnaParam$trait, read.TraitEvolMod)
#' 
#' cfg[[1]]
#' 
#' cfg[[2]]
#' 
#' cfg[[3]]
#' 
#' cfg[[4]]
#' 
#' @export
read.TraitEvolMod <- function(cfg) {
  if(is.null(cfg$name))
    stop("All traits must be named!")
  if(is.null(cfg$sigma))
    stop("No value for sigma!")
  if(cfg$sigma <= 0)
    stop("Sigma value must be non-zero and positive!")
  if(is.null(cfg$step))
    stop("No value for step!")
  if(cfg$step < 0)
    stop("Step value must be positive!")
  if(is.null(cfg$alpha)) {
    warning("No alpha value provided: assumed to be 0 (neutral).")
    cfg$alpha <- 0
  }
  if(cfg$alpha < 0)
    stop("Alpha value must be positive!")
  if(is.null(cfg$optima)) {
    if(cfg$alpha != 0)
      stop("No optimal value provided for a non-neutral process!")
  }
  if(length(cfg$optima) > 1L) {
    if(is.null(cfg$transition))
      stop("Many optima provided without a transition intensity matrix!")
    matrix(
      unlist(
        lapply(
          cfg$transition,
          function(x)
            lapply(
              x,
              FUN=function(x) if(is.null(x)) NA else x
            )
        )
      ),
      length(cfg$optima),
      length(cfg$optima),
      byrow = TRUE
    ) -> cfg$transition
    diag(cfg$transition) <- -rowSums(cfg$transition, na.rm = TRUE)
  }
  class(cfg) <- "TEMcfg"
  cfg
}
#' 
#' @rdname traitEvolMod
#' 
#' @examples ## The four trait evolution models can also be obtained at once as
#' ## follows:
#' tem <- lapply(cfg, traitEvolMod)
#' 
#' ## The first model is non-neutral with three optima:
#' tem[[1]]
#' 
#' ## The second model is neutral:
#' tem[[2]]
#' 
#' ## The third model is non-neutral with a single optimum:
#' tem[[3]]
#' 
#' ## The fourth trait is also non-neutral with a single optimum (but having a
#' ## different value than the latter):
#' tem[[4]]
#' 
#' ## Each model has a set of embedded member functions that can be called
#' ## directly.
#' 
#' ## Member $getName() returns the name of the trait as follows:
#' unlist(lapply(tem, function(x) x$getName()))
#' 
#' ## Member $getStep() returns the step sizes for the traits as follows:
#' unlist(lapply(tem, function(x) x$getStep()))
#' 
#' ## Member $setStep() sets the step sizes as follows:
#' lapply(tem, function(x) x$setStep(0.1))
#' 
#' ## and returns the recalculated transition probability matrix of models
#' ## having a transition intensity matrix (i.e., for non-neutral models with
#' ## multiple trait optima).
#' 
#' ## This is the modified step sizes:
#' unlist(lapply(tem, function(x) x$getStep()))
#' 
#' ## Member $getOptima() returns the model optimum (if available) or a vector
#' ## of model optima (if multiple optima are available) as follows:
#' lapply(tem, function(x) x$getOptima())
#' 
#' ## Member $getTransition() returns the transition intensity matrix, when
#' ## available (NULL otherwise) as follows:
#' lapply(tem, function(x) x$getTransition())
#' 
#' ## Member $getProb() returns the transition intensity matrix, whenever
#' ## available (NULL otherwise) as follows:
#' lapply(tem, function(x) x$getProb())
#' 
#' ## When multiple optima are available, member $updateState() enables to
#' ## simulate the transition from one optimal trait state (the one given as the
#' ## argument) to another (the one which is returned by the function):
#' state <- 1
#' newstate <- tem[[1]]$updateState(state)
#' newstate
#' 
#' ## Member $updateValue() simulates the evolution of the trait from one value
#' ## to another. For a non-neutral with multiple optima, The trait state is
#' ## provided using argument state (default: 1, which is the only applicable
#' ## value for a single optimum).
#' oldvalue <- 31.5
#' newvalue <- tem[[1]]$updateValue(oldvalue, state=1)
#' newvalue
#' 
#' ## Member $dumpConfig() returns the configuration list, which can be used
#' 
#' cfg2 <- lapply(tem, function(x) x$dumpConfig())
#' tem2 <- lapply(cfg2, traitEvolMod)
#' tem2
#' 
#' ## Clean up:
#' rm(cfg2, tem2)
#' 
#' ## Simulate trait evolution using the four models described previously:
#' 
#' ## Set step size to 0.05
#' lapply(tem, function(x, a) x$setStep(a), a = 0.05)
#' 
#' ## Results list:
#' res <- NULL
#' trNms <- lapply(tem, function(x) x$getName())
#' res$state <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
#' res$optim <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
#' res$trait <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
#' rm(trNms)
#' 
#' ## res$state contains the trait state at the simulation time
#' ## res$optim contains the trait's optimum value at the simulation time
#' ## res$trait contains the trait value at the simulation time
#' 
#' ## Setting the optimal state of the four traits at the beginning of the
#' ## simulation period:
#' res$state[1L,] <- c(2,NA,1,1)  ## NBL trait #2 evolves neutrally.
#' 
#' ## Getting the trait optima at the beginning of the simulation period:
#' unlist(
#'   lapply(
#'     tem,
#'     function(x)
#'       x$getOptima()[res$state[1L,x$getName()]]
#'   )
#' ) -> res$optim[1L,]
#' 
#' ## Setting the initial trait values:
#' res$trait[1L,] <- c(50,0,15,-25)
#' 
#' ## The state of the simulation at the beginning of the simulation:
#' head(res$state)
#' head(res$optim)
#' head(res$trait)
#' 
#' ## Setting RNG state to obtain 
#' set.seed(1234567)
#' 
#' ## This loop simulates time steps #2 through #1001
#' for(i in 2L:1001L) {
#' 
#'   ## Simulate the evolution of the trait states (if relevant, which it is
#'   ## only for the first trait):
#'   unlist(
#'     lapply(
#'       tem,
#'       function(x)
#'         x$updateState(res$state[i - 1L,x$getName()])
#'     ) 
#'   )-> res$state[i,]
#'   
#'   ## Obtain the optimal trait value (relevant for all traits but the second,
#'   ## trait, which evolves neutrally):
#'   unlist(
#'     lapply(
#'       tem,
#'       function(x)
#'         x$getOptima()[res$state[i,x$getName()]]
#'     )
#'   ) -> res$optim[i,]
#'   
#'   ## Simulate the evolution of the trait value:
#'   unlist(
#'     lapply(
#'       tem,
#'       function(x)
#'         x$updateValue(
#'           state = res$state[i,x$getName()],
#'           value = res$trait[i - 1L,x$getName()]
#'         )
#'     )
#'   ) -> res$trait[i,]
#' }
#' 
#' ## Plot the results:
#' par(mar=c(4,4,1,1))
#' plot(NA, xlim=c(0,0.05*1000), ylim=range(res$trait), xlab="Time",
#'      ylab="Trait value", las=1L)
#' lines(x=0.05*(0:1000), y=res$trait[,1L], col="black")
#' if(!is.na(tem[[1L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,1L], col="black", lty=3L)
#' lines(x=0.05*(0:1000), y=res$trait[,2L], col="red")
#' if(!is.na(tem[[2L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,2L], col="red", lty=3L)
#' lines(x=0.05*(0:1000), y=res$trait[,3L], col="blue")
#' if(!is.na(tem[[3L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,3L], col="blue", lty=3L)
#' lines(x=0.05*(0:1000), y=res$trait[,4L], col="green")
#' if(!is.na(tem[[4L]]$getOptima())[1L])
#'   lines(x=0.05*(0:1000), y=res$optim[,4L], col="green", lty=3L)
#' 
#' ## Save the figure to a PNG file:
#' dev.copy(png, filename="Trait value simulation (linear).png", width=600,
#' height=400)
#' dev.off()
#' 
#' @export
traitEvolMod <- function(par) {
  if(!inherits(par, "TEMcfg"))
    stop("The argument must be a of class TEMcfg!")
  name <- par$name
  sigma <- par$sigma
  step <- par$step
  alpha <- par$alpha
  optima <- par$optima
  transition <- par$transition
  nstate <- length(optima)
  
  if(!is.null(transition)) {
    prob <- eigen(transition)
    prob$vectors %*%
      diag(exp(step*prob$values)) %*%
      t(prob$vectors) -> prob
  } else
    prob <- NULL
  
  updateState <- function(state) {
    if(missing(state))
      stop("An initial state must be provided!")
    if(length(optima) > 1L)
      state <- sample(x=nstate, size=1L, prob=prob[state,])
    state
  }
  
  updateValue <- function(value, state = 1L) {
    if(missing(value))
      stop("No initial value!")
    if(alpha) {
      w <- exp(-alpha*step)
      s <- sigma*sqrt((1 - exp(-2*alpha*step))/(2*alpha))
      value <- rnorm(1L, w*value + (1 - w)*optima[state], s)
    } else {
      s <- sigma*sqrt(step)
      value <- rnorm(1L, value, s)
    }
    value
  }
  
  structure(
    list(
      getName = function() name,
      getStep = function() step,
      setStep = function(step) {
        step <<- step
        if(!is.null(transition)) {
          prob <- eigen(transition)
          prob$vectors %*%
            diag(exp(step*prob$values)) %*%
            t(prob$vectors) ->> prob
        }
      },
      getOptima = function() if(is.null(optima)) NA else optima,
      getTransition = function() transition,
      getProb = function() prob,
      updateState = updateState,
      updateValue = updateValue,
      dumpConfig = function()
        structure(
          list(
            name = name,
            sigma = sigma,
            step = step,
            alpha = alpha,
            optima = optima,
            transition = transition
          ),
          class = "TEMcfg"
        ),
      print = function() {
        cat("\nQuantitative trait evolution model (Ornstein-Uhlenbeck process)",
            "\n----------------------------------------------\n")
        cat("Name:", name, "\n")
        if(alpha) {
          cat("Trait evolving non-neutrally (alpha = ",alpha,")\n",sep="")
        } else
          cat("Trait evolving neutrally (alpha = 0)\n")
        cat("Sigma = ",sigma,"\n",sep="")
        if(!is.null(optima))
          if(nstate > 1L) {
            cat("Trait optima:",paste(optima,collapse=", "),"\n")
          } else
            cat("Trait optimum:",optima,"\n")
        if(!is.null(transition)) {
          cat("Transition intensity matrix:\n")
          print(transition)
        }
        cat("\n")
        invisible(NULL)
      }
    ),
    class = "TEM"
  )
}
#' 
#' @rdname traitEvolMod
#' 
#' @method print TEM
#' 
#' @export
print.TEM <- function(x, ...)
  x$print()
#' 
