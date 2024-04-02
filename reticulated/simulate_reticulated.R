
## rm(list=ls())
library(yaml)
library(Rcpp)
library(PhyloSim)  ## detach("package:PhyloSim", unload=TRUE)
library(ape)

source("aux.R")

## Reticulated sequence simulation

## Load the example of a configuration file provided with the
## package:
yaml.load_file(
  system.file(
    package = "PhyloSim",
    "extdata",
    "evolMod.yml"
  )
) -> dnaParam

## The configuration list is located in member `$DNA`
read.mutationMat(
  dnaParam$DNA 
) -> I

## The transition intensity matrix:
I

## The trait evolution model parameters:
dnaParam$trait %>%
  lapply(read.TraitEvolMod) -> traitMod

simulate_network(
  I = I,
  NS = 100,
  NN = 200,
  traitMod = traitMod,
  initState = c(2,NA,1,1),
  initValue = c(50,0,15,-25),
  NC = function() 1L + rpois(1L, 1.25),
  NP = function() 1L + rpois(1L, 0.25),
  gamma.shape = 5,
  gamma.scale = 5e-4,
  timestep = 1,
  maxDst = 4,
  prob = c(0.3,rep(0.175,4L)),
  removeGapOnly = TRUE,
  verbose = TRUE
) -> res

## res$trait
## res$rate
## res$em
## res$sqn
## res$cld
## res$dst
## res$edge
## res$cont

## Get the concatenated sequences
show.sequences(
  concatenate(res$sqn),
  cex=0.75,
  font=2L
)

show.sequences(
  concatenate(res$sqn, discard = "-"),
  cex=0.75,
  font=2L
)

## Display the network as a round diagram:
network_round_plot(res)

## These are the trait evolution models that were used:
res$trait$tem

## Obtain the nodes that are tips:
res$edge %>%
  lapply(function(x) !length(x)) %>%
  unlist %>%
  which -> tips

## Plot the trait simulation results for the first trail:
i <- getDescTrail(res, tips[1L])

par(mar=c(4,4,1,1))
plot(NA, xlim=c(0,(length(i) - 1)), ylim=range(res$trait$traitLog[i,]),
     xlab="Time", ylab="Trait value", las=1L)
lines(x=0:(length(i) - 1), y=res$trait$traitLog[i,1L], col="black")
if(!is.na(res$trait$tem[[1L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res$trait$optimLog[i,1L], col="black", lty=3L)
lines(x=0:(length(i) - 1), y=res$trait$traitLog[i,2L], col="red")
if(!is.na(res$trait$tem[[2L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res$trait$optimLog[i,2L], col="red", lty=3L)
lines(x=0:(length(i) - 1), y=res$trait$traitLog[i,3L], col="blue")
if(!is.na(res$trait$tem[[3L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res$trait$optimLog[i,3L], col="blue", lty=3L)
lines(x=0:(length(i) - 1), y=res$trait$traitLog[i,4L], col="green")
if(!is.na(res$trait$tem[[4L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res$trait$optimLog[i,4L], col="green", lty=3L)

## A larger example:
simulate_network(
  I = I,
  NS = 25000,
  NN = 250,
  traitMod = traitMod,
  initState = c(2,NA,1,1),
  initValue = c(50,0,15,-25),
  NC = function() 1L + rpois(1L, 0.75),
  NP = function() 1L + rpois(1L, 0.25),
  gamma.shape = 5,
  gamma.scale = 5e-4,
  timestep = 1,
  maxDst = 4,
  prob = c(0.3,rep(0.175,4L)),
  removeGapOnly = TRUE,
  verbose = TRUE
) -> res2

res2$edge %>%
  lapply(function(x) length(x)) %>%
  unlist -> desc

## Obtain the tips (node without children):
res2$edge %>%
  lapply(function(x) length(x) == 0L) %>%
  unlist %>%
  which -> tips

## any(desc[tips] != 0L)

## Plot the trait simulation results for the first trail:
i <- getDescTrail(res2, tips[100L])

par(mar=c(4,4,1,1))
plot(NA, xlim=c(0,(length(i) - 1)), ylim=range(res2$trait$traitLog[i,]),
     xlab="Time", ylab="Trait value", las=1L)
lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,1L], col="black")
if(!is.na(res2$trait$tem[[1L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,1L], col="black", lty=3L)
lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,2L], col="red")
if(!is.na(res2$trait$tem[[2L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,2L], col="red", lty=3L)
lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,3L], col="blue")
if(!is.na(res2$trait$tem[[3L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,3L], col="blue", lty=3L)
lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,4L], col="green")
if(!is.na(res2$trait$tem[[4L]]$getOptima())[1L])
  lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,4L], col="green", lty=3L)

## Clean up:
rm(dnaParam,I,traitMod,res,res2,show.sequences,getDescTrail,i,tips,desc,
   network_round_plot)
