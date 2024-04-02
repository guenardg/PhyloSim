
##

## rm(list=ls())
library(yaml)
library(Rcpp)
library(PhyloSim)  ## detach("package:PhyloSim", unload=TRUE)
library(ape)

## source("PhyloNet-aux.R")
## sourceCpp("PhyloNet-dev.cpp")
## file.copy("evolMod.yml", "package/PhyloSim/inst/extdata/", overwrite=TRUE)

## Trait Evolution Simulator

if(FALSE) {
  
  ## Load the example of a configuration file provided with the
  ## package:
  yaml.load_file(
    system.file(
      package = "PhyloSim",
      "extdata",
      "evolMod.yml"
    )
  ) -> dnaParam
  
  ## This file contains four trait configurations.
  
  ## The first trait file appears as follows:
  
  dnaParam$trait[[1]]
  
  ## The second, third, fourth appears as follows:
  
  dnaParam$trait[[2]]
  
  dnaParam$trait[[3]]
  
  ## and
  
  dnaParam$trait[[4]]
  
  ## The configuration list is obtained from the raw file as follows:
  
  cfg <- read.TraitEvolMod(dnaParam$trait[[1]])
  cfg
  
  ## All four lists are obtained as follows:
  
  cfg <- lapply(dnaParam$trait, read.TraitEvolMod)
  
  cfg[[1]]
  
  cfg[[2]]
  
  cfg[[3]]
  
  cfg[[4]]
  
  ## The four trait evolution models can also be obtained at once as
  ## follows:
  tem <- lapply(cfg, traitEvolMod)
  
  ## The first model is non-neutral with three optima:
  tem[[1]]
  
  ## The second model is neutral:
  tem[[2]]
  
  ## The third model is non-neutral with a single optimum:
  tem[[3]]
  
  ## The fourth trait is also non-neutral with a single optimum (but having a
  ## different value than the latter):
  tem[[4]]
  
  ## Each model has a set of embedded member functions that can be called
  ## directly.
  
  ## Member $getName() returns the name of the trait as follows:
  unlist(lapply(tem, function(x) x$getName()))
  
  ## Member $getStep() returns the step sizes for the traits as follows:
  unlist(lapply(tem, function(x) x$getStep()))
  
  ## Member $setStep() sets the step sizes as follows:
  lapply(tem, function(x) x$setStep(0.1))
  
  ## and returns the recalculated transition probability matrix of models
  ## having a transition intensity matrix (i.e., for non-neutral models with
  ## multiple trait optima).
  
  ## This is the modified step sizes:
  unlist(lapply(tem, function(x) x$getStep()))
  
  ## Member $getOptima() returns the model optimum (if available) or a vector
  ## of model optima (if multiple optima are available) as follows:
  lapply(tem, function(x) x$getOptima())
  
  ## Member $getTransition() returns the transition intensity matrix, when
  ## available (NULL otherwise) as follows:
  lapply(tem, function(x) x$getTransition())
  
  ## Member $getProb() returns the transition intensity matrix, whenever
  ## available (NULL otherwise) as follows:
  lapply(tem, function(x) x$getProb())
  
  ## When multiple optima are available, member $updateState() enables to
  ## simulate the transition from one optimal trait state (the one given as the
  ## argument) to another (the one which is returned by the function):
  state <- 1
  newstate <- tem[[1]]$updateState(state)
  newstate
  
  ## Member $updateValue() simulates the evolution of the trait from one value
  ## to another. For a non-neutral with multiple optima, The trait state is
  ## provided using argument state (default: 1, which is the only applicable
  ## value for a single optimum).
  oldvalue <- 31.5
  newvalue <- tem[[1]]$updateValue(oldvalue, state=1)
  newvalue
  
  ## Member $dumpConfig() returns the configuration list, which can be used
  
  cfg2 <- lapply(tem, function(x) x$dumpConfig())
  tem2 <- lapply(cfg2, traitEvolMod)
  tem2
  
  ## Clean up:
  rm(cfg2, tem2)
  
  ## Simulate trait evolution using the four models described previously:
  
  ## Set step size to 0.05
  lapply(tem, function(x, a) x$setStep(a), a = 0.05)
  
  ## Results list:
  res <- NULL
  trNms <- lapply(tem, function(x) x$getName())
  res$state <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
  res$optim <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
  res$trait <- matrix(NA, 1001L, length(tem), dimnames = list(NULL, trNms))
  rm(trNms)
  
  ## res$state contains the trait state at the simulation time
  ## res$optim contains the trait's optimum value at the simulation time
  ## res$trait contains the trait value at the simulation time
  
  ## Setting the optimal state of the four traits at the beginning of the
  ## simulation period:
  res$state[1L,] <- c(2,NA,1,1)  ## NBL trait #2 evolves neutrally.
  
  ## Getting the trait optima at the beginning of the simulation period:
  unlist(
    lapply(
      tem,
      function(x)
        x$getOptima()[res$state[1L,x$getName()]]
    )
  ) -> res$optim[1L,]
  
  ## Setting the initial trait values:
  res$trait[1L,] <- c(50,0,15,-25)
  
  ## The state of the simulation at the beginning of the simulation:
  head(res$state)
  head(res$optim)
  head(res$trait)
  
  ## Setting RNG state to obtain 
  set.seed(1234567)
  
  ## This loop simulates time steps #2 through #1001
  for(i in 2L:1001L) {
    
    ## Simulate the evolution of the trait states (if relevant, which it is
    ## only for the first trait):
    unlist(
      lapply(
        tem,
        function(x)
          x$updateState(res$state[i - 1L,x$getName()])
      ) 
    )-> res$state[i,]
    
    ## Obtain the optimal trait value (relevant for all traits but the second,
    ## trait, which evolves neutrally):
    unlist(
      lapply(
        tem,
        function(x)
          x$getOptima()[res$state[i,x$getName()]]
      )
    ) -> res$optim[i,]
    
    ## Simulate the evolution of the trait value:
    unlist(
      lapply(
        tem,
        function(x)
          x$updateValue(
            state = res$state[i,x$getName()],
            value = res$trait[i - 1L,x$getName()]
          )
      )
    ) -> res$trait[i,]
  }
  
  ## Plot the results:
  par(mar=c(4,4,1,1))
  plot(NA, xlim=c(0,0.05*1000), ylim=range(res$trait), xlab="Time",
       ylab="Trait value", las=1L)
  lines(x=0.05*(0:1000), y=res$trait[,1L], col="black")
  if(!is.na(tem[[1L]]$getOptima())[1L])
    lines(x=0.05*(0:1000), y=res$optim[,1L], col="black", lty=3L)
  lines(x=0.05*(0:1000), y=res$trait[,2L], col="red")
  if(!is.na(tem[[2L]]$getOptima())[1L])
    lines(x=0.05*(0:1000), y=res$optim[,2L], col="red", lty=3L)
  lines(x=0.05*(0:1000), y=res$trait[,3L], col="blue")
  if(!is.na(tem[[3L]]$getOptima())[1L])
    lines(x=0.05*(0:1000), y=res$optim[,3L], col="blue", lty=3L)
  lines(x=0.05*(0:1000), y=res$trait[,4L], col="green")
  if(!is.na(tem[[4L]]$getOptima())[1L])
    lines(x=0.05*(0:1000), y=res$optim[,4L], col="green", lty=3L)
  
  ## Save the figure to a PNG file:
  dev.copy(png, filename="Trait value simulation (linear).png", width=600,
           height=400)
  dev.off()
  
  ## Clean up
  rm(cfg,dnaParam,res,tem,i,newstate,newvalue,oldvalue,state)
  
}

if(FALSE) {
  ## Define a raw vector for storing nuceotide values:
  c(Sequence_1 = "ATCG-TTTCG--CCCCA--TTA--TTAA-GTAA-GTAATCTTTCA",
    Sequence_2 = "TTGGCTTCC--TC--CAT-TTC--TTCA-GT-ACG-ATTCTTTTA",
    Sequence_3 = "TTCG-TACC-T-T---A-ATAA--T-AA-GTCTTGTAATCGTTCA") %>%
    sapply(charToRaw) %>%
    t -> sqn
  
  ## Display the raw sequence:
  sqn
  
  ## Transforming the sequence to character strings
  concatenate(sqn)

  ## Transforming the sequence to character strings without the gaps:
  concatenate(sqn, discard="-")
  
  ## Clean-up:
  rm(sqn)
}

if(FALSE) {
  ## Define a raw vector for storing nuceotide values:
  c(Sequence_1 = "ATCG-TTTCG--CCCCA--TTA--TTAA-GTAA-GTAATCTTTCA",
    Sequence_2 = "TTGGCTTCC--TC--CAT-TTC--TTCA-GT-ACG-ATTCTTTTA",
    Sequence_3 = "TTCG-TACC-T-T---A-ATAA--T-AA-GTCTTGTAATCGTTCA") %>%
    sapply(charToRaw) %>%
    t -> sqn
  
  ## Display the raw sequence:
  sqn
  
  ## Transforming the sequence to character strings
  tmp <- concatenate(sqn)
  tmp
  
  ## Transforming the sequence to character strings without the gaps:
  concatenate(sqn, discard="-")
  
  ## Write the sequences into a text file:
  write.fasta(tmp, file="Sequences.fst", linebreak=15)
  
  ## Clean-up:
  file.remove("Sequences.fst")
  rm(sqn, tmp)
}

if(FALSE) {
  
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
  
  ## Implement the molecular evolution model for a single nucleotide:
  molEvolMod(I, 1, 1) -> em1
  
  ## Get the transition probability matrix as follows:
  em1$getMt()
  
  ## A vector of raw as examples of initial traits:
  tr <- charToRaw("-ACGT")
  
  ## Simulate molecular evolution from:
  rawToChar(em1$evolve(tr[1L]))    ## a gap.
  rawToChar(em1$evolve(tr[2L]))    ## an adenine base.
  rawToChar(em1$evolve(tr[3L]))    ## a cytosine base.
  rawToChar(em1$evolve(tr[4L]))    ## a guanine base.
  rawToChar(em1$evolve(tr[5L]))    ## a thymine base.
  
  ## Recalculate the probabilities for a lower mean evolution rate (one tenth
  ## the previous one):
  em1$recalculate(1, 0.1)
  
  em1$getMt()        ## The recalculated transition probability matrix.
  
  ## Simulate molecular evolution from:
  rawToChar(em1$evolve(tr[1L]))    ## a gap.
  rawToChar(em1$evolve(tr[2L]))    ## an adenine base.
  rawToChar(em1$evolve(tr[3L]))    ## a cytosine base.
  rawToChar(em1$evolve(tr[4L]))    ## a guanine base.
  rawToChar(em1$evolve(tr[5L]))    ## a thymine base.
  
  ## Base changes are now less probable.
  
  ## Clean up:
  rm(dnaParam, I, em1, tr)
  
}

if(FALSE) {
  
  ## Simple sequence simulation
  
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
  
  ## Setting RNG seed:
  set.seed(12345)
  
  ## Generate the sequences:
  simulate_linear(
    I = I,
    NS = 35,
    NN = 150,
    gamma.shape = 5,
    gamma.scale = 5e-4,
    timestep = 1
  ) -> res
  
  ## Get the concatenated sequences
  original <- concatenate(res$sqn)
  nogaps <-  concatenate(res$sqn, discard = "-")
  
  ## To save the sequences:
  ## write.fasta(original, "Simulated sequences - aligned.fst",
  ##             linebreak = 70L)
  ## write.fasta(nogaps, "Simulated sequences - raw.fst", linebreak = 70L)
  
  ## Function for showing sequence:
  show.sequences <- function(x, xlim, ylim, xlab="Position", text = FALSE,
    col=c(A="red",`T`="blue",C="green",G="yellow", `-`="grey"),
    ...) {
    
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
  
  ## Showing the resulting sequences:
  show.sequences(original, cex=0.75, font=2L)   ## Natively aligned
  show.sequences(nogaps, cex=0.75, font=2L)     ## Raw
  
  ## Sequence simulation with traits
  
  ## Resetting RNG seed:
  set.seed(12345)
  
  ## The trait evolution model parameters:
  dnaParam$trait %>%
    lapply(read.TraitEvolMod) -> traitMod
  
  ## Generate the sequences:
  simulate_linear(
    I = I,
    NS = 300,
    NN = 150,
    traitMod = traitMod,
    initState = c(2,NA,1,1),
    initValue = c(50,0,15,-25),
    gamma.shape = 5,
    gamma.scale = 5e-4,
    timestep = 1
  ) -> res2
  
  ## These are the trait evolution models that were used:
  res2$trait$tem
  
  ## Plot the trait simulation results:
  par(mar=c(4,4,1,1))
  plot(NA, xlim=c(0,300), ylim=range(res2$trait$traitLog), xlab="Time",
       ylab="Trait value", las=1L)
  lines(x=0:299, y=res2$trait$traitLog[,1L], col="black")
  if(!is.na(res2$trait$tem[[1L]]$getOptima())[1L])
    lines(x=0:299, y=res2$trait$optimLog[,1L], col="black", lty=3L)
  lines(x=0:299, y=res2$trait$traitLog[,2L], col="red")
  if(!is.na(res2$trait$tem[[2L]]$getOptima())[1L])
    lines(x=0:299, y=res2$trait$optimLog[,2L], col="red", lty=3L)
  lines(x=0:299, y=res2$trait$traitLog[,3L], col="blue")
  if(!is.na(res2$trait$tem[[3L]]$getOptima())[1L])
    lines(x=0:299, y=res2$trait$optimLog[,3L], col="blue", lty=3L)
  lines(x=0:299, y=res2$trait$traitLog[,4L], col="green")
  if(!is.na(res2$trait$tem[[4L]]$getOptima())[1L])
    lines(x=0:299, y=res2$trait$optimLog[,4L], col="green", lty=3L)
  
  original2 <- concatenate(res2$sqn)
  
  show.sequences(original2, cex=0.75, font=2L)
  
  ## Clean up:
  rm(dnaParam,I,traitMod,res,res2,nogaps,original,original2,show.sequences)
  
}

if(FALSE) {
  
  ## Tree sequence simulation
  
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
  
  ## Resetting RNG seed:
  set.seed(12345)
  
  ## Generate the sequences:
  simulate_tree(
    I = I,
    NS = 300,
    NN = 150
  ) -> res
  
  ## Get the concatenated sequences
  original <- concatenate(res$sqn)
  nogaps <-  concatenate(res$sqn, discard = "-")
  
  ## To save the sequences:
  ## write.fasta(original, "Simulated sequences - aligned.fst",
  ##             linebreak = 70L)
  ## write.fasta(nogaps, "Simulated sequences - raw.fst", linebreak = 70L)
  
  ## Function for showing sequence:
  show.sequences <- function(
    x, xlim, ylim, xlab = "Position", text = FALSE,
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
  
  ## Showing the resulting sequences:
  show.sequences(original, cex=0.75, font=2L)   ## Natively aligned
  show.sequences(nogaps, cex=0.75, font=2L)     ## Raw
  
  ## Sequence simulation with traits
  
  ## Resetting RNG seed:
  set.seed(12345)
  
  ## The trait evolution model parameters:
  dnaParam$trait %>%
    lapply(read.TraitEvolMod) -> traitMod
  
  ## Generate the sequences with the traits:
  simulate_tree(
    I = I,
    NS = 300,
    NN = 150,
    traitMod = traitMod,
    initState = c(2,NA,1,1),
    initValue = c(50,0,15,-25),
    NC = function() 1L + rpois(1L, 0.25)
  ) -> res2
  
  ## These are the trait evolution models that were used:
  res2$trait$tem
  
  ## Obtain the nodes that are tips:
  res2$edge %>%
    {unique(.$to)[!(unique(.$to) %in% unique(.$from))]} %>%
    {.[order(.)]} -> tips
  
  ## A function to obtain all the descendents from a tip:
  getDescendents <- function(x, y) {
    while(head(y,1L) != 1)
      y <- c(x$edge$from[head(y,1L) == x$edge$to],y)
    y
  }
  
  ## Get the concatenated sequences
  original2 <- concatenate(res2$sqn)
  
  ## Get the trail from the first tip:
  show.sequences(
    original2[getDescendents(res2, tips[1L])],
    cex=0.75,
    font=2L
  )
  
  ## Get the trail from the tenth tip:
  show.sequences(
    original2[getDescendents(res2, tips[10L])],
    cex=0.75,
    font=2L
  )
  
  ## Plot the trait simulation results for the first trail:
  i <- getDescendents(res2, tips[1L])
  
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
  
  ## Plot the trait simulation results for the tenth trail:
  i <- getDescendents(res2, tips[10L])
  
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
  rm(dnaParam,I,traitMod,res,res2,nogaps,original,original2,show.sequences,
     getDescendents,i,tips)
  
}

## rm(l)

if(FALSE) {
  
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
  
  ## Setting seed:
  set.seed(12345)
  
  ## Generate the sequences:
  simulate_network(
    I = I,
    NS = 300,
    NN = 150,
    NC = function() 1L + rpois(1L, 1.5),
    NP = function() 1L + rpois(1L, 0.5),
    verbose = FALSE
  ) -> res
  
  ## Get the concatenated sequences
  original <- concatenate(res$sqn)
  nogaps <-  concatenate(res$sqn, discard = "-")
  
  ## To save the sequences:
  ## write.fasta(original, "Simulated sequences - aligned.fst",
  ##             linebreak = 70L)
  ## write.fasta(nogaps, "Simulated sequences - raw.fst", linebreak = 70L)
  
  ## Function for showing sequence:
  show.sequences <- function(
    x, xlim, ylim, xlab = "Position", text = FALSE,
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
  
  ## Showing the resulting sequences:
  show.sequences(original, cex=0.75, font=2L)   ## Natively aligned
  show.sequences(nogaps, cex=0.75, font=2L)     ## Raw
  
  ## Sequence simulation with traits
  
  ## Resetting RNG seed:
  set.seed(12345)
  
    ## The trait evolution model parameters:
  dnaParam$trait %>%
    lapply(read.TraitEvolMod) -> traitMod
  
  ## Generate the sequences with the traits:
  simulate_network(
    I = I,
    NS = 300,
    NN = 150,
    traitMod = traitMod,
    initState = c(2,NA,1,1),
    initValue = c(50,0,15,-25),
    NC = function() 1L + rpois(1L, 1.5),
    NP = function() 1L + rpois(1L, 0.5)
  ) -> res2
  
  ## These are the trait evolution models that were used:
  res2$trait$tem
  
  ## Obtain the nodes that are tips:
  res2$edge %>%
    lapply(function(x) !length(x)) %>%
    unlist %>%
    which -> tips
  
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
  
  ## Get the concatenated sequences
  original2 <- concatenate(res2$sqn)
  
  ## Get the trail from the first tip:
  show.sequences(
    original2[getDescTrail(res2, tips[1L])],
    cex=0.75,
    font=2L
  )
  
  ## Get the trail from the tenth tip:
  show.sequences(
    original2[getDescTrail(res2, tips[10L])],
    cex=0.75,
    font=2L
  )
  
  ## Plot the trait simulation results for the first trail:
  i <- getDescTrail(res2, tips[1L])
  
  par(mar=c(4,4,1,1))
  plot(NA, xlim=c(0,(length(i) - 1)), ylim=range(res2$trait$traitLog[i,]),
       xlab="Time", ylab="Trait value", las=1L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,1L], col="black")
  if(!is.na(res2$trait$tem[[1L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,1L], col="black",
          lty=3L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,2L], col="red")
  if(!is.na(res2$trait$tem[[2L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,2L], col="red", lty=3L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,3L], col="blue")
  if(!is.na(res2$trait$tem[[3L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,3L], col="blue", lty=3L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,4L], col="green")
  if(!is.na(res2$trait$tem[[4L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,4L], col="green",
          lty=3L)
  
  ## Plot the trait simulation results for the tenth trail:
  i <- getDescTrail(res2, tips[10L])
  
  par(mar=c(4,4,1,1))
  plot(NA, xlim=c(0,(length(i) - 1)), ylim=range(res2$trait$traitLog[i,]),
       xlab="Time", ylab="Trait value", las=1L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,1L], col="black")
  if(!is.na(res2$trait$tem[[1L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,1L], col="black",
          lty=3L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,2L], col="red")
  if(!is.na(res2$trait$tem[[2L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,2L], col="red", lty=3L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,3L], col="blue")
  if(!is.na(res2$trait$tem[[3L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,3L], col="blue", lty=3L)
  lines(x=0:(length(i) - 1), y=res2$trait$traitLog[i,4L], col="green")
  if(!is.na(res2$trait$tem[[4L]]$getOptima())[1L])
    lines(x=0:(length(i) - 1), y=res2$trait$optimLog[i,4L], col="green",
          lty=3L)
  
  ## Clean up:
  rm(dnaParam,I,traitMod,res,res2,nogaps,original,original2,show.sequences,
     getDescTrail,i,tips)
  
}

## Garbage from here...












if(FALSE) {
  
  
  
  
  
  
  res$sqn %>%
    concatenate(width=70L) -> out
  names(out) <- sprintf("sequence%04d", 1L:length(out))
  write.fasta(out, "Simulated sequences.fst")
  sprintf(
    "clustalo --force -i %s -o %s",
    "'Simulated sequences.fst'",
    "'Aligned sequences - clustalo.fst'"
  ) %>%
    system
  
  sprintf(
    "muscle -in %s -out %s",
    "'Simulated sequences.fst'",
    "'Aligned sequences - muscle.fst'"
  ) %>%
    system
  
  sqn %>%
    concatenate(gap=TRUE, width=70L) -> out
  names(out) <- sprintf("sequence%04d", 1L:length(out))
  write.fasta(out, "Simulated sequences - aligned.fst")
  
  list(
    aligned = read.FASTA("Simulated sequences - aligned.fst"),
    clustalo = read.FASTA("Aligned sequences - clustalo.fst"),
    muscle = read.FASTA("Aligned sequences - muscle.fst")
  ) -> simseqs
  
  dnaDst <- lapply(simseqs, dist.dna, model="K80")
  
  data.frame(
    time = as.numeric(dst),
    aligned = as.numeric(dnaDst$aligned),
    clustalo = as.numeric(dnaDst$clustalo),
    muscle = as.numeric(dnaDst$muscle)
  ) -> dat
  
  par(mar = c(4,4,1,1))
  boxplot(
    aligned ~ time,
    data = dat,
    xlim = range(dat$time),
    ylim = range(dat$aligned),
    xlab = "Simulation time",
    ylab = "Phylogenetic distance",
    las = 1L
  )
  lmd <- lm(aligned~time, data=dat)
  abline(lmd)
  summary(lmd)
  
  par(mar = c(4,4,1,1))
  boxplot(
    clustalo ~ time,
    data = dat,
    xlim = range(dat$time),
    ylim = range(dat$clustalo),
    xlab = "Simulation time",
    ylab = "Phylogenetic distance",
    las = 1L
  )
  lmd <- lm(clustalo~time-1, data=dat)
  abline(lmd)
  summary(lmd)
  
  par(mar = c(4,4,1,1))
  boxplot(
    muscle ~ time,
    data = dat,
    xlim = range(dat$time),
    ylim = range(dat$muscle),
    xlab = "Simulation time",
    ylab = "Phylogenetic distance",
    las = 1L
  )
  lmd <- lm(muscle~time-1, data=dat)
  abline(lmd)
  summary(lmd)
  
  rm(dat,lmd)
  
  ## max(dst)
  nrow(sqn) %>%
    {sapply(1L:., function(x) max(dst[dst_idx(., x)]))} %>%
    which.max -> i
  
  (1L:nrow(sqn))[-i][which.max(dst[dst_idx(nrow(sqn), i)])] -> j
  
  dst[dst_idx(nrow(sqn), i, j)]
  
  sqn[c(i,j),] %>%
    concatenate(gap=TRUE, width=70L) -> out2
  
  names(out2) <- sprintf("sequence%04d", c(i,j))
  
  write.fasta(out2, "Simulated sequence - most different.fst")
  
  rm(i,j,out,out2)
  
  edge %>%
    {unique(.$to)[!(unique(.$to) %in% unique(.$from))]} %>%
    {.[order(.)]} -> tips
  
  ## max(dst)
  ## dst[dst_idx(nrow(sqn), 1, tips)]
  x <- tips[3L]
  while(head(x,1L) != 1)
    x <- c(edge$from[head(x,1L) == edge$to],x)
  stateLog[x,]
  traitLog[x,]
  optimLog[x,]
  
  
  
  
}
