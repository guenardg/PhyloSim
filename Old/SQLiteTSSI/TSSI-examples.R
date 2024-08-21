
### TSSI examples

if(recalculateExamples) {
  
  system.file(
    package = "PhyloSim",
    "extdata",
    "evolMod.yml"
  ) %>%
    yaml.load_file -> dnaParam
  
  dnaParam$DNA %>%
    read.mutationMat -> I
  
  dnaParam$trait %>%
    lapply(read.TraitEvolMod) -> traitMod
  
  set.seed(162745L)
  
  net <- openNetwork()
  
  makeReticulatedNetwork(
    net,
    NS = 400,
    NC = function() 1 + rpois(1,2),
    NP = function() 1 + rpois(1,1),
    timestep = function() exp(rnorm(1,0,0.5)),
    maxDiss = function() runif(1, 2, 3),
    contrib = function(x) rep(1/x,x),
    root = TRUE,
    verbose = TRUE
  )
  
  data.frame(
    name = c("SEQ1","SEQ2","SEQ3","SEQ4","SEQ5"),
    NN = c(500,800,120,700,230)
  ) -> cond
  
  for(i in 1:NROW(cond))
    simulateSequence(
      net,
      I,
      name = cond$name[i],
      sqn = drawDNASequence(cond$NN[i]),
      rate = drawEvolRate(cond$NN[i])
    )
  
  data.frame(
    name = sprintf("trait%d",1:4),
    state = c(2,NA,1,1),
    value = c(50,0,15,-25),
    note = c(
      "Non-neutral (OU) with multiple optima",
      "Neutral",
      "Non-neutral (OU) with single optimum",
      "Non-neutral (OU) with single optimum"
    )
  ) -> cond
  
  for(i in 1:NROW(cond))
    simulateTrait(
      net,
      tem = traitEvolMod(traitMod[[i]]),
      name = cond$name[i],
      state = cond$state[i],
      value = cond$value[i],
      note = cond$note[i]
    )
  
  closeNetwork(net, TRUE, TRUE, "reticulated-ex1")
  
  rm(net)
  
}

## net1 <- openNetwork("reticulated-ex1", load=TRUE)
## closeNetwork(net1, FALSE)

if(recalculateExamples) {
  
  system.file(
    package = "PhyloSim",
    "extdata",
    "evolMod.yml"
  ) %>%
    yaml.load_file -> dnaParam
  
  dnaParam$DNA %>%
    read.mutationMat -> I
  
  dnaParam$trait %>%
    lapply(read.TraitEvolMod) -> traitMod
  
  set.seed(263747L)
  
  net <- openNetwork()
  
  makeReticulatedNetwork(
    net,
    NS = 100,
    NC = function() 1 + rpois(1,1),
    NP = function() 1 + rpois(1,0.5),
    timestep = function() exp(rnorm(1,0.5,0.25)),
    maxDiss = function() runif(1, 2, 3),
    contrib = function(x) rep(1/x,x),
    root = TRUE,
    verbose = TRUE
  )
  
  data.frame(
    name = c("SEQ1","SEQ2","SEQ3"),
    NN = c(50,80,70)
  ) -> cond
  
  for(i in 1:NROW(cond))
    simulateSequence(
      net,
      I,
      name = cond$name[i],
      sqn = drawDNASequence(cond$NN[i]),
      rate = drawEvolRate(cond$NN[i])
    )
  
  data.frame(
    name = sprintf("trait%d",1:4),
    state = c(2,NA,1,1),
    value = c(50,0,15,-25),
    note = c(
      "Non-neutral (OU) with multiple optima",
      "Neutral",
      "Non-neutral (OU) with single optimum",
      "Non-neutral (OU) with single optimum"
    )
  ) -> cond
  
  for(i in 1:NROW(cond))
    simulateTrait(
      net,
      tem = traitEvolMod(traitMod[[i]]),
      name = cond$name[i],
      state = cond$state[i],
      value = cond$value[i],
      note = cond$note[i]
    )
  
  closeNetwork(net, TRUE, TRUE, "reticulated-ex2")
  
  rm(net)
  
}
## net2 <- openNetwork("reticulated-ex2", load=TRUE)
## closeNetwork(net2, FALSE)

if(recalculateExamples) {
  
  system.file(
    package = "PhyloSim",
    "extdata",
    "evolMod.yml"
  ) %>%
    yaml.load_file -> dnaParam
  
  dnaParam$DNA %>%
    read.mutationMat -> I
  
  dnaParam$trait %>%
    lapply(read.TraitEvolMod) -> traitMod
  
  set.seed(263747L)
  
  net <- openNetwork("reticulated-ex3.sim")
  
  makeReticulatedNetwork(
    net,
    NS = 10000,
    NC = function() 1 + rpois(1,1.5),
    NP = function() 1 + rpois(1,0.75),
    timestep = function() exp(rnorm(1,0.75,0.5)),
    maxDiss = function() runif(1, 2, 3),
    contrib = function(x) rep(1/x,x),
    root = TRUE,
    verbose = TRUE
  )
  
  closeNetwork(net)
  
  rm(net)
  
}

