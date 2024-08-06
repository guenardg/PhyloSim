
## SQLite-based reticulated evolution API

## rm(list=ls())
## library(yaml)
## library(Rcpp)
library(PhyloSim)  ## detach("package:PhyloSim", unload=TRUE)
## library(RSQLite)
## library(rgl)

source("aux.R")

## Load the example of a configuration file provided with the
## package:
system.file(
  package = "PhyloSim",
  "extdata",
  "evolMod.yml"
) %>%
  yaml.load_file -> dnaParam

## The configuration list is located in member `$DNA`
dnaParam$DNA %>%
  read.mutationMat -> I

## The transition intensity matrix:
I

## The trait evolution model parameters:
dnaParam$trait %>%
  lapply(read.TraitEvolMod) -> traitMod


### Linear case:

## Initializing a new data base to store the network:
net <- connectNetwork(name = "net_linear", init = TRUE)

## Setting RNG seed to have a consistent example:
set.seed(162745)

## These are the tables that are created in the database:
dbListTables(net$con)

## The following code enables one to extract the SQLite table information:
lapply(
  dbListTables(net$con),
  function(x, con)
    dbGetQuery(con, sprintf("PRAGMA table_info(%s)", x))
  ,
  con = net$con
) -> table_info

names(table_info) <- dbListTables(net$con)

table_info

## The following code retrieves the number of records in all the tables.
## All of them are empty at the beginning of the simulation:
sapply(
  dbListTables(net$con),
  function(x, con)
    unlist(
      dbGetQuery(
        con,
        sprintf(
          "SELECT COUNT(*) AS %s
         FROM %s",
          x,x),
      )
    ),
  con = net$con
)

## Create a linear network with 50 species having 1 descendant each with a
## constant time step of 1:
makeLinearNetwork(
  net = net,
  NS = 50,
  NC = function() 1,
  timestep = function() 1,
  root = TRUE,
  verbose = TRUE
)

## After the linear phylogenetic network has been created, tables diss, edge,
## and node conntain information. The other tables are used for storing the
## sequence and trait information.
sapply(
  dbListTables(net$con),
  function(x, con)
    unlist(
      dbGetQuery(
        con,
        sprintf("SELECT COUNT(*) FROM %s", x),
      )
    ),
  con = net$con
)

## Simulate two sequences, one named 'SEQ1' having 100 locations (i.e., one of
## the four nucleotide or a gap) and a second sequence named 'SEQ2' having 250
## locations.
data.frame(
  name = c("SEQ1","SEQ2"),
  NN = c(100,250)
) -> cond

## Calling function simulateSequence once for each of the simulated sequences,
## each time drawing a random root sequence (using drawDNASequence) and a set
## of evolution rate (using function drawEvolRate).
for(i in 1:NROW(cond))
  simulateSequence(
    net,
    I,
    name = cond$name[i],
    sqn = drawDNASequence(cond$NN[i]),
    rate = drawEvolRate(cond$NN[i])
  )

## Now, tables seq and seq_content contain information:
sapply(
  dbListTables(net$con),
  function(x, con)
    unlist(
      dbGetQuery(
        con,
        sprintf("SELECT COUNT(*) FROM %s", x),
      )
    ),
  con = net$con
)

## Show the first sequence (here, the y-axis corresponds to the time step):
net %>%
  getSequence(name = "SEQ1") %>%
  concatenate(discard = "-") %>%
  show.sequence

## Show the second sequence:
net %>%
  getSequence(name = "SEQ2") %>%
  concatenate(discard = "-") %>%
  show.sequence

## Set four traits to be simulated.
data.frame(
  name = sprintf("trait%d",1:4),
  state = c(2,NA,1,1),
  value = c(50,0,15,-25),
  note = c("Non-neutral (OU) with multiple optima","Neutral",
           "Non-neutral (OU) with single optimum",
           "Non-neutral (OU) with single optimum")
) -> cond

## The four traits are simulated as follows:
for(i in 1:NROW(cond))
  simulateTrait(
    net,
    tem = traitEvolMod(traitMod[[i]]),
    name = cond$name[i],
    state = cond$state[i],
    value = cond$value[i],
    note = cond$note[i]
  )

## Now, tables trait and trait_content also contain information:
sapply(
  dbListTables(net$con),
  function(x, con)
    unlist(
      dbGetQuery(
        con,
        sprintf("SELECT COUNT(*) FROM %s", x),
      )
    ),
  con = net$con
)

## To obtain the trait values at the nodes, one has to use an SQL query (in the
## SQLite dialect) as follows:
lapply(
  sprintf("trait%d",1:4),
  function(x, con)
    dbGetQuery(
      con,
      sprintf(
        "SELECT trait_content.node, trait_content.optim, trait_content.value
         FROM trait_content
         JOIN trait ON trait.rowid = trait_content.trait
         WHERE trait.name = '%s'
         ORDER BY trait_content.node",
        x
      )
    ),
  con = net$con
) -> sim_trait

## Plot the trait simulation results:
par(mar=c(4,4,1,1))
plot(NA, xlim=c(0,50), xlab="Time", ylab="Trait value", las=1,
     ylim=range(-25,sapply(sim_trait, function(x) x$value),80))
col <- c("black","red","blue","green")

## The solid line is the trait value, whereas the dotted line is the effective
## trait optimum:
for(i in 1:4) {
  lines(x=0:49, y=sim_trait[[i]]$value, col=col[i])
  if(!all(is.na(sim_trait[[i]]$optim)))
    lines(x=0:49, y=sim_trait[[i]]$optim, col=col[i], lty=3)
}
## Note: there are no effective optimum for a trait evolving neutrally.

## Clear the simulated linear network:
clearNetwork(net, TRUE)


### Tree case:

## Initializing a new data base to store the network:
net <- connectNetwork(name = "net_linear", init = TRUE)

## Setting RNG seed to have a consistent example:
set.seed(162745)

## Create a tree network (phylogenetic tree) with 100 species having
## a fix number of 2 descendants (dichotomic) and a random log-normal time step
## drawn at each time step:
makeTreeNetwork(
  net = net,
  NS = 100,
  NC = function() 2,
  timestep = function() exp(rnorm(1,0,0.5)),
  root = TRUE,
  verbose = TRUE
)

## Simulate two sequences, one named 'SEQ1' having 50 locations (i.e., one of
## the four nucleotide or a gap) and a second sequence named 'SEQ2' having 80
## locations.
data.frame(
  name = c("SEQ1","SEQ2"),
  NN = c(50,80)
) -> cond

## Calling function simulateSequence once for each of the simulated sequences,
## each time drawing a random root sequence (using drawDNASequence) and a set
## of evolution rate (using function drawEvolRate).
for(i in 1:NROW(cond))
  simulateSequence(
    net,
    I,
    name = cond$name[i],
    sqn = drawDNASequence(cond$NN[i]),
    rate = 10*drawEvolRate(cond$NN[i])  ## Accelerate evolution by 10 folds
  )

## Set four traits to be simulated.
data.frame(
  name = sprintf("trait%d",1:4),
  state = c(2,NA,1,1),
  value = c(50,0,15,-25),
  note = c("Non-neutral (OU) with multiple optima","Neutral",
           "Non-neutral (OU) with single optimum",
           "Non-neutral (OU) with single optimum")
) -> cond

## The four traits are simulated as follows:
for(i in 1:NROW(cond))
  simulateTrait(
    net,
    tem = traitEvolMod(traitMod[[i]]),
    name = cond$name[i],
    state = cond$state[i],
    value = cond$value[i],
    note = cond$note[i]
  )

## To obtain the trait values at the nodes, one has to use an SQL query (in the
## SQLite dialect) as follows:
lapply(
  sprintf("trait%d",1:4),
  function(x, con)
    dbGetQuery(
      con,
      sprintf(
        "SELECT trait_content.node, trait_content.optim, trait_content.value
         FROM trait_content
         JOIN trait ON trait.rowid = trait_content.trait
         WHERE trait.name = '%s'
         ORDER BY trait_content.node",
        x
      )
    ),
  con = net$con
) -> sim_trait

## Showing the sequences and traits evolving is not as straightforward for a
## phylogenetic tree as it is for a linear sequence of descendants.

## We first need to retrieve the edge list as follows:
dbGetQuery(
  net$con,
  "SELECT i,j FROM edge ORDER BY rowid"
) -> edge

## From that list, we can determine which node is a leaf:
leaf <- edge$j[!(edge$j %in% unique(edge$i))]
leaf

## We can use this simple function to obtain lineage of any node:
getLineage <- function(x, e) {
  while(length(wh <- which(e$j == head(x,1))))
    x <- c(e$i[wh], x)
  unname(x)
}

## Here is the lineage of the first leaf:
getLineage(leaf[1], edge)

## Show the first sequence (here, the y-axis corresponds to the time step):
net %>%
  getSequence(name = "SEQ1") %>%
  .[getLineage(leaf[1], edge),] %>%
  concatenate(discard = "-") %>%
  show.sequence

## Show the second sequence:
net %>%
  getSequence(name = "SEQ2") %>%
  .[getLineage(leaf[1], edge),] %>%
  concatenate(discard = "-") %>%
  show.sequence

## Obtain the distances between the root node and all the nodes:
rbind(
  list(1L,0),
  dbGetQuery(
    net$con,
    "SELECT j, d FROM diss WHERE i = 1 ORDER BY j"
  )
) -> dst

## Function to show every lineages 'l':
plotit <- function(l) {
  nn <- getLineage(leaf[l], edge)
  dd <- dst$d[match(nn, dst$j)]
  
  ## Plot the trait simulation results:
  par(mar=c(4,4,1,1))
  plot(NA, xlim=range(dd), xlab="Time", ylab="Trait value", las=1,
       ylim=range(-25,sapply(sim_trait, function(x) x$value),80))
  col <- c("black","red","blue","green")
  
  ## The solid line is the trait value, whereas the dotted line is the effective
  ## trait optimum:
  for(i in 1:4) {
    lines(x=dd, y=sim_trait[[i]]$value[match(nn, dst$j)], col=col[i])
    if(!all(is.na(sim_trait[[i]]$optim)))
      lines(x=dd, y=sim_trait[[i]]$optim[match(nn, dst$j)], col=col[i], lty=3)
  }
  ## Note: there are no effective optimum for a trait evolving neutrally.
}

## Plot the lineages of some of the leaves:
plotit(1)
plotit(2)
plotit(3)
plotit(6)
plotit(7)

## Clear the simulated linear network:
clearNetwork(net, TRUE)


### The reticulated case:

## Initializing a new data base to store the network:
net <- connectNetwork(name = "net_reticulated", init = TRUE)
## clearNetwork(net, TRUE)

## Setting RNG seed to have a consistent example:
set.seed(162745)

## Create a reticulated network with 100 species having a random number of
## descendants and parents per node, a random log-normal time step, a random
## maximum hybridization distance, and equal parental contributions on
## hybridization events.
makeReticulatedNetwork(
  net,
  NS = 100,
  NC = function() 1 + rpois(1,2),
  NP = function() 1 + rpois(1,1),
  timestep = function() exp(rnorm(1,0,0.5)),
  maxDiss = function() runif(1, 2, 3),
  contrib = function(x) rep(1/x,x),
  root = TRUE,
  verbose = TRUE
)

## Simulate one sequence named 'SEQ1' having 50 locations (i.e., one of
## the four nucleotide or a gap). Calling function simulateSequence drawing a
## random root sequence (using drawDNASequence) and a set of evolution rate
## (using function drawEvolRate).
simulateSequence(
  net,
  I,
  name = "SEQ1",
  sqn = drawDNASequence(50),
  rate = 10*drawEvolRate(50)  ## Accelerate evolution by 10 folds
)

## Set four traits to be simulated.
data.frame(
  name = sprintf("trait%d",1:4),
  state = c(2,NA,1,1),
  value = c(50,0,15,-25),
  note = c("Non-neutral (OU) with multiple optima","Neutral",
           "Non-neutral (OU) with single optimum",
           "Non-neutral (OU) with single optimum")
) -> cond

## The four traits are simulated as follows:
for(i in 1:NROW(cond))
  simulateTrait(
    net,
    tem = traitEvolMod(traitMod[[i]]),
    name = cond$name[i],
    state = cond$state[i],
    value = cond$value[i],
    note = cond$note[i]
  )

## To obtain the trait values at the nodes, one has to use an SQL query (in the
## SQLite dialect) as follows:
lapply(
  sprintf("trait%d",1:4),
  function(x, con)
    dbGetQuery(
      con,
      sprintf(
        "SELECT trait_content.node, trait_content.optim, trait_content.value
         FROM trait_content
         JOIN trait ON trait.rowid = trait_content.trait
         WHERE trait.name = '%s'
         ORDER BY trait_content.node",
        x
      )
    ),
  con = net$con
) -> sim_trait

## Showing the sequences and traits evolving is even less straightforward for a
## reticulated phylogenetic network as it is for a phylogenetic tree, let alone
## for a linear sequence of descendants.

## We first need to retrieve the edge list as follows (same as for a tree):
dbGetQuery(
  net$con,
  "SELECT i,j FROM edge ORDER BY rowid"
) -> edge

## From that list, we can determine which node is a leaf (also same as for a
## tree):
leaf <- unique(edge$j[!(edge$j %in% unique(edge$i))])
leaf

## Obtain the distances between the root node and all the nodes:
rbind(
  list(1L,0),
  dbGetQuery(
    net$con,
    "SELECT j, d FROM diss WHERE i = 1 ORDER BY j"
  )
) -> dst

## To obtain simple lineage, we need a criterion to choose which parent to
## follow at hybridization events. Here, the parent that is the closest to the
## root is the one chosen:
getLineageDR <- function(x, e, d) {
  while(length(wh <- which(e$j == head(x,1)))) {
    wh <- wh[which.min(d$d[match(e$i[wh],d$j)])]
    x <- c(e$i[wh], x)
  }
  unname(x)
}

## Here is the lineage of the first leaf of the network:
getLineageDR(leaf[1], edge, dst)

## Show the sequences of the lineage of the first leaf:
net %>%
  getSequence(name = "SEQ1") %>%
  .[getLineageDR(leaf[1], edge, dst),] %>%
  concatenate(discard = "-") %>%
  show.sequence

## Show the sequences of the lineage of the second leaf:
net %>%
  getSequence(name = "SEQ1") %>%
  .[getLineageDR(leaf[2], edge, dst),] %>%
  concatenate(discard = "-") %>%
  show.sequence

## Show the sequences of the lineage of the third leaf:
net %>%
  getSequence(name = "SEQ1") %>%
  .[getLineageDR(leaf[6], edge, dst),] %>%
  concatenate(discard = "-") %>%
  show.sequence


## Function to show every lineages 'l':
plotit <- function(l) {
  nn <- getLineageDR(leaf[l], edge, dst)
  dd <- dst$d[match(nn, dst$j)]
  
  ## Plot the trait simulation results:
  par(mar=c(4,4,1,1))
  plot(NA, xlim=range(dd), xlab="Time", ylab="Trait value", las=1,
       ylim=range(-25,sapply(sim_trait, function(x) x$value),80))
  col <- c("black","red","blue","green")
  
  ## The solid line is the trait value, whereas the dotted line is the effective
  ## trait optimum:
  for(i in 1:4) {
    lines(x=dd, y=sim_trait[[i]]$value[match(nn, dst$j)], col=col[i])
    if(!all(is.na(sim_trait[[i]]$optim)))
      lines(x=dd, y=sim_trait[[i]]$optim[match(nn, dst$j)], col=col[i], lty=3)
  }
  ## Note: there are no effective optimum for a trait evolving neutrally.
}

## Plot the lineages of some of the leaves:
plotit(7)
plotit(15)
plotit(22)
plotit(length(leaf))

## Clear the simulated linear network:
clearNetwork(net, TRUE)
