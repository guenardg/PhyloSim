## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Trait Evolution Simulation - SQLite TSI Functions **
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
#' SQLite Trait Simulation Interface
#' 
#' A database interface for simulating traits along evolutionaly time on linear
#' evolutionary series, phylogenerit trees, or phylogenetic networks.
#' 
#' @name SQLite-TSI
#' 
#' @aliases connectNetwork ClearNetwork makeLinearNetwork makeTreeNetwork
#' makeReticulatedNetwork
#' 
#' @param name The name of the SQLite database, sequence, or trait to be created
#' (character).
#' @param init Whether or not to initialize the SQLite database upon creation
#' (logical).
#' @param net A two elements list containing the name of the database an a
#' SQLite database connection from function \code{\link[DBI]{dbConnect}}.
#' @param delete Whether to delete the database file upon closing the database.
#' @param NS The number of sequences to generate (integer).
#' @param NC A function without arguments returning the simulated number of
#' children for each ancestor (see details).
#' @param timestep A function (without arguments) returning the value (numeric)
#' of the time step during the simulation (see details).
#' @param root Whether to root the series or not (logical; default:
#' \code{TRUE}). It is carried out on new data series (a single lineage, tree,
#' or network) and set to \code{FALSE} when extending an existing data series
#' (see details).
#' @param verbose Whether to print details about the simulation process during
#' the calculations (logical; default: \code{FALSE}).
#' @param NP A function without arguments returning the number of parents during
#' a hybridization event.
#' @param maxDiss The maximum distance (in term of the value of argument
#' \code{timestep}.)
#' @param contrib A function returning the relative contribution of the
#' different parents during an hybridization event.
#' @param NN The number of locations to be generated (integer). A locations is
#' either one of the four nucleotides or a gap).
#' @param prob Probabilities for drawing a gap (-) or one of the for DNA bases
#' (A, C, G, or T), in that order (A length 5 numeric vector). Default:
#' \code{c(0.3,0.175,0.175,0.175,0.175)}. Probabilities not summing to 1 will be
#' made to sum to 1.
#' @param gamma.shape Shape parameter of the beta distribution used to draw the
#' nucleotide (and gaps) evolution rates (numeric).
#' @param gamma.scale Scale parameter of the beta distribution used to draw the
#' nucleotide (and gaps) evolution rates (numeric).
#' @param I A 5 x 5 transition intensity matrix (see
#' \code{\link[PhyloSim]{molEvolMod}} for details).
#' @param sqn A raw vector containing the seed sequence used at the root of the
#' network.
#' @param rate A vector of nucleotide evolution rates (numeric).
#' @param tem A trait evolution model (see \code{\link[PhyloSim]{traitEvolMod}}
#' for the details.)
#' @param state For a trait evolving according to an Ornstein Uhlenbeck process,
#' trait state at the onset of the simulation.
#' @param value Value of the trait at the onset of the simulation.
#' @param note Note about the trait to be stored in the SQLite database.
#' @param removeGapOnly Remove positions that are all gaps from the output
#' (logical; default: \code{TRUE}).
#' 
#' @details Details here...
#' 
#' @return
#' \describe{
#'   \item{connectNetwork}{A two elements list containing the name of the
#'   database an a SQLite database connection from function
#'   \code{\link[DBI]{dbConnect}}.}
#'   \item{clearNetwork}{\code{NULL} (invisibly).}
#'   \item{makeLinearNetwork}{The number of the last node that has been
#'   created.}
#'   \item{makeTreeNetwork}{The number of the last node that has been created.}
#'   \item{makeReticulatedNetwork}{The number of the last node that has been
#'   created.}
#'   \item{drawDNASequence}{A raw vector of length \code{NN} containing the
#'   ASCII values for characters '-' (0x2d), 'A' (0x41), 'C' (0x43), 'G' (0x47),
#'   and 'T' (0x54) representing random nucleotides to be used as the seed
#'   sequence for a DNA fragment.}
#'   \item{drawEvolRate}{A numeric vectors of length \code{NN} containing the
#'   evolution rate for each of the nucleotides in the sequence.}
#'   \item{simulateSequence}{\code{NULL} (invisibly).}
#'   \item{simulateTrait}{\code{NULL} (invisibly).}
#'   \item{getSequence}{A vector of type character.}
#' }
#' 
#' @examples ## Load the example of a configuration file provided with the
#' ## package:
#' system.file(
#'   package = "PhyloSim",
#'   "extdata",
#'   "evolMod.yml"
#' ) %>%
#'   yaml.load_file -> dnaParam
#' 
#' ## The configuration list is located in member `$DNA`
#' dnaParam$DNA %>%
#'   read.mutationMat -> I
#' 
#' ## The transition intensity matrix:
#' I
#' 
#' ## The trait evolution model parameters:
#' dnaParam$trait %>%
#'   lapply(read.TraitEvolMod) -> traitMod
#' 
#' 
#' ### Linear case:
#' 
#' ## Initializing a new data base to store the network:
#' net <- connectNetwork(name = "net_linear", init = TRUE)
#' 
#' ## Setting RNG seed to have a consistent example:
#' set.seed(162745)
#' 
#' ## These are the tables that are created in the database:
#' dbListTables(net$con)
#' 
#' ## The following code enables one to extract the SQLite table information:
#' lapply(
#'   dbListTables(net$con),
#'   function(x, con)
#'     dbGetQuery(con, sprintf("PRAGMA table_info(%s)", x)),
#'   con = net$con
#' ) -> table_info
#' 
#' names(table_info) <- dbListTables(net$con)
#' 
#' table_info
#' 
#' ## The following code retrieves the number of records in all the tables.
#' ## All of them are empty at the beginning of the simulation:
#' sapply(
#'   dbListTables(net$con),
#'   function(x, con)
#'     unlist(
#'       dbGetQuery(
#'         con,
#'         sprintf(
#'           "SELECT COUNT(*) AS %s
#'          FROM %s",
#'           x,x),
#'       )
#'     ),
#'   con = net$con
#' )
#' 
#' ## Create a linear network with 50 species having 1 descendant each with a
#' ## constant time step of 1:
#' makeLinearNetwork(
#'   net = net,
#'   NS = 50,
#'   NC = function() 1,
#'   timestep = function() 1,
#'   root = TRUE,
#'   verbose = TRUE
#' )
#' 
#' ## After the linear phylogenetic network has been created, tables diss, edge,
#' ## and node conntain information. The other tables are used for storing the
#' ## sequence and trait information.
#' sapply(
#'   dbListTables(net$con),
#'   function(x, con)
#'     unlist(
#'       dbGetQuery(
#'         con,
#'         sprintf("SELECT COUNT(*) FROM %s", x),
#'       )
#'     ),
#'   con = net$con
#' )
#' 
#' ## Simulate two sequences, one named 'SEQ1' having 100 locations (i.e., one
#' ## of the four nucleotide or a gap) and a second sequence named 'SEQ2' having
#' ## 250 locations.
#' data.frame(
#'   name = c("SEQ1","SEQ2"),
#'   NN = c(100,250)
#' ) -> cond
#' 
#' ## Calling function simulateSequence once for each of the simulated
#' ## sequences, each time drawing a random root sequence (using
#' ## drawDNASequence) and a set of evolution rate (using function
#' ## drawEvolRate).
#' for(i in 1:NROW(cond))
#'   simulateSequence(
#'     net,
#'     I,
#'     name = cond$name[i],
#'     sqn = drawDNASequence(cond$NN[i]),
#'     rate = drawEvolRate(cond$NN[i])
#'   )
#' 
#' ## Now, tables seq and seq_content contain information:
#' sapply(
#'   dbListTables(net$con),
#'   function(x, con)
#'     unlist(
#'       dbGetQuery(
#'         con,
#'         sprintf("SELECT COUNT(*) FROM %s", x),
#'       )
#'     ),
#'   con = net$con
#' )
#' 
#' ## Show the first sequence (here, the y-axis corresponds to the time step):
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Show the second sequence:
#' net %>%
#'   getSequence(name = "SEQ2") %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Set four traits to be simulated.
#' data.frame(
#'   name = sprintf("trait%d",1:4),
#'   state = c(2,NA,1,1),
#'   value = c(50,0,15,-25),
#'   note = c("Non-neutral (OU) with multiple optima","Neutral",
#'            "Non-neutral (OU) with single optimum",
#'            "Non-neutral (OU) with single optimum")
#' ) -> cond
#' 
#' ## The four traits are simulated as follows:
#' for(i in 1:NROW(cond))
#'   simulateTrait(
#'     net,
#'     tem = traitEvolMod(traitMod[[i]]),
#'     name = cond$name[i],
#'     state = cond$state[i],
#'     value = cond$value[i],
#'     note = cond$note[i]
#'   )
#' 
#' ## Now, tables trait and trait_content also contain information:
#' sapply(
#'   dbListTables(net$con),
#'   function(x, con)
#'     unlist(
#'       dbGetQuery(
#'         con,
#'         sprintf("SELECT COUNT(*) FROM %s", x),
#'       )
#'     ),
#'   con = net$con
#' )
#' 
#' ## To obtain the trait values at the nodes, one has to use an SQL query (in
#' ## the SQLite dialect) as follows:
#' lapply(
#'   sprintf("trait%d",1:4),
#'   function(x, con)
#'     dbGetQuery(
#'       con,
#'       sprintf(
#'         "SELECT trait_content.node, trait_content.optim, trait_content.value
#'          FROM trait_content
#'          JOIN trait ON trait.rowid = trait_content.trait
#'          WHERE trait.name = '%s'
#'          ORDER BY trait_content.node",
#'         x
#'       )
#'     ),
#'   con = net$con
#' ) -> sim_trait
#' 
#' ## Plot the trait simulation results:
#' par(mar=c(4,4,1,1))
#' plot(NA, xlim=c(0,50), xlab="Time", ylab="Trait value", las=1,
#'      ylim=range(-25,sapply(sim_trait, function(x) x$value),80))
#' col <- c("black","red","blue","green")
#' 
#' ## The solid line is the trait value, whereas the dotted line is the
#' ## effective trait optimum:
#' for(i in 1:4) {
#'   lines(x=0:49, y=sim_trait[[i]]$value, col=col[i])
#'   if(!all(is.na(sim_trait[[i]]$optim)))
#'     lines(x=0:49, y=sim_trait[[i]]$optim, col=col[i], lty=3)
#' }
#' ## Note: there are no effective optimum for a trait evolving neutrally.
#' 
#' ## Clear the simulated linear network:
#' clearNetwork(net, TRUE)
#' 
#' 
#' ### Tree case:
#' 
#' ## Initializing a new data base to store the network:
#' net <- connectNetwork(name = "net_linear", init = TRUE)
#' 
#' ## Setting RNG seed to have a consistent example:
#' set.seed(162745)
#' 
#' ## Create a tree network (phylogenetic tree) with 100 species having
#' ## a fix number of 2 descendants (dichotomic) and a random log-normal time
#' ## step drawn at each time step:
#' makeTreeNetwork(
#'   net = net,
#'   NS = 100,
#'   NC = function() 2,
#'   timestep = function() exp(rnorm(1,0,0.5)),
#'   root = TRUE,
#'   verbose = TRUE
#' )
#' 
#' ## Simulate two sequences, one named 'SEQ1' having 50 locations (i.e., one of
#' ## the four nucleotide or a gap) and a second sequence named 'SEQ2' having 80
#' ## locations.
#' data.frame(
#'   name = c("SEQ1","SEQ2"),
#'   NN = c(50,80)
#' ) -> cond
#' 
#' ## Calling function simulateSequence once for each of the simulated
#' ## sequences, each time drawing a random root sequence (using
#' ## drawDNASequence) and a set of evolution rate (using function
#' ## drawEvolRate).
#' for(i in 1:NROW(cond))
#'   simulateSequence(
#'     net,
#'     I,
#'     name = cond$name[i],
#'     sqn = drawDNASequence(cond$NN[i]),
#'     rate = 10*drawEvolRate(cond$NN[i])  ## Accelerate evolution by 10 folds
#'   )
#' 
#' ## Set four traits to be simulated.
#' data.frame(
#'   name = sprintf("trait%d",1:4),
#'   state = c(2,NA,1,1),
#'   value = c(50,0,15,-25),
#'   note = c("Non-neutral (OU) with multiple optima","Neutral",
#'            "Non-neutral (OU) with single optimum",
#'            "Non-neutral (OU) with single optimum")
#' ) -> cond
#' 
#' ## The four traits are simulated as follows:
#' for(i in 1:NROW(cond))
#'   simulateTrait(
#'     net,
#'     tem = traitEvolMod(traitMod[[i]]),
#'     name = cond$name[i],
#'     state = cond$state[i],
#'     value = cond$value[i],
#'     note = cond$note[i]
#'   )
#' 
#' ## To obtain the trait values at the nodes, one has to use an SQL query (in
#' ## the SQLite dialect) as follows:
#' lapply(
#'   sprintf("trait%d",1:4),
#'   function(x, con)
#'     dbGetQuery(
#'       con,
#'       sprintf(
#'         "SELECT trait_content.node, trait_content.optim, trait_content.value
#'          FROM trait_content
#'          JOIN trait ON trait.rowid = trait_content.trait
#'          WHERE trait.name = '%s'
#'          ORDER BY trait_content.node",
#'         x
#'       )
#'     ),
#'   con = net$con
#' ) -> sim_trait
#' 
#' ## Showing the sequences and traits evolving is not as straightforward for a
#' ## phylogenetic tree as it is for a linear sequence of descendants.
#' 
#' ## We first need to retrieve the edge list as follows:
#' dbGetQuery(
#'   net$con,
#'   "SELECT i,j FROM edge ORDER BY rowid"
#' ) -> edge
#' 
#' ## From that list, we can determine which node is a leaf:
#' leaf <- edge$j[!(edge$j %in% unique(edge$i))]
#' leaf
#' 
#' ## We can use this simple function to obtain lineage of any node:
#' getLineage <- function(x, e) {
#'   while(length(wh <- which(e$j == head(x,1))))
#'     x <- c(e$i[wh], x)
#'   unname(x)
#' }
#' 
#' ## Here is the lineage of the first leaf:
#' getLineage(leaf[1], edge)
#' 
#' ## Show the first sequence (here, the y-axis corresponds to the time step):
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   .[getLineage(leaf[1], edge),] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Show the second sequence:
#' net %>%
#'   getSequence(name = "SEQ2") %>%
#'   .[getLineage(leaf[1], edge),] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Obtain the distances between the root node and all the nodes:
#' rbind(
#'   list(1L,0),
#'   dbGetQuery(
#'     net$con,
#'     "SELECT j, d FROM diss WHERE i = 1 ORDER BY j"
#'   )
#' ) -> dst
#' 
#' ## Function to show every lineages 'l':
#' plotit <- function(l) {
#'   nn <- getLineage(leaf[l], edge)
#'   dd <- dst$d[match(nn, dst$j)]
#'   
#'   ## Plot the trait simulation results:
#'   par(mar=c(4,4,1,1))
#'   plot(NA, xlim=range(dd), xlab="Time", ylab="Trait value", las=1,
#'        ylim=range(-25,sapply(sim_trait, function(x) x$value),80))
#'   col <- c("black","red","blue","green")
#'   
#'   ## The solid line is the trait value, whereas the dotted line is the
#'   ## effective trait optimum:
#'   for(i in 1:4) {
#'     lines(x=dd, y=sim_trait[[i]]$value[match(nn, dst$j)], col=col[i])
#'     if(!all(is.na(sim_trait[[i]]$optim)))
#'       lines(x=dd, y=sim_trait[[i]]$optim[match(nn, dst$j)], col=col[i],
#'             lty=3)
#'   }
#'   ## Note: there are no effective optimum for a trait evolving neutrally.
#' }
#' 
#' ## Plot the lineages of some of the leaves:
#' plotit(1)
#' plotit(2)
#' plotit(3)
#' plotit(6)
#' plotit(7)
#' 
#' ## Clear the simulated linear network:
#' clearNetwork(net, TRUE)
#' 
#' 
#' ### The reticulated case:
#' 
#' ## Initializing a new data base to store the network:
#' net <- connectNetwork(name = "net_reticulated", init = TRUE)
#' ## clearNetwork(net, TRUE)
#' 
#' ## Setting RNG seed to have a consistent example:
#' set.seed(162745)
#' 
#' ## Create a reticulated network with 100 species having a random number of
#' ## descendants and parents per node, a random log-normal time step, a random
#' ## maximum hybridization distance, and equal parental contributions on
#' ## hybridization events.
#' makeReticulatedNetwork(
#'   net,
#'   NS = 100,
#'   NC = function() 1 + rpois(1,2),
#'   NP = function() 1 + rpois(1,1),
#'   timestep = function() exp(rnorm(1,0,0.5)),
#'   maxDiss = function() runif(1, 2, 3),
#'   contrib = function(x) rep(1/x,x),
#'   root = TRUE,
#'   verbose = TRUE
#' )
#' 
#' ## Simulate one sequence named 'SEQ1' having 50 locations (i.e., one of
#' ## the four nucleotide or a gap). Calling function simulateSequence drawing a
#' ## random root sequence (using drawDNASequence) and a set of evolution rate
#' ## (using function drawEvolRate).
#' simulateSequence(
#'   net,
#'   I,
#'   name = "SEQ1",
#'   sqn = drawDNASequence(50),
#'   rate = 10*drawEvolRate(50)  ## Accelerate evolution by 10 folds
#' )
#' 
#' ## Set four traits to be simulated.
#' data.frame(
#'   name = sprintf("trait%d",1:4),
#'   state = c(2,NA,1,1),
#'   value = c(50,0,15,-25),
#'   note = c("Non-neutral (OU) with multiple optima","Neutral",
#'            "Non-neutral (OU) with single optimum",
#'            "Non-neutral (OU) with single optimum")
#' ) -> cond
#' 
#' ## The four traits are simulated as follows:
#' for(i in 1:NROW(cond))
#'   simulateTrait(
#'     net,
#'     tem = traitEvolMod(traitMod[[i]]),
#'     name = cond$name[i],
#'     state = cond$state[i],
#'     value = cond$value[i],
#'     note = cond$note[i]
#'   )
#' 
#' ## To obtain the trait values at the nodes, one has to use an SQL query (in
#' ## the SQLite dialect) as follows:
#' lapply(
#'   sprintf("trait%d",1:4),
#'   function(x, con)
#'     dbGetQuery(
#'       con,
#'       sprintf(
#'         "SELECT trait_content.node, trait_content.optim, trait_content.value
#'          FROM trait_content
#'          JOIN trait ON trait.rowid = trait_content.trait
#'          WHERE trait.name = '%s'
#'          ORDER BY trait_content.node",
#'         x
#'       )
#'     ),
#'   con = net$con
#' ) -> sim_trait
#' 
#' ## Showing the sequences and traits evolving is even less straightforward for
#' ## a reticulated phylogenetic network as it is for a phylogenetic tree, let
#' ## alone for a linear sequence of descendants.
#' 
#' ## We first need to retrieve the edge list as follows (same as for a tree):
#' dbGetQuery(
#'   net$con,
#'   "SELECT i,j FROM edge ORDER BY rowid"
#' ) -> edge
#' 
#' ## From that list, we can determine which node is a leaf (also same as for a
#' ## tree):
#' leaf <- unique(edge$j[!(edge$j %in% unique(edge$i))])
#' leaf
#' 
#' ## Obtain the distances between the root node and all the nodes:
#' rbind(
#'   list(1L,0),
#'   dbGetQuery(
#'     net$con,
#'     "SELECT j, d FROM diss WHERE i = 1 ORDER BY j"
#'   )
#' ) -> dst
#' 
#' ## To obtain simple lineage, we need a criterion to choose which parent to
#' ## follow at hybridization events. Here, the parent that is the closest to
#' ## the root is the one chosen:
#' getLineageDR <- function(x, e, d) {
#'   while(length(wh <- which(e$j == head(x,1)))) {
#'     wh <- wh[which.min(d$d[match(e$i[wh],d$j)])]
#'     x <- c(e$i[wh], x)
#'   }
#'   unname(x)
#' }
#' 
#' ## Here is the lineage of the first leaf of the network:
#' getLineageDR(leaf[1], edge, dst)
#' 
#' ## Show the sequences of the lineage of the first leaf:
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   .[getLineageDR(leaf[1], edge, dst),] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Show the sequences of the lineage of the second leaf:
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   .[getLineageDR(leaf[2], edge, dst),] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Show the sequences of the lineage of the third leaf:
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   .[getLineageDR(leaf[6], edge, dst),] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Function to show every lineages 'l':
#' plotit <- function(l) {
#'   nn <- getLineageDR(leaf[l], edge, dst)
#'   dd <- dst$d[match(nn, dst$j)]
#'   
#'   ## Plot the trait simulation results:
#'   par(mar=c(4,4,1,1))
#'   plot(NA, xlim=range(dd), xlab="Time", ylab="Trait value", las=1,
#'        ylim=range(-25,sapply(sim_trait, function(x) x$value),80))
#'   col <- c("black","red","blue","green")
#'   
#'   ## The solid line is the trait value, whereas the dotted line is the
#'   ## effective trait optimum:
#'   for(i in 1:4) {
#'     lines(x=dd, y=sim_trait[[i]]$value[match(nn, dst$j)], col=col[i])
#'     if(!all(is.na(sim_trait[[i]]$optim)))
#'       lines(x=dd, y=sim_trait[[i]]$optim[match(nn, dst$j)], col=col[i],
#'             lty=3)
#'   }
#'   ## Note: there are no effective optimum for a trait evolving neutrally.
#' }
#' 
#' ## Plot the lineages of some of the leaves:
#' plotit(7)
#' plotit(15)
#' plotit(22)
#' plotit(length(leaf))
#' 
#' ## Clear the simulated linear network:
#' clearNetwork(net, TRUE)
#' 
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom RSQLite SQLite
#' @importFrom stats rgamma
#' 
NULL
#' 
#' @describeIn SQLite-TSI
#' 
#' Create Network
#' 
#' Create an SQLite database in which to simulate various evolutionary series.
#' 
#' @export
connectNetwork <- function(name, init = FALSE) {
  
  net <- list(name=sprintf("%s.sqlite",name))
  
  net$con <- dbConnect(SQLite(), net$name)
  
  if(init) {
    yaml.load_file(
      system.file(
        package = "PhyloSim",
        "extdata",
        "SQLiteAPI.yml"
      )
    ) -> cfg
    
    cfg %>%
      lapply(
        function(x, net) dbExecute(net$con, x$create),
        net = net
      ) %>%
      unlist
  }
  
  net
}
#' 
#' @describeIn SQLite-TSI
#' 
#' Clear Network
#' 
#' Close an SQLite database, possibly deleting the database file.
#' 
#' @export
clearNetwork <- function(net, delete = FALSE) {
  
  dbDisconnect(net$con)
  
  if(delete)
    file.remove(net$name)
  
  invisible(NULL)
}
#' 
#' @describeIn SQLite-TSI
#' 
#' Linear Evolutionary Series
#' 
#' Simulates a linear evolutionary series.
#' 
#' @export
makeLinearNetwork <- function(net, NS, NC, timestep, root = TRUE,
                              verbose = FALSE) {
  
  if(root)
    node <- addNode(net, cld=NC())
  
  for(k in (if(root) 2L else 1L):NS) {
    asc <- tail(getAsc(net), 1L)
    node <- addNode(net, cld=NC())
    addEdge(net, asc, node, contrib=1L)
    decrementChildren(net, asc)
    tstp <- timestep()
    others <- (1L:node)[-c(asc,node)]
    if(length(others)) {
      dd <- getDiss(net, i=asc, j=others)
      addDiss(net, i=node, j=others, dd + tstp)
    }
    addDiss(net, i=asc, j=node, value=tstp)
    
    if(verbose)
      cat(sprintf("Link: %d -> %d (%f)\n", asc, node, tstp))
    
  }
  
  node
}
#' 
#' @describeIn SQLite-TSI
#' 
#' Phylogenetic Tree
#' 
#' Simulates a phylogenetic tree.
#' 
#' @export
makeTreeNetwork <- function(net, NS, NC, timestep, root = TRUE,
                            verbose = FALSE) {
  
  if(root)
    node <- addNode(net, cld=NC())
  
  for(k in (if(root) 2L else 1L):NS) {
    asc <- getAsc(net)
    asc <- if(length(asc) == 1L) asc else sample(asc,1L)
    node <- addNode(net, cld=NC())
    addEdge(net, asc, node, contrib=1L)
    decrementChildren(net, asc)
    tstp <- timestep()
    others <- (1L:node)[-c(asc,node)]
    if(length(others)) {
      dd <- getDiss(net, i=asc, j=others)
      addDiss(net, i=node, j=others, dd + tstp)
    }
    addDiss(net, i=asc, j=node, value=tstp)
    
    if(verbose)
      cat(sprintf("Link: %d -> %d (%f)\n", asc, node, tstp))
    
  }
  
  node
}
#' 
#' @describeIn SQLite-TSI
#' 
#' Phylogenetic Network
#' 
#' Simulates a phylogenetic network.
#' 
#' @export
makeReticulatedNetwork <- function(net, NS, NC, NP, timestep, maxDiss,
                                   contrib, root = TRUE, verbose = FALSE) {
  
  if(root)
    node <- addNode(net, cld=NC())
  
  for(k in (if(root) 2L else 1L):NS) {
    asc <- getAsc(net)
    if(length(asc) > 1L) {
      np <- NP()
      md <- maxDiss()
      s <- sample(asc,1L)
      r <- asc[!(asc %in% s)]
      r <- r[getDiss(net,s,r) <= md]
      asc <- c(s,if(length(r) > (np - 1L)) sample(r, np - 1L) else r)
    }
    np <- length(asc)
    node <- addNode(net, cld=NC())
    pasc <- contrib(np)
    tstp <- timestep()
    addEdge(net, asc, node, contrib=pasc)
    decrementChildren(net, asc)
    others <- (1L:node)[-c(asc,node)]
    if(length(others)) {
      dd <- pasc[1L]*getDiss(net,asc[1L],others)
      if(length(asc) > 1L)
        for(i in 2L:length(asc))
          dd <- dd + pasc[i]*getDiss(net,asc[i],others)
      addDiss(net, i=node, j=others, dd + tstp)
    }
    addDiss(net, i=asc, j=node, value=tstp)
    
    if(verbose) {
      cat(sprintf("%d links, ts = %f\n",length(asc), tstp))
      cat(sprintf("\tLink: %d -> %d (contribution = %f)\n", asc, node, pasc))
      cat("-----------------------\n")
    }
    
  }
  
  node
}
#' 
#' @describeIn SQLite-TSI
#' 
#' Random Sequence Generator
#' 
#' Generates a random sequence of gaps or DNA bases.
#' 
#' @export
drawDNASequence <- function(NN, prob = c(0.3, rep(0.175, 4L)))
  sample(
    x = charToRaw("-ACGT"),
    size = NN,
    prob = prob,
    replace = TRUE
  )
#' 
#' @describeIn SQLite-TSI
#' 
#' Random Evolution Rate Generator
#' 
#' Generates a set of random nucleotide evolution rate from a gamma
#' distribution.
#' 
#' @export
drawEvolRate <- function(NN, gamma.shape = 5, gamma.scale = 5e-04)
  rgamma(
    NN,
    shape = gamma.shape,
    scale = gamma.scale
  )
#' 
#' @describeIn SQLite-TSI
#' 
#' Sequence Evolution Simulator
#' 
#' Generates a set of DNA sequences along the edge of a phylogenetic network by
#' evolving filial sequences from parental ones following a random Markov
#' process.
#' 
#' @export
simulateSequence <- function(net, I, name, sqn, rate, note = "") {
  
  NN <- length(sqn)
  if(length(sqn) != length(rate)) {
    warning(
      "Length of 'rate' (",length(rate),") does not match that of 'sqn'",
      length(sqn),".")
    
    rate <- rep(rate, length.out=length(sqn))
  }
  
  
  dbGetQuery(
    net$con,
    sprintf("SELECT true FROM seq WHERE name = '%s'", name)
  ) %>%
    unlist -> q
  
  if(length(q))
    stop("Sequence name '",name,"' aleady exists.")
  
  dbExecute(
    net$con,
    sprintf("INSERT INTO seq(name,note) VALUES ('%s','%s')", name, note)
  )
  
  dbGetQuery(
    net$con,
    sprintf("SELECT rowid FROM seq WHERE name = '%s'", name)
  ) %>%
    unlist -> seqId
  
  dbGetQuery(
    net$con,
    "SELECT rowid FROM node"
  ) %>% unlist -> node
  
  timestep <- 1
  
  em <- list()
  for(i in 1L:NN)
    em[[i]] <- molEvolMod(I, timestep, rate[i])
  
  
  for(k in node) {
    if(k == 1L) {
      
      dbExecute(
        net$con,
        sprintf(
          "INSERT INTO seq_content(seq,node,content) VALUES (%d,%d,%s)",
          seqId, k, sprintf("X'%s'", paste(sqn, collapse=""))
        )
      )
      
    } else {
      
      dbGetQuery(
        net$con,
        sprintf("SELECT i, contrib FROM edge WHERE j = %d",k)
      ) -> asc
      
      ## getDiss(net,asc$i,k)
      nts <- getDiss(net,asc$i[1L],k)
      if(nts != timestep) {
        timestep <- nts
        for(i in 1L:NN)
          em[[i]]$recalculate(timestep, rate[i])
      }
      
      ## Select one of the ascendants:
      asc <- if(NROW(asc) == 1L) asc$i else sample(asc$i, 1L, prob=asc$contrib)
      
      dbGetQuery(
        net$con,
        sprintf(
          "SELECT content FROM seq_content
           JOIN seq ON seq.rowid = seq_content.seq
           WHERE seq.name = '%s' AND seq_content.node = %d",
          name, asc
        )
      ) %>%
        unlist -> sqn
      
      for(i in 1L:NN)
        sqn[i] <- em[[i]]$evolve(sqn[i])
      
      dbExecute(
        net$con,
        sprintf(
          "INSERT INTO seq_content (seq,node,content) VALUES (%d,%d,X'%s')",
          seqId, k, paste(sqn,collapse="")
        )
      )
      
    }
  }
  
  invisible(NULL)
}
#' 
#' @describeIn SQLite-TSI
#' 
#' Trait Evolution Simulator
#' 
#' Simulates the evolution of a quantitative traits following a random walk
#' process.
#' 
#' @export
simulateTrait <- function(net, tem, name, state, value, note = "") {
  
  dbGetQuery(
    net$con,
    sprintf("SELECT true FROM trait WHERE name = '%s'", name)
  ) %>%
    unlist -> q
  
  if(length(q))
    stop("Trait name '",name,"' aleady exists.")
  
  dbExecute(
    net$con,
    sprintf("INSERT INTO trait(name,note) VALUES ('%s','%s')", name, note)
  )
  
  dbGetQuery(
    net$con,
    sprintf("SELECT rowid FROM trait WHERE name = '%s'", name)
  ) %>%
    unlist -> traitId
  
  dbGetQuery(
    net$con,
    "SELECT rowid FROM node"
  ) %>% unlist -> node
  
  timestep <- 1
  
  for(k in node) {
    if(k == 1L) {
      
      optval <- tem$getOptima()[state]
      
      dbExecute(
        net$con,
        sprintf(
          "INSERT INTO trait_content(trait,node,state,optim,value)
           VALUES (%d,%d,%s,%s,%f)",
          traitId, k,
          if(is.na(state)) "null" else state,
          if(is.na(optval)) "null" else optval,
          value
        )
      )
      
    } else {
      
      dbGetQuery(
        net$con,
        sprintf("SELECT i, contrib FROM edge WHERE j = %d ORDER BY i", k)
      ) -> asc
      
      pasc <- asc$contrib
      
      timestep <- getDiss(net,asc$i[1L],k)
      
      if(timestep != tem$getStep())
        tem$setStep(timestep)
      
      dbGetQuery(
        net$con,
        sprintf(
          "SELECT node AS asc, state, optim, value FROM trait_content
           JOIN trait ON trait.rowid = trait_content.trait
           WHERE trait.name = '%s' AND trait_content.node IN (%s)
           ORDER BY node",
          name, paste(asc$i, collapse=",")
        )
      ) -> asc
      
      if(NROW(asc) > 1L) {
        state <- tem$updateState(sample(asc$state, 1L, prob = pasc))
      } else
        state <- tem$updateState(asc$state)
      
      optval <- tem$getOptima()[state]
      
      value <- tem$updateValue(sum(pasc*asc$value), state)
      
      dbExecute(
        net$con,
        sprintf(
          "INSERT INTO trait_content(trait,node,state,optim,value)
           VALUES (%d,%d,%s,%s,%f)",
          traitId, k,
          if(is.na(state)) "null" else state,
          if(is.na(optval)) "null" else optval,
          value
        )
      )
      
    }
  }
  
  invisible(NULL)
}
#' 
#' 
#' @describeIn SQLite-TSI
#' 
#' Sequence Exractor
#' 
#' Extracts all sequences of a given fragment from a phylogenetic network.
#' 
#' @export
getSequence <- function(net, name, removeGapOnly = TRUE) {
  dbGetQuery(
    net$con,
    sprintf(
      "SELECT content FROM seq_content
       JOIN seq ON seq.rowid = seq_content.seq
       WHERE seq.name = '%s'
       ORDER BY seq_content.node",
      name
    )
  ) -> sqn
  
  if(!NROW(sqn)) {
    warning("No content found for sequence name '",name,"'.")
    return(NULL)
  }
  
  sqn %<>% {matrix(unlist(.), nrow=length(.[[1L]]), byrow=TRUE)}
  
  if(removeGapOnly)
    sqn %<>% .[,apply(., 2L, function(x) !all(x == charToRaw("-")))]
  
  sqn
}
