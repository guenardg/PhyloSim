## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Trait Evolution Simulation - SQLite TSSI Functions **
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
#' A database interface for simulating traits along evolutionary time on linear
#' evolutionary series, phylogenerit trees, or phylogenetic networks.
#' 
#' @name SQLite-TSSI
#' 
#' @aliases openNetwork closeNetwork makeLinearNetwork makeTreeNetwork
#' makeReticulatedNetwork simulateSequence simulateTrait drawDNASequence
#' drawEvolRate
#' 
#' @param filename A character string; the name of the SQLite database file,
#' without the ".sim" suffix. The default depends on the function (see details).
#' @param load A logical; should the database specified using argument
#' \code{filename} be loaded in the computer memory instead of accessed from the
#' file. The default is \code{FALSE}: access from the file.
#' @param save A logical; when the SQLite database has been created in the
#' computer memory, whether or not to save it before to close it. It is ignored
#' for an SQLite database stored in a file. The default is \code{TRUE}
#' @param overwrite A logical; whether or not to erase an existing database file
#' with the same name as the one specified by argument \code{filename} before
#' saving an SQLite database that has been created in the computer memory. The
#' default is \code{FALSE}. See detail section for more details.
#' @param net A two elements list: the name of the network and an SQLite
#' database connection from function \code{\link[DBI]{dbConnect}}.
#' @param NS An integer; the number of sequences to generate.
#' @param NC A function without arguments returning the simulated number of
#' children for each ancestor (see details).
#' @param timestep A function without arguments returning the value (numeric)
#' of the time step during the simulation (see details).
#' @param root A logical; whether or not to root the series (default:
#' \code{TRUE}). It is carried out on new data series (a single lineage, tree,
#' or network) and set to \code{FALSE} when extending an existing data series
#' (see details).
#' @param verbose A logical; whether or not to print details about the
#' simulation process during the calculations (default: \code{FALSE}).
#' @param NP A function without arguments returning the number of parents during
#' a hybridization event.
#' @param maxDiss A numeric; the maximum distance, in term of the value of
#' argument \code{timestep}, beyond which hybridization is disallowed to occur.
#' @param contrib A function taking the number of parents as its one argument,
#' and returning a numeric vector of the same length as the number of parents
#' giving the relative contribution of these parents during an hybridization
#' event (see details).
#' @param name A character string; the name of the DNA sequence or trait to be
#' simulated or retrieved.
#' @param NN An integer; the number of locations to be generated. A locations is
#' either one of the four nucleotides or a gap.
#' @param prob A length 5 numeric vector with a sum of 1; probabilities for
#' drawing a gap (-) or one of the for DNA bases (A, C, G, or T), in that order.
#' Default is \code{c(0.3,0.175,0.175,0.175,0.175)}. Probabilities not summing
#' to 1 are made to sum to 1.
#' @param gamma.shape A numeric; shape parameter of the beta distribution used
#' to draw the nucleotide (and gaps) evolution rates.
#' @param gamma.scale A numeric; scale parameter of the beta distribution used
#' to draw the nucleotide (and gaps) evolution rates.
#' @param I A 5 x 5 numeric matrix; transition intensity matrix (see
#' \code{\link[PhyloSim]{molEvolMod}} for details).
#' @param sqn A raw vector; the seed sequence used at the root of the network.
#' @param rate A numeric vector of length \code{NN}; values of nucleotide
#' evolution rate.
#' @param tem A trait evolution model; see \code{\link[PhyloSim]{traitEvolMod}}
#' for the details.
#' @param state An integer; for a trait evolving according to an
#' Ornstein-Uhlenbeck process, index of the trait state at the onset of the
#' simulation.
#' @param value A numeric; value of the trait at the onset of the simulation.
#' @param note A character string; note about the trait to be stored in the
#' SQLite database.
#' @param removeGapOnly A logical; whether or not to remove positions that are
#' all gaps from the output (default: \code{TRUE}).
#' 
#' @details The simulation workflow involves the following xxx steps:
#' 
#' Firstly, a network is opened, either in the computer memory or into a file.
#' Opening a database in the computer memory is much faster than opening it into
#' a file, but uses computer memory more heavily. A memory database is preferred
#' whenever possible. It is also possible to make part of the process in the
#' memory, then close the connection while saving into a file, and then opening
#' the same database as a file and perform the more memory-intensive tasks
#' directly into a file. That task is performed by function
#' \code{openNetwork()}. That function opens a database connection and creates
#' the necessary tables for a new database. When argument
#' \code{filename = NULL}, the database created resides in the computer memory,
#' whereas when argument \code{filename} is a valid and accessible file path,
#' the database created resides in the specified file. If the database file
#' already exists, it is opened and any of the tables required by the simulation
#' process is only created if it is not already present.
#' 
#' Secondly, a single simulated network is generated in the database using one
#' of three available function: \code{makeLinearNetwork},
#' \code{makeTreeNetwork}, and \code{makeReticulatedNetwork}, which are shown
#' here in increasing order of complexity. These functions share arguments
#' \code{net}, \code{NS}, \code{NC}, \code{timestep}, \code{root}, and
#' \code{verbose}.
#' 
#' Function \code{makeLinearNetwork()} creates a simple linear network by taking
#' the last node created as the parent of the new descendant node(s). While
#' more than one descendant per not is allowed, the function will produce a
#' consistent linear backbone of ancestors, somewhat similar to a minimum
#' spanning tree, albeit not formally one.
#' 
#' Function \code{makeTreeNetwork()} has the very same arguments as
#' \code{makeLinearNetwork()}, but will extend the network randomly from any
#' ancestors that still can bear a descendant, resulting in a structure with
#' greater arborescence than a single linear backbone.
#' 
#' Function \code{makeReticulatedNetwork()} has three more arguments, namely
#' \code{NP}, \code{maxDiss}, and \code{contrib}, and also perform a more
#' complex task. It creates a reticulated network where descendants are allowed
#' to have multiple parents. The number of parents is drawn randomly by the
#' function given to argument \code{NP}. When the number of parent drawn is
#' larger than 1, a single ancestor is selected randomly as the first parent.
#' Than the remaining parents are drawn among the available ancestors within a
#' certain maximum hybridization distance (drawn by the function given to
#' argument \code{maxDiss}) from the first parent. Then, the relative
#' contributions of the parents are drawn by the function given to argument
#' \code{contrib}. This process results in a random directed graph with a
#' more or less reticulated structure depending on the choice of parameters.
#' 
#' Thirdly, it is possible to populate the network with any number of simulated
#' DNA sequences or simulated traits using functions \code{simulateSequence()}
#' and \code{simulateTrait()}, respectively.
#' 
#' Function \code{simulateSequence()} assigns a random sequence at the root node
#' of the network through its argument \code{sqn}, and also takes a set of
#' evolution rate (one for each location) through its argument \code{rate}.
#' The sequence also needs to be given a name (unique for each network; argument
#' \code{name}) and, optionally, a textual note about it can be added (argument
#' \code{note}). The DNA sequence evolution is simulated as a Markov process
#' described by the transition intensity matrix given as argument \code{I}. This
#' process is described in \code{\link[PhyloSim]{molEvolMod}}. 
#' 
#' Function \code{simulateTrait()} assigns a value to the trait (and,
#' optionally, the optimal trait state) at the root node of the network through
#' its argument \code{value} (and, optionally, \code{state}). As for
#' \code{simulateSequence()}, the trait also needs to be given a name (unique
#' for each network; argument \code{name}) and, optionally, a textual note about
#' the it can be added (argument \code{note}). Trait evolution is simulated as
#' a combination of the Ornstein-Uhlenbeck and Markov processes described in
#' \code{\link[PhyloSim]{traitEvolMod}}.
#' 
#' Fourthly, the sequences and traits are generated can be retrieved using
#' functions \code{getSequence()} and \code{getTrait()}.
#' 
#' Finally, the network can be closed using function \code{closeNetwork}. If it
#' has been residing in the computer memory, it is the occasion to copy it to a
#' file by providing a value to argument \code{filename}. It the file already
#' exists, \code{closeNetwork} will have to be allowed to overwrite it by
#' passing argument \code{overwrite = TRUE}. Otherwise, the function will end in
#' an error. If argument \code{save = FALSE} (the default is \code{TRUE}), the
#' the network will be discarded if it resides in the computer memory. Arguments
#' \code{filename}, \code{save}, and \code{overwrite} are not used when closing
#' a network residing in a file.
#' 
#' @return
#' \describe{
#'   \item{openNetwork}{A two elements list containing the name of the database
#'   and an SQLite database connection from function
#'   \code{\link[DBI]{dbConnect}}.}
#'   \item{closeNetwork}{\code{NULL} (invisibly).}
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
#'   \item{getTrait}{A three-column data frame with the node numbers, optimum
#'   values and trait values.}
#' }
#' 
#' @author \packageAuthor{PhyloSim}
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
#' ## Open an empty data base in a file to store the network (slow):
#' net <- openNetwork(filename = "net_linear.sim")
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
#' ## One can obtain the trait values at the nodes as follows:
#' lapply(
#'   sprintf("trait%d",1:4),
#'   function(x)
#'     getTrait(net, x)
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
#' ## Close the simulated linear network while saving the data:
#' closeNetwork(net, filename = "net_linear.sim")
#' 
#' ## Reopen in the computer memory (load = TRUE):
#' net <- openNetwork(filename = "net_linear.sim", load = TRUE)
#' 
#' ## Close the network for good:
#' closeNetwork(net, save = FALSE)
#' 
#' 
#' ### Tree case:
#' 
#' ## Initializing an empty data base in memory to store the network (fast):
#' net <- openNetwork()
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
#' ## One can obtain the trait values at the nodes as follows:
#' lapply(
#'   sprintf("trait%d",1:4),
#'   function(x)
#'     getTrait(net, x)
#' ) -> sim_trait
#' 
#' ## Showing the sequences and traits evolving is not as straightforward for a
#' ## phylogenetic tree as it is for a linear sequence of descendants.
#' 
#' ## We can determine which node is a leaf as follows:
#' leaf <- getLeaves(net)
#' leaf
#' 
#' ## The lineage of each leaf is obtained as follows:
#' lineages <- getLineage(net, leaf)
#' 
#' ## Here is the lineage of the first leaf:
#' lineages[[1L]]
#' 
#' ## Show the first sequence (here, the y-axis corresponds to the time step):
#' net %>%
#' getSequence(name = "SEQ1") %>%
#'   .[lineages[[1L]],] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#'   
#' ## Show the second sequence:
#' net %>%
#' getSequence(name = "SEQ2") %>%
#'   .[lineages[[2L]],] %>%
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
#' plotit <- function(nn) {
#'   
#'   dd <- dst$d[match(nn, dst$j)]
#'   
#'   ## Plot the trait simulation results:
#'   par(mar=c(4,4,1,1))
#'   plot(NA, xlim=range(dd), xlab="Time", ylab="Trait value", las=1,
#'        ylim=range(-25,sapply(sim_trait, function(x) x$value),80))
#'   col <- c("black","red","blue","green")
#'   
#'   ## The solid line is the trait value, whereas the dotted line is the effective
#'   ## trait optimum:
#'   for(i in 1:4) {
#'     lines(x=dd, y=sim_trait[[i]]$value[match(nn, dst$j)], col=col[i])
#'     if(!all(is.na(sim_trait[[i]]$optim)))
#'       lines(x=dd, y=sim_trait[[i]]$optim[match(nn, dst$j)], col=col[i], lty=3)
#'   }
#'   ## Note: there are no effective optimum for a trait evolving neutrally.
#' }
#' 
#' ## Plot the lineages of some of the leaves:
#' plotit(lineages[[1]])
#' plotit(lineages[[2]])
#' plotit(lineages[[3]])
#' plotit(lineages[[6]])
#' plotit(lineages[[7]])
#' 
#' ## Clear the simulated linear network:
#' closeNetwork(net)
#' 
#' 
#' ### The reticulated case:
#' 
#' ## Open an empty data base in computer memory to store the network (fast):
#' net <- openNetwork()
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
#' ## One can obtain the trait values at the nodes as follows:
#' lapply(
#'   sprintf("trait%d",1:4),
#'   function(x)
#'     getTrait(net, x)
#' ) -> sim_trait
#' 
#' ## Showing the sequences and traits evolving is even less straightforward for
#' ## a reticulated phylogenetic network as it is for a phylogenetic tree, let
#' ## alone for a linear sequence of descendants.
#' 
#' ## We can determine which node is a leaf as follows:
#' leaf <- getLeaves(net)
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
#' ## The lineage of each leaf is obtained as follows:
#' lineages <- getLineage(net, leaf)
#' 
#' ## Here is the lineage of the first leaf of the network:
#' lineages[[1L]]
#' 
#' ## Show the sequences of the lineage of the first leaf:
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   .[lineages[[1L]],] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Show the sequences of the lineage of the second leaf:
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   .[lineages[[2L]],] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Show the sequences of the lineage of the third leaf:
#' net %>%
#'   getSequence(name = "SEQ1") %>%
#'   .[lineages[[3L]],] %>%
#'   concatenate(discard = "-") %>%
#'   show.sequence
#' 
#' ## Function to show every lineages 'l':
#' plotit <- function(nn) {
#'   
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
#' plotit(lineages[[7]])
#' plotit(lineages[[15]])
#' plotit(lineages[[22]])
#' plotit(unlist(tail(lineages,1)))
#' 
#' ## Close the simulated reticulated network, discarding it:
#' closeNetwork(net, save = FALSE)
#' 
#' ## Erase the database files:
#' file.remove(sprintf("net_%s.sim", c("linear","tree","reticulated")))
#' 
#' @importFrom DBI dbConnect dbDisconnect dbListTables
#' @importFrom RSQLite SQLite sqliteCopyDatabase
#' @importFrom stats rgamma
#' 
NULL
#' 
#' @describeIn SQLite-TSSI
#' 
#' Open Network
#' 
#' Opens an SQLite database in which to simulate various evolutionary series.
#' 
#' @export
openNetwork <- function(filename = NULL, load = FALSE) {
  
  yaml.load_file(
    system.file(
      package = "PhyloSim",
      "extdata",
      "SQLiteAPI.yml"
    )
  ) -> cfg
  
  net <- list(name = if(!load) filename else NULL)
  
  if(is.null(filename)) {
    net$con <- dbConnect(SQLite(), ":memory:")
  } else if(load) {
    con <- dbConnect(SQLite(), filename)
    net$con <- dbConnect(SQLite(), ":memory:")
    sqliteCopyDatabase(con, net$con)
    dbDisconnect(con)
  } else
    net$con <- dbConnect(SQLite(), filename)
  
  tbl <- dbListTables(net$con)
  
  lapply(
    cfg[!(names(cfg) %in% tbl)],
    function(x) dbExecute(net$con, x$create)
  )
  
  net
}
#' 
#' @describeIn SQLite-TSSI
#' 
#' Close Network
#' 
#' Closes an SQLite database, with some saving options.
#' 
#' @export
closeNetwork <- function(net, save = TRUE, overwrite = FALSE,
                         filename = "default") {
  
  if(is.null(net$name)) {
    
    if(save) {
      
      if(file.exists(filename)) {
        if(overwrite) {
          file.remove(filename)
        } else
          stop(
            "Cannot save because file '",filename,
            "' exists and argument overwrite = FALSE"
          )
      }
      
      con <- dbConnect(SQLite(), filename)
      sqliteCopyDatabase(net$con, con)
      dbDisconnect(con)
      dbDisconnect(net$con)
      
    } else
      dbDisconnect(net$con)
    
  } else
    dbDisconnect(net$con)
  
  invisible(NULL)
}
#' 
#' @describeIn SQLite-TSSI
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
#' @describeIn SQLite-TSSI
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
#' @describeIn SQLite-TSSI
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
#' @describeIn SQLite-TSSI
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
#' @describeIn SQLite-TSSI
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
#' @describeIn SQLite-TSSI
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
#' @describeIn SQLite-TSSI
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
#' @describeIn SQLite-TSSI
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
#' 
#' @describeIn SQLite-TSSI
#' 
#' Trait Exractor
#' 
#' Extracts all the nodal optima and values of a given trait from a phylogenetic
#' network.
#' 
#' @export
getTrait <- function(net, name) {
  
  if(!(name %in% unlist(dbGetQuery(net$con,"SELECT name FROM trait"))))
    stop("There is no trait named '", name, "' associated with this network")
  
  dbGetQuery(
    net$con,
    sprintf(
      "SELECT trait_content.node, trait_content.optim, trait_content.value
       FROM trait_content
       JOIN trait ON trait_content.trait = trait.rowid
       WHERE trait.name = '%s'
       ORDER BY trait_content.node",
      name
    )
  )
}
