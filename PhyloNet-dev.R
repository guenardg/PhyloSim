
##

library(magrittr)
library(yaml)
library(Rcpp)
library(PhyloSim)

## rm(list=ls())

source("PhyloNet-aux.R")
## sourceCpp("PhyloNet-dev.cpp")

I <- read.mutationMat("DNAEvol.yml")
## rowSums(I)

if(FALSE) {
  molEvolMod(I, 1, 0.001) -> em1
  
  ## em1$evolve("Z")
  
  em1$getMt()
  
  em1$evolve("-")
  em1$evolve("A")
  em1$evolve("C")
  em1$evolve("G")
  em1$evolve("T")
  
  em1$recalculate(1,0.01)
  
  em1$getMt()
  
  em1$evolve("-")
  em1$evolve("A")
  em1$evolve("C")
  em1$evolve("G")
  em1$evolve("T")
  
  rm(em1)
}

## Simulation

## Number of sequences, sequences length, time step, maximum number of children
## per parent node, and maximum patristic distance for hybridization (currently
## not implemented):
nseq <- 5000L
seql <- 20000L
tm <- 1
maxCld <- 2L
## maxHyb <- 10

## Draw the evolution rates for each nucleotides:
rgamma(
  seql,
  shape = 5,
  scale = 5e-4
) -> rate
## plot(density(rate))

## Create the nucleotides evolution models from the mutation intensity matrix
## 'I', the time step 'tm' and the evolution rates drawn previously:
em <- list()
for(i in 1L:seql)
  em[[i]] <- molEvolMod(I, tm, rate[i])
rm(i)

## Create the sequence storage matrix 'sqn' and the allowed child numbers 'cld':
sqn <- matrix(NA, nseq, seql)
cld <- integer(nseq)

## Generate the root sequence (with hidden gaps) and initiate its number of
## childs:
sample(
  x = c("-","A","C","G","T"),
  size = seql,
  prob = c(0.3,rep(0.175,4L)),
  replace = TRUE
) -> sqn[1L,]
cld[1L] <- maxCld

## Create the patristic distance storage vector and the edge list:
dst <- rep(NA,nseq*(nseq - 1L)/2)
edge <- list(from=c(), to=c())

## Evolving the sequences:
## i=2L # i=3L # i=4L # i=5L # i=11L
for(i in 2L:nseq) {
  ## 'i' is the new child node to be created
  
  ## 'wh' is the indices of the nodes that are allowed to be parent.
  wh <- which(as.logical(cld))
  
  ## This is the patristic distance of the potential parents from the root.
  ## Parent nodes that exist at the same time are allowed to hybridize:
  ## dst[dst_idx(nseq, 1L, wh)]
  
  ## For the time being, no hybridization is simulated.
  
  ## Selects one parent 'k' among the nodes that are allowed to reproduce:
  k <- if(length(wh)==1L) wh else sample(wh, 1L)
  ## If there is only a single parent available, select this one.
  
  ## Generate the new sequence by evolving parent 'k' using the nucleotide
  ## evolution model defined previously:
  for(j in 1L:seql)
    sqn[i,j] <- em[[j]]$evolve(sqn[k,j])
  
  ## Update the patristic distance matrix:
  
  ## 1. Insert the time in the locations:
  dst[dst_idx(nseq, i, k)] <- tm
  
  ## 2. Add the time to all other distances for active nodes:
  if(i>2L) {
    idx_k_others <- dst_idx(nseq, k, (1L:i)[-c(i,k)])
    idx_i_others <- dst_idx(nseq, i, (1L:i)[-c(i,k)])
    dst[idx_k_others] + tm -> dst[idx_i_others]
  }
  
  ## Initialize the parental status of the newly created node:
  cld[i] <- maxCld
  
  ## Decrease the parental status of the parent node:
  cld[k] <- cld[k] - 1L
  
  ## Insert the edge representing the parent-child link:
  edge$from <- c(edge$from, k)
  edge$to <- c(edge$to, i)
  
  cat(sprintf("%04d -> %04d (%d)\n",k,i,dst[dst_idx(nseq, 1L, i)]))
}



## Add the attributes to the distance matrix to make it loo as though it comes
## from function 'dist()'
list(
  Size = nseq,
  Diag = FALSE,
  Upper = FALSE,
  method = "patristic",
  call = match.call(),
  class = "dist"
) -> attributes(dst)

## Cleanup
rm(i,j,k,wh,idx_k_others,idx_i_others)

## Purge gap-only sequences:
sqn %<>% .[,apply(., 2L, function(x) !all(x == "-"))]

## as.data.frame(edge)

## concatenate(head(sqn), gap = TRUE, width = 70L)
## concatenate(head(sqn), gap = FALSE, width = 70L)

sqn %>%
  concatenate(width=70L) -> out

names(out) <- sprintf("sequence%03d", 1L:length(out))

write.fasta(out, "Simulated DNA sequence without the gaps before alignment.fst")

sqn %>%
  concatenate(gap=TRUE, width=70L) -> out2

names(out2) <- sprintf("sequence%03d", 1L:length(out2))

write.fasta(out2, "Simulated DNA sequence with the gaps.fst")

rm(I,em,nseq,seql,tm,maxCld,rate,cld,out,out2)

## max(dst)
nrow(sqn) %>%
  {sapply(1L:., function(x) max(dst[dst_idx(., x)]))} %>%
  which.max -> i

(1L:nrow(sqn))[-i][which.max(dst[dst_idx(nrow(sqn), i)])] -> j

dst[dst_idx(nrow(sqn), i, j)]

sqn[c(i,j),] %>%
  concatenate(gap=TRUE, width=70L) -> out3

names(out3) <- sprintf("sequence%03d", 1L:length(out3))

write.fasta(out3, "Simulated DNA sequence - most different.fst")

rm(out3)

## max(dst[dst_idx(nrow(sqn), 1L)])

edge %>%
  {unique(.$to)[!(unique(.$to) %in% unique(.$from))]} %>%
  {.[order(.)]} -> tips

(1L:nrow(sqn))[-tips] -> nodes







max(dst[dst_idx(nseq, 1L)])

## For convergence, pairs of nodes at the same (or a similar) distance from the
## root can be selected with a probability. That probability can be
## conditioned by the distance between the nodes

(1L:nseq)[-1L][which(dst[dst_idx(nseq, 1L)] == 10)] -> sel
dst[dst_idx(nseq, sel[2L], sel[3L])]




## and merged, with the distances
## between the mergers and any other initialized nodes 




sapply(sqn[1L,], function(x,y) as.numeric(x==y), y = sqn[1L,]) -> tmp
## dim(tmp)
image(tmp, col=grey(c(1,0)), asp=1)


tmp[1L:10L,1L:10L]


## Rendu ici...




base64enc::base64encode(serialize(I,NULL)) -> tmp
unserialize(base64enc::base64decode(tmp))
rm(tmp)

