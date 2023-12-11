
##

## rm(list=ls())

library(magrittr)

c(NA,1,1,1,1, 1,NA,1,2,3, 1,1,NA,1,2, 1,2,1,NA,1, 1,3,2,1,NA) %>%
  matrix(
    5L, 5L,
    dimnames=list(
      c("-","A","C","G","T"),
      c("-","A","C","G","T")
    )
  ) -> I
diag(I) <- -rowSums(I, na.rm = TRUE)
## rowSums(I)

genem <- function(I, t, rho) {
  expMat <- function(x) {
    eig <- eigen(x)
    y <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
    dimnames(y) <- dimnames(x)
    y
  }
  Mt <- expMat(t*rho*I)
  list(
    recalculate = function(t, rho)
      Mt <<- expMat(t*rho*I),
    evolve = function(N)
      sample(x = c("-","A","C","G","T"), size = 1L, prob=Mt[N,]),
    getMt = function() Mt
  )
}

genem(I, 1, 0.001) -> em1

## gen$evolve("Z")

gen$getMt()

gen$evolve("-")
gen$evolve("A")
gen$evolve("C")
gen$evolve("G")
gen$evolve("T")

gen$recalculate(1,0.01)

gen$getMt()

gen$evolve("-")
gen$evolve("A")
gen$evolve("C")
gen$evolve("G")
gen$evolve("T")

rgamma(
  1500L,
  shape = 5,
  scale = 0.0001
) -> rate
## plot(density(rate))

em <- list()

for(i in 1L:1500L)
  em[[i]] <- genem(I, 1, rate[i])

sample(
  x = c("-","A","C","G","T"),
  size = 1500L,
  prob = rep(0.2,5L),
  replace = TRUE
) -> initial

final <- character(1500L)

## i=1L
for(i in 1L:1500L)
  final[i] <- em[[i]]$evolve(initial[i])

## paste(initial[initial!="-"], collapse="")
## paste(final[final!="-"], collapse="")

data.frame(initial, final, compare = initial!=final) -> cmp
cmp[cmp$compare,-3L]


