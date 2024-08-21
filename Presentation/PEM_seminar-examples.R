
### Code generating the examples

## Linear
if(recalculateExamples) {
  
  list(
    N = 7L,
    len = c(1,1,1,1,1,1)
  ) -> linear
  
  linear %>%
    {pop.graph(
      n = linear$N,
      vertex = list(
        x = c(0,cumsum(.$len)/sum(.$len)),
        y = rep(0.5, linear$N)
      ),
      label = LETTERS[1L:linear$N]
    )} -> linear$graph
  
  linear$graph %<>%
    add.edge(
      from = 1L:(linear$N - 1L),
      to = 2L:linear$N,
      edge = list(len=linear$len),
      label = sprintf("%02d",1L:(linear$N - 1L))
    )
  
  linear$Infl <- PEMInfluence(linear$graph)
  
  save(linear, file="linear.rda")
  
} else
  load(file="linear.rda")

## Tree
if(recalculateExamples) {
  
  list(
    N = 9L
  ) -> tree
  
  tree %>%
    {pop.graph(
      n = tree$N,
      vertex = list(
        s = c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
        x = c(0.00,0.50,0.70,0.60,1.00,1.00,1.00,1.00,1.00),
        y = c(0.50,0.75,0.87,0.12,1.00,0.75,0.50,0.25,0.00)
      ),
      label = LETTERS[1L:tree$N]
    )} -> tree$graph
  
  tree$graph %<>%
    add.edge(
      from = c(1,1,2,3,3,2,4,4),
      to   = c(2,4,3,5,6,7,8,9),
      edge = list(len=c(0.5,0.6,0.2,0.3,0.3,0.5,0.4,0.4)),
      label = sprintf("%02d",1L:8L)
    )
  
  tree$Infl <- PEMInfluence(tree$graph)
  
  tree$PEM <- PEM.build(tree$graph, d="len", sp="s")
  
  save(tree, file="tree.rda")
  
} else
  load(file="tree.rda")

## Reticulate
if(recalculateExamples) {
  
  list(
    N = 10L
  ) -> reticulate
  
  reticulate %>%
    {pop.graph(
      n = reticulate$N,
      vertex = list(
        s = c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
        x = c(0.00,0.20,0.30,0.50,0.70,1.00,1.00,1.00,1.00,1.00),
        y = c(0.50,0.70,0.20,0.45,0.87,1.00,0.75,0.50,0.25,0.00)
      ),
      label = LETTERS[1L:reticulate$N]
    )} -> reticulate$graph
  
  reticulate$graph %<>%
    add.edge(
      from = c(1,1,2,3,2,5,5,4,4,3),
      to   = c(2,3,4,4,5,6,7,8,9,10),
      edge = list(len=c(0.2,0.3,0.3,0.2,0.4,0.4,0.4,0.5,0.5,0.7)),
      label = sprintf("%02d",1L:10L)
    )
  
  reticulate$Infl <- PEMInfluence(reticulate$graph)
  
  reticulate$PEM <- PEM.build(reticulate$graph, d="len", sp="s")
  
  save(reticulate, file="reticulate.rda")
  
} else
  load(file="reticulate.rda")
