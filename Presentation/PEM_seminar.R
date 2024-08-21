
### Animations

library(MPSEM)
library(magrittr)

## rm(list=ls())
source("PEM_seminar-aux.R")

recalculateExamples <- FALSE
## recalculateExamples <- TRUE
source("PEM_seminar-examples.R")

### Helps explain the calculation of the influence matrix

## The linear case:
plot1_graph(linear, "\"")           ## Display the graph
plot1_vertex(linear, 1L)            ## Pinpoint a vertex
plot1_edge(linear, 1L)              ## Pinpoint an edge

plot1_influence(linear, 3L, 4L)     ## Displays the partial influence matrix
plot1_influence(linear, 4L)         ## Displays the partial influence matrix

plot1_graph(linear)                 ## Displays the complete influence matrix


## The tree case:
plot1_graph(tree, "\"")             ## Display the graph
plot1_vertex(tree, 1L)              ## Pinpoint a vertex
plot1_edge(tree, 6L)                ## Pinpoint an edge

plot1_influence(tree, 4L, 5L)       ## Displays the partial influence matrix
plot1_influence(tree, 2L)           ## Displays the partial influence matrix

plot1_graph(tree)                   ## Displays the complete influence matrix

tree$PEM$u



## The reticulated case:

plot1_graph(reticulate, "\"")       ## Display the graph
plot1_vertex(reticulate, 3L)        ## Pinpoint a vertex
plot1_edge(reticulate, 6L)          ## Pinpoint an edge

plot1_influence(reticulate, 6L, 7L) ## Displays the partial influence matrix
plot1_influence(reticulate, 7L)     ## Displays the partial influence matrix


plot1_graph(reticulate)             ## Displays the complete influence matrix

reticulate$PEM$u


### 

PEM.build(Linear, d="len", sp = "")




