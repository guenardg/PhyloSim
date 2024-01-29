
## 

dMatChk <- function(object) {
  
  errors <- character()
  n <- object@size
  l <- length(object@re)
  
  if(l != n*(n - 1L)/2)
    c(
      errors,
      sprintf("Object size does not match the number of data!")
    ) -> errors
  
  if(length(object@im))
    if(l != length(object@im))
      c(
        errors,
        sprintf("Number of real and imaginary values do not match!")
      ) -> errors
  
  if(length(object@names))
    if(n != length(object@names))
      c(
        errors,
        sprintf("Length of 'names' does not equal size!")
      ) -> errors
  
  if(length(errors) == 0) TRUE else errors
}

setClass(
  Class = "dMat",
  representation(
    re = "numeric",
    im = "numeric",
    size = "integer",
    names = "character",
    method = "character"),
  validity = dMatChk
)

dMat <- function(data, size, names = character(), method = character()) {
  
  if(missing(data)) {
    if(missing(size)) {
      data <- numeric()
      size <- 0L
    } else
      data <- rep(NA_real_, size*(size - 1L)/2L)
  } else {
    if(is.matrix(data)) {
      if(nrow(data) != ncol(data))
        stop("Matrix given as data is not square!")
      size <- nrow(data)
      data <- data[lower.tri(data)]
    } else {
      if(missing(size))
        size <- as.integer(0.5 + 0.5*sqrt(1L + 8L*length(data)))
    }
  }
  
  if(is.complex(data)) {
    new("dMat", re=Re(data), im=Im(data), size=size, names=names, method=method)
  } else {
    new("dMat", re=data, im=numeric(), size=size, names=names, method=method)
  }
}

dMat()
dMat(method = "patristic")
dMat(size = 8L, method = "patristic")
dMat(1:10, method = "patristic")
dMat(matrix(1:9,3L,3L), method = "patristic")
dMat(matrix(1:9,3L,3L), names = LETTERS[1L:3L], method = "patristic")
dMat(matrix(1:25,5L,5L), names = LETTERS[1L:5L], method = "patristic")
dMat(matrix(1:49,7L,7L), names = LETTERS[1L:7L], method = "patristic")
dMat(matrix(complex(real = 1:49, imaginary = 4:52),7L,7L),
     names = LETTERS[1L:7L], method = "patristic") -> tmp
tmp@im
rm(tmp)

dMat(
  matrix(1:100,10L,10L),
  names = LETTERS[1L:10L],
  method = "patristic"
) -> object








setMethod("[",signature(object="dMat"),
          function(object, i, j, drop=TRUE) {
            
          })

setMethod("show",signature(object="dMat"),
          function(object) {
            cat("-----------\nA dMat-class object\n")
            cat("-----------\n")
            cat("Distance method:",
                if(length(object@method)) object@method else "Unknown",
                "\n")
            cat("Number of objects:", object@size, "\n")
            cat("Complex-valued:", !!length(object@im), "\n")
            if(length(object@names)) {
              if(length(object@names) <= 8L) {
                nms <- paste(object@names,collapse=", ")
              } else
                nms <- paste(c(head(object@names, 7L),"..."),collapse=", ")
              cat("Names:",nms,"\n")
            }
            if(length(object@re)) {
              mx <- if(object@size <= 8L) object@size else 7L
              lab <- if(length(object@names)) object@names[1L:mx] else as.character(1L:mx)
              mm <- matrix("", mx - 1L, mx - 1L, dimnames = list(lab[-1L],lab[-mx]))
              for(j in 2L:mx) {
                .Call(
                  "PhyloSim_dstIdx",
                  PACKAGE="PhyloSim",
                  object@size, 1L:(j - 1L), j) -> idx
                if(length(object@im)) {
                  mm[j - 1L,1L:(j - 1L)] <- complex(
                    real = object@re[idx],
                    imaginary = object@im[idx]
                  )
                } else
                  mm[j - 1L,1L:(j - 1L)] <- object@re[idx]
              }
              mm[is.na(mm)] <- "NA"
              print(as.data.frame(mm))
              if(object@size > 8L) cat("...\n")
            }
            invisible(NULL)
          })





## Read mutation intensity matrix:
read.mutationMat <- function(file) {
  
  yaml.load_file(file) -> tmp
  
  matrix(
    NA, 5L, 5L,
    dimnames=list(
      c("-","A","C","G","T"),
      c("-","A","C","G","T")
    )
  ) -> I
  
  I["-",names(tmp$insertion)] <- unlist(tmp$insertion)
  I[names(tmp$deletion),"-"] <- unlist(tmp$deletion)
  
  unlist(tmp$transition)[c("AG","GA","CT","TC")] ->
    I[rbind(c("A","G"),c("G","A"),c("C","T"),c("T","C"))]
  
  for(i in names(tmp$transversion))
    I[i,names(tmp$transversion[[i]])] <- unlist(tmp$transversion[[i]])
  
  diag(I) <- -rowSums(I, na.rm = TRUE)
  
  I
}

molEvolMod <- function(I, t, rho) {
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

if(FALSE) {
  dst_idx_R <- function(n, i, j) {
    if(missing(j))
      j <- (1L:n)[-i]
    ni <- length(i)
    nj <- length(j)
    nn <- max(ni,nj)
    idx <- cbind(rep(i,length.out=nn),rep(j,length.out = nn))
    sel <- which(idx[,1L]>idx[,2L])
    if(length(sel)) {
      idx[sel,1L] -> ii
      idx[sel,2L] -> idx[sel,1L]
      idx[sel,2L] <- ii
    }
    sel <- which(idx[,1L]!=idx[,2L])
    y <- rep(NA,nn)
    y[sel] <- idx[sel,2L] + (idx[sel,1L] - 1L)*n - idx[sel,1L]*(idx[sel,1L] + 1)/2
    y
  }
  
  ## dst_idx_R(10, 1L)
  ## dst_idx_R(10, 2L)
  ## dst_idx_R(10, 5L)
  ## dst_idx_R(10, 5L, c(1L,3L:7L,9L))
}

## dst_idx(10, 1L)
## dst_idx(10, 2L)
## dst_idx(10, 5L)
## dst_idx(10, 5L, c(1L,3L:7L,9L))

concatenate <- function(x, gap = FALSE, width = NULL) {
  apply(
    x,
    1L,
    function(x, gap, width) {
      if(!gap) x <- x[x!="-"]
      x <- paste(x, collapse="")
      if(!is.null(width)) {
        ss <- seq(1L,nchar(x),width)
        paste(
          paste(substring(x,ss[-length(ss)],ss[-1L] - 1L),collapse = "\n"),
          substring(x,ss[length(ss)],nchar(x)),sep="\n"
        ) -> x
      }
      x
    },
    gap = gap,
    width = width
  )
}

write.fasta <- function(x, file, ...) {
  con <- file(file, open="w", blocking=TRUE, ...)
  
  if(is.null(names(x)))
    names(x) <- sprintf("SEQ%d",1L:length(x))
  
  for(i in 1L:length(x)) {
    cat(sprintf(">%s\n", names(x)[i]), file=con, append=TRUE)
    cat(sprintf("%s\n", x[i]), file=con, append=TRUE)
  }
  
  close(con)
  
  invisible(NULL)
}
