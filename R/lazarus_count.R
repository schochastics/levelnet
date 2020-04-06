#' @title Lazarus Count
#' @description  Calculates the Lazarus count of a matrx/graph.
#'
#' @param g either an igraph object or a (0,1)-matrix
#' @param perm permutation or NULL
#' @param mode one of "mcl" (clique vertex matrix), "cols" (Lazarus count of columns) "rows" (Lazarus count of rows) or "sym" (Lazarus count of both columns and rows).
#' @details The Lazarus count of a matrix is the number of "holes" in each column.
#' A hole is a number of zero entries surrounded by ones. For an interval graph, this count is zero for the [clique_vertex_mat].
#' If `perm` is NULL, a row permutation based on the Fiedler vector of the Laplacian is calculated.
#' @return Lazarus count of g
#' @author David Schoch
#' @examples
#' set.seed(10)
#' #the lazarus count of an interval graph is zero
#' g <- graph_interval(n = 10,r = 1)
#' lazarus_count(g, mode = "mcl")
#' @export
lazarus_count <- function(g,perm = NULL,mode = "cols"){
  #cols,rows,sym,mcl
  if(is.null(perm)){
    perm <- order(fiedler_order(g,mode))
  }
  if(mode=="mcl" & !igraph::is.igraph(g)){
    stop("g must be an igraph object if mode=mcl")
  }

  if(igraph::is.igraph(g)){
    A <- igraph::get.adjacency(g,sparse = F)
    if(mode=="mcl"){
      if(igraph::is.directed(g)){
        warning("graph is directed, using undirected version instead")
        g <- igraph::as.undirected(g)
      }
      mcl <- igraph::max_cliques(g)
      M <- matrix(0,length(mcl),igraph::vcount(g))
      for(i in 1:length(mcl)){
        M[i,mcl[[i]]] <- 1
      }
      A <- M[perm,]
    }
  } else{
    A <- g
  }
  if(mode=="sym"){
    diag(A) <- NA
  }
  if(mode=="cols"){
    A <- A[perm,]
  } else if(mode=="rows"){
    A <- t(A[,perm])
  } else if(mode=="sym"){
    A <- A[perm,perm]
  }

  id1 <- t(apply(A,2,function(x){
    idx <- which(x==1)
    if(length(idx)==0){
      idx <- 0
    }
    idx[c(1,length(idx))]
  }))
  lcount <- sum(sapply(1:ncol(A),function(x) sum(A[id1[x,1]:id1[x,2],x]==0,na.rm=T)),na.rm = T)
  lcount
}

#' @title Check whether graph is interval graph
#' @description  Check whether graph is interval graph.
#'
#' @param g igraph object
#' @return Logical scalar, whether graph is an interval graph
#' @details This function is not very efficient since it relies on the clique vertex matrix. More efficient linear time algorithms will be implemented in the future.
#' @author David Schoch
#' @export
is_interval <- function(g){
  perm <- order(fiedler_order(g,mode="mcl"))
  lz <- lazarus_count(g,perm,mode="mcl")
  lz==0
}

#' @title Permutation induced by Fiedler vector
#' @description  Returns the permutation induced by sorting the Fiedler vector of the Laplacian matrix of a graph
#'
#' @param g an igraph object or a (0,1)-matrix
#' @param mode one of "mcl" (clique vertex matrix), "cols" (Lazarus count of columns) "rows" (Lazarus count of rows) or "sym" (Lazarus count of both columns and rows).
#' @return numeric vector
#' @author David Schoch
#' @export
fiedler_order <- function(g,mode="cols"){
  #cols, rows, mcl
  if(mode=="mcl" & !igraph::is.igraph(g)){
    stop("g must be an igraph object if mode=mcl")
  }

  if(igraph::is.igraph(g)){
    A <- igraph::get.adjacency(g,sparse = F)
  } else{
    A <- g
  }

  if(mode=="mcl"){
    if(igraph::is.directed(g)){
      warning("graph is directed, using undirected version instead")
      g <- igraph::as.undirected(g)
    }
    mcl <- igraph::max_cliques(g)
    M <- matrix(0,length(mcl),igraph::vcount(g))
    for(i in 1:length(mcl)){
      M[i,mcl[[i]]] <- 1
    }
    A <- M%*%t(M)
  }
  if(mode=="cols"){
    A <- t(A)%*%A
  } else if(mode=="rows"){
    A <- A%*%t(A)
  }
  L <- diag(rowSums(A))-A
  sL <- eigen(L)
  fiedler <- which(round(sL$values,8)==0)[1]-1
  fv <- sL$vectors[,fiedler]

  fv
}
