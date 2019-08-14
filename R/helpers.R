#' @title helper function
#' @description small helper functions
#'
#' @param g igraph object.
#' @name helpers
#' @return igraph object
#' @author David Schoch

NULL

#' @rdname helpers
#' @export
clique_vertex_mat <- function(g){
  if(!igraph::is.igraph(g)){
    stop("g must be an igraph object")
  }
  if(igraph::is.directed(g)){
    warning("g is directed. Underlying undirected graph is used")
    g <- igraph::as.undirected(g)
  }
  mcl <- igraph::max_cliques(g)
  # mcl <- lapply(mcl,sort)
  # mcl <- mcl[order(unlist(lapply(mcl,sum)))]
  M <- matrix(0,length(mcl),igraph::vcount(g))
  for(i in 1:length(mcl)){
    M[i,mcl[[i]]] <- 1
  }
  M
}

#' @rdname helpers
#' @export
is_bipartite1 <- function(g){
  adj <- igraph::get.adjacency(g,"both",sparse=F)
  comps <- igraph::components(g,"weak")
  if(comps$no==1){
    return(isBipartite(adj,igraph::vcount(g),1))
  } else{
    bip_bool <- rep(F,comps$no)
    for(i in 1:comps$no){
      g1 <- igraph::induced_subgraph(g,which(comps$membership==i))
      adj <- igraph::get.adjacency(g1,"both",sparse=F)
      bip_bool[i] <- isBipartite(adj,igraph::vcount(g1),1)
    }
    return(all(bip_bool))
  }
}

robinson_min <- function(A){
  n <- nrow(A)
  Anew <- matrix(NA,n,n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(i!=j){
        a <- A[i,j]
        b <- Anew[i,j-1]
        c <- Anew[i+1,j]
        Anew[i,j] <- Anew[j,i] <- min(c(a,b,c),na.rm = T)
      }
    }
  }
  diag(Anew) <- diag(A)
  Anew
}


robinson_max <- function(A){
  n <- nrow(A)
  Anew <- matrix(NA,n,n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(i!=j){
        a <- A[i,j]
        b <- ifelse(j==n,NA,Anew[i,j+1])
        c <- Anew[i-1,j]
        Anew[i,j] <- Anew[j,i] <- min(c(a,b,c),na.rm = T)
      }
    }
  }
  diag(Anew) <- diag(A)
  Anew
}

#' @useDynLib levelnet
#' @importFrom Rcpp sourceCpp
NULL


