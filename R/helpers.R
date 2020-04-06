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
  M <- matrix(0,length(mcl),igraph::vcount(g))
  for(i in 1:length(mcl)){
    M[i,mcl[[i]]] <- 1
  }
  M
}

#' @rdname helpers
#' @export
is_bipartite1 <- function(g){
  adj <- igraph::as_adj(g,"both",sparse=F)
  comps <- igraph::components(g,"weak")
  if(comps$no==1){
    return(isBipartite(adj,igraph::vcount(g),1))
  } else{
    bip_bool <- rep(F,comps$no)
    for(i in 1:comps$no){
      g1 <- igraph::induced_subgraph(g,which(comps$membership==i))
      adj <- igraph::as_adj(g1,"both",sparse=F)
      bip_bool[i] <- isBipartite(adj,igraph::vcount(g1),1)
    }
    return(all(bip_bool))
  }
}

#' @useDynLib levelnet
#' @importFrom Rcpp sourceCpp
NULL


