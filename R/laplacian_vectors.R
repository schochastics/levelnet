#' @title Laplacian Vectors
#' @description  Returns Laplacian eigenvectors associated with the k smallest positive eigenvalues
#'
#' @param g igraph object
#' @param k number of vectors to return
#' @return data.frame of vectors
#' @author David Schoch
#' @export

laplacian_vectors <- function(g,k=2){
  if(igraph::components(g)$no>1){
    stop("g is not connected")
  }
  sL <- eigen(igraph::laplacian_matrix(g))
  n <- igraph::vcount(g)
  as.data.frame(sL$vectors[,(n-k):(n-1)])
}
