#' @title hypergeometric filter
#' @description Extract significant edges with a hypergeometric distribution.
#' @import igraph
#' @param g igraph object. The two-mode network
#' @param proj which mode to project on
#' @param alpha significants level
#' @return backbone of weighted network
#' @author David Schoch
#' @export
#'

hypergeometric_filter <- function(g,proj="true",alpha=0.05){

  gb <- bipartite_projection(g,proj = type)
  P <- get.adjacency(gb,"both",attr="weight",sparse = FALSE,)
  degs <- as.integer(degree(g)[V(g)$type==as.logical(type)])
  m <- vcount(g) - vcount(gb)
  H <- hypergeom(P,degs,m)
  hg <- graph_from_adjacency_matrix(H<=alpha,"undirected",diag = FALSE)
  hg
}
