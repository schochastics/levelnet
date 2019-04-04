#' @title hypergeometric filter
#' @description Extract significant edges with a hypergeometric distribution.
#' @import igraph
#' @param g igraph object. The two-mode network
#' @param type which mode to project on
#' @param alpha significants level
#' @return backbone of weighted network
#' @author David Schoch
#' @export
#'

hypergeometric_filter <- function(g,type="true",alpha=0.05){
  cat("calculating projection\n")
  gb <- bipartite_projection(g,which = type)
  cat("getting matrix\n")
  P <- get.adjacency(gb,"both",attr="weight",sparse = FALSE,)

  degs <- as.integer(degree(g)[V(g)$type==as.logical(type)])
  m <- vcount(g) - vcount(gb)

  cat("calculate hypergeometric\n")
  H <- hypergeom(P,degs,m)
  hg <- graph_from_adjacency_matrix(H<=alpha,"undirected",diag = FALSE)
  hg
}
