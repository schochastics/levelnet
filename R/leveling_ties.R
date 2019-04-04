#' @title Extract leveling ties
#' @description Extracting leveling ties.
#'
#' @param g igraph object. The two-mode network
#' @param proj string indicating mode to project on
#' @param thresh weight threshold
#' @param alpha se addition
#' @return network of leveling ties
#' @author David Schoch
#' @export
#'

leveling_ties <- function(g,proj="true",thresh=0,alpha=2){
  gs <- igraph::bipartite.projection(g,which=proj)

  proj_logi <- as.logical(proj)

  df <- igraph::as_data_frame(gs,"edges")
  deg <- igraph::degree(g)[igraph::V(g)$type==proj_logi]
  tot <- sum(igraph::V(g)$type!=proj_logi)

  df[["nA"]] <- deg[match(df[["from"]],names(deg))]
  df[["nB"]] <- deg[match(df[["to"]],names(deg))]
  df[["expn"]] <- df[["nA"]]/tot * df[["nB"]]/tot * tot
  df[["se"]] <- sqrt(df[["weight"]]/tot * (1-df[["weight"]]/tot) / tot)
  df[["bound"]] <- df[["expn"]] + alpha * df[["se"]]

  idx <- (df[["bound"]] >= df[["weight"]] & df[["weight"]] > thresh)

  igraph::graph_from_edgelist(as.matrix(df[idx,1:2]),directed = FALSE)
}
