#' @title Multisweep Lex-BFS
#' @description  Multisweep lexicograpical BFS to recoginze interval graphs.
#'
#' @param g igraph object
#' @param k number of sweeps
#' @return permutation
#' @author David Schoch
#' @export
#'

multiLexBFS <- function(g,k=4){
  perm_init <- lexBFS(g)
  for(i in 1:k){
    perm <- lexBFS_plus(g,perm_init)
    perm_init <- perm
  }
  perm
}

lexBFS <- function(g){
  label <- vector("list",igraph::vcount(g))
  visited <- rep(F,igraph::vcount(g))
  perm <- 1:igraph::vcount(g)
  adj <- igraph::get.adjlist(g)
  for(i in 1:igraph::vcount(g)){
    S <- max_label(label,visited)
    p <- sample(S,1)
    perm[p] <- i
    visited[p] <- T
    Np <- adj[[p]]
    Np <- Np[!visited[Np]]
    for(x in Np){
      label[[x]] <- c(label[[x]],igraph::vcount(g)-i)
    }
  }
  order(perm)
}
lexBFS_plus <- function(g,perm_init){
  label <- vector("list",igraph::vcount(g))
  visited <- rep(F,igraph::vcount(g))
  perm <- 1:igraph::vcount(g)
  adj <- igraph::get.adjlist(g)
  for(i in 1:igraph::vcount(g)){
    S <- max_label(label,visited)
    p <- S[which.max(match(S,perm_init))]
    perm[p] <- i
    visited[p] <- T
    Np <- adj[[p]]
    Np <- Np[!visited[Np]]
    for(x in Np){
      label[[x]] <- c(label[[x]],igraph::vcount(g)-i)
    }
  }
  order(perm)
}

max_label <- function(label,visited){
  lex <- unlist(lapply(label,function(x) paste0("0",x,collapse="")))
  idx <- which(!visited)
  lex <- lex[idx]
  idmax <- which(lex==max(lex))
  idx[idmax]
}
