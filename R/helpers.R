#' @title helper function
#' @description small functions to deal with typical network problems
#'
#' @param g igraph object.
#' @param eattr edge attribute that contains 1/-1
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
count_triangle_type <- function(g,eattr="type"){
  if(!eattr%in%igraph::edge_attr_names(g)){
    stop("edge attribute not found")
  }
  eattrV <- get.edge.attribute(g,eattr)
  if(!all(eattrV%in%c(-1,1))){
    stop("edge attribute may only contain -1 and 1")
  }
  tmat <- t(matrix(igraph::triangles(g),nrow=3))
  # igraph::E(g)$typeI <- ifelse(igraph::E(g)$type=="pos",1,-1)

  emat <- t(apply(tmat,1,function(x) c(igraph::get.edge.ids(g,x[1:2]),
                                       igraph::get.edge.ids(g,x[2:3]),
                                       igraph::get.edge.ids(g,x[c(3,1)]))))


  emat[,1] <- eattrV[emat[,1]]
  emat[,2] <- eattrV[emat[,2]]
  emat[,3] <- eattrV[emat[,3]]
  emat <- t(apply(emat,1,sort))
  emat_df <- as.data.frame(emat)
  res <- by(emat_df,list(emat_df[["V1"]],emat_df[["V2"]],emat_df[["V3"]]),
            function(x) c(E1=mean(x$V1),E2=mean(x$V2),E3=mean(x$V3),count=nrow(x)))
  do.call(rbind,res)
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


