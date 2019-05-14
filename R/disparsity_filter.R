#' @title Disparity Filter
#' @description Extract significant edges with disparsity filter.
#'
#' @param g igraph object. either two-mode or weighted network
#' @param proj string. Which mode to project on ("true"/"false")
#' @param alpha significants level
#' @param cut_mode 'and' or 'or'
#' @return backbone of weighted network
#' @author David Schoch
#' @references Serrano et al. (2009). Extracting the multiscale backbone of complex weighted networks
#' @export
#'

disparsity_filter <- function(g,proj="true",alpha=0.05,cut_mode="or"){
  if(!any(igraph::vertex_attr_names(g)=="name")){
    igraph::V(g)$name <- 1:igraph::vcount(g)
  }
  if(igraph::is.bipartite(g)){
    g <- igraph::bipartite_projection(g,which=proj)
  }
  cut_mode <- match.arg(cut_mode,c("and","or"))
  A <- igraph::get.adjacency(g,sparse=F,attr="weight")
  n <- nrow(A)
  strength <- rowSums(A)
  B <- diag(1/strength) %*% A
  degs <- rowSums(A>0)

  # df <- as.data.frame(t(utils::combn(1:n,2)))
  df <- as.data.frame(expand.grid(1:n,1:n))
  df <- df[df[["Var1"]]!=df[["Var2"]],]
  names(df) <- c("V1","V2")
  df[["d"]] <- degs[df[["V1"]]]
  df[["pij"]] <- apply(df[,1:2],1,function(x) B[x[1],x[2]])
  df[["aij"]] <- apply(df[,3:4],1,function(x) tryCatch(1-(x[1]-1)*stats::integrate(alpha_func,lower=0,upper=x[2],k=x[1])$value,error=function(e) 1))
  df[["V1"]] <- rownames(A)[df[["V1"]]]
  df[["V2"]] <- rownames(A)[df[["V2"]]]


  gf <- igraph::graph_from_data_frame(df[which(df[["aij"]]<alpha),1:2],directed = T)
  if(cut_mode=="or"){
    d <- igraph::as.undirected(gf,mode="collapse")
  } else if(cut_mode == "and"){
    d <- igraph::as.undirected(gf,mode="mutual")
  }
  d
}


alpha_func <- function(x,k){
  (1-x)^(k-2)
}
