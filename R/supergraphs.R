#' @title Supergraph with given boxicity
#' @description Create a supergraph with given boxicity using simulated annealing (SA)
#'
#' @param g igraph object
#' @param dim integer. target boxicity
#' @param perm starting permutation for SA. If NULL, a random permutation is created
#' @param iter integer. number of iterations for SA
#' @param temp integer. starting temperature for SA
#' @param tmax integer. number of function evaluations at each temperature for SA
#' @param verbose logical. print report during SA (defaults to FALSE)
#' @return a list with entries
#' \item{perm}{permutation vector. All permutations are concatenated to one long vector}
#' \item{ged}{graph edit distance from original graph}
#' \item{A}{adjacency matrix of supergraph with given boxicity}
#' @importFrom igraph vcount as_adj delete.vertices neighborhood
#' @importFrom stats optim
#' @references
#' Chandran, L. S., Francis, M. C. & Sivadasan, N. Geometric representation of graphs in low dimension using axis parallel boxes. Algorithmica 56, 129.
#' @export

superbox_graph <- function(l,dim=1,perm=NULL,iter=15000,temp=10,tmax=5,verbose=FALSE){
  B <- as_adj(l,sparse = FALSE)
  MSE <- struct_equi(l)
  l1 <- delete.vertices(l,which(duplicated(MSE)))
  adj <- neighborhood(l1)
  adj <- lapply(adj,function(x) as.integer(x))
  A <- as_adj(l1,sparse=FALSE)
  if(is.null(perm)){
    perm <- c()
    for(i in 1:dim){
      perm <- c(perm,sample(1:vcount(l1)))
    }
  } else{
    perm <- perm[!duplicated(MSE)]
    if(length(perm)!=vcount(l1)/dim){
      stop("perm has the wrong length")
    }
  }
  res_sa <- optim(par = perm, fn = gedn, A = A,adj = adj,dim = dim,gr = genpermn,method = "SANN",
                  control = list(maxit = iter, temp = temp, tmax = tmax, trace = verbose))

  perm <- res_sa$par
  Ared <- perm2adj(adj,perm,dim)

  Afull <- Ared[MSE,MSE]
  for(grp in unique(MSE)){
    Afull[MSE==grp,MSE==grp] <- B[MSE==grp,MSE==grp]
  }

  # res <- list(perm=perm,ged=res_sa$value,A_MSE = Ared,A=Afull,MSE = MSE)
  n <- vcount(l1)
  permMSE <- c()
  for(i in 1:dim){
    perm_cur <- perm[((i-1)*n+1):(i*n)][MSE]
    permMSE <- c(permMSE,perm_cur)
  }
  res <- list(perm=permMSE,ged=res_sa$value,A=Afull)
  res
}

#' @title Box representation from permutations
#' @description Create a box representation from permutations
#' @param g igraph object.
#' @param perm integer vector of length n times dim
#' @param dim integer. dimensionality of boxes
#' @importFrom igraph vcount as_adj delete.vertices neighborhood
#' @return coordinates
#' @references
#' Chandran, L. S., Francis, M. C. & Sivadasan, N. Geometric representation of graphs in low dimension using axis parallel boxes. Algorithmica 56, 129.
#' @export

perm2box <- function(g,perm,dim){
  adj <- neighborhood(g)
  adj <- lapply(adj,function(x) as.integer(x))
  n <- length(adj)
  if(length(perm)!=n*dim){
    stop("perm must have the same length as dim*nodes")
  }
  coords <- vector("list",dim)
  for(i in 1:dim){
    perm_cur <- perm[((i-1)*n+1):(i*n)]
    Aperm <- lapply(adj,function(x) perm_cur[x])
    coords[[i]] <- cbind(sapply(Aperm,function(x) c(min(x))),perm_cur)

  }
  coords
}

# helper ----
perm2adj <- function(N,perm,dim){
  n <- length(N)
  As <- vector("list",dim)
  for(i in 1:dim){
    perm_cur <- perm[((i-1)*n+1):(i*n)]
    Nperm <- lapply(N,function(x) perm_cur[x])
    # xy1 <- cbind(sapply(Nperm,function(x) c(min(x))),perm_cur)
    # a <- xy1[,1]
    # b <- xy1[,2]
    a <- sapply(Nperm,function(x) c(min(x)))
    b <- perm_cur
    As[[i]] <- getA_cpp(a,b)
    diag(As[[i]]) <- 0
  }
  B <- As[[1]]==1
  if(dim>1){
    for(i in 2:dim){
      B <- B & (As[[i]]==As[[i-1]])
    }
  }
  B <- B + 0
  B
}

gedn <- function(perm,A,adj,dim){
  Anew <- perm2adj(adj,perm,dim)
  sum(Anew-A)/2
}

genpermn <- function(perm,A,adj,dim) {
  n <- length(adj)
  d <- length(perm)/n
  for(i in 1:d){
    perm_cur <- perm[((i-1)*n+1):(i*n)]
    changepoints <- sample(1:n,2)
    tmp <- perm_cur[changepoints[1]]
    perm_cur[changepoints[1]] <- perm_cur[changepoints[2]]
    perm_cur[changepoints[2]] <- tmp
    perm[((i-1)*n+1):(i*n)] <- perm_cur
  }
  perm
}


struct_equi <- function(g) {
  # adj <- lapply(igraph::get.adjlist(g), function(x) x - 1)
  adj <- lapply(igraph::neighborhood(g,mindist = 1), function(x) x - 1)
  deg <- igraph::degree(g)
  P <- mse(adj, deg)
  MSE <- which((P + t(P)) == 2, arr.ind = T)
  if (length(MSE) >= 1) {
    MSE <- t(apply(MSE, 1, sort))
    MSE <- MSE[!duplicated(MSE), ]
    g <- igraph::graph.empty()
    g <- igraph::add.vertices(g, nrow(P))
    g <- igraph::add.edges(g, c(t(MSE)))
    g <- igraph::as.undirected(g)
    MSE <- igraph::components(g,"weak")$membership
  } else {
    MSE <- 1:nrow(P)
  }
  return(MSE)
}
