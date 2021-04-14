#' @title Supergraph with low boxicity
#' @description Create supergraph with low boxicity
#'
#' @param g igraph object
#' @param dim integer. target boxicity (1-5)
#' @param iter number of iterations for simulated annealing step
#' @importFrom igraph vcount as_adj delete.vertices neighborhood
#' @importFrom stats optim
#' @export

superbox_graph <- function(g,dim,iter=50000){
  if(!igraph::is.igraph(g)){
    stop("g must be an igraph object")
  }
  if(igraph::is.directed(g)){
    stop("g must be undirected")
  }
  if(dim==1){
    return(create_superbox1(g,iter))
  } else if(dim==2){
    return(create_superbox2(g,iter))
  } else if(dim==3){
    return(create_superbox3(g,iter))
  } else if(dim==4){
    return(create_superbox4(g,iter))
  } else if(dim==5){
    return(create_superbox5(g,iter))
  } else {
    stop("dim must be less equal 5")
  }

}



create_superbox1 <- function(l,iter=500000){
  B <- as_adj(l,sparse = FALSE)
  MSE <- struct_equi(l)
  l1 <- delete.vertices(l,which(duplicated(MSE)))
  adj <- neighborhood(l1)
  adj <- lapply(adj,function(x) as.integer(x))
  A <- as_adj(l1,sparse=FALSE)
  perm <- sample(1:vcount(l1))
  res_sa <- optim(par = perm, fn = ged1, A = A,adj = adj,gr = genperm_sb1,method = "SANN",
                  control = list(maxit = iter, temp = 100, tmax = iter/1000, trace = FALSE,
                                 REPORT = 5))

  res_sa <- optim(par = res_sa$par, fn = ged1, A = A,adj = adj,gr = genperm_sb1,method = "SANN",
                  control = list(maxit = iter/10, temp = 10, tmax = iter/1000, trace = FALSE,
                                 REPORT = 5))

  Ared <- perm2int(adj,res_sa$par)

  Afull <- Ared[MSE,MSE]
  for(grp in unique(MSE)){
    Afull[MSE==grp,MSE==grp] <- B[MSE==grp,MSE==grp]
  }

  res <- list(perm1 = res_sa$par,ged=res_sa$value,A_MSE = Ared,A=Afull,MSE = MSE)
  res
}

create_superbox2 <- function(l,iter=500000){
  B <- as_adj(l,sparse = FALSE)
  MSE <- struct_equi(l)
  l1 <- delete.vertices(l,which(duplicated(MSE)))
  adj <- neighborhood(l1)
  adj <- lapply(adj,function(x) as.integer(x))
  A <- as_adj(l1,sparse=FALSE)
  perm <- c(sample(1:vcount(l1)),sample(1:vcount(l1)))
  res_sa <- optim(par = perm, fn = ged2, A = A,adj = adj,gr = genperm_sb2,method = "SANN",
                  control = list(maxit = iter, temp = 100, tmax = iter/1000, trace = FALSE,
                                 REPORT = 5))
  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):length(res_sa$par)]

  res_sa <- optim(par = c(perm1,perm2), fn = ged2, A = A,adj = adj,gr = genperm_sb2,method = "SANN",
                  control = list(maxit = iter/10, temp = 10, tmax = iter/1000, trace = FALSE,
                                 REPORT = 5))

  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):length(res_sa$par)]

  Ared <- perm2box(adj,perm1,perm2)

  Afull <- Ared[MSE,MSE]
  for(grp in unique(MSE)){
    Afull[MSE==grp,MSE==grp] <- B[MSE==grp,MSE==grp]
  }

  res <- list(perm1 = perm1,perm2 = perm2,ged=res_sa$value,A_MSE = Ared,A=Afull,MSE = MSE)
  res
}

create_superbox3 <- function(l,iter=500000){
  B <- as_adj(l,sparse = FALSE)
  MSE <- struct_equi(l)
  l1 <- delete.vertices(l,which(duplicated(MSE)))
  adj <- neighborhood(l1)
  adj <- lapply(adj,function(x) as.integer(x))
  A <- as_adj(l1,sparse=FALSE)
  perm <- c(sample(1:vcount(l1)),sample(1:vcount(l1)),sample(1:vcount(l1)))
  res_sa <- optim(par = perm, fn = ged3, A = A,adj = adj,gr = genperm_sb3,method = "SANN",
                  control = list(maxit = iter, temp = 100, tmax = iter/1000, trace = FALSE,
                                 REPORT = 5))

  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):(2*vcount(l1))]
  perm3 <- res_sa$par[(2*vcount(l1)+1):(3*vcount(l1))]

  res_sa <- optim(par = c(perm1,perm2,perm3), fn = ged3, A = A,adj = adj,gr = genperm_sb3,method = "SANN",
                  control = list(maxit = iter/10, temp = 10, tmax = iter/1000, trace = FALSE, REPORT = 5))

  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):(2*vcount(l1))]
  perm3 <- res_sa$par[(2*vcount(l1)+1):(3*vcount(l1))]

  Ared <- perm2cube(adj,perm1,perm2,perm3)

  Afull <- Ared[MSE,MSE]
  for(grp in unique(MSE)){
    Afull[MSE==grp,MSE==grp] <- B[MSE==grp,MSE==grp]
  }

  res <- list(perm1 = perm1,perm2 = perm2,perm3=perm3,ged=res_sa$value,A_MSE = Ared,A=Afull,MSE = MSE)
  res
}

create_superbox4 <- function(l,iter=500000){
  B <- as_adj(l,sparse = FALSE)
  MSE <- struct_equi(l)
  l1 <- delete.vertices(l,which(duplicated(MSE)))
  adj <- neighborhood(l1)
  adj <- lapply(adj,function(x) as.integer(x))
  A <- as_adj(l1,sparse=FALSE)
  perm <- c(sample(1:vcount(l1)),sample(1:vcount(l1)),sample(1:vcount(l1)),sample(1:vcount(l1)))
  res_sa <- optim(par = perm, fn = ged4, A = A,adj = adj,gr = genperm_sb4,method = "SANN",
                  control = list(maxit = iter, temp = 100, tmax = iter/1000, trace = FALSE,
                                 REPORT = 5))

  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):(2*vcount(l1))]
  perm3 <- res_sa$par[(2*vcount(l1)+1):(3*vcount(l1))]
  perm4 <- res_sa$par[(3*vcount(l1)+1):(4*vcount(l1))]

  res_sa <- optim(par = c(perm1,perm2,perm3,perm4), fn = ged4, A = A,adj = adj,gr = genperm_sb4,method = "SANN",
                  control = list(maxit = iter/10, temp = 10, tmax = iter/1000, trace = FALSE, REPORT = 5))

  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):(2*vcount(l1))]
  perm3 <- res_sa$par[(2*vcount(l1)+1):(3*vcount(l1))]
  perm4 <- res_sa$par[(3*vcount(l1)+1):(4*vcount(l1))]

  Ared <- perm2hyper(adj,perm1,perm2,perm3,perm4)

  Afull <- Ared[MSE,MSE]
  for(grp in unique(MSE)){
    Afull[MSE==grp,MSE==grp] <- B[MSE==grp,MSE==grp]
  }

  res <- list(perm1 = perm1,perm2 = perm2,perm3=perm3,perm4=perm4,ged=res_sa$value,A_MSE = Ared,A=Afull,MSE = MSE)
  res
}

create_superbox5 <- function(l){
  B <- as_adj(l,sparse = FALSE)
  MSE <- struct_equi(l)
  l1 <- delete.vertices(l,which(duplicated(MSE)))
  adj <- neighborhood(l1)
  adj <- lapply(adj,function(x) as.integer(x))
  A <- as_adj(l1,sparse=FALSE)
  perm <- c(sample(1:vcount(l1)),sample(1:vcount(l1)),
            sample(1:vcount(l1)),sample(1:vcount(l1)),sample(1:vcount(l1)))
  res_sa <- optim(par = perm, fn = ged5, A = A,adj = adj,gr = genperm_sb5,method = "SANN",
                  control = list(maxit = iter, temp = 100, tmax = iter/1000, trace = FALSE,
                                 REPORT = 5))

  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):(2*vcount(l1))]
  perm3 <- res_sa$par[(2*vcount(l1)+1):(3*vcount(l1))]
  perm4 <- res_sa$par[(3*vcount(l1)+1):(4*vcount(l1))]
  perm5 <- res_sa$par[(4*vcount(l1)+1):(5*vcount(l1))]

  res_sa <- optim(par = c(perm1,perm2,perm3,perm4,perm5), fn = ged5, A = A,adj = adj,gr = genperm_sb5,method = "SANN",
                  control = list(maxit = iter/10, temp = 10, tmax = iter/1000, trace = FALSE, REPORT = 5))

  perm1 <- res_sa$par[1:vcount(l1)]
  perm2 <- res_sa$par[(vcount(l1)+1):(2*vcount(l1))]
  perm3 <- res_sa$par[(2*vcount(l1)+1):(3*vcount(l1))]
  perm4 <- res_sa$par[(3*vcount(l1)+1):(4*vcount(l1))]
  perm5 <- res_sa$par[(4*vcount(l1)+1):(5*vcount(l1))]

  Ared <- perm2pente(adj,perm1,perm2,perm3,perm4,perm5)

  Afull <- Ared[MSE,MSE]
  for(grp in unique(MSE)){
    Afull[MSE==grp,MSE==grp] <- B[MSE==grp,MSE==grp]
  }

  res <- list(perm1 = perm1,perm2 = perm2,perm3=perm3,perm4=perm4,perm5=perm5,
              ged=res_sa$value,A_MSE = Ared,A=Afull,MSE = MSE)
  res
}

# create interval graph from permutation
perm2int <- function(N,perm){

  # xy <- getxy_cpp(N,perm)
  Nperm <- lapply(N,function(x) perm[x])
  xy <- cbind(sapply(Nperm,function(x) c(min(x))),perm)
  a <- xy[,1]
  b <- xy[,2]

  A <- getA_cpp(a,b)
  diag(A) <- 0
  A
}

perm2box <- function(N,perm1,perm2){

  # xy1 <- getxy_cpp(N,perm1)
  # xy2 <- getxy_cpp(N,perm2)
  Nperm <- lapply(N,function(x) perm1[x])
  xy1 <- cbind(sapply(Nperm,function(x) c(min(x))),perm1)

  Nperm <- lapply(N,function(x) perm2[x])
  xy2   <- cbind(sapply(Nperm,function(x) c(min(x))),perm2)

  a <- xy1[,1]
  b <- xy1[,2]

  A1 <- getA_cpp(a,b)

  a <- xy2[,1]
  b <- xy2[,2]

  A2 <- getA_cpp(a,b)

  diag(A1) <- diag(A2) <- 0

  (A1==1 & A1==A2)+0
}

perm2cube <- function(N,perm1,perm2,perm3){

  # xy1 <- getxy_cpp(N,perm1)
  # xy2 <- getxy_cpp(N,perm2)
  # xy3 <- getxy_cpp(N,perm3)
  Nperm <- lapply(N,function(x) perm1[x])
  xy1 <- cbind(sapply(Nperm,function(x) c(min(x))),perm1)

  Nperm <- lapply(N,function(x) perm2[x])
  xy2   <- cbind(sapply(Nperm,function(x) c(min(x))),perm2)

  Nperm <- lapply(N,function(x) perm3[x])
  xy3   <- cbind(sapply(Nperm,function(x) c(min(x))),perm3)

  a <- xy1[,1]
  b <- xy1[,2]
  A1 <- getA_cpp(a,b)

  a <- xy2[,1]
  b <- xy2[,2]
  A2 <- getA_cpp(a,b)

  a <- xy3[,1]
  b <- xy3[,2]
  A3 <- getA_cpp(a,b)

  diag(A1) <- diag(A2) <- diag(A3) <- 0

  (A1==1 & A1==A2 & A2==A3)+0
}

perm2hyper <- function(N,perm1,perm2,perm3,perm4){

  # xy1 <- getxy_cpp(N,perm1)
  # xy2 <- getxy_cpp(N,perm2)
  # xy3 <- getxy_cpp(N,perm3)
  # xy4 <- getxy_cpp(N,perm4)

  Nperm <- lapply(N,function(x) perm1[x])
  xy1 <- cbind(sapply(Nperm,function(x) c(min(x))),perm1)

  Nperm <- lapply(N,function(x) perm2[x])
  xy2   <- cbind(sapply(Nperm,function(x) c(min(x))),perm2)

  Nperm <- lapply(N,function(x) perm3[x])
  xy3   <- cbind(sapply(Nperm,function(x) c(min(x))),perm3)

  Nperm <- lapply(N,function(x) perm4[x])
  xy4   <- cbind(sapply(Nperm,function(x) c(min(x))),perm4)


  a <- xy1[,1]
  b <- xy1[,2]
  A1 <- getA_cpp(a,b)

  a <- xy2[,1]
  b <- xy2[,2]
  A2 <- getA_cpp(a,b)

  a <- xy3[,1]
  b <- xy3[,2]
  A3 <- getA_cpp(a,b)

  a <- xy4[,1]
  b <- xy4[,2]
  A4 <- getA_cpp(a,b)

  diag(A1) <- diag(A2) <- diag(A3) <- diag(A4) <- 0

  (A1==1 & A1==A2 & A2==A3 & A3==A4)+0
}

perm2pente <- function(N,perm1,perm2,perm3,perm4,perm5){

  # xy1 <- getxy_cpp(N,perm1)
  # xy2 <- getxy_cpp(N,perm2)
  # xy3 <- getxy_cpp(N,perm3)
  # xy4 <- getxy_cpp(N,perm4)
  # xy5 <- getxy_cpp(N,perm5)

  Nperm <- lapply(N,function(x) perm1[x])
  xy1 <- cbind(sapply(Nperm,function(x) c(min(x))),perm1)

  Nperm <- lapply(N,function(x) perm2[x])
  xy2   <- cbind(sapply(Nperm,function(x) c(min(x))),perm2)

  Nperm <- lapply(N,function(x) perm3[x])
  xy3   <- cbind(sapply(Nperm,function(x) c(min(x))),perm3)

  Nperm <- lapply(N,function(x) perm4[x])
  xy4   <- cbind(sapply(Nperm,function(x) c(min(x))),perm4)

  Nperm <- lapply(N,function(x) perm4[x])
  xy5   <- cbind(sapply(Nperm,function(x) c(min(x))),perm5)

  a <- xy1[,1]
  b <- xy1[,2]
  A1 <- getA_cpp(a,b)

  a <- xy2[,1]
  b <- xy2[,2]
  A2 <- getA_cpp(a,b)

  a <- xy3[,1]
  b <- xy3[,2]
  A3 <- getA_cpp(a,b)

  a <- xy4[,1]
  b <- xy4[,2]
  A4 <- getA_cpp(a,b)

  a <- xy5[,1]
  b <- xy5[,2]
  A5 <- getA_cpp(a,b)


  diag(A1) <- diag(A2) <- diag(A3) <- diag(A4) <- diag(A5) <- 0


  (A1==1 & A1==A2 & A2==A3 & A3==A4 & A4==A5)+0
}


# get the interval representation of a permutation
perm2rep <- function(N,perm){
  Nperm <- lapply(N,function(x) perm[x])
  xy <- cbind(sapply(Nperm,function(x) c(min(x))),perm)
  colnames(xy) <- NULL
  xy
}

#generate new permutation for lazarus count problem
genperm_lz <- function(sq,A) {
  changepoints <- sample(sq,2)
  tmp <- sq[changepoints[1]]
  sq[changepoints[1]] <- sq[changepoints[2]]
  sq[changepoints[2]] <- tmp
  sq
}


#generate new permutation for super interval problem
genperm_sb1 <- function(sq,A,adj) {
  changepoints <- sample(sq,2)
  tmp <- sq[changepoints[1]]
  sq[changepoints[1]] <- sq[changepoints[2]]
  sq[changepoints[2]] <- tmp
  sq
}

#generate new permutation for super box2 problem
genperm_sb2 <- function(sq,A,adj) {
  n <- length(sq)/2
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):length(sq)]

  changepoints <- sample(perm1,2)
  tmp <- perm1[changepoints[1]]
  perm1[changepoints[1]] <- perm1[changepoints[2]]
  perm1[changepoints[2]] <- tmp

  changepoints <- sample(perm2,2)
  tmp <- perm2[changepoints[1]]
  perm2[changepoints[1]] <- perm2[changepoints[2]]
  perm2[changepoints[2]] <- tmp
  sq <- c(perm1,perm2)
  sq
}

#generate new permutation for super box3 problem
genperm_sb3 <- function(sq,A,adj) {
  n <- length(sq)/3
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):(2*n)]
  perm3 <- sq[(2*n+1):length(sq)]

  changepoints <- sample(perm1,2)
  tmp <- perm1[changepoints[1]]
  perm1[changepoints[1]] <- perm1[changepoints[2]]
  perm1[changepoints[2]] <- tmp

  changepoints <- sample(perm2,2)
  tmp <- perm2[changepoints[1]]
  perm2[changepoints[1]] <- perm2[changepoints[2]]
  perm2[changepoints[2]] <- tmp

  changepoints <- sample(perm3,2)
  tmp <- perm3[changepoints[1]]
  perm3[changepoints[1]] <- perm3[changepoints[2]]
  perm3[changepoints[2]] <- tmp

  sq <- c(perm1,perm2,perm3)
  sq
}

#generate new permutation for super box4 problem
genperm_sb4 <- function(sq,A,adj) {
  n <- length(sq)/4
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):(2*n)]
  perm3 <- sq[(2*n+1):(3*n)]
  perm4 <- sq[(3*n+1):(4*n)]

  changepoints <- sample(perm1,2)
  tmp <- perm1[changepoints[1]]
  perm1[changepoints[1]] <- perm1[changepoints[2]]
  perm1[changepoints[2]] <- tmp

  changepoints <- sample(perm2,2)
  tmp <- perm2[changepoints[1]]
  perm2[changepoints[1]] <- perm2[changepoints[2]]
  perm2[changepoints[2]] <- tmp

  changepoints <- sample(perm3,2)
  tmp <- perm3[changepoints[1]]
  perm3[changepoints[1]] <- perm3[changepoints[2]]
  perm3[changepoints[2]] <- tmp

  changepoints <- sample(perm4,2)
  tmp <- perm4[changepoints[1]]
  perm4[changepoints[1]] <- perm4[changepoints[2]]
  perm4[changepoints[2]] <- tmp

  sq <- c(perm1,perm2,perm3,perm4)
  sq
}

#generate new permutation for super box5 problem
genperm_sb5 <- function(sq,A,adj) {
  n <- length(sq)/5
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):(2*n)]
  perm3 <- sq[(2*n+1):(3*n)]
  perm4 <- sq[(3*n+1):(4*n)]
  perm5 <- sq[(4*n+1):(5*n)]

  changepoints <- sample(perm1,2)
  tmp <- perm1[changepoints[1]]
  perm1[changepoints[1]] <- perm1[changepoints[2]]
  perm1[changepoints[2]] <- tmp

  changepoints <- sample(perm2,2)
  tmp <- perm2[changepoints[1]]
  perm2[changepoints[1]] <- perm2[changepoints[2]]
  perm2[changepoints[2]] <- tmp

  changepoints <- sample(perm3,2)
  tmp <- perm3[changepoints[1]]
  perm3[changepoints[1]] <- perm3[changepoints[2]]
  perm3[changepoints[2]] <- tmp

  changepoints <- sample(perm4,2)
  tmp <- perm4[changepoints[1]]
  perm4[changepoints[1]] <- perm4[changepoints[2]]
  perm4[changepoints[2]] <- tmp

  changepoints <- sample(perm5,2)
  tmp <- perm5[changepoints[1]]
  perm5[changepoints[1]] <- perm5[changepoints[2]]
  perm5[changepoints[2]] <- tmp

  sq <- c(perm1,perm2,perm3,perm4,perm5)
  sq
}


#graph edit distance
ged1 <- function(sq,A,adj){
  Anew <- perm2int(adj,sq)
  sum(Anew-A)/2
}

ged2 <- function(sq,A,adj){
  n <- length(sq)/2
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):(2*n)]
  Anew <- perm2box(adj,perm1,perm2)
  sum(Anew-A)/2
}

ged3 <- function(sq,A,adj){
  n <- length(sq)/3
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):(2*n)]
  perm3 <- sq[(2*n+1):length(sq)]
  Anew <- perm2cube(adj,perm1,perm2,perm3)
  sum(Anew-A)/2
}

ged4 <- function(sq,A,adj){
  n <- length(sq)/4
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):(2*n)]
  perm3 <- sq[(2*n+1):(3*n)]
  perm4 <- sq[(3*n+1):(4*n)]
  Anew <- perm2hyper(adj,perm1,perm2,perm3,perm4)
  sum(Anew-A)/2
}

ged5 <- function(sq,A,adj){
  n <- length(sq)/5
  perm1 <- sq[1:n]
  perm2 <- sq[(n+1):(2*n)]
  perm3 <- sq[(2*n+1):(3*n)]
  perm4 <- sq[(3*n+1):(4*n)]
  perm5 <- sq[(4*n+1):(5*n)]
  Anew <- perm2pente(adj,perm1,perm2,perm3,perm4,perm5)
  sum(Anew-A)/2
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
