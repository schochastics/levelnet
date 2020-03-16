#' @title Random Interval Graph
#' @description  Create a random interval graph. In an interval graph, each node is characterized by an
#' interval on the real line. Two nodes are connected, if their intervals overlap.
#' @param n number of nodes
#' @param r radius (see details)
#' @param sd standard deviation (see details)
#' @details Interval graphs are created as follows. First, n random points x are created uniformly at random between 0 and `r`.
#' For each point, a value Y is created from a normal distribution with mean X and standard deviation is `sd`. In this way, it is
#' possible to control the density of the network. The larger `r` and the larger `sd` the more likely do intervals overlap.
#' @return interval graph as igraph object and interval representation as node attribute (a,b)
#' @seealso [graph_indifference,graph_tolerance]
#' @author David Schoch
#' @examples
#' graph_interval(n = 10)
#' @export

graph_interval <- function(n,r = 2,sd = 0.5){
  if(missing(n)){
    stop('argument "n" is missing with no default')
  }
  x <- stats::runif(n,0,r)
  y <- stats::rnorm(n,x,sd)
  xy <- matrix(c(pmin(x,y),pmax(x,y)),ncol=2)
  a <- xy[,1]
  b <- xy[,2]

  A <- ((outer(a,a,">=") & outer(a,b,"<=")) | (outer(a,b,"<=") & outer(b,b,">="))) + 0

  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A,"undirected")
  igraph::V(g)$a <- a
  igraph::V(g)$b <- b
  g$name <- "Interval Graph"
  return(g)
}

#' @title Random Indifference Graph
#' @description  Create a random indifference graph. An indifference graph is an interval graph where intervals
#' have length 1.
#' @param n number of nodes
#' @param r radius
#' @details `n` points (x) are sampled uniformly at random between 0 and `r`. The interval is then given by (x,x+1)
#' @return indifference graph as igraph object and interval representation (a,b)
#' @seealso [graph_interval,graph_tolerance]
#' @author David Schoch
#' @examples
#' graph_indifference(n = 10)
#' @export

graph_indifference <- function(n,r = 2){
  if(missing(n)){
    stop('argument "n" is missing with no default')
  }
  x <- stats::runif(n,0,r)
  A <- outer(x,x,function(x,y) abs(x-y)<=1)+0
  xy <-  cbind(x,x+rep(1,n))
  colnames(xy) <- NULL
  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A,"undirected")
  igraph::V(g)$a <- xy[,1]
  igraph::V(g)$b <- xy[,2]
  g$name <- "Indifference Graph"
  return(g)
}

#' @title Random Tolerance Graph
#' @description  Create a random tolerance graph. A tolerance graph is an interval graph, where nodes are
#' only connected if the overlap is larger than a nodes tolerance level. These graphs are directed.
#'
#' @param n number of nodes
#' @param r radius (see details)
#' @param sd standard deviation (see details)
#' @param tol tolerance
#' @details Tolerance graphs are created as follows. First, n random points x are created uniformly at random between 0 and `r`.
#' For each point, a value Y is created from a normal distribution with mean X and standard deviation is `sd`. In this way, it is
#' possible to control the density of the network. The larger `r` and the larger `sd` the more likely do intervals overlap. When overlaps
#' are calculated, it is checked whether the overlap is larger than the tolerance of the node. If so, the edge is included.
#' @return tolerance graph as igraph object and interval representation and tolerance as node attributes
#' @seealso [graph_interval,graph_indifference]
#' @author David Schoch
#' @examples
#' graph_tolerance(n = 10)
#' @export
#
graph_tolerance <- function(n,r = 2,sd = 0.5,tol = 0.5){
  if(missing(n)){
    stop('argument "n" is missing with no default')
  }
  x <- stats::runif(n,0,r)
  y <- stats::rnorm(n,x,sd)
  tol <- abs(stats::rnorm(n,0,tol))
  xy <- matrix(c(pmin(x,y),pmax(x,y)),ncol=2)
  colnames(xy) <- NULL
  a <- xy[,1]
  b <- xy[,2]
  A <- matrix(0,n,n)
  # TODO: properly vectorize
  # (outer(a,a,">=") & outer(a,b,"<=") & outer(-a+b,tol,">")) +
  #   (outer(a,a,"<=") & outer(b,a,">=") & outer(b,b,"<=") & outer(b-a,tol,">")) +
  #   (outer(a,a,">=") & outer(b,b,"<=") & b-a>=tol) +
  #   (outer(a,a,"<=") & outer(b,b,">=") & b-a>=tol)
  for(i in 1:n){
    for(j in 1:n){
      if(a[i]>=a[j] & a[i]<=b[j] ){
        if(b[j]-a[i]>tol[i]){
          A[i,j] <- 1
        }
      } else if(a[i]<=a[j] & b[i]>=a[j] & b[i]<=b[j]){
        if(b[i]-a[j]>tol[i]){
          A[i,j] <- 1
        }
      } else if(a[i]>=a[j] & b[i]<=b[j]){
        if(b[i]-a[i]>=tol[i]){
          A[i,j] <- 1
        }
      } else if(a[i]<=a[j] & b[i]>=b[j]){
        if(b[j]-a[j]>=tol[i]){
          A[i,j] <- 1
        }
      }
    }
  }
  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A,"directed")
  igraph::V(g)$a <- xy[,1]
  igraph::V(g)$b <- xy[,2]
  igraph::V(g)$tolerance <- tol

  g$name <- "Tolerance Graph"
  return(g)
}

#' @title Boxicity 2 graph
#' @description  Create a random graph with boxicity 2.
#'
#' @param n number of nodes
#' @param r radius
#' @param sd standard deviation
#' @return Boxicity 2 graph as igraph object
#' @author David Schoch
#' @export
#

graph_rectangle <- function(n,r=2,sd=0.5){
  x1 <- stats::runif(n,0,r)
  x2 <- stats::rnorm(n,x1,sd)
  y1 <- stats::runif(n,0,r)
  y2 <- stats::rnorm(n,y1,sd)
  x <-  t(apply(cbind(x1,x2),1,sort))
  y <-  t(apply(cbind(y1,y2),1,sort))
  A <- matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      XA1 <- x[i,1]
      XA2 <- x[i,2]
      YA1 <- y[i,1]
      YA2 <- y[i,2]
      XB1 <- x[j,1]
      XB2 <- x[j,2]
      YB1 <- y[j,1]
      YB2 <- y[j,2]
      SI <-  max(0, min(XA2, XB2) - max(XA1, XB1)) * max(0, min(YA2, YB2) - max(YA1, YB1))
      if(SI>0){
        A[i,j] <- 1
      }
    }
  }
  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A,"undirected")
  rect <- cbind(x,y)
  igraph::V(g)$x1 <- rect[,1]
  igraph::V(g)$x2 <- rect[,2]
  igraph::V(g)$y1 <- rect[,3]
  igraph::V(g)$y2 <- rect[,4]
  return(g)
}

#' @title two-mode network from a data.frame
#' @description Create a two-mode network from a data.frame
#'
#' @param df data.frame
#' @param type1 column name of mode 1
#' @param type2 column name of mode 2
#' @return two mode network as igraph object
#' @author David Schoch
#' @export

bipartite_from_data_frame <- function(df,type1,type2){
  mode1 <- unique(df[[type1]])
  mode2 <- unique(df[[type2]])
  el <- cbind(df[[type1]],df[[type2]])

  g <- igraph::graph.empty(directed=F)
  g <- igraph::add_vertices(g,nv=length(mode1),attr = list(name=mode1,type=T))
  g <- igraph::add_vertices(g,nv=length(mode2),attr = list(name=mode2,type=F))
  g <- igraph::add_edges(g,c(t(el)))
  if(igraph::any_multiple(g)){
    igraph::E(g)$weight <- 1
    g <- igraph::simplify(g,remove.multiple = TRUE,
                          remove.loops = TRUE,
                          edge.attr.comb = "sum")
  }
  g
}

#' @title generate random roll-call votes based on ideology space
#' @param M number of members
#' @param D distance between means
#' @param p dimension of space
#' @param pd dimensioms where distributions are separated
#' @param beta scaling parameter for probabilistic voting
#' @param r radius of hypersphere for random generation
#' @param noprob probabilit of non voting
#' @param Nrand number of randomly generated votes
#' @param N number of votes to sample from randomly generated votes
#' @return list with random votes and ideologies
#' @author David Schoch
#' @references Aldrich, John H., and Montgomery, Jacob M., and Sparks, David B. (2014). Polarization and Ideology: Partisan Sources of Low Dimensionality in Scaled Roll Call Analyses. Political Analysis 22:435-456
#' @export

graph_random_vote <- function(M=101,D=1,p=4,pd=2,beta=1,r=9,noprob=0.05,Nrand=1000,N=525){
  #sanity checks
  if(pd>p){
    stop("pd must be smaller than p")
  }
  if(N>Nrand){
    stop("N must be smaller than Nrand")
  }
  if(M<2){
    stop("M should be larger than 1")
  }

  n1 <- ceiling(M/2)
  n2 <- M-n1

  # create random ideologies according to parameters
  xd <- matrix(stats::rnorm(n1*pd,D/2,1),nrow=n1,ncol=pd)
  yd <- matrix(stats::rnorm(n2*pd,-D/2,1),nrow=n2,ncol=pd)
  x <- matrix(stats::rnorm(n1*(p-pd),0,1),nrow=n1,ncol=p-pd)
  y <- matrix(stats::rnorm(n2*(p-pd),0,1),nrow=n2,ncol=p-pd)
  Mnames <- paste0("M",1:M)
  x_tbl <- data.frame(cbind(xd,x),stringsAsFactors = FALSE)
  x_tbl[["party"]] <- "D"
  x_tbl[["name"]] <- Mnames[1:n1]
  y_tbl <- data.frame(cbind(yd,y),stringsAsFactors = FALSE)
  y_tbl[["party"]] <- "R"
  y_tbl[["name"]] <- Mnames[(n1+1):M]

  ideo_tbl <- rbind(x_tbl,y_tbl)
  ideo_mat <-  as.matrix(ideo_tbl[,1:p])
  ideo_tbl <- ideo_tbl[,c(p+2,p+1,1:p)]
  names(ideo_tbl)[3:(p+2)] <- paste0("Dim",1:p)

  sim <- replicate(Nrand,prob_vote(r=r,p=p,ideo_mat))

  idx <- which(colSums(sim)<=(M*0.02) | colSums(sim) >= (M*0.98))
  sim <- sim[,-idx]
  sim <- sim[,sample(1:ncol(sim),N)]

  res <- data.frame(name=rep(paste0("M",1:M),N),
                    vote=c(sim), bill=rep(paste0("B",1:N),each=M),stringsAsFactors = FALSE)
  not_vote <- which(stats::runif(nrow(res))<=noprob)
  res <- res[-not_vote,]
  res[["vote_ident"]] <- paste0(res[["bill"]],"_",res[["vote"]])
  return(list(votes=res,ideo=ideo_tbl))
}


#helper
prob_vote <- function(r,p,x,beta=1){
  a <- random_proposal(r,p)
  b <- random_proposal(r,p)
  wp <- 1/length(a)
  dist_a <- t(apply(x,1,function(y) sqrt((y-a)^2)))
  dist_b <- t(apply(x,1,function(y) sqrt((y-b)^2)))
  if(p==1){
    dist_a <- t(dist_a)
    dist_b <- t(dist_b)
  }
  val <- beta * (rowSums(wp * dist_b^2)-rowSums(wp * dist_a^2))
  probs <- stats::pnorm(val)
  sapply(probs,function(x) sample(c(0,1),1,prob=c(x,1-x)))
}

random_proposal <- function(r,p){
  pts <- stats::runif(p,-r,r)
  while(sqrt(sum(pts^2))>r){
    pts <- stats::runif(p,-r,r)
  }
  pts
}
