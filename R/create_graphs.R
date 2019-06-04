#' @title Interval graph
#' @description  Create a random interval graph.
#'
#' @param n number of nodes
#' @param r radius
#' @param sd standard deviation
#' @return interval graph as igraph object and interval representation
#' @author David Schoch
#' @export

graph_interval <- function(n,r=2,sd=0.5){
  x <- stats::runif(n,0,r)
  y <- stats::rnorm(n,x,sd)
  xy <-  t(apply(cbind(x,y),1,sort))
  colnames(xy) <- NULL
  a <- xy[,1]
  b <- xy[,2]
  A <- matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      if((a[i]>=a[j] & a[i]<=b[j]) | (b[j]>=a[i] & b[j]<=b[i])){
        A[i,j] <- 1
      }
    }
  }
  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A,"undirected")
  igraph::V(g)$a <- xy[,1]
  igraph::V(g)$b <- xy[,2]
  return(g)
}

#' @title Indifference graph
#' @description  Create a random indifference graph.
#'
#' @param n number of nodes
#' @param r radius
#' @return indifference graph as igraph object and interval representation
#' @author David Schoch
#' @export

graph_indifference <- function(n,r=2){
  x <- stats::runif(n,0,r)
  A <- outer(x,x,function(x,y) abs(x-y)<=1)+0
  xy <-  t(apply(cbind(x,x+rep(1,n)),1,sort))
  colnames(xy) <- NULL
  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A,"undirected")
  igraph::V(g)$a <- xy[,1]
  igraph::V(g)$b <- xy[,2]
  return(g)
}

#' @title Tolerance graph
#' @description  Create a random interval graph.
#'
#' @param n number of nodes
#' @param r radius
#' @param sd standard deviation
#' @param tol tolerance
#' @return interval graph as igraph object and interval representation
#' @author David Schoch
#' @export
#
graph_tolerance <- function(n,r=2,sd=0.5,tol=0.5){
  x <- stats::runif(n,0,r)
  y <- stats::rnorm(n,x,sd)
  tol <- abs(stats::rnorm(n,0,tol))
  xy <-  t(apply(cbind(x,y),1,sort))
  colnames(xy) <- NULL
  a <- xy[,1]
  b <- xy[,2]
  A <- matrix(0,n,n)
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

#' @title circular arc graph
#' @description  circular arc graph with positive and negative ties.
#'
#' @param n number of nodes
#' @param r radius
#' @param pos distance fraction between positive edges
#' @param neg distance fraction between negative edges
#' @param skew create clustered graph
#' @return circular arc graph as igraph object
#' @author David Schoch
#' @export
#
graph_circ_arc <- function(n,r,pos=1/10,neg=1/10,skew=0){
  pts <- circleFun(r=r,npoints=n,skew = skew)
  # D1 <- matrix(0,n,n)
  # for(i in 1:n){
  #   for(j in 1:n){
  #     D1[i,j] <- c(arc_dist(pts[i,],pts[j,],r))$x
  #   }
  # }
  D <- arcDistMat(as.matrix(pts),r)
  # print(D1-D)
  thr <- (2*pi*r)*pos
  anti <- arc_dist(c(0,r),c(0,-r),r)*(1-neg)
  P <- (D<=thr & D!=0)+0
  N <- (D>=anti & D!=0)+0

  A <- P-N
  g <- igraph::graph_from_adjacency_matrix(A,mode="undirected",weighted = T)
  igraph::E(g)$type <- ifelse(igraph::E(g)$weight==1,"pos","neg")
  g <- igraph::delete_edge_attr(g,"weight")
  igraph::V(g)$x <- pts$x
  igraph::V(g)$y <- pts$y
  g
}

#' @title k partite graphs
#' @description  Create a random k-partite graph.
#'
#' @param n number of nodes
#' @param grp vector of partition sizes
#' @return kpartite graph
#' @author David Schoch
#' @export

graph_kpartite <- function(n=10,grp=c(5,5)){
  g <- igraph::graph.empty(n=n,directed=FALSE)
  cur_node <- 1
  nodes <- 1:n
  for(i in 1:(length(grp)-1)){
    add_nodes <- cur_node:(cur_node + grp[i] - 1)
    add_edges <- c(t(expand.grid(add_nodes,nodes[nodes>max(add_nodes)])))
    g <- igraph::add.edges(g,add_edges)
    cur_node <- cur_node + grp[i]
  }
  return(g)
}

#' @title convert igraph object to sage format
#' @description  convert igraph object to sage format
#'
#' @param g igraph object
#' @details
#' sage code
#'
#' gis = g.is_interval(certificate=true)
#'
#'gis
#'
#'gisstr=str(gis)
#'
#'o = open('interval_raw.txt','w')
#'
#'o.write(gisstr)
#'
#'o.close()
#'
#'system('grep -oe "[0-9]\\+,\\s[0-9]\\+" interval_raw.txt > intervals.txt')
#' @return sage string
#' @author David Schoch
#' @export

graph_to_sage <- function(g){
  igraph::V(g)$name <- 1:igraph::vcount(g)
  tst <- igraph::get.adjlist(g)
  gstr <- sapply(1:igraph::vcount(g),function(x)paste0(x,":","[",paste(tst[[x]],collapse = ","),"]"))
  gstr <- paste(gstr,collapse = ",")
  gstr <- paste0("g=Graph({",gstr,"})")
  gstr
}

#' @title two-mode network from a data.frame
#' @description Create a two-mode network from a data.frame
#'
#' @param df data.frame
#' @param type1 column name of mode 1
#' @param type2 column name of mode 2
#' @return two mode network
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
circleFun <- function(center = c(0,0),r = 1, npoints = 20,skew = 0){
  ttseq <- seq(0,2*pi,length.out = npoints*100)
  if(skew==0){
    tt <- sample(ttseq,npoints)
  } else{
    probs <- abs(sin(ttseq/skew))
    probs <- probs/sum(probs)
    tt <- sample(ttseq[order(probs)[1:(npoints*5)]],npoints)
  }

  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

arc_dist <- function(x,y,r){
  c <- sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)
  theta <- acos((2*r^2-c^2)/(2*r^2))
  2*pi*r*theta/(2*pi)
}

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
