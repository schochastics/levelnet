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
graph_circ_arc <- function(n,r,pos=1/10,neg=1/10,skew=1){
  pts <- circleFun(r=r,npoints=n)
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

#helper
circleFun <- function(center = c(0,0),r = 1, npoints = 20,skew = 0){
  ttseq <- seq(0,2*pi,length.out = npoints*100)
  if(skew==0){
    tt <- sample(ttseq,npoints)
  } else{
    # probs <- abs(ttseq)/skew
    probs <- abs(sin(ttseq/skew))
    probs <- probs/sum(probs)
    tt <- sample(ttseq,npoints,prob = probs)
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
