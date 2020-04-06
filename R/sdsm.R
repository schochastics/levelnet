#' @title Stochastic Degree Sequence Model
#' @description Stochastic Degree Sequence Model.

#' @param g igraph object. The two-mode network
#' @param proj string. Which mode to project on ("true"/"false")
#' @param model string. which link to be used ('logit','probit','cloglog' or 'scobit')
#' @param max_iter number of randomly sampled networks
#' @param alpha significance level
#' @param params named parameter list for scobit model
#' @param verbose print status during execution

#' @return backbone of one-mode projection
#' @author David Schoch
#' @references Neal, Zachary (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors
#' @export
#'
sdsm <- function(g,proj="true",model="logit",max_iter=1000,alpha=0.05,
                 params=list(b0=0.1,b1=0.00005,b2=0.00005,b3=0.00005,a=0.01),
                 verbose = FALSE){
  if(!igraph::is_bipartite(g)){
    stop("network is not bipartite")
  }
  if(!proj%in%c("true","false")){
    stop("proj must be one of 'true' or 'false'")
  }
  if(any(igraph::edge_attr_names(g)=="weight")){
    stop("sdsm does not work with weighted two-mode networks")
  }
  if(!model%in%c("logit","probit","cloglog","scobit")){
    stop("model must be one of 'logit','probit', 'cloglog' or 'scobit'")
  }
  if(!any(igraph::vertex_attr_names(g)=="name")){
    igraph::V(g)$name <- 1:igraph::vcount(g)
  }
  proj_logi <- as.logical(proj)

  bip <- igraph::bipartite_projection(g,which=proj)
  P <- as.matrix(igraph::get.adjacency(bip,attr = "weight",sparse=T))

  deg_artif <- unname(igraph::degree(g)[igraph::V(g)$type!=proj_logi])
  deg_agent <- unname(igraph::degree(g)[igraph::V(g)$type==proj_logi])
  A <- igraph::as_incidence_matrix(g,sparse = F)
  if(nrow(A)!=sum(igraph::V(g)$type==proj_logi)){
    A <- t(A)
  }
  D_agent <- rep(deg_agent,length(deg_artif))
  D_artif <- rep(deg_artif,each=length(deg_agent))
  resp <- c(A)
  df <- data.frame(resp,D_agent,D_artif)
  if(verbose){
    message(paste0("    Fitting ",model," model\n"))
  }
  if(model!="scobit"){
    model.fit <- stats::glm(resp ~ D_agent + D_artif + D_agent * D_artif,
                     family = stats::binomial(link = model),data = df)
    df[["prob"]] <- stats::predict(model.fit, newdata = df, type = "response")
  } else{
    y  <- resp
    x1 <- D_agent
    x2 <- D_artif

    model.fit <- stats::optim(params,scobit_loglike_cpp,
                              gr=scobit_loglike_gr_cpp,
                              method="BFGS",
                              x1=x1,x2=x2,y=y)

    pars <- c(model.fit$par[1],model.fit$par[2],model.fit$par[3],model.fit$par[4])

    df[["prob"]] <- scobit_fct(D_agent,D_artif,pars,model.fit$par[5])
  }

  P_test <- matrix(0,length(deg_agent),length(deg_agent))
  if(verbose){
    message("    Simulating random networks\n")
    pb <- utils::txtProgressBar(min = 1, max = max_iter, style = 3)
  }
  for(i in 1:max_iter){
    if(verbose){
      utils::setTxtProgressBar(pb,i)
    }
    b_vec <- stats::runif(length(deg_agent)*length(deg_artif))
    B <- Matrix::Matrix((b_vec<=df[["prob"]])+0,length(deg_agent),length(deg_artif),sparse=T)
    # P_rand <- eigenMatMult(B)
    P_rand <- Matrix::tcrossprod(B)
    P_test <- P_test + (P_rand >= P)+0
  }
  if(verbose){
    close(pb)
  }
  A_new <- as.matrix((P_test<=alpha*max_iter)+0)
  l <- igraph::graph_from_adjacency_matrix(A_new,mode = "undirected",diag = F)
  igraph::V(l)$name <- igraph::V(bip)$name
  return(l)
}

#' @title sdsm model diagnostics
#' @description check which binary outcome model fits the data best
#'
#' @param g igraph object. The two-mode network
#' @param proj string. Which mode to project on
#' @param iter number of fits per model
#' @param verbose logical. print additional information (default: FALSE)
#' @param params named parameter list for scobit model
#' @return rmse and runtime of various models
#' @author David Schoch
#' @importFrom graphics abline par plot
#' @export
#'
sdsm_diagnostic <- function(g,proj="true",iter=10,verbose=FALSE,
                            params=list(b0=0.1,b1=0.00005,b2=0.00005,b3=0.00005,a=0.01)){

  if(!igraph::is_bipartite(g)){
    stop("network is not bipartite")
  }
  if(!proj%in%c("true","false")){
    stop("proj must be one of 'true' or 'false'")
  }
  if(any(igraph::edge_attr_names(g)=="weight")){
    stop("sdsm does not work with weighted two-mode networks")
  }
  if(!any(igraph::vertex_attr_names(g)=="name")){
    igraph::V(g)$name <- 1:igraph::vcount(g)
  }
  proj_logi <- as.logical(proj)

  bip <- igraph::bipartite_projection(g,which=proj)
  P <- as.matrix(igraph::get.adjacency(bip,attr = "weight",sparse=T))


  deg_artif <- unname(igraph::degree(g)[igraph::V(g)$type!=proj_logi])
  deg_agent <- unname(igraph::degree(g)[igraph::V(g)$type==proj_logi])
  A <- igraph::as_incidence_matrix(g,sparse = F)
  if(nrow(A)!=sum(igraph::V(g)$type==proj_logi)){
    A <- t(A)
  }
  D_agent <- rep(deg_agent,length(deg_artif))
  D_artif <- rep(deg_artif,each=length(deg_agent))
  resp <- c(A)
  df <- data.frame(resp,D_agent,D_artif)
  models <- c("logit","probit","cloglog")
  df_mod <- data.frame(name=c(models,"scobit"),rmse_row=rep(0,length(models)+1),
                       rmse_col=rep(0,length(models)+1),time=rep(0,length(models)+1))
  bkp <- par()$mfcol
  par(mfcol=c(2,4))
  for(m in seq_along(models)){
    if(verbose){
      cat(paste0("fitting ",models[m]," model\n"))
    }
    df_mod$time[m] <- system.time(model.fit <- stats::glm(resp ~ D_agent + D_artif + D_agent * D_artif,
                            family = stats::binomial(link = models[m]),data = df))[[3]]

    df[["prob"]] <- stats::predict(model.fit, newdata = df, type = "response")


    for(i in 1:iter){
      b_vec <- stats::runif(length(deg_agent)*length(deg_artif))
      B <- Matrix::Matrix((b_vec<=df[["prob"]])+0,length(deg_agent),length(deg_artif),sparse=T)
      res <- c(sqrt(sum((rowSums(A)-Matrix::rowSums(B))^2)/nrow(A)),
        sqrt(sum((colSums(A)-Matrix::colSums(B))^2)/ncol(A)))
      df_mod$rmse_row[m] <- df_mod$rmse_row[m]+res[1]
      df_mod$rmse_col[m] <- df_mod$rmse_col[m]+res[2]
    }
    plot(rowSums(A),Matrix::rowSums(B),main=paste0(models[m], "(rows)"),xlab="data",ylab="sample")
    abline(0,1)
    plot(colSums(A),Matrix::colSums(B),main=paste0(models[m], "(cols)"),xlab="data",ylab="sample")
    abline(0,1)
  }
  if(verbose){
    message("fitting scobit model\n")
  }
  y  <- resp
  x1 <- D_agent
  x2 <- D_artif

  df_mod$time[m+1] <- system.time(model.fit <- stats::optim(params,scobit_loglike_cpp,gr=scobit_loglike_gr_cpp,method="BFGS",x1=x1,x2=x2,y=y))[[3]]
  pars <- c(model.fit$par[1],model.fit$par[2],model.fit$par[3],model.fit$par[4])
  df[["prob"]] <- scobit_fct(D_agent,D_artif,pars,model.fit$par[5])
  for(i in 1:iter){
    b_vec <- stats::runif(length(deg_agent)*length(deg_artif))
    B <- Matrix::Matrix((b_vec<=df[["prob"]])+0,length(deg_agent),length(deg_artif),sparse=T)
    res <- c(sqrt(sum((rowSums(A)-Matrix::rowSums(B))^2)/nrow(A)),
             sqrt(sum((colSums(A)-Matrix::colSums(B))^2)/ncol(A)))
    df_mod$rmse_row[m+1] <- df_mod$rmse_row[m+1]+res[1]
    df_mod$rmse_col[m+1] <- df_mod$rmse_col[m+1]+res[2]
  }
  plot(rowSums(A),Matrix::rowSums(B),main=paste0("scobit", "(rows)"),xlab="data",ylab="sample")
  abline(0,1)
  plot(colSums(A),Matrix::colSums(B),main=paste0("scobit", "(cols)"),xlab="data",ylab="sample")
  abline(0,1)
  par(mfcol=bkp)
  df_mod$rmse_row <- df_mod$rmse_row/iter
  df_mod$rmse_col <- df_mod$rmse_col/iter
  df_mod
}


###########################################
# scobit helper ----
###########################################
scobit_fct <- function(x1,x2,beta,alpha){
  fct <- 1-1/(1+exp(beta[1]+beta[2]*x1+beta[3]*x2+beta[4]*x1*x2))^alpha
  fct
}


