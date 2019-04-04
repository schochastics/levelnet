#' @title Stochastic Degree Sequence Model
#' @description Stochastic Degree Sequence Model.

#' @param g igraph object. The two-mode network
#' @param proj string. Which mode to project on
#' @param model string. which link to be used ('logit','probit','cloglog' or 'scobit')
#' @param max_iter number of randomly sampled networks
#' @param alpha significance level
#' @param method optimization method used for scobit
#' @param params named parameter list for scobit model
#' @return backbone of one-mode projection
#' @author David Schoch
#' @references Neal, Zachary (2014). The backbone of bipartite projections: Inferring relationships from co-authorship, co-sponsorship, co-attendance and other co-behaviors
#' @export
#'
sdsm <- function(g,proj="true",model="logit",max_iter=1000,alpha=0.05,
                 method = "BFGS",
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
  if(!model%in%c("logit","probit","cloglog","scobit")){
    stop("model must be one of 'logit','probit', 'cloglog' or 'scobit'")
  }
  if(!any(igraph::vertex_attr_names(g)=="name")){
    igraph::V(g)$name <- 1:igraph::vcount(g)
  }
  proj_logi <- as.logical(proj)

  bip <- igraph::bipartite_projection(g,which=proj)
  P <- igraph::get.adjacency(bip,attr = "weight",sparse=F)


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
  cat(paste0("    Fitting ",model," model\n"))
  if(model!="scobit"){
    model.fit <- stats::glm(resp ~ D_agent + D_artif + D_agent * D_artif,
                     family = stats::binomial(link = model),data = df)
    df[["prob"]] <- stats::predict(model.fit, newdata = df, type = "response")
  } else{
    y  <- resp
    x1 <- D_agent
    x2 <- D_artif

    scobit_loglike <- function(b0,b1,b2,b3,a){
      -sum((1 - y)*log((1 + exp(b0+b1*x1+b2*x2+b3*x1*x2))^(-a)) +
             y*log(1-(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))^(-a)))
    }

    scobit_mle <- function(y,x1,x2){
      mle.res <- stats4::mle(scobit_loglike, start = params,
                             method=method)

      beta <- stats4::coef(mle.res)[1:4]
      alpha <- stats4::coef(mle.res)[5]
      list(beta=beta,alpha=alpha)
    }

    model.fit <- scobit_mle(y,x1,x2)

    df[["prob"]] <- scobit_fct(D_agent,D_artif,model.fit$beta,model.fit$alpha)
  }

  # print(summary(model.fit))
  P_test <- matrix(0,length(deg_agent),length(deg_agent))

  cat("    Simulating random networks\n")
  pb <- utils::txtProgressBar(min = 1, max = max_iter, style = 3)
  for(i in 1:max_iter){
    utils::setTxtProgressBar(pb,i)
    b_vec <- stats::runif(length(deg_agent)*length(deg_artif))
    B <- Matrix::Matrix((b_vec<=df[["prob"]])+0,length(deg_agent),length(deg_artif),sparse=T)
    P_rand <- B%*%Matrix::t(B)
    P_test <- P_test + (P_rand >= P)+0
  }
  close(pb)
  A_new <- as.matrix((P_test<=alpha*max_iter)+0)
  l <- igraph::graph_from_adjacency_matrix(A_new,mode = "undirected",diag = F)
  igraph::V(l)$name <- igraph::V(bip)$name
  return(l)
}

#' @title sdsm model diagnostics
#' @description check which link function to use
#'
#' @param g igraph object. The two-mode network
#' @param proj string. Which mode to project on
#' @param iter number of fits per model
#' @param verbose logical. print additional information (default: FALSE)
#' @param method optimization method used for scobit
#' @param params named parameter list for scobit model
#' @return rmse and runtime of various models
#' @author David Schoch
#' @importFrom graphics abline par plot
#' @export
#'
sdsm_diagnostic <- function(g,proj="true",iter=10,verbose=FALSE,
                            method = "BFGS",
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
  P <- igraph::get.adjacency(bip,attr = "weight",sparse=F)


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
  par(mfcol=c(2,4))
  for(m in seq_along(models)){
    if(verbose){
      cat(paste0("fitting ",models[m]," model\n"))
    }
    model.fit <- stats::glm(resp ~ D_agent + D_artif + D_agent * D_artif,
                            family = stats::binomial(link = models[m]),data = df)
    df_mod$time[m] <- system.time(
      df[["prob"]] <- stats::predict(model.fit, newdata = df, type = "response")
    )[[3]]

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
    cat("fitting scobit model\n")
  }
  y  <- resp
  x1 <- D_agent
  x2 <- D_artif

  scobit_loglike <- function(b0,b1,b2,b3,a){
    -sum((1 - y)*log((1 + exp(b0+b1*x1+b2*x2+b3*x1*x2))^(-a)) +
           y*log(1-(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))^(-a)))
  }

  scobit_mle <- function(y,x1,x2){
    mle.res <- stats4::mle(scobit_loglike, start = params,method = method)

    beta <- stats4::coef(mle.res)[1:4]
    alpha <- stats4::coef(mle.res)[5]
    list(beta=beta,alpha=alpha)
  }

  df_mod$time[m+1] <- system.time(
    model.fit <- scobit_mle(y,x1,x2)
  )[[3]]

  df[["prob"]] <- scobit_fct(D_agent,D_artif,model.fit$beta,model.fit$alpha)
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



