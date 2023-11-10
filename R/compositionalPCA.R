#' @export onedimconvexprojection
onedimconvexprojection<-function(c,x,v){
  # one-dimensional convex projection function
  
  # c,x is an element of S^{p} and v is a unit vector perpendicular with 1=(1,1,...,1)\in \mathbb{R}^{p}
  
  # direct projection score
  t=sum((x-c)*v)
  
  # define index set
  indexminus<-which(v<0)
  indexplus<-which(v>0)
  
  lminus<-length(indexminus)
  lplus<-length(indexplus)
  
  # calculate convec projection score case by case
  if (lminus>0 & lplus>0){
    m=max( -(c/v)[indexplus])
    M=min( -(c/v)[indexminus])
    return(min(max(m,t),M))
  }
  
  if (lminus>0 & lplus==0){
    m=max( -(c/v)[indexplus])
    return(max(m,t))
  }
  
  if (lminus==0 & lplus>0){
    M=min( -(c/v)[indexminus])
    return(min(t,M))
  }
  
  if (lminus==0 & lplus==0){
    return(t)
  }
  
}

# # example of onedim convex projection
# {
#   p<-100
#   c<-rep(1,p)/p
#   x<-runif(p,0,1)
#   x<-x/sum(x)
#   v<-rnorm(p,0,1)
#   v<-v-sum(v*1)/p
#   v<-v/sqrt(sum(v^2))
#   
#   t<-onedimconvexprojection(c,x,v)
# }




#' @export multidimconvexprojection
multidimconvexprojection<-function(c,x,V){
  # c,x is an element of S^{p} and V=(V1, V2, ..., Vr) is r number of orthonormal vectors perpendicular with 1=(1,1,...,1)\in \mathbb{R}^{p}, given by V\in \mathbb{R}^{p \times r}
  
  if(FALSE){
    set.seed(1)
    n=1000; p=3; r=2
    X <- sim.simplex2(n=n, p=p, r=r)$X2
    
    
    mu <- colMeans(X)
    vhat <- prcomp(X)$rotation[,1:(ncol(X)-1)]
    x <- X[1,]
    
    mu + tcrossprod( multidimconvexprojection(mu,x,vhat[,1:2]), vhat[,1:2] )
    
    c=C[i,]; x=X[i,]; V=Vold
  }
  
  
  require(quadprog)
  
  r<-ncol(V)
  U<-solve.QP(Dmat=diag(rep(1,r)), dvec=t(V) %*% (x-c), Amat=t(V), bvec=-c )$solution
  return(U)
  
}











#' @export png.pca
png.pca <- function(X, nrank=2){
  d0 <- sum(svd(X)$d>1e-10)
  # if( nrank > (d0-1) ) stop("nrank should be smaller than or equal to rank(X)-1")
  
  n=nrow(X); p=ncol(X); mu=colMeans(X)
  uhat <- prcomp(X)$x[,1:nrank,drop=F]
  vhat <- prcomp(X)$rotation[,1:nrank,drop=F]
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat, X=X, method="pca") )
}




#' @export png.ppca
png.ppca <- function(X, nrank=2){
  
  if(FALSE){
    set.seed(2)
    n=500; p=4; r=2
    X <- sim.simplex2(n=n, p=p, r=r, snr=1, d=1, d0=10)$X2
    png.quaternary3d(X,vhat=png.ppca(X,2)$vhat, xhat=png.ppca(X,2)$xhat)
    
    png.ppca(X,2)$xhat %>% {sum(.<0)}
    
    png.ppca(X,2)$xhat[png.ppca(X,2)$xhat<0]
  }
  
  
  d0 <- sum(svd(X)$d>1e-10)
  
  if( nrank > (d0-1) ) stop("nrank should be smaller than rank(X)-1")
  
  
  n=nrow(X); p=ncol(X)
  mu <- colMeans(X)
  vhat <- prcomp(X)$rotation[,1:nrank,drop=F]
  
  
  uhat <- NULL
  for( j in 1:nrank ){
    if(j==1) Xnew=X
    v <- vhat[,j,drop=F]
    uhat_j <- apply(Xnew,1,function(x) onedimconvexprojection(mu, x, v))
    uhat <- cbind(uhat, uhat_j)
    Xnew <- Xnew #- tcrossprod(uhat,vhat[,1:j])
  }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(uhat,vhat)
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat, X=X, method="ppca") )
}




#' @export png.gppca
png.gppca <- function(X, nrank=2, V=prcomp(X)$rotation[,1:nrank,drop=F]){
  
  if(FALSE){
    set.seed(3)
    n=500; p=4; r=2
    X <- sim.simplex2(n=n, p=p, r=r, delta=runif(1,0.8,1.2))$X2
    png.ppca(X,2)$xhat %>% {sum(.<(-1e-10))}
    png.gppca(X,2)$xhat %>% {sum(.<(-1e-10))}
    nrank=1
  }
  
  
  d0 <- sum(svd(X)$d>1e-10)
  
  if( nrank > (d0-1) ) stop("nrank should be smaller than rank(X)-1")
  
  
  mu <- colMeans(X)
  n <- nrow(X)
  
  uhat <- apply(X,1,function(x){
    multidimconvexprojection(mu,x,V)
  }) %>% {matrix(., n, nrank, byrow=T)}
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(uhat, V)
  
  return( list(mu=mu, uhat=uhat, vhat=V, xhat=xhat, X=X, method="gppca") )
}




#' @export png.fit_all
png.fit_all <- function(X, nrank, ...){
  # delta1 <- 1e-12
  # delta2 <- 1e-8
  # delta.seq <- c(1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
  X <- matrix( as.numeric(X), nrow(X), ncol(X) )
  delta.seq <- paste0("min", c(-2,-1,0,1,2))
  
  fit1 <- purrr::map(delta.seq, ~ png.lrpca(X, nrank=nrank, zero.replace="simple", delta=.x))
  names(fit1) <- paste0("lrPCA (", gsub("min","",delta.seq), ")")
  
  fit2 <- png.gppca(X, nrank=nrank)
  
  gamma.seq <- 0.1 # 10^(-seq(1, 5, 2))
  fit3 <- purrr::map(gamma.seq, ~try(png.ppca_qp(X, nrank=nrank, gamma=.x, ...)))
  fit4 <- purrr::map(gamma.seq, ~try(png.gppca_qp(X, nrank=nrank, gamma=.x, ...)))
  
  names(fit3) <- paste0("aCPCA (", format(gamma.seq,digits=2), ")" )
  names(fit4) <- paste0("CPCA (", format(gamma.seq,digits=2), ")" )
  
  # fit5 %>% sapply(function(fit){
  #   sqrt(mean((data$X0 - fit$xhat)^2))
  # })
  # fit6 %>% sapply(function(fit){
  #   sqrt(mean((data$X0 - fit$xhat)^2))
  # })
  # fit6[[3]] %>% png.crit.path()
  
  fit.list <- list(fit1,list(crPCA=fit2),fit3,fit4) %>% Reduce(append, .)
  fit.list
}







#' @export png.pca.rmse.linear
png.pca.rmse.linear <- function(fit.list, data, n.test="10x"){
  X0=data$X0
  params=data$params
  
  if(n.test == "10x"){
    params_test <- png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.simplex.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    # mean(rowSums((X0 - fit$xhat)^2))
    sqrt(mean((X0 - fit$xhat)^2))
  }))
  rmspe.list <- try(sapply(fit.list, function(fit){
    xhat_test <- png.projection(Xtest, fit, method=fit$method)$xhat
    # mean(rowSums((Xtest - xhat_test)^2))
    sqrt(mean((Xtest - xhat_test)^2))
  }))
  
  cbind(rmse=rmse.list, rmspe=rmspe.list)
}




#' @export png.pca.rmse.LogNormal
png.pca.rmse.LogNormal <- function(fit.list, data, n.test="10x"){
  X0=data$X0
  params=data$params
  
  if(n.test == "10x"){
    params_test <- png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.LogNormal.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    # mean(rowSums((X0 - fit$xhat)^2))
    sqrt(mean((X0 - fit$xhat)^2))
  }))
  rmspe.list <- try(sapply(fit.list, function(fit){
    xhat_test <- png.projection(Xtest, fit, method=fit$method)$xhat
    # mean(rowSums((Xtest - xhat_test)^2))
    sqrt(mean((Xtest - xhat_test)^2))
  }))
  
  cbind(rmse=rmse.list, rmspe=rmspe.list)
}



#' @export png.pca.rmse.LogNormal_rank
png.pca.rmse.LogNormal_rank <- function(fit.list, data, n.test="10x"){
  X0=data$X0
  params=data$params
  
  if(n.test == "10x"){
    params_test <- png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.LogNormal.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- png.projection(data$X2, fit=fit, nrank=r, method=fit$method)
      sqrt(mean((data$X0 - Xhat)^2))
    }) %>% unlist
    
  }))
  
  
  # rmse.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(method, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=method)) + png.utils::png.ggplot.scale_y_log10()
  
  
  rmspe.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- png.projection(Xtest, fit=fit, nrank=r, method=fit$method)
      sqrt(mean((Xtest - Xhat)^2))
    }) %>% unlist
    
  }))
  
  # rmspe.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(method, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=method)) + png.utils::png.ggplot.scale_y_log10()
  
  
  list(rmse=rmse.list, rmspe=rmspe.list)
}







#' @export png.pca.rmse.linear_rank
png.pca.rmse.linear_rank <- function(fit.list, data, n.test="10x"){
  X0=data$X0
  params=data$params
  
  if(n.test == "10x"){
    params_test <- png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.simplex.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- png.projection(data$X2, fit=fit, nrank=r, method=fit$method)
      sqrt(mean((data$X0 - Xhat)^2))
    }) %>% unlist
    
  }))
  
  
  # rmse.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(method, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=method)) + png.utils::png.ggplot.scale_y_log10()
  
  
  rmspe.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- png.projection(Xtest, fit=fit, nrank=r, method=fit$method)
      sqrt(mean((Xtest - Xhat)^2))
    }) %>% unlist
    
  }))
  
  # rmspe.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(method, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=method)) + png.utils::png.ggplot.scale_y_log10()
  
  
  list(rmse=rmse.list, rmspe=rmspe.list)
}







#' @export png.fit_all.cv
png.fit_all.cv <- function(X, nfold=5, nrank=5, ...){
  set.seed(123);  # nfold=nrow(X)
  foldid <- sample( rep(1:nfold, length=nrow(X)) )
  
  fit.list <- NULL
  for( i in 1:nfold ){
    cat(i, " / ", nfold, "\n")
    X.train <- X[foldid != i,]
    if(FALSE) X=X.train; nrank=nrank
    fit.list[[i]] <- try(png.fit_all(X=X.train, nrank=nrank, ...))
  }
  return(list(fit.list=fit.list, foldid=foldid))
}







png.adist <- function(x,y){
  sqrt(sum((png.clr(x) - png.clr(x))^2))
}

png.STRESS <- function(X){
  if(FALSE){
    set.seed(1); 
    X=rdirichlet(10,c(1,2,3)) %>% {ifelse(.<0.2,0,.) %>% apply(1,function(x)x/sum(x)) %>% t} 
    Y=rdirichlet(10,c(1,2,3)) %>% {ifelse(.<0.2,0,.) %>% apply(1,function(x)x/sum(x)) %>% t} 
  }
  
  MIN <-  min(X[X>0])
  delta.seq <- MIN * pmin(1,10^seq(-2,2,1))
  
  delta <- delta.seq[1]
  
  Xnew <- png.ZeroReplace.simple(X, delta=delta)
  
  for( i in 1:nrow(X) ){
    png.adist( X[i,], Xnew[i,] )
  }
  
  
  sqrt(sum((png.clr(x1) - png.clr(x2))^2))
}


