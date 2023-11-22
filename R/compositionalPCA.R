
#' @export PCA
PCA <- function(X, nrank=2){
  d0 <- sum(svd(X)$d>1e-10)
  # if( nrank > (d0-1) ) stop("nrank should be smaller than or equal to rank(X)-1")
  
  n=nrow(X); p=ncol(X); mu=colMeans(X)
  uhat <- prcomp(X)$x[,1:nrank,drop=F]
  vhat <- prcomp(X)$rotation[,1:nrank,drop=F]
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat, X=X, type.projection="proj") )
}




#' @export lrPCA
lrPCA <- function(X, nrank=2, zero.replace=NULL, delta="min0", delta.multiple=FALSE){
  
  if(grepl("min",delta)){
    delta_num <- as.numeric(gsub("min","",delta))
    
    if( delta.multiple ){
      delta <- apply(X,2,function(x) ifelse(all(x>0),1,min(x[x>0]) * 10^(delta_num)))
    } else {
      delta <- rep( min(X[X>0]) * 10^(delta_num), ncol(X) )
    }
    
  }
  
  
  n=nrow(X); p=ncol(X);
  
  Xnew <- X
  if(!is.null(zero.replace)){
    for(i in 1:nrow(X)){
      x <- X[i,]
      
      f <- switch(zero.replace,
                  "simple"=ZeroReplace.simple,
                  "additive"=ZeroReplace.additive,
                  "multiplicative"=ZeroReplace.multiplicative)
      
      Xnew[i,] <- f(x, delta)
    }
  }
  
  
  mu = iclr( colMeans( log(Xnew) ) )
  
  Xclr <- t(apply(Xnew,1,clr))
  fit <- PCA(Xclr, nrank=nrank)
  uhat <- fit$uhat
  vhat <- lapply( seq(-20,20,0.5), function(z){
    z*fit$vhat
  })
  xhat <- fit$xhat %>% {t(apply(.,1,iclr))}
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, logmu=colMeans( log(Xnew) ), logvhat=fit$vhat, xhat=xhat, X=X, Xnew=Xnew, type.projection="logratio", zero.replace=zero.replace, delta=delta) )
}








# Deprecated
pPCA <- function(X, nrank=2){
  
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
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat, X=X, type.projection="1dim") )
}





#' @export crPCA
crPCA <- function(X, nrank=2, V=prcomp(X)$rotation[,1:nrank,drop=F]){
  
  if(FALSE){
    set.seed(3)
    n=500; p=4; r=2
    X <- sim.simplex2(n=n, p=p, r=r, delta=runif(1,0.8,1.2))$X2
    ppca(X,2)$xhat %>% {sum(.<(-1e-10))}
    crPCA(X,2)$xhat %>% {sum(.<(-1e-10))}
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
  
  return( list(mu=mu, uhat=uhat, vhat=V, xhat=xhat, X=X, type.projection="mdim") )
}










#' @export aCPCA
aCPCA <- function(X, nrank=2, maxit=500, eps=1e-8, kappa_v=1e-8, gamma=1e-1, phi=0.01, V.init=c("PC","random"), verbose=FALSE){
  
  if(length(V.init)>1) V.init <- "PC"
  
  
  if(FALSE){
    X; nrank=2; maxit=500; eps=1e-8; kappa_v=1e-8; gamma=1e-2; phi=0.01; V.init=c("PC","random"); verbose=FALSE
  }
  
  library(quadprog)
  
  if(is.null(V.init)){
    V.init <- c("PC","random")[1]
  }
  
  n=nrow(X); p=ncol(X); r=nrank
  mu=colMeans(X);
  
  
  Uhat <- Vhat <- NULL
  fit.path <- NULL
  for( iter in 1:nrank ){
    if( iter == 1 ){
      fit.rank1 <- rank12(X, maxit=maxit, eps=eps, kappa_v=kappa_v, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
      Uhat <- cbind(Uhat, fit.rank1$uhat)
      Vhat <- cbind(Vhat, fit.rank1$vhat)
      
      fit.path[[iter]] <- fit.rank1
    } else {
      
      fit.UkVk <- update_UkVk(X=X, Uhat=Uhat, Vhat=Vhat, maxit=maxit, eps=eps, kappa_v=kappa_v, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
      Uhat <- cbind(Uhat, fit.UkVk$uhat)
      Vhat <- cbind(Vhat, fit.UkVk$vhat)
      
      # if(!save.est.path){
      #   fit.UkVk$est.path <- NULL
      # }
      
      fit.path[[iter]] <- fit.UkVk
    }
  }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Uhat, Vhat)
  
  params=list(nrank=nrank, 
              maxit=maxit, 
              eps=eps, 
              kappa_v=kappa_v, 
              gamma=gamma)
  
  return(list(mu=mu, uhat=Uhat, vhat=Vhat, xhat=xhat, 
              X=X, fit.path=fit.path, maxit=maxit, type.projection="1dim", params=params))
  
}











#' @export CPCA
CPCA <- function(X, nrank=2, maxit=500, eps=1e-8, kappa_v=1e-8, gamma=1e-1, phi=0.01, kappa_u=1e-8, V.init=c("PC","random"), verbose=FALSE ){
  
  if(length(V.init)>1) V.init <- "PC"
  
  if(is.null(V.init)){
    V.init <- c("PC","random")[1]
  }
  
  
  n=nrow(X); p=ncol(X); r=nrank
  mu=colMeans(X);
  
  
  if( p < r ) stop("r should be less than or equal to p !")
  
  
  start <- proc.time()
  
  U_total <- V_total <- NULL
  fit.path <- NULL
  for( iter in 1:nrank ){
    if( iter == 1 ){
      # fit.rank1 <- rank1(X, maxit=maxit, eps=eps, kappa_v=kappa_v, gamma=gamma)
      fit.rank1 <- rank12(X, maxit=maxit, eps=eps, kappa_v=kappa_v, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
      U_total <- cbind(U_total, fit.rank1$uhat)
      V_total <- cbind(V_total, fit.rank1$vhat)
      
      fit.path[[iter]] <- fit.rank1
    } else {
      # fit.UV <- UV_update(X=X, Vhat=V_total, maxit=maxit, eps=eps, kappa_v=kappa_v, gamma=gamma)
      
      fit.UV <- UV_update2(X, U_total, V_total, maxit=maxit, eps=eps, kappa_v=kappa_v, gamma=gamma, phi=phi, kappa_u=1e-8, V.init=V.init, verbose=verbose)
      
      
      U_total <- fit.UV$uhat
      V_total <- fit.UV$vhat
      
      # if(!save.est.path){
      #   fit.UV$est.path <- NULL
      # }
      
      fit.path[[iter]] <- fit.UV
    }
  }
  
  end <- proc.time()
  
  colnames(U_total) <- colnames(V_total) <- NULL
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(U_total, V_total)
  xhat <- ifelse(abs(xhat)<1e-12,0,xhat)
  
  
  params=list(nrank=nrank, 
              maxit=maxit, 
              eps=eps, 
              kappa_v=kappa_v, 
              gamma=gamma)
  
  return(list(mu=mu, uhat=U_total, vhat=V_total, 
              xhat=xhat, X=X, fit.path=fit.path, 
              fit.rank1=fit.rank1, maxit=maxit, 
              time=end-start, params=params, type.projection="mdim"))
  
}












#' @export cv.CPCA
cv.CPCA <- function(X, nrank=2, nfold=5){
  if(FALSE){
    nrank <- 5
    nfold <- 5
    X <- pseq_list_total$urine$Phylum %>% .@otu_table %>% t
  }
  
  foldid <- sample( rep(1:nfold, length=nrow(X)) )
  
  fit.list <- NULL
  for( i in 1:nfold ){
    X.train <- X[foldid != i,]
    
    fit.PCA <- PCA(X, nrank=nrank)
    fit.lrPCA <- lrPCA( X.train, nrank=nrank, zero.replace="simple" )
    fit.aCPCA <- try(aCPCA( X.train, nrank=nrank, kappa_v=1e-6, gamma=0 ))
    fit.CPCA <- try(CPCA( X.train, nrank=nrank, kappa_v=1e-6, gamma=0 ))
    
    fit.list[[i]] <- list(fit.PCA=fit.PCA, 
                          fit.lrPCA=fit.lrPCA, 
                          fit.aCPCA=fit.aCPCA, 
                          fit.CPCA=fit.CPCA)
  }
  
  out <- NULL
  for( i in 1:nfold ){
    X.test <- X[foldid == i,]
    out[[i]] <- purrr::map(fit.list[[i]], ~{
      try( utils::rmse(X.test, projection(X.test, .x, type.projection=.x$type.projection)) )
    })
  }
  
  list(out=out, fit.list=fit.list, foldid=foldid)
}

















#' @export fit_all
fit_all <- function(X, nrank, ...){
  
  X <- matrix( as.numeric(X), nrow(X), ncol(X) )
  delta.seq <- paste0("min", c(-2,-1,0,1,2))
  
  fit1 <- purrr::map(delta.seq, ~ lrPCA(X, nrank=nrank, zero.replace="simple", delta=.x))
  names(fit1) <- paste0("lrPCA (", gsub("min","",delta.seq), ")")
  
  fit2 <- crPCA(X, nrank=nrank)
  
  gamma.seq <- 0.1 # 10^(-seq(1, 5, 2))
  fit3 <- purrr::map(gamma.seq, ~try(aCPCA(X, nrank=nrank, gamma=.x, ...)))
  fit4 <- purrr::map(gamma.seq, ~try(CPCA(X, nrank=nrank, gamma=.x, ...)))
  
  names(fit3) <- paste0("aCPCA (", format(gamma.seq,digits=2), ")" )
  names(fit4) <- paste0("CPCA (", format(gamma.seq,digits=2), ")" )
  
  fit.list <- list(fit1,list(crPCA=fit2),fit3,fit4) %>% Reduce(append, .)
  fit.list
}





#' @export fit_all.cv
fit_all.cv <- function(X, nfold=5, nrank=5, ...){
  set.seed(123);  # nfold=nrow(X)
  foldid <- sample( rep(1:nfold, length=nrow(X)) )
  
  fit.list <- NULL
  for( i in 1:nfold ){
    cat(i, " / ", nfold, "\n")
    X.train <- X[foldid != i,]
    if(FALSE) X=X.train; nrank=nrank
    fit.list[[i]] <- try(fit_all(X=X.train, nrank=nrank, ...))
  }
  return(list(fit.list=fit.list, foldid=foldid))
}







