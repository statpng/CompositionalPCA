#' @export png.gppca_qp
png.gppca_qp <- function(X, nrank=2, maxit=500, eps=1e-8, kappa=1e-8, gamma=1e-1, phi=0.01, nu=1e-8, verbose=FALSE, V.init=NULL){
  if(FALSE){
    X=data$X2; nrank=r; V.init="PC"
    nrank=2; maxit=500; eps=1e-6; kappa=1e-4; gamma=0; verbose=TRUE; V.init=c("PC","random")[1]
  }
  
  if(FALSE){
    X <- sim.LogNormal(n=100, p=4, r=2, snr=2, d=2, seed=1, verbose=TRUE)$X2
    
    png.quaternary3d(X)
    
    fit <- png.gppca_qp(X, nrank=1, V.init="PC")
    png.quaternary3d(X, vhat=fit$vhat, xhat=fit$xhat)
  }
  
  if(FALSE){
    library(png.compositionalPCA)
    X <- R.utils::loadToEnv("./data/X.train.RData")$X.train
    X <- R.utils::loadToEnv("./data/X.RData")$X
    
    fit <- png.gppca_qp(X, nrank=5, kappa=1e-8, gamma=0.9, V.init="PC")
    
    fit$vhat
    mean(fit$X<1e-8)
    mean(fit$xhat<1e-4)
    fit %>% png.crit.path()
    
    
    fit$fit.path[[1]]$crit.path
    
    nrank=5; maxit=500; eps=1e-6; 
    kappa=1e-4; gamma=0; save.est.path=FALSE
    kappa.seq=10^(-seq(2,8,2))
    
    fit.ppca <- purrr::map(kappa.seq, ~try(png.ppca_qp(X, nrank=nrank, kappa=.x, gamma=0.5)))
    fit.gppca <- purrr::map(kappa.seq, ~try(png.gppca_qp(X, nrank=nrank, kappa=.x, gamma=0.2)))
    
    fit.pca <- png.pca(X,nrank)
    fit.ppca <- png.ppca_qp(X, nrank=nrank, kappa=1e-4)
    fit.gppca <- png.gppca_qp(X, nrank=nrank, kappa=1e-4)
    png.crit.path(fit.gppca[[6]])
    
    fit.ppca[[6]]$uhat %>% head
    fit.gppca$uhat %>% head
    
    
    fit.ppca[[6]] %>% png.crit.path
    fit.gppca$fit.path[[2]]$crit.path
    
    png.crit.path(fit.ppca[[1]])
    png.crit.path(fit.ppca[[2]])
    png.crit.path(fit.ppca[[3]])
    png.crit.path(fit.ppca[[4]])
  }
  
  if(FALSE){
    nrank=5; kappa=1e-8; gamma=0.2
  }
  if(FALSE){
    nrank=2; epsilon=1e-4; maxit=100; kappa=1e-8
  }
  
  
  
  
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
      # fit.rank1 <- png.rank1(X, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma)
      fit.rank1 <- png.rank12(X, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
      U_total <- cbind(U_total, fit.rank1$uhat)
      V_total <- cbind(V_total, fit.rank1$vhat)
      
      fit.path[[iter]] <- fit.rank1
    } else {
      # fit.UV <- UV_update(X=X, Vhat=V_total, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma)
      
      fit.UV <- UV_update2(X, U_total, V_total, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma, phi=phi, nu=1e-8, V.init=V.init, verbose=verbose)
      
      
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
              kappa=kappa, 
              gamma=gamma)
  
  return(list(mu=mu, uhat=U_total, vhat=V_total, 
              xhat=xhat, X=X, fit.path=fit.path, 
              fit.rank1=fit.rank1, maxit=maxit, 
              time=end-start, params=params, method="gppca"))
  
}












#' @export png.cv.cpca
png.cv.cpca <- function(X, nrank=2, nfold=5){
  if(FALSE){
    nrank <- 5
    nfold <- 5
    X <- pseq_list_total$urine$Phylum %>% .@otu_table %>% t
  }
  
  foldid <- sample( rep(1:nfold, length=nrow(X)) )
  
  fit.list <- NULL
  for( i in 1:nfold ){
    X.train <- X[foldid != i,]
    
    fit.pca <- png.pca(X, nrank=nrank)
    fit.lrpca <- png.lrpca( X.train, nrank=nrank, zero.replace="simple" )
    fit.ppca <- try(png.ppca_qp( X.train, nrank=nrank, kappa=1e-6, gamma=0 ))
    fit.gppca <- try(png.gppca_qp( X.train, nrank=nrank, kappa=1e-6, gamma=0 ))
    
    fit.list[[i]] <- list(fit.pca=fit.pca, 
                          fit.lrpca=fit.lrpca, 
                          fit.ppca=fit.ppca, 
                          fit.gppca=fit.gppca)
  }
  
  out <- NULL
  for( i in 1:nfold ){
    X.test <- X[foldid == i,]
    out[[i]] <- purrr::map(fit.list[[i]], ~{
      try( png.utils::png.rmse(X.test, png.projection(X.test, .x, method=.x$method)) )
    })
  }
  
  list(out=out, fit.list=fit.list, foldid=foldid)
}
