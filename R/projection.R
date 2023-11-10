#' @export png.projection_old
png.projection_old <- function(X, fit, nrank=NULL, method=c("pca","ppca", "gppca", "lrpca")){
  if(FALSE){
    X=Xtrain; method=fit$method
  }
  
  n=nrow(X); p=ncol(X); 
  
  if(is.null(nrank)){
    r=NCOL(fit$uhat)
  } else {
    r=nrank
    if(method == "lrpca"){
      fit$logvhat <- fit$logvhat[,1:r,drop=F]
    } else {
      fit$vhat <- fit$vhat[,1:r,drop=F]
    }
  }
  
  
  mu=fit$mu
  vhat=fit$vhat
  
  if( method == "pca" ){
    
    xhat2 <- (X-tcrossprod(rep(1,n),fit$mu)) %*% tcrossprod(fit$vhat)
    xhat <- t(apply(xhat2,1,png.proj2simplex))
    
    return(xhat)
    
  } else if( method == "ppca" ){
    
    uhat <- matrix(0,n,r)
    for( k in 1:r ){
      if(k-1 == 0){
        chat <- tcrossprod(rep(1,n), mu)
      } else {
        chat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat[,1:(k-1)], vhat[,1:(k-1)])
      }
      
      for( i in 1:n ){
        # uhat[i,k] <- onedimconvexprojection(chat[i,], as.vector(X[i,]), vhat[,k])
        uhat[i,k] <- Solve_U_SP(as.vector(X[i,]), chat[i,], vhat[,k], gamma=0)
      }
      # it=fit$fit.path[[k]]$it;  gamma=fit$params$gamma
      # uhat[,k] <- uhat[,k] * (1-gamma/it)
    }
    
  } else if( method == "gppca" ){
    
    uhat <- matrix(0,n,r)
    for( i in 1:n ){
      # uhat[i,] <- multidimconvexprojection(mu, as.vector(X[i,]), vhat)
      uhat[i,] <- Solve_U_GP(as.vector(X[i,]), mu, vhat, gamma=0)
    }
    # it=fit$fit.path[[r]]$it;  gamma=fit$params$gamma
    # uhat <- uhat * (1-gamma/it)
    
  } else if( method == "lrpca" ){
    
    f <- switch(fit$zero.replace, 
                "simple"=png.ZeroReplace.simple,
                "additive"=png.ZeroReplace.additive,
                "multiplicative"=png.ZeroReplace.multiplicative)
    
    Xnew <- t(apply(X, 1, f, delta=fit$delta))
    Xclr <- t(apply(Xnew, 1, png.clr))
    
    mu <- t(apply(fit$Xnew, 1, png.clr)) %>% colMeans()
    vhat <- fit$logvhat
    uhat <- (Xclr - tcrossprod(rep(1,n),mu)) %*% vhat
    
  }
  
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  if( method == "lrpca" ){
    xhat <- xhat %>% {t(apply(.,1,png.iclr))}
  }
  
  return(xhat)
}













#' @export png.projection
png.projection <- function(X, fit, nrank=NULL, method=c("pca_proj","ppca", "gppca", "lrpca"), nu=1e-8){
  if(FALSE){
    X=Xtrain; method=fit$method
    
    data <- sim.LogNormal(n=500, p=100, r=2, snr=0.001, d=d, seed=1, verbose=TRUE); 
    Xtest <- sim.LogNormal.test(data$params)$X2
    
    fit1 <- png.ppca_qp(X, nrank=2)
    fit2 <- png.gppca_qp(X, nrank=2)
    
    sum( abs( fit1$uhat[,1] - fit2$uhat[,1] ) )
    sum( abs( fit1$xhat[,1] - fit2$xhat[,1] ) )
    png.projection(X, fit1, method=fit1$method)$xhat[1:5,1:5]
    proj1 <- fit1 %>% { png.projection3(Xtest, ., nrank=1, method=.$method, nu=1e-16) }
    proj2 <- fit2 %>% { png.projection3(Xtest, ., nrank=1, method=.$method, nu=1e-16) }
    
    sum(abs(proj1$xhat - proj2$xhat ))
    
    rmspe.list <- try(sapply(fit.list, function(fit){
      xhat_test <- png.projection(Xtest, fit, method=fit$method)
      # mean(rowSums((Xtest - xhat_test)^2))
      sqrt(mean((Xtest - xhat_test)^2))
    }))
  }
  
  n=nrow(X); p=ncol(X); 
  
  if(is.null(nrank)){
    r=NCOL(fit$uhat)
  } else {
    r=nrank
    if(method == "lrpca"){
      fit$logvhat <- fit$logvhat[,1:r,drop=F]
    } else {
      fit$vhat <- fit$vhat[,1:r,drop=F]
    }
  }
  
  
  mu=fit$mu
  vhat=fit$vhat
  
  if( method == "pca_proj" ){
    
    xhat2 <- (X-tcrossprod(rep(1,n),fit$mu)) %*% tcrossprod(fit$vhat)
    xhat <- t(apply(xhat2,1,png.proj2simplex))
    
    return(list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat))
    
  } else if( method == "ppca" ){
    
    uhat <- matrix(0,n,r)
    for( k in 1:r ){
      if(k-1 == 0){
        chat <- tcrossprod(rep(1,n), mu)
      } else {
        chat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat[,1:(k-1)], vhat[,1:(k-1)])
      }
      
      for( i in 1:n ){
        # uhat[i,k] <- onedimconvexprojection(chat[i,], as.vector(X[i,]), vhat[,k])
        uhat[i,k] <- Solve_U_SP(as.vector(X[i,]), chat[i,], vhat[,k], gamma=0)
      }
      # it=fit$fit.path[[k]]$it;  gamma=fit$params$gamma
      # uhat[,k] <- uhat[,k] * (1-gamma/it)
    }
    
  } else if( method == "gppca" ){
    
    uhat <- matrix(0,n,r)
    for( i in 1:n ){
      # uhat[i,] <- multidimconvexprojection(mu, as.vector(X[i,]), vhat)
      if( r > 1 ){
        uhat[i,] <- Solve_U_GP(as.vector(X[i,]), mu, vhat, gamma=0, nu=nu)
      } else if( r==1 ) {
        uhat[i,1] <- Solve_U_SP(as.vector(X[i,]), mu, vhat[,1], gamma=0)
      }
    }
    # it=fit$fit.path[[r]]$it;  gamma=fit$params$gamma
    # uhat <- uhat * (1-gamma/it)
    
  } else if( method == "lrpca" ){
    
    if( length(fit$delta) == 1 ){
      fit$delta <- rep(fit$delta, ncol(X))
    }
    f <- switch(fit$zero.replace, 
                "simple"=png.ZeroReplace.simple,
                "additive"=png.ZeroReplace.additive,
                "multiplicative"=png.ZeroReplace.multiplicative)
    
    Xnew <- t(apply(X, 1, f, delta=fit$delta))
    Xclr <- t(apply(Xnew, 1, png.clr))
    
    mu <- t(apply(fit$Xnew, 1, png.clr)) %>% colMeans()
    vhat <- fit$logvhat
    uhat <- (Xclr - tcrossprod(rep(1,n),mu)) %*% vhat
    
  }
  
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  if( method == "lrpca" ){
    xhat <- xhat %>% {t(apply(.,1,png.iclr))}
  }
  
  return(list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat))
}
