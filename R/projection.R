#' @export projection
projection <- function(X, fit, nrank=NULL, kappa_u=1e-8){
  
  
  # type.projection=c("proj","1dim","mdim","logratio")
  type.projection <- fit$type.projection
  
  
  if(FALSE){
    X=Xtrain;
    
    data <- sim.LogNormal(n=500, p=100, r=2, snr=0.001, d=d, seed=1, verbose=TRUE); 
    Xtest <- sim.LogNormal.test(data$params)$X2
    
    fit1 <- aCPCA(X, nrank=2)
    fit2 <- CPCA(X, nrank=2)
    
    sum( abs( fit1$uhat[,1] - fit2$uhat[,1] ) )
    sum( abs( fit1$xhat[,1] - fit2$xhat[,1] ) )
    projection(X, fit1)$xhat[1:5,1:5]
    proj1 <- fit1 %>% { projection(Xtest, ., nrank=1, kappa_u=1e-16) }
    proj2 <- fit2 %>% { projection(Xtest, ., nrank=1, kappa_u=1e-16) }
    
    sum(abs(proj1$xhat - proj2$xhat ))
    
    rmspe.list <- try(sapply(fit.list, function(fit){
      xhat_test <- projection(Xtest, fit)
      # mean(rowSums((Xtest - xhat_test)^2))
      sqrt(mean((Xtest - xhat_test)^2))
    }))
  }
  
  n=nrow(X); p=ncol(X); 
  
  if(is.null(nrank)){
    r=NCOL(fit$uhat)
  } else {
    r=nrank
    if(type.projection == "logratio"){
      fit$logvhat <- fit$logvhat[,1:r,drop=F]
    } else {
      fit$vhat <- fit$vhat[,1:r,drop=F]
    }
  }
  
  
  mu <- fit$mu
  uhat <- fit$uhat
  vhat <- fit$vhat
  
  if( type.projection == "proj" ){
    
    xhat2 <- (X-tcrossprod(rep(1,n),fit$mu)) %*% tcrossprod(fit$vhat)
    xhat <- t(apply(xhat2,1,proj2simplex))
    
    return(list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat))
    
  } else if( type.projection == "1dim" ){
    
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
    
  } else if( type.projection == "mdim" ){
    
    uhat <- matrix(0,n,r)
    for( i in 1:n ){
      # uhat[i,] <- multidimconvexprojection(mu, as.vector(X[i,]), vhat)
      if( r > 1 ){
        uhat[i,] <- Solve_U_GP(as.vector(X[i,]), mu, vhat, gamma=0, kappa_u=kappa_u)
      } else if( r==1 ) {
        uhat[i,1] <- Solve_U_SP(as.vector(X[i,]), mu, vhat[,1], gamma=0)
      }
    }
    # it=fit$fit.path[[r]]$it;  gamma=fit$params$gamma
    # uhat <- uhat * (1-gamma/it)
    
  } else if( type.projection == "logratio" ){
    
    if( length(fit$delta) == 1 ){
      fit$delta <- rep(fit$delta, ncol(X))
    }
    f <- switch(fit$zero.replace, 
                "simple"=ZeroReplace.simple,
                "additive"=ZeroReplace.additive,
                "multiplicative"=ZeroReplace.multiplicative)
    
    Xnew <- t(apply(X, 1, f, delta=fit$delta))
    Xclr <- t(apply(Xnew, 1, clr))
    
    mu <- t(apply(fit$Xnew, 1, clr)) %>% colMeans()
    vhat <- fit$logvhat
    uhat <- (Xclr - tcrossprod(rep(1,n),mu)) %*% vhat
    
  }
  
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  if( type.projection == "logratio" ){
    xhat <- xhat %>% {t(apply(.,1,iclr))}
  }
  
  return(list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat))
}





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



