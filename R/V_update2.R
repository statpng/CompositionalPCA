#' @export V_update2
V_update2 <- function(X, Uhat, Vhat, Uk, kappa=1e-8){
  if(FALSE){
    X=X; Uhat=Unew[,1:(r-1)]; Vhat=Vhat; Uk=Unew[,r]; kappa=kappa
  }
  
  n=nrow(X); p=ncol(X); r=NCOL(Uhat)
  mu=colMeans(X)
  
  C <- tcrossprod(rep(1,n),mu) + tcrossprod(Uhat, Vhat)
  
  Vnew <- rep(0,p)
  Dmat <- diag(rep(1,p))
  dvec <- t(X-tcrossprod(rep(1,n),mu)) %*% Uk/sum(Uk^2)
  dvec <- t(X-C) %*% Uk/sum(Uk^2)
  
  solve_V <- function(Vhat, lbmat=NULL, ubmat=NULL, kappa=0) {
    p=NROW(Vhat);  r=NCOL(Vhat)
    lb <- if(!is.null(lbmat)) apply(lbmat,2,max) else NULL
    ub <- if(!is.null(ubmat)) apply(ubmat,2,min) else NULL
    
    bvec <- c(rep(0,r+1), lb-kappa, -(ub+kappa))
    Amat <- cbind(rep(1,p), Vhat, if (!is.null(lbmat)) diag(rep(1,p)) else NULL, if (!is.null(ubmat)) -diag(rep(1,p)) else NULL)
    Vnew <- try( solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=r+1, factorized=FALSE)$solution, silent=TRUE )
    if(inherits(Vnew, "try-error")){
      return( "No solution" )
    } else {
      Vnew <- Vnew/norm(Vnew,"2")
      Vnew <- ifelse(abs(Vnew)<1e-8, 0, Vnew)
      return( Vnew/norm(Vnew,"2") )
    }
    # bvec <- c(rep(0,1), rep(-1e-12,2*r), lb-kappa, -(ub+kappa))
    # Amat <- cbind(rep(1,p), Vhat, -Vhat, diag(rep(1,p)), -diag(rep(1,p)))
    # 
    # solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=1, factorized=FALSE)$solution
  }
  
  positiveindex = which(Uk>0); negativeindex = which(Uk<0)
  n.pos=length(positiveindex); n.neg=length(negativeindex)
  
  if(n.pos > 0 && n.neg > 0) {
    Vnew <- solve_V(Vhat=Vhat,
                    lbmat=-diag(1/Uk[positiveindex],n.pos,n.pos) %*% C[positiveindex,],
                    ubmat=-diag(1/Uk[negativeindex],n.neg,n.neg) %*% C[negativeindex,],
                    kappa=kappa)
  }
  
  if (n.pos>0 && n.neg==0) {
    Vnew <- solve_V(Vhat=Vhat,
                    lbmat=-diag(1/Uk[positiveindex],n.pos,n.pos) %*% C[positiveindex,], 
                    kappa=kappa)
  }
  
  if (n.pos==0 && n.neg>0) {
    Vnew <- solve_V(Vhat=Vhat,
                    ubmat=-diag(1/Uk[negativeindex],n.neg,n.neg) %*% C[negativeindex,], 
                    kappa=kappa)
  }
  
  return(Vnew)
}













#' @export Solve_U_SP
Solve_U_SP <- function(x, mu, v, gamma=0){
  if(FALSE){
    x=X[i,]; mu=C[i,]; v=Vold; gamma=gamma/(it^(1/2))
  }
  
  t=sum((x-mu)*v)
  
  indexminus=which(v<0);  indexplus=which(v>0)
  lminus=length(indexminus);  lplus=length(indexplus)
  
  if (lminus>0 & lplus>0){
    m <- max( -(mu/v)[indexplus]);  M=min( -(mu/v)[indexminus])
    U <- min(max(m,t),M)
  }
  
  if (lminus>0 & lplus==0){
    m <- max( -(mu/v)[indexminus])
    U <- max(m,t)
  }
  
  if (lminus==0 & lplus>0){
    M <- min( -(mu/v)[indexplus])
    U <- min(t,M)
  }
  
  if (lminus==0 & lplus==0){
    U <- t
  }
  
  return(U * (1-gamma))
}

#' @export Solve_U_GP
Solve_U_GP <- function(x, mu, V, gamma=0, nu=1e-8){
  if(FALSE){
    x=X[i,]; mu=C[i,]; V=cbind(Vhat, Vk_old); gamma=gamma/(it^(1/2))
  }
  
  x=as.numeric(x)
  mu=as.numeric(mu)
  
  U <- solve.QP(Dmat=crossprod(V), dvec=t(V) %*% (x - mu), Amat=t(V), bvec=-mu - nu)$solution
  # U <- solve.QP(Dmat=crossprod(V), dvec=t(V) %*% (x - mu), Amat=t(V), bvec=-mu)$solution
  return(U * (1-gamma))
}

