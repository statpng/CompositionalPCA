
#' @export update_UkVk
update_UkVk <- function(X, Uhat, Vhat, maxit=500, eps=1e-8, kappa_v=1e-4, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  
  require(quadprog)
  
  n=nrow(X); p=ncol(X); d=NCOL(Uhat)
  mu=colMeans(X)
  
  C <- tcrossprod(rep(1,n), mu) + tcrossprod(Uhat, Vhat)
  
  # Initial V0
  V0.PC <- prcomp(X-(tcrossprod(rep(1,n), mu) + X %*% tcrossprod(Vhat)))$rotation[,1]
  if(V.init == "PC"){
    V0 <- V0.PC
  } else {
    V0 <- (V0.PC + qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]*phi) %>% {./norm(.,"2")}
  }
  
  # angle(V0.PC, V0)[[1]] %>% print
  
  
  # est.path <- 
  crit.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vold <- V0
    } else {
      Vold <- Vnew
    }
    
    Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vold, gamma=gamma/(it^(1/2))))
    Vnew <- V_update2(X, Uhat=Uhat, Vhat=as.matrix(Vhat), Uk=as.matrix(Unew), kappa_v=kappa_v)
    if( Vnew[1] == "No solution" ){
      Vnew <- Vold
      NoSolution <- TRUE
    } else {
      NoSolution <- FALSE
    }
    
    # est.path[[it]] <- list(uhat=cbind(Uhat, Unew), vhat=cbind(Vhat, Vnew))
    crit.path[it] <- min( sum((Vnew-Vold)^2), sum((Vnew+Vold)^2) )
    
    if( crit.path[it] < eps ) break
    
  }
  
  if(NoSolution) crit.path[it] <- 0.1
  
  
  Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vnew, gamma=0))
  
  xhat <- C + tcrossprod(Unew,Vnew)
  
  loss <- mean(rowSums((X - xhat)^2))
  
  # result <- list(uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path)
  result <- list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss)
  
  return(result)
  
}













#' @export rank12
rank12 <- function(X, maxit=500, eps=1e-8, kappa_v=1e-8, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  
  if(FALSE){
    data <- sim.Linear(n=500, p=100, r=2, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE); 
    X <- data$X2
  }
  require(quadprog)
  
  n=nrow(X); p=ncol(X);
  mu=colMeans(X); 
  C=tcrossprod(rep(1,n), mu)
  
  # Initial V0
  V0.PC <- prcomp(X-C)$rotation[,1]
  if(V.init == "PC"){
    V0 <- V0.PC
  } else {
    V0 <- (V0.PC + qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]*phi) %>% {./norm(.,"2")}
  }
  
  # angle(V0.PC, V0)[[1]] %>% print
  
  
  est.path <- crit.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vold <- V0
    } else {
      Vold <- Vnew
    }
    
    
    Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vold, gamma=gamma/(it^(1/2))))
    
    Vnew <- V_update2(X, Uhat=matrix(0,n,1), Vhat=matrix(0,p,1), Uk=as.matrix(Unew), kappa_v=kappa_v)
    
    crit <- min( sum((Vnew-Vold)^2), sum((Vnew+Vold)^2) )
    
    crit.path[it] <- crit
    est.path[[it]] <- list(uhat=Unew, vhat=Vnew)
    
    if( crit < eps ) break
  }
  
  
  Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vnew, gamma=0))
  
  # Unew1 <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vnew, gamma=0))
  # Unew2 <- sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=Vnew, gamma=0))
  # 
  # Unew1[ ( abs(Unew1 - Unew2) > 1e-5 ) ]
  # 
  # sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vnew, gamma=0)) %>% head
  # sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=Vnew, gamma=0)) %>% head
  
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew,Vnew)
  
  loss <- mean(rowSums((X - xhat)^2))
  
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path) )
  # return( list(# xhat=xhat, 
  #   mu=mu, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss) )
  
}











#' @export UV_update2
UV_update2 <- function(X, Uhat, Vhat, maxit=500, eps=1e-6, kappa_v=1e-4, gamma=0, kappa_u=1e-8, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  if(FALSE){
    Uhat=U_total; Vhat=V_total
  }
  
  
  n=nrow(X); p=ncol(X); r=NCOL(Vhat) + 1
  mu=colMeans(X)
  
  C <- tcrossprod(rep(1,n), mu)
  # C <- tcrossprod(rep(1,n), mu) + X %*% tcrossprod(Vhat)
  
  
  # Initial V0
  V0.PC <- prcomp(X-(tcrossprod(rep(1,n), mu) + X %*% tcrossprod(Vhat)))$rotation[,1]
  if(V.init == "PC"){
    V0 <- V0.PC
  } else {
    V0 <- (V0.PC + qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]*phi) %>% {./norm(.,"2")}
  }
  
  # angle(V0.PC, V0)[[1]] %>% print
  
  crit.path <- est.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vk_old <- V0
    } else {
      Vk_old <- Vk_new
    }
    
    
    
    # U-update
    Unew <- t(sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=cbind(Vhat, Vk_old), gamma=gamma/(it^(1/2)), kappa_u=kappa_u)))
    
    # V-update
    Vk_new <- V_update2(X, Uhat=Unew[,1:(r-1)], Vhat=Vhat, Uk=Unew[,r], kappa_v=kappa_v)
    if( Vk_new[1] == "No solution" ){
      Vk_new <- Vk_old
      NoSolution <- TRUE
    } else {
      NoSolution <- FALSE
    }
    
    crit <- min( (sum((Vk_old-Vk_new)^2)), (sum((Vk_old+Vk_new)^2)) )
    
    crit.path[it] <- crit
    est.path[[it]] <- list(uhat=Unew, vhat=cbind(Vhat,Vk_new))
    
    if( crit < eps ) break
  }
  
  if(NoSolution) crit.path[it] <- 0.1
  
  Vnew <- cbind(Vhat, Vk_new)
  
  # U-update for Vnew
  Unew <- t(sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=Vnew, gamma=0, kappa_u=kappa_u)))
  # Unew <- t(apply(X, 1, function(x) Solve_U_GP(x=x, mu=mu, V=Vnew, gamma=0)))
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew, Vnew)
  
  loss <- sqrt(mean((X - xhat)^2))
  
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path) )
  # return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss) )
}















#' @export V_update2
V_update2 <- function(X, Uhat, Vhat, Uk, kappa_v=1e-8){
  if(FALSE){
    X=X; Uhat=Unew[,1:(r-1)]; Vhat=Vhat; Uk=Unew[,r]; kappa_v=kappa_v
  }
  
  n=nrow(X); p=ncol(X); r=NCOL(Uhat)
  mu=colMeans(X)
  
  C <- tcrossprod(rep(1,n),mu) + tcrossprod(Uhat, Vhat)
  
  Vnew <- rep(0,p)
  Dmat <- diag(rep(1,p))
  dvec <- t(X-tcrossprod(rep(1,n),mu)) %*% Uk/sum(Uk^2)
  dvec <- t(X-C) %*% Uk/sum(Uk^2)
  
  solve_V <- function(Vhat, lbmat=NULL, ubmat=NULL, kappa_v=0) {
    p=NROW(Vhat);  r=NCOL(Vhat)
    lb <- if(!is.null(lbmat)) apply(lbmat,2,max) else NULL
    ub <- if(!is.null(ubmat)) apply(ubmat,2,min) else NULL
    
    bvec <- c(rep(0,r+1), lb-kappa_v, -(ub+kappa_v))
    Amat <- cbind(rep(1,p), Vhat, if (!is.null(lbmat)) diag(rep(1,p)) else NULL, if (!is.null(ubmat)) -diag(rep(1,p)) else NULL)
    Vnew <- try( solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=r+1, factorized=FALSE)$solution, silent=TRUE )
    if(inherits(Vnew, "try-error")){
      return( "No solution" )
    } else {
      Vnew <- Vnew/norm(Vnew,"2")
      Vnew <- ifelse(abs(Vnew)<1e-8, 0, Vnew)
      return( Vnew/norm(Vnew,"2") )
    }
    # bvec <- c(rep(0,1), rep(-1e-12,2*r), lb-kappa_v, -(ub+kappa_v))
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
                    kappa_v=kappa_v)
  }
  
  if (n.pos>0 && n.neg==0) {
    Vnew <- solve_V(Vhat=Vhat,
                    lbmat=-diag(1/Uk[positiveindex],n.pos,n.pos) %*% C[positiveindex,], 
                    kappa_v=kappa_v)
  }
  
  if (n.pos==0 && n.neg>0) {
    Vnew <- solve_V(Vhat=Vhat,
                    ubmat=-diag(1/Uk[negativeindex],n.neg,n.neg) %*% C[negativeindex,], 
                    kappa_v=kappa_v)
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
Solve_U_GP <- function(x, mu, V, gamma=0, kappa_u=1e-8){
  if(FALSE){
    x=X[i,]; mu=C[i,]; V=cbind(Vhat, Vk_old); gamma=gamma/(it^(1/2))
  }
  
  x=as.numeric(x)
  mu=as.numeric(mu)
  
  U <- solve.QP(Dmat=crossprod(V), dvec=t(V) %*% (x - mu), Amat=t(V), bvec=-mu - kappa_u)$solution
  
  return(U * (1-gamma))
}

