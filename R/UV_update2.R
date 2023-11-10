#' @export png.rank12
png.rank12 <- function(X, maxit=500, eps=1e-8, kappa=1e-4, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  
  if(FALSE){
    data <- sim.simplex(n=500, p=100, r=2, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE); 
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
  
  # png.angle(V0.PC, V0)[[1]] %>% print
  
  
  est.path <- crit.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vold <- V0
    } else {
      Vold <- Vnew
    }
    
    
    Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vold, gamma=gamma/(it^(1/2))))
    
    Vnew <- V_update2(X, Uhat=matrix(0,n,1), Vhat=matrix(0,p,1), Uk=as.matrix(Unew), kappa=kappa)
    
    crit <- min( sum((Vnew-Vold)^2), sum((Vnew+Vold)^2) )
    
    crit.path[it] <- crit
    # est.path[[it]] <- list(uhat=Unew, vhat=Vnew)
    
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
  
  # return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path) )
  return( list(# xhat=xhat, 
    mu=mu, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss) )
  
}











#' @export UV_update2
UV_update2 <- function(X, Uhat, Vhat, maxit=500, eps=1e-6, kappa=1e-4, gamma=0, nu=1e-8, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
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
  
  # png.angle(V0.PC, V0)[[1]] %>% print
  
  crit.path <- est.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vk_old <- V0
    } else {
      Vk_old <- Vk_new
    }
    
    
    
    # U-update
    Unew <- t(sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=cbind(Vhat, Vk_old), gamma=gamma/(it^(1/2)), nu=nu)))
    
    # V-update
    Vk_new <- V_update2(X, Uhat=Unew[,1:(r-1)], Vhat=Vhat, Uk=Unew[,r], kappa=kappa)
    if( Vk_new[1] == "No solution" ){
      Vk_new <- Vk_old
      NoSolution <- TRUE
    } else {
      NoSolution <- FALSE
    }
    
    crit <- min( (sum((Vk_old-Vk_new)^2)), (sum((Vk_old+Vk_new)^2)) )
    
    crit.path[it] <- crit
    # est.path[[it]] <- list(uhat=Unew, vhat=cbind(Vhat,Vk_new))
    
    if( crit < eps ) break
  }
  
  if(NoSolution) crit.path[it] <- 0.1
  
  Vnew <- cbind(Vhat, Vk_new)
  
  # U-update for Vnew
  Unew <- t(sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=Vnew, gamma=0, nu=nu)))
  # Unew <- t(apply(X, 1, function(x) Solve_U_GP(x=x, mu=mu, V=Vnew, gamma=0)))
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew, Vnew)
  
  loss <- sqrt(mean((X - xhat)^2))
  
  # return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path) )
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss) )
}
