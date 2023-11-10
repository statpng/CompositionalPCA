
#' @export png.ppca_qp
png.ppca_qp <- function(X, nrank=2, maxit=500, eps=1e-8, kappa=1e-8, gamma=1e-1, phi=0.01, V.init=NULL, verbose=FALSE){
  if(FALSE){
    X; nrank=2; maxit=500; eps=1e-8; kappa=1e-8; gamma=1e-2; phi=0.01; V.init=c("PC","random"); verbose=FALSE
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
      fit.rank1 <- png.rank12(X, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
      Uhat <- cbind(Uhat, fit.rank1$uhat)
      Vhat <- cbind(Vhat, fit.rank1$vhat)
      
      fit.path[[iter]] <- fit.rank1
    } else {
      
      fit.UkVk <- update_UkVk(X=X, Uhat=Uhat, Vhat=Vhat, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
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
              kappa=kappa, 
              gamma=gamma)
  
  return(list(mu=mu, uhat=Uhat, vhat=Vhat, xhat=xhat, 
              X=X, fit.path=fit.path, maxit=maxit, method="ppca", params=params))
  
}













#' @export update_UkVk
update_UkVk <- function(X, Uhat, Vhat, maxit=500, eps=1e-8, kappa=1e-4, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  
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
  
  # png.angle(V0.PC, V0)[[1]] %>% print
  
  
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
    Vnew <- V_update2(X, Uhat=Uhat, Vhat=as.matrix(Vhat), Uk=as.matrix(Unew), kappa=kappa)
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

