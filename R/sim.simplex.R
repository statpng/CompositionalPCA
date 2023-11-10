#' @export calc_snr
calc_snr <- function(delta, mu,U,V, seed1){
  n=nrow(U); p=nrow(V)
  set.seed(seed1)
  E <- matrix(runif(n*p,-delta,delta),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- tcrossprod(rep(1,n),mu) + tcrossprod(U,V) + E
  sum(tcrossprod(U,V)^2) / sum(E^2)
}

#' @export calc_snr2
calc_snr2 <- function(delta, mu, X0, seed1){
  n=nrow(X0); p=ncol(X0)
  UV <- X0 - tcrossprod(rep(1,n),mu) # scale(X0, center=TRUE, scale=FALSE)
  
  set.seed(seed1)
  E <- matrix(runif(n*p,-delta,delta),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- UV + E
  sum(UV^2) / sum(E^2)
}



#' @export sim.simplex
sim.simplex <- function(n, p, r, snr=2, d=10, d0=0, seed=1, seed.U=seed, seed.V=seed, alpha=NULL, eta=0, verbose=FALSE){
  # Simulate compositional data
  
  if(FALSE){
    snr=2; d=10; d0=0.1; seed=1; seed.U=seed; seed.V=seed; alpha=NULL
  }
  
  if(FALSE){
    n=100; p=4; r=1; d=10
    
    line=-4.5
    setEPS()
    postscript(file="Figure-sim.simplex.eps", width=10, height=6)
    # pdf(file="Figure-sim.simplex.pdf", width=10, height=6)
    par(mfrow=c(2,3), mar = rep(0.2, 4), mai=c(0.4,0,0,0), omi=c(0,0,0,0))
    snr=2; eta=-0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    
    snr=1; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=5; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    
    mtext(text=sprintf("[n, p, r, d] = [%d, %d, %d, %d]", n,p,r,d), side=3, outer=TRUE, line=-2) # adj=0.02
    
    dev.off()
  }
  
  
  
  
  
  
  if(is.null(alpha)){
    alpha <- rep(10,p)
  }
  
  if(sum(alpha)==1){
    mu <- alpha
  } else {
    set.seed(seed.U)
    mu <- as.vector(dirmult::rdirichlet(1, alpha=alpha))
  }
  
  
  set.seed(seed.V)
  V <- qr.Q(qr(cbind(1, do.call("cbind", lapply(1:r, function(x) rnorm(p))) )))[,-1,drop=F]
  
  
  D <- sapply(1:r, function(k) d/k + d0)
  
  set.seed(seed.U)
  U <- matrix(0,n,r)
  for( k in 1:r ){
    v = V[,k];
    start=onedimconvexprojection(mu,-100*v,v); end=onedimconvexprojection(mu,100*v,v)
    U[,k]=truncnorm::rtruncnorm(n, a=start-eta, b=end+eta, mean=0, sd=D[k])
  }
  
  
  delta_snr <- optimize(function(x, mu=mu,U=U,V=V, seed2) norm(calc_snr(delta=x, mu=mu,U=U,V=V, seed1=seed2)-snr, "2"), c(0,1), tol=1e-10, mu=mu,U=U,V=V, seed2=seed.U)$minimum
  
  
  set.seed(seed.U)
  E <- matrix(runif(n*p,-delta_snr,delta_snr),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- tcrossprod(rep(1,n),mu) + tcrossprod(U,V) + E
  X2 <- t(apply(X, 1, png.proj2simplex))
  X0 <- (tcrossprod(rep(1,n),mu) + tcrossprod(U,V)) %>% {
    t(apply(., 1, png.proj2simplex))
  }
  
  snr_out <- sum(tcrossprod(U,V)^2) / sum(E^2)
  if(verbose) print(paste0("delta=",round(delta_snr,4), "; snr=", round(snr_out,4)))
  
  idx_not_sum_to_unit <- which(apply(X2,1,sum)>0.9999 & apply(X2,1,sum)<1)
  X2[idx_not_sum_to_unit,] <- t(apply(X2[idx_not_sum_to_unit,,drop=F],1,function(x) x/sum(x)))
  
  if(verbose) print(paste0("seed=", seed))
  
  
  params <- list(n=n, p=p, r=r, snr=snr, d=d, d0=d0, seed=seed, seed.U=seed.U, seed.V=seed.V, alpha=alpha, eta=eta)
  
  result <- list(mu=mu, U=U, D=D, V=V, E=E, X0=X0, X=X, X2=X2, params=params)
  
    
  return( result )
  
}






#' @export sim.simplex.test
sim.simplex.test <- function(params){
  TestData <- params %>% 
    { sim.simplex(n=.[["n"]], 
                  p=.[["p"]], 
                  r=.[["r"]], 
                  snr=.[["snr"]], 
                  d=.[["d"]],
                  d0=.[["d0"]],
                  seed.U=.[["seed.U"]]*123,
                  seed.V=.[["seed.V"]],
                  alpha=.[["alpha"]],
                  eta=.[["eta"]]) }
  
  TestData
}







#' @export sim.LogNormal
sim.LogNormal <- function(n, p, r, snr=5, d=3, d0=0, zero.prop=0.1, seed=1, seed.U=seed, seed.V=seed, verbose=FALSE){
  # Simulate compositional data
  
  if(FALSE){
    n=500; p=4; r=1; snr=2; d=10; d0=0; seed=1; seed.U=seed; seed.V=seed; verbose=FALSE
    
    data <- sim.LogNormal(n=n, p=p, r=r, snr=10, d=50, d0=0, zero.prop=0.1, seed=1, verbose=TRUE)
    mean(data$X2==0)
    data$X2 %>% apply(1,function(x) any(x==0)) %>% mean
    
    
    lapply(1:100, function(i) sim.LogNormal(n=50, p=100, r=5, snr=5, d=10, d0=0, seed=i, verbose=TRUE)) %>% 
      sapply(function(x) x$X2 %>% apply(1,function(x) any(x==0)) %>% mean ) %>% mean
    
    lapply(1:100, function(i) sim.LogNormal(n=50, p=200, r=5, snr=5, d=10, d0=0, seed=i, verbose=TRUE)) %>% 
      sapply(function(x) x$X2 %>% apply(1,function(x) any(x==0)) %>% mean ) %>% mean
    
    
    lapply(1:100, function(i) sim.LogNormal(n=50, p=100, r=5, snr=5, d=3, d0=0, seed=i, verbose=TRUE)) %>% 
      sapply(function(x) mean(x$X2==0) ) %>% mean
    
    lapply(1:100, function(i) sim.LogNormal(n=50, p=200, r=5, snr=5, d=3, d0=0, seed=i, verbose=TRUE)) %>% 
      sapply(function(x) mean(x$X2==0) ) %>% mean
    
    
    data %>% {
      Xnew <- png.lrpca(.$X2, nrank=2, zero.replace="simple", delta=1e-10)$Xnew
      png.quaternary3d(Xnew, vhat=.$Vlist, mu=.$mu)
    }
    
    data$X2 %>% apply(1,function(x) any(x==0)) %>% mean
    
    if(FALSE){
      n=100; p=4; r=1;
      
      line=-4.5
      setEPS()
      postscript(file="Figure-sim.LogNormal.eps", width=10, height=6)
      # pdf(file="Figure-sim.LogNormal.pdf", width=10, height=6)
      par(mfrow=c(2,3), mar = rep(0.2, 4), mai=c(0.4,0,0,0), omi=c(0,0,0,0))
      snr=5; d=1
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=5; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=5; d=5
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      
      snr=2; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=5; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=10; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      
      mtext(text=sprintf("[n, p, r] = [%d, %d, %d]", n,p,r), side=3, outer=TRUE, line=-2) # adj=0.02
      
      dev.off()
    }
    
    
  }
  
  
  
  mu <- rep(0,p)
  
  set.seed(seed.V)
  V <- qr.Q(qr( do.call("cbind", lapply(1:r, function(x) rnorm(p))) ))
  D <- sapply(1:r, function(k) (d/k + d0))
  
  set.seed(seed.U)
  U <- matrix(0,n,r)
  for( k in 1:r ){
    U[,k]=rnorm(n, mean=0, sd=D[k])
  }
  
  X0 <- tcrossprod(rep(1,n),mu) + tcrossprod(U,V)
  X02 <- png.ICLR(X0)
  
  E <- matrix(rnorm(n*p,0,1),n,p)
  sigma <- sqrt(sum(as.numeric(X0)^2)/sum(as.numeric(E)^2)/snr)
  E <- E * sigma
  
  # Z <- matrix(0,n,p)
  # for(j in 1:p){
  #   eta_j <- runif(1,0,zero.prop)
  #   Z[,j] <- rbinom(n, size=1, prob=eta_j)
  # }
  
  X <- X0 + E
  X2 <- X %>% 
    png.ICLR() %>% {
      # .[Z==1] <- 0
      .[.<1e-2/log(p)] <- 0
      .
    } %>% 
    apply(1, function(x){
      x/sum(x)
    }) %>% t()
  
  
  # png.quaternary3d(X)
  
  mu0 <- png.iclr(mu)
  
  idx_not_sum_to_unit <- which(apply(X2,1,sum)>0.9999 & apply(X2,1,sum)<1)
  X2[idx_not_sum_to_unit,] <- t(apply(X2[idx_not_sum_to_unit,,drop=F],1,function(x) x/sum(x)))
  #
  #
  params <- list(n=n, p=p, r=r, snr=snr, d=d, d0=d0, seed=seed, seed.U=seed.U, seed.V=seed.V)
  #
  result <- list(mu=mu, U=U, D=D, V=V, Vlist=lapply(seq(-20,20,1), function(z) z*V), E=E, X0=X02, X=X, X2=X2, params=params)
  #
  return( result )
  
}



#' @export sim.LogNormal.test
sim.LogNormal.test <- function(params){
  TestData <- params %>% 
    { sim.LogNormal(n=.[["n"]], 
                    p=.[["p"]], 
                    r=.[["r"]], 
                    snr=.[["snr"]], 
                    d=.[["d"]],
                    d0=.[["d0"]],
                    seed.U=.[["seed.U"]]*123,
                    seed.V=.[["seed.V"]] ) }
  
  TestData
}






