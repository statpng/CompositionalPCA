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



#' @export sim.Linear
sim.Linear <- function(n, p, r, snr=2, d=10, d0=0, seed=1, seed.U=seed, seed.V=seed, alpha=NULL, eta=0, verbose=FALSE){
  # Simulate compositional data
  
  if(FALSE){
    snr=2; d=10; d0=0.1; seed=1; seed.U=seed; seed.V=seed; alpha=NULL
  }
  
  if(FALSE){
    n=100; p=4; r=1; d=10
    
    line=-4.5
    setEPS()
    postscript(file="Figure-sim.Linear.eps", width=10, height=6)
    # pdf(file="Figure-sim.Linear.pdf", width=10, height=6)
    par(mfrow=c(2,3), mar = rep(0.2, 4), mai=c(0.4,0,0,0), omi=c(0,0,0,0))
    snr=2; eta=-0.1
    sim.Linear(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0
    sim.Linear(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0.1
    sim.Linear(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    
    snr=1; eta=0.1
    sim.Linear(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0.1
    sim.Linear(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=5; eta=0.1
    sim.Linear(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
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
  X2 <- t(apply(X, 1, proj2simplex))
  X0 <- (tcrossprod(rep(1,n),mu) + tcrossprod(U,V)) %>% {
    t(apply(., 1, proj2simplex))
  }
  
  snr_out <- sum(tcrossprod(U,V)^2) / sum(E^2)
  if(verbose) print(paste0("delta=",round(delta_snr,4), "; snr=", round(snr_out,4)))
  
  idx_not_sum_to_unit <- which(apply(X2,1,sum)>0.9999 & apply(X2,1,sum)<1)
  X2[idx_not_sum_to_unit,] <- t(apply(X2[idx_not_sum_to_unit,,drop=F],1,function(x) x/sum(x)))
  
  if(verbose) print(paste0("seed=", seed))
  
  
  params <- list(n=n, p=p, r=r, snr=snr, d=d, d0=d0, seed=seed, seed.U=seed.U, seed.V=seed.V, alpha=alpha, eta=eta)
  
  result <- list(mu=mu, U=U, D=D, V=V, E=E, X0=X0, X=X, X2=X2, params=params, type="Linear")
  
    
  return( result )
  
}






#' @export sim.Linear.test
sim.Linear.test <- function(params){
  TestData <- params %>% 
    { sim.Linear(n=.[["n"]], 
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
      Xnew <- lrPCA(.$X2, nrank=2, zero.replace="simple", delta=1e-10)$Xnew
      quaternary3d(Xnew, vhat=.$Vlist, mu=.$mu)
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
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=5; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=5; d=5
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      
      snr=2; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=5; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[snr, d] = [%d, %d]", snr, d), side=3, line=line)
      snr=10; d=3
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
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
  X02 <- ICLR(X0)
  
  E <- matrix(rnorm(n*p,0,1),n,p)
  sigma <- sqrt(sum(as.numeric(X0)^2)/sum(as.numeric(E)^2)/snr)
  E <- E * sigma
  
  snr_out <- sum(tcrossprod(U,V)^2) / sum(E^2)
  if(verbose) print(paste0("sigma=",round(sigma,4), "; snr=", round(snr_out,4)))
  
  # Z <- matrix(0,n,p)
  # for(j in 1:p){
  #   eta_j <- runif(1,0,zero.prop)
  #   Z[,j] <- rbinom(n, size=1, prob=eta_j)
  # }
  
  X <- X0 + E
  X2 <- X %>% 
    ICLR() %>% {
      # .[Z==1] <- 0
      .[.<1e-2/log(p)] <- 0
      .
    } %>% 
    apply(1, function(x){
      x/sum(x)
    }) %>% t()
  
  
  # quaternary3d(X)
  
  mu0 <- iclr(mu)
  
  idx_not_sum_to_unit <- which(apply(X2,1,sum)>0.9999 & apply(X2,1,sum)<1)
  X2[idx_not_sum_to_unit,] <- t(apply(X2[idx_not_sum_to_unit,,drop=F],1,function(x) x/sum(x)))
  #
  #
  params <- list(n=n, p=p, r=r, snr=snr, d=d, d0=d0, seed=seed, seed.U=seed.U, seed.V=seed.V)
  #
  result <- list(mu=mu, U=U, D=D, V=V, Vlist=lapply(seq(-20,20,1), function(z) z*V), E=E, X0=X02, X=X, X2=X2, params=params, type="LogNormal")
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




















# Other Papers --------------------------------------------------------------




#' @export rmixnorm
rmixnorm <- function(n,mean,varcov,prop){
  if(FALSE){
    n=1
    mean=list(c(-1,1),c(2,1.5),c(0.5,-1.5))
    varcov=lapply(1:3,function(x)0.5*diag(2,2))
    prop=c(0.4,0.3,0.3)
  }
  
  cumprop <- cumsum(prop)
  p <- unique(sapply(mean,length))
  
  x <- matrix(NA,n,p)
  for( i in 1:n ){
    mixture_prop <- runif(1,0,1)
    if(mixture_prop<cumprop[1]){
      x[i,] <- mnormt::rmnorm(1, c(-1,1), 0.5*diag(2,2))
    } else if(mixture_prop<cumprop[2]) {
      x[i,] <- mnormt::rmnorm(1, c(2,1.5), 0.5*diag(2,2))
    } else if(mixture_prop<cumprop[3]) {
      x[i,] <- mnormt::rmnorm(1, c(0.5,-1.5), 0.5*diag(2,2))
    }
  }
  
  x
}



#' @export sim.pierson
sim.pierson <- function(n=100,p=50,r=2){
  if(FALSE){
    n = 100; p = 50; r = 2;
  }
  
  f_i <- matrix(0,nrow = n, ncol = r)
  beta_j <- matrix(0,nrow = p, ncol = r)
  phi_j <- 1
  
  # Pierson and Yau (2015): zero-inflated Gaussian model
  #
  # x_ij ~ N(mu_ij, sigma2_j)  if z_ij = 0
  #      ~ 0  if z_ij = 1
  #
  # f_i ~ N(0, 1)
  # beta_j ~ U(-4,4)
  # alpha0, beta0 ~ U(3, 5)
  # eta_j ~ U(0, 1)
  
  alpha <- rep(0,n)
  beta0 <- runif(p,2.7,3.3)
  sigma2_j <- runif(p, 0.27, 3.3)
  
  for(i in 1:n){
    f_i[i,] <- rnorm(r, mean = 0, sd = 1)
  }
  
  for(j in 1:p){
    beta_j[j,] <- runif(r, -0.5, 0.5)
  }
  
  
  mu <- matrix(alpha,n,p) + matrix(beta0,n,p,byrow=TRUE) + f_i %*% t(beta_j)
  X_star <- matrix(0,n,p,byrow=TRUE)
  for(i in 1:n){
    for(j in 1:p){
      X_star[i,j] <- rnorm(1, mu[i,j], sqrt(sigma2_j[j]))
    }
  }
  eta_ij <- exp(-0.1*X_star^2)
  
  z <- matrix(0,n,p)
  for(i in 1:n){
    for(j in 1:p){
      z[i,j] <- rbinom(1, size=phi_j, prob=eta_ij[i,j])
    }
  }
  
  X <- ifelse(z==1, 0, round(exp(X_star)))
  X[z==1]=0
  
  
  zerorow <- which(rowSums(X)==0)
  if(length(zerorow) >0 ){
    X <- X[-zerorow,];f_i <- f_i[-zerorow,]
  }
  zerocol <- which(colSums(X)==0)
  if(length(zerocol) >0 ){
    X <- X[,-zerocol]; beta_j <- t(t(beta_j)[,-zerocol])
  }
  
  list(X=X,
       X2=t(apply(X,1,function(x) x/sum(x))),
       f_i=f_i,
       beta_j=beta_j,
       beta0=beta0,
       alpha=alpha,
       mu=mu,
       eta_ij=eta_ij,
       z=z)
}



#' @export sim.niku
sim.niku <- function(n=100,p=50,r=2,type=c(1,2,3)){
  if(FALSE){
    n = 100; p = 50; r = 2;
  }
  
  f_i <- matrix(0,nrow = n, ncol = r)
  beta_j <- matrix(0,nrow = p, ncol = r)
  phi_j <- 1
  
  
  if(type == 1){
    # no zero inflation with % of zero = 50%
    # f_i ~ 0.4*N((-1,1),Sigma) + 0.3*N((2,1.5),Sigma) + 0.3*N((0.5,-1.5),Sigma), where Sigma = 0.5*diag(2)
    # beta_j ~ U(-2,2)
    # alpha0, beta0 ~ U(-1, 1)
    # phi_j = 1
    
    eta_j <- runif(p,0,0)
    
    alpha <- runif(n,-1,1)
    beta0 <- runif(p,-1,1)
    
    for(i in 1:n){
      # mixture bivariate normal with prop (0.4,0.3,0.3)
      f_i[i,] <- rmixnorm(n=1,
                          mean=list(c(-1,1),c(2,1.5),c(0.5,-1.5)),
                          varcov=lapply(1:3,function(x)0.5*diag(2,2)),
                          prop=c(0.4,0.3,0.3) )
    }
    
    for(j in 1:p){
      beta_j[j,] <- runif(r, -2, 2)
    }
    
  } else if(type == 2){
    # true absence of species (eta_j): structural zero percentage = 50%
    # f_i ~ 0.4*N((-1,1),Sigma) + 0.3*N((2,1.5),Sigma) + 0.3*N((0.5,-1.5),Sigma), where Sigma = 0.5*diag(2)
    # beta_j ~ U(-5,5)
    # alpha0, beta0 ~ U(0, 2)
    # phi_j = 1
    # eta_j ~ U(0, 0.5)
    
    eta_j <- runif(p,0,0.5)
    
    alpha <- runif(n,0,2)
    beta0 <- runif(p,0,2)
    
    for(i in 1:n){
      # mixture bivariate normal with prop (0.4,0.3,0.3)
      f_i[i,] <- rmixnorm(n=1,
                          mean=list(c(-1,1),c(2,1.5),c(0.5,-1.5)),
                          varcov=lapply(1:3,function(x)0.5*diag(2,2)),
                          prop=c(0.4,0.3,0.3) )
    }
    
    for(j in 1:p){
      beta_j[j,] <- runif(r, -5, 5)
    }
    
  } else if(type == 3){
    # high structural zero percentage = 90%
    # f_i ~ 0.4*N((-1,1),Sigma) + 0.3*N((2,1.5),Sigma) + 0.3*N((0.5,-1.5),Sigma), where Sigma = 0.5*diag(2)
    # beta_j ~ U(-4,4)
    # alpha0, beta0 ~ U(3, 5)
    # phi_j = 1
    # eta_j ~ U(0, 1)
    
    eta_j <- runif(p,0,1)
    
    alpha <- runif(n,3,5)
    beta0 <- runif(p,3,5)
    
    for(i in 1:n){
      # mixture bivariate normal with prop (0.4,0.3,0.3)
      f_i[i,] <- rmixnorm(n=1,
                          mean=list(c(-1,1),c(2,1.5),c(0.5,-1.5)),
                          varcov=lapply(1:3,function(x)0.5*diag(2,2)),
                          prop=c(0.4,0.3,0.3) )
    }
    
    for(j in 1:p){
      beta_j[j,] <- runif(r, -4, 4)
    }
    
  }
  
  
  
  
  log_mu <- matrix(alpha,n,p) + matrix(beta0,n,p,byrow=TRUE) + f_i %*% t(beta_j)
  mu <- exp(log_mu)
  
  z <- matrix(0,n,p)
  for(i in 1:n){
    z[i,] <- rbinom(p, size=1, prob=eta_j)
  }
  
  X <- matrix(0,n,p,byrow = TRUE)
  for(i in 1:n){
    for(j in 1:p){
      X[i,j] <- rnbinom(n=1, size=phi_j, mu=mu[i,j])
    }
  }
  X[z==1]=0
  
  
  zerorow <- which(rowSums(X)==0)
  if(length(zerorow) >0 ){
    X <- X[-zerorow,];f_i <- f_i[-zerorow,]
  }
  zerocol <- which(colSums(X)==0)
  if(length(zerocol) >0 ){
    X <- X[,-zerocol]; beta_j <- t(t(beta_j)[,-zerocol])
  }
  
  list(X=X,
       X2=t(apply(X,1,function(x) x/sum(x))),
       f_i=f_i,
       beta_j=beta_j,
       beta0=beta0,
       alpha=alpha,
       mu=mu,
       eta_j=eta_j,
       z=z)
}



#' @export sim.zippca
sim.zippca <- function(n=100,p=50,r=2, phi_j=10, sigma=1){
  if(FALSE){
    n = 100; p = 50; r = 2;
  }
  
  # Zero-Inflated Negative Binomial Distributions
  # f_i ~ N(0,1)
  # beta_j ~ U(-2,2)
  # alpha0, beta0 ~ U(-1, 1)
  # phi_j = 1
  
  
  f_i <- matrix(0,nrow = n, ncol = r)
  beta_j <- matrix(0,nrow = p, ncol = r)
  
  beta0 <- rep(1,p)
  alpha <- rep(0,n)
  eta_j <- runif(p,0,1)
  
  for(i in 1:n){
    f_i[i,] <- rnorm(r, mean = 0, sd = sigma)
  }
  
  for(j in 1:p){
    beta_j[j,] <- rnorm(r, mean =0, sd = sigma)
  }
  
  
  log_mu <- matrix(alpha,n,p) + matrix(beta0,n,p,byrow=TRUE) + f_i %*% t(beta_j)
  mu <- exp(log_mu)
  
  z <- matrix(0,n,p)
  for(i in 1:n){
    z[i,] <- rbinom(p, size=1, prob=eta_j)
  }
  
  X <- matrix(0,n,p,byrow = TRUE)
  for(i in 1:n){
    for(j in 1:p){
      X[i,j] <- rnbinom(n=1, size=phi_j, mu=mu[i,j])
    }
  }
  X[z==1]=0
  
  zerorow <- which(rowSums(X)==0)
  if(length(zerorow) >0 ){
    X <- X[-zerorow,]; f_i <- f_i[-zerorow,]
  }
  zerocol <- which(colSums(X)==0)
  if(length(zerocol) >0 ){
    X <- X[,-zerocol]; beta_j <- t(t(beta_j)[,-zerocol])
  }
  
  list(X=X,
       X2=t(apply(X,1,function(x) x/sum(x))),
       f_i=f_i,
       beta_j=beta_j,
       beta0=beta0,
       alpha=alpha,
       mu=mu,
       eta_j=eta_j,
       z=z)
}





#' @export simplot.rank2
simplot.rank2 <- function(n,p,r,f=sim.zippca, xlim=NULL, ylim=NULL, ...){
  
  zero_prop <- mean(f(n=n, p=p, r=r, ...)$X2<1e-12)
  
  f(n=n, p=p, r=r, ...)$X2 %>%
    pca.approx(nrank=2) %>%
    {
      . <- ifelse( abs(.)<1e-10, 0, .)
      . <- ifelse( abs(.-1)<1e-10, 1, .)
      outliers <- which( apply(.,1,function(x) any(x<0|x>1)) )
      Color <- ifelse(1:n %in% outliers, "red", "black")
      zero_prop <- apply(.,1,function(x) any(x<0)) %>% mean()
      
      {
        # print(order(apply(.,2,sd),decreasing=TRUE)[1:2])
        # .[,order(apply(.,2,sd),decreasing=TRUE)[1:2]]
        prcomp(.)$x[,1:2]
      } %>%
        plot(cex=1, pch=18, col=Color, xlab="PC 1", ylab="PC 2",
             main=paste0("(n, p, r) = (",n,", ",p,", ",r,")"),
             sub=paste0("Outside the simplex = ",round(zero_prop*100,1),"%"), xlim=xlim, ylim=ylim)
    }
  
  
  
}


