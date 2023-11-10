#' @export png.rmixnorm
png.rmixnorm <- function(n,mean,varcov,prop){
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
      f_i[i,] <- png.rmixnorm(n=1,
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
      f_i[i,] <- png.rmixnorm(n=1,
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
      f_i[i,] <- png.rmixnorm(n=1,
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
    png.pca.approx(nrank=2) %>%
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

