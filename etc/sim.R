detach("package:png.compositionalPCA",unload=TRUE)
library(png.compositionalPCA)


{
  library(parallel)

  numCores <- detectCores() - 1  # leave one core free
  numCores <- 4

  {
    n.seq <- c(50,100,500,1000)
    p.seq <- c(100,500)
    eta.seq <- c(0, 0.1)
    seed.seq <- c(1:10)
    param.list <- list(n.seq=n.seq,
                       p.seq=p.seq,
                       eta.seq=eta.seq,
                       seed.seq=seed.seq)

    # Create the indices
    grid <- expand.grid(param.list)
  }
}



# consistency
function(){
  {
    #1: eta=0
    eta=0

    Fnorm_V <- function(V,vhat){
      if(FALSE){
        V=data$V; vhat=fit$vhat
      }
      r=NCOL(V)
      out <- NULL
      for( k in 1:r ){
        out[k] <- min(sum((V[,k]-vhat[,k])^2),
                      sum((V[,k]+vhat[,k])^2))
      }
      sqrt(sum(out))
    }



    n.seq <- c(10,100,1000,5000)
    eta.seq <- c(-0.1, 0, 0.1)
    seed.seq <- 1:10

    out <- array(NA, dim=c(4,3,2,3,10), dimnames=list( paste0("n=",c(10,100,1000,5000)), c("V vs PC", "V vs GPC", "PC vs GPC"), c("angle", "fnorm"), paste0("eta=", c(-0.1, 0, 0.1)), paste0("seed=", seed.seq) ))

    for( k in 1:length(seed.seq) ){
      seed <- seed.seq[k]
      for( j in 1:length(eta.seq) ){
        eta=eta.seq[j]

        for( i in 1:length(n.seq) ){
          n=n.seq[i]

          p=10; r=2; d=1; d0=0.1; snr=5; seed=seed
          data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=d,d0=d0,seed=seed,eta=eta )
          X <- data$X2
          fit <- png.gppca_qp(X, nrank=r, kappa=1e-4)

          out[i,1,1,j,k] <- png.utils::png.angle(data$V, prcomp(X)$rot[,1:2])$max
          out[i,2,1,j,k] <- png.utils::png.angle(data$V, fit$vhat)$max
          out[i,3,1,j,k] <- png.utils::png.angle(prcomp(X)$rot[,1:2], fit$vhat)$max

          out[i,1,2,j,k] <- Fnorm_V(data$V, prcomp(X)$rot[,1:2])
          out[i,2,2,j,k] <- Fnorm_V(data$V, fit$vhat)
          out[i,3,2,j,k] <- Fnorm_V(prcomp(X)$rot[,1:2], fit$vhat)
        }
      }
    }

    apply(out,1:4,mean)
    #
    #

  }

}








function(){

  {


    {
      n=100; p=4; r=2; snr=5; seed=1; eta=0
      X <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta)$X2
      fit1_a <- png.lrpca(X, nrank=r, zero.replace="simple")
      fit1_b <- png.lrpca(X, nrank=r, zero.replace="additive")
      fit1_c <- png.lrpca(X, nrank=r, zero.replace="multiplicative")
      # png.quaternary3d(fit1$Xnew, vhat=fit1$vhat, xhat=fit1$xhat)
      fit2 <- png.ppca(X, nrank=r)
      fit3 <- png.gppca(X, nrank=r)
      kappa=1e-4; gamma=0.2
      fit4 <- png.ppca_qp(X, nrank=r, kappa=kappa, gamma=gamma)
      fit5 <- png.gppca_qp(X, nrank=r, kappa=kappa, gamma=gamma)

      png.ternary(fit1_a$Xnew, vhat=fit1_a$vhat)
      png.ternary(X, vhat=fit4$vhat)
      png.ternary(X, vhat=fit5$vhat)

      png.ternary(fit1_a$Xnew, vhat=fit1_a$vhat)
      png.ternary(X, vhat=fit4$vhat)
      png.quaternary3d(X, vhat=fit4$vhat)
      png.quaternary3d(X, vhat=fit5$vhat)
    }




    # main function
    sim1_convergence <- function(n, p, r, snr, eta, seed){
      if(FALSE){
        n=100; p=4; r=2; snr=2; eta=0.1/log(p); seed=1
        n=100; p=50; r=5; snr=5; eta=0.0; seed=1
      }

      data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta/log(p))
      X <- data$X2

      kappa=1e-6; gamma=0.5
      fit4 <- png.ppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma)
      fit5 <- png.gppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma)

      # png.pca.plot_convergence(fit4)
      # png.pca.plot_convergence(fit5)

      return( list(fit4, fit5) %>% map(~png.pca.convergence(.x)) )

    }




    res1_convergence <- function(){
      # n.seq eta.seq seed.seq
      # 1     10     0.0        1
      # 2    100     0.0        1
      # 3    500     0.0        1
      # 4   1000     0.0        1
      # 5     10     0.1        1
      # 6    100     0.1        1
      # 7    500     0.1        1
      # 8   1000     0.1        1
      # 9     10     0.0        2
      # 10   100     0.0        2
      # 11   500     0.0        2
      # 12  1000     0.0        2
      # 13    10     0.1        2
      # 14   100     0.1        2
      # 15   500     0.1        2
      # 16  1000     0.1        2

      {
        n.seq <- c(50,100,500)
        p.seq <- c(100,200)
        eta.seq <- c(0, 0.1)
        seed.seq <- c(1:10)
        param.list <- list(n.seq=n.seq,
                           p.seq=p.seq,
                           eta.seq=eta.seq,
                           seed.seq=seed.seq)
        grid <- expand.grid(param.list)

        res1 <- mcmapply(function(n, p, r, snr, eta, seed) sim1_convergence(n=n, p=p, r=r, snr=snr, eta=eta/log(p), seed=seed),
                         n=grid$n.seq, p=grid$p.seq, eta=grid$eta.seq, seed=grid$seed.seq, r=5, snr=2,
                         mc.cores = 4)

        save(res1, file="./res1.RData")
      }




      # Reshape the res1 to a matrix form
      n.out <- length(res1) / prod(sapply(param.list,length))
      res1_array <- array(res1,
                             dim=c(n.out, sapply(param.list,length)),
                             dimnames=append(list(out=paste0("V", 1:n.out)), param.list) )

      dimnames(res1_array)

      # eta=0
      # ppca
      ## n=10
      res1_array[1,1,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[1,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[1,3,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=1000
      res1_array[1,4,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      # gppca
      ## n=10
      res1_array[2,1,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[2,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[2,3,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=1000
      res1_array[2,4,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)


      # eta=0
      # ppca
      ## n=10
      res1_array[1,1,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[1,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[1,3,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=1000
      res1_array[1,4,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      # gppca
      ## n=10
      res1_array[2,1,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[2,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[2,3,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=1000
      res1_array[2,4,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)

    }







    sim2_kappa <- function(n, p, r, snr, eta, seed){

      data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta/log(p))
      X <- data$X2

      gamma=0.0
      fitA_0.0 <- png.ppca_qp(X, nrank=r, kappa=1e-0, maxit=500, eps=1e-6, gamma=gamma)
      fitA_0.1 <- png.ppca_qp(X, nrank=r, kappa=1e-1, maxit=500, eps=1e-6, gamma=gamma)
      fitA_0.2 <- png.ppca_qp(X, nrank=r, kappa=1e-2, maxit=500, eps=1e-6, gamma=gamma)
      fitA_0.3 <- png.ppca_qp(X, nrank=r, kappa=1e-3, maxit=500, eps=1e-6, gamma=gamma)
      fitA_0.4 <- png.ppca_qp(X, nrank=r, kappa=1e-4, maxit=500, eps=1e-6, gamma=gamma)

      gamma=0.5
      fitA_5.0 <- png.ppca_qp(X, nrank=r, kappa=1e-0, maxit=500, eps=1e-6, gamma=gamma)
      fitA_5.1 <- png.ppca_qp(X, nrank=r, kappa=1e-1, maxit=500, eps=1e-6, gamma=gamma)
      fitA_5.2 <- png.ppca_qp(X, nrank=r, kappa=1e-2, maxit=500, eps=1e-6, gamma=gamma)
      fitA_5.3 <- png.ppca_qp(X, nrank=r, kappa=1e-3, maxit=500, eps=1e-6, gamma=gamma)
      fitA_5.4 <- png.ppca_qp(X, nrank=r, kappa=1e-4, maxit=500, eps=1e-6, gamma=gamma)

      gamma=0.0
      fitB_0.0 <- png.gppca_qp(X, nrank=r, kappa=1e-0, maxit=500, eps=1e-6, gamma=gamma)
      fitB_0.1 <- png.gppca_qp(X, nrank=r, kappa=1e-1, maxit=500, eps=1e-6, gamma=gamma)
      fitB_0.2 <- png.gppca_qp(X, nrank=r, kappa=1e-2, maxit=500, eps=1e-6, gamma=gamma)
      fitB_0.3 <- png.gppca_qp(X, nrank=r, kappa=1e-3, maxit=500, eps=1e-6, gamma=gamma)
      fitB_0.4 <- png.gppca_qp(X, nrank=r, kappa=1e-4, maxit=500, eps=1e-6, gamma=gamma)

      gamma=0.5
      fitB_5.0 <- png.gppca_qp(X, nrank=r, kappa=1e-0, maxit=500, eps=1e-6, gamma=gamma)
      fitB_5.1 <- png.gppca_qp(X, nrank=r, kappa=1e-1, maxit=500, eps=1e-6, gamma=gamma)
      fitB_5.2 <- png.gppca_qp(X, nrank=r, kappa=1e-2, maxit=500, eps=1e-6, gamma=gamma)
      fitB_5.3 <- png.gppca_qp(X, nrank=r, kappa=1e-3, maxit=500, eps=1e-6, gamma=gamma)
      fitB_5.4 <- png.gppca_qp(X, nrank=r, kappa=1e-4, maxit=500, eps=1e-6, gamma=gamma)


      prcomp(X)$rot[,1:r] %>% head
      fit5_0$vhat %>% head
      fit5_1$vhat %>% head
      fit5_2$vhat %>% head
      fit5_3$vhat %>% head
      fit5_4$vhat %>% head

      png.quaternary3d(X, vhat=png.pca(X)$vhat)
      png.quaternary3d(X, vhat=fit5_0$vhat)

      fit5_0$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_1$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_2$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_3$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_4$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")

      fit5_0$fit.path[[2]]$crit.path %>% tail
      fit5_1$fit.path[[2]]$crit.path %>% tail
      fit5_2$fit.path[[2]]$crit.path %>% tail
      fit5_3$fit.path[[2]]$crit.path %>% tail
      fit5_4$fit.path[[2]]$crit.path %>% tail

      norm(fit5_1$xhat - fit5_2$xhat)/length(fit5_2$xhat)
      norm(X - fit5_0$xhat)/length(fit5_2$xhat)
      norm(X - fit5_1$xhat)/length(fit5_2$xhat)
      norm(X - fit5_2$xhat)/length(fit5_2$xhat)
      norm(X - fit5_3$xhat)/length(fit5_2$xhat)
      norm(X - fit5_4$xhat)/length(fit5_2$xhat)

      png.utils::png.angle(prcomp(X)$rot[,1:r], fit5_0$vhat)
      png.utils::png.angle(data$V, prcomp(X)$rot[,1:r])
      png.utils::png.angle(data$V, fit5_0$vhat)
      png.utils::png.angle(data$V, fit5_1$vhat)
      png.utils::png.angle(data$V, fit5_2$vhat)
      png.utils::png.angle(data$V, fit5_3$vhat)
      png.utils::png.angle(data$V, fit5_4$vhat)
    }



    sim.gamma <- function(){
      n=50; p=100; r=5; snr=2; eta=0.5; seed=1

      data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta)
      X <- data$X2

      kappa=10; maxit=500; eps=1e-6
      fit5_0 <- png.gppca_qp(X, nrank=r, maxit=maxit, eps=eps, kappa=kappa, gamma=0)
      fit5_1 <- png.gppca_qp(X, nrank=r, maxit=maxit, eps=eps, kappa=kappa, gamma=0.1)
      fit5_2 <- png.gppca_qp(X, nrank=r, maxit=maxit, eps=eps, kappa=kappa, gamma=0.2)
      fit5_3 <- png.gppca_qp(X, nrank=r, maxit=maxit, eps=eps, kappa=kappa, gamma=0.3)
      fit5_4 <- png.gppca_qp(X, nrank=r, maxit=maxit, eps=eps, kappa=kappa, gamma=0.9)

      fit6_0 <- png.ppca_qp(X, nrank=r, maxit=maxit, eps=eps, kappa=1e-1, gamma=0.0)

      fit5_0$time
      fit5_1$time
      fit5_2$time
      fit5_3$time
      fit5_4$time

      X[1:10,]
      fit5_0$xhat[1:10,]
      fit5_1$xhat[1:10,]
      fit5_2$xhat[1:10,]
      fit5_3$xhat[1:10,]
      fit5_4$xhat[1:10,]

      fit5_0$vhat %>% head
      fit5_1$vhat %>% head
      fit5_2$vhat %>% head
      fit5_3$vhat %>% head
      fit5_4$vhat %>% head

      fit5_0$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_1$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_2$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_3$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")
      fit5_4$fit.path[[2]]$crit.path[-(1:9)] %>% plot(type="l")


      fit5_0$fit.path[[2]]$crit.path %>% tail
      fit5_1$fit.path[[2]]$crit.path %>% tail
      fit5_2$fit.path[[2]]$crit.path %>% tail
      fit5_3$fit.path[[2]]$crit.path %>% tail
      fit5_4$fit.path[[2]]$crit.path %>% tail
      #

      norm(fit5_1$xhat - fit5_2$xhat)/length(fit5_2$xhat)
      norm(X - fit5_0$xhat)/length(fit5_2$xhat)
      norm(X - fit5_1$xhat)/length(fit5_2$xhat)
      norm(X - fit5_2$xhat)/length(fit5_2$xhat)
      norm(X - fit5_3$xhat)/length(fit5_2$xhat)
      norm(X - fit5_4$xhat)/length(fit5_2$xhat)

      png.utils::png.angle(data$V, fit5_0$vhat)
      png.utils::png.angle(data$V, fit5_1$vhat)
      png.utils::png.angle(data$V, fit5_2$vhat)
      png.utils::png.angle(data$V, fit5_3$vhat)
      png.utils::png.angle(data$V, fit5_4$vhat)

    }




    sim2 <- function(n, p, r, snr, eta, seed){
      if(FALSE){
        n=100; p=5; r=2; snr=2; eta=0.1; seed=1
        n=50; p=4; r=1; snr=2; eta=0.5; seed=1
      }
      data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta)
      X <- data$X2

      fit1_a <- png.lrpca(X, nrank=r, zero.replace="simple")
      fit1_b <- png.lrpca(X, nrank=r, zero.replace="additive")
      fit1_c <- png.lrpca(X, nrank=r, zero.replace="multiplicative")
      # fit1_a %>% { png.quaternary3d(.$Xnew, vhat=.$vhat, xhat=.$xhat) }
      fit2 <- png.ppca(X, nrank=r)
      fit3 <- png.gppca(X, nrank=r)
      kappa=1e-2; gamma=0.5; maxit=500; eps=1e-8
      fit4 <- png.ppca_qp(X, nrank=r, kappa=kappa, maxit=maxit, eps=eps, gamma=gamma)
      fit5 <- png.gppca_qp(X, nrank=r, kappa=kappa, maxit=maxit, eps=eps, gamma=gamma)

      eps=1e-8
      fit5_0 <- png.gppca_qp(X, nrank=2, kappa=10, eps=eps, maxit=500, gamma=gamma)
      fit5_1 <- png.gppca_qp(X, nrank=2, kappa=1, eps=eps, maxit=500, gamma=gamma)
      fit5_2 <- png.gppca_qp(X, nrank=2, kappa=1e-2, eps=eps, maxit=500, gamma=gamma)
      fit5_3 <- png.gppca_qp(X, nrank=2, kappa=1e-4, eps=eps, maxit=500, gamma=gamma)
      fit5_4 <- png.gppca_qp(X, nrank=2, kappa=1e-8, eps=eps, maxit=500, gamma=gamma)
      fit5_5 <- png.gppca_qp(X, nrank=2, kappa=1e-10, eps=eps, maxit=500, gamma=gamma)

      X[1:10,]
      fit5_0$xhat[1:10,]
      fit5_1$xhat[1:10,]
      fit5_2$xhat[1:10,]
      fit5_3$xhat[1:10,]
      fit5_4$xhat[1:10,]
      fit5_5$xhat[1:10,]

      fit5_0$vhat
      fit5_1$vhat
      fit5_2$vhat
      fit5_3$vhat
      fit5_4$vhat
      fit5_5$vhat

      png.quaternary(X, vhat=fit5_3$vhat)
      png.quaternary(X, vhat=fit5_4$vhat)

      png.utils::png.angle(data$V, fit5_3$vhat[,1:2])
      #

      fit5_0$fit.path[[r]]$crit.path %>% plot(type="l")
      fit5_1$fit.path[[r]]$crit.path %>% plot(type="l")
      fit5_2$fit.path[[r]]$crit.path %>% plot(type="l")
      fit5_3$fit.path[[r]]$crit.path %>% plot(type="l")
      fit5_4$fit.path[[r]]$crit.path %>% plot(type="l")
      fit5_5$fit.path[[r]]$crit.path %>% plot(type="l")
      #
      norm(fit5_1$xhat - fit5_2$xhat)/length(fit5_2$xhat)
      norm(X - fit5_0$xhat)/length(fit5_2$xhat)
      norm(X - fit5_1$xhat)/length(fit5_2$xhat)
      norm(X - fit5_2$xhat)/length(fit5_2$xhat)
      norm(X - fit5_3$xhat)/length(fit5_2$xhat)
      norm(X - fit5_4$xhat)/length(fit5_2$xhat)
      norm(X - fit5_5$xhat)/length(fit5_2$xhat)
      #


      list(fit1_a, fit1_b, fit1_c, fit2, fit3, fit4, fit5) %>%
        lapply(function(fit) sum((data$X2-fit$xhat)^2)/(n*p) ) %>%
        Reduce("c", .)

      list(fit2, fit3, fit4, fit5) %>%
        lapply(function(fit) png.utils::png.angle(data$V, fit$vhat)$max ) %>%
        Reduce("c", .)
    }



    result_eta0 <- function(eta0){
      # Use mcmapply function
      results_vector <- mcmapply(function(n, p, r, snr, seed) sim1(n=n, p=p, r=r, snr=snr, eta=eta0/log(p), seed=seed),
                                 n=grid$n.seq, seed=grid$seed.seq,
                                 p=100, r=5, snr=5,
                                 mc.cores = numCores)

      # Reshape the results_vector to a matrix form
      n.out <- length(results_vector) / prod(sapply(param.list,length))
      results_array <- array(results_vector,
                             dim=c(n.out, sapply(param.list,length)),
                             dimnames=append(list(out=paste0("V", 1:n.out)), param.list) )

      apply(results_array,1:2,mean)
    }

    res_m0.15 <- result_eta0(eta0=-0.1)
    res_0.0 <- result_eta0(eta0=0.0)
    res_0.1 <- result_eta0(eta0=0.1)
    res_0.2 <- result_eta0(eta0=0.2)


  }

  n.seq <- c(10,100,1000,5000)
  eta.seq <- c(-0.1, 0, 0.1)
  seed.seq <- 1:10

  out <- array(NA, dim=c(4,3,2,3,10), dimnames=list( paste0("n=",c(10,100,1000,5000)), c("V vs PC", "V vs GPC", "PC vs GPC"), c("angle", "fnorm"), paste0("eta=", c(-0.1, 0, 0.1)), paste0("seed=", seed.seq) ))

  for( k in 1:length(seed.seq) ){
    seed <- seed.seq[k]
    for( j in 1:length(eta.seq) ){
      eta=eta.seq[j]

      for( i in 1:length(n.seq) ){
        n=n.seq[i]

        p=10; r=2; d=1; d0=0.1; snr=5; seed=seed
        data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=d,d0=d0,seed=seed,eta=eta )
        X <- data$X2
        fit <- png.gppca_qp(X, nrank=r, kappa=1e-4)

        out[i,1,1,j,k] <- png.utils::png.angle(data$V, prcomp(X)$rot[,1:2])$max
        out[i,2,1,j,k] <- png.utils::png.angle(data$V, fit$vhat)$max
        out[i,3,1,j,k] <- png.utils::png.angle(prcomp(X)$rot[,1:2], fit$vhat)$max

        out[i,1,2,j,k] <- Fnorm_V(data$V, prcomp(X)$rot[,1:2])
        out[i,2,2,j,k] <- Fnorm_V(data$V, fit$vhat)
        out[i,3,2,j,k] <- Fnorm_V(prcomp(X)$rot[,1:2], fit$vhat)
      }
    }
  }

}








function(){
  {
    #1: eta=0
    eta=0

    n=200; p=4; r=2; d=10; d0=0.1; snr=5; seed=2
    X <- sim.simplex(n=n,p=p,r=r,snr=snr,d=d,d0=d0,seed=seed,eta=eta)$X2 #1, 5, 10, 17

    delta=1e-6
    f <- png.ZeroReplace.simple
    X.zero_replaced <- t(apply(X, 1, f, delta=delta))
    X.zero_replaced_logratio <- t(apply(X.zero_replaced,1,png.clr))
    fit.pca <- png.pca( X.zero_replaced_logratio, nrank=2 )

    png.quaternary3d(X, xhat=fit.pca$xhat %>% {t(apply(.,1,png.iclr))})
  }

  {
    #2: eta=0.2
    eta=-0.2

    n=200; p=4; r=2; d=10; d0=0.01; snr=100; seed=2
    X <- sim.simplex(n=n,p=p,r=r,snr=snr,d=d,d0=d0,seed=seed,eta=eta)$X2 #1, 5, 10, 17
    X.zero_replaced <- t(apply(X, 1, f, delta=delta))

    png.quaternary3d(X.zero_replaced, vhat=png.lrpca(X.zero_replaced)$vhat)
    png.quaternary3d(X.zero_replaced, vhat=png.gppca_qp(X.zero_replaced)$vhat)

    png.ternary(X.zero_replaced, vhat=png.lrpca(X.zero_replaced)$vhat)
    png.ternary(X.zero_replaced, vhat=png.gppca_qp(X.zero_replaced)$vhat)

    norm(X - png.lrpca(X.zero_replaced)$xhat, "F")
    norm(X - png.gppca_qp(X.zero_replaced)$xhat, "F")


    #
    #
    #


    delta=1e-6
    f <- png.ZeroReplace.simple
    fit.pca0 <- png.pca( X, nrank=2 )
    fit.lrpca <- png.lrpca( X.zero_replaced, nrank=2 )
    fit.gppca_qp <- png.gppca_qp( X, nrank=2 )
    fit.ppca_qp <- png.ppca_qp( X, nrank=2 )


    #

    A=runif(100,0,1); B=runif(100,0,1);
    C=exp(log(A)+0.2*log(B)+0.1)
    X=cbind(A,B,C) %>% apply(1,function(x) x/sum(x)) %>% t
    png.ternary(X, vhat=png.lrpca(X)$vhat)

    log()

    library(ToolsForCoDa)
    data(bentonites)
    Ben <- bentonites[,1:8]
    Ben.com <- Ben/rowSums(Ben)
    out.lrpca <- lrpca(Ben.com)
    plot(out.lrpca)

    #

    png.quaternary3d(X.zero_replaced, vhat=fit.lrpca$vhat)
    png.quaternary(X, vhat=fit.lrpca$vhat)
    #
    #

    list(X=X, xhat=fit.pca0$xhat) %>% { norm(.$X-.$xhat,"F")/length(.$xhat) }
    list(X=X, xhat=fit.pca0$xhat) %>% { norm(.$X-.$xhat,"F")/length(.$xhat) }
    list(X=X, xhat=fit.pca0$xhat) %>% { norm(.$X-.$xhat,"F")/length(.$xhat) }

    png.quaternary3d(X, xhat=fit.lrpca$xhat)
    png.quaternary3d(X, xhat=fit.pca0$xhat)
    png.quaternary3d(X, xhat=fit.gppca_qp$xhat)
    png.quaternary3d(X, xhat=fit.ppca_qp$xhat)

    fit.pca0$vhat
    fit.gppca_qp$vhat
    fit.ppca_qp$vhat

    png.ternary(X, xhat=fit.pca0$xhat)


  }

}




function(){
  # 1.zero-replacement methods
  # criteria
  #


  #

  png.quaternary(X.zero_replaced, vhat=fit.pca$vhat)
  png.quaternary(X.zero_replaced, vhat=png.pca(X.zero_replaced)$vhat)
  fit.pca$vhat
#

}

{
  # Install the necessary packages if not already installed
  if (!require(compositions)) install.packages("compositions")
  if (!require(FactoMineR)) install.packages("FactoMineR")

  # Load the required libraries
  library(compositions)
  library(FactoMineR)

  # Create an example matrix with compositional data
  # Usually, your data would come from a dataset
  data <- matrix(runif(30), nrow=10, ncol=3) %>% {t(apply(.,1,function(x) x/sum(x)))}
  data[1,] <- (c(1,1,0)+1e-10) %>% {./sum(.)}
  clr(data)

  #

  data[1,]
  clr_data[1,]
  t(apply(data,1,png.clr))[1,]


  # Apply PCA on clr-transformed data
  pca_result <- PCA(clr_data)

  # Print summary of PCA result
  summary(pca_result)



  png.pca(data,3)$uhat %>% apply(2,var) %>% {cumsum(.)/sum(.)}
  summary(prcomp(data))

  ilr(data[1,])

}




function(){



}

x = c(0,0.1,0.8,0.1,0);  delta=0.01

png.ZeroReplace.simple(x, delta=0.01)
png.ZeroReplace.additive(x, delta=0.01)
png.ZeroReplace.multiplicative(x, delta=0.01)


