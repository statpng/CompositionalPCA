# #' @export png.lrpca
# png.lrpca <- function(X, nrank=2, zero.replace=NULL, delta="min0"){
#   
#   if(grepl("min",delta)){
#     delta_num <- as.numeric(gsub("min","",delta))
#     
#     delta <- min(X[X>0]) * 10^(delta_num)
#   }
#   
#   
#   n=nrow(X); p=ncol(X); 
#   
#   
#   if(!is.null(zero.replace)){
#     f <- switch(zero.replace, 
#                 "simple"=png.ZeroReplace.simple,
#                 "additive"=png.ZeroReplace.additive,
#                 "multiplicative"=png.ZeroReplace.multiplicative)
#     Xnew <- t(apply(X, 1, f, delta=delta))
#   } else {
#     Xnew <- X
#   }
#   
#   mu = png.iclr( colMeans( log(Xnew) ) )
#   
#   Xclr <- t(apply(Xnew,1,png.clr))
#   fit <- png.pca(Xclr, nrank=nrank)
#   uhat <- fit$uhat
#   vhat <- lapply( seq(-20,20,0.5), function(z){
#     z*fit$vhat
#   })
#   xhat <- fit$xhat %>% {t(apply(.,1,png.iclr))}
#   
#   return( list(mu=mu, uhat=uhat, vhat=vhat, logmu=colMeans( log(Xnew) ), logvhat=fit$vhat, xhat=xhat, X=X, Xnew=Xnew, method="lrpca", zero.replace=zero.replace, delta=delta) )
# }


#' @export png.lrpca
png.lrpca <- function(X, nrank=2, zero.replace=NULL, delta="min0", delta.multiple=FALSE){

  if(grepl("min",delta)){
    delta_num <- as.numeric(gsub("min","",delta))

    if( delta.multiple ){
      delta <- apply(X,2,function(x) ifelse(all(x>0),1,min(x[x>0]) * 10^(delta_num)))
    } else {
      delta <- rep( min(X[X>0]) * 10^(delta_num), ncol(X) )
    }
    
  }


  n=nrow(X); p=ncol(X);

  Xnew <- X
  if(!is.null(zero.replace)){
    for(i in 1:nrow(X)){
      x <- X[i,]

      f <- switch(zero.replace,
                  "simple"=png.ZeroReplace.simple,
                  "additive"=png.ZeroReplace.additive,
                  "multiplicative"=png.ZeroReplace.multiplicative)

      Xnew[i,] <- f(x, delta)
    }
  }


  mu = png.iclr( colMeans( log(Xnew) ) )

  Xclr <- t(apply(Xnew,1,png.clr))
  fit <- png.pca(Xclr, nrank=nrank)
  uhat <- fit$uhat
  vhat <- lapply( seq(-20,20,0.5), function(z){
    z*fit$vhat
  })
  xhat <- fit$xhat %>% {t(apply(.,1,png.iclr))}

  return( list(mu=mu, uhat=uhat, vhat=vhat, logmu=colMeans( log(Xnew) ), logvhat=fit$vhat, xhat=xhat, X=X, Xnew=Xnew, method="lrpca", zero.replace=zero.replace, delta=delta) )
}





# #' @export png.ZeroReplace.simple
# png.ZeroReplace.simple <- function(x, delta=1/2*min(x[x>0])){
# x <- ifelse( abs(x)<1e-12, 0, x)
# x[x==0] <- delta
# x/sum(x)
# }


#' @export png.ZeroReplace.simple
png.ZeroReplace.simple <- function(x, delta){
  if(length(x) != length(delta)) stop("length(x) != length(delta)")
  x <- ifelse( abs(x)<1e-12, 0, x)

  for( j in which(x==0) ){
    x[j] <- delta[j]
  }

  x/sum(x)
}


#' @export png.ZeroReplace.additive
png.ZeroReplace.additive <- function(x, delta){
  x <- ifelse( abs(x)<1e-12, 0, x)

  D=length(x);  Z=sum(abs(x==0))
  
  for( j in 1:length(x) ){
    
    if(x[j]==0){
      x[j] = delta[j]*(Z+1)*(D-Z) / D^2
    } else {
      if( any(x[j] < delta*(Z+1)*Z / D^2) ) stop("Reduce the delta value!")
      x[j] = x[j] - delta[j]*(Z+1)*Z / D^2
    }
  }
  
  return(x)
}


#' @export png.ZeroReplace.multiplicative
png.ZeroReplace.multiplicative <- function(x, delta){
  x <- ifelse( abs(x)<1e-12, 0, x)
  
  c=sum(x); Z=sum(x==0)
  
  for( j in 1:length(x) ){
    
    if(x[j]==0){
      x[j] <- delta[j]
    } else {
      x[j] <- (1-delta[j]*Z/c)*x[j]
    }
    
  }
  

  return(x)
}



