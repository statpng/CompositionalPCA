#' @export reconstruction_error
reconstruction_error <- function(fit, data, n.test="10x"){
  
  if("type.projection" %in% names(fit)){
    fit.list <- list(fit)
  } else {
    fit.list <- fit
  }
  
  
  X0 <- data$X0
  params <- data$params
  
  if(n.test == "10x"){
    params_test <- png.utils::png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.utils::png.list.replace(params, list(n=n.test))
  }
  
  if( data$type == "Linear" ){
    Xtest <- sim.Linear.test(params_test)$X2
  } else if( data$type == "LogNormal" ) {
    Xtest <- sim.LogNormal.test(params_test)$X2
  }
  
  
  rmse.list <- try(sapply(fit.list, function(fit){
    sqrt(mean((X0 - fit$xhat)^2))
  }))
  
  rmspe.list <- try(sapply(fit.list, function(fit){
    xhat_test <- projection(Xtest, fit)$xhat
    # mean(rowSums((Xtest - xhat_test)^2))
    sqrt(mean((Xtest - xhat_test)^2))
  }))
  
  cbind("within-sample"=rmse.list, "out-of-sample"=rmspe.list)
}




#' @export pca.rmse.linear
pca.rmse.linear <- function(fit.list, data, n.test="10x", relative=FALSE){
  X0 <- data$X0
  params <- data$params
  
  if(n.test == "10x"){
    params_test <- png.utils::png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.utils::png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.Linear.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    if(relative){
      sqrt(1/length(X0)) * norm(X0 - fit$xhat, "F") / norm(X0, "F")
    } else {
      sqrt(mean((X0 - fit$xhat)^2))
    }
  }))
  
  rmspe.list <- try(sapply(fit.list, function(fit){
    xhat_test <- projection(Xtest, fit)$xhat
    # mean(rowSums((Xtest - xhat_test)^2))
    
    if(relative){
      norm(Xtest - xhat_test, "F") / norm(Xtest, "F")
    } else {
      sqrt(mean((Xtest - xhat_test)^2))
    }
    
  }))
  
  cbind(rmse=rmse.list, rmspe=rmspe.list)
}




#' @export pca.rmse.LogNormal
pca.rmse.LogNormal <- function(fit.list, data, n.test="10x", relative=FALSE){
  X0=data$X0
  params=data$params
  
  if(n.test == "10x"){
    params_test <- png.utils::png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.utils::png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.LogNormal.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    if(relative){
      sqrt(1/length(X0)) * norm(X0 - fit$xhat, "F") / norm(X0, "F")
    } else {
      sqrt(mean((X0 - fit$xhat)^2))
    }
  }))
  rmspe.list <- try(sapply(fit.list, function(fit){
    xhat_test <- projection(Xtest, fit)$xhat
    # mean(rowSums((Xtest - xhat_test)^2))
    
    if(relative){
      norm(Xtest - xhat_test, "F") / norm(Xtest, "F")
    } else {
      sqrt(mean((Xtest - xhat_test)^2))
    }
    
  }))
  
  cbind(rmse=rmse.list, rmspe=rmspe.list)
}





#' @export pca.rmse.linear_rank
pca.rmse.linear_rank <- function(fit.list, data, n.test="10x", relative=FALSE){
  X0=data$X0
  params=data$params
  
  if(n.test == "10x"){
    params_test <- png.utils::png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.utils::png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.Linear.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- projection(data$X2, fit=fit, nrank=r)$xhat
      
      if(relative){
        sqrt(1/length(data$X0)) * norm(data$X0 - Xhat, "F") / norm(data$X0, "F")
      } else {
        sqrt(mean((data$X0 - Xhat)^2))
      }
      
    }) %>% unlist
    
  }))
  
  
  # rmse.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(type.projection, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=type.projection)) + utils::ggplot.scale_y_log10()
  
  
  rmspe.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- projection(Xtest, fit=fit, nrank=r)$xhat
      
      if(relative){
        norm(Xtest - Xhat, "F") / norm(Xtest, "F")
      } else {
        sqrt(mean((Xtest - Xhat)^2))
      }
      
    }) %>% unlist
    
  }))
  
  # rmspe.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(type.projection, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=type.projection)) + utils::ggplot.scale_y_log10()
  
  
  list(rmse=rmse.list, rmspe=rmspe.list)
}



#' @export pca.rmse.LogNormal_rank
pca.rmse.LogNormal_rank <- function(fit.list, data, n.test="10x", relative=FALSE){
  X0=data$X0
  params=data$params
  
  if(n.test == "10x"){
    params_test <- png.utils::png.list.replace(params, list(n=10*params[["n"]]))
  } else {
    params_test <- png.utils::png.list.replace(params, list(n=n.test))
  }
  
  Xtest <- sim.LogNormal.test(params_test)$X2
  
  rmse.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- projection(data$X2, fit=fit, nrank=r)$xhat
      
      if(relative){
        sqrt(1/length(data$X0)) * norm(data$X0 - Xhat, "F") / norm(data$X0, "F")
      } else {
        sqrt(mean((data$X0 - Xhat)^2))
      }
      
    }) %>% unlist
    
  }))
  
  
  # rmse.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(type.projection, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=type.projection)) + utils::ggplot.scale_y_log10()
  
  
  rmspe.list <- try(sapply(fit.list, function(fit){
    
    lapply(1:data$params$r, function(r){
      Xhat <- projection(Xtest, fit=fit, nrank=r)$xhat
      
      if(relative){
        norm(Xtest - Xhat, "F") / norm(Xtest, "F")
      } else {
        sqrt(mean((Xtest - Xhat)^2))
      }
      
    }) %>% unlist
    
  }))
  
  # rmspe.list %>% {cbind.data.frame(rank=1:nrow(.), .)} %>% as.data.frame %>% gather(type.projection, value, -rank) %>% ggplot() + geom_line(aes(rank, value, color=type.projection)) + utils::ggplot.scale_y_log10()
  
  
  list(rmse=rmse.list, rmspe=rmspe.list)
}















adist <- function(x,y){
  sqrt(sum((clr(x) - clr(x))^2))
}

STRESS <- function(X){
  if(FALSE){
    set.seed(1); 
    X=rdirichlet(10,c(1,2,3)) %>% {ifelse(.<0.2,0,.) %>% apply(1,function(x)x/sum(x)) %>% t} 
    Y=rdirichlet(10,c(1,2,3)) %>% {ifelse(.<0.2,0,.) %>% apply(1,function(x)x/sum(x)) %>% t} 
  }
  
  MIN <-  min(X[X>0])
  delta.seq <- MIN * pmin(1,10^seq(-2,2,1))
  
  delta <- delta.seq[1]
  
  Xnew <- ZeroReplace.simple(X, delta=delta)
  
  for( i in 1:nrow(X) ){
    adist( X[i,], Xnew[i,] )
  }
  
  
  sqrt(sum((clr(x1) - clr(x2))^2))
}

