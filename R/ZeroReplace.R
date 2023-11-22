




#' @export ZeroReplace.simple
ZeroReplace.simple <- function(x, delta, zero.eps=1e-12){
  if(length(x) != length(delta)) stop("length(x) != length(delta)")
  x <- ifelse( abs(x)<zero.eps, 0, x)

  for( j in which(x==0) ){
    x[j] <- delta[j]
  }

  x/sum(x)
}


#' @export ZeroReplace.additive
ZeroReplace.additive <- function(x, delta, zero.eps=1e-12){
  x <- ifelse( abs(x)<zero.eps, 0, x)

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


#' @export ZeroReplace.multiplicative
ZeroReplace.multiplicative <- function(x, delta, zero.eps=1e-12){
  x <- ifelse( abs(x)<zero.eps, 0, x)
  
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



