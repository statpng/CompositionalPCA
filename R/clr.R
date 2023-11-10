#' @export png.iclr
png.iclr <- function(x){
  exp(x)/sum(exp(x))
}

#' @export png.clr
png.clr <- function(x){
  gx <- sum(log(x))/length(x)
  log(x)-gx
}

#' @export png.ICLR
png.ICLR <- function(X){
  t(apply(X,1,png.iclr))
}

#' @export png.CLR
png.CLR <- function(X){
  t(apply(X,1,png.clr))
}


