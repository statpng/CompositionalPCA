#' @importFrom mnormt rmnorm
#' @export png.example.ppca
png.example.ppca <- function(X, newx=X[1,]){
  
  library(Ternary)
  
  x <- newx # X[idx,]
  mu <- colMeans(X)
  r <- 2
  vhat <- prcomp(X)$rotation[,1:r]
  
  png.ternary.init(X)
  png.ternary.PC_axis(mu,vhat,color="darkblue")
  
  
  uhat_total <- NULL
  for( j in 1:r ){
    vj <- vhat[,j]
    uhat <- onedimconvexprojection(mu, x, v=vj)
    x.pc1 <- mu + crossprod(x-mu, vj) %*% vj
    x.cpc1 <- mu + uhat*vj
    AddToTernary(text, x, "o", cex = 2, font = 2, col="purple")
    AddToTernary(text, x.pc1, "x", cex = 1.5, font = 2, col="red")
    AddToTernary(text, x.cpc1, "+", cex = 1.5, font = 2, col="blue")
    
    uhat_total <- c(uhat_total, uhat)
  }
  
  xhat <- mu + tcrossprod( uhat_total, vhat[,1:r] )
  xhat.pc <- mu + crossprod(x-mu,vhat[,1:r])%*%t(vhat[,1:r])
  AddToTernary(text, xhat.pc, "x", cex = 1.5, font = 2, col="red")
  AddToTernary(text, xhat, "+", cex = 1.5, font = 2, col="blue")
  
  text(x = 0.0, y = 1.0, "o: data points;  x: PCA;  +: Projected PCA", col = "black", cex=0.9)
  
  AddToTernary(text, X, "o", cex = 0.8, font = 2, col="grey50")
  
  
}











#' @export png.example.gppca
png.example.gppca <- function(X, newx=X[1,]){
  
  if(FALSE){
    set.seed(2)
    n=50; p=3; r=2
    X <- sim.simplex2(n=n, p=p, r=r, delta=0.1)$X
    newx=c(0.35,-0.05,0.7)
  }
  
  library(Ternary)
  
  x <- newx # X[idx,]
  mu <- colMeans(X)
  vhat <- prcomp(X)$rotation %>% {.[,1:(ncol(.)-1)]}
  uhat <- multidimconvexprojection(mu, x, vhat)
  xhat <- mu + tcrossprod(uhat,vhat)
  
  png.ternary.init(X)
  png.ternary.PC_axis(mu,vhat,color="darkblue")
  
  x.pc <- mu + crossprod(x-mu,vhat)%*%t(vhat)
  x.cpc <- mu + tcrossprod(uhat,vhat)
  
  AddToTernary(text, x, "o", cex = 2, font = 2, col="purple")
  AddToTernary(text, x.pc, "x", cex = 1.5, font = 1.5, col="red")
  AddToTernary(text, x.cpc, "+", cex = 1.5, font = 1.5, col="blue")
  
  
  for( j in 1:ncol(vhat) ){
    v <- vhat[,j]
    x.pc.uni <- mu + crossprod(x-mu,vhat[,j])%*%vhat[,j]
    x.cpc.uni <- mu + tcrossprod(uhat[j],vhat[,j])
    
    png.ternary.point(x.pc.uni, "x", cex=1.5, color="red")
    png.ternary.point(x.cpc.uni, "+", cex=1.5, color="blue")
  }
  
  
  text(x = 0.0, y = 1.0, "o: data points;  x: PCA;  +: Projected PCA", col = "black", cex=0.9)
  
  
  
}








#' @export png.example.compare.gppca
png.example.compare.gppca <- function(X){
  if(FALSE){
    set.seed(2)
    n=50; p=3; r=2
    X <- sim.simplex2(n=n, p=p, r=r, delta=0.1)$X
  }
  
  n=nrow(X); p=ncol(X); r=2
  
  library(Ternary)
  png.ternary.init(X)
  
  x <- X
  mu <- colMeans(X)
  vhat <- prcomp(X)$rotation[,1:r]
  uhat <- cbind(apply(x,1,function(xx)onedimconvexprojection(mu, xx, vhat[,1])),
                apply(x,1,function(xx)onedimconvexprojection(mu, xx, vhat[,2])))
  uhat <- apply(x,1,function(xx)multidimconvexprojection(mu, xx, vhat)) %>% {matrix(.,n,r,byrow=T)}
  xhat <- mu + tcrossprod(uhat,vhat)
  
  png.ternary.PC_axis(mu,vhat,color="darkblue")
  
  x.pc <- mu + tcrossprod(uhat,vhat)
  x.cpc <- mu + tcrossprod(uhat,vhat)
  
  AddToTernary(text, x, "o", cex = 2, font = 2, col="purple")
  AddToTernary(text, x.pc, "x", cex = 1.5, font = 1.5, col="red")
  AddToTernary(text, x.cpc, "+", cex = 1.5, font = 1.5, col="blue")
  
  
  for( j in 1:ncol(vhat) ){
    v <- vhat[,j]
    x.pc.uni <- mu + crossprod(x-mu,vhat[,j])%*%vhat[,j]
    x.cpc.uni <- mu + tcrossprod(uhat[j],vhat[,j])
    
    AddToTernary(text, x.pc.uni, "x", cex = 1.5, font = 2, col="red")
    AddToTernary(text, x.cpc.uni, "+", cex = 1.5, font = 2, col="blue")
    
  }
  
  text(x = 0.0, y = 1.0, "o: data points;  x: PCA;  +: Projected PCA", col = "black", cex=0.9)
  
  
}







#' @export png.example_SQP
png.example_SQP <- function(X, sample=FALSE){
  library(Ternary)
  
  if(FALSE){
    set.seed(2)
    n=50; p=3; r=2
    X <- sim.simplex2(n=n, p=p, r=r, delta=1.0)$X2
  }
  
  
  fit.qppca <- png.ppca_qp(X, eta=1e-4)
  mu <- fit.qppca$mu
  uhat <- fit.qppca$uhat
  vhat <- fit.qppca$vhat
  xhat <- fit.qppca$xhat
  r <- ncol(uhat)
  
  png.ternary.init(X)
  
  # QP
  png.ternary.PC_axis(mu,vhat)
  
  for( k in 1:r ){
    xhat <- tcrossprod(rep(1,nrow(X)),mu) + tcrossprod(uhat[,k],vhat[,k])
    AddToTernary(text, xhat, "+", cex = 0.8, font = 2, col="blue")
  }
  
  xhat <- tcrossprod(rep(1,nrow(X)),mu) + tcrossprod(uhat[,k],vhat[,k])
  AddToTernary(text, fit.qppca$xhat, "*", cex = 1.0, font = 2, col="blue")
  
  
  
  # Ordinary
  vhat <- prcomp(X)$rot[,-ncol(X)]
  png.ternary.PC_axis(mu,vhat,color="darkred")
  for( k in 1:r ){
    xhat <- tcrossprod(rep(1,n),mu) + with(prcomp(X), tcrossprod(x[,k],rotation[,k]))
    AddToTernary(text, xhat, "+", cex = 0.8, font = 2, col="red")
  }
  
  legend("topright", 
         legend = c("PCA", "Sequential QP"),
         cex = 1.0, bty = "n", pch = "+", #pt.cex = 1.0,
         col=c("red", "blue")
  )
  
}









#' @export png.example_SQP2
png.example_SQP2 <- function(X, newx=NULL, sample=FALSE){
  library(Ternary)
  
  if(FALSE){
    n=10; p=3; r=2; seed=2
    X <- sim.simplex2(n=n, p=p, r=r, seed=seed)$X2
  }
  
  
  fit.qppca <- png.ppca_qp(X, eta=1e-4)
  mu <- fit.qppca$mu
  uhat <- fit.qppca$uhat
  vhat <- fit.qppca$vhat
  xhat <- fit.qppca$xhat
  
  
  png.ternary.init(X)
  # QP
  GetStartEnd <- png.ternary.PC_axis(mu,vhat,color="grey20")
  # QP 2
  {
    boundary1 <- (mu+GetStartEnd$Start[1]*vhat[,1])
    point1 <- boundary1 + c(0.1,0,-0.1)
    point2 <- (((png.ternary.ternary2coord(point1)-png.ternary.ternary2coord(boundary1)) %*% png.rotmat(320)) + png.ternary.ternary2coord(boundary1)) %>% {png.ternary.coord2ternary(.) } %>% c
    point3 <- point1 + 0.05*vhat[,1]
    
    X1 <- rbind(X, point1);  fit.qppca1 <- png.ppca_qp(X1, eta=1e-4)
    X2 <- rbind(X, point2);  fit.qppca2 <- png.ppca_qp(X2, eta=1e-4)
    X3 <- rbind(X, point3);  fit.qppca3 <- png.ppca_qp(X3, eta=1e-4)
    
    png.ternary.point(boundary1, "x", cex=1.0, color="black")
    png.ternary.point(point1, "x", cex=1.0, color="red")
    png.ternary.point(point2, "x", cex=1.0, color="blue")
    png.ternary.point(point3, "x", cex=1.0, color="green")
    
    png.ternary.PC_axis(colMeans(X1), fit.qppca1$vhat, color="red")
    png.ternary.PC_axis(colMeans(X2), fit.qppca2$vhat, color="blue")
    png.ternary.PC_axis(colMeans(X3), fit.qppca3$vhat, color="green")
  }
  
}






#' @export png.example_NoSolutionInHDLSS
png.example_NoSolutionInHDLSS <- function(){
  library(Ternary)
  library(dplyr)
  
  if(FALSE){
    set.seed(2)
    n=4; p=3; r=3
    
    X <- sim.simplex2(n=n, p=p, r=r, snr=1)$X2
    
    png.quaternary3d(X)
    png.gppca_qp(X, nrank=3)
  }
  
  
  fit.qppca <- png.ppca_qp(X, eta=1e-4)
  mu <- fit.qppca$mu
  uhat <- fit.qppca$uhat
  vhat <- fit.qppca$vhat
  xhat <- fit.qppca$xhat
  
  
  png.ternary.init(X)
  # QP
  GetStartEnd <- png.ternary.PC_axis(mu,vhat,color="grey20")
  # QP 2
  {
    boundary1 <- (mu+GetStartEnd$Start[1]*vhat[,1])
    point1 <- boundary1 + c(0.1,0,-0.1)
    point2 <- (((png.ternary.png.ternary.ternary2coord(point1)-png.ternary.ternary2coord(boundary1)) %*% png.rotmat(320)) + png.ternary.ternary2coord(boundary1)) %>% {png.ternary.coord2ternary(.) } %>% c
    point3 <- point1 + 0.05*vhat[,1]
    
    X1 <- rbind(X, point1);  fit.qppca1 <- png.ppca_qp(X1, eta=1e-4)
    X2 <- rbind(X, point2);  fit.qppca2 <- png.ppca_qp(X2, eta=1e-4)
    X3 <- rbind(X, point3);  fit.qppca3 <- png.ppca_qp(X3, eta=1e-4)
    
    png.ternary.point(boundary1, "x", cex=1.0, color="black")
    png.ternary.point(point1, "x", cex=1.0, color="red")
    png.ternary.point(point2, "x", cex=1.0, color="blue")
    png.ternary.point(point3, "x", cex=1.0, color="green")
    
    png.ternary.PC_axis(colMeans(X1), fit.qppca1$vhat, color="red")
    png.ternary.PC_axis(colMeans(X2), fit.qppca2$vhat, color="blue")
    png.ternary.PC_axis(colMeans(X3), fit.qppca3$vhat, color="green")
  }
  
}

