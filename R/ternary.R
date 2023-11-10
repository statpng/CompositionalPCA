#' @export png.loading2StartEnd
png.loading2StartEnd <- function(mu,vhat){
  dir_list <- NULL
  for( k in 1:ncol(vhat) ){
    v <- vhat[,k]
    dir_list[[k]] <- rbind(start=mu+onedimconvexprojection(mu, -10*v, v)*v,
                           end=mu+onedimconvexprojection(mu, +10*v, v)*v )
  }
  dir_list
}




#' @export png.ternary.init
png.ternary.init <- function(X, cex=0.8, cex.axis=1.0, color="grey50"){
  
  # old
  # {
  #   par(mar = rep(0.5, 4))
  #   TernaryPlot("A","B","C", grid.lines = 5, grid.lty = "dotted",
  #               grid.minor.lines = 1, grid.minor.lty = "dotted")
  #   AddToTernary(text, X, "o", cex = cex, font = 2, col=color)
  # }
  
  C1 <- colnames(X)[1]
  C2 <- colnames(X)[2]
  C3 <- colnames(X)[3]
  
  par(mar = rep(0.5, 4))
  TernaryPlot(C1,C2,C3, grid.lines = 5, grid.lty = "dotted",
              lab.offset=c(0.16,0.2)[2],
              grid.minor.lines = 1, grid.minor.lty = "dotted", lab.cex=cex.axis, axis.cex=cex.axis)
  AddToTernary(text, X, "o", cex = cex, font = 2, col=color)
}




#' @export png.ternary.PC_axis
png.ternary.PC_axis <- function(mu, vhat, color="darkblue", cex=1, print.PCtext=TRUE, adjust.pc=1.1){
  
  Start <- apply(vhat,2,function(v) Solve_U_SP(v*-10, mu, v, gamma=0))
  End <- apply(vhat,2,function(v) Solve_U_SP(v*10, mu, v, gamma=0))
  for( k in 1:length(Start) ){
    TernaryLines(rbind((mu+Start[k]*vhat[,k]), (mu+End[k]*vhat[,k])), col = color)
    if(print.PCtext){
      TernaryText((mu+(End[k]*adjust.pc)*vhat[,k]), paste0("PC",k), cex = cex, col = color, font = 2)
    }
    
  }
  
  return(list(Start=Start, End=End))
}

#' @export png.ternary.point
png.ternary.point <- function(tuple, shape="+", cex, color="red"){
  AddToTernary(text, tuple, shape, cex = cex, font = 1, col=color)
}


#' @export png.ternary
png.ternary <- function(X, vhat=NULL, xhat=NULL, mu=NULL, cex=1.0, cex.axis=1.0, adjust.pc=1.1, print.caption=TRUE, print.PCtext=FALSE, color.pc="darkblue"){
  if(FALSE){
    cex=0.5
  }
  
  
  png.ternary.init(X, color="grey40", cex=cex*0.6, cex.axis=cex.axis)
  
  if(!is.null(vhat)){
    
    if(is.list(vhat)){
      
      if(is.null(mu)){
        mu <- (colMeans(log(X)))
      }
      
      for( jj in 1:NCOL(vhat[[1]]) ){
        for( ii in 1:(length(vhat)-1) ){
          vhat.i1 = vhat[[ii]][,jj]
          vhat.i2 = vhat[[ii+1]][,jj]
          
          vhat.start <- png.iclr( mu+vhat.i1 )
          vhat.end <- png.iclr( mu+vhat.i2 )
          
          
          TernaryLines(rbind(vhat.start,vhat.end), col = color.pc)
        }
        if(print.PCtext){
          TernaryText(vhat.end*adjust.pc, paste0("PC",jj), cex = cex*0.8, col = color.pc, font = 2)
        }
        
      }
      
    } else {
      
      if(is.null(mu)){
        mu <- colMeans(X)
      }
      
      png.ternary.PC_axis(mu, vhat, color=color.pc, cex=cex*0.8, print.PCtext=print.PCtext, adjust.pc=adjust.pc)
      
    }
  }
  
  if(!is.null(xhat)){
    AddToTernary(text, xhat, "+", font = 1.5, col="darkred", cex=cex*0.8)
    
    if(print.caption){
      text(x = 0, y = -0.15, "o: data points;  x: fitted values", col = "black", cex=cex*0.8)
    }
    
  }
  
  
}




# Unused -----


#' @export png.ternary.ternary2coord
png.ternary.ternary2coord <- function(x){
  Ternary::TernaryToXY(x)
  # a=x[1]; b=x[2]; c=x[3]
  # c(1/2*(2*b+c)/(a+b+c), sqrt(3)/2*(c)/(a+b+c))
}


#' @export png.ternary.coord2ternary
png.ternary.coord2ternary <- function(x){
  Ternary::XYToTernary(x[1], x[2])
}


#' @export png.rotmat
png.rotmat <- function(angle){
  theta = angle * 180/pi
  matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
}