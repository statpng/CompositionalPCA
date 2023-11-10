#' @export tmp.simplex
tmp.simplex <- function(n) {
  qr.Q(qr(matrix(1, nrow=n)) ,complete = TRUE)[,-1]
}



#' @export png.quaternary
png.quaternary <- function(X, vhat=NULL, xhat=NULL, mu=NULL, xhat.col="darkred", cex=0.5, cex.text=1.1, pch=1, col="grey50", theta = 675, phi=-15, print.surface=TRUE, alpha=0.05, n.grid=100, use.par=FALSE){
  if(FALSE){
    vhat=NULL; xhat=NULL; xhat.col="darkred"; cex=0.5; theta = 675; phi=-15; print.surface=TRUE; alpha=0.05; n.grid=100
    cex=0.5;  theta = 675;  phi=-15
  }
  
  library(geometry)
  
  # Compute tetrahedron coordinates according to https://mathoverflow.net/a/184585
  tetra <- tmp.simplex(4)
  
  if( is.null(colnames(X)) ){
    colnames(X) <- LETTERS[1:4]
  }
  df3D <- bary2cart(tetra, X)
  vertex <- bary2cart(tetra, diag(1,4,4))
  
  
  
  pm <- par("mfrow")
  # pmar <- par("mar")
  # pmai <- par("mai")
  # pomi <- par("omi")
  
  
  # Plot data
  library(plot3D)
  if(use.par){
    par(mar = rep(0.5, 4), mai=c(0.4,0,0,0), omi=c(0,0,0,0))
  }
  scatter3D(df3D[,1], df3D[,2], df3D[,3], cex=cex,
            xlim = c(-0.55,0.5), 
            ylim = c(-0.55,0.5), 
            zlim = c(-0.55,0.9), 
            col = col, pch = pch, box = FALSE, theta = theta, phi=phi)
  lines3D(tetra[c(1,2,3,4,1,3,1,2,4),1],
          tetra[c(1,2,3,4,1,3,1,2,4),2],
          tetra[c(1,2,3,4,1,3,1,2,4),3],
          col = "grey10", add = TRUE)
  
  
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
          
          dir.xyz <- rbind(vhat.start, vhat.end) %>% bary2cart(tetra, .)
          
          lines3D(dir.xyz[,1],
                  dir.xyz[,2],
                  dir.xyz[,3], col = "blue", add=TRUE)
          
        }
      }
      
    } else {
      
      if(is.null(mu)){
        mu <- colMeans(X)
      }
      
      dir <- png.loading2StartEnd(mu, vhat)
      
      dir.xyz <- NULL
      for( i in 1:length(dir) ){
        dir.xyz[[i]] <- bary2cart(tetra, dir[[i]])
      }
      
      for( i in 1:length(dir.xyz) ){
        lines3D(dir.xyz[[i]][,1],
                dir.xyz[[i]][,2],
                dir.xyz[[i]][,3], col = "blue", add=TRUE)
      }
      
      
      if(print.surface & ncol(vhat)==2){
        
        u.StartEnd <- lapply( png.loading2StartEnd(colMeans(X), vhat), function(x){
          c(-norm(x[1,],"2"),norm(x[2,],"2"))
        } )
        
        n.grid <- n.grid
        u.grid <- as.matrix(expand.grid(seq(u.StartEnd[[1]][1]*1.2,u.StartEnd[[1]][2]*1.2,length.out=n.grid),
                                        seq(u.StartEnd[[2]][1]*1.2,u.StartEnd[[2]][2]*1.2,length.out=n.grid)))
        grid <- (tcrossprod(rep(1,n.grid^2),mu)+tcrossprod(u.grid,vhat))
        wh.outliers <- which( apply(grid,1,function(x) any(x < -1e-5)) )
        
        surface.xyz <- t(apply(grid, 1, function(surf) bary2cart(tetra, surf)))
        surface.xyz[wh.outliers,] <- 0
        colnames(surface.xyz) <- letters[24:26]
        
        
        surf3D(x = matrix(surface.xyz[,1],n.grid,n.grid),
               y = matrix(surface.xyz[,2],n.grid,n.grid),
               z = matrix(surface.xyz[,3],n.grid,n.grid), 
               col="blue", colkey=FALSE, add=TRUE,
               alpha = alpha)
        
      }
      
      
    }
    
    
  }
  
  
  
  
  if(!is.null(xhat)){
    df_xhat <- bary2cart(tetra, xhat)
    
    df_xhat %>% {
      scatter3D(.[,1], .[,2], .[,3], cex=cex,
                xlim = c(-0.55,0.5), 
                ylim = c(-0.55,0.5), 
                zlim = c(-0.55,0.9), 
                col = xhat.col, pch = 3, add=TRUE)
    }
  }
  
  # tetra[1,1] <- tetra[1,1]+0.1
  tetra[1,3] <- tetra[1,3]-0.1
  tetra[2,1] <- tetra[2,1]+0.2
  tetra[3,1] <- tetra[3,1]-0.2
  tetra[4,3] <- tetra[4,3]+0.10
  
  text3D(tetra[,1], tetra[,2], tetra[,3], colnames(X), add = TRUE, cex=cex.text)
  
  # par(mfrow = pm)
  # par(mar = pmar)
  # par(mai = pmai)
  # par(omi = pomi)
  
  # return(df3D)
}




#' @export png.quaternary3d
png.quaternary3d <- function(X, vhat=NULL, xhat=NULL, mu=NULL, xhat.size=2, size = 2, print.surface=TRUE, ...) {
  if(FALSE){
    vhat=NULL; xhat=NULL; mu=NULL;
    xhat.size=2; size = 2; print.surface=TRUE; ...=NULL
    size=2; ...=NULL
  }
  
  library(plotly)
  library(dplyr)
  library(geometry)
  
  
  tetra <- tmp.simplex(4)
  
  df <- bary2cart(tetra, X)
  colnames(df) <- c("x", "y", "z")
  
  p <- plot_ly(as.data.frame(df),
               x = ~x, y = ~y, z = ~z, 
               type = 'scatter3d', mode="markers", alpha=1.0,
               marker = list(symbol = 3, size = size, color = "grey10",
                             line = list(color="grey10", width=1), ... ) ) %>%
    plotly::layout( scene = list(xaxis = list( range=c(-1,1) ),
                         yaxis = list( range=c(-1,1) ),
                         zaxis = list( range=c(-1,1) ))
    ) %>% 
    add_trace(x=tetra[c(1,2,3,4,1,3,1,2,4),1],
              y=tetra[c(1,2,3,4,1,3,1,2,4),2],
              z=tetra[c(1,2,3,4,1,3,1,2,4),3],
              type="scatter3d", mode="lines",
              name="lines", color=I("grey20"), 
              showlegend=F, inherit=F)
  
  
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
          
          dir.xyz <- rbind(vhat.start, vhat.end) %>% bary2cart(tetra, .)
          
          p <- p %>% 
            add_trace(x = dir.xyz[,1],
                      y = dir.xyz[,2],
                      z = dir.xyz[,3],
                      type = "scatter3d", mode = "lines", 
                      name = "lines", line=list(color="darkred"), 
                      showlegend=FALSE, inherit=FALSE)
          
        }
      }
      
    } else {
      
      if(is.null(mu)){
        mu <- colMeans(X)
      }
      
      dir <- png.loading2StartEnd(mu, vhat)
      
      for( i in 1:length(dir) ){
        dir.xyz <- bary2cart(tetra, dir[[i]])
        
        p <- p %>% 
          add_trace(x = dir.xyz[,1],
                    y = dir.xyz[,2],
                    z = dir.xyz[,3],
                    type = "scatter3d", mode = "lines", 
                    name = "lines", line=list(color="darkred"), 
                    showlegend=FALSE, inherit=FALSE)
      }
      
      
      
      if(print.surface & ncol(vhat)==2){
        
        n.grid <- 100
        u.grid <- as.matrix(expand.grid(seq(-1,1,length.out=n.grid),
                                        seq(-1,1,length.out=n.grid)))
        grid <- (tcrossprod(rep(1,n.grid^2),mu)+tcrossprod(u.grid,vhat)) %>% 
          {.[!apply(.,1,function(x) any(x < -1e-5)),]}
        
        surface.xyz <- t(apply(grid, 1, function(surf) bary2cart(tetra, surf)))
        colnames(surface.xyz) <- letters[24:26]
        
        
        p <- p %>%
          add_trace(x = surface.xyz[,1],
                    y = surface.xyz[,2],
                    z = surface.xyz[,3],
                    type = "mesh3d", 
                    inherit=F,
                    colors=c("red"),
                    name = "surface", showlegend=FALSE, 
                    opacity=0.3)
        
      }
    }
  }
  
  
  
  
  if(!is.null(xhat)){
    df_xhat <- bary2cart(tetra, xhat)
    
    p <- p %>% 
      add_trace(x = df_xhat[,1],
                y = df_xhat[,2],
                z = df_xhat[,3],
                type = "scatter3d", mode = "markers", 
                name = "xhat", color=I("darkred"), 
                showlegend = FALSE,
                marker = list(symbol="cross", 
                              size=xhat.size) )
    
  }
  
  
  return(p)
  
}




