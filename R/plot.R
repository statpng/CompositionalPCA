#' @export png.pca.plot
png.pca.plot <- function(X, ...){
  X %>% {
    PCs <- prcomp(.)$x
    ExplainedVariance <- round(apply(PCs,2,var) %>% {(.)/sum(.)} %>% head(2),4)*100
    
    pch=18
    PCs %>% plot(pch=pch, 
                 xlab=paste0("PC1 (", format(ExplainedVariance[1], digits=1, nsmall=1), "%)"),
                 ylab=paste0("PC2 (", format(ExplainedVariance[2], digits=1, nsmall=1), "%)"), ...)
    # mtext(side=3, text=site, adj=0)
    # mtext(side=3, text=paste(paste0(c("PC1: ", "PC2: "), round(PCs %>% apply(2,var) %>% {(.)/sum(.)} %>% head(2),3)*100,"%"), collapse="; "), padj=0, adj=1, cex=0.6)
  }
  
}


#' @export png.pca.convergence
png.pca.convergence <- function(fit){
  if(FALSE){
    fit <- fit5_1
  }
  n=nrow(fit$X); p=ncol(fit$X); r=NCOL(fit$vhat)
  mu=fit$mu
  maxit=fit$maxit
  
  crit.array <- array(NA, 
                      dim=c(maxit-1, r, 3), 
                      dimnames=list(paste0("iter=",1:(maxit-1)), 
                                    paste0("rank=",1:r), 
                                    c("uhat", "vhat", "xhat") ))
  for( k in 1:r ){
    est.path <- fit$fit.path[[k]]$est.path
    
    if( length(est.path) == 1 ){ 
      crit.array[1,k,1] <- 1e-16
      crit.array[1,k,2] <- 1e-16
      crit.array[1,k,3] <- 1e-16
      
      next;
    }
    
    for( kk in 1:(length(est.path)-1) ){
      # kk
      uhat_old <- as.matrix(est.path[[kk]]$uhat)
      vhat_old <- as.matrix(est.path[[kk]]$vhat)
      xhat_old <- tcrossprod(rep(1,n),mu)+tcrossprod(uhat_old,vhat_old)
      
      # kk+1
      uhat_new <- as.matrix(est.path[[kk+1]]$uhat)
      vhat_new <- as.matrix(est.path[[kk+1]]$vhat)
      xhat_new <- tcrossprod(rep(1,n),mu)+tcrossprod(uhat_new,vhat_new)
      
      crit.array[kk,k,1] <- sqrt(sum((uhat_new-uhat_old)^2)/(n*r))
      crit.array[kk,k,2] <- sqrt(sum((vhat_new-vhat_old)^2)/(n*r))
      crit.array[kk,k,3] <- sqrt(sum((xhat_new-xhat_old)^2)/(n*r))
    }
  }
  
  out <- crit.array[,,] %>% plyr::adply(1:3)
  colnames(out) <- c("iter", "rank", "est", "value")
  
  out2 <- subset(out, !is.na(value))
  
  print( out2 %>% group_by(est, rank) %>% slice_tail(n=1) )
  
  return( out2 )
}


#' @export png.pca.plot_convergence
png.pca.plot_convergence <- function(fit, maxit=100){
  if(is.data.frame(fit)){
    df <- fit
  } else {
    df <- png.pca.convergence(fit)
  }
  
  df[,1] <- as.numeric(df[,1])
  
  df %>%
    # filter(est=="xhat") %>% 
    filter(iter>1, iter<=maxit) %>% 
    ggplot() +
    geom_line(aes(iter, value, color=rank)) + 
    ggsci::scale_color_lancet() +
    # scale_y_sqrt() +
    # scale_y_continuous(trans='log2') +
    xlab("Iteration") + ylab("Changes") +
    facet_grid(est~., scales="free") +
    theme_bw()
}













#' @export png.angle
png.angle <- function(true, est){
  # true: n x p; est: n x p
  
  # The largest principal angle
  qt <- qr.Q(qr(true))
  qe <- qr.Q(qr(est))
  fit.svd <- svd( crossprod(qe, qt) )
  theta <- acos(fit.svd$d |> round(12))
  
  # theta[1] * 180 / pi (in degree)
  list( max = theta[1] * 180 / pi, Grassmannian = norm( theta, "2" ) * 180 / pi )
}




#' @export png.pca.criteria
png.pca.criteria <- function(fit, data, n.test){
  if(FALSE){
    n=50; p=50; r=5; snr=2; eta=0.1/log(p); seed=1
    data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta)
    
    fit1 <- png.ppca_qp(data$X2, nrank=r, kappa=1e-6, maxit=2000, eps=1e-6, gamma=0.5, save.est.path = TRUE)
    fit2 <- png.gppca_qp(data$X2, nrank=r, kappa=1e-6, maxit=2000, eps=1e-6, gamma=0.5, save.est.path = TRUE)
    
    fit1 %>% png.pca.criteria(data, n.test=n*5)
    fit2 %>% png.pca.criteria(data, n.test=n*5)
    
  }
  
  if(FALSE){
    fit=fit_ppca
    data
    n.test=200
  }
  

  data$params["n"] <- n.test
  true <- data %>% { list(Xtrain=.$X2,
                          Xtest=sim.simplex.test(.$params)$X2,
                          V=.$V) }
  
  Xtrain=true$Xtrain
  Xtest=true$Xtest
  Vtrue=true$V
  
  vhat <- fit$vhat
  xhat_train <- fit$xhat
  xhat_test <- png.projection(Xtest, fit, method=fit$method)
  
  # png.projection(Xtrain, fit, method=fit$method)[1:5,1:5]
  # xhat_train[1:5,1:5]
  
  obj.Xtrain <- sqrt(mean((Xtrain-xhat_train)^2))
  obj.Xtest <- sqrt(mean((Xtest-xhat_test)^2))
  
  rmse.Xtrain <- sqrt(mean((Xtrain-xhat_train)^2))
  rmse.Xtest <- sqrt(mean((Xtest-xhat_test)^2))
  Pangle.V <- png.angle(Vtrue, vhat)$max
  Gangle.V <- png.angle(Vtrue, vhat)$Grassmannian
  # Out-of-simplex Sample Percentage
  OutOfSimplex <- mean(apply(xhat_train,1,function(x) any(x < -1e-8)))
  Sparsity <- mean( abs(xhat_train) < 1e-12 )
  Convergence <- sapply(fit$fit.path, function(fit_rank) fit_rank$it < fit$maxit)
  # fit$fit.path[[4]]$crit.path[-(1:100)] %>% plot(type="l")
  # fit$fit.path[[4]]$crit.path %>% tail
  
  
  c(rmse.Xtrain=rmse.Xtrain,
    rmse.Xtest=rmse.Xtest,
    Pangle.V=Pangle.V,
    Gangle.V=Gangle.V,
    OutOfSimplex=OutOfSimplex,
    Sparsity=Sparsity,
    Convergence=Convergence)
  
}





#' @export png.CompositionalPlot
png.CompositionalPlot <- function(pseq, xhat=NULL, uhat=NULL, title="", df.color=NULL, legend.ncol=10){
  if(FALSE){
    pseq <- pseq_list_total_xhat$urine$Phylum
    xhat <- fit4$xhat
  }
  
  if(!is.null(xhat)){
    pseq@otu_table <- t(xhat) %>% otu_table(taxa_are_rows = TRUE )
  }
  rownames(pseq@otu_table) <- pseq@tax_table[,"unique"]
  
  n.taxa <- length(taxa(pseq))
  # otu.sort <- "abundance"
  otu.sort <- top_taxa(pseq, n=n.taxa)
  # otu.sort <- rev(c(rev(names(sort(rowSums(abu)))[seq(1, nrow(abu), 2)]), 
  # names(sort(rowSums(abu)))[seq(2, nrow(abu), 2)]))
  
  
  if(is.null(uhat)){
    TopTaxa <- pseq@tax_table[ top_taxa(pseq, n=1), "unique" ]
    sample.sort <- (sample_names(pseq)[order(abundances(pseq)[TopTaxa, 
    ])])
  } else {
    sample.sort <- (sample_names(pseq)[order(uhat)])
  }
  
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% 
    # {.[c("Set1","Pastel1","Dark2","Set3"),]} %>% 
    {.[c("Set1","Set3","Pastel1","Dark2","Pastel2", "Set2", "Accent"),]} %>%
    { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }
  
  # top_taxa(pseq, n=n.taxa)
  # 
  # adjustcolor( "red", alpha.f = 0.2)
  
  # remotes::install_github("KarstensLab/microshades")
  # library(microshades)
  # microshades::create_color_dfs()
  # microshades_palettes
  # library(randomcoloR)
  # set.seed(1)
  # Taxa_cols <- distinctColorPalette(k=n.taxa)
  
  
  
  
  
  abu <- abundances(pseq)
  dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
  names(dfm) <- c("Tax", "Sample", "Abundance")
  
  dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels=otu.sort)
  
  {
    Tax <- Sample <- Abundance <- NULL
    dfm <- dfm %>% arrange(Tax)
    dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
    p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) + 
      geom_bar(position="stack", stat="identity") + 
      scale_x_discrete(labels = dfm$xlabel, 
                       breaks = dfm$Sample)
    p <- p + labs(y = "Abundance")
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                              hjust = 0))
    p <- p + guides(fill = guide_legend(reverse = FALSE))
    # if (!is.null(group_by)) {
    #   p <- p + facet_grid(. ~ Group, drop = TRUE, space = "free", 
    #                       scales = "free")
    # }
  }
  
  
  
  if( is.null(df.color) ){
    p +
      # scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
      scale_fill_manual("", values = Taxa_cols, na.value = "black") +
      # scale_fill_manual("", values = df.color$Color, breaks=df.color$Tax, na.value = "black") +
      scale_y_continuous(label = scales::percent) +
      # hrbrthemes::theme_ipsum(grid="Y") +
      theme_bw(base_size=14) +
      labs(x = "Samples", y = "Relative Abundance",
           subtitle = title) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "bottom") +
      guides(fill = guide_legend(ncol = legend.ncol))
  } else {
    p +
      # scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
      # scale_fill_manual("", values = Taxa_cols, na.value = "black") +
      scale_fill_manual("", values = df.color$Color, breaks=df.color$Tax, na.value = "black") +
      scale_y_continuous(label = scales::percent) +
      # hrbrthemes::theme_ipsum(grid="Y") +
      theme_bw(base_size=14) +
      labs(x = "Samples", y = "Relative Abundance",
           subtitle = title) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "bottom") +
      guides(fill = guide_legend(ncol = legend.ncol))
  }
  
  
}




#' @export png.MultiCompositionalPlot
png.MultiCompositionalPlot <- function(p.list, title="", legend.ncol=10){
  if(FALSE){
    p.list <- list(ppca=p1, gppca=p2, lrpca=p3)
  }
  
  
  p.names <- names(p.list)
  data.list <- lapply(p.list, function(x) x$data)
  
  
  TopTaxa.list <- lapply(data.list, function(x) x %>% group_by(Tax) %>% summarise(avg=mean(Abundance)) %>% arrange(desc(avg)) %>% slice_head(n=1) %>% .$Tax)
  
  
  OrderedSample.list <- lapply(1:length(data.list), function(i) data.list[[i]] %>% filter(Tax == TopTaxa.list[[i]]) %>% arrange(desc(Abundance)) %>% .$Sample)
  
  tmp <- gtools::mixedsort(as.character(unique(OrderedSample.list[[1]])))
  sample.order <- lapply(OrderedSample.list, function(y){
    purrr::map_int( tmp, ~which( y %in% .x ) )
  }) %>% {Reduce("+", .)/length(.)} %>% order
  #
  
  for( i in 1:length(data.list) ){
    data.list[[i]]$Sample <- factor(data.list[[i]]$Sample, levels=tmp[sample.order])
  }
  
  
  dfm <- NULL
  for(i in 1:length(data.list) ){
    dfm <- rbind.data.frame(dfm, cbind.data.frame(Type=p.names[i], data.list[[i]]))
  }
  
  
  n.taxa <- length(unique(p.list[[1]]$data$Tax))
  TopTaxa <- dfm %>% group_by(Tax) %>% summarise(avg=mean(Abundance)) %>% slice_head(n=n.taxa) %>% .$Tax
  otu.sort <- TopTaxa
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  
  dfm$Type <- factor(dfm$Type, levels=p.names)
  dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels=otu.sort)
  
  {
    Tax <- Sample <- Abundance <- NULL
    dfm <- dfm %>% arrange(Tax)
    dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
    p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) + 
      geom_bar(position="stack", stat="identity") + 
      scale_x_discrete(labels = dfm$xlabel, 
                       breaks = dfm$Sample)
    p <- p + labs(y = "Abundance")
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                              hjust = 0))
    p <- p + guides(fill = guide_legend(reverse = FALSE))
    p <- p + facet_grid( Type ~ ., drop = TRUE, space = "free",
                        scales = "free")
    
    p +
      # scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
      scale_fill_manual("", values = Taxa_cols, na.value = "black") +
      scale_y_continuous(label = scales::percent) +
      # hrbrthemes::theme_ipsum(grid="Y") +
      theme_bw(base_size=14) +
      labs(x = "Samples", y = "Relative Abundance",
           subtitle = title) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "bottom") +
      guides(fill = guide_legend(ncol = legend.ncol))
  }
  
  
}



#' @export png.CheckComposition
png.CheckComposition <- function(X, eps=1e-8){
  if(FALSE){
    n=10; p=500
    set.seed(2)
    
    X <- matrix(runif(n*p,0,1),n,p) %>% {t(apply(.,1,function(x)x/sum(x)))}
    X <- png.pca(X,nrank=1)$xhat
    sum(Xhat>1|Xhat<0)
  }
  
  wh <- t( apply(X,1,function(x) {
    a=(abs(sum(x)-1)>eps)
    b=any(x<0-eps)
    c=any(x>1+eps)
    c(a,b,c)
  }) )
  
  out <- apply(wh,2,which)
  
  if(length(out)==0){
    return(list(`=1`=NA,
                `<0`=NA, 
                `>1`=NA, 
                any=NA))
  }
  names(out) <- c("=1", "<0", ">1")
  out$any <- unique( unlist(out) )
  out
}

#' @export png.CompositionalViolation
png.CompositionalViolation <- function(X){
  if(FALSE){
    n=10; p=100
    X <- matrix(runif(n*p,0,1),n,p)
  }
  
  png.CheckComposition( png.pca(X)$xhat )
  
  X
  
}



#' @export png.crit.path
png.crit.path <- function(fit, remove=NULL){
  fit.path <- fit$fit.path
  N <- length(fit.path)
  
  par(mfrow=c(N,1), 
      mai=c(.2,.4,.2,.2), 
      omi=c(.2,.2,.1,.2))
  
  if(is.null(remove)){
    purrr::map( fit.path, ~{
      .x$crit.path %>% plot(type="l")
      .x$crit.path %>% tail
    } )
  } else {
    purrr::map( fit.path, ~{
      .x$crit.path[-remove] %>% plot(type="l")
      .x$crit.path %>% tail
    } )
  }
  
}



































# Multi-source -----------------------------------------------------------------------



#' @export make_barplot2
make_barplot2 <- function (dfm, group_by, direction="vertical") {
  
  Tax <- Sample <- Abundance <- NULL
  
  # Provide barplot
  dfm <- dfm %>% arrange(Tax)  # Show Taxs always in the same order
  dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
  
  p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) +
    geom_bar(position="stack", stat="identity") +
    scale_x_discrete(labels=dfm$xlabel, breaks=dfm$Sample)
  
  # Name appropriately
  p <- p + labs(y = "Abundance")
  
  # Rotate horizontal axis labels, and adjust
  p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5,
                                          hjust=0))
  p <- p + guides(fill=guide_legend(reverse=FALSE))
  
  if (!is.null(group_by)) {
    if( direction == "vertical" ){
      p <- p + facet_grid(Group1~Group2)
    } else if( direction == "horizontal" ) {
      p <- p + facet_grid(Group2~Group1, drop = TRUE,
                          space = "free", scales = "free") 
    }
    
    
  }
  list(df=dfm, plot=p)
}








#' @export png.top_taxa_list
png.top_taxa_list <- function(LIST, cutoff=0.0001, n=NULL){
  
  if( is.null(n) ){
    n <- 0.5*length(taxa(LIST[[1]]))
  }
  
  lapply( LIST, function(x) taxa_sums(x) / ncol(x@otu_table) ) %>% do.call("c", .) %>% sort(decreasing=TRUE) %>% {.[.>cutoff]} %>% {names(.)[!duplicated(names(.))]} %>% head(n)
  
}




#' @export png.PlotComposition
png.PlotComposition <- function(pseq, taxa=NULL, group_by="Site", group.level=c("urine", "serum", "stool", "stoolp", "stools"), sample.pattern="[a-z]+\\_", filename="./plot.pdf", height=5, width=10, legend.ncol=10){
  
  
  # tmp.list[[idx]] %>% { prune_taxa( taxa_sums(.)/N > 0.001, . ) } %>% taxa
  # tmp.list[[idx]] %>% top_taxa(14)
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  TOP_taxa <- pseq@tax_table[top_taxa(pseq, n = 1), "unique"]
  
  fit.comp <- pseq %>% 
    plot_composition(sample.sort = TOP_taxa, otu.sort = "abundance", group_by=group_by, group.level=group.level, direction="vertical", sample.pattern=sample.pattern)
  
  fit.comp$df %>% head %>% print
  
  if( is.null(taxa) ){
    taxa <- levels(fit.comp$df$Tax)
  }
  
  if( length(taxa) > 74 ) warnings("The maximum number of colors is 74.")
  
  
  p <- fit.comp$plot +
    scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
    scale_y_continuous(label = scales::percent) + 
    # hrbrthemes::theme_ipsum(grid="Y") +
    theme_bw(base_size=16) +
    labs(x = "Samples", y = "Relative Abundance",
         title = "Relative Abundance data") + 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          legend.position = "bottom") + 
    guides(fill = guide_legend(ncol = legend.ncol))
  
  
  
  pdf(file=filename, height=height, width=width)
  print(p)
  dev.off()
  
  
  invisible(fit.comp)
  
}









#' @export png.PlotComposition2
png.PlotComposition2 <- function(pseq, taxa=NULL, group_by=c("Structure", "Site"), group.level=list( c("data", "joint", "individual"), c("urine", "serum", "stools") ), filename="./plot.pdf", sample.sort="top", height=5, width=10, legend.ncol=10){
  
  if(FALSE){
    pseq <- decomp_ajive
    level="Phylum"
    filename="./plot.pdf"
    height=5
    width=10
    legend.ncol=10
    sample.sort="top"
  }
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  
  if(sample.sort == "top"){
    TOP_taxa <- pseq@tax_table[top_taxa(pseq, n = 1), "unique"]
  } else {
    TOP_taxa <- NULL
  }
  
  
  
  if( is.null(taxa) ){
    fit.comp0 <- pseq %>% phyloseq::subset_samples(Structure == "data") %>% 
      plot_composition(sample.sort = TOP_taxa, otu.sort = "abundance", group_by=c("Site"), group.level=c("urine", "serum", "stool", "stoolp", "stools"), group.type="vertical", sample.pattern="[a-z]+\\_[a-z]+\\_")
    
    taxa <- levels(fit.comp0$df$Tax)
  }
  
  
  fit.comp <- pseq %>% 
    plot_composition2(sample.sort = TOP_taxa, otu.sort = "abundance", group_by=group_by, group.level=group.level, group.type="vertical", sample.pattern="[a-z]+\\_[a-z]+\\_")
  
  
  fit.comp$df %>% head %>% print
  
  
  if( length(taxa) > 74 ) warnings("The maximum number of colors is 74.")
  
  
  p <- fit.comp$plot +
    scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
    scale_y_continuous(label = scales::percent) + 
    # hrbrthemes::theme_ipsum(grid="Y") +
    theme_bw(base_size=16) +
    labs(x = "Samples", y = "Relative Abundance",
         title = "Relative Abundance data") + 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          legend.position = "bottom") + 
    guides(fill = guide_legend(ncol = legend.ncol))
  
  
  
  pdf(file=filename, height=height, width=width)
  print(p)
  dev.off()
  
  
  invisible(fit.comp)
  
}










#' @export plot_composition2
plot_composition2 <- function(x,
                              sample.sort=NULL,
                              otu.sort=NULL,
                              x.label="sample",
                              plot.type="barplot",
                              verbose=FALSE, 
                              average_by=NULL,
                              group_by = NULL, 
                              group.level = NULL,
                              direction = "vertical", 
                              sample.pattern = NULL, ...) {
  
  if(FALSE){
    x <- decomp_jive
    sample.sort = TOP_taxa
    otu.sort = "abundance"
    group_by = c("Structure", "Site")
    group.level = list( c("data", "joint", "individual"), c("urine", "serum", "stool", "stoolp", "stools") )
    direction = "vertical"
    sample.pattern="[a-z]+\\_"
    
    x.label="sample"
    plot.type="barplot"
    verbose=FALSE
    average_by=NULL
    
  }
  
  
  # Avoid warnings
  Sample <- Abundance <- Taxon <- Group <- Tax <-
    horiz <- value <- scales <- ID <- 
    meta <- OTU <- taxic <- otu.df <- taxmat <-  new.tax<- NULL
  if (!is.null(x@phy_tree)){
    x@phy_tree <- NULL
  }
  
  xorig <- x
  if (verbose) {message("Pick the abundance matrix taxa x samples")}
  abu <- abundances(x)
  if (verbose) {message("Average the samples by group")}
  group <- NULL
  if (!is.null(average_by)) {
    dff <- as.data.frame(t(abu))
    dff$group <- sample_data(x)[[average_by]]
    if (is.numeric(dff$group) || is.character(dff$group)) {
      dff$group <- factor(dff$group, levels=sort(unique(dff$group)))
    }
    # Remove samples with no group info
    dff <- dff %>% filter(!is.na(group))
    dff$group <- droplevels(dff$group)
    # av <- ddply(dff, "group", colwise(mean))
    av <- aggregate(. ~ group, data = dff, mean)
    rownames(av) <- as.character(av$group)
    av$group <- NULL
    abu <- t(av)  # taxa x groups
  }
  
  if (verbose) {message("Sort samples")}
  if (is.null(sample.sort) || sample.sort == "none" ||
      !is.null(average_by)) {
    # No sorting sample.sort <- sample_names(x)
    sample.sort <- colnames(abu)
  } else if (length(sample.sort) == 1 && sample.sort %in% taxa(xorig)) {
    tax <- sample.sort
    sample.sort <- rev(sample_names(x)[order(abundances(x)[tax,])])
  } else if (length(sample.sort) == 1 &&
             sample.sort %in% names(sample_data(x)) && 
             is.null(average_by)) {
    # Sort by metadata field
    sample.sort <- rownames(sample_data(x))[order(
      sample_data(x)[[sample.sort]])]
  } else if (all(sample.sort %in% sample_names(x)) & is.null(average_by)) {
    # Use predefined order
    sample.sort <- sample.sort
  } else if (length(sample.sort) == 1 && sample.sort == "neatmap") {
    sample.sort <- neatsort(x, method="NMDS", distance="bray",
                            target="sites", first=NULL)
  } else if (is.vector(sample.sort) && length(sample.sort) > 1) {
    sample.sort <- sample_names(x)[sample.sort]            
  } else if (!sample.sort %in% names(sample_data(x))) {
    warning(paste("The sample.sort argument", sample.sort,
                  "is not included in sample_data(x). 
            Using original sample ordering."))
    sample.sort <- sample_names(x)
  }
  
  # Sort taxa
  if (is.null(otu.sort) || otu.sort == "none") {
    # No sorting
    otu.sort <- taxa(x)
  } else if (length(otu.sort) == 1 && otu.sort == "abundance2") {
    otu.sort <- rev(c(rev(names(sort(rowSums(abu)))[seq(1, nrow(abu), 2)]),
                      names(sort(rowSums(abu)))[seq(2, nrow(abu), 2)]))
  } else if (length(otu.sort) == 1 && otu.sort == "abundance") {
    otu.sort <- rev(names(sort(rowSums(abu))))
  } else if (length(otu.sort) == 1 && otu.sort %in%
             colnames(tax_table(x))) {
    otu.sort <- rownames(sample_data(x))[order(tax_table(x)[[otu.sort]])]
  } else if (all(otu.sort %in% taxa(x))) {
    # Use predefined order
    otu.sort <- otu.sort
  } else if (length(otu.sort) == 1 && otu.sort == "neatmap") {
    otu.sort <- neatsort(x, method="NMDS", distance="bray",
                         target="species", first=NULL)
  }
  
  # Abundances as data.frame dfm <- psmelt(x)
  dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
  names(dfm) <- c("Tax", "Sample", "Abundance")
  dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels=otu.sort)
  
  if (!is.null(group_by)) {
    if (!is.null(average_by)) {
      dfm$Group1 <- meta(x)[[group_by[1]]][match(as.character(dfm$Sample),
                                                 meta(x)[[average_by]])]
      dfm$Group2 <- meta(x)[[group_by[2]]][match(as.character(dfm$Sample),
                                                 meta(x)[[average_by]])]
    }else{
      dfm$Group1 <- meta(x)[[group_by[1]]][match(as.character(dfm$Sample),
                                                 sample_names(x))]
      dfm$Group2 <- meta(x)[[group_by[2]]][match(as.character(dfm$Sample),
                                                 sample_names(x))]
    }
  }
  
  if(!is.null(group.level)){
    dfm$Group1 <- factor(dfm$Group1, levels=group.level[[1]])
    dfm$Group2 <- factor(dfm$Group2, levels=group.level[[2]])
  }
  
  # SampleIDs for plotting
  if (x.label %in% colnames(sample_data(x)) & is.null(average_by)) {
    
    meta <- sample_data(x)
    dfm$xlabel <- as.vector(unlist(meta[as.character(dfm$Sample), x.label]))
    
    # Sort the levels as in the original metadata
    if (is.factor(meta[, x.label])) {            
      lev <- levels(meta[, x.label])            
    } else {            
      lev <- unique(as.character(unname(unlist(meta[, x.label]))))
    }       
    dfm$xlabel <- factor(dfm$xlabel, levels=lev)        
  } else {       
    dfm$xlabel <- dfm$Sample        
  }
  
  
  if( !is.null(sample.pattern) ){
    dfm$Sample <- gsub(sample.pattern, "", dfm$Sample)
    dfm$xlabel <- gsub(sample.pattern, "", dfm$xlabel)
  }
  
  
  if (verbose) {message("Construct the plots")}   
  if (plot.type == "barplot") {
    p <- make_barplot2(dfm, group_by, direction=direction)
  } else if (plot.type == "heatmap") {
    p <- make_heatmap1(x, otu.sort, sample.sort, verbose)        
  } else if (plot.type == "lineplot") {
    p <- make_lineplot1(dfm) 
  } else {
    stop("plot.type argument not recognized")
  }
  
  p
  
}











#' @export png.PlotComposition3
png.PlotComposition3 <- function(pseq, taxa=NULL, sample.sort="unique", filename="./plot.pdf", height=5, width=10, legend.ncol=10){
  
  # png.PlotComposition2(pseq_comp, taxa=TopN.Taxa, sample.sort=level,
  #                      filename=paste0("./tmp_", idx, ".",level,".pdf"),
  #                      height=10, width=21, legend.ncol=10)
  
  # tmp.list[[idx]] %>% { prune_taxa( taxa_sums(.)/N > 0.001, . ) } %>% taxa
  # tmp.list[[idx]] %>% top_taxa(14)
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  TOP_taxa <- pseq@tax_table[top_taxa(pseq, n = 1), sample.sort]
  
  fit.comp <- pseq %>%
    plot_composition(sample.sort = TOP_taxa, otu.sort = "abundance")
  
  # fit.comp$df %>% head %>% print
  
  if( is.null(taxa) ){
    taxa <- levels(fit.comp$df$Tax)
  }
  
  if( length(taxa) > 74 ) warnings("The maximum number of colors is 74.")
  
  
  p <- fit.comp$plot +
    scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
    scale_y_continuous(label = scales::percent) +
    # hrbrthemes::theme_ipsum(grid="Y") +
    theme_bw(base_size=16) +
    labs(x = "Samples", y = "Relative Abundance",
         title = "Relative Abundance data") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_legend(ncol = legend.ncol))
  
  
  
  pdf(file=filename, height=height, width=width)
  print(p)
  dev.off()
  
  
  invisible(fit.comp)
  
}

