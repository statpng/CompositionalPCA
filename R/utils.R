#' @export png.merge_pseq_list
png.merge_pseq_list <- function(pseq_list, label="Site"){
  
  if(FALSE){
    load("./data/pseq_list.rda")
    # pseq_list: list(data=pseq_total_phylum, joint=jpseq.jive, individual=apseq.jive)
    label="Site"
  }
  
  
  N <- length(pseq_list)
  
  if( !is.null(names(pseq_list)) ){
    Name <- names(pseq_list)
  } else {
    Name <- LETTERS[seq_len(N)]
  }
  
  
  otu_list <- lapply(pseq_list, otu_table)
  tax_list <- lapply(pseq_list, tax_table)
  sample_list <- lapply(pseq_list, sample_data)
  
  # taxa_all <- unique(do.call("c", lapply(otu_list, taxa_names)))
  
  
  # start <- proc.time()
  otu_new <- lapply(1:N, function(idx){
    x_new <- otu_list[[idx]] %>% { cbind.data.frame(Taxa=rownames(.), .) }
    colnames(x_new)[-1] <- paste0(Name[idx],"_",colnames(x_new)[-1])
    x_new
  })
  tax_new <- lapply(1:N, function(idx){
    tax_list[[idx]] %>% { cbind.data.frame(Taxa=rownames(.), .) }
  })
  
  otu_all <- otu_new %>% reduce(full_join, by="Taxa")
  otu_all[is.na(otu_all)] <- 0
  tax_all <- tax_new %>% reduce(full_join)
  tax_all[is.na(tax_all)] <- 0
  
  rownames(otu_all) <- otu_all[,1]
  otu_all <- as.matrix(otu_all[,-1])
  # colnames(otu_all) <- gsub("[a-z]+\\_", "", colnames(otu_all))
  
  # end <- proc.time()
  # end - start
  
  
  # sample_new %>% apply(2, function(x) table(x) %>% {c(head(.), tail(.))})
  
  sample_all <- lapply(1:N, function(idx){
    sample_new <- cbind.data.frame(Site=names(sample_list)[idx], sample_list[[idx]])
    colnames(sample_new) <- c(label, colnames(sample_list[[idx]]))
    sample_new
  }) %>% dplyr::bind_rows()
  
  
  my_sample_data <- sample_data( sample_all )
  my_taxonomyTable <- tax_table( as.matrix(tax_all[,-1]) )
  my_otu_table <- otu_table( otu_all, taxa_are_rows = TRUE )
  
  rownames(my_taxonomyTable) <- rownames(my_otu_table)
  my_sample_data@row.names <- colnames(my_otu_table)
  
  my_phyloseq <- phyloseq(my_otu_table, my_taxonomyTable, my_sample_data)
  
  my_phyloseq
}






#' @export png.df.split_taxa
png.df.split_taxa <- function(df, col="taxonomy", sep="; ", col_names=c("Kingdom","Phylum","Class","Order","Family", "Genus", "Species")){
  if(FALSE){
    df <- data_list[[1]]
    col="taxonomy"
    sep="; "
    col_names=c("Kingdom","Phylum","Class","Order","Family", "Genus", "Species")

    split_taxa(data_list[[1]])[1:5,1:10]
  }

  # df %>% tidyr::separate("taxonomy", c("Kingdom","Phylum","Class","Order","Family", "Genus", "Species"), sep="; ") %>% .[5504,1:10]

  # taxa_df <- df$taxonomy %>% strsplit("; ") %>% do.call("rbind", .)
  taxa_df <- df[[col]] %>% strsplit(sep) %>% do.call("rbind", .)
  colnames(taxa_df) <- col_names
  cbind.data.frame( taxa_df, df %>% select(-all_of(col)) )
}


#' @export png.otu2pseq
png.otu2pseq <- function(otu, taxa, sample){

  if(FALSE){
    otu <- data_list3[[1]][,-(1:7)]
    taxa <- data_list3[[1]][,1:7]
    sample <- info
  }

  otu[1:2,1:10]
  taxa %>% head(2)
  sample %>% head(2)
  # >     otu[1:2,1:10]
  # UR2199 UR2203 UR2208 UR2211 UR2212 UR2213 UR2214 UR2215 UR2218 UR2219
  # 114785      0      0      0      0      0      0      0      0      0      0
  # 110008      0      0      0      0      0      0      0      0      0      0
  # >     taxa %>% head(2)
  # Kingdom         Phylum          Class            Order
  # 114785 Bacteria Planctomycetes Planctomycetia Planctomycetales
  # 110008 Bacteria  Bacteroidetes  Flavobacteria Flavobacteriales
  # Family          Genus    Species
  # 114785 Planctomycetaceae Rhodopirellula FJ624356_s
  # 110008 Crocinitomicaceae     JN874375_g EU652656_s
  # >     sample %>% head(2)
  # urine  serum  STool Age Sex     Date AST ALT Glucose Total.Cholesterol
  # 1 UR2199 SR0622 ST0307  52   M 20150811  30  42    147                215
  # 2 UR2203 SR0626 ST0308  56   M 20150818  18  13    131                136
  # Creatinine vGTP Triglyceride HDL.Cholesterol  Hgb
  # 1       0.86   45           62              54 14.9
  # 2       1.14   40           60              36 16.1


  my_otu_table <- otu_table( otu, taxa_are_rows = TRUE )
  my_taxonomyTable <- tax_table( as.matrix(taxa) )
  my_sample_data <- sample_data( sample )

  #check
  # identical( my_sample_data[[info_name]], colnames(my_otu_table) )

  rownames(my_otu_table) <- rownames(my_taxonomyTable)
  colnames(my_otu_table) <- my_sample_data@row.names

  my_phyloseq <- phyloseq(my_otu_table, my_taxonomyTable, my_sample_data)

  return(my_phyloseq)

}






#' @export Remove_All_Zeros
Remove_All_Zeros <- function(df, axis){
  df_taxa <- df
  df <- df %>% select(-taxonomy)
  AllZeros <- apply( df, axis, function(x) all(x==0) )
  if( sum(AllZeros) > 0 ){
    if(axis == 1){
      return( df_taxa[-which(AllZeros),] )
    } else {
      return( df_taxa[,-which(AllZeros)] )
    }
  }
  return( df_taxa )
}



#' @export png.make_zero
png.make_zero <- function(x, cutoff=0.01){
  # x: compositional data vector with zero values
  if(FALSE){
    x = c(0.01,0.01,0.5,0.48)
    cutoff=0.01
  }

  x[x<=cutoff] <- 0
  x[x>cutoff] <- x[x>cutoff] / sum(x[x>cutoff])

  x
}


#' @export png.replace_simple
png.replace_simple <- function(x, delta=0.01){
  # x: compositional data vector with zero values
  if(FALSE){
    x = c(0,0.1,0.8,0.1,0)
    delta=0.01
  }
  x[x==0] <- delta
  x/sum(x)
}




#' @export png.pca.approx
png.pca.approx <- function(X, nrank){
  # X: n x p dimensional matrix in real space
  n=nrow(X)
  mu=colMeans(X)

  with( prcomp(X),
        tcrossprod(rep(1,n),mu)+tcrossprod(x[,1:nrank], rotation[,1:nrank]) )
}



#' @export png.pca.reconstruction
png.pca.reconstruction <- function(fit, nrank){
  mu <- fit$mu
  uhat <- matrix(fit$uhat[,1:nrank], ncol=nrank)
  vhat <- matrix(fit$vhat[,1:nrank], ncol=nrank)

  n=nrow(uhat); p=nrow(vhat)
  return( tcrossprod(rep(1,n),mu) + tcrossprod(uhat,vhat) )
}






#' @export png.cpca.conv
png.cpca.conv <- function(LIST, eps=1e-6, accept0.1=FALSE){
  
  if( accept0.1 ){
    conv.list <- LIST %>% purrr::map(function(.x){
      out <- try(purrr::map_dbl(.x$fit.path, function(.y){
        .y$crit.path %>% tail(1) %>% {.<=eps | .==0.1}
      }))
    })
  } else {
    conv.list <- LIST %>% purrr::map(function(.x){
      out <- try(purrr::map_dbl(.x$fit.path, function(.y){
        .y$crit.path %>% tail(1) %>% {.<=eps}
      }))
    })
  }
  
  out <- do.call("rbind",conv.list)
  colnames(out) <- paste0("rank=",1:ncol(out))
  out
}




#' @export png.proj2simplex
png.proj2simplex <- function(v){
  ord = order(v, decreasing = TRUE)
  v = v[ord]
  s = 0
  for (i in 1:length(v)) {
    s = s + v[i]
    if (i == length(v)) {
      delta = (s - 1)/i
      break
    }
    d = s - v[i + 1] * i
    if (d >= 1) {
      delta = (s - 1)/i
      break
    }
  }
  v[1:i] = v[1:i] - delta
  v[-(1:i)] = 0
  w = rep(0, length(v))
  w[ord] = v
  w
}





# Depreciated -------------------------------------------------------------

#' @export png.NormalVector
png.NormalVector <- function(x,v){
  x-sum(x*v)/norm(v,"2") * v/norm(v,"2")
}




