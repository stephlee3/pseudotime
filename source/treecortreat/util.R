library(Matrix)
library(preprocessCore)

get_cluster_pb = function(s, clu, pt,
                          topclu = NULL, batch = NULL, 
                          HVG = T, qnorm = F, combat = F){
  if (is.factor(clu))
    uc = levels(clu)
  else uc = unique(clu) 
  m <- sapply(uc,function(sc) {
    tmp <- s[,which(clu==sc), drop = F]
    p <- pt[which(clu == sc)]
    sapply(unique(p),function(sp) rowSums(tmp[,p==sp,drop=F]))
  })
  d <- lapply(m,function(i) {
    rc <- colSums(i)/1e6
    log2(t(t(i)/rc + 1))
  })
  names(d) = uc
  print(str(d))
  
  ## combat
  if(is.null(topclu))
    topclu = 1:length(s)
  s = s[topclu]
  d <- sapply(1:length(s),function(i) {
    d <- s[[i]][!grepl('MT-',rownames(s[[i]])),]
    d <- d[rowMeans(d > 0.01) > 0.1,]
    
    ## HVG
    if (HVG){
      cm <- rowMeans(d)
      csd <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
      mod <- loess(csd~cm)
      rid <- which(resid(mod) > 0)
    } else{
      rid = 1:nrow(d)
    }
    
    
    ## qnorm
    if (qnorm){
      dn<- dimnames(d)
      d<- normalize.quantiles(d); dimnames(d) = dn;
    }
    
    ## combat
    if (combat){
      match_id = match(colnames(d),names(batch))
      d = ComBat(d,batch = batch[match_id])
    }
    
    
    d = d[rid,]
    rownames(d) <- paste0(rownames(d),':',names(s)[i])
    return(d)
  }) 
  n <- table(unlist(sapply(d,colnames)))
  n <- names(n)[n==length(d)]
  d <- do.call(rbind,lapply(d,function(i) i[,n]))
  print(str(d))
  return(d)
}

get_sample_pb = function(s, pt, HVG = T, qnorm = F, combat = F, batch = NULL){
  m <- sapply(unique(pt),function(sp){
    rowSums(s[,pt ==sp,drop=F])
  })
  rc <- colSums(m)/1e6
  d <- log2(t(t(m)/rc + 1))
  
  ## pre-filtering
  d = d[!grepl('MT-', rownames(d)),]
  d = d[rowMeans(d > 0.01) > 0.1, ]
  
  ## hvg
  if (HVG){
    cm <- rowMeans(d)
    csd <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
    mod <- loess(csd~cm)
    rid <- which(resid(mod) > 0)
  } else{
    rid = 1:nrow(d)
  }
  
  ## qnorm
  if (qnorm){
    dn<- dimnames(d)
    d<- normalize.quantiles(d); dimnames(d) = dn;
  }
  
  ## combat
  if (combat){
    match_id = match(colnames(d),names(batch))
    d = ComBat(d,batch = batch[match_id])
  }
  
  d = d[rid,]
  return(d)
}


clr_transform = function(X,pseudocount = 1){
  # X: each row is a cluster, each column is a sample
  X = X + pseudocount # pseudocount
  ans = apply(X, 2, function(x){log(x/(exp(mean(log(x)))))})
  return(ans)
}

remove_batch_effect = function(dat, batch, adj_cov = NULL){
  batchmod = model.matrix(~ -1 + batch)
  nbatch = ncol(batchmod)
  design = cbind(batchmod, adj_cov)
  
  dat.adj = t(sapply(1:nrow(dat),function(r){
    y = dat[r,]
    beta.hat = solve(crossprod(design), crossprod(design, y))
    gamma.hat = batchmod %*% beta.hat[1:nbatch,]
    grand.mean = mean(gamma.hat)
    gamma.hat = gamma.hat-grand.mean
    y - gamma.hat
  }))
  dimnames(dat.adj) = dimnames(dat)
  return(dat.adj)
}


get_ctp = function(clu, pt, clr = F, batch_correction = F, batch = NULL, adj_cov = NULL){
  cell_type_count = table(clu, pt) %>% unclass()
  cell_type_prop = apply(cell_type_count,2,function(x) x/sum(x))
  
  if(clr){
    cell_type_prop = clr_transform(cell_type_count, 1)
  }
  print(str(cell_type_prop))
  if(batch_correction){
    match_id = match(colnames(cell_type_prop), names(batch))
    print(str(batch[match_id]))
    cell_type_prop = remove_batch_effect(cell_type_prop, batch[match_id], adj_cov =  adj_cov)
  }
  print(str(cell_type_prop))
  return(cell_type_prop)
}
