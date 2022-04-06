library(tidyverse)
library(Matrix)
library(sva)
library(preprocessCore)


# s: A list of normalized gene expression data. Each element of the list represents the gene expression (gene \times samples) within a cell cluster. The number of elements equal to the number of cell clusters.
# cluster: indicate the assigned cell cluster for each cell.
# topclu: indicate the top clusters to preserve for pseudobulk construction


get_pseudobulk = function(s, cluster, topclu = NULL, batch = NULL, 
                          HVG = T, qnorm = T, combat = T){
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
  return(d)
}


## run pseudobulk
s = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/pb_norm.rds")
clu = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/cluster.rds")
clu = factor(clu,levels = paste0("c",1:length(unique(clu))))
s = s[match(levels(clu),names(s))]
sampsize = sapply(1:length(s),function(i) ncol(s[[i]])); names(sampsize) = names(s)
topclu = which(sampsize > 400)
batch = sub('.*-','',colnames(s[[1]]))
names(batch) = colnames(s[[1]])

pb = get_pseudobulk(s = s, cluster = clu, topclu = topclu, batch = batch, 
                    HVG = T, qnorm = T, combat = T)

print(str(pb))