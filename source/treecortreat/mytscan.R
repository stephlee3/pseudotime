library(TSCAN)
library(igraph)
library(patchwork)
mclustobj = exprmclust(t(rd), clusternum = 4, modelNames = 'VII',reduce = F)
p1 = plotmclust(mclustobj)
TSCANproj = myTSCANorder(mclustobj, listbranch = T)
backbone = data.frame(TSCANproj$orth_proj$`backbone 1,2,3`)
backbone = data.frame(TSCANproj$orth_proj$`branch: 1,2,4`)
backbone = data.frame(TSCANproj$orth_proj$`branch: 2,3,4`)

orth_dist = do.call(c,TSCANproj$orth_dist) 
names(orth_dist) = do.call(c,lapply(TSCANproj$orth_dist, names))
orth_dist = orth_dist[which(!duplicated(names(orth_dist)))]
max(orth_dist)
sum(orth_dist)
R2 = 1- sum(orth_dist)/tot_var


colnames(backbone) = c('PCA_dimension_1','PCA_dimension_2')

p1 +
  geom_point(aes(x = PCA_dimension_1, y = PCA_dimension_2), data = backbone)







myTSCANorder <- function(mclustobj,startcluster=NULL,MSTorder = NULL,orderonly=F,flip=F,listbranch=F,divide=T) {
  if (!is.null(MSTorder) & length(MSTorder) == 1) {
    stop("MSTorder is not a path!")
  }
  set.seed(12345)
  clucenter <- mclustobj$clucenter
  row.names(clucenter) <- paste0("clu",1:nrow(clucenter))
  clusterid <- mclustobj$clusterid             
  pcareduceres <- mclustobj$pcareduceres            
  adjmat <- as_adjacency_matrix(mclustobj$MSTtree,sparse=FALSE)
  if (is.null(MSTorder)) {
    orderinMST <- 1
    clutable <- table(mclustobj$clusterid)
    alldeg <- degree(mclustobj$MSTtree)
    if (is.null(startcluster)) {
      allcomb <- expand.grid(as.numeric(names(alldeg)[alldeg==1]),as.numeric(names(alldeg)[alldeg==1]))
      allcomb <- allcomb[allcomb[,1] < allcomb[,2],]
      numres <- t(apply(allcomb, 1, function(i) {
        tmp <- as.vector(get.shortest.paths(mclustobj$MSTtree,i[1],i[2])$vpath[[1]])
        c(length(tmp), sum(clutable[tmp]))
      }))
      optcomb <- allcomb[order(numres[,1],numres[,2],decreasing = T)[1],]
      branchcomb <- allcomb[-order(numres[,1],numres[,2],decreasing = T)[1],,drop=F]
      MSTorder <- get.shortest.paths(mclustobj$MSTtree,optcomb[1],optcomb[2])$vpath[[1]]       
    } else {
      allcomb <- cbind(startcluster,setdiff(as.numeric(names(alldeg)[alldeg==1]),startcluster))
      numres <- t(apply(allcomb, 1, function(i) {
        tmp <- as.vector(get.shortest.paths(mclustobj$MSTtree,i[1],i[2])$vpath[[1]])
        c(length(tmp), sum(clutable[tmp]))
      }))
      optcomb <- allcomb[order(numres[,1],numres[,2],decreasing = T)[1],]
      branchcomb <- allcomb[-order(numres[,1],numres[,2],decreasing = T)[1],,drop=F]
      MSTorder <- get.shortest.paths(mclustobj$MSTtree,optcomb[1],optcomb[2])$vpath[[1]] 
    }
    if (flip) MSTorder <- rev(MSTorder)
  } else {
    edgeinMST <- sapply(1:(length(MSTorder)-1),function(i) {
      adjmat[MSTorder[i],MSTorder[i+1]]
    })
    if (divide) {
      if (sum(edgeinMST==0) > 0) {
        orderinMST <- 0
      } else {
        orderinMST <- 1
      }      
    } else {
      orderinMST <- 0
    }
  }          
  internalorderfunc <- function(internalorder,MSTinout) {
    TSCANorder <- NULL
    orth_proj<- NULL
    orth_dist<- NULL
    
    for (i in 1:(length(internalorder)-1)) {                  
      currentcluid <- internalorder[i]
      nextcluid <- internalorder[i + 1]
      currentclucenter <- clucenter[currentcluid,]
      nextclucenter <- clucenter[nextcluid,]
      
      currentreduceres <- pcareduceres[clusterid==currentcluid,]
      if (MSTinout) {
        connectcluid <- as.numeric(names(which(adjmat[currentcluid,] == 1)))      
      } else {
        if (i == 1) {
          connectcluid <- nextcluid      
        } else {
          connectcluid <- c(nextcluid,internalorder[i - 1])
        }                        
      }
      
      cludist <- sapply(connectcluid, function(x) {                              
        rowSums(sweep(currentreduceres,2,clucenter[x,],"-")^2)
      })
      mindistid <- apply(cludist,1,which.min)
      
      edgecell <- names(which(mindistid == which(connectcluid == nextcluid)))
      
      difvec <- nextclucenter - currentclucenter
      tmppos <- pcareduceres[edgecell,,drop=F] %*% difvec
      pos <- as.vector(tmppos)
      names(pos) <- row.names(tmppos)
      print('current cluster projection')
      print(str(pos))
      
      P <- difvec %*% t(difvec)/as.numeric(t(difvec) %*% difvec)
      edgecell_coord <- pcareduceres[edgecell,,drop=F]
      orth_proj_tmp <- t(sapply(1:nrow(edgecell_coord), function(r){
        P %*% edgecell_coord[r, ] + (diag(nrow(P))-P) %*% currentclucenter
      }))
      rownames(orth_proj_tmp) = rownames(edgecell_coord) 
      print(str(orth_proj_tmp))
      orth_dist_tmp = sapply(1:nrow(orth_proj_tmp), function(r){
        sum((edgecell_coord[r,]-orth_proj_tmp[r,])^2)
      })
      names(orth_dist_tmp) = rownames(edgecell_coord)

      orth_proj = rbind(orth_proj, orth_proj_tmp)
      orth_dist = c(orth_dist, orth_dist_tmp)
          
      
      
      TSCANorder <- c(TSCANorder,names(sort(pos)))  
      
      nextreduceres <- pcareduceres[clusterid==nextcluid,,drop=F]     
      if (MSTinout) {
        connectcluid <- as.numeric(names(which(adjmat[nextcluid,] == 1)))
      } else {
        if (i == length(internalorder)-1) {
          connectcluid <- currentcluid      
        } else {
          connectcluid <- c(currentcluid,internalorder[i + 2])
        }                        
      }
      
      cludist <- sapply(connectcluid, function(x) { 
        rowSums(sweep(nextreduceres,2,clucenter[x,],"-")^2)
      })
      if (length(cludist)==1) {
        mindistid <- 1
      } else {
        mindistid <- apply(cludist,1,which.min)
      }
      
      edgecell <- names(which(mindistid == which(connectcluid == currentcluid)))
      
      difvec <- nextclucenter - currentclucenter
      tmppos <- pcareduceres[edgecell,,drop=F] %*% difvec
      pos <- as.vector(tmppos)
      names(pos) <- row.names(tmppos)
      print('next cluster projection')
      print(str(pos))
      
      P <- difvec %*% t(difvec)/as.numeric(t(difvec) %*% difvec)
      edgecell_coord <- pcareduceres[edgecell,,drop=F]
      orth_proj_tmp <- t(sapply(1:nrow(edgecell_coord), function(r){
        P %*% edgecell_coord[r, ] + (diag(nrow(P))-P) %*% currentclucenter
      }))
      rownames(orth_proj_tmp) = rownames(edgecell_coord) 
      print(str(orth_proj_tmp))
      orth_dist_tmp = sapply(1:nrow(orth_proj_tmp), function(r){
        sum((edgecell_coord[r,]-orth_proj_tmp[r,])^2)
      })
      names(orth_dist_tmp) = rownames(edgecell_coord)
      
      orth_proj = rbind(orth_proj, orth_proj_tmp)
      orth_dist = c(orth_dist, orth_dist_tmp)
      
      TSCANorder <- c(TSCANorder,names(sort(pos)))
     
      
    }
    
    if (orderonly) {
      TSCANorder = TSCANorder      
    } else {
      #                   datadist <- dist(mclustobj$pcareduceres)
      #                   distmat <- as.matrix(datadist)
      #                   alldist <- sapply(1:(length(TSCANorder)-1), function(x) {
      #                         distmat[TSCANorder[x],TSCANorder[x+1]]
      #                   })
      #                   ptime <- c(0,cumsum(alldist))
      #                   ptime <- ptime/max(ptime) * 100
      TSCANorder = data.frame(sample_name=TSCANorder,State=clusterid[TSCANorder],Pseudotime=1:length(TSCANorder),stringsAsFactors = F)            
    }
    
    return(list(TSCANorder = TSCANorder, orth_proj = orth_proj, orth_dist = orth_dist))
  }
  if (!orderinMST) {
    internalorderfunc(MSTorder,0)            
  } else {
    if (exists("branchcomb") & listbranch) {                  
      allres <- list()
      backbone_result = internalorderfunc(MSTorder,1)
      allres["TSCANorder"] = list();allres[["orth_proj"]] = list(); allres[["orth_dist"]] = list()
      
      allres[["TSCAN_order"]][[paste("backbone",paste(MSTorder,collapse = ','))]] <- backbone_result$TSCANorder
      allres[["orth_proj"]][[paste("backbone",paste(MSTorder,collapse = ','))]] = backbone_result$orth_proj
      allres[["orth_dist"]][[paste("backbone",paste(MSTorder,collapse = ','))]] = backbone_result$orth_dist
      
      for (tmpcombid in 1:nrow(branchcomb)) {
        tmporder <- get.shortest.paths(mclustobj$MSTtree,branchcomb[tmpcombid,1],branchcomb[tmpcombid,2])$vpath[[1]] 
        if (flip) tmporder <- rev(tmporder)
        branch_result = internalorderfunc(tmporder,1)
        allres[["TSCAN_order"]][[paste("branch:",paste(tmporder,collapse = ','))]] <- branch_result$TSCANorder 
        allres[["orth_proj"]][[paste("branch:",paste(tmporder,collapse = ','))]] <- branch_result$orth_proj
        allres[["orth_dist"]][[paste("branch:",paste(tmporder,collapse = ','))]] <- branch_result$orth_dist
      }
      allres
    } else {
      internalorderfunc(MSTorder,1)            
    }      
  }
}
