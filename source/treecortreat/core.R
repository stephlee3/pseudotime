# `x`: the data matrix of features $\times$ samples, i.e. each row is a feature of gene expression or cell type proportion, each column is a patient sample.

# `sample_metadata`: patient-level metadata, including several important variables.

# `sample_name` and `severity_name`: the variable name in `sample_metadata` that indicates the sample and severity. The default setting of `sample_name` is "Patient" and `severity_name` is "severity". 

# `method`: users now have two options, either using TSCAN for sample-level pseudotime construction or use CCA for optimal severity axis detection. 

# `cluster = NULL, reduce = F, clusternum = 2, startcluster = 1`: parameters for TSCAN.

# `severity_num_name`: the variable name in 'sample_metadata' that indicates the numerical valaue of severity severity. 

get_sample_pseudotime = function(x, sample_metadata = NULL, 
                                 sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                 method = c('tscan','cca'),
                                 cluster = NULL, reduce = F, clusternum = 2, startcluster = 1){
  method = match.arg(method)
  
  ## arrange sample metadata to have the same order as x
  match_idx = match(colnames(x), sample_metadata[,sample_name])
  sample_metadata = sample_metadata[match_idx, ]
  
  ## tscan
  if (method == 'tscan'){
    sample_cluster = exprmclust(x, 
                                reduce = reduce,
                                cluster = cluster,
                                clusternum = clusternum)
    sample_metadata$cluster = sample_cluster$clusterid
    
    sample_order = TSCANorder(sample_cluster,
                              startcluster = startcluster,
                              orderonly = F,
                              listbranch = F)
    
    sample_info = sample_metadata %>%
      inner_join(sample_order, by = setNames("sample_name",sample_name)) %>%
      arrange(Pseudotime)
  }
  if (method == 'cca'){
    if (reduce) {
      pc = prcomp(t(x), scale = T)$x[1:min(20, nrow(x))]
    }
    else {
      pc = t(x)
    }
    severity_num = sample_metadata[,severity_num_name]
    res = cancor(pc,severity_num)
    proj = (res$xcoef["PC1",1]*(pc[,"PC1"]-res$xcenter["PC1"]))+ (res$xcoef["PC2",1]*(pc[,"PC2"]-res$xcenter["PC2"]))  
    pseudotime = rank(proj)
    if (cor(pseudotime, severity_num)<0)
      pseudotime = length(pseudotime)+1 -pseudotime
    slope = res$xcoef["PC2",1]/res$xcoef["PC1",1]
    print(paste("slope is", slope))
    print(paste("corr is ", res$cor))
    
    sample_info = sample_metadata %>%
      mutate(pseudotime = pseudotime) %>%
      arrange(pseudotime)
  }
  
  return(sample_info)
}

run_pseudo_diff_gene = function(x, sample_metadata = NULL, cluster = 'all',
                                sample_name = 'sample', pseudotime_name = 'pseudotime', 
                                cv_quantile_cutoff = 0.5, postcurve_cutoff = 0.5, filter_method = 'cv'
){
  ptime = sample_metadata[,pseudotime_name]
  sample_metadata = sample_metadata[order(ptime),]
  
  
  genes = rownames(x)
  x_reorder = x[,match(sample_metadata[,sample_name],colnames(x)) %>% na.omit()]
  
  if(filter_method == 'cv'){
    cv = apply(x_reorder, 1, function(x) sd(x)/mean(x))
    top_genes = which(cv > quantile(cv, cv_quantile_cutoff))
    eff = cv[top_genes]
    x_reorder = x_reorder[top_genes,]
  }
  
  
  
  res = lapply(1:nrow(x_reorder), function(i){
    model = mgcv::gam(x_reorder[i,]~s(c(1:ncol(x_reorder)),k=3))
    smooth = model$fitted.values
    eff = (max(smooth)-min(smooth))/sqrt(sum((x_reorder[i,]-smooth)^2))
    pval = pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F) 
    list(smooth = smooth, pval = pval, eff = eff)
  })
  pval = sapply(res, function(i) i$pval)
  fdr = p.adjust(pval, method = 'fdr')
  xsmooth = lapply(res, function(i) i$smooth) %>% do.call(rbind, .)
  eff = sapply(res, function(i) i$eff)
  
  if (filter_method == 'postcurve'){
    top_genes = which(eff > quantile(eff, postcurve_cutoff))
    eff = eff[top_genes]
    x_reorder = x_reorder[top_genes,]
    xsmooth = xsmooth[top_genes,]
    pval = pval[top_genes]
    fdr = fdr[top_genes]
  }
  
  result = data.frame(gene = rownames(x_reorder), pval = pval, fdr = fdr, 
                      effect_size = eff, filter_method = filter_method ) %>%
    mutate(cell_cluster = cluster,
           signif = (fdr < 0.05)) %>%
    arrange(fdr, effect_size)
  
  gene_symbols = result$gene
  gene_ids = AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns= c("ENSEMBL", "ENTREZID", "SYMBOL"))
  result = result %>%
    left_join(gene_ids, by = c('gene' = 'SYMBOL'))
  return(result)
}

get_tree_node_feature = function(leaves_info, features = c('expr','prop'), 
                                 raw_count, cell_meta, filter_pct = 0.9, 
                                 HVG = F, qnorm = F, batch_correction = F, batch = NULL, 
                                 clr = F, adj_cov = NULL, parallel = T, num_cores = 10, 
){
  unq_id <- unique(leaves_info$label)
  unq_y <- setdiff(unique(leaves_info$y),max(leaves_info$y)) #remove root
  unq_y_root = unique(leaves_info$y)
  tot_sample_size = length(unique(cell_meta$sample))
  
  
  if (features == 'expr'){
    if (parallel){
      pb.ls <- mclapply(unq_id, function(tid){
        print(tid)
        node_info <- leaves_info %>% filter(label==tid)
        print(node_info)
        # note: cell_meta must contain these columns: barcode, sample, celltype
        sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
        sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]
        unq_sample_size = length(unique(sub_meta$sample))
        if (unq_sample_size < filter_pct * tot_sample_size){
          print(paste("insufficient sample size for node", tid))
          return(NULL)
        }
        # pseudobulk construction (please fill in)
        pb <- get_sample_pb(s = sub_count, pt = sub_meta$sample, 
                            HVG = HVG, qnorm = qnorm, combat = batch_correction, batch = batch)
        print(str(pb))
        return(pb)
      }, mc.cores = num_cores)
    } else{
      pb.ls <- lapply(unq_id,function(tid){
        print(tid)
        node_info <- leaves_info %>% filter(label==tid)
        print(node_info)
        # note: cell_meta must contain these columns: barcode, sample, celltype
        sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
        sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]
        unq_sample_size = length(unique(sub_meta$sample))
        if (unq_sample_size < filter_pct * tot_sample_size){
          print(paste("insufficient sample size for node", tid))
          return(NULL)
        }
        # pseudobulk construction (please fill in)
        pb <- get_sample_pb(s = sub_count, pt = sub_meta$sample, 
                            HVG = HVG, qnorm = qnorm, combat = batch_correction, batch = batch)
        print(str(pb))
        return(pb)
      })
    }
    
    
    names(pb.ls) = unq_id
    pb.ls <- pb.ls[-which(sapply(pb.ls, is.null))]
    
    
    pb.ls.agg = lapply(unq_y_root, function(ty){
      label_names = leaves_info %>% filter(y == ty) %>% pull(label) %>% unique()
      pb.ls.sub = pb.ls[which(names(pb.ls) %in% label_names)]
      pb.ls.agg = lapply(1:length(pb.ls.sub), function(i){
        pb = pb.ls.sub[[i]]
        rownames(pb) = paste0(rownames(pb),':', names(pb.ls.sub)[i])
        return(pb)
      })  
      samp <- table(unlist(sapply(pb.ls.agg,colnames)))
      samp <- names(samp)[samp==length(pb.ls.agg)]
      pb.ls.agg <- do.call(rbind,lapply(pb.ls.agg,function(pb) pb[,samp]))
    })
    names(pb.ls.agg) = unq_y_root
    print(str(pb.ls.agg))
    
    pb.all = list(node = pb.ls, agg = pb.ls.agg)
    return(pb.all)
    
  }
  
  if(features == 'prop'){
    cell_clu = cell_meta$celltype
    cell_sample = cell_sample$sample
    
    ctp.ls <- lapply(unq_y, function(ty){
      #print(ty)
      node_info <- leaves_info %>% filter(y==ty)
      #print(node_info)
      cell_clu = cell_clu[cell_clu %in% node_info$children]
      cell_sample = cell_sample[cell_clu %in% node_info$children]
      
      match_idx = match(cell_clu, node_info$children)
      cell_clu= node_info$label[match_idx]
      ctp = get_ctp(clu = cell_clu, pt = cell_sample, 
                    clr = clr, batch_correction = batch_correction, 
                    batch = batch, adj_cov = NULL)
      print(str(ctp))
      return(ctp)
    })
    names(ctp.ls) = unq_y
    return(ctp.ls)
  }
  
}
pseudotime_tree_node_feature = function(leaves_info, features = c('expr','prop'), 
                                        raw_count, cell_meta, pb.ls, ctp.ls, 
                                        sample_metadata = NULL,
                                        sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                        ptime_method = c('tscan','cca'), cluster = NULL, reduce = F, clusternum = 2,
                                        startcluster = 1, pc_scale = F
){
  unq_id <- unique(leaves_info$label)
  unq_y <- setdiff(unique(leaves_info$y),max(leaves_info$y)) #remove root
  unq_y_root = unique(leaves_info$y)
  
  
  if (features == 'expr'){
    pb.ls.agg = pb.ls$agg
    ptime.ls = lapply(1:length(pb.ls.agg), function(i){
      pb = pb.ls.agg[[i]]
      rvar = apply(pb, 1, var)
      pb = pb[which(rvar > 0),] # 
      pbpca = try(prcomp(t(pb),scale. = pc_scale)$x[,1:min(20, nrow(pb))])
      if ('try-error' %in% class(pbpca))
        return(NULL)
      pbpca = pbpca[rownames(pbpca) %in% sample_metadata[,sample_name],]
      sample_meta_ptime = get_sample_pseudotime(x = t(pbpca), sample_metadata = sample_metadata, 
                                                sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                method = ptime_method, cluster = cluster, reduce = reduce, clusternum = clusternum)
    })
    names(ptime.ls) = paste0('expr_', unq_y_root)
    return(ptime.ls)
  }
  
  if(features == 'prop'){
    ptime.ls = lapply(1:length(ctp.ls), function(i){
      ctp = ctp.ls[[i]]
      rvar = apply(ctp, 1, var)
      ctp = ctp[which(rvar > 0),]
      ctpca = prcomp(t(ctp),scale. = pc_scale)$x[,1:min(20, nrow(ctp))]                 
      ctpca = ctpca[rownames(ctpca) %in% sample_metadata[,sample_name],]
      sample_meta_ptime = get_sample_pseudotime(x = t(ctpca), sample_metadata = sample_metadata, 
                                                sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                method = ptime_method, cluster = cluster, reduce = reduce, clusternum = clusternum)
    })
    names(ptime.ls) = paste0('prop_', unq_y)
    return(ptime.ls)
  }
}
