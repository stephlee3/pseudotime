library(tidyverse)
library(mgcv)
library(org.Hs.eg.db)
library(parallel)
library(Matrix)
library(TreeCorTreat)



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

## main
sample_meta_ptime = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/sample_meta_ptime.rds")
pb.ls = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/pb.ls.rds")

deg_result = lapply(1:length(pb.ls), function(i){
  print(names(pb.ls)[i])
  res = run_pseudo_diff_gene(x = pb.ls[[i]], sample_metadata = sample_meta_ptime, cluster = names(pb.ls)[i],
          cv_quantile_cutoff = 0.5, postcurve_cutoff = 0.5, filter_method = 'cv')
  print(head(res))
  return(res)
})
names(deg_result) = names(pb.ls)
saveRDS(deg_result, "/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/deg_result.rds")

## deg summary
deg_result = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/deg_result.rds")
s <- paste0('@All(', 
            '@T(',paste0('T_c',c(1,2,3,4,8), collapse = ','),'),',
            '@Mono(',paste0('Mono_c',c(13,14), collapse = ','),'),',
            '@B(',paste0('B_c',c(24,25), collapse = ','),'))')

print("Construct hierarchical tree structure...")
pdf('tree.pdf')
hierarchy_list <- extract_hrchy_string(s,'@') # tree is plotted in bottom-up fashion.
dev.off()

hierarchy_list$layout
hierarchy_list$leaves_info


label_info <- hierarchy_list$layout
if(sum(duplicated(label_info$label))) stop('Duplicated cell type name.')
leaves_info <- hierarchy_list$leaves_info

deg_summary = lapply(1:length(deg_result), function(i){
  num_DEGs = length(which(deg_result[[i]]$signif))
  df = leaves_info %>%
    dplyr::select(label, x, y ,id, leaf) %>%
    unique() %>%
    filter(label == names(deg_result)[i]) %>%
    mutate(num_DEGs = num_DEGs)
})
deg_summary = deg_summary %>% do.call(rbind, .)
saveRDS(deg_summary, "/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/deg_summary.rds")
