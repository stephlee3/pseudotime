library(tidyverse)
library(TreeCorTreat)

source("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/util.R")
source("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/core.R")

get_eval_metric = function(sample_metadata,
                           sample_name = 'sample',
                           pseudotime_name = 'pseudotime',
                           severity_num_name = 'sev.level',
                           eval_method = c('corr','pairwise_order')){
  eval_method = match.arg(eval_method)
  
  if (eval_method == 'corr'){
    ans = cor(sample_metadata[,severity_num_name], sample_metadata[,pseudotime_name],
              use = 'complete')
  }
  
  if (eval_method == 'pairwise_order'){
    pt = sample_metadata[, sample_name]
    submeta = sample_metadata %>% select(!!sample_name, !!pseudotime_name, !! severity_num_name)
    df = expand.grid(x1 = pt, x2 = pt, stringsAsFactors = F) %>%
      filter(x1 != x2)
    df = df %>%
      inner_join(submeta, by = c('x1' = sample_name)) %>%
      inner_join(submeta, by = c('x2' = sample_name)) 
    colnames(df) = c('x1','x2','ptime1','sev1','ptime2','sev2')
    df = df %>%
      mutate(correct_order = ((sev1 <= sev2) == (ptime1 <= ptime2)))
    ans = mean(df$correct_order)
    
  }
  return(ans)
}

## leaves info

pseudotime_tree_node = function(leaves_info, features = c('expr','prop'), pb.ls = NULL, ctp.ls = NULL,
                                raw_count, cell_clu, cell_sample,
                                sample_metadata = NULL,
                                sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                ptime_method = c('tscan','cca'), cluster = NULL, reduce = F, clusternum = 2,
                                startcluster = 1, HVG = T, qnrom = F, batch_correction = F, batch = NULL, 
                                clr = F, adj_cov = F, pc_scale = F
                                ){
  unq_id <- unique(leaves_info$label)
  unq_y <- setdiff(unique(leaves_info$y),max(leaves_info$y)) #remove root
  unq_y_root = unique(leaves_info$y)
  
  
  if (features == 'expr'){
    if (is.null(pb.ls)){
      pb.ls <- lapply(unq_id,function(tid){
        print(tid)
        node_info <- leaves_info %>% filter(label==tid)
        print(node_info)
        # note: cell_meta must contain these columns: barcode, sample, celltype
        sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
        sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]
        
        
        # pseudobulk construction (please fill in)
        pb <- get_sample_pb(s = sub_count, pt = sub_meta$sample, 
                            HVG = HVG, qnorm = qnorm, combat = batch_coorection, batch = batch)
        print(str(pb))
        return(pb)
      })
    }
    
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
    

    ptime.ls = lapply(1:length(pb.ls.agg), function(i){
      pb = pb.ls.agg[[i]]
      pbpca = prcomp(t(pb),scale. = pc_scale)$x[,1:min(20, nrow(pb))]                  
      pbpca = pbpca[rownames(pbpca) %in% sample_metadata[,sample_name],]
      sample_meta_ptime = get_sample_pseudotime(x = t(pbpca), sample_metadata = sample_metadata, 
                                                sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                method = ptime_method)
    })
    names(ptime.ls) = paste0('expr_', unq_y_root)
   
  }
  
  if(features == 'prop'){
      if (is.null(ctp.ls)){
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
        
        ptime.ls = lapply(1:length(ctp.ls), function(i){
          ctp = ctp.ls[[i]]
          ctpca = prcomp(t(ctp),scale. = pc_scale)$x[,1:min(20, nrow(ctp))]                 
          ctpca = ctpca[rownames(ctpca) %in% sample_metadata[,sample_name],]
          sample_meta_ptime = get_sample_pseudotime(x = t(ctpca), sample_metadata = sample_metadata, 
                                                    sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                    method = ptime_method)
        })
        names(ptime.ls) = paste0('prop_', unq_y, '_clr', clr, '_batch', batch_correction)
      }
  }
  
  return(ptime.ls)
}



## automatic optimal pseudotime selection
auto_pseudotime_sel = function(ptime.ls, 
                               sample_name = 'sample', pseudotime_name = 'pseudotime', severity_num_name = 'sev.level',
                               eval_method = c('corr', 'pairwise_order')
                               ){
  eval_method = match.arg(eval_method)
  result = sapply(1:length(ptime.ls), function(i){
    ptime = ptime.ls[[i]]
    ans = get_eval_metric(sample_metadata = ptime, sample_name = sample_name, pseudotime_name = pseudotime_name, severity_num_name = severity_num_name,
                          eval_method = eval_method)
  })
  
  optime = ptime.ls[which.max(result)]
  names(result) = names(ptime.ls)
  result = sort(result, decreasing = T)
  return(list(opt = optime, eval = result))
}


## main
## load data
s <- paste0('@All(', 
            '@T(',paste0('T_c',c(1,2,3,4,8), collapse = ','),'),',
            '@Mono(',paste0('Mono_c',c(13,14), collapse = ','),'),',
            '@B(',paste0('B_c',c(24,25), collapse = ','),'))')

s <- paste0('@All(', 
            '@T(',paste0('T_c',c(1:12), collapse = ','),'),',
            '@Mono(',paste0('Mono_c',c(13:23), collapse = ','),'),',
            '@B(',paste0('B_c',c(24:28), collapse = ','),'),',
            '@Neutrophil(',paste0('Neutrophil_c',c(29:30), collapse = ','),'),',
            '@Megakaryocyte(',paste0('Megakaryocyte_c',c(31), collapse = ','),'))')
hierarchy_list <- extract_hrchy_string(s,'@') # tree is plotted in bottom-up fashion.
leaves_info <- hierarchy_list$leaves_info

sample_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/sample_meta.rds")
sample_meta = sample_meta %>%
  filter(severity %in% c('HD','Mi','Mod','Se')) %>%
  mutate(sev.level = case_when(
    severity == 'HD' ~ 1,
    severity == 'Mi' ~ 2,
    severity == 'Mod' ~ 3,
    severity == 'Se' ~ 4
  ))
batch = sample_meta$batch; names(batch) = sample_meta$sample

cell_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/cell_meta.rds")
cell_meta = cell_meta %>%
  mutate(large = gsub('(.*)_(.*)','\\1', celltype)) %>%
  filter(sample %in% sample_meta$sample)
clu = cell_meta$celltype
pt = cell_meta$sample

pb.ls = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/pb.ls.rds")

## run pipeline
ptime_expr = pseudotime_tree_node(leaves_info, features = 'expr', pb.ls = pb.ls,
                                  sample_metadata = sample_meta,
                                  ptime_method = 'cca', pc_scale = F)
combo = expand.grid(clr = c(T,F), 
                    batch_correction = c(T,F),
                    stringsAsFactors = F)
  
ptime_prop = lapply(1:nrow(combo), function(r){
  pseudotime_tree_node(leaves_info, features = 'prop', 
                       cell_clu = clu , cell_sample = pt,
                       sample_metadata = sample_meta,
                       ptime_method = 'cca', batch = batch, 
                       clr = combo$clr[r], batch_correction =  combo$batch_correction[r],
                       pc_scale = T)
  
})
ptime_prop = do.call(c, ptime_prop)


ptime.ls = c(ptime_expr, ptime_prop)
ptime.result = auto_pseudotime_sel(ptime.ls)
  
  


                                             
  
