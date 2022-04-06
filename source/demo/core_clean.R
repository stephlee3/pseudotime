## input: 
# raw count 
# cell meta
# sample meta
# tree structure

## output:
# pb.ls, ctp.ls
# sample meta
runPseudotime = function(raw_count, 
                         cell_meta, 
                         sample_meta, 
                         leaves_info,
                         features = 'expr_and_prop', 
                         batch = NULL, 
                         pb_filter_pct = 0.9, 
                         pb_HVG = T, 
                         pb_qnorm = T, 
                         pb_batch_correction = F, 
                         pb_parallel = F,
                         ctp_clr = F, 
                         ctp_adj_cov = NULL, 
                         ctp_batch_correction = F, 
                         sample_name = 'sample', 
                         severity_name = 'severity', 
                         severity_num_name = 'sev.level',
                         ptime_method = 'cca', 
                         clsutermethod = 'louvain', 
                         cluster_id = NULL, 
                         reduce = F, 
                         clusternum = 2, 
                         startcluster = 1, 
                         near_nbr_num = 3, 
                         modelNames = 'VVV',
                         pc_scale = T,
                         eval_method = 'corr'
                         
){
  
  if (features %in% c('expr','expr_and_prop')){
    pb.ls = get_tree_node_feature(
      leaves_info, features = 'expr',
      raw_count = raw_count, cell_meta = cell_meta, filter_pct = pb_filter_pct,
      HVG = pb_HVG, qnorm = pb_qnorm, batch_correction = pb_batch_correction, batch = batch,
      parallel = pb_parallel
    )
    
    if (ptime_method == 'cca'){
      expr.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'expr',
        pb.ls = pb.ls, sample_metadata = sample_meta, ptime_method = 'cca',
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        pc_scale = pc_scale
      )
    }
    if (ptime_method == 'tscan'){
      expr.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'expr',
        pb.ls = pb.ls, sample_metadata = sample_meta, 
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        ptime_method = 'tscan', clustermethod = cluster_method, modelNames = modelNames,
        cluster = cluster_id, reduce = F, clusternum = clusternum, startcluster = startcluster, near_nbr_num = near_nbr_num, 
        pc_scale = T
      )
    }
    
  }
  
  if (features %in% c('prop','expr_and_prop')){
    ctp.ls = get_tree_node_feature(
      leaves_info = leaves_info, features = 'prop',
      raw_count = NULL, cell_meta = cell_meta, filter_pct = NULL,
      clr = ctp_clr, adj_cov = ctp_adj_cov, batch_correction = ctp_batch_correction, batch = batch
    )
    
    if (ptime_method == 'cca'){
      prop.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'prop',
        ctp.ls = ctp.ls, sample_metadata = sample_meta, ptime_method = 'cca',
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        pc_scale = pc_scale
      )
    }
    if (ptime_method == 'tscan'){
      prop.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'prop',
        ctp.ls, sample_metadata = sample_meta, 
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        ptime_method = 'tscan', clustermethod = cluster_method, modelNames = modelNames,
        cluster = cluster_id, reduce = F, clusternum = clusternum, startcluster = startcluster, near_nbr_num = near_nbr_num, 
        pc_scale = T
      )
    }
  }
  
  if (features == 'expr'){
    ptime.ls = expr.ptime
  } else if (features == 'prop'){
    ptime.ls = prop.ptime
  } else if (features == 'expr_and_prop'){
    ptime.ls = c(expr.ptime, prop.ptime)
  }
  
  #print(str(ptime.ls))
  auto_ptime = auto_pseudotime_sel(ptime.ls, ptime_method = ptime_method, eval_method = eval_method)
  sample_meta_ptime = auto_ptime$opt
  
  
  if (features == 'expr'){
    ans = list(
      pb.ls = pb.ls,
      ptime.ls = ptime.ls, 
      sample_meta_ptime = sample_meta_ptime
    )
  } else if (features == 'prop'){
    ans = list(
      ctp.ls = ctp.ls,
      ptime.ls = ptime.ls, 
      sample_meta_ptime = sample_meta_ptime
    )
  } else if (features == 'expr_and_prop'){
    ans = list(
      pb.ls = pb.ls, 
      ctp.ls = ctp.ls,
      ptime.ls = ptime.ls,
      sample_meta_ptime = sample_meta_ptime
    )
  }
}

## input: 
# sample meta
# leaves_info: tree structure

## output:
# deg result
# deg summary

runPseudoDEG = function(sample_meta, 
                        pb.ls,
                        leaves_info, 
                        phenotype_name = 'severity',
                        sample_name = 'sample', 
                        pseudotime_name = 'pseudotime', 
                        cv_quantile_cutoff = 0.5, 
                        postcurve_cutoff = 0.5, 
                        filter_method = 'cv')
{
  pb.ls.node = pb.ls$node
  
  deg_result = lapply(1:length(pb.ls.node), function(i){
    #print(names(pb.ls.node)[i])
    res = run_pseudo_diff_gene(x = pb.ls.node[[i]], sample_metadata = sample_meta_ptime, cluster = names(pb.ls.node)[i],
    pseudotime_name = pseudotime_name, cv_quantile_cutoff = cv_quantile_cutoff, filter_method = filter_method)
    #print(head(res))
    return(res)
  })
  names(deg_result) = names(pb.ls.node)
  
  num_DEGs = sapply(1:length(deg_result), function(i){
    num_DEGs = length(which(deg_result[[i]]$signif))
  })
  names(num_DEGs) = names(deg_result)
  #print(num_DEGs)
  num_DEGs_col_name = paste0(phenotype_name, '.num_DEGs')
  
  
  deg_summary = leaves_info %>%
    dplyr::select(label, x, y ,id, leaf) %>%
    unique() %>%
    mutate(!!num_DEGs_col_name := 0)
  match_idx = match(names(num_DEGs), deg_summary$label)
  deg_summary[match_idx, num_DEGs_col_name] = num_DEGs
  
  return(list(tab = deg_result, num_DEG = deg_summary))
}


