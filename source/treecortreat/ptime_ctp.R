library(tidyverse)
cell_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/cell_meta.rds")
cell_meta = cell_meta %>%
  mutate(large = gsub('(.*)_(.*)','\\1', celltype))
clu = cell_meta$large
clu = cell_meta$celltype
pt = cell_meta$sample
sample_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/sample_meta.rds")
batch = sample_meta$batch; names(batch) = sample_meta$sample



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
  if(batch_correction){
    match_id = match(colnames(cell_type_prop), names(batch))
    cell_type_prop = remove_batch_effect(cell_type_prop, batch[match_id], adj_cov =  NULL)
  }
  return(cell_type_prop)
}

sample_meta = sample_meta %>%
  filter(severity %in% c('HD','Mi','Mod','Se')) %>%
  #filter(sample %in% colnames(ctp)) %>%
  mutate(sev.level = case_when(
    severity == 'HD' ~ 1,
    severity == 'Mi' ~ 2,
    severity == 'Mod' ~ 3,
    severity == 'Se' ~ 4
  ))

combo = expand.grid(clr = c(T,F), 
                  batch_correction = c(T,F),
                  stringsAsFactors = F)

cpt_result = lapply(1:nrow(combo), function(r){
  print(combo[r,])
  ctp = get_ctp(clu, pt, clr = combo$clr[r], batch_correction = combo$batch_correction[r], batch = batch, adj_cov = NULL)
  ctpca = prcomp(t(ctp), scale. = T)$x[,1:min(20, nrow(ctp))]
  ctpca = ctpca[rownames(ctpca) %in% sample_meta$sample,]
  sample_meta_ptime = get_sample_pseudotime(x = t(ctpca), sample_metadata = sample_meta, 
                                            sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                            method = c('cca'))
})


sample_meta_pc = ctpca %>% 
  as.data.frame() %>%
  select(PC1:PC2) %>%
  mutate(sample = rownames(.)) %>%
  inner_join(sample_meta_ptime %>% select(sample, severity, pseudotime, study), by = 'sample')

study.color = c(brewer.pal(5,'Set1'),brewer.pal(7,'Set2')); names(study.color) = sort(unique(sample_meta$study))
severity.color = c("HD" = '#4DAF4A', "Mi" = '#377EB8' ,"Mod" ='orange', "Se" = '#E41A1C',
                   "Rec" = 'grey', "Flu" = 'purple')

theme_paper = function(){
  theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(legend.position = "top", legend.key.size = unit(0.3, "in"),legend.text = element_text(size = 30),legend.title=element_text(size = 30),legend.box = "vertical") + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) +
    theme(axis.text.x = element_text(size=30,color="black"),
          axis.text.y = element_text(size=30,color='black'),
          axis.title.x = element_text(size=40,vjust=-1),
          axis.title.y = element_text(size=40,vjust=1),
          plot.margin=unit(c(1,1,1,1),"cm"))+
    theme(plot.title = element_text(size = 40, hjust = 0.5))
}
create_pseudotime_plot = function(sample_info, title, cluster = F, opt_axis = F, slope = NULL, plot_study = F){
  p1 =  sample_info %>%
    ggplot(aes(x = PC1, y = PC2, color = severity))+
    geom_point(size = 3)+
    coord_fixed(1)+
    scale_colour_manual(values = severity.color)+
    guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
    labs(color = "Severity") +
    ggtitle(title)+
    theme_paper()
  
  if(cluster){
    p1 =  sample_info %>%
      mutate(cluster = factor(cluster)) %>%
      ggplot(aes(x = PC1, y = PC2, color = severity, alpha = cluster))+
      geom_point(size = 3)+
      coord_fixed(1)+
      scale_colour_manual(values = severity.color)+
      guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
      labs(color = "Severity") +
      ggtitle(title)+
      theme_paper()
  }
  
  if(opt_axis){
    p1 =  sample_info %>%
      mutate(cluster = factor(cluster)) %>%
      ggplot(aes(x = PC1, y = PC2, color = severity))+
      geom_point(size = 3)+
      geom_abline(intercept = 0, slope = slope, linetype = 'dashed')+
      coord_fixed(1)+
      scale_colour_manual(values = severity.color)+
      guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
      labs(color = "Severity") +
      ggtitle(title)+
      theme_paper()
  }
  
  if(plot_study){
    p1 =  sample_info %>%
      mutate(cluster = factor(cluster)) %>%
      ggplot(aes(x = PC1, y = PC2, color = study))+
      geom_point(size = 3)+
      geom_abline(intercept = 0, slope = slope, linetype = 'dashed')+
      coord_fixed(1)+
      scale_colour_manual(values = study.color)+
      guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
      labs(color = "Severity") +
      ggtitle(title)+
      theme_paper()
  }
  
  print(p1)
}



pdf('/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/figure/sample_pc.pdf', width =  16, height = 16)
create_pseudotime_plot(sample_meta_pc, "Gene Expression", cluster = F, opt_axis = T, slope = -0.4586)
create_pseudotime_plot(sample_meta_pc, "Gene Expression", cluster = F, opt_axis = T, slope = -0.4586, plot_study = T)
dev.off()