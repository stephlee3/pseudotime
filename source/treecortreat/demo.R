library(tidyverse)
library(TreeCorTreat)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
severity.color = readRDS("/dcl02/hongkai/data/covid/data/current/palette/severity.rds")
study.color = readRDS("/dcl02/hongkai/data/covid/data/current/palette/study.rds")

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



source("/dcl02/hongkai/data/rli/covid/wrapup/demo/util.R")
source("/dcl02/hongkai/data/rli/covid/wrapup/demo/core.R")
source("/dcl02/hongkai/data/rli/covid/wrapup/demo/core_clean.R")

## main
## tree structure
s <- paste0('@All(', 
            '@T(',paste0('T_c',c(1:12), collapse = ','),'),',
            '@Mono(',paste0('Mono_c',c(13:23), collapse = ','),'),',
            '@B(',paste0('B_c',c(24:28), collapse = ','),'),',
            '@Neutrophil(',paste0('Neutrophil_c',c(29:30), collapse = ','),'),',
            '@Megakaryocyte(',paste0('Megakaryocyte_c',c(31), collapse = ','),'))')

pdf('/dcl02/hongkai/data/rli/covid/wrapup/demo/figure/tree_struc.pdf')
hierarchy_list <- extract_hrchy_string(s,'@') # tree is plotted in bottom-up fashion.
dev.off() 
leaves_info <- hierarchy_list$leaves_info
saveRDS(leaves_info, "/dcl02/hongkai/data/rli/covid/wrapup/demo/data/leaves_info.rds")

## sample meta
sample_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/sample_meta.rds")
sample_meta = sample_meta %>%
  filter(severity %in% c('HD','Mi','Mod','Se')) %>%
  mutate(sev.level = case_when(
    severity == 'HD' ~ 1,
    severity == 'Mi' ~ 2,
    severity == 'Mod' ~ 3,
    severity == 'Se' ~ 4
  ))
batch = sample_meta$batch; names(batch) = sample_meta$sample ## use names of batch to indicate sample

set.seed(2022)
sample_meta = sample_meta %>%
  filter(study == 'Su',
         severity %in% c('Mi','Mod')) %>%
  mutate(id = 1:nrow(.)) %>%
  filter(id %in% sample(1:nrow(.), 10))

## cell meta
cell_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/cell_meta.rds")
cell_meta = cell_meta %>%
  mutate(large = gsub('(.*)_(.*)','\\1', celltype)) %>%
  filter(sample %in% sample_meta$sample)



## run pipeline
## list of pseudobulk 
raw_count = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/count.rds")
raw_count = raw_count[,colnames(raw_count) %in% cell_meta$barcode]

ptime_result = runPseudotime(raw_count = raw_count,
                             cell_meta = cell_meta, 
                             sample_meta = sample_meta,
                             leaves_info = leaves_info,
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
                             eval_method = 'corr')

## DEG
sample_meta_ptime = ptime_result$sample_meta_ptime
pb.ls = ptime_result$pb.ls
pb.ls.node = pb.ls$node

deg_result = runPseudoDEG(sample_meta = sample_meta_ptime, pb.ls = pb.ls,
                          leaves_info = leaves_info, phenotype_name = 'severity')
deg_table = deg_result$tab
deg_summary = deg_result$num_DEG

## plot 
create_pseudotime_plot = function(sample_info, rd, title, slope = NULL){
  
  sample_info = rd %>% 
    as.data.frame() %>%
    mutate(sample = rownames(.)) %>%
    inner_join(sample_info, by = 'sample') %>%
    filter(study == 'Su')

  p1 = ggplot(aes(x = PC1, y = PC2, color = severity, label = sample), data= sample_info)+
    geom_point(size = 5)+
    scale_colour_manual(values = severity.color)+
    geom_text_repel()+
    geom_abline(slope = slope, intercept = 0, linetype = 'dashed', col = 'red', size = 3)+
    guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
    labs(color = "Severity") +
    ggtitle(title)+
    theme_paper()
  
  return(p1)
}


ptime.ls = ptime_result$ptime.ls
pdf('/dcl02/hongkai/data/rli/covid/wrapup/demo/figure/pseudotime.pdf', width =  16, height = 16)
for (i in 1:length(ptime.ls)){
  print(i)
  p1 = create_pseudotime_plot(sample_info = ptime.ls[[i]]$ptime$sample_info, rd = ptime.ls[[i]]$rd,
                              title = paste0("CCA", names(ptime.ls)[i]),
                              slope = ptime.ls[[i]]$ptime$slope)
  print(p1)
  
}
dev.off()