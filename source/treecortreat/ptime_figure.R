library(tidyverse)
library(ggpubr)
library(RColorBrewer)

sample_meta_pc = pbpca %>% 
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
