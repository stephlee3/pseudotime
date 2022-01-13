library(TSCAN)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

# x: features * samples
## sample_metadata

get_sample_pseudotime = function(x, sample_metadata = NULL, 
                                 sample_name = 'Patient', phenotype_name = 'Phenotype',
                                 method = c('tscan','cca'),
                                 cluster = NULL, reduce = F, clusternum = 2, startcluster = 1,
                                 phenotype_num = NULL){
  method = match.arg(method)
  
  ## match x and meta
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
    res = cancor(pc,phenotype_num)
    proj = (res$xcoef["PC1",1]*(pc[,"PC1"]-res$xcenter["PC1"]))+ (res$xcoef["PC2",1]*(pc[,"PC1"]-res$xcenter["PC2"]))
    pseudotime = order(proj)
    slope = res$xcoef["PC2",1]/res$xcoef["PC1",1]
    print(slope)
    
    sample_info = sample_metadata %>%
      mutate(Pseudotime = pseudotime) %>%
      arrange(Pseudotime)
  }
  
  return(sample_info)
}


## examples
pbpca = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/pbpca.rds")
cell_prop_pca = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/cell_prop_pca.rds")
patient_info = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/patient_info.rds")
sample_info_pb_center = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/sample_info_pb_center.rds")

sev.level = patient_info[match(rownames(pbpca),patient_info$Patient),]$sev.level
precluster = sample_info_pb_center$cluster; names(precluster) = sample_info_pb_center$Patient



## exampel 1: pseudobulk pca + tscan + model-based clustering with 2 clusters
sample_info1 = get_sample_pseudotime(x = t(pbpca[,1:2]), sample_metadata = patient_info, 
                                    sample_name = 'Patient', phenotype_name = 'Phenotype',
                                    method = 'tscan',
                                    reduce = F, clusternum = 2, startcluster = 1)

## exampel 2: pseudobulk pca + tscan + predefined 2 clusters connecting HD center and COVID center
sample_info2 = get_sample_pseudotime(x = t(pbpca[,1:2]), sample_metadata = patient_info, 
                                    sample_name = 'Patient', phenotype_name = 'Phenotype',
                                    method = 'tscan', cluster = precluster,
                                    reduce = F, startcluster = 1)


## exampel 3: cell proportion pca + tscan + model-based clustering with 2 clusters
sample_info3 = get_sample_pseudotime(x = t(cell_prop_pca[,1:2]), sample_metadata = patient_info, 
                                    sample_name = 'Patient', phenotype_name = 'Phenotype',
                                    method = 'tscan',
                                    reduce = F, clusternum = 2, startcluster = 1)

## example 4: pseudobulk pca + cca 
sample_info4 = get_sample_pseudotime(x = t(pbpca[,1:2]), sample_metadata = patient_info, 
                                     sample_name = 'Patient', phenotype_name = 'Phenotype',
                                     method = 'cca', phenotype_num = sev.level)


## plots
theme_paper = function(){
  theme_classic()+
    theme(plot.title = element_text(size = 40, hjust = 0.5)) +
    theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(legend.position = "top", legend.key.size = unit(0.3, "in"),legend.text = element_text(size = 30),legend.title=element_text(size = 30),legend.box = "vertical") + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) +
    theme(axis.text.x = element_text(size=30,color="black"),
          axis.text.y = element_text(size=30,color='black'),
          axis.title.x = element_text(size=40,vjust=-1),
          axis.title.y = element_text(size=40,vjust=1),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
} 

Study.color = c(brewer.pal(5,'Set1'),brewer.pal(7,'Set2')); names(Study.color) = sort(unique(patient_info$Study))
Phenotype.color = c("HD" = '#4DAF4A', "Mi" = '#377EB8' ,"Mod" ='orange', "Se" = '#E41A1C',
                    "Rec" = 'grey', "Flu" = 'purple')

create_pseudotime_plot = function(sample_info, title){
  p1 =  sample_info %>%
    ggplot(aes(x = PC1, y = PC2, color = Phenotype))+
    geom_point(size = 3 )+
    coord_fixed(1)+
    scale_colour_manual(values = Phenotype.color)+
    guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
    labs(color = "Severity") +
    ggtitle(title)+
    theme_paper()
  print(p1)
}


sample_info_pc = pbpca %>% 
  as.data.frame() %>%
  select(PC1:PC2) %>%
  mutate(Patient = rownames(.)) %>%
  inner_join(sample_info4 %>% select(Patient, Phenotype, Pseudotime), by = 'Patient')

create_pseudotime_plot(sample_info_pc, "Gene Expression")



#### DO NOT RUN ###
DO_NOT_RUN = function(){
  pbpca = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/lctqnorm_combat/m1/data/pbpca_m1.rds")
  cell_prop_pca = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/lctqnorm_combat/m1/data/cell_prop_pca_m1.rds")
  sample_info_pb_center = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/lctqnorm_combat/m1/data/sample_info_pb_center_m1.rds")
  
  
  ## 400 samples after some filtering
  patient_info = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/lctqnorm_combat/m1/data/patient_info_m1.rds")
  patient_info = patient_info %>% filter(Patient %in% rownames(pbpca),
                                         Phenotype %in% c("HD","Mi","Mod","Se"))
  
  
  sev.level = patient_info$Phenotype
  sev.level = case_when(
    sev.level == "HD" ~ 1,
    sev.level == "Mi" ~ 2,
    sev.level == "Mod" ~3,
    sev.level == "Se"~ 4
  )
  patient_info$sev.level = sev.level
  
  pbpca = pbpca[rownames(pbpca) %in% patient_info$Patient,]
  cell_prop_pca = cell_prop_pca[rownames(cell_prop_pca) %in% patient_info$Patient,]
  
  
  sample_info_pb_center = sample_info_pb_center %>% filter(Patient %in% patient_info$Patient)
  sample_info_pb_center = sample_info_pb_center[match(rownames(pbpca), sample_info_pb_center$Patient),]
  
  saveRDS(pbpca, "/dcl02/hongkai/data/rli/covid/wrapup/data/pbpca.rds")
  saveRDS(cell_prop_pca, "/dcl02/hongkai/data/rli/covid/wrapup/data/cell_prop_pca.rds")
  saveRDS(patient_info, "/dcl02/hongkai/data/rli/covid/wrapup/data/patient_info.rds")
  saveRDS(sample_info_pb_center,"/dcl02/hongkai/data/rli/covid/wrapup/data/sample_info_pb_center.rds")
  
}

