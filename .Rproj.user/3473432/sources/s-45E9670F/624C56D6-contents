library(tidyverse)
library(Matrix)

get_cell_prop_traj = function(cell_type_prop, cluster_list = NULL, celltype_list = NULL, 
                              sample_metadata = NULL, 
                              sample_name = 'Patient', pseudotime_name = 'Pseudotime'){
  celltype_label = paste0('c', cluster_list, '\n', celltype_list)
  patient_order = sample_metadata[, sample_name]
  cp = cell_type_prop[cluster_list, match(patient_order,colnames(cell_type_prop)),drop = F]
  
  cp_df = cp %>%
    as.data.frame() %>%
    mutate(celltype = celltype_label) %>%
    gather(Patient, value, -celltype) %>%
    inner_join(sample_metadata %>% select({{sample_name}}, {{pseudotime_name}}), by = setNames("Patient",sample_name))
  
  
  p = cp_df %>%
    mutate(celltype = factor(celltype, levels = celltype_label)) %>%
    ggplot(aes(x = Pseudotime, y = value))+
    geom_point()+ geom_smooth()+ 
    facet_wrap(vars(celltype),scales = "free_y")+
    ylab("Cell Type Proportion")+
    theme_classic()
  
  return(p)
    
}

## load data
cell_type_prop = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/cell_type_prop.rds")
celltype_info = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/celltype_info.rds")
sample_info_opt = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/sample_info_pb_opt.rds")

## top clusters
topclu = c(1,2,3,4,8,13,14,24,25); names(topclu) = paste0("c",topclu)
topct = celltype_info$celltype[topclu]
p = get_cell_prop_traj(cell_type_prop,cluster_list = topclu, celltype_list = topct,
                       sample_metadata = sample_info_opt)

p = get_cell_prop_traj(cell_type_prop,cluster_list = topclu[1], celltype_list = topct[1],
                       sample_metadata = sample_info_opt)





## DO_NOT_RUN
DO_NOT_RUN = function(){
  
  ## load data
  sample_info = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/lctoptaxis/m1/data/sample_info_pb_opt.rds")
  cell_type_prop = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/qnorm_combat/m1/data/cell_type_prop_m1.rds")
  cell_type_count = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/qnorm_combat/m1/data/cell_type_count_m1.rds")
  patient_info = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/patient_info.rds")
  
  ## 
  sample_info = sample_info %>% filter(Patient %in% patient_info$Patient)
  cell_type_prop = cell_type_prop[,which(colnames(cell_type_prop) %in% patient_info$Patient)]
  cell_type_count = cell_type_count[,which(colnames(cell_type_count) %in% patient_info$Patient)]
  
  
  celltype = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/celltype.rds")
  ctlct_info = data.frame(celltype = celltype, 
                          cluster = names(celltype),
                          lct = gsub('(.*):(.*)','\\2',celltype),
                          label =  paste0(names(celltype),"\n",gsub('(.*)_(.*)', '\\1', gsub('(.*):(.*)','\\1',celltype)))) 
  
  
  saveRDS(sample_info, "/dcl02/hongkai/data/rli/covid/wrapup/data/sample_info_pb_opt.rds")
  saveRDS(cell_type_prop, "/dcl02/hongkai/data/rli/covid/wrapup/data/cell_type_prop.rds")
  saveRDS(cell_type_count, "/dcl02/hongkai/data/rli/covid/wrapup/data/cell_type_count.rds")
  saveRDS(ctlct_info, "/dcl02/hongkai/data/rli/covid/wrapup/data/celltype_info.rds")
  
}