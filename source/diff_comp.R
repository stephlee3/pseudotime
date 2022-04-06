library(tidyverse)
topclu = c(1,2,3,4,8,13,14,24,25)
cv_result = lapply(topclu, function(i){
  result = read_csv(paste0("/dcl02/hongkai/data/rli/covid/wrapup/data/diff0201/cv/c",
                          i, '.csv'))
})
cv_result = do.call(rbind, cv_result) 

eff_result = lapply(topclu, function(i){
  result = read_csv(paste0("/dcl02/hongkai/data/rli/covid/wrapup/data/diff0201/postcurve/c",
                           i, '.csv'))
})
eff_result = do.call(rbind, eff_result) 

result = lapply(topclu, function(i){
  cv_gene = cv_result %>% filter(cell_cluster == i & signif) %>% pull(gene)
  eff_gene = eff_result %>% filter(cell_cluster == i & signif) %>% pull(gene)
  res = data.frame(cv_num = length(cv_gene),
             eff_num = length(eff_gene),
             common_num = length(intersect(cv_gene, eff_gene)),
             cell_cluster = paste0('c',i))
  res
}) 
result = do.call(rbind, result)