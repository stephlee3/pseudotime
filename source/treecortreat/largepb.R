
library(tidyverse)
library(Matrix)
library(sva)
library(preprocessCore)



source("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/util.R")

cell_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/cell_meta.rds")
raw_count = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/count.rds")
sample_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/sample_meta.rds")

cell_meta = cell_meta %>%
  mutate(large = gsub('(.*)_(.*)','\\1',celltype))
largect = cell_meta$large



batch = sample_meta$study; names(batch) = sample_meta$sample
largepb = get_cluster_pb(s = raw_count, clu = cell_meta$large, pt = cell_meta$sample,
                         topclu = c('T','Mono','B'), batch = batch, 
                         HVG = T, qnorm = T, combat = T)
saveRDS(largepb,'/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/largepb.rds')