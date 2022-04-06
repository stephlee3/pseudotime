######################
## treecortreat
######################
library(harmony)
library(TreeCorTreat)
library(tidyverse)
library(sva)

source("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/util.R")

#s <- '@A(@B(1,2),@C(3,4))' # intermediate node use '@' to distinguish

cell_meta = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/cell_meta.rds")
raw_count = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/count.rds")

                       
s <- paste0('@All(', 
            '@T(',paste0('T_c',1:12, collapse = ','),'),',
            '@Mono(',paste0('Mono_c',13:23, collapse = ','),'),',
            '@B(',paste0('B',24:28, collapse = ','),'))')

s <- paste0('@All(', 
            '@T(',paste0('T_c',c(1,2,3,4,8), collapse = ','),'),',
            '@Mono(',paste0('Mono_c',c(13,14), collapse = ','),'),',
            '@B(',paste0('B',c(24,25), collapse = ','),'))')

print("Construct hierarchical tree structure...")
pdf('tree.pdf')
hierarchy_list <- extract_hrchy_string(s,'@') # tree is plotted in bottom-up fashion.
dev.off()

hierarchy_list$layout
hierarchy_list$leaves_info

#######################
## column names (must have) - to be consistent in the pipeline
#######################
## sample_meta: sample
## cell_meta: barcode, sample, celltype

######################
## pseudobulk
######################
## extract info
label_info <- hierarchy_list$layout
if(sum(duplicated(label_info$label))) stop('Duplicated cell type name.')
leaves_info <- hierarchy_list$leaves_info

## construct pseudobulk
print("---------------------------------")
print("Construct pseudobulk...")
unq_id <- unique(leaves_info$label)
pb.ls <- lapply(unq_id,function(tid){
  print(tid)
  node_info <- leaves_info %>% filter(label==tid)
  print(node_info)
  # note: cell_meta must contain these columns: barcode, sample, celltype
  sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
  sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]
  
  
  # pseudobulk construction (please fill in)
  pb <- get_sample_pb(s = sub_count, pt = sub_meta$sample, 
                      HVG = T, qnorm = T, combat = T, batch = cell_meta$study)
  print(str(pb))
  return(pb)
})
names(pb.ls) <- unq_id
saveRDS(pb.ls, '/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/pb.ls.rds')


##############
## DEG summary
##############
# format: data frame (column names: label,x,y,id,leaf, (num_DEGs - can use any name you'd like));
# first 5 columns are contained in hierarchy_list$layout