######################
## treecortreat
######################
library(harmony)
library(TreeCorTreat)
library(tidyverse)
s <- '@A(@B(1,2),@C(3,4))' # intermediate node use '@' to distinguish
hierarchy_list <- extract_hrchy_string(s,'@') # tree is plotted in bottom-up fashion.
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
unq_id <- unique(leaves_info$label)
pb.ls <- lapply(unq_id,function(tid){
    print(tid)
    node_info <- leaves_info %>% filter(label==tid)
    print(node_info)
    # note: cell_meta must contain these columns: barcode, sample, celltype
    #sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)

    # pseudobulk construction (please fill in)
    #pb <- {}
})
names(pb.ls) <- unq_id


##############
## DEG summary
##############
# format: data frame (column names: label,x,y,id,leaf, (num_DEGs - can use any name you'd like));
# first 5 columns are contained in hierarchy_list$layout



