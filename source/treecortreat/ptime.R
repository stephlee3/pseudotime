library(TSCAN)
library(tidyverse)
library(Matrix)
library(sva)
library(preprocessCore)


source("/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/util.R")


# `x`: the data matrix of features $\times$ samples, i.e. each row is a feature of gene expression or cell type proportion, each column is a patient sample.

# `sample_metadata`: patient-level metadata, including several important variables.

# `sample_name` and `severity_name`: the variable name in `sample_metadata` that indicates the sample and severity. The default setting of `sample_name` is "Patient" and `severity_name` is "severity". 

# `method`: users now have two options, either using TSCAN for sample-level pseudotime construction or use CCA for optimal severity axis detection. 

# `cluster = NULL, reduce = F, clusternum = 2, startcluster = 1`: parameters for TSCAN.

# `severity_num_name`: the variable name in 'sample_metadata' that indicates the numerical valaue of severity severity. 

get_sample_pseudotime = function(x, sample_metadata = NULL, 
                                 sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                 method = c('tscan','cca'),
                                 cluster = NULL, reduce = F, clusternum = 2, startcluster = 1){
  method = match.arg(method)
  
  ## arrange sample metadata to have the same order as x
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
    severity_num = sample_metadata[,severity_num_name]
    res = cancor(pc,severity_num)
    proj = (res$xcoef["PC1",1]*(pc[,"PC1"]-res$xcenter["PC1"]))+ (res$xcoef["PC2",1]*(pc[,"PC2"]-res$xcenter["PC2"]))  
    pseudotime = rank(proj)
    if (cor(pseudotime, severity_num)<0)
      pseudotime = length(pseudotime)+1 -pseudotime
    slope = res$xcoef["PC2",1]/res$xcoef["PC1",1]
    print(paste("slope is", slope))
    print(paste("corr is ", res$cor))
    
    sample_info = sample_metadata %>%
      mutate(pseudotime = pseudotime) %>%
      arrange(pseudotime)
  }
  
  return(sample_info)
}

## main ##
largepb = readRDS('/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/largepb.rds')
sample_meta = sample_meta %>%
  filter(severity %in% c('HD','Mi','Mod','Se')) %>%
  filter(sample %in% colnames(largepb)) %>%
  mutate(sev.level = case_when(
    severity == 'HD' ~ 1,
    severity == 'Mi' ~ 2,
    severity == 'Mod' ~ 3,
    severity == 'Se' ~ 4
  ))
pbpca = prcomp(t(largepb))$x[,1:20]
pbpca = pbpca[rownames(pbpca) %in% sample_meta$sample,]
sample_meta_ptime = get_sample_pseudotime(x = t(pbpca), sample_metadata = sample_meta, 
                                 sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                 method = c('cca'))

print(head(sample_meta_ptime))

saveRDS(sample_meta_ptime,"/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/sample_meta_ptime.rds")
