library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(mgcv)
library(parallel)

sample_info_opt = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/sample_info_pb_opt.rds")
pb_raw = readRDS("/dcl02/hongkai/data/rli/covid/wrapup/data/pb_no_scale.rds")
common_pt = intersect(colnames(pb_raw), sample_info_opt$Patient)

sample_info = sample_info_opt %>%
  filter(Patient %in% common_pt)

pb_raw = pb_raw[,colnames(pb_raw) %in% common_pt]
pb_reorder = pb_raw[,match(common_pt, colnames(pb_raw))]

example_genes = c('GZMA','RAP1GAP2','FXYD7','KLRB1','RPS27','EEF1A1')
pb_example = pb_reorder[paste0(example_genes,':c1'),]
saveRDS(pb_example, "/dcl02/hongkai/data/rli/covid/wrapup/data/pb_example_c1.rds")


get_gene_expr_traj = function(x, gene_list, cluster, sample_metadata, 
                              pseudotime_name = 'Pseudotime'){
  ptime = sample_metadata[,pseudotime_name]
  sample_metadata = sample_metadata[order(ptime),]
  xg = x[paste0(gene_list,":c",cluster),,drop = F]
  colnames(xg) = 1:ncol(xg)
  df = xg %>%
    as.data.frame() %>%
    mutate(gene=rownames(.)) %>%
    pivot_longer(!gene, names_to = 'Pseudotime', values_to = 'Expression') %>%
    mutate(Pseudotime = as.numeric(Pseudotime))
  str(df)
  p = df%>%
    ggplot(aes(x = Pseudotime, y = Expression))+
    geom_point()+ geom_smooth()+
    facet_wrap(vars(gene)) + theme_classic()
  return(p)
}
p = get_gene_expr_traj(x = pb_reorder, 
                       gene_list = c('ATXN1','ADGRG1','C18orf25','CCDC9','ACVR2A'), 
                       cluster = 1, 
                       sample_metadata = sample_info)
#print(p)

run_pseudo_diff_gene = function(x, sample_metadata = NULL, cluster = 1,
                         sample_name = 'Patient', pseudotime_name = 'Pseudotime',
                         save = F, csv_path = './data', plot = F
                         ){
  print(cluster)
  ptime = sample_metadata[,pseudotime_name]
  sample_metadata = sample_metadata[order(ptime),]
  
  fn = rownames(x)
  gn = gsub('(.*):(.*)','\\1',fn)
  cn = gsub('(.*):(.*)','\\2',fn)
  
  xc = x[which(cn == paste0("c",cluster)),]
  xc_reorder = xc[,match(sample_metadata$Patient,colnames(x))]
  
  res = lapply(1:nrow(xc_reorder), function(i){
    model = mgcv::gam(xc_reorder[i,]~s(c(1:ncol(xc_reorder)),k=3))
    smooth = model$fitted.values
    pval = pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F) 
    list(smooth = smooth, pval = pval)
  })
  pval = sapply(res, function(i) i$pval)
  xsmooth = lapply(res, function(i) i$smooth) %>% do.call(rbind, .)
  rownames(xsmooth) = rownames(xc_reorder)
  print(str(xsmooth))
  
  fdr = p.adjust(pval, method = 'fdr')
  result = data.frame(gene = gsub('(.*):(.*)','\\1',rownames(xc_reorder)), pval = pval, fdr = fdr) %>%
    mutate(cluster = cluster,
           signif = (fdr < 0.05)) %>%
    arrange(fdr)
  
  sig = result %>% filter(signif) %>% pull(gene) 
  sig = paste0(sig,':c',cluster)
  
  xsig = xsmooth[rownames(xsmooth) %in% sig,]
  xsig = t(apply(xsig, 1, scale))
  rownames(xsig) = gsub('(.*):(.*)','\\1',rownames(xsig))
  colnames(xsig) = colnames(xsmooth)
  
  
  if (save){
    dir.create(file.path(csv_path, "diff"))
    write_csv(result, file.path(csv_path, "diff", paste0("c", cluster,".csv")))
  }
  if (plot){
    Study.color = readRDS("/dcl02/hongkai/data/covid/data/current/palette/study.rds")
    Phenotype.color = readRDS("/dcl02/hongkai/data/covid/data/current/palette/severity.rds")
    patient_annot = HeatmapAnnotation(Severity = sample_metadata$Phenotype,
                                      Study = sample_metadata$Study,
                                      #Sex = sample_info$Sex,
                                      #Age = sample_info$Age,
                                      col = list(Severity = Phenotype.color,
                                                 Study = Study.color),
                                      #Sex = sex.color,
                                      #Age = col_fun),
                                      show_legend = F, simple_anno_size = unit(10,"mm"),
                                      annotation_name_gp = gpar(fontsize = 25)
    )
    p = Heatmap(xsig, name = "Expression", 
                show_row_dend = F, cluster_rows = T,
                column_title = "Pseudotime -->", column_title_side = "bottom",
                top_annotation = patient_annot, 
                cluster_columns = F, column_dend_reorder=F,
                row_names_gp = gpar(fontsize = 1), row_names_centered = T,
                column_names_gp = gpar(fontsize = 10),
                column_title_gp = gpar(fontsize = 40),
                show_heatmap_legend = F,
                show_column_names = F,
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 25),
                  labels_gp = gpar(fontsize = 25),
                  legend_height = unit(8, "cm"),
                  grid_width = unit(1,"cm")
                ),
                use_raster = T, raster_device = "tiff")
    pdf(file.path(csv_path, "diff", paste0("c", cluster,".pdf")), width = 16, height = 80)
    print(p)
    dev.off()
  }
  return(xsig)
}

topclu = c(1,2,3,4,8,13,14,24,25); 
diff_result = mclapply(topclu, run_pseudo_diff_gene, mc.cores = 9,
                     x = pb_reorder, sample_metadata = sample_info,
                     sample_name = 'Patient', pseudotime_name = 'Pseudotime',
                     save = T, csv_path = '/dcl02/hongkai/data/rli/covid/wrapup/data', plot = T)

saveRDS(diff_result, "/dcl02/hongkai/data/rli/covid/wrapup/data/diff_smooth.rds")


## DO NOT RUN ## 
DO_NOT_RUN = function(){
  pb = readRDS("/dcl02/hongkai/data/rli/covid/nb1221/qnorm_combat/m1/data/pbhv_m1.rds")
  saveRDS(pb,"/dcl02/hongkai/data/rli/covid/wrapup/data/pb.rds" )
}

## STIP
STIP_ht <- function(fit,fit_raw, scale = F, patient_annot_l = NULL, patient_annot_r = NULL,plot.title = "Cytokine/Receptor Signatures", font_size = 10, row_names = F, 
                    row_annot = NULL, row_annot_size = 25, row_annot_col = "black", link_width = 20, direction = "all") {
  fit.colname = colnames(fit)
  colnames(fit) = NULL
  
  dn <- dimnames(fit)
  if (scale == T) {
    fit <- t(apply(fit,1,scale))
  }
  dimnames(fit) <- dn
  gene <- row.names(fit)
  
  zpdirection <- fit[,1] < fit[,ncol(fit)]
  
  zp <- apply(fit,1,function(sf) {
    which(sapply(1:(length(sf)-1),function(i) sf[i]*sf[i+1] < 0))
  })
  zpnum <- sapply(zp,length)
  inczp <- names(which(zpdirection[zpnum==1]))
  deczp <- names(which(!zpdirection[zpnum==1]))
  multipoint <- names(zpnum)[zpnum > 1]
  
  geneorder <- NULL
  if (length(deczp) > 0) {
    geneorder <- c(geneorder,names(sort(unlist(zp[deczp]),decreasing = F)))
  }
  if (length(multipoint) > 0) {
    geneorder <- c(geneorder,names(sort(sapply(zp[multipoint],function(i) i[1]))))
  }
  if (length(inczp) > 0) {
    geneorder <- c(geneorder,names(sort(unlist(zp[inczp]))))
  }
  
  if (direction == "dec"){
    geneorder = geneorder[1:length(deczp)]
  } else if (direction == "inc"){
    geneorder = geneorder[tail(1:length(geneorder), length(inczp))]
  }
  
  fit_reorder =  fit[geneorder,]
  colnames(fit_reorder) = fit.colname
  print(str(fit_reorder))
  
  
  fit_raw = fit_raw[geneorder,]
  print(str(fit_raw))
  
  ## create heatmap
  p = Heatmap(fit_reorder, name = "Smooth", 
              show_row_names = row_names, show_row_dend = F, cluster_rows = F,
              column_title = "Pseudotime -->", column_title_side = "bottom",
              top_annotation = patient_annot_l, 
              cluster_columns = F, column_dend_reorder=F,
              row_names_gp = gpar(fontsize = font_size), row_names_centered = T,
              column_names_gp = gpar(fontsize = 10),
              column_title_gp = gpar(fontsize = 40),
              show_heatmap_legend = F,
              show_column_names = F,
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 25),
                labels_gp = gpar(fontsize = 25),
                legend_height = unit(8, "cm"),
                grid_width = unit(1,"cm")
              ),
              use_raster = T, raster_device = "tiff")
  
  p_raw = Heatmap(fit_raw, name = "Raw", 
                  show_row_names = row_names, show_row_dend = F, cluster_rows = F,
                  column_title = "Pseudotime -->", column_title_side = "bottom",
                  top_annotation = patient_annot_r, 
                  cluster_columns = F, column_dend_reorder=F,
                  row_names_gp = gpar(fontsize = font_size), row_names_centered = T,
                  column_names_gp = gpar(fontsize = 10), 
                  column_title_gp = gpar(fontsize = 40),
                  show_column_names = F,
                  show_heatmap_legend = F,
                  heatmap_legend_param = list(
                    title_gp = gpar(fontsize = 25),
                    labels_gp = gpar(fontsize = 25),
                    legend_height = unit(8, "cm"),
                    grid_width = unit(1,"cm")
                  ),
                  use_raster = T, raster_device = "tiff")
  
  if (!is.null(row_annot)){
    ha = rowAnnotation(gene = anno_mark(at = match(row_annot,geneorder), labels = row_annot, labels_gp = gpar(fontsize = row_annot_size, col = row_annot_col),link_width = unit(link_width, "mm")))
    p_raw = Heatmap(fit_raw, name = "Raw", 
                    show_row_names = F, show_row_dend = F, cluster_rows = F,
                    right_annotation = ha,
                    column_title = "Pseudotime -->", column_title_side = "bottom",
                    top_annotation = patient_annot_r, 
                    cluster_columns = F, column_dend_reorder=F,
                    row_names_gp = gpar(fontsize = font_size), row_names_centered = T,
                    column_names_gp = gpar(fontsize = 10),
                    column_title_gp = gpar(fontsize = 40),
                    show_heatmap_legend = F,
                    show_column_names = F,
                    heatmap_legend_param = list(
                      title_gp = gpar(fontsize = 25),
                      labels_gp = gpar(fontsize = 25),
                      legend_height = unit(8, "cm"),
                      grid_width = unit(1,"cm")
                    ),
                    use_raster = T, raster_device = "tiff")
  }
  
  draw(p + p_raw, column_title = plot.title, column_title_gp = gpar(fontsize = 20))
  return(geneorder)
  
}
