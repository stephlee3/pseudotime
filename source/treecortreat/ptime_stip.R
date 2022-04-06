sig = result %>% filter(signif) %>% pull(gene) 
sig = paste0(sig,':c',cluster)

xsig = xsmooth[rownames(xsmooth) %in% sig,]
xsig = t(apply(xsig, 1, scale))
rownames(xsig) = gsub('(.*):(.*)','\\1',rownames(xsig))
colnames(xsig) = colnames(xsmooth)


if (gene_clustering){
  if (length(cluster_num) == 1){
    k = cluster_num
  } else{
    wss = sapply(cluster_num, function(k){
      kmeans(xsig, centers = k)$tot.withinss
    })
    k = find_cluster_num(wss, cluster_num, method = cluster_sel_method)
  }
  gene_clusters = data.frame(gene = rownames(xsig),
                             gene_cluster = kmeans(xsig, centers = k)$cluster)
  result = result %>%
    left_join(gene_clusters, by = 'gene')
}

print(head(result))

if (save){
  #dir.create(file.path(csv_path, "difftest"))
  #write_csv(result, file.path(csv_path, "difftest", paste0("c", cluster,".csv")))
  write_csv(result, file.path(csv_path, filter_method, paste0("opt_kmeans_", cluster_num, "_c", cluster,".csv")))
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
  row_clusters = result$gene_cluster[match(rownames(xsig),result$gene)]
  p = Heatmap(xsig, name = "Expression", 
              show_row_dend = F, cluster_rows = T,
              column_title = "Pseudotime -->", column_title_side = "bottom",
              top_annotation = patient_annot, 
              cluster_columns = F, column_dend_reorder=F,
              row_split = row_clusters)
  xsig_raw = x_reorder[match(paste0(rownames(xsig),":c",cluster), rownames(x_reorder)),]
  xsig_raw = t(apply(xsig_raw, 1, scale))
  rownames(xsig_raw) = gsub('(.*):(.*)','\\1',rownames(xsig_raw))
  colnames(xsig_raw) = colnames(xsig)
  
  gene_order = row_order(p)
  order_res = lapply(1:length(gene_order), function(i){
    genes = gene_order[[i]]
    max_id = sapply(genes, function(i) which.max(xsig[i,]))
    gene_reorder = genes[order(max_id)]
    avg_max_id = mean(max_id)
    return(list(gene_reorder = gene_reorder, avg_max_id = avg_max_id))
  })
  names(order_res) = names(gene_order)
  row_max_id = sapply(order_res, function(i) i$avg_max_id)
  slice_order = order(row_max_id)
  
  gene_reorder = lapply(slice_order, function(i){
    order_res[[i]]$gene_reorder
  })
  gene_reorder = lapply(1:length(gene_reorder), function(i){
    if (i == length(gene_reorder)){
      return(rev(gene_reorder[[i]]))
    } else{
      return(gene_reorder[[i]])
    }
  })
  gene_reorder = do.call(c,gene_reorder)
  
  p = Heatmap(xsig[gene_reorder,], name = "Expression", 
              show_row_dend = F, cluster_rows = F,
              column_title = "Pseudotime -->", column_title_side = "bottom",
              top_annotation = patient_annot, 
              cluster_columns = F, column_dend_reorder=F,
              row_split = factor(row_clusters[gene_reorder], levels = names(gene_order)[slice_order]), row_gap = unit(5, "mm"),
              row_names_gp = gpar(fontsize = 3), row_names_centered = T,
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
  
  
  p_raw = Heatmap(xsig_raw[gene_reorder,], name = "Expression", 
                  show_row_dend = F, cluster_rows = F,
                  column_title = "Pseudotime -->", column_title_side = "bottom",
                  top_annotation = patient_annot, 
                  cluster_columns = F, column_dend_reorder=F,
                  row_split = factor(row_clusters[gene_reorder], levels = names(gene_order)[slice_order]), row_gap = unit(5, "mm"),
                  row_names_gp = gpar(fontsize = 3), row_names_centered = T,
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
  
  
  
  
  pdf(file.path(csv_path, filter_method, paste0("opt_kmeans_", cluster_num, "_c", cluster,".pdf")), width = 32, height = 80)
  draw(p + p_raw)
  dev.off()
}
return(xsig)