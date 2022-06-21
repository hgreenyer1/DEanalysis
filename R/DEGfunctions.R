#Haley Greenyer 
#Updated: 6/16/2022

#' Get DEGs for one Condition
#'
#' @param dds DESeqDataSet object
#' @param c1 condition 1, e.g. 'control'
#' @param c2 condition 2, e.g. 'dTAG'
#' @param organism e.g. 'hsapiens', for conversion of ENSG IDs to gene symbols
#' @param cell option to specify cell line for output naming 
#' @param out_dir option to specify output directory
#'
#' @return list of 2, DE matrix and table of DEGs
#' @export 
#'
#' @examples get_condition_DEGs(kasumi_dds_object, 'control', 'dTAG', cell='Kasumi' )
get_condition_DEGs <- function(dds, c1, c2, organism='hsapiens', cell=NULL, out_dir=NULL){
  
  resA <- results(dds, contrast=c('condition',c1,c2))
  summary_deseq(resA)
  resA <- resA[order(resA$padj), ]
  resA_data <- merge(as.data.frame(resA), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resA_data)[1] <- "Gene"
  
  #remove the version id
  resA_data2 <- resA_data %>%  
    mutate(Gene_symbol = str_replace(Gene, ".[0-9]+$", ""))
  colnames(resA_data2)
  resA_symbol <- gconvert(resA_data2$Gene, organism = "hsapiens", target = "ENSG", numeric_ns = "", mthreshold = Inf, filter_na = T)
  names(resA_symbol)[2] <- "Gene"
  data <- merge(resA_data2, resA_symbol, by="Gene", all.x=T, sort=F)
  
  c1_v_c2 <- data %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
  tableDEGs <- c1_v_c2[order(c1_v_c2$log2FoldChange), ] %>% dplyr::select("Gene_symbol", "name", "log2FoldChange", "padj")
  rownames(tableDEGs) = NULL
  
  #for title spacing
  cell <- paste(' ',cell, sep = '')

  DEG_tbl <- knitr::kable(tableDEGs, caption = paste(c1,' vs ', c2, ' DEGs', sep = ''), digits = 3) %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    scroll_box(width = "500px", height = "200px")
  
  show(DEG_tbl)
  
  vsVolcano(
    x = c2, y = c1, 
    data = dds, d.factor = 'condition', type = 'deseq', 
    padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
    legend = TRUE, grid = TRUE, data.return = FALSE)
   
  #export DEGs and figures
  if(!is.null(out_dir)){
    write_csv(c1_v_c2, file = paste(out_dir, '/', cell, c1, '_v_', c2,'.csv',sep=''))
    write_csv(tableDEGs, file = paste(out_dir, '/', cell, c1, '_v_', c2,'DEGtable.csv',sep=''))
    
    pdf(file=paste(out_dir,'/', cell, c1,'_',c2,'_volcano.pdf',sep=''))
    vsVolcano(
      x = c2, y = c1, 
      data = dds, d.factor = 'condition', type = 'deseq', 
      padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
      legend = TRUE, grid = TRUE, data.return = FALSE)
    dev.off()
  }

    
  return(list('DE_df'=c1_v_c2, 'DEG_table'=tableDEGs))
  
}

#function compare DEGs: compare get_DEGs output across conditions
#args to_compare: list containing DE_dfs produced by get_DEGs
#
#' Get DEG lists
#'
#' @param to_compare list of DE_df objects to compare (produced by get_DEGs) 
#' @param vst_mat matrix of VST counts 
#' @param annotation sample annotation df for the appropriate cell line 
#' @param out_dir option to specify output directory 
#' @param cell optional specification of cell line 
#'
#' @return list of DE dfs
#' @export
#'
#' @examples get_DEG_lists(control_dTAG, NO_NOdTAG,...)
compare_DEGs <- function(to_compare, vst_mat, annotation, out_dir=NULL, cell=NULL){
  
  compare_list <- list()
  for_heatmap <- NULL
  #get names for visualization 
  for(df in to_compare){
    #DEG list
    df_deg <- as_tibble(df) %>% 
      magrittr::use_series(Gene) %>% 
      as.character()
    compare_list <- c(compare_list, list(name=df_deg))
    for_heatmap <- c(for_heatmap,df_deg)
  }
  #set names
  names(compare_list) <- names(to_compare)

  #visualize
  show(upset(fromList(compare_list)))
  
  #show(ggVennDiagram(compare_list))
  
  #export upset and comparison list 
  if(!is.null(out_dir)){
    write.xlsx(compare_list, file = paste(out_dir, '/', cell, 'DEG_lists.xlsx', sep=''))
    pdf(file=paste(out_dir,'/', cell, 'DEGupset.pdf',sep=''))
    upset(fromList(compare_list))
    dev.off()
    pdf(file=paste(out_dir,'/', cell, 'DEGvennDiagram.pdf',sep=''))
    ggVennDiagram(compare_list)
    dev.off()
    
  }

  #get heatmaps 
  #DEG_heatmaps(for_heatmap, vst_mat = vst_mat, out_dir=out_dir, cell=cell, annotation=annotation )
  
  return(list('cl'=compare_list, 'heat'=for_heatmap))
}


#' DEG Heatmaps
#'
#' @param DEG_lists concatenated lits of DEG_df data
#' @param vst_mat matrix of VST counts 
#' @param annotation sample annotation df for the appropriate cell line  
#' @param out_dir option to specify output directory 
#' @param cell option to specify cell line 
#'
#' @return heatmaps
#' @export
#'
#' @examples DEG_heatmaps(DEG_compare_list, vst_counts)
DEG_heatmaps <- function(DEG_lists, vst_mat, annotation, out_dir=NULL, cell=NULL){
  
  dTAG_degs = unique(DEG_lists)
  
  length(dTAG_degs)
  
  #heatmaps of all DEGs
  breaksList = seq(-2, 2, by = 0.05)
  
  DEGZ = as.character(dTAG_degs)
  
  #with clustering
  DEGs_vst <- vst_mat[intersect(row.names(vst_mat), DEGZ),]
  heat_cluster <- pheatmap(DEGs_vst, border_color = NA, scale = "row", 
           cluster_rows = TRUE, cluster_cols = TRUE, 
           clustering_distance_rows = "correlation", 
           clustering_distance_cols = "correlation",
           annotation_col = annotation,  main = "dTAG DEGs", 
           breaks=breaksList, show_rownames = FALSE,
           fontsize_row = 7)
  
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell, 'Heatmap_cluster.pdf',sep=''))
    show(heat_cluster)
    dev.off()
  }
  else{
    show(heat_cluster)
  }

  
  #without clustering
  heat_noCluster <- pheatmap(DEGs_vst,  scale = "row", 
           cluster_rows = TRUE, cluster_cols = FALSE, 
           clustering_distance_rows = "correlation", 
           clustering_distance_cols = "correlation",
           annotation_col = annotation,  main = "dTAG DEGs no column clustering", 
           breaks=breaksList, show_rownames = FALSE,
           fontsize_row = 7)
  
  if(!is.null(out_dir)){
    pdf(file=paste(out_dir,'/', cell, 'Heatmap_NOcluster.pdf',sep=''))
    show(heat_noCluster)
    dev.off()
  }
  else{
    show(heat_noCluster)
  }
  
  
}



