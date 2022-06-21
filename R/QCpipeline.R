#Haley Greenyer 
#Updated: 6/16/2022

#upper pipeline

#' Quality control and DESeq
#'
#' @param gtf_dir Directory where salmon output file folders are located
#' @param salmon_dir Data frame containing sample annotations
#' @param out_dir option to specify directory to which figures and tables will be exported
#' @param cell optional name of cell line to subset
#' @param remove_batch option to remove batch effects (assume replicate=batch)
#' @param grp_by attribute to group by, e.g. 'condition', 'replicate', etc.
#'
#' @return list of DESeqDataSet, VST DESeqDataSet, VST counts, and annotation df 
#' @export
#'
#' @examples QC_DESeq('~/refDirectory', '~/salmonDirectory', cell='Kasumi')
QC_DESeq <- function(gtf_dir, salmon_dir, out_dir=NULL, cell=NULL, remove_batch=TRUE, grp_by='condition'){
  
  #get reference metadata
  tx2gene <- get_ref(gtf_dir=gtf_dir)
  
  s_import <- import_salmon_data(salmon_dir=salmon_dir, tx2gene=tx2gene, out_dir=out_dir)
  
  #get full or sub-setting txi_salmon list and annotation
  if(!is.null(cell)){
    
    txi_subset <- get_cell_subset(cell=cell, txi_salmon=s_import$txi_salmon, out_dir=out_dir)
    annotation <- txi_subset$annotation
    salmon_data <- txi_subset$abundances
  }
  
  else{
    
    annotation <- s_import$annotation
    salmon_data <- s_import$txi_salmon
  }
  
  # run DESeq 
  DESeq_out <- run_DESeq(salmon_data=salmon_data, annotation=annotation, batch_effect=remove_batch, out_dir=out_dir, cell=cell)
  
  vst <- DESeq_out$vst
  vst_mat <- DESeq_out$vst_mat
  dds <- DESeq_out$dds
  
  #Quality Control 
  salmon_QC(salmon_data, dds, vst, annotation, grp_by = grp_by, out_dir=out_dir)
  
  return(list('dds'=dds, 'vst'=vst, 'vst_mat'=vst, 'annotation'=annotation))
  
}