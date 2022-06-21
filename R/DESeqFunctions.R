#Haley Greenyer 
#Updated: 6/15/2022


run_DESeq <- function(salmon_data, annotation, cell=NULL, batch_effect=TRUE, out_dir=NULL) UseMethod('run_DESeq', salmon_data)


#' DESeq with double input 
#'
#' @param salmon_data count matrix 
#' @param annotation sample annotation df
#' @param cell option to specify cell line 
#' @param batch_effect option to remove batch effects (assume replicate=batch)
#' @param out_dir option to specify output directory 
#'
#' @return list containing DESeqDataSet object, VST version of theDESeqDataSet, and VST counts
#' @export
#'
#' @examples run_DESeq(counts, annotation_df)
setMethod("run_DESeq",
          signature(salmon_data = "double"),
          function (salmon_data, annotation, cell=NULL, batch_effect=TRUE, out_dir=NULL) 
          {
            #remove batch effects
            if(batch_effect){
              dds <- DESeqDataSetFromMatrix(countData = round(salmon_data), 
                                            colData = annotation, 
                                            design = ~replicate + condition)
              dds <- DESeq(dds)
              vst <- varianceStabilizingTransformation(dds, blind=TRUE)
              #assume replicate=batch
              assay(vst) <- limma::removeBatchEffect(assay(vst), vst$replicate)
              vst_mat <- assay(vst) 
            }
            
            
            else{
              dds <- DESeqDataSetFromMatrix(countData = round(salmon_data), 
                                            colData = annotation, 
                                            design = ~condition)
              dds <- DESeq(dds)
              vst <- varianceStabilizingTransformation(dds, blind=TRUE)
              vst_mat <- assay(vst) 
            }
            
            #option to return annotation and VST counts 
            if(!is.null(out_dir)){
              filepath1 = paste(out_dir,'/',cell,'Annotation.csv',sep='')
              filepath2 = paste(out_dir,'/',cell,'VSTcounts.csv',sep='')
              write(as.matrix(annotation), file=filepath1)
              write(as.matrix(vst_mat), file=filepath2)
            }
            
            return(list('vst'=vst,'vst_mat'=vst_mat, 'dds'=dds))
          }
)


#' DESeq with matrix input 
#'
#' @param salmon_data count matrix 
#' @param annotation sample annotation df
#' @param cell option to specify cell line 
#' @param batch_effect option to remove batch effects (assume replicate=batch)
#' @param out_dir option to specify output directory 
#'
#' @return list containing DESeqDataSet object, VST version of theDESeqDataSet, and VST counts
#' @export
#'
#' @examples run_DESeq(count_matrix, annotation_df)
setMethod("run_DESeq",
          signature(salmon_data = "matrix"),
          function (salmon_data, annotation, cell=NULL, batch_effect=TRUE, out_dir=NULL) 
          {
            #remove batch effects
            if(batch_effect){
              dds <- DESeqDataSetFromMatrix(countData = round(salmon_data), 
                                            colData = annotation, 
                                            design = ~replicate + condition)
              dds <- DESeq(dds)
              vst <- varianceStabilizingTransformation(dds, blind=TRUE)
              #assume replicate=batch
              assay(vst) <- limma::removeBatchEffect(assay(vst), vst$replicate)
              vst_mat <- assay(vst) 
            }
            
            
            else{
              dds <- DESeqDataSetFromMatrix(countData = round(salmon_data), 
                                            colData = annotation, 
                                            design = ~condition)
              dds <- DESeq(dds)
              vst <- varianceStabilizingTransformation(dds, blind=TRUE)
              vst_mat <- assay(vst) 
            }
            
            #option to return annotation and VST counts 
            if(!is.null(out_dir)){
              filepath1 = paste(out_dir,'/',cell,'Annotation.csv',sep='')
              filepath2 = paste(out_dir,'/',cell,'VSTcounts.csv',sep='')
              write(as.matrix(annotation), file=filepath1)
              write(as.matrix(vst_mat), file=filepath2)
            }
            
            
            return(list('vst'=vst,'vst_mat'=vst_mat, 'dds'=dds))
          }
)



#' DESeq with list input 
#'
#' @param salmon_data salmon tximport list 
#' @param annotation sample annotation df
#' @param cell option to specify cell line 
#' @param batch_effect option to remove batch effects (assume replicate=batch)
#' @param out_dir option to specify output directory 
#'
#' @return list containing DESeqDataSet object, VST version of theDESeqDataSet, and VST counts
#' @export
#'
#' @examples run_DESeq(txi_list, annotation_df)
setMethod("run_DESeq",
          signature(salmon_data = "list"),
          function (salmon_data, annotation, cell=NULL, batch_effect=TRUE, out_dir=NULL) 
          {
            
            #remove batch effects 
            if(batch_effect){
              dds <- DESeqDataSetFromTximport(txi = salmon_data, 
                                              colData = annotation, 
                                              design = ~replicate + condition)
              dds <- DESeq(dds)
              vst <- varianceStabilizingTransformation(dds, blind=TRUE)
              #assume replicate=batch
              assay(vst) <- limma::removeBatchEffect(assay(vst), vst$replicate)
              vst_mat <- assay(vst) 
              
              print(type(annotation))
              print(type(vst_mat))
            }
            
            
            else{
              dds <- DESeqDataSetFromTximport(txi = salmon_data, 
                                              colData = annotation, 
                                              design = ~condition)
              dds <- DESeq(dds)
              vst <- varianceStabilizingTransformation(dds, blind=TRUE)
              vst_mat <- assay(vst) 
            }
            
            #option to return annotation and VST counts 
            if(!is.null(out_dir)){
              filepath1 = paste(out_dir,'/',cell,'Annotation.csv',sep='')
              filepath2 = paste(out_dir,'/',cell,'VSTcounts.csv',sep='')
              write(as.matrix(annotation), file=filepath1)
              write(as.matrix(vst_mat), file=filepath2)
            }
            
            return(list('vst'=vst,'vst_mat'=vst_mat, 'dds'=dds))
          }
)


